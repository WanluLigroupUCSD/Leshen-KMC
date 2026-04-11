"""
SuperCache — manages multiple superbasins with dynamic barrier tolerance.

Ported from openFLY kinetic/cache.hpp + cache.cpp.

Decides whether to expand the active superbasin (low barrier transition)
or create a new one (high barrier transition). Caches inactive superbasins
for reuse, and dynamically adjusts the barrier tolerance threshold.
"""

import numpy as np
from .basin import Basin
from .superbasin import SuperBasin


class SuperCache:
    """
    SuperBasin cache with dynamic barrier tolerance.

    Manages the active superbasin and a cache of inactive ones.
    When following a low-barrier mechanism, the destination basin is
    added to the current superbasin. When following a high-barrier
    mechanism, the current superbasin is cached and a new one is
    created (or a matching cached one is restored).

    Parameters
    ----------
    temperature : float
        Temperature in Kelvin.
    barrier_tol : float
        Barrier threshold in eV. Mechanisms below this are "low barrier"
        and trigger superbasin expansion.
    state_tol : float
        Position tolerance for basin state matching.
    cache_size : int
        Maximum number of cached superbasins.
    max_superbasin_size : int
        Maximum number of basins per superbasin before shrinking tolerance.
    dynamic_tol : bool
        If True, automatically adjust barrier_tol.
    tol_grow : float
        Factor to grow barrier_tol when cache is frequently reused.
    tol_shrink : float
        Factor to shrink barrier_tol when superbasin overflows.
    max_barrier : float
        Maximum barrier for basin construction.
    debug : bool
    """

    def __init__(self, temperature, barrier_tol=0.3, state_tol=0.5,
                 cache_size=10, max_superbasin_size=20,
                 dynamic_tol=True, tol_grow=1.5, tol_shrink=0.7,
                 max_barrier=5.0, debug=False):
        self.temperature = temperature
        self.barrier_tol = barrier_tol
        self.state_tol = state_tol
        self.cache_size = cache_size
        self.max_superbasin_size = max_superbasin_size
        self.dynamic_tol = dynamic_tol
        self.tol_grow = tol_grow
        self.tol_shrink = tol_shrink
        self.max_barrier = max_barrier
        self.debug = debug

        self._sb = None  # Active SuperBasin
        self._cache = []  # List of cached SuperBasins
        self._in_cache_count = 0

    @property
    def active_superbasin(self):
        return self._sb

    def initialize(self, positions, catalogue):
        """
        Initialize with the first basin.

        Parameters
        ----------
        positions : ndarray, shape (n_atoms, 3)
        catalogue : Catalogue
        """
        basin = Basin(positions, catalogue, self.temperature,
                      max_barrier=self.max_barrier, debug=self.debug)
        self._sb = SuperBasin(basin, debug=self.debug)

    def kmc_choice(self, rng=None):
        """
        Select a mechanism from the active superbasin.

        Returns (mechanism, atom_index, dt, source_basin, switched).
        """
        if self._sb is None:
            raise RuntimeError("SuperCache not initialized")
        return self._sb.kmc_choice(rng)

    def connect_from(self, source_basin, atom_index, mechanism,
                     positions, catalogue):
        """
        Process a completed mechanism transition.

        Determines whether the destination is:
          1. An existing basin in the current superbasin (internal jump)
          2. A new low-barrier destination (expand superbasin)
          3. A high-barrier destination (cache current, start/restore new)

        Parameters
        ----------
        source_basin : int
            Source basin index in the superbasin.
        atom_index : int
        mechanism : Mechanism
        positions : ndarray, shape (n_atoms, 3)
            Current (relaxed) positions after the transition.
        catalogue : Catalogue
        """
        # Compute state hash for the new configuration
        from hashlib import sha256
        n_atoms = len(positions)
        hash_indices = []
        for i in range(n_atoms):
            entry = catalogue.get_entry(i)
            if entry is not None:
                hash_indices.append(entry.cat_index)
            else:
                hash_indices.append(-1)
        h = sha256(np.array(hash_indices, dtype=np.int32).tobytes())
        state_hash = int(h.hexdigest()[:16], 16)

        # Case 1: Existing basin in current superbasin
        prev = self._sb.find_occupy(state_hash, positions, self.state_tol)
        if prev is not None:
            self._sb.connect_from(source_basin, atom_index, mechanism)
            if self.debug:
                print(f"SuperCache: Internal jump, SB size={self._sb.size}")
            return

        # Check if this is a low-barrier transition
        max_barrier = max(mechanism.barrier,
                          mechanism.barrier - mechanism.delta)

        if max_barrier < self.barrier_tol:
            # Case 2: Low barrier → expand superbasin
            if (self.dynamic_tol and
                    self._sb.size >= self.max_superbasin_size):
                # Overflow: shrink tolerance and reset
                prev_tol = self.barrier_tol
                self.barrier_tol = max(0.0, self.barrier_tol * self.tol_shrink)
                if self.debug:
                    print(f"SuperCache: Shrinking barrier_tol: "
                          f"{prev_tol:.3f} → {self.barrier_tol:.3f}")
                self._reset_with(positions, catalogue)
            else:
                new_basin = Basin(positions, catalogue, self.temperature,
                                  max_barrier=self.max_barrier,
                                  debug=self.debug)
                self._sb.expand_occupy(new_basin)
                self._sb.connect_from(source_basin, atom_index, mechanism)
                if self.debug:
                    print(f"SuperCache: Expanded SB to {self._sb.size}")
            return

        # Case 3: High barrier → cache current SB, try restore or create new
        cached_sb = self._try_restore(state_hash, positions)

        if cached_sb is not None:
            if self.debug:
                print(f"SuperCache: Restored cached SB with "
                      f"{cached_sb.size} basins, "
                      f"cache size={len(self._cache)}")
            self._cache_current()
            self._sb = cached_sb
            self._in_cache_count += 1
        else:
            if self.debug:
                print(f"SuperCache: New SB, cache size={len(self._cache)}")
            self._cache_current()
            basin = Basin(positions, catalogue, self.temperature,
                          max_barrier=self.max_barrier, debug=self.debug)
            self._sb = SuperBasin(basin, debug=self.debug)
            self._in_cache_count = 0

        # Dynamic tolerance: grow if cache is frequently reused
        if self.dynamic_tol and self._in_cache_count > self.cache_size:
            prev_tol = self.barrier_tol
            self.barrier_tol *= self.tol_grow
            if self.debug:
                print(f"SuperCache: Growing barrier_tol: "
                      f"{prev_tol:.3f} → {self.barrier_tol:.3f}")
            self._reset_with(positions, catalogue)

    def _cache_current(self):
        """Move current superbasin to cache."""
        if self._sb is not None:
            self._cache.append(self._sb)
            # Evict oldest if over capacity
            while len(self._cache) > self.cache_size:
                self._cache.pop(0)

    def _try_restore(self, state_hash, positions):
        """Try to find and restore a cached superbasin."""
        for i, sb in enumerate(self._cache):
            prev = sb.find_occupy(state_hash, positions, self.state_tol)
            if prev is not None:
                return self._cache.pop(i)
        return None

    def _reset_with(self, positions, catalogue):
        """Reset: create fresh SB from current state."""
        basin = Basin(positions, catalogue, self.temperature,
                      max_barrier=self.max_barrier, debug=self.debug)
        self._sb = SuperBasin(basin, debug=self.debug)

    def reset(self, positions=None, catalogue=None):
        """
        Full reset: clear cache and optionally re-initialize.
        """
        self._cache.clear()
        self._in_cache_count = 0
        if positions is not None and catalogue is not None:
            self._reset_with(positions, catalogue)
        else:
            self._sb = None

    def summary(self):
        """Print cache summary."""
        sb_size = self._sb.size if self._sb else 0
        n_cached = len(self._cache)
        print(f"SuperCache: active SB has {sb_size} basins, "
              f"{n_cached} cached SBs, "
              f"barrier_tol={self.barrier_tol:.3f} eV")
