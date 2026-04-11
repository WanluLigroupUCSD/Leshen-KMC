"""
Environment catalogue — maps local environments to their mechanisms.

Ported from openFLY env/catalogue.hpp + catalogue.cpp.

Three-step matching pipeline:
  1. Graph hash → bucket lookup (O(1))
  2. Fingerprint equivalence → fast reject (O(n))
  3. Full permutation matching → Kabsch + greedy (O(n!), pruned)
"""

import json
import pickle
import numpy as np
from collections import defaultdict

from .mechanism import Mechanism
from .environment import Geometry, Fingerprint, build_geometry, build_neighbor_list


class CatalogueEntry:
    """
    A catalogued local environment with its mechanisms.

    Parameters
    ----------
    geometry : Geometry
        Reference geometry (canonical form).
    fingerprint : Fingerprint
        Distance fingerprint for fast matching.
    cat_index : int
        Unique index in the catalogue.
    delta_max : float
        Maximum RMSD tolerance for matching.
    """

    def __init__(self, geometry, fingerprint, cat_index, delta_max):
        self.geometry = geometry
        self.fingerprint = fingerprint
        self.cat_index = cat_index
        self.delta_max = delta_max
        self.mechanisms = []
        self.freq = 1
        self.false_pos = 0

    @property
    def n_mechs(self):
        return len(self.mechanisms)

    def __repr__(self):
        return (f"CatalogueEntry(idx={self.cat_index}, "
                f"n_atoms={self.geometry.n_atoms}, "
                f"n_mechs={self.n_mechs}, freq={self.freq})")


class _AtomEnv:
    """Per-atom runtime environment data (rebuilt each catalogue.rebuild())."""
    __slots__ = ('geo', 'fingerprint', 'graph_hash', 'entry')

    def __init__(self):
        self.geo = None
        self.fingerprint = None
        self.graph_hash = 0
        self.entry = None  # Reference to CatalogueEntry or None if new


class Catalogue:
    """
    Environment-to-mechanism catalogue with 3-step matching.

    Manages the mapping from local atomic environments (defined by
    geometry, element types, and connectivity) to known transition
    mechanisms. New environments trigger saddle point searches;
    known environments reuse cached mechanisms.

    Parameters
    ----------
    r_env : float
        Radius for local environment cutoff in Angstrom.
    r_edge : float
        Bonding distance for graph construction in Angstrom.
    delta_max : float
        Initial RMSD tolerance for environment matching.
    overfuzz : float
        Fingerprint fuzz factor (multiplied onto delta for fingerprint check).
    min_delta_max : float
        Minimum allowed delta_max after refinement.
    debug : bool
        Print debug information.
    """

    def __init__(self, r_env=5.0, r_edge=3.0, delta_max=0.15,
                 overfuzz=1.5, min_delta_max=0.01, debug=False):
        self.r_env = r_env
        self.r_edge = r_edge
        self.delta_max = delta_max
        self.overfuzz = overfuzz
        self.min_delta_max = min_delta_max
        self.debug = debug

        # Hash table: graph_hash → list[CatalogueEntry]
        self._buckets = defaultdict(list)
        self._size = 0

        # Per-atom runtime data (rebuilt on each rebuild())
        self._atoms = []

        # Neighbour list cache
        self._nl_cache = None

    @property
    def n_environments(self):
        """Total number of unique catalogued environments."""
        return self._size

    @property
    def n_mechanisms_total(self):
        """Total number of mechanisms across all environments."""
        return sum(e.n_mechs for bucket in self._buckets.values()
                   for e in bucket)

    def rebuild(self, positions, types, frozen, cell, pbc, r_nl=None):
        """
        Rebuild the catalogue for the current atomic configuration.

        Detects local environments around each atom, matches them against
        known catalogue entries, and returns indices of new (unknown)
        environments that need saddle point searches.

        Parameters
        ----------
        positions : ndarray, shape (n_atoms, 3)
        types : ndarray of int, shape (n_atoms,)
            Species type IDs (0, 1, 2, ...).
        frozen : ndarray of bool, shape (n_atoms,)
        cell : ndarray, shape (3, 3)
            Unit cell matrix.
        pbc : array-like of bool
            Periodic boundary conditions.
        r_nl : float, optional
            Neighbor list cutoff. Defaults to r_env + skin.

        Returns
        -------
        new_indices : list of int
            Atom indices with new (unknown) environments needing search.
        """
        n_atoms = len(positions)
        if r_nl is None:
            r_nl = self.r_env + 1.0  # skin

        # Build neighbor list
        neighbors = build_neighbor_list(positions, cell, pbc, r_nl)

        # Build per-atom environments
        self._atoms = []
        for i in range(n_atoms):
            ae = _AtomEnv()
            ae.geo = build_geometry(i, positions, types, frozen,
                                    neighbors[i], self.r_env)
            ae.fingerprint = ae.geo.fingerprint()
            ae.graph_hash = ae.geo.graph_hash(self.r_edge)
            ae.entry = None
            self._atoms.append(ae)

        # Ensure all hashes exist in bucket map
        for ae in self._atoms:
            if ae.graph_hash not in self._buckets:
                self._buckets[ae.graph_hash] = []

        # Find matches (can be parallelized in the future)
        for ae in self._atoms:
            ae.entry = self._find(ae)
            if ae.entry is not None:
                ae.entry.freq += 1

        # Collect new indices
        new_indices = []
        # Track newly inserted entries to avoid duplicates among new envs
        new_entries = []

        for i, ae in enumerate(self._atoms):
            if ae.entry is not None:
                continue

            # Check against other new environments first
            found_dup = False
            for ne_hash, ne_entry in new_entries:
                if ae.graph_hash == ne_hash and self._equiv(ae, ne_entry):
                    ae.entry = ne_entry
                    found_dup = True
                    break

            if not found_dup:
                entry = self._insert(ae)
                ae.entry = entry
                new_indices.append(i)
                new_entries.append((ae.graph_hash, entry))
                if self.debug:
                    print(f"CAT: New environment at atom {i}, "
                          f"n_atoms={ae.geo.n_atoms}")

        if self.debug and not new_indices:
            print("CAT: No new environments")

        return new_indices

    def get_entry(self, atom_index):
        """Get the CatalogueEntry for an atom."""
        return self._atoms[atom_index].entry

    def get_geometry(self, atom_index):
        """Get the local Geometry for an atom."""
        return self._atoms[atom_index].geo

    def get_mechanisms(self, atom_index):
        """Get all mechanisms for the environment around an atom."""
        entry = self._atoms[atom_index].entry
        if entry is None:
            return []
        return entry.mechanisms

    def set_mechanisms(self, atom_index, mechanisms):
        """
        Set the mechanisms for a newly discovered environment.

        Parameters
        ----------
        atom_index : int
            Atom index whose environment to set mechanisms for.
        mechanisms : list of Mechanism
        """
        entry = self._atoms[atom_index].entry
        if entry is None:
            raise ValueError(f"No catalogue entry for atom {atom_index}")
        if entry.mechanisms:
            raise ValueError(f"Entry already has {len(entry.mechanisms)} mechanisms")
        entry.mechanisms = list(mechanisms)

    def reconstruct(self, mechanism, atom_index, positions):
        """
        Reconstruct a mechanism on the current configuration.

        Applies the local displacements from the mechanism to the current
        positions, producing a new configuration (un-relaxed).

        Parameters
        ----------
        mechanism : Mechanism
            The mechanism to apply.
        atom_index : int
            Index of the central atom.
        positions : ndarray, shape (n_atoms, 3)
            Current positions.

        Returns
        -------
        new_positions : ndarray, shape (n_atoms, 3)
            Positions after applying mechanism displacements.
        O : ndarray, shape (3, 3)
            Rotation matrix used for alignment.
        """
        ae = self._atoms[atom_index]
        entry = ae.entry

        if entry is None:
            raise ValueError(f"No catalogue entry for atom {atom_index}")

        # Compute alignment: current geo → reference geo
        delta = self._calc_delta(ae.fingerprint, entry)
        result = ae.geo.permute_onto(entry.geometry, delta)

        if result is None:
            raise RuntimeError(
                f"Failed to align environment at atom {atom_index}")

        O = result['O'].T  # Transpose to go from ref frame to current frame

        # Apply displacements
        new_positions = positions.copy()
        geo = ae.geo

        for j in range(geo.n_atoms):
            global_idx = int(geo.indices[j])
            new_positions[global_idx] += O @ mechanism.delta_fwd[j]

        return new_positions, O

    def reconstruct_sp(self, mechanism, atom_index, positions):
        """
        Reconstruct the saddle point configuration.

        Same as reconstruct() but uses delta_sp instead of delta_fwd.
        """
        ae = self._atoms[atom_index]
        entry = ae.entry

        if entry is None:
            raise ValueError(f"No catalogue entry for atom {atom_index}")

        delta = self._calc_delta(ae.fingerprint, entry)
        result = ae.geo.permute_onto(entry.geometry, delta)

        if result is None:
            raise RuntimeError(
                f"Failed to align environment at atom {atom_index}")

        O = result['O'].T

        new_positions = positions.copy()
        geo = ae.geo

        for j in range(geo.n_atoms):
            global_idx = int(geo.indices[j])
            new_positions[global_idx] += O @ mechanism.delta_sp[j]

        return new_positions, O

    def refine_tolerance(self, atom_index, min_delta=None):
        """
        Tighten the matching tolerance for an atom's environment.

        Called when mechanism reconstruction fails — the tolerance was too
        loose, causing a false positive match.

        Parameters
        ----------
        atom_index : int
        min_delta : float, optional
            Minimum allowed delta_max.

        Returns
        -------
        new_delta : float
            The refined tolerance.
        """
        if min_delta is None:
            min_delta = self.min_delta_max

        ae = self._atoms[atom_index]
        entry = ae.entry

        delta = self._calc_delta(ae.fingerprint, entry)
        result = ae.geo.permute_onto(entry.geometry, delta)

        if result is None:
            raise RuntimeError("Cannot refine: permutation failed")

        new_delta = max(min_delta, result['rmsd'] / 1.5)

        if self.debug:
            print(f"CAT: Refining delta_max at atom {atom_index}: "
                  f"{entry.delta_max:.5f} → {new_delta:.5f}")

        if new_delta < self.min_delta_max:
            raise RuntimeError(
                f"delta_max too small: {entry.delta_max} → {new_delta}")

        entry.delta_max = new_delta
        entry.freq = 1
        entry.false_pos = 0

        return new_delta

    def optimize(self):
        """Re-order catalogue buckets by frequency (most common first)."""
        for h, bucket in self._buckets.items():
            if len(bucket) > 1:
                bucket.sort(
                    key=lambda e: e.freq / (e.false_pos + 1.0),
                    reverse=True,
                )

    # ---- Internal matching methods ----

    def _calc_delta(self, fingerprint, entry):
        """Compute matching tolerance from fingerprint r_min and entry delta_max."""
        return min(fingerprint.r_min() * 0.4, entry.delta_max)

    def _find(self, ae):
        """Find a matching CatalogueEntry for the atom environment, or None."""
        bucket = self._buckets.get(ae.graph_hash)
        if bucket is None:
            return None

        for entry in bucket:
            if self._equiv(ae, entry):
                return entry
        return None

    def _equiv(self, ae, entry):
        """Check if ae matches entry using fingerprint + geometry."""
        delta = self._calc_delta(ae.fingerprint, entry)

        # Step 2: Fingerprint check (fast reject)
        if not entry.fingerprint.equiv(ae.fingerprint, delta * self.overfuzz):
            return False

        # Step 3: Full permutation matching
        result = ae.geo.permute_onto(entry.geometry, delta)
        if result is not None:
            return True
        else:
            entry.false_pos += 1
            return False

    def _insert(self, ae):
        """Insert a new environment into the catalogue."""
        fp = ae.fingerprint
        delta_max = min(fp.r_min() * 0.4, self.delta_max)
        entry = CatalogueEntry(
            geometry=ae.geo,
            fingerprint=fp,
            cat_index=self._size,
            delta_max=delta_max,
        )
        self._size += 1
        self._buckets[ae.graph_hash].append(entry)
        return entry

    # ---- Serialization ----

    def save(self, filepath):
        """
        Save catalogue to file (pickle format).

        Saves all catalogue entries and their mechanisms for restart.
        """
        data = {
            'r_env': self.r_env,
            'r_edge': self.r_edge,
            'delta_max': self.delta_max,
            'overfuzz': self.overfuzz,
            'min_delta_max': self.min_delta_max,
            'entries': [],
        }

        for h, bucket in self._buckets.items():
            for entry in bucket:
                entry_data = {
                    'hash': h,
                    'cat_index': entry.cat_index,
                    'delta_max': entry.delta_max,
                    'freq': entry.freq,
                    'positions': entry.geometry.positions.tolist(),
                    'colours': entry.geometry.colours.tolist(),
                    'indices': entry.geometry.indices.tolist(),
                    'mechanisms': [m.to_dict() for m in entry.mechanisms],
                }
                data['entries'].append(entry_data)

        with open(filepath, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

        if self.debug:
            print(f"CAT: Saved {self._size} environments to {filepath}")

    @classmethod
    def load(cls, filepath, debug=False):
        """Load catalogue from file."""
        with open(filepath, 'rb') as f:
            data = pickle.load(f)

        cat = cls(
            r_env=data['r_env'],
            r_edge=data['r_edge'],
            delta_max=data['delta_max'],
            overfuzz=data.get('overfuzz', 1.5),
            min_delta_max=data.get('min_delta_max', 0.01),
            debug=debug,
        )

        for ed in data['entries']:
            geo = Geometry(
                np.array(ed['positions']),
                np.array(ed['colours']),
                np.array(ed['indices']),
            )
            fp = geo.fingerprint()
            entry = CatalogueEntry(
                geometry=geo,
                fingerprint=fp,
                cat_index=ed['cat_index'],
                delta_max=ed['delta_max'],
            )
            entry.freq = ed.get('freq', 1)
            entry.mechanisms = [Mechanism.from_dict(m)
                                for m in ed['mechanisms']]
            cat._buckets[ed['hash']].append(entry)
            cat._size = max(cat._size, entry.cat_index + 1)

        if debug:
            print(f"CAT: Loaded {cat._size} environments from {filepath}")

        return cat

    def summary(self):
        """Print catalogue summary."""
        print(f"Catalogue: {self.n_environments} unique environments, "
              f"{self.n_mechanisms_total} total mechanisms")
        for h, bucket in self._buckets.items():
            for entry in bucket:
                print(f"  Entry #{entry.cat_index}: "
                      f"n_atoms={entry.geometry.n_atoms}, "
                      f"n_mechs={entry.n_mechs}, "
                      f"freq={entry.freq}")

    def __repr__(self):
        return (f"Catalogue(n_envs={self.n_environments}, "
                f"n_mechs={self.n_mechanisms_total}, "
                f"r_env={self.r_env})")
