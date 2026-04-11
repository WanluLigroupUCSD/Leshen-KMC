"""
Basin — a local energy minimum with all accessible mechanisms and rates.

Ported from openFLY kinetic/basin.hpp + basin.cpp.
"""

import numpy as np
from hashlib import sha256


# Inverse Boltzmann constant in eV^-1, matching openFLY exactly
INV_BOLTZ = 16021766340.0 / 1380649.0


class LocalisedMech:
    """A mechanism localized to a specific atom with its rate."""
    __slots__ = ('atom_index', 'rate', 'mechanism', 'exit_mech')

    def __init__(self, atom_index, rate, mechanism):
        self.atom_index = atom_index
        self.rate = rate
        self.mechanism = mechanism
        self.exit_mech = True  # True = leads outside superbasin


class Basin:
    """
    A local energy minimum in the PES with all accessible mechanisms.

    Collects all mechanisms from the catalogue for the current configuration,
    computes their Arrhenius rates, and provides KMC event selection.

    Parameters
    ----------
    positions : ndarray, shape (n_atoms, 3)
        Atomic positions at this minimum.
    catalogue : Catalogue
        Environment catalogue with mechanisms.
    temperature : float
        Temperature in Kelvin.
    max_barrier : float
        Maximum barrier to include (eV). Mechanisms with higher barriers
        are excluded for efficiency.
    debug : bool
    """

    def __init__(self, positions, catalogue, temperature,
                 max_barrier=5.0, debug=False):
        self.positions = positions.copy()
        self.temperature = temperature
        self.debug = debug

        self.mechs = []
        self.rate_sum = 0.0
        self.connected = False  # True if any internal connection made

        # State hash: identifies unique configuration by catalogue indices
        n_atoms = len(positions)
        hash_indices = []
        for i in range(n_atoms):
            entry = catalogue.get_entry(i)
            if entry is not None:
                hash_indices.append(entry.cat_index)
            else:
                hash_indices.append(-1)

        h = sha256(np.array(hash_indices, dtype=np.int32).tobytes())
        self.state_hash = int(h.hexdigest()[:16], 16)

        # Collect all mechanisms with rates
        for i in range(n_atoms):
            for mech in catalogue.get_mechanisms(i):
                if mech.barrier < max_barrier and mech.barrier > 0:
                    rate = mech.kinetic_pre * np.exp(
                        -mech.barrier * INV_BOLTZ / temperature)
                    lm = LocalisedMech(i, rate, mech)
                    self.mechs.append(lm)
                    self.rate_sum += rate

        if debug:
            print(f"Basin: {len(self.mechs)} mechanisms, "
                  f"rate_sum={self.rate_sum:.4e}")

    def kmc_choice(self, rng=None):
        """
        Select a mechanism via KMC (n-fold way).

        Parameters
        ----------
        rng : numpy.random.Generator, optional
            Random number generator.

        Returns
        -------
        mechanism : Mechanism
            Selected mechanism.
        atom_index : int
            Central atom for the selected mechanism.
        dt : float
            KMC time step.
        """
        if rng is None:
            rng = np.random.default_rng()

        if self.rate_sum <= 0 or not self.mechs:
            raise RuntimeError("Basin has no mechanisms (rate_sum=0)")

        # Select mechanism proportional to rate
        r = rng.random()
        target = r * self.rate_sum
        cumul = 0.0
        selected = self.mechs[-1]

        for lm in self.mechs:
            cumul += lm.rate
            if cumul >= target:
                selected = lm
                break

        # Time step
        dt = -np.log(rng.random()) / self.rate_sum

        if self.debug:
            pct = selected.rate / self.rate_sum * 100
            print(f"Basin: KMC choice @atom={selected.atom_index}, "
                  f"barrier={selected.mechanism.barrier:.3f} eV, "
                  f"{pct:.1f}% of {len(self.mechs)} choices")

        return selected.mechanism, selected.atom_index, dt

    def most_likely(self, tol=0.01):
        """
        Get mechanisms with rate > tol * rate_sum.

        Returns list of (mechanism, atom_index, rate_fraction).
        """
        result = []
        for lm in self.mechs:
            frac = lm.rate / self.rate_sum if self.rate_sum > 0 else 0
            if frac > tol:
                result.append((lm.mechanism, lm.atom_index, frac))
        return result

    def state(self):
        """Return positions at this minimum."""
        return self.positions
