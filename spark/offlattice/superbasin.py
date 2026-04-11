"""
SuperBasin — group of basins connected by low barriers, accelerated via bac-MRM.

Ported from openFLY kinetic/superbasin.hpp + superbasin.cpp.

The mean-rate method (MRM) analytically treats transitions between
low-barrier basins, avoiding wasted KMC steps flickering between
degenerate minima. Only exit mechanisms (to outside the superbasin)
are selected.
"""

import numpy as np


class SuperBasin:
    """
    A group of interconnected basins (superbasin) with mean-rate-method
    acceleration.

    When basins are connected by low barriers, KMC wastes time flickering
    between them. The superbasin groups these basins and analytically
    computes the effective exit rate, selecting only exit mechanisms.

    Parameters
    ----------
    initial_basin : Basin
        The first basin in the superbasin.
    debug : bool
    """

    def __init__(self, initial_basin, debug=False):
        self.basins = [initial_basin]
        self.occupied = 0  # index of currently occupied basin
        self.debug = debug
        # Transition probability matrix (grows dynamically)
        self._prob = np.zeros((1, 1))

    @property
    def size(self):
        return len(self.basins)

    def expand_occupy(self, new_basin):
        """
        Add a new basin to the superbasin and occupy it.

        Parameters
        ----------
        new_basin : Basin
        """
        n = self.size
        self.basins.append(new_basin)

        # Expand probability matrix
        new_prob = np.zeros((n + 1, n + 1))
        new_prob[:n, :n] = self._prob
        self._prob = new_prob

        self.occupied = n  # Occupy the new basin

        if self.debug:
            print(f"SuperBasin: expanded to {self.size} basins")

    def find_occupy(self, state_hash, positions, tol):
        """
        Find a basin matching the given state and occupy it.

        Parameters
        ----------
        state_hash : int
            Hash of the configuration.
        positions : ndarray, shape (n_atoms, 3)
            Current positions.
        tol : float
            Position matching tolerance.

        Returns
        -------
        prev_occupied : int or None
            Index of previously occupied basin, or None if not found.
        """
        centroid_x = np.mean(positions, axis=0)

        for i, basin in enumerate(self.basins):
            if state_hash == basin.state_hash:
                # Compare positions (with centroid correction for drift)
                delta = centroid_x - np.mean(basin.positions, axis=0)
                diff = positions - basin.positions - delta
                d = np.sqrt(np.sum(diff**2))

                if self.debug:
                    print(f"SuperBasin: hash match at basin {i}, "
                          f"dr={d:.6f}, drift={np.linalg.norm(delta):.6f}")

                if d < tol:
                    prev = self.occupied
                    self.occupied = i
                    return prev

        return None

    def connect_from(self, basin_idx, atom_index, mechanism):
        """
        Mark a mechanism as internal (connecting two basins in the SB).

        Parameters
        ----------
        basin_idx : int
            Index of the source basin.
        atom_index : int
            Atom index of the mechanism.
        mechanism : Mechanism
            The mechanism that connected the basins.
        """
        # Find the matching mechanism in the source basin
        source = self.basins[basin_idx]
        for lm in source.mechs:
            if lm.atom_index == atom_index and lm.mechanism is mechanism:
                # Mark as internal (not exit)
                lm.exit_mech = False
                source.connected = True

                # Update transition probability matrix
                to_idx = self.occupied
                from_idx = basin_idx
                self._prob[to_idx, from_idx] = (
                    lm.rate / source.rate_sum)

                if self.debug:
                    print(f"SuperBasin: connected basin {from_idx} → "
                          f"{to_idx}")
                return

        raise RuntimeError(
            f"Mechanism not found in basin {basin_idx} "
            f"for atom {atom_index}")

    def kmc_choice(self, rng=None):
        """
        Select an exit mechanism from the superbasin.

        If the occupied basin has no internal connections, falls back
        to normal basin KMC selection. Otherwise uses the mean-rate
        method to select among exit mechanisms weighted by mean
        residence time.

        Parameters
        ----------
        rng : numpy.random.Generator, optional

        Returns
        -------
        mechanism : Mechanism
        atom_index : int
        dt : float
        source_basin : int
            Basin index from which the mechanism exits.
        switched : bool
            True if the selected mechanism is from a different basin
            than the currently occupied one.
        """
        if rng is None:
            rng = np.random.default_rng()

        # If no connections, use normal KMC
        if not self.basins[self.occupied].connected:
            mech, atom, dt = self.basins[self.occupied].kmc_choice(rng)
            return mech, atom, dt, self.occupied, False

        # Mean-rate method
        tau = self._compute_tau()
        tau_sum = np.sum(tau)
        if tau_sum <= 0:
            # Fallback
            mech, atom, dt = self.basins[self.occupied].kmc_choice(rng)
            return mech, atom, dt, self.occupied, False

        tau_norm = tau / tau_sum

        # Compute total exit rate weighted by residence probability
        exit_rates = []
        for i, basin in enumerate(self.basins):
            for lm in basin.mechs:
                if lm.exit_mech:
                    eff_rate = tau_norm[i] * lm.rate
                    exit_rates.append((i, lm, eff_rate))

        if not exit_rates:
            # No exit mechanisms — all internal
            mech, atom, dt = self.basins[self.occupied].kmc_choice(rng)
            return mech, atom, dt, self.occupied, False

        r_sum = sum(er for _, _, er in exit_rates)

        # Select exit mechanism
        target = rng.random() * r_sum
        cumul = 0.0
        selected_basin = 0
        selected_lm = exit_rates[0][1]

        for basin_idx, lm, eff_rate in exit_rates:
            cumul += eff_rate
            if cumul >= target:
                selected_basin = basin_idx
                selected_lm = lm
                break

        dt = -np.log(rng.random()) / r_sum
        switched = (self.occupied != selected_basin)

        if self.debug:
            pct = tau_norm[selected_basin] * selected_lm.rate / r_sum * 100
            print(f"SuperBasin: SKMC choice @atom={selected_lm.atom_index}, "
                  f"barrier={selected_lm.mechanism.barrier:.3f} eV, "
                  f"{pct:.1f}% of {len(exit_rates)} exit mechs")

        return (selected_lm.mechanism, selected_lm.atom_index,
                dt, selected_basin, switched)

    def _compute_tau(self):
        """
        Compute mean residence times in each basin.

        Solves: tau = (I - P)^{-1} @ theta
        where theta is the Kronecker delta at the occupied basin.

        Returns
        -------
        tau : ndarray, shape (n_basins,)
            Mean residence time ratios.
        """
        n = self.size
        theta = np.zeros(n)
        theta[self.occupied] = 1.0

        I = np.eye(n)
        A = I - self._prob

        try:
            tau = np.linalg.solve(A, theta)
        except np.linalg.LinAlgError:
            tau = theta.copy()

        # Convert to actual times: tau_i / rate_sum_i
        for i in range(n):
            if self.basins[i].rate_sum > 0:
                tau[i] /= self.basins[i].rate_sum

        return tau
