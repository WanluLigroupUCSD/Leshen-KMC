"""
SKMC Engine — on-the-fly off-lattice kinetic Monte Carlo simulation.

Ported from openFLY kinetic/skmc.hpp.

Main simulation loop:
  1. Relax initial structure
  2. Detect local environments
  3. For unknown environments: run saddle point searches
  4. Build basin with all mechanisms
  5. KMC step: select mechanism → reconstruct → relax → validate
  6. Update catalogue and superbasin connectivity
  7. Repeat from 5
"""

import numpy as np
import time as _time
import json
from pathlib import Path

from .potential import CalculatorAdapter
from .minimize import Minimizer
from .saddle import DimerSearch, SaddleMaster
from .catalogue import Catalogue
from .cache import SuperCache
from .basin import Basin


class SKMCEngine:
    """
    On-the-fly off-lattice kinetic Monte Carlo engine.

    Orchestrates the complete OLKMC simulation: environment detection,
    on-the-fly saddle point searches, mechanism catalogue management,
    and superbasin-accelerated KMC event selection.

    Parameters
    ----------
    atoms : ase.Atoms
        Initial atomic configuration with a calculator attached.
    temperature : float
        Simulation temperature in Kelvin.
    r_env : float
        Environment radius in Angstrom.
    r_edge : float
        Graph edge distance for environment hashing.
    delta_max : float
        RMSD tolerance for environment matching.
    frozen_indices : list of int, optional
        Indices of atoms to freeze during simulation.
    catalogue_file : str, optional
        Path to load/save catalogue for restart.
    barrier_tol : float
        SuperBasin barrier tolerance in eV.
    max_barrier : float
        Maximum mechanism barrier to consider in eV.
    recon_e_tol : float
        Energy tolerance for reconstruction validation (eV).
    recon_r_tol : float
        Position tolerance for reconstruction validation (Angstrom).
    debug : bool
    """

    def __init__(self, atoms, temperature=300.0,
                 r_env=5.0, r_edge=3.0, delta_max=0.15,
                 frozen_indices=None,
                 catalogue_file=None,
                 # Dimer options
                 dimer_f_tol=0.05, dimer_max_steps=300,
                 # Minimizer options
                 min_f_tol=0.01, min_max_steps=500,
                 # Master options
                 max_searches=50, max_failed=20,
                 r_pert=3.0, stddev=0.3,
                 # Basin options
                 barrier_tol=0.3, max_barrier=5.0,
                 # Reconstruction tolerances
                 recon_e_tol=0.1, recon_r_tol=0.5,
                 debug=False):

        self.atoms = atoms.copy()
        self.temperature = temperature
        self.debug = debug
        self.recon_e_tol = recon_e_tol
        self.recon_r_tol = recon_r_tol

        n = len(atoms)

        # Frozen mask
        self.frozen_mask = np.zeros(n, dtype=bool)
        if frozen_indices is not None:
            self.frozen_mask[frozen_indices] = True

        # Calculator adapter
        self.calc = CalculatorAdapter(atoms.calc, self.frozen_mask)

        # Minimizer
        self.minimizer = Minimizer(
            self.calc, f_tol=min_f_tol, max_steps=min_max_steps,
            debug=debug)

        # Dimer
        self.dimer = DimerSearch(
            self.calc, f_tol=dimer_f_tol, max_steps=dimer_max_steps,
            debug=debug)

        # Saddle master
        self.master = SaddleMaster(
            self.dimer, self.minimizer, self.calc,
            r_pert=r_pert, stddev=stddev,
            max_searches=max_searches, max_failed=max_failed,
            debug=debug)

        # Catalogue
        if catalogue_file and Path(catalogue_file).exists():
            self.catalogue = Catalogue.load(catalogue_file, debug=debug)
        else:
            self.catalogue = Catalogue(
                r_env=r_env, r_edge=r_edge, delta_max=delta_max,
                debug=debug)
        self._catalogue_file = catalogue_file

        # SuperCache
        self.super_cache = SuperCache(
            temperature=temperature,
            barrier_tol=barrier_tol,
            max_barrier=max_barrier,
            debug=debug)

        # Simulation state
        self.kmc_time = 0.0
        self.kmc_step = 0
        self.total_searches = 0
        self.total_mechanisms = 0

        # Random number generator
        self.rng = np.random.default_rng()

        # Trajectory
        self._trajectory = []

    def run(self, max_steps=1000, max_time=None, callback=None):
        """
        Run the SKMC simulation.

        Parameters
        ----------
        max_steps : int
            Maximum number of KMC steps.
        max_time : float, optional
            Maximum simulation time in seconds.
        callback : callable, optional
            Called after each step with signature:
                callback(step, time, energy_before, atom_index, mechanism,
                         energy_after, positions)
            Return True from callback to stop the simulation.

        Returns
        -------
        trajectory : list of dict
            Trajectory data for each step.
        """
        t_start = _time.time()

        # Step 1: Initial minimization
        if self.debug:
            print("SKMC: Minimizing initial structure...")
        self.atoms, e_init, conv, nit = self.minimizer.minimize(self.atoms)
        if not conv:
            print(f"WARNING: Initial minimization did not converge "
                  f"({nit} steps)")

        positions = self.atoms.get_positions()
        types = self._get_types()
        cell = np.array(self.atoms.get_cell())
        pbc = self.atoms.get_pbc()

        if self.debug:
            print(f"SKMC: Initial E={e_init:.6f} eV")

        # Step 2: Initial catalogue rebuild
        self._update_catalogue(positions, types, cell, pbc)

        # Step 3: Initialize SuperCache
        self.super_cache.initialize(positions, self.catalogue)

        # Main loop
        for step in range(max_steps):
            # Check time limit
            if max_time is not None:
                elapsed = _time.time() - t_start
                if elapsed > max_time:
                    if self.debug:
                        print(f"SKMC: Time limit reached ({elapsed:.1f}s)")
                    break

            # Step 4: KMC selection
            try:
                mech, atom_idx, dt, source_basin, switched = \
                    self.super_cache.kmc_choice(self.rng)
            except RuntimeError as e:
                print(f"SKMC: KMC selection failed: {e}")
                break

            e_before = self.calc.energy(self.atoms)

            # Step 5: Reconstruct mechanism
            try:
                new_positions, O = self.catalogue.reconstruct(
                    mech, atom_idx, positions)
            except RuntimeError as e:
                if self.debug:
                    print(f"SKMC: Reconstruction failed: {e}")
                # Refine tolerance and rebuild
                self.catalogue.refine_tolerance(atom_idx)
                self._update_catalogue(positions, types, cell, pbc)
                self.super_cache.reset(positions, self.catalogue)
                continue

            # Step 6: Relax reconstructed configuration
            atoms_trial = self.atoms.copy()
            atoms_trial.calc = self.calc.calculator
            atoms_trial.set_positions(new_positions)
            atoms_trial, e_after, conv, _ = self.minimizer.minimize(atoms_trial)

            if not conv:
                if self.debug:
                    print(f"SKMC: Post-reconstruction relaxation "
                          f"did not converge")
                # Refine and retry
                self.catalogue.refine_tolerance(atom_idx)
                self._update_catalogue(positions, types, cell, pbc)
                self.super_cache.reset(positions, self.catalogue)
                continue

            new_positions = atoms_trial.get_positions()

            # Step 7: Validate reconstruction
            dE = abs(e_after - e_before - mech.delta)
            dR = np.sqrt(np.mean(np.sum(
                (new_positions - positions)**2, axis=1)))

            if dE > self.recon_e_tol and dR > self.recon_r_tol:
                if self.debug:
                    print(f"SKMC: Reconstruction validation failed: "
                          f"dE={dE:.4f}, dR={dR:.4f}")
                self.catalogue.refine_tolerance(atom_idx)
                self._update_catalogue(positions, types, cell, pbc)
                self.super_cache.reset(positions, self.catalogue)
                continue

            # Step 8: Accept step
            self.kmc_time += dt
            self.kmc_step += 1
            self.atoms.set_positions(new_positions)
            positions = new_positions

            # Record trajectory
            step_data = {
                'step': self.kmc_step,
                'time': self.kmc_time,
                'dt': dt,
                'atom_index': atom_idx,
                'barrier': mech.barrier,
                'delta': mech.delta,
                'energy_before': e_before,
                'energy_after': e_after,
            }
            self._trajectory.append(step_data)

            if self.debug:
                print(f"SKMC: Step {self.kmc_step}, t={self.kmc_time:.6e}s, "
                      f"atom={atom_idx}, barrier={mech.barrier:.4f}eV, "
                      f"dE={mech.delta:.4f}eV")

            # Step 9: Update catalogue with new environments
            self._update_catalogue(positions, types, cell, pbc)

            # Step 10: Update superbasin connectivity
            self.super_cache.connect_from(
                source_basin, atom_idx, mech, positions, self.catalogue)

            # Callback
            if callback is not None:
                stop = callback(
                    self.kmc_step, self.kmc_time,
                    e_before, atom_idx, mech,
                    e_after, positions)
                if stop:
                    if self.debug:
                        print("SKMC: Stopped by callback")
                    break

        # Save catalogue
        if self._catalogue_file:
            self.catalogue.save(self._catalogue_file)

        return self._trajectory

    def _update_catalogue(self, positions, types, cell, pbc):
        """Rebuild catalogue and run searches for new environments."""
        new_indices = self.catalogue.rebuild(
            positions, types, self.frozen_mask, cell, pbc)

        if not new_indices:
            return

        # Run saddle point searches for each new environment
        for idx in new_indices:
            geo = self.catalogue.get_geometry(idx)
            if geo is None:
                continue

            if self.debug:
                print(f"SKMC: Searching mechanisms for atom {idx} "
                      f"(n_local={geo.n_atoms})")

            mechanisms = self.master.find_mechanisms(
                self.atoms, idx, geo, self.catalogue)

            self.catalogue.set_mechanisms(idx, mechanisms)
            self.total_searches += 1
            self.total_mechanisms += len(mechanisms)

    def _get_types(self):
        """Get integer type array from atoms."""
        numbers = self.atoms.get_atomic_numbers()
        unique = np.unique(numbers)
        type_map = {z: i for i, z in enumerate(unique)}
        return np.array([type_map[z] for z in numbers], dtype=int)

    # ---- Observables ----

    def get_trajectory(self):
        """Get simulation trajectory as list of step dictionaries."""
        return self._trajectory

    def get_positions(self):
        """Get current atomic positions."""
        return self.atoms.get_positions().copy()

    def get_displacements(self):
        """Get total displacement from initial positions."""
        if not self._trajectory:
            return np.zeros_like(self.atoms.get_positions())
        # Would need to store initial positions
        return self.atoms.get_positions()

    def summary(self):
        """Print simulation summary."""
        print(f"\nSKMC Simulation Summary")
        print(f"  Steps: {self.kmc_step}")
        print(f"  Time: {self.kmc_time:.6e} s")
        print(f"  Temperature: {self.temperature} K")
        print(f"  Atoms: {len(self.atoms)} "
              f"({np.sum(~self.frozen_mask)} free)")
        print(f"  Catalogue: {self.catalogue.n_environments} environments, "
              f"{self.catalogue.n_mechanisms_total} mechanisms")
        print(f"  Searches: {self.total_searches}")
        if self.super_cache.active_superbasin:
            print(f"  Active SB: {self.super_cache.active_superbasin.size} "
                  f"basins")
        self.super_cache.summary()

    def save_trajectory(self, filepath):
        """Save trajectory to JSON."""
        with open(filepath, 'w') as f:
            json.dump(self._trajectory, f, indent=2)
        print(f"Trajectory saved to {filepath} ({len(self._trajectory)} steps)")

    def __repr__(self):
        return (f"SKMCEngine(n_atoms={len(self.atoms)}, T={self.temperature}K, "
                f"step={self.kmc_step}, time={self.kmc_time:.6e}s)")
