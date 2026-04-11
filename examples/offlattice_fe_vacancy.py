"""
Off-Lattice KMC Example: Vacancy diffusion in BCC Fe.

Demonstrates the full SKMC pipeline:
  1. Build a BCC Fe supercell with a vacancy
  2. Attach an EMT calculator (fast, for demonstration)
  3. Run on-the-fly off-lattice KMC
  4. Track vacancy migration path and MSD

Usage:
    python examples/offlattice_fe_vacancy.py
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from spark.offlattice import SKMCEngine


def main():
    # ---- Step 1: Build BCC Fe supercell with vacancy ----
    try:
        from ase.build import bulk
        from ase.calculators.emt import EMT
    except ImportError:
        print("This example requires ASE. Install with: pip install ase")
        return

    # Create 4x4x4 BCC Fe supercell
    atoms = bulk('Fe', 'bcc', a=2.87, cubic=True) * (4, 4, 4)
    print(f"Created {len(atoms)} atom Fe supercell")

    # Create vacancy by removing center atom
    center = atoms.get_positions().mean(axis=0)
    dists = np.linalg.norm(atoms.get_positions() - center, axis=1)
    vac_idx = np.argmin(dists)
    print(f"Removing atom {vac_idx} to create vacancy")
    del atoms[vac_idx]
    print(f"System: {len(atoms)} atoms with 1 vacancy")

    # ---- Step 2: Attach calculator ----
    atoms.calc = EMT()

    # Freeze edge atoms (2 layers from each boundary)
    positions = atoms.get_positions()
    cell_diag = np.diag(atoms.get_cell())
    margin = 2.87 * 2  # 2 lattice constants
    frozen = []
    for i in range(len(atoms)):
        pos = positions[i]
        if (pos[0] < margin or pos[0] > cell_diag[0] - margin or
            pos[1] < margin or pos[1] > cell_diag[1] - margin or
            pos[2] < margin or pos[2] > cell_diag[2] - margin):
            frozen.append(i)
    print(f"Frozen {len(frozen)} edge atoms, "
          f"{len(atoms) - len(frozen)} free atoms")

    # ---- Step 3: Run SKMC ----
    print("\nStarting SKMC simulation...")

    engine = SKMCEngine(
        atoms,
        temperature=500.0,
        r_env=4.0,
        r_edge=3.2,
        frozen_indices=frozen,
        catalogue_file='fe_vacancy_catalogue.pkl',
        max_searches=30,
        max_failed=15,
        barrier_tol=0.2,
        debug=True,
    )

    # Track vacancy position
    vacancy_positions = []

    def callback(step, time, e_before, atom_idx, mech, e_after, positions):
        vacancy_positions.append({
            'step': step,
            'time': time,
            'atom': atom_idx,
            'barrier': mech.barrier,
            'delta': mech.delta,
        })

        if step % 10 == 0:
            print(f"\n--- Step {step}, t={time:.6e}s ---")
            print(f"  E: {e_before:.4f} → {e_after:.4f} eV")
            print(f"  Atom {atom_idx}: barrier={mech.barrier:.4f} eV")

        return step >= 100  # Stop after 100 steps

    trajectory = engine.run(
        max_steps=100,
        max_time=600,  # 10 minute wall-time limit
        callback=callback,
    )

    # ---- Step 4: Print results ----
    engine.summary()

    if trajectory:
        barriers = [t['barrier'] for t in trajectory]
        print(f"\nBarrier statistics:")
        print(f"  Mean: {np.mean(barriers):.4f} eV")
        print(f"  Min:  {np.min(barriers):.4f} eV")
        print(f"  Max:  {np.max(barriers):.4f} eV")

        times = [t['time'] for t in trajectory]
        if len(times) > 1:
            dt_avg = np.diff(times).mean()
            print(f"\nAverage time step: {dt_avg:.6e} s")
            print(f"Total simulated time: {times[-1]:.6e} s")

    engine.save_trajectory('fe_vacancy_trajectory.json')


if __name__ == '__main__':
    main()
