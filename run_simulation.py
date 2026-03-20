#!/usr/bin/env python3
"""
Main simulation script for N2 reduction microkinetic calculations on Mo surface.

Runs both:
  1. Mean-field microkinetic model (ODE) - fast, for parameter scanning
  2. Lattice KMC simulation - stochastic, for spatial effects

Usage:
    python run_simulation.py                  # Run both models
    python run_simulation.py --mkm-only       # Mean-field only
    python run_simulation.py --kmc-only       # KMC only
    python run_simulation.py --scan-potential  # Scan applied potential
    python run_simulation.py --scan-temp      # Scan temperature
"""

import sys
import os
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from models.n2_reduction_Mo import (
    build_kmc_model, build_microkinetic_model, build_dft_microkinetic_model,
)
from mykmc.engine import KMCEngine
from mykmc.analysis import TrajectoryRecorder, run_to_steady_state
from mykmc.rates import tst_rate, electrochemical_rate
from mykmc.polarization import PolarizationCurve, compute_polarization_curve


def run_microkinetic(T=300.0, p_N2=1.0, U=-0.5):
    """Run mean-field microkinetic model and print results."""
    print("=" * 65)
    print("  Mean-Field Microkinetic Model: N2 Reduction on Mo")
    print("=" * 65)

    mkm = build_microkinetic_model()
    mkm.parameters['T'] = T
    mkm.parameters['p_N2'] = p_N2
    mkm.parameters['U'] = U

    print(f"\nConditions: T={T} K, p_N2={p_N2} bar, U={U} V vs RHE")

    # Solve steady state
    print("\nSolving steady state...")
    ss = mkm.solve_steady_state()
    mkm.print_summary(ss)

    # Time-dependent solution
    print("\nSolving time-dependent ODE...")
    t_span = (0, 1e3)  # 0 to 1000 seconds
    t_eval = np.logspace(-6, 3, 500)
    sol = mkm.solve_ode(t_span, t_eval=t_eval)

    if sol.success:
        print(f"ODE solved: {len(sol.t)} time points, "
              f"t_final = {sol.t[-1]:.2e} s")

        # Final coverages from ODE
        print("\nFinal coverages (ODE):")
        theta_final = sol.y[:, -1]
        theta_empty = 1.0 - np.sum(theta_final)
        print(f"  {'empty':<20s} {theta_empty:.6e}")
        for i, sp in enumerate(mkm.species):
            if theta_final[i] > 1e-12:
                print(f"  {sp:<20s} {theta_final[i]:.6e}")

        # TOF at final state
        tof = mkm.get_tof(
            {sp: theta_final[i] for i, sp in enumerate(mkm.species)})
        print(f"\nTOF at t={sol.t[-1]:.2e} s:")
        for obs, val in tof.items():
            print(f"  {obs}: {val:.6e} s^-1")
    else:
        print(f"ODE solver failed: {sol.message}")

    return mkm, ss


def run_kmc(T=300.0, p_N2=1.0, U=-0.5, lattice_size=20, n_steps=1000000):
    """Run lattice KMC simulation."""
    print("=" * 65)
    print("  Lattice KMC Simulation: N2 Reduction on Mo")
    print("=" * 65)

    pt = build_kmc_model()
    engine = KMCEngine(pt, size=[lattice_size, lattice_size],
                       print_rates=True, banner=True)

    # Set parameters
    engine.parameters.T = T
    engine.parameters.p_N2 = p_N2
    engine.parameters.U = U

    print(f"\nConditions: T={T} K, p_N2={p_N2} bar, U={U} V vs RHE")
    print(f"Lattice: {lattice_size}x{lattice_size} = "
          f"{lattice_size**2} sites")
    print(f"Running {n_steps} KMC steps...")

    # Run with trajectory recording
    recorder = TrajectoryRecorder(engine)
    steps_per_sample = max(1, n_steps // 200)

    for i in range(200):
        engine.do_steps(steps_per_sample)
        recorder.record()

        if (i + 1) % 20 == 0:
            cov = engine.get_coverage()
            non_empty = {k: v for k, v in cov.items()
                         if k != 'empty' and v > 1e-6}
            print(f"  [{100*(i+1)/200:5.1f}%] "
                  f"time={engine.kmc_time:.4e} s, "
                  f"coverages: {non_empty}")

    # Final results
    print(f"\nFinal state:")
    print(f"  KMC steps: {engine.kmc_step}")
    print(f"  KMC time:  {engine.kmc_time:.6e} s")
    engine.print_coverages()

    tof = engine.get_tof()
    if tof:
        print(f"\nTurn-over frequencies:")
        for obs, val in tof.items():
            print(f"  {obs}: {val:.6e} site^-1 s^-1")

    print(f"\nProcess execution counts:")
    stats = engine.get_process_stats()
    for name, count in sorted(stats.items(), key=lambda x: -x[1]):
        if count > 0:
            print(f"  {name}: {count}")

    return engine, recorder


def scan_potential(T=300.0, p_N2=1.0, U_range=None):
    """Scan applied potential and compute TOFs."""
    if U_range is None:
        U_range = np.linspace(-2.0, 0.0, 41)

    print("=" * 65)
    print("  Potential Scan: N2 Reduction on Mo (Mean-Field)")
    print("=" * 65)
    print(f"T = {T} K, p_N2 = {p_N2} bar")
    print(f"U range: {U_range[0]:.2f} to {U_range[-1]:.2f} V vs RHE\n")

    mkm = build_microkinetic_model()
    mkm.parameters['T'] = T
    mkm.parameters['p_N2'] = p_N2

    results = mkm.scan_parameter('U', U_range, observable='NH3_production')

    print(f"{'U [V]':>8s} {'TOF_NH3 [s^-1]':>15s} "
          f"{'theta_empty':>12s} {'theta_N2':>12s} "
          f"{'theta_N':>12s} {'theta_NH3':>12s}")
    print("-" * 72)

    for i, U in enumerate(results['values']):
        tof = results['tofs'][i]
        cov = results['coverages'][i]
        theta_empty = 1.0 - sum(cov.values())
        print(f"{U:>8.2f} {tof:>15.6e} "
              f"{theta_empty:>12.6e} {cov.get('N2', 0):>12.6e} "
              f"{cov.get('N', 0):>12.6e} {cov.get('NH3', 0):>12.6e}")

    # Save results
    outfile = 'potential_scan_results.dat'
    with open(outfile, 'w') as f:
        header = ['U_V', 'TOF_NH3'] + list(mkm.species) + ['empty']
        f.write('# ' + '\t'.join(header) + '\n')
        for i, U in enumerate(results['values']):
            tof = results['tofs'][i]
            cov = results['coverages'][i]
            theta_empty = 1.0 - sum(cov.values())
            vals = [f"{U:.4f}", f"{tof:.8e}"]
            vals += [f"{cov.get(sp, 0):.8e}" for sp in mkm.species]
            vals.append(f"{theta_empty:.8e}")
            f.write('\t'.join(vals) + '\n')
    print(f"\nResults saved to {outfile}")

    return results


def scan_temperature(p_N2=1.0, U=-0.5, T_range=None):
    """Scan temperature and compute TOFs."""
    if T_range is None:
        T_range = np.linspace(200, 500, 31)

    print("=" * 65)
    print("  Temperature Scan: N2 Reduction on Mo (Mean-Field)")
    print("=" * 65)
    print(f"p_N2 = {p_N2} bar, U = {U} V vs RHE")
    print(f"T range: {T_range[0]:.0f} to {T_range[-1]:.0f} K\n")

    mkm = build_microkinetic_model()
    mkm.parameters['p_N2'] = p_N2
    mkm.parameters['U'] = U

    results = mkm.scan_parameter('T', T_range, observable='NH3_production')

    print(f"{'T [K]':>8s} {'TOF_NH3 [s^-1]':>15s} "
          f"{'theta_empty':>12s} {'theta_N':>12s}")
    print("-" * 50)

    for i, T in enumerate(results['values']):
        tof = results['tofs'][i]
        cov = results['coverages'][i]
        theta_empty = 1.0 - sum(cov.values())
        print(f"{T:>8.1f} {tof:>15.6e} "
              f"{theta_empty:>12.6e} {cov.get('N', 0):>12.6e}")

    # Save results
    outfile = 'temperature_scan_results.dat'
    with open(outfile, 'w') as f:
        header = ['T_K', 'TOF_NH3'] + list(mkm.species) + ['empty']
        f.write('# ' + '\t'.join(header) + '\n')
        for i, T in enumerate(results['values']):
            tof = results['tofs'][i]
            cov = results['coverages'][i]
            theta_empty = 1.0 - sum(cov.values())
            vals = [f"{T:.2f}", f"{tof:.8e}"]
            vals += [f"{cov.get(sp, 0):.8e}" for sp in mkm.species]
            vals.append(f"{theta_empty:.8e}")
            f.write('\t'.join(vals) + '\n')
    print(f"\nResults saved to {outfile}")

    return results


def run_polarization(T=300.0, p_N2=1.0, U_range=None, dft_data=None):
    """
    Compute polarization curve (j-U) and TOF vs potential.

    If dft_data is provided, uses interpolated DFT barriers.
    Otherwise, uses Butler-Volmer model as demonstration.
    """
    print("=" * 65)
    print("  Polarization Curve: N2 Reduction on Mo")
    print("=" * 65)

    if dft_data is not None:
        # DFT data mode: potential-dependent barriers from constant-potential DFT
        print(f"\nUsing DFT energy data: {dft_data}")
        pc = build_dft_microkinetic_model(dft_data)
        pc.mkm.parameters['p_N2'] = p_N2
    else:
        # Butler-Volmer mode: use existing model as demo
        print("\nUsing Butler-Volmer model (no DFT data file provided)")
        mkm = build_microkinetic_model()
        mkm.parameters['p_N2'] = p_N2

        # N2 + 6H+ + 6e- -> 2NH3
        # NH3_production counts each NH3 molecule: 3 e- per count
        pc = PolarizationCurve(
            mkm=mkm,
            n_electrons={'NH3_production': 3},
            A_site=(3.2e-10) ** 2,
        )

    if U_range is None:
        U_range = np.linspace(-2.0, 0.0, 41)

    print(f"T = {T} K, p_N2 = {p_N2} bar")
    print(f"U range: {U_range[0]:.2f} to {U_range[-1]:.2f} V vs RHE\n")

    results = pc.compute(U_range, T=T, verbose=True)

    # Print summary table
    pc.print_results()

    # Save results
    outfile = 'polarization_curve.dat'
    pc.save(outfile)

    return results


def main():
    parser = argparse.ArgumentParser(
        description='N2 Reduction Microkinetic Simulation on Mo Surface')
    parser.add_argument('--mkm-only', action='store_true',
                        help='Run mean-field model only')
    parser.add_argument('--kmc-only', action='store_true',
                        help='Run KMC simulation only')
    parser.add_argument('--scan-potential', action='store_true',
                        help='Scan applied potential')
    parser.add_argument('--scan-temp', action='store_true',
                        help='Scan temperature')
    parser.add_argument('--polarization', action='store_true',
                        help='Compute polarization curve (j-U)')
    parser.add_argument('--dft-data', type=str, default=None,
                        help='Path to DFT energy data JSON file '
                             'for potential-dependent barriers')
    parser.add_argument('--T', type=float, default=300.0,
                        help='Temperature [K] (default: 300)')
    parser.add_argument('--p_N2', type=float, default=1.0,
                        help='N2 partial pressure [bar] (default: 1.0)')
    parser.add_argument('--U', type=float, default=-0.5,
                        help='Applied potential [V vs RHE] (default: -0.5)')
    parser.add_argument('--lattice-size', type=int, default=20,
                        help='KMC lattice size (default: 20)')
    parser.add_argument('--kmc-steps', type=int, default=1000000,
                        help='Number of KMC steps (default: 1000000)')

    args = parser.parse_args()

    if args.polarization:
        run_polarization(T=args.T, p_N2=args.p_N2, dft_data=args.dft_data)
    elif args.scan_potential:
        scan_potential(T=args.T, p_N2=args.p_N2)
    elif args.scan_temp:
        scan_temperature(p_N2=args.p_N2, U=args.U)
    elif args.mkm_only:
        run_microkinetic(T=args.T, p_N2=args.p_N2, U=args.U)
    elif args.kmc_only:
        run_kmc(T=args.T, p_N2=args.p_N2, U=args.U,
                lattice_size=args.lattice_size, n_steps=args.kmc_steps)
    else:
        # Run both
        run_microkinetic(T=args.T, p_N2=args.p_N2, U=args.U)
        print("\n\n")
        run_kmc(T=args.T, p_N2=args.p_N2, U=args.U,
                lattice_size=args.lattice_size, n_steps=args.kmc_steps)


if __name__ == '__main__':
    main()
