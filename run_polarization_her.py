#!/usr/bin/env python3
"""
T4: HER on Pt(111) — Polarization Curves & Experimental Comparison

Generates:
  1. MKM polarization curve (j-V) over U = [0, -0.5] V vs RHE
  2. KMC polarization at selected potentials for comparison
  3. Tafel plot (log|j| vs η)
  4. Coverage vs potential
  5. Comparison with experimental exchange current density j₀ ≈ 1-10 mA/cm²

Output:
  results/her_mkm_polarization.dat   — MKM data
  results/her_kmc_polarization.dat   — KMC data
  results/fig_polarization.png       — j-V curve
  results/fig_tafel.png              — Tafel plot
  results/fig_coverage.png           — θ_H vs U
  results/fig_mkm_vs_kmc.png         — MKM vs KMC comparison
  results/fig_summary.png            — 4-panel summary figure
"""

import sys
import os
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from spark.microkinetic import MicroKineticModel
from spark.rates import tst_rate, electrochemical_rate
from spark.units import kB, h, eV, e_charge
from spark.polarization import tof_to_current_density

# ============================================================================
#  DFT Parameters
# ============================================================================
EA_VOLMER_FWD = 0.67
EA_VOLMER_REV = 0.62
EA_TAFEL = 0.85
EA_HEYROVSKY = 0.70
EA_DIFFUSION = 0.10
EPS_HH = 0.10
ALPHA_BEP = 0.5
BETA_BV = 0.5
U0 = 0.0
A_PT = 2.775  # Angstrom
A_SITE = (A_PT * 1e-10) ** 2  # m^2
T_DEFAULT = 298.0

# Experimental reference: Pt(111) HER exchange current density
# Markovic et al., Surface Science 499, 149 (2002)
# Sheng et al., Nature Energy 2, 17070 (2017)
J0_EXP_RANGE = (0.5, 10.0)  # mA/cm^2


# ============================================================================
#  MKM Model
# ============================================================================

def build_mkm():
    """Build the mean-field MKM for HER on Pt(111)."""
    mkm = MicroKineticModel()
    mkm.add_species('H')
    mkm.parameters = {'T': T_DEFAULT, 'U': 0.0, 'beta_BV': BETA_BV}

    def k_volmer_fwd(p):
        return electrochemical_rate(EA_VOLMER_FWD, p['T'], p['U'],
                                    U0=0.0, beta_bv=p.get('beta_BV', 0.5))

    def k_volmer_rev(p):
        T, U = p['T'], p['U']
        beta = p.get('beta_BV', 0.5)
        Ea_eff = EA_VOLMER_REV - (1.0 - beta) * (U - U0)
        Ea_eff = max(Ea_eff, 0.0)
        return (kB * T / h) * np.exp(-Ea_eff * eV / (kB * T))

    def k_tafel(p):
        return tst_rate(EA_TAFEL, p['T'])

    def k_heyrovsky_fwd(p):
        return electrochemical_rate(EA_HEYROVSKY, p['T'], p['U'],
                                    U0=0.0, beta_bv=p.get('beta_BV', 0.5))

    mkm.add_reaction(name='Volmer', reactants={'empty': 1}, products={'H': 1},
                     rate_fwd=k_volmer_fwd, rate_rev=k_volmer_rev,
                     tof_count={'Volmer_net': 1})
    mkm.add_reaction(name='Tafel', reactants={'H': 2}, products={'empty': 2},
                     rate_fwd=k_tafel, tof_count={'H2_Tafel': 1})
    mkm.add_reaction(name='Heyrovsky', reactants={'H': 1}, products={'empty': 1},
                     rate_fwd=k_heyrovsky_fwd, tof_count={'H2_Heyrovsky': 1})
    return mkm


# ============================================================================
#  MKM Polarization Curve
# ============================================================================

def compute_mkm_polarization(U_range, T=T_DEFAULT, verbose=True):
    """Compute MKM polarization curve."""
    mkm = build_mkm()
    mkm.parameters['T'] = T

    results = {
        'U': np.array(U_range),
        'theta_H': np.zeros(len(U_range)),
        'tof_Tafel': np.zeros(len(U_range)),
        'tof_Heyrovsky': np.zeros(len(U_range)),
        'tof_H2_total': np.zeros(len(U_range)),
        'j_Tafel': np.zeros(len(U_range)),
        'j_Heyrovsky': np.zeros(len(U_range)),
        'j_total': np.zeros(len(U_range)),
    }

    theta_guess = None
    for i, U in enumerate(U_range):
        mkm.parameters['U'] = U
        if theta_guess is not None:
            ss = mkm.solve_steady_state(theta0=theta_guess)
        else:
            ss = mkm.solve_steady_state()
        theta_guess = [ss.get(sp, 0.0) for sp in mkm.species]

        tof = mkm.get_tof(ss)
        theta_H = ss.get('H', 0.0)
        tof_T = tof.get('H2_Tafel', 0.0)
        tof_Hey = tof.get('H2_Heyrovsky', 0.0)

        results['theta_H'][i] = theta_H
        results['tof_Tafel'][i] = tof_T
        results['tof_Heyrovsky'][i] = tof_Hey
        results['tof_H2_total'][i] = tof_T + tof_Hey
        results['j_Tafel'][i] = tof_to_current_density(tof_T, 2, A_SITE)
        results['j_Heyrovsky'][i] = tof_to_current_density(tof_Hey, 2, A_SITE)
        results['j_total'][i] = results['j_Tafel'][i] + results['j_Heyrovsky'][i]

        if verbose and (i % 10 == 0 or i == len(U_range) - 1):
            print(f"  U = {U:+.3f} V | θ_H = {theta_H:.4f} | "
                  f"j = {results['j_total'][i]:.4e} mA/cm²")

    return results


# ============================================================================
#  KMC Polarization Curve
# ============================================================================

def compute_kmc_polarization(U_range, T=T_DEFAULT, lattice_size=20,
                              equil_steps=50000, prod_steps=50000,
                              with_lateral=True, verbose=True):
    """Compute KMC polarization curve (without diffusion for speed)."""
    from spark.types import Project, Site, Condition, Action, Coord
    from spark.engine import KMCEngine

    U_range = np.asarray(U_range)
    results = {
        'U': U_range,
        'theta_H': np.zeros(len(U_range)),
        'tof_Tafel': np.zeros(len(U_range)),
        'tof_Heyrovsky': np.zeros(len(U_range)),
        'j_Tafel': np.zeros(len(U_range)),
        'j_Heyrovsky': np.zeros(len(U_range)),
        'j_total': np.zeros(len(U_range)),
    }

    for i, U in enumerate(U_range):
        pt = Project()
        pt.set_meta(model_name='HER_Pt111', model_dimension=2)
        pt.add_species(name='empty', color='#ffffff')
        pt.add_species(name='H', color='#ff4444')

        layer = pt.add_layer(name='Pt111')
        layer.sites.append(
            Site(name='top', pos=(0.5, 0.5, 0.5), default_species='empty'))
        pt.lattice.cell = np.diag([A_PT, A_PT, 15.0])

        pt.add_parameter(name='T', value=T)
        pt.add_parameter(name='U', value=U)
        pt.add_parameter(name='beta_BV', value=BETA_BV)
        pt.add_parameter(name='Ea_Volmer_fwd', value=EA_VOLMER_FWD)
        pt.add_parameter(name='Ea_Volmer_rev', value=EA_VOLMER_REV)
        pt.add_parameter(name='Ea_Tafel', value=EA_TAFEL)
        pt.add_parameter(name='Ea_Heyrovsky', value=EA_HEYROVSKY)

        coord = pt.lattice.generate_coord('top')
        coord_nn = Coord(offset=(1, 0, 0), site='top', layer='Pt111')

        pt.add_process(
            name='Volmer_fwd',
            conditions=[Condition(coord, 'empty')],
            actions=[Action(coord, 'H')],
            rate_constant='kB*T/h*exp(-(Ea_Volmer_fwd + beta_BV*U)*eV/(kB*T))',
            tof_count={'Volmer_fwd': 1},
        )
        pt.add_process(
            name='Volmer_rev',
            conditions=[Condition(coord, 'H')],
            actions=[Action(coord, 'empty')],
            rate_constant='kB*T/h*exp(-(Ea_Volmer_rev - (1-beta_BV)*U)*eV/(kB*T))',
            tof_count={'Volmer_rev': 1},
        )
        pt.add_process(
            name='Tafel',
            conditions=[Condition(coord, 'H'), Condition(coord_nn, 'H')],
            actions=[Action(coord, 'empty'), Action(coord_nn, 'empty')],
            rate_constant='kB*T/h*exp(-Ea_Tafel*eV/(kB*T))',
            tof_count={'H2_Tafel': 1},
        )
        pt.add_process(
            name='Heyrovsky_fwd',
            conditions=[Condition(coord, 'H')],
            actions=[Action(coord, 'empty')],
            rate_constant='kB*T/h*exp(-(Ea_Heyrovsky + beta_BV*U)*eV/(kB*T))',
            tof_count={'H2_Heyrovsky': 1},
        )

        if with_lateral:
            pt.add_lateral_interaction('H', 'H', energy=EPS_HH)
            pt.add_bep_relation('Volmer_fwd', alpha=ALPHA_BEP)
            pt.add_bep_relation('Volmer_rev', alpha=1.0 - ALPHA_BEP)

        engine = KMCEngine(pt, size=[lattice_size, lattice_size],
                           print_rates=False, banner=False)

        # Equilibrate
        engine.do_steps(equil_steps)
        # Reset TOF counters by calling get_tof (updates internal baseline)
        engine.get_tof()

        # Production
        engine.do_steps(prod_steps)
        cov = engine.get_coverage()
        tof = engine.get_tof()

        theta_H = cov.get('H', 0.0)
        tof_T = tof.get('H2_Tafel', 0.0)
        tof_Hey = tof.get('H2_Heyrovsky', 0.0)

        results['theta_H'][i] = theta_H
        results['tof_Tafel'][i] = tof_T
        results['tof_Heyrovsky'][i] = tof_Hey
        results['j_Tafel'][i] = tof_to_current_density(tof_T, 2, A_SITE)
        results['j_Heyrovsky'][i] = tof_to_current_density(tof_Hey, 2, A_SITE)
        results['j_total'][i] = results['j_Tafel'][i] + results['j_Heyrovsky'][i]

        if verbose:
            print(f"  U = {U:+.3f} V | θ_H = {theta_H:.4f} | "
                  f"j = {results['j_total'][i]:.4e} mA/cm² "
                  f"(Tafel: {results['j_Tafel'][i]:.4e}, "
                  f"Hey: {results['j_Heyrovsky'][i]:.4e})")

    return results


# ============================================================================
#  Analysis: Tafel Slope
# ============================================================================

def compute_tafel_slope(U, j, U_range=(-0.15, -0.02)):
    """
    Compute Tafel slope from polarization data.

    Fits log10|j| vs η (η = -U for HER at RHE) in the specified range.
    Returns slope in mV/decade.
    """
    eta = -np.array(U)  # overpotential (positive)
    log_j = np.log10(np.abs(np.array(j)) + 1e-30)

    mask = (np.array(U) >= U_range[0]) & (np.array(U) <= U_range[1])
    if np.sum(mask) < 3:
        return None, None, None

    eta_fit = eta[mask]
    log_j_fit = log_j[mask]

    # Remove points with very low current
    valid = log_j_fit > -20
    if np.sum(valid) < 3:
        return None, None, None

    eta_fit = eta_fit[valid]
    log_j_fit = log_j_fit[valid]

    # Linear fit: η = a + b * log10|j|  =>  b = Tafel slope
    # Or equivalently: log10|j| = c + d * η  =>  Tafel slope = 1/d
    coeffs = np.polyfit(log_j_fit, eta_fit * 1000, 1)  # mV
    tafel_slope = coeffs[0]  # mV/decade

    return tafel_slope, eta_fit, log_j_fit


# ============================================================================
#  Data Saving
# ============================================================================

def save_data(results, filename, method='MKM'):
    """Save polarization data to tab-separated file."""
    with open(filename, 'w') as f:
        f.write(f"# HER on Pt(111) - {method} Polarization Curve\n")
        f.write(f"# T = {T_DEFAULT} K, A_site = {A_SITE:.4e} m^2\n")
        f.write("# U_V\ttheta_H\tTOF_Tafel\tTOF_Heyrovsky\t"
                "j_Tafel_mA_cm2\tj_Heyrovsky_mA_cm2\tj_total_mA_cm2\n")
        for i in range(len(results['U'])):
            f.write(f"{results['U'][i]:.4f}\t"
                    f"{results['theta_H'][i]:.6e}\t"
                    f"{results['tof_Tafel'][i]:.6e}\t"
                    f"{results['tof_Heyrovsky'][i]:.6e}\t"
                    f"{results['j_Tafel'][i]:.6e}\t"
                    f"{results['j_Heyrovsky'][i]:.6e}\t"
                    f"{results['j_total'][i]:.6e}\n")
    print(f"  Data saved to {filename}")


# ============================================================================
#  Plotting
# ============================================================================

def make_plots(mkm_results, kmc_results=None, outdir='results'):
    """Generate publication-quality figures."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 14,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 10,
        'figure.dpi': 150,
        'savefig.dpi': 200,
        'savefig.bbox': 'tight',
    })

    U_mkm = mkm_results['U']
    j_mkm = mkm_results['j_total']
    theta_mkm = mkm_results['theta_H']

    # ---- Figure 1: Polarization curve (j vs U) ----
    fig1, ax1 = plt.subplots(figsize=(7, 5))
    ax1.plot(U_mkm, j_mkm, 'b-', lw=2, label='MKM (mean-field)')
    if kmc_results is not None:
        ax1.plot(kmc_results['U'], kmc_results['j_total'],
                 'ro', ms=8, label='KMC (spatial)')
    # Experimental j₀ band
    ax1.axhspan(J0_EXP_RANGE[0], J0_EXP_RANGE[1], alpha=0.15, color='green',
                label=f'Exp. j₀ = {J0_EXP_RANGE[0]}-{J0_EXP_RANGE[1]} mA/cm²')
    ax1.set_xlabel('U [V vs RHE]')
    ax1.set_ylabel('j [mA/cm²]')
    ax1.set_title('HER on Pt(111) — Polarization Curve')
    ax1.legend()
    ax1.set_xlim([U_mkm[-1] - 0.02, U_mkm[0] + 0.02])
    ax1.grid(True, alpha=0.3)
    fig1.savefig(os.path.join(outdir, 'fig_polarization.png'))
    plt.close(fig1)
    print(f"  Saved {outdir}/fig_polarization.png")

    # ---- Figure 2: Tafel plot (log|j| vs η) ----
    fig2, ax2 = plt.subplots(figsize=(7, 5))
    eta_mkm = -U_mkm * 1000  # mV
    log_j_mkm = np.log10(np.abs(j_mkm) + 1e-30)

    ax2.plot(log_j_mkm, eta_mkm, 'b-', lw=2, label='MKM')
    if kmc_results is not None:
        eta_kmc = -kmc_results['U'] * 1000
        log_j_kmc = np.log10(np.abs(kmc_results['j_total']) + 1e-30)
        ax2.plot(log_j_kmc, eta_kmc, 'ro', ms=8, label='KMC')

    # Compute and annotate Tafel slope
    tafel_slope, _, _ = compute_tafel_slope(U_mkm, j_mkm)
    if tafel_slope is not None:
        ax2.annotate(f'Tafel slope ≈ {tafel_slope:.0f} mV/dec',
                     xy=(0.05, 0.95), xycoords='axes fraction',
                     fontsize=12, va='top',
                     bbox=dict(boxstyle='round', fc='wheat', alpha=0.8))

    # Reference lines for 30 and 120 mV/dec
    log_j_ref = np.linspace(-2, 4, 50)
    ax2.plot(log_j_ref, 30 * log_j_ref + 30, 'g--', alpha=0.4,
             label='30 mV/dec (Tafel-limited)')
    ax2.plot(log_j_ref, 120 * log_j_ref - 100, 'm--', alpha=0.4,
             label='120 mV/dec (Volmer-limited)')

    ax2.set_xlabel('log₁₀|j| [mA/cm²]')
    ax2.set_ylabel('η [mV]')
    ax2.set_title('HER on Pt(111) — Tafel Plot')
    ax2.set_ylim([0, 500])
    ax2.set_xlim([-2, 4])
    ax2.legend(loc='lower right')
    ax2.grid(True, alpha=0.3)
    fig2.savefig(os.path.join(outdir, 'fig_tafel.png'))
    plt.close(fig2)
    print(f"  Saved {outdir}/fig_tafel.png")

    # ---- Figure 3: Coverage vs potential ----
    fig3, ax3 = plt.subplots(figsize=(7, 5))
    ax3.plot(U_mkm, theta_mkm, 'b-', lw=2, label='θ_H (MKM)')
    if kmc_results is not None:
        ax3.plot(kmc_results['U'], kmc_results['theta_H'],
                 'ro', ms=8, label='θ_H (KMC)')
    ax3.set_xlabel('U [V vs RHE]')
    ax3.set_ylabel('θ_H')
    ax3.set_title('HER on Pt(111) — H Coverage vs Potential')
    ax3.set_ylim([-0.02, 1.02])
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    fig3.savefig(os.path.join(outdir, 'fig_coverage.png'))
    plt.close(fig3)
    print(f"  Saved {outdir}/fig_coverage.png")

    # ---- Figure 4: MKM vs KMC comparison ----
    if kmc_results is not None:
        fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(12, 5))

        # Coverage comparison
        ax4a.plot(U_mkm, theta_mkm, 'b-', lw=2, label='MKM')
        ax4a.plot(kmc_results['U'], kmc_results['theta_H'],
                  'ro', ms=8, label='KMC')
        ax4a.set_xlabel('U [V vs RHE]')
        ax4a.set_ylabel('θ_H')
        ax4a.set_title('Coverage Comparison')
        ax4a.legend()
        ax4a.grid(True, alpha=0.3)

        # Current comparison
        ax4b.semilogy(U_mkm, np.abs(j_mkm) + 1e-30, 'b-', lw=2, label='MKM')
        ax4b.semilogy(kmc_results['U'],
                      np.abs(kmc_results['j_total']) + 1e-30,
                      'ro', ms=8, label='KMC')
        ax4b.set_xlabel('U [V vs RHE]')
        ax4b.set_ylabel('|j| [mA/cm²]')
        ax4b.set_title('Current Density Comparison')
        ax4b.legend()
        ax4b.grid(True, alpha=0.3)

        fig4.tight_layout()
        fig4.savefig(os.path.join(outdir, 'fig_mkm_vs_kmc.png'))
        plt.close(fig4)
        print(f"  Saved {outdir}/fig_mkm_vs_kmc.png")

    # ---- Figure 5: Summary (4-panel) ----
    fig5, axes = plt.subplots(2, 2, figsize=(12, 10))

    # (a) j vs U
    ax = axes[0, 0]
    ax.plot(U_mkm, j_mkm, 'b-', lw=2, label='MKM total')
    ax.plot(U_mkm, mkm_results['j_Tafel'], 'g--', lw=1.5, label='Tafel')
    ax.plot(U_mkm, mkm_results['j_Heyrovsky'], 'r--', lw=1.5, label='Heyrovsky')
    if kmc_results is not None:
        ax.plot(kmc_results['U'], kmc_results['j_total'], 'ko', ms=6,
                label='KMC total')
    ax.axhspan(J0_EXP_RANGE[0], J0_EXP_RANGE[1], alpha=0.12, color='green')
    ax.set_xlabel('U [V vs RHE]')
    ax.set_ylabel('j [mA/cm²]')
    ax.set_title('(a) Polarization Curve')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (b) Tafel plot
    ax = axes[0, 1]
    eta_mkm_v = -U_mkm  # in V
    log_j_abs = np.log10(np.abs(j_mkm) + 1e-30)
    ax.plot(log_j_abs, eta_mkm_v * 1000, 'b-', lw=2, label='MKM')
    if kmc_results is not None:
        log_j_kmc = np.log10(np.abs(kmc_results['j_total']) + 1e-30)
        ax.plot(log_j_kmc, -kmc_results['U'] * 1000, 'ro', ms=6, label='KMC')
    if tafel_slope is not None:
        ax.annotate(f'Tafel slope ≈ {tafel_slope:.0f} mV/dec',
                    xy=(0.05, 0.95), xycoords='axes fraction',
                    fontsize=11, va='top',
                    bbox=dict(boxstyle='round', fc='wheat', alpha=0.8))
    ax.set_xlabel('log₁₀|j| [mA/cm²]')
    ax.set_ylabel('η [mV]')
    ax.set_title('(b) Tafel Plot')
    ax.set_ylim([0, 500])
    ax.set_xlim([-2, 4])
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (c) Coverage
    ax = axes[1, 0]
    ax.plot(U_mkm, theta_mkm, 'b-', lw=2, label='MKM')
    if kmc_results is not None:
        ax.plot(kmc_results['U'], kmc_results['theta_H'], 'ro', ms=6,
                label='KMC')
    ax.set_xlabel('U [V vs RHE]')
    ax.set_ylabel('θ_H')
    ax.set_title('(c) H Coverage')
    ax.set_ylim([-0.02, 1.02])
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (d) Pathway selectivity
    ax = axes[1, 1]
    total = mkm_results['j_Tafel'] + mkm_results['j_Heyrovsky'] + 1e-30
    frac_tafel = mkm_results['j_Tafel'] / total
    frac_hey = mkm_results['j_Heyrovsky'] / total
    ax.fill_between(U_mkm, 0, frac_tafel, alpha=0.5, color='green',
                    label='Tafel pathway')
    ax.fill_between(U_mkm, frac_tafel, 1, alpha=0.5, color='red',
                    label='Heyrovsky pathway')
    ax.set_xlabel('U [V vs RHE]')
    ax.set_ylabel('Fractional contribution')
    ax.set_title('(d) Pathway Selectivity')
    ax.set_ylim([0, 1])
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig5.suptitle('HER on Pt(111) — SPARK Validation (T4)', fontsize=15)
    fig5.tight_layout(rect=[0, 0, 1, 0.96])
    fig5.savefig(os.path.join(outdir, 'fig_summary.png'))
    plt.close(fig5)
    print(f"  Saved {outdir}/fig_summary.png")


# ============================================================================
#  Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='T4: HER Pt(111) polarization curves')
    parser.add_argument('--no-kmc', action='store_true',
                        help='Skip KMC (only MKM)')
    parser.add_argument('--kmc-size', type=int, default=20,
                        help='KMC lattice size (default: 20)')
    parser.add_argument('--kmc-equil', type=int, default=50000,
                        help='KMC equilibration steps (default: 50000)')
    parser.add_argument('--kmc-prod', type=int, default=50000,
                        help='KMC production steps (default: 50000)')
    parser.add_argument('--kmc-points', type=int, default=11,
                        help='Number of KMC potential points (default: 11)')
    parser.add_argument('--no-lateral', action='store_true',
                        help='Disable lateral interactions in KMC')
    parser.add_argument('--T', type=float, default=T_DEFAULT,
                        help='Temperature [K]')
    parser.add_argument('--outdir', default='results',
                        help='Output directory')
    args = parser.parse_args()

    outdir = os.path.join(os.path.dirname(__file__), args.outdir)
    os.makedirs(outdir, exist_ok=True)

    print("=" * 65)
    print("  T4: HER on Pt(111) — Polarization Curves")
    print("  SPARK Validation")
    print("=" * 65)

    # ---- 1. MKM Polarization ----
    print("\n--- MKM Polarization Curve ---")
    t0 = time.time()
    U_mkm = np.linspace(0.0, -0.5, 101)
    mkm_results = compute_mkm_polarization(U_mkm, T=args.T)
    t_mkm = time.time() - t0
    print(f"  MKM done in {t_mkm:.1f} s")

    # Tafel slope analysis
    tafel_slope, _, _ = compute_tafel_slope(U_mkm, mkm_results['j_total'])
    if tafel_slope is not None:
        print(f"  Tafel slope = {tafel_slope:.1f} mV/decade")

    # Exchange current density (j at U=0 for HER at RHE)
    j0 = mkm_results['j_total'][0]
    print(f"  Exchange current density j₀ = {j0:.4f} mA/cm²")
    if J0_EXP_RANGE[0] <= j0 <= J0_EXP_RANGE[1]:
        print(f"  ✓ j₀ within experimental range "
              f"({J0_EXP_RANGE[0]}-{J0_EXP_RANGE[1]} mA/cm²)")
    else:
        print(f"  Note: j₀ outside experimental range "
              f"({J0_EXP_RANGE[0]}-{J0_EXP_RANGE[1]} mA/cm²)")
        print(f"  (Model uses simplified BV kinetics; "
              f"DFT barriers may need refinement)")

    # Rate-determining step
    tof_T_max = np.max(mkm_results['j_Tafel'])
    tof_H_max = np.max(mkm_results['j_Heyrovsky'])
    dominant = 'Tafel' if tof_T_max > tof_H_max else 'Heyrovsky'
    print(f"  Dominant pathway: {dominant} "
          f"(j_Tafel_max = {tof_T_max:.2e}, j_Hey_max = {tof_H_max:.2e})")

    save_data(mkm_results, os.path.join(outdir, 'her_mkm_polarization.dat'),
              'MKM')

    # ---- 2. KMC Polarization ----
    kmc_results = None
    if not args.no_kmc:
        print(f"\n--- KMC Polarization Curve ---")
        print(f"  Lattice: {args.kmc_size}×{args.kmc_size}, "
              f"equil: {args.kmc_equil}, prod: {args.kmc_prod} steps")
        print(f"  Lateral interactions: {'ON' if not args.no_lateral else 'OFF'}")

        U_kmc = np.linspace(0.0, -0.5, args.kmc_points)
        t0 = time.time()
        kmc_results = compute_kmc_polarization(
            U_kmc, T=args.T,
            lattice_size=args.kmc_size,
            equil_steps=args.kmc_equil,
            prod_steps=args.kmc_prod,
            with_lateral=not args.no_lateral,
        )
        t_kmc = time.time() - t0
        print(f"  KMC done in {t_kmc:.1f} s")

        save_data(kmc_results,
                  os.path.join(outdir, 'her_kmc_polarization.dat'), 'KMC')

        # ---- 3. MKM vs KMC comparison table ----
        print(f"\n--- MKM vs KMC Comparison ---")
        print(f"  {'U [V]':>8s}  {'θ_MKM':>8s}  {'θ_KMC':>8s}  "
              f"{'Δθ/θ%':>8s}  {'j_MKM':>12s}  {'j_KMC':>12s}  {'Δj/j%':>8s}")
        print(f"  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}  "
              f"{'-'*12}  {'-'*12}  {'-'*8}")
        for i, U in enumerate(U_kmc):
            # Find closest MKM point
            idx = np.argmin(np.abs(U_mkm - U))
            th_m = mkm_results['theta_H'][idx]
            th_k = kmc_results['theta_H'][i]
            j_m = mkm_results['j_total'][idx]
            j_k = kmc_results['j_total'][i]

            dth = abs(th_m - th_k) / max(th_m, 1e-10) * 100
            dj = abs(j_m - j_k) / max(abs(j_m), 1e-10) * 100

            print(f"  {U:>+8.3f}  {th_m:>8.4f}  {th_k:>8.4f}  "
                  f"{dth:>7.1f}%  {j_m:>12.4e}  {j_k:>12.4e}  {dj:>7.1f}%")

    # ---- 4. Plots ----
    print(f"\n--- Generating Figures ---")
    make_plots(mkm_results, kmc_results, outdir=outdir)

    # ---- Summary ----
    print("\n" + "=" * 65)
    print("  T4 Summary")
    print("=" * 65)
    print(f"  MKM polarization: 101 points, U = [0, -0.5] V")
    print(f"  Exchange current density: j₀ = {j0:.4f} mA/cm²")
    if tafel_slope is not None:
        print(f"  Tafel slope: {tafel_slope:.1f} mV/decade")
    print(f"  Dominant pathway: {dominant}")
    if kmc_results is not None:
        print(f"  KMC: {args.kmc_points} points, "
              f"{args.kmc_size}×{args.kmc_size} lattice")
    print(f"  Output directory: {outdir}/")
    print("=" * 65)


if __name__ == '__main__':
    main()
