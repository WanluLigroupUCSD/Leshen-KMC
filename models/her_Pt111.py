"""
HER on Pt(111): KMC + MKM model for Leshen-KMC validation.

Reaction network (acidic media):
  1. Volmer (fwd):  H+(aq) + e- + *  -> H*         (PCET, Ea=0.67 eV)
  2. Volmer (rev):  H*  -> H+(aq) + e- + *          (PCET, Ea=0.62 eV)
  3. Tafel:         H* + H*(adj) -> H2(g) + 2*      (thermal, Ea=0.85 eV)
  4. Heyrovsky:     H* + H+(aq) + e- -> H2(g) + *   (PCET, Ea=0.70 eV)
  5. Diffusion:     H*@s1 + *@s2 -> *@s1 + H*@s2    (thermal, Ea=0.10 eV)

Lateral: eps(H*-H*) = +0.10 eV (repulsive, NN)
BEP:     alpha = 0.5 for Volmer step

DFT sources:
  Li et al., ACS Catal. 14, 2696 (2024)
  Skulason et al., J. Phys. Chem. C 114, 18182 (2010)
  Karlberg et al., Phys. Rev. Lett. 99, 126101 (2007)
  Greeley & Mavrikakis, J. Phys. Chem. B 109, 3460 (2005)
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from mykmc.types import Project, Site, Layer, Condition, Action, Coord
from mykmc.rates import tst_rate, electrochemical_rate
from mykmc.units import kB, h, eV, e_charge

# ============================================================================
#  DFT Parameters
# ============================================================================

# Activation energies [eV] at zero coverage
EA_VOLMER_FWD = 0.67    # Li 2024
EA_VOLMER_REV = 0.62    # Li 2024
EA_TAFEL = 0.85         # Skulason 2010
EA_HEYROVSKY = 0.70     # Li 2024
EA_DIFFUSION = 0.10     # Greeley 2004

# Lateral interaction [eV]
EPS_HH = 0.10           # Karlberg 2007, repulsive NN

# BEP
ALPHA_BEP = 0.5         # Volmer step

# Electrochemistry
BETA_BV = 0.5           # Butler-Volmer symmetry factor
U0 = 0.0                # Equilibrium potential [V vs RHE]

# Lattice
A_PT = 2.775            # Pt(111) nearest-neighbor distance [Angstrom]
A_SITE = (A_PT * 1e-10) ** 2  # site area [m^2]


# ============================================================================
#  KMC Model (spatial, lattice-based)
# ============================================================================

def build_kmc_model():
    """
    Build the lattice KMC model for HER on Pt(111).

    Returns a Project object with:
      - 2 species (empty, H)
      - 5 process types (Volmer fwd/rev, Tafel, Heyrovsky, diffusion)
      - H-H lateral interaction
      - BEP relation for Volmer
    """
    pt = Project()
    pt.set_meta(
        author='Leshen-KMC',
        model_name='HER_Pt111',
        model_dimension=2,
    )

    # --- Species ---
    pt.add_species(name='empty', color='#ffffff')
    pt.add_species(name='H', color='#ff4444')

    # --- Lattice ---
    layer = pt.add_layer(name='Pt111')
    layer.sites.append(
        Site(name='top', pos=(0.5, 0.5, 0.5), default_species='empty'))
    pt.lattice.cell = np.diag([A_PT, A_PT, 15.0])

    # --- Parameters ---
    pt.add_parameter(name='T', value=298.0, adjustable=True, min=200, max=500)
    pt.add_parameter(name='U', value=-0.2, adjustable=True, min=-0.6, max=0.1)
    pt.add_parameter(name='beta_BV', value=BETA_BV)

    # Activation energies as adjustable parameters
    pt.add_parameter(name='Ea_Volmer_fwd', value=EA_VOLMER_FWD)
    pt.add_parameter(name='Ea_Volmer_rev', value=EA_VOLMER_REV)
    pt.add_parameter(name='Ea_Tafel', value=EA_TAFEL)
    pt.add_parameter(name='Ea_Heyrovsky', value=EA_HEYROVSKY)
    pt.add_parameter(name='Ea_diff', value=EA_DIFFUSION)

    # --- Coordinates ---
    coord = pt.lattice.generate_coord('top')
    # Neighbor coordinate for Tafel (adjacent site)
    coord_nn = Coord(offset=(1, 0, 0), site='top', layer='Pt111')

    # --- Processes ---

    # 1. Volmer forward: * -> H* (PCET)
    #    k = (kBT/h) * exp(-(Ea + beta_BV*(U - U0)) / (kBT))
    pt.add_process(
        name='Volmer_fwd',
        conditions=[Condition(coord, 'empty')],
        actions=[Action(coord, 'H')],
        rate_constant='kB*T/h*exp(-(Ea_Volmer_fwd + beta_BV*U)*eV/(kB*T))',
        tof_count={'Volmer_fwd': 1},
    )

    # 2. Volmer reverse: H* -> * (PCET)
    #    k = (kBT/h) * exp(-(Ea_rev - (1-beta_BV)*(U - U0)) / (kBT))
    #    = (kBT/h) * exp(-(Ea_rev - (1-beta_BV)*U) / (kBT))
    pt.add_process(
        name='Volmer_rev',
        conditions=[Condition(coord, 'H')],
        actions=[Action(coord, 'empty')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_Volmer_rev - (1-beta_BV)*U)*eV/(kB*T))'
        ),
        tof_count={'Volmer_rev': 1},
        reverse_of='Volmer_fwd',
    )

    # 3. Tafel: H* + H*(adj) -> 2* + H2(g) (thermal, bimolecular surface)
    pt.add_process(
        name='Tafel',
        conditions=[
            Condition(coord, 'H'),
            Condition(coord_nn, 'H'),
        ],
        actions=[
            Action(coord, 'empty'),
            Action(coord_nn, 'empty'),
        ],
        rate_constant='kB*T/h*exp(-Ea_Tafel*eV/(kB*T))',
        tof_count={'H2_Tafel': 1},
    )

    # 4. Heyrovsky: H* + H+(aq) + e- -> * + H2(g) (PCET)
    pt.add_process(
        name='Heyrovsky_fwd',
        conditions=[Condition(coord, 'H')],
        actions=[Action(coord, 'empty')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_Heyrovsky + beta_BV*U)*eV/(kB*T))'
        ),
        tof_count={'H2_Heyrovsky': 1},
    )

    # 5. H* diffusion (auto-generated for all NN directions)
    pt.add_diffusion(
        species='H',
        rate_constant='kB*T/h*exp(-Ea_diff*eV/(kB*T))',
        tof_count={'H_diff': 1},
    )

    # --- Lateral Interactions ---
    pt.add_lateral_interaction('H', 'H', energy=EPS_HH)

    # --- BEP Relations ---
    pt.add_bep_relation('Volmer_fwd', alpha=ALPHA_BEP)
    pt.add_bep_relation('Volmer_rev', alpha=1.0 - ALPHA_BEP)

    return pt


# ============================================================================
#  MKM Model (mean-field microkinetic)
# ============================================================================

def build_microkinetic_model():
    """
    Build the mean-field microkinetic model for HER on Pt(111).

    Same reaction network as KMC but without spatial resolution.
    Returns a MicroKineticModel object.
    """
    from mykmc.microkinetic import MicroKineticModel

    mkm = MicroKineticModel()
    mkm.add_species('H')

    mkm.parameters = {
        'T': 298.0,
        'U': -0.2,
        'beta_BV': BETA_BV,
    }

    # --- Rate functions ---

    def k_volmer_fwd(p):
        """Volmer forward: * -> H* (PCET)"""
        Ea = EA_VOLMER_FWD
        return electrochemical_rate(Ea, p['T'], p['U'], U0=0.0,
                                    beta_bv=p.get('beta_BV', 0.5))

    def k_volmer_rev(p):
        """Volmer reverse: H* -> * (PCET)"""
        Ea = EA_VOLMER_REV
        T, U = p['T'], p['U']
        beta_bv = p.get('beta_BV', 0.5)
        # Reverse PCET: Ea_eff = Ea_rev - (1-beta)*(U-U0)
        effective_Ea = Ea - (1.0 - beta_bv) * (U - U0)
        if effective_Ea < 0:
            effective_Ea = 0.0
        return (kB * T / h) * np.exp(-effective_Ea * eV / (kB * T))

    def k_tafel(p):
        """Tafel: H* + H* -> H2 (thermal)"""
        return tst_rate(EA_TAFEL, p['T'])

    def k_heyrovsky_fwd(p):
        """Heyrovsky: H* + H+ + e- -> H2 (PCET)"""
        Ea = EA_HEYROVSKY
        return electrochemical_rate(Ea, p['T'], p['U'], U0=0.0,
                                    beta_bv=p.get('beta_BV', 0.5))

    # --- Reactions ---

    # Volmer: empty + H+(aq) + e- <-> H*
    mkm.add_reaction(
        name='Volmer',
        reactants={'empty': 1},
        products={'H': 1},
        rate_fwd=k_volmer_fwd,
        rate_rev=k_volmer_rev,
        tof_count={'Volmer_net': 1},
    )

    # Tafel: 2H* -> H2(g) + 2*
    # Mean-field: rate = k_Tafel * theta_H^2
    mkm.add_reaction(
        name='Tafel',
        reactants={'H': 2},
        products={'empty': 2},
        rate_fwd=k_tafel,
        tof_count={'H2_Tafel': 1},
    )

    # Heyrovsky: H* + H+(aq) + e- -> H2(g) + *
    mkm.add_reaction(
        name='Heyrovsky',
        reactants={'H': 1},
        products={'empty': 1},
        rate_fwd=k_heyrovsky_fwd,
        tof_count={'H2_Heyrovsky': 1},
    )

    return mkm


# ============================================================================
#  Polarization Curve Builder
# ============================================================================

def build_polarization_curve(method='mkm'):
    """
    Build a PolarizationCurve object for HER on Pt(111).

    Parameters
    ----------
    method : str
        'mkm' for mean-field, 'kmc' for spatial KMC.

    Returns
    -------
    PolarizationCurve (if method='mkm')
    """
    from mykmc.polarization import PolarizationCurve

    if method == 'mkm':
        mkm = build_microkinetic_model()
        # HER: 2H+ + 2e- -> H2
        # Each H2 produced = 2 electrons transferred
        # TOF H2_Tafel and H2_Heyrovsky each produce one H2
        pc = PolarizationCurve(
            mkm=mkm,
            n_electrons={
                'H2_Tafel': 2,
                'H2_Heyrovsky': 2,
            },
            A_site=A_SITE,
        )
        return pc
    else:
        raise ValueError(
            "KMC polarization requires running the engine; "
            "use run_kmc_polarization() instead.")


def run_kmc_polarization(U_range, T=298.0, lattice_size=50,
                         equil_steps=1_000_000, prod_steps=1_000_000,
                         n_runs=1, verbose=True):
    """
    Compute polarization curve from KMC simulations.

    Runs KMC at each potential, extracts steady-state TOF, converts to j.

    Parameters
    ----------
    U_range : array-like
        Potentials [V vs RHE].
    T : float
        Temperature [K].
    lattice_size : int
        Lattice dimension (lattice_size x lattice_size).
    equil_steps : int
        Equilibration KMC steps (discarded).
    prod_steps : int
        Production KMC steps.
    n_runs : int
        Independent runs for error bars.
    verbose : bool
        Print progress.

    Returns
    -------
    dict with keys:
        'U'           : potential array
        'theta_H'     : H coverage array (mean over runs)
        'tof_Tafel'   : TOF_Tafel array
        'tof_Heyrovsky': TOF_Heyrovsky array
        'j_total'     : total current density [mA/cm^2]
    """
    from mykmc.engine import KMCEngine
    from mykmc.polarization import tof_to_current_density

    U_range = np.asarray(U_range, dtype=float)
    results = {
        'U': U_range,
        'theta_H': np.zeros(len(U_range)),
        'tof_Tafel': np.zeros(len(U_range)),
        'tof_Heyrovsky': np.zeros(len(U_range)),
        'j_Tafel': np.zeros(len(U_range)),
        'j_Heyrovsky': np.zeros(len(U_range)),
        'j_total': np.zeros(len(U_range)),
        'theta_H_std': np.zeros(len(U_range)),
    }

    for i, U in enumerate(U_range):
        theta_runs = []
        tof_T_runs = []
        tof_H_runs = []

        for run in range(n_runs):
            pt = build_kmc_model()
            engine = KMCEngine(pt, size=[lattice_size, lattice_size],
                               print_rates=False, banner=False)
            engine.parameters.T = T
            engine.parameters.U = U

            # Equilibrate
            engine.do_steps(equil_steps)

            # Reset TOF counters after equilibration
            engine.reset_tof()

            # Production
            engine.do_steps(prod_steps)

            cov = engine.get_coverage()
            tof = engine.get_tof()

            theta_runs.append(cov.get('H', 0.0))
            tof_T_runs.append(tof.get('H2_Tafel', 0.0))
            tof_H_runs.append(tof.get('H2_Heyrovsky', 0.0))

        # Average over runs
        results['theta_H'][i] = np.mean(theta_runs)
        results['theta_H_std'][i] = np.std(theta_runs) if n_runs > 1 else 0.0
        results['tof_Tafel'][i] = np.mean(tof_T_runs)
        results['tof_Heyrovsky'][i] = np.mean(tof_H_runs)

        # Convert to current density
        results['j_Tafel'][i] = tof_to_current_density(
            results['tof_Tafel'][i], n_electrons=2, A_site=A_SITE)
        results['j_Heyrovsky'][i] = tof_to_current_density(
            results['tof_Heyrovsky'][i], n_electrons=2, A_site=A_SITE)
        results['j_total'][i] = (results['j_Tafel'][i]
                                 + results['j_Heyrovsky'][i])

        if verbose:
            print(f"  U = {U:+.3f} V | theta_H = {results['theta_H'][i]:.4f}"
                  f" | j = {results['j_total'][i]:.4e} mA/cm^2"
                  f" (Tafel: {results['j_Tafel'][i]:.4e},"
                  f" Heyrovsky: {results['j_Heyrovsky'][i]:.4e})")

    return results


# ============================================================================
#  Quick Test / Demo
# ============================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='HER on Pt(111) simulation')
    parser.add_argument('--mode', choices=['kmc', 'mkm', 'both', 'test'],
                        default='test', help='Run mode')
    parser.add_argument('--T', type=float, default=298.0)
    parser.add_argument('--U', type=float, default=-0.2)
    parser.add_argument('--lattice-size', type=int, default=30)
    parser.add_argument('--steps', type=int, default=500000)
    args = parser.parse_args()

    if args.mode == 'test':
        # Quick sanity check
        print("=" * 60)
        print("  HER on Pt(111) — Model Sanity Check")
        print("=" * 60)

        # 1. Build and summarize KMC model
        print("\n--- KMC Model ---")
        pt = build_kmc_model()
        pt.summary()

        # 2. Build MKM and compute steady state
        print("\n--- MKM Steady State ---")
        mkm = build_microkinetic_model()

        for U in [0.0, -0.1, -0.2, -0.3, -0.4, -0.5]:
            mkm.parameters['U'] = U
            ss = mkm.solve_steady_state()
            tof = mkm.get_tof(ss)
            theta_H = ss.get('H', 0.0)
            tof_total = tof.get('H2_Tafel', 0.0) + tof.get('H2_Heyrovsky', 0.0)
            j = tof_total * 2 * e_charge / A_SITE * 0.1  # mA/cm^2
            print(f"  U = {U:+.2f} V | theta_H = {theta_H:.4f}"
                  f" | TOF_H2 = {tof_total:.4e} s^-1"
                  f" | j = {j:.4e} mA/cm^2")

        # 3. Quick KMC run
        print(f"\n--- Quick KMC (U={args.U} V, {args.lattice_size}x"
              f"{args.lattice_size}, {args.steps} steps) ---")
        from mykmc.engine import KMCEngine

        pt = build_kmc_model()
        engine = KMCEngine(pt, size=[args.lattice_size, args.lattice_size],
                           print_rates=True, banner=True)
        engine.parameters.T = args.T
        engine.parameters.U = args.U

        engine.do_steps(args.steps, progress=True)

        cov = engine.get_coverage()
        tof = engine.get_tof()
        print(f"\n  Coverage: theta_H = {cov.get('H', 0.0):.4f}")
        print(f"  TOF: {tof}")
        print(f"  Time: {engine.kmc_time:.4e} s")
        print(f"  Process stats: {engine.get_process_stats()}")

    elif args.mode == 'mkm':
        print("=" * 60)
        print("  HER Pt(111) — MKM Polarization Curve")
        print("=" * 60)

        pc = build_polarization_curve(method='mkm')
        U_range = np.linspace(-0.5, 0.0, 51)
        results = pc.compute(U_range, T=args.T)
        pc.print_results()
        pc.save('her_mkm_polarization.dat')

    elif args.mode == 'kmc':
        print("=" * 60)
        print("  HER Pt(111) — KMC Polarization Curve")
        print("=" * 60)

        U_range = np.linspace(-0.5, 0.0, 11)
        results = run_kmc_polarization(
            U_range, T=args.T,
            lattice_size=args.lattice_size,
            equil_steps=args.steps,
            prod_steps=args.steps,
        )

    elif args.mode == 'both':
        print("=" * 60)
        print("  HER Pt(111) — MKM vs KMC Comparison")
        print("=" * 60)

        # MKM
        print("\n--- MKM ---")
        pc = build_polarization_curve(method='mkm')
        U_range = np.linspace(-0.5, 0.0, 11)
        mkm_results = pc.compute(U_range, T=args.T)

        # KMC
        print("\n--- KMC ---")
        kmc_results = run_kmc_polarization(
            U_range, T=args.T,
            lattice_size=args.lattice_size,
            equil_steps=args.steps,
            prod_steps=args.steps,
        )

        # Comparison
        print("\n--- Comparison ---")
        print(f"{'U [V]':>8s}  {'theta_MKM':>10s}  {'theta_KMC':>10s}"
              f"  {'j_MKM':>12s}  {'j_KMC':>12s}")
        print("-" * 60)
        for i, U in enumerate(U_range):
            theta_mkm = mkm_results['coverages'][i].get('H', 0.0)
            theta_kmc = kmc_results['theta_H'][i]
            j_mkm = mkm_results['total_current'][i]
            j_kmc = kmc_results['j_total'][i]
            print(f"{U:>8.3f}  {theta_mkm:>10.4f}  {theta_kmc:>10.4f}"
                  f"  {j_mkm:>12.4e}  {j_kmc:>12.4e}")
