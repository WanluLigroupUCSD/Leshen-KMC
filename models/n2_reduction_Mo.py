"""
Kinetic Monte Carlo model for electrochemical N2 reduction on Mo surface.

Reaction network (from DFT transition state calculations):
  Distal pathway:   N2 -> NNH -> NNH2 -> N + NH3 -> NH -> NH2 -> NH3
  Alternating path:  N2 -> NNH -> HNNH -> HNNH2 -> H2NNH2 -> 2NH2 -> 2NH3

Elementary steps (all proton-coupled electron transfer unless noted):
  1.  N2(g) + * -> N2*                    (adsorption)
  2.  N2* + H+/e- -> NNH*                 (Ea = 2.82 eV)
  3a. NNH* + H+/e- -> HNNH*              (Ea = 1.14 eV, alternating)
  3b. NNH* + H+/e- -> NNH2*              (Ea = 3.22 eV, distal)
  4a. HNNH* + H+/e- -> HNNH2*            (Ea = 3.32 eV)
  4b. NNH2* -> N* + NH3(g)               (Ea = 2.94 eV, N-N cleavage)
  5a. HNNH2* + H+/e- -> H2NNH2*          (Ea = 3.36 eV)
  5b. N* + H+/e- -> NH*                  (Ea = 2.01 eV)
  6a. H2NNH2* -> 2NH2* (or NH2* + NH2(g)) (Ea = 5.14 eV, N-N cleavage)
  6b. NH* + H+/e- -> NH2*                (Ea = 2.68 eV)
  7.  NH2* + H+/e- -> NH3*               (Ea = 4.35 eV)
  8.  NH3* -> * + NH3(g)                  (desorption)
  9.  N2* -> * + N2(g)                    (desorption)

Activation barriers from DFT calculations on Mo surface.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from spark.types import Project, Site, Layer, Condition, Action, Coord
from spark.rates import tst_rate, electrochemical_rate, hertz_knudsen
from spark.units import bar, angstrom, pi, umass, kB, eV


def build_kmc_model():
    """
    Build the KMC model for N2 reduction on Mo surface.

    Returns a Project object ready for simulation.
    """
    pt = Project()
    pt.set_meta(
        author='Project4',
        model_name='N2_reduction_Mo',
        model_dimension=2,
    )

    # --- Species ---
    # First species is the default (empty site)
    pt.add_species(name='empty', color='#ffffff')
    pt.add_species(name='N2', color='#0000ff')
    pt.add_species(name='NNH', color='#0044ff')
    pt.add_species(name='HNNH', color='#0088ff')
    pt.add_species(name='NNH2', color='#00ccff')
    pt.add_species(name='HNNH2', color='#00ffcc')
    pt.add_species(name='H2NNH2', color='#00ff88')
    pt.add_species(name='N', color='#ff0000')
    pt.add_species(name='NH', color='#ff4400')
    pt.add_species(name='NH2', color='#ff8800')
    pt.add_species(name='NH3', color='#ffcc00')

    # --- Lattice ---
    layer = pt.add_layer(name='Mo_surface')
    layer.sites.append(Site(name='top', pos=(0.5, 0.5, 0.5),
                            default_species='empty'))
    pt.lattice.cell = np.diag([3.2, 3.2, 15.0])  # Mo surface unit cell

    # --- Parameters ---
    pt.add_parameter(name='T', value=300.0, adjustable=True,
                     min=200, max=500)
    pt.add_parameter(name='p_N2', value=1.0, adjustable=True,
                     min=1e-3, max=10.0)
    pt.add_parameter(name='U', value=-0.5, adjustable=True,
                     min=-2.0, max=0.5)  # Applied potential [V vs RHE]
    pt.add_parameter(name='A_site', value='(3.2*angstrom)**2')
    pt.add_parameter(name='E_bind_N2', value=0.8)   # N2 binding energy [eV]
    pt.add_parameter(name='E_bind_NH3', value=0.5)  # NH3 binding energy [eV]

    # --- DFT activation barriers [eV] ---
    pt.add_parameter(name='Ea_N2_to_NNH', value=2.82)
    pt.add_parameter(name='Ea_NNH_to_HNNH', value=1.14)
    pt.add_parameter(name='Ea_NNH_to_NNH2', value=3.22)
    pt.add_parameter(name='Ea_HNNH_to_HNNH2', value=3.32)
    pt.add_parameter(name='Ea_NNH2_to_N', value=2.94)
    pt.add_parameter(name='Ea_HNNH2_to_H2NNH2', value=3.36)
    pt.add_parameter(name='Ea_H2NNH2_to_NH2', value=5.14)
    pt.add_parameter(name='Ea_N_to_NH', value=2.01)
    pt.add_parameter(name='Ea_NH_to_NH2', value=2.68)
    pt.add_parameter(name='Ea_NH2_to_NH3', value=4.35)

    # Symmetry factor for Butler-Volmer
    pt.add_parameter(name='beta_BV', value=0.5)

    # --- Coordinate (single-site processes) ---
    coord = pt.lattice.generate_coord('top')

    # --- Processes ---

    # 1. N2 adsorption: * -> N2*
    # Rate: Hertz-Knudsen impingement
    pt.add_process(
        name='N2_adsorption',
        conditions=[Condition(coord, 'empty')],
        actions=[Action(coord, 'N2')],
        rate_constant='p_N2*bar*A_site/sqrt(2*pi*m_N2*umass/beta)',
    )

    # 2. N2 desorption: N2* -> *
    pt.add_process(
        name='N2_desorption',
        conditions=[Condition(coord, 'N2')],
        actions=[Action(coord, 'empty')],
        rate_constant=(
            'p_N2*bar*A_site/sqrt(2*pi*m_N2*umass/beta)'
            '*exp(-beta*E_bind_N2*eV)'
        ),
    )

    # 3. N2* + H+/e- -> NNH* (first hydrogenation)
    pt.add_process(
        name='N2_to_NNH',
        conditions=[Condition(coord, 'N2')],
        actions=[Action(coord, 'NNH')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_N2_to_NNH + beta_BV*U)*eV/(kB*T))'
        ),
        tof_count={'N2_consumption': 1},
    )

    # 4a. NNH* + H+/e- -> HNNH* (alternating pathway)
    pt.add_process(
        name='NNH_to_HNNH',
        conditions=[Condition(coord, 'NNH')],
        actions=[Action(coord, 'HNNH')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_NNH_to_HNNH + beta_BV*U)*eV/(kB*T))'
        ),
    )

    # 4b. NNH* + H+/e- -> NNH2* (distal pathway)
    pt.add_process(
        name='NNH_to_NNH2',
        conditions=[Condition(coord, 'NNH')],
        actions=[Action(coord, 'NNH2')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_NNH_to_NNH2 + beta_BV*U)*eV/(kB*T))'
        ),
    )

    # 5a. HNNH* + H+/e- -> HNNH2*
    pt.add_process(
        name='HNNH_to_HNNH2',
        conditions=[Condition(coord, 'HNNH')],
        actions=[Action(coord, 'HNNH2')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_HNNH_to_HNNH2 + beta_BV*U)*eV/(kB*T))'
        ),
    )

    # 5b. NNH2* -> N* + NH3(g) (N-N bond cleavage, releases NH3)
    pt.add_process(
        name='NNH2_to_N_NH3',
        conditions=[Condition(coord, 'NNH2')],
        actions=[Action(coord, 'N')],
        rate_constant='kB*T/h*exp(-Ea_NNH2_to_N*eV/(kB*T))',
        tof_count={'NH3_production': 1},
    )

    # 6a. HNNH2* + H+/e- -> H2NNH2*
    pt.add_process(
        name='HNNH2_to_H2NNH2',
        conditions=[Condition(coord, 'HNNH2')],
        actions=[Action(coord, 'H2NNH2')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_HNNH2_to_H2NNH2 + beta_BV*U)*eV/(kB*T))'
        ),
    )

    # 6b. H2NNH2* -> NH2* + NH2(g) (N-N bond cleavage)
    # Simplified: one NH2 stays on surface, other released
    pt.add_process(
        name='H2NNH2_to_NH2',
        conditions=[Condition(coord, 'H2NNH2')],
        actions=[Action(coord, 'NH2')],
        rate_constant='kB*T/h*exp(-Ea_H2NNH2_to_NH2*eV/(kB*T))',
    )

    # 7. N* + H+/e- -> NH*
    pt.add_process(
        name='N_to_NH',
        conditions=[Condition(coord, 'N')],
        actions=[Action(coord, 'NH')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_N_to_NH + beta_BV*U)*eV/(kB*T))'
        ),
    )

    # 8. NH* + H+/e- -> NH2*
    pt.add_process(
        name='NH_to_NH2',
        conditions=[Condition(coord, 'NH')],
        actions=[Action(coord, 'NH2')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_NH_to_NH2 + beta_BV*U)*eV/(kB*T))'
        ),
    )

    # 9. NH2* + H+/e- -> NH3*
    pt.add_process(
        name='NH2_to_NH3',
        conditions=[Condition(coord, 'NH2')],
        actions=[Action(coord, 'NH3')],
        rate_constant=(
            'kB*T/h*exp(-(Ea_NH2_to_NH3 + beta_BV*U)*eV/(kB*T))'
        ),
    )

    # 10. NH3* -> * + NH3(g) (product desorption)
    pt.add_process(
        name='NH3_desorption',
        conditions=[Condition(coord, 'NH3')],
        actions=[Action(coord, 'empty')],
        rate_constant='kB*T/h*exp(-E_bind_NH3*eV/(kB*T))',
        tof_count={'NH3_production': 1},
    )

    return pt


def build_microkinetic_model():
    """
    Build the mean-field microkinetic model for N2 reduction on Mo surface.

    Returns a MicroKineticModel object.
    """
    from spark.microkinetic import MicroKineticModel
    from spark.rates import tst_rate, electrochemical_rate

    mkm = MicroKineticModel()

    # Surface species
    for sp in ['N2', 'NNH', 'HNNH', 'NNH2', 'HNNH2',
               'H2NNH2', 'N', 'NH', 'NH2', 'NH3']:
        mkm.add_species(sp)

    # Default parameters
    mkm.parameters = {
        'T': 300.0,
        'p_N2': 1.0,       # bar
        'U': -0.5,         # V vs RHE
        'beta_BV': 0.5,    # Butler-Volmer symmetry factor
        'E_bind_N2': 0.8,  # eV
        'E_bind_NH3': 0.5, # eV
    }

    # DFT barriers [eV]
    barriers = {
        'N2_to_NNH': 2.82,
        'NNH_to_HNNH': 1.14,
        'NNH_to_NNH2': 3.22,
        'HNNH_to_HNNH2': 3.32,
        'NNH2_to_N': 2.94,
        'HNNH2_to_H2NNH2': 3.36,
        'H2NNH2_to_NH2': 5.14,
        'N_to_NH': 2.01,
        'NH_to_NH2': 2.68,
        'NH2_to_NH3': 4.35,
    }

    # Helper: electrochemical rate with potential dependence
    def echem_rate(Ea_key):
        def rate_func(p):
            Ea = barriers[Ea_key]
            T = p['T']
            U = p['U']
            beta_bv = p.get('beta_BV', 0.5)
            return electrochemical_rate(Ea, T, U, U0=0.0, beta_bv=beta_bv)
        return rate_func

    # Helper: thermal rate (no potential dependence)
    def thermal_rate(Ea_key):
        def rate_func(p):
            return tst_rate(barriers[Ea_key], p['T'])
        return rate_func

    # N2 adsorption rate
    def k_ads_N2(p):
        T = p['T']
        p_N2 = p['p_N2']
        A = (3.2e-10) ** 2
        m_N2 = 28.014 * 1.66054e-27
        return hertz_knudsen(p_N2 * 1e5, T, m_N2, A)

    # N2 desorption rate
    def k_des_N2(p):
        T = p['T']
        E_bind = p.get('E_bind_N2', 0.8)
        return tst_rate(E_bind, T)

    # NH3 desorption rate
    def k_des_NH3(p):
        T = p['T']
        E_bind = p.get('E_bind_NH3', 0.5)
        return tst_rate(E_bind, T)

    # --- Reactions ---

    # 1. N2 adsorption/desorption
    mkm.add_reaction(
        name='N2_adsorption',
        reactants={'empty': 1},
        products={'N2': 1},
        rate_fwd=k_ads_N2,
        rate_rev=k_des_N2,
    )

    # 2. N2* -> NNH* (first PCET)
    mkm.add_reaction(
        name='N2_to_NNH',
        reactants={'N2': 1},
        products={'NNH': 1},
        rate_fwd=echem_rate('N2_to_NNH'),
        tof_count={'N2_consumption': 1},
    )

    # 3a. NNH* -> HNNH* (alternating)
    mkm.add_reaction(
        name='NNH_to_HNNH',
        reactants={'NNH': 1},
        products={'HNNH': 1},
        rate_fwd=echem_rate('NNH_to_HNNH'),
    )

    # 3b. NNH* -> NNH2* (distal)
    mkm.add_reaction(
        name='NNH_to_NNH2',
        reactants={'NNH': 1},
        products={'NNH2': 1},
        rate_fwd=echem_rate('NNH_to_NNH2'),
    )

    # 4a. HNNH* -> HNNH2*
    mkm.add_reaction(
        name='HNNH_to_HNNH2',
        reactants={'HNNH': 1},
        products={'HNNH2': 1},
        rate_fwd=echem_rate('HNNH_to_HNNH2'),
    )

    # 4b. NNH2* -> N* + NH3(g) (N-N cleavage, releases NH3)
    mkm.add_reaction(
        name='NNH2_to_N',
        reactants={'NNH2': 1},
        products={'N': 1},
        rate_fwd=thermal_rate('NNH2_to_N'),
        tof_count={'NH3_production': 1},
    )

    # 5a. HNNH2* -> H2NNH2*
    mkm.add_reaction(
        name='HNNH2_to_H2NNH2',
        reactants={'HNNH2': 1},
        products={'H2NNH2': 1},
        rate_fwd=echem_rate('HNNH2_to_H2NNH2'),
    )

    # 5b. H2NNH2* -> NH2* (+ NH2 gas, N-N cleavage)
    mkm.add_reaction(
        name='H2NNH2_to_NH2',
        reactants={'H2NNH2': 1},
        products={'NH2': 1},
        rate_fwd=thermal_rate('H2NNH2_to_NH2'),
    )

    # 6. N* -> NH*
    mkm.add_reaction(
        name='N_to_NH',
        reactants={'N': 1},
        products={'NH': 1},
        rate_fwd=echem_rate('N_to_NH'),
    )

    # 7. NH* -> NH2*
    mkm.add_reaction(
        name='NH_to_NH2',
        reactants={'NH': 1},
        products={'NH2': 1},
        rate_fwd=echem_rate('NH_to_NH2'),
    )

    # 8. NH2* -> NH3*
    mkm.add_reaction(
        name='NH2_to_NH3',
        reactants={'NH2': 1},
        products={'NH3': 1},
        rate_fwd=echem_rate('NH2_to_NH3'),
    )

    # 9. NH3 desorption: NH3* -> *
    mkm.add_reaction(
        name='NH3_desorption',
        reactants={'NH3': 1},
        products={'empty': 1},
        rate_fwd=k_des_NH3,
        tof_count={'NH3_production': 1},
    )

    return mkm


def build_dft_microkinetic_model(dft_data):
    """
    Build microkinetic model using potential-dependent barriers from
    constant-potential DFT calculations.

    Parameters
    ----------
    dft_data : dict or str
        If str, path to a JSON file (see load_energy_data).
        If dict, must contain:
          'potentials'     : list of U values [V vs RHE]
          'state_energies' : {state_name: [E(U1), E(U2), ...]} in eV
          'ts_energies'    : {(IS, FS): [E_TS(U1), ...]} in eV
                             or {"IS->FS": [...]} in JSON format

    Returns
    -------
    PolarizationCurve
        Ready to call .compute(U_range, T).

    Example JSON file (n2_dft_energies.json)::

        {
          "potentials": [-2.0, -1.5, -1.0, -0.5, 0.0],
          "state_energies": {
            "clean":   [0.00, 0.00, 0.00, 0.00, 0.00],
            "N2*":     [-0.80, -0.72, -0.60, -0.50, -0.38],
            "NNH*":    [1.50, 1.65, 1.82, 2.00, 2.20],
            "HNNH*":   [1.20, 1.30, 1.45, 1.62, 1.80],
            "NNH2*":   [1.80, 1.95, 2.15, 2.38, 2.60],
            "N*":      [0.10, 0.10, 0.10, 0.10, 0.10],
            "NH*":     [0.50, 0.55, 0.62, 0.70, 0.80],
            "NH2*":    [0.80, 0.88, 0.98, 1.10, 1.25],
            "NH3*":    [0.30, 0.32, 0.35, 0.38, 0.42]
          },
          "ts_energies": {
            "N2*->NNH*":    [1.80, 2.00, 2.30, 2.60, 3.00],
            "NNH*->HNNH*":  [2.20, 2.40, 2.65, 2.90, 3.15],
            "NNH*->NNH2*":  [2.50, 2.72, 3.00, 3.30, 3.65],
            "NNH2*->N*":    [3.00, 3.10, 3.20, 3.35, 3.50],
            "N*->NH*":      [1.50, 1.65, 1.85, 2.10, 2.40],
            "NH*->NH2*":    [2.00, 2.18, 2.40, 2.65, 2.95],
            "NH2*->NH3*":   [2.80, 3.00, 3.28, 3.60, 3.98]
          }
        }
    """
    from spark.polarization import (
        EnergyLandscape, PolarizationCurve, load_energy_data,
    )
    from spark.microkinetic import MicroKineticModel
    from spark.rates import tst_rate, hertz_knudsen

    # Load data
    if isinstance(dft_data, str):
        dft_data = load_energy_data(dft_data)

    # Convert string keys "IS->FS" to tuples if needed
    ts_energies = {}
    for key, vals in dft_data['ts_energies'].items():
        if isinstance(key, str) and '->' in key:
            is_name, fs_name = key.split('->')
            ts_energies[(is_name.strip(), fs_name.strip())] = vals
        else:
            ts_energies[key] = vals

    landscape = EnergyLandscape(
        potentials=dft_data['potentials'],
        state_energies=dft_data['state_energies'],
        ts_energies=ts_energies,
    )

    # Build microkinetic model
    mkm = MicroKineticModel()
    species_list = ['N2', 'NNH', 'HNNH', 'NNH2', 'N', 'NH', 'NH2', 'NH3']
    for sp in species_list:
        mkm.add_species(sp)

    mkm.parameters = {
        'T': 300.0,
        'U': -0.5,
        'p_N2': 1.0,
        'E_bind_N2': 0.8,
        'E_bind_NH3': 0.5,
    }

    # N2 adsorption/desorption (gas-phase kinetics, not from DFT landscape)
    def k_ads_N2(p):
        T = p['T']
        p_N2 = p['p_N2']
        A = (3.2e-10) ** 2
        m_N2 = 28.014 * 1.66054e-27
        return hertz_knudsen(p_N2 * 1e5, T, m_N2, A)

    def k_des_N2(p):
        return tst_rate(p.get('E_bind_N2', 0.8), p['T'])

    def k_des_NH3(p):
        return tst_rate(p.get('E_bind_NH3', 0.5), p['T'])

    mkm.add_reaction(
        name='N2_adsorption',
        reactants={'empty': 1}, products={'N2': 1},
        rate_fwd=k_ads_N2, rate_rev=k_des_N2,
    )

    # Electrochemical steps with interpolated DFT barriers
    echem_steps = []
    step_map = {
        ('N2*', 'NNH*'):   ('N2_to_NNH',   {'N2': 1},   {'NNH': 1},
                             {'N2_consumption': 1}),
        ('NNH*', 'HNNH*'): ('NNH_to_HNNH',  {'NNH': 1},  {'HNNH': 1}, {}),
        ('NNH*', 'NNH2*'): ('NNH_to_NNH2',  {'NNH': 1},  {'NNH2': 1}, {}),
        ('NNH2*', 'N*'):    ('NNH2_to_N',    {'NNH2': 1}, {'N': 1},
                             {'NH3_production': 1}),
        ('N*', 'NH*'):      ('N_to_NH',      {'N': 1},    {'NH': 1},   {}),
        ('NH*', 'NH2*'):    ('NH_to_NH2',    {'NH': 1},   {'NH2': 1},  {}),
        ('NH2*', 'NH3*'):   ('NH2_to_NH3',   {'NH2': 1},  {'NH3': 1},  {}),
    }

    for key, (name, reactants, products, tof_count) in step_map.items():
        if key in ts_energies:
            rate_fwd, rate_rev = landscape.get_rate_funcs(*key)
            mkm.add_reaction(
                name=name,
                reactants=reactants, products=products,
                rate_fwd=rate_fwd, rate_rev=rate_rev,
                tof_count=tof_count,
            )

    # NH3 desorption
    mkm.add_reaction(
        name='NH3_desorption',
        reactants={'NH3': 1}, products={'empty': 1},
        rate_fwd=k_des_NH3,
        tof_count={'NH3_production': 1},
    )

    # N2 overall: N2 + 6H+ + 6e- -> 2NH3, so 6 electrons per NH3 pair
    # But our TOF counts per NH3 molecule, so 3 electrons per count
    # (depends on how tof_count is defined; here NH3_production counts each NH3)
    pc = PolarizationCurve(
        mkm=mkm,
        n_electrons={'NH3_production': 3},  # 3 e- per NH3 molecule
        A_site=(3.2e-10) ** 2,
    )

    return pc


if __name__ == '__main__':
    print("Building KMC model...")
    pt = build_kmc_model()
    pt.summary()

    print("\nBuilding microkinetic model...")
    mkm = build_microkinetic_model()
    print(f"Species: {mkm.species}")
    print(f"Reactions: {[r['name'] for r in mkm.reactions]}")
