#!/usr/bin/env python3
"""
Validation test: CO adsorption/desorption on Pd(100).

Replicates the kmos tutorial "A First kMC Model" to verify
that spark produces correct results.

Expected behavior:
  - At T=600K, deltaG=-0.5eV: surface mostly covered with CO
  - At T=600K, deltaG=+0.1eV: surface mostly empty
  - Coverage follows Langmuir isotherm: theta = K*p / (1 + K*p)
    where K = exp(-deltaG/(kBT))
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from spark import Project, Site, Condition, Action, KMCEngine
from spark.units import kB, eV


def build_co_model():
    """Build the CO adsorption/desorption model from kmos tutorial."""
    pt = Project()
    pt.set_meta(author='Tutorial', model_name='CO_on_Pd100',
                model_dimension=2)

    # Species (first = default)
    pt.add_species(name='empty')
    pt.add_species(name='CO', color='#ff0000')

    # Lattice
    layer = pt.add_layer(name='simple_cubic')
    layer.sites.append(Site(name='hollow', pos=(0.5, 0.5, 0.5),
                            default_species='empty'))
    pt.lattice.cell = np.diag([3.5, 3.5, 10])

    # Parameters
    pt.add_parameter(name='T', value=600.0, adjustable=True,
                     min=400, max=800)
    pt.add_parameter(name='p_CO', value=1.0, adjustable=True,
                     min=1e-10, max=100)
    pt.add_parameter(name='A', value='(3.5*angstrom)**2')
    pt.add_parameter(name='deltaG', value=-0.5, adjustable=True,
                     min=-1.3, max=0.3)

    # Coordinate
    coord = pt.lattice.generate_coord('hollow')

    # CO adsorption
    pt.add_process(
        name='CO_adsorption',
        conditions=[Condition(coord, 'empty')],
        actions=[Action(coord, 'CO')],
        rate_constant='p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)',
    )

    # CO desorption (kmos convention: deltaG < 0 = exothermic adsorption)
    pt.add_process(
        name='CO_desorption',
        conditions=[Condition(coord, 'CO')],
        actions=[Action(coord, 'empty')],
        rate_constant=(
            'p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)'
            '*exp(beta*deltaG*eV)'
        ),
    )

    return pt


def langmuir_coverage(T, deltaG, p_CO=1.0):
    """Analytical Langmuir isotherm coverage."""
    K = np.exp(-deltaG * eV / (kB * T))
    return K * p_CO / (1.0 + K * p_CO)


def test_co_model():
    """Run validation test."""
    print("=" * 60)
    print("  Validation: CO adsorption/desorption on Pd(100)")
    print("=" * 60)

    pt = build_co_model()

    # Test at several deltaG values
    T = 600.0
    test_cases = [
        (-0.5, "exothermic, should have high CO coverage"),
        (-0.3, "moderate binding"),
        (-0.1, "weak binding"),
        (0.0,  "zero free energy"),
        (0.1,  "endothermic, should have low CO coverage"),
    ]

    print(f"\nT = {T} K, p_CO = 1 bar, Lattice = 30x30")
    print(f"{'deltaG [eV]':>12s} {'theta_KMC':>10s} {'theta_Langmuir':>15s} "
          f"{'Relative Error':>15s}")
    print("-" * 55)

    all_pass = True
    for deltaG, desc in test_cases:
        # Analytical
        theta_exact = langmuir_coverage(T, deltaG)

        # KMC
        engine = KMCEngine(pt, size=[30, 30], print_rates=False, banner=False)
        engine.parameters.T = T
        engine.parameters.deltaG = deltaG

        # Equilibrate
        engine.do_steps(500000)
        # Sample
        cov_samples = []
        for _ in range(20):
            engine.do_steps(10000)
            cov_samples.append(engine.get_coverage()['CO'])
        theta_kmc = np.mean(cov_samples)

        rel_err = abs(theta_kmc - theta_exact) / max(theta_exact, 1e-10)
        status = "PASS" if rel_err < 0.05 else "FAIL"
        if status == "FAIL":
            all_pass = False

        print(f"{deltaG:>12.1f} {theta_kmc:>10.4f} {theta_exact:>15.4f} "
              f"{rel_err:>14.2%}  {status}")

    print()
    if all_pass:
        print("All tests PASSED - KMC matches Langmuir isotherm!")
    else:
        print("Some tests FAILED - check implementation.")

    return all_pass


if __name__ == '__main__':
    test_co_model()
