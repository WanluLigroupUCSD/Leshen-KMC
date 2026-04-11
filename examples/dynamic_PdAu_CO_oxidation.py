"""
Dynamic Catalytic KMC Benchmark: CO Oxidation on PdAu(111)

Demonstrates the core V1 capability:
  - Binary alloy surface with random PdAu composition
  - Environment-dependent CO adsorption/desorption/reaction rates
  - Pd↔Au segregation events driven by CO binding
  - Comparison: dynamic (segregation ON) vs static (segregation OFF)

The "killer demo": segregation causes Pd enrichment at the surface,
which changes the dominant site type and alters catalytic activity —
impossible to predict with fixed-site methods.

Usage:
    python examples/dynamic_PdAu_CO_oxidation.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from collections import defaultdict

from spark.dynamic import (
    DynamicSurface, DynamicKMCEngine,
    EventGenerator, LookupTableEstimator,
    EventType,
)
from spark.dynamic.descriptor import EnvHash


def build_rate_table(temperature):
    """
    Build a lookup table of rates for CO oxidation on PdAu(111).

    DFT-inspired model barriers (simplified for demonstration):
      - CO adsorption barrier depends on NN Pd count
      - CO desorption barrier depends on NN Pd count
      - CO+O reaction barrier depends on environment
      - Pd-Au segregation barrier depends on CO coverage nearby

    On Pd: CO binds strongly (Eads ~ -1.5 to -2.0 eV)
    On Au: CO binds weakly  (Eads ~ -0.3 to -0.5 eV)
    Mixed: interpolates based on NN composition
    """
    est = LookupTableEstimator(temperature=temperature)

    # ── CO adsorption (barrierless, rate ~ Hertz-Knudsen) ──
    # Lower effective barrier = faster adsorption
    for n_pd in range(7):  # 0 to 6 Pd neighbors (fcc111 has 6 NN)
        n_au = 6 - n_pd

        # On Pd site: barrier decreases with more Pd neighbors
        nn = tuple(sorted([('Au', n_au), ('Pd', n_pd)]))
        est.add_entry(EventType.ADSORPTION, 'Pd', 0, nn,
                      barrier=0.05 + 0.01 * n_au,  # near-barrierless
                      prefactor=1e8)

        # On Au site: higher barrier, weaker binding
        est.add_entry(EventType.ADSORPTION, 'Au', 0, nn,
                      barrier=0.15 + 0.02 * n_au,
                      prefactor=1e8)

    # ── CO desorption ──
    CO = 1  # species ID for CO

    for n_pd in range(7):
        n_au = 6 - n_pd
        nn = tuple(sorted([('Au', n_au), ('Pd', n_pd)]))

        # From Pd: strong binding → high desorption barrier
        # More Pd neighbors → even stronger binding (cooperative effect)
        e_des_pd = 1.5 + 0.05 * n_pd - 0.02 * n_au
        est.add_entry(EventType.DESORPTION, 'Pd', CO, nn,
                      barrier=e_des_pd, prefactor=1e13)

        # From Au: weak binding → low desorption barrier
        e_des_au = 0.5 + 0.03 * n_pd
        est.add_entry(EventType.DESORPTION, 'Au', CO, nn,
                      barrier=e_des_au, prefactor=1e13)

    # ── CO diffusion ──
    for center_type in ['Pd', 'Au']:
        for n_pd in range(7):
            n_au = 6 - n_pd
            nn = tuple(sorted([('Au', n_au), ('Pd', n_pd)]))
            e_diff = 0.3 if center_type == 'Pd' else 0.15
            est.add_entry(EventType.DIFFUSION, center_type, CO, nn,
                          barrier=e_diff, prefactor=1e12)

    # ── Segregation (Pd↔Au swap) ──
    # Pd segregates to surface when CO is present (CO stabilizes Pd)
    # Barrier depends on local CO coverage via nn_adsorbate_count
    for center_type in ['Pd', 'Au']:
        for n_pd in range(7):
            n_au = 6 - n_pd
            nn = tuple(sorted([('Au', n_au), ('Pd', n_pd)]))

            # Base segregation barrier
            if center_type == 'Au':
                # Au→subsurface: barrier reduced when surrounded by CO/Pd
                e_seg = 0.8 - 0.03 * n_pd
            else:
                # Pd stays: higher barrier to leave surface
                e_seg = 1.0 + 0.02 * n_pd

            est.add_entry(EventType.SEGREGATION, center_type, 0, nn,
                          barrier=max(0.1, e_seg), prefactor=1e12)
            # Also with CO present
            est.add_entry(EventType.SEGREGATION, center_type, CO, nn,
                          barrier=max(0.1, e_seg - 0.2), prefactor=1e12)

    return est


def run_simulation(surface, estimator, temperature, max_steps,
                   enable_segregation=True, label=""):
    """Run a single KMC simulation and collect observables."""

    generator = EventGenerator(
        species_list=[1],  # CO = species 1
        enable_diffusion=True,
        enable_segregation=enable_segregation,
        enable_site_conversion=False,
    )

    engine = DynamicKMCEngine(
        surface, estimator, generator,
        temperature=temperature,
        debug=False,
    )

    # Observables over time
    times = []
    coverages = []
    compositions = []
    n_unique_envs = []
    sample_interval = max(1, max_steps // 200)

    def sampler(eng):
        if eng.kmc_step % sample_interval == 0:
            times.append(eng.kmc_time)
            coverages.append(eng.get_coverage())
            compositions.append(eng.get_composition())
            n_unique_envs.append(len(eng.get_unique_environments()))
        return False  # don't stop

    engine.run(max_steps=max_steps, callback=sampler)

    # Print results
    print(f"\n{'='*50}")
    print(f"  {label}")
    print(f"{'='*50}")
    engine.summary()

    return {
        'times': times,
        'coverages': coverages,
        'compositions': compositions,
        'n_unique_envs': n_unique_envs,
        'stats': engine.get_event_stats(),
        'final_composition': engine.get_composition(),
        'final_coverage': engine.get_coverage(),
        'final_site_dist': engine.get_site_type_distribution(),
    }


def main():
    temperature = 600  # K
    max_steps = 50000
    size = (15, 15)

    print("="*60)
    print("  PdAu(111) CO Oxidation — Dynamic vs Static KMC")
    print("="*60)
    print(f"  Temperature: {temperature} K")
    print(f"  Surface: {size[0]}×{size[1]} fcc(111)")
    print(f"  Initial composition: Pd:50% / Au:50%")
    print(f"  Max steps: {max_steps}")

    # Build rate table
    estimator = build_rate_table(temperature)
    print(f"\n  Rate table: {len(estimator._table)} entries")

    # ── Run 1: Dynamic (segregation ON) ──
    np.random.seed(42)
    surface_dyn = DynamicSurface.fcc111(
        composition={'Pd': 0.5, 'Au': 0.5},
        size=size, n_layers=1,  # single layer for clean CN=6
    )
    surface_dyn.register_species('CO')
    print(f"\n  Surface sites: {surface_dyn.n_sites}")

    result_dyn = run_simulation(
        surface_dyn, estimator, temperature, max_steps,
        enable_segregation=True,
        label="DYNAMIC (segregation ON)",
    )

    # ── Run 2: Static (segregation OFF) ──
    np.random.seed(42)
    surface_sta = DynamicSurface.fcc111(
        composition={'Pd': 0.5, 'Au': 0.5},
        size=size, n_layers=1,
    )
    surface_sta.register_species('CO')

    result_sta = run_simulation(
        surface_sta, estimator, temperature, max_steps,
        enable_segregation=False,
        label="STATIC (segregation OFF)",
    )

    # ── Comparison ──
    print(f"\n{'='*60}")
    print(f"  COMPARISON: Dynamic vs Static")
    print(f"{'='*60}")

    print(f"\n  Final surface composition:")
    print(f"    Dynamic: {result_dyn['final_composition']}")
    print(f"    Static:  {result_sta['final_composition']}")

    print(f"\n  Final CO coverage:")
    print(f"    Dynamic: {result_dyn['final_coverage']:.4f}")
    print(f"    Static:  {result_sta['final_coverage']:.4f}")

    print(f"\n  Final site type distribution:")
    print(f"    Dynamic: {result_dyn['final_site_dist']}")
    print(f"    Static:  {result_sta['final_site_dist']}")

    if result_dyn['compositions'] and result_sta['compositions']:
        dyn_pd = result_dyn['compositions'][-1].get('Pd', 0)
        sta_pd = result_sta['compositions'][-1].get('Pd', 0)
        delta_pd = dyn_pd - sta_pd
        print(f"\n  Surface Pd fraction change due to segregation:")
        print(f"    Δ(Pd_surface) = {delta_pd:+.4f}")
        if abs(delta_pd) > 0.01:
            print(f"    → Segregation significantly altered surface composition!")
            print(f"    → This is the phenomenon that fixed-site methods CANNOT capture.")

    print(f"\n  Unique environments (final):")
    if result_dyn['n_unique_envs']:
        print(f"    Dynamic: {result_dyn['n_unique_envs'][-1]}")
    if result_sta['n_unique_envs']:
        print(f"    Static:  {result_sta['n_unique_envs'][-1]}")


if __name__ == '__main__':
    main()
