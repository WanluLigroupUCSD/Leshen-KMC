"""SPARK multi-lattice example — Pd(100) / sqrt5-PdO surface oxide reduction.

A reduced-scope reproduction of:
  Hoffmann, Scheffler, Reuter — Multi-Lattice Kinetic Monte Carlo Simulations
  from First Principles: Reduction of the Pd(100) Surface Oxide by CO.
  ACS Catal. 2015, 5, 1199-1209. DOI 10.1021/cs501352t

This script demonstrates SPARK's new multi-lattice capability (Hoffmann et al.
algorithm, paper §2.1-2.3): two commensurate sub-lattices coexisting on a
shared super-cell, with cross-layer "destruct" elementary processes that
transition local regions from the oxide phase to the metal phase.

Algorithmic correspondence with the paper:
  - paper "super-lattice" -> SPARK ``Lattice`` containing two ``Layer`` objects
    sharing one ``cell``
  - paper "null" species occupying disabled sites -> ordinary species
    ('frozen') marking the not-yet-reduced metal layer (kmcos convention,
    see docs/multi_lattice_design.md §1)
  - paper "lattice-swap elementary process" -> the ``destruct`` and
    ``oxidize`` processes below, whose conditions/actions cross layer
    boundaries
  - paper §2.3 lateral interactions baked into rate constants -> SPARK
    ``LateralInteraction`` (handled in ``spark.engine``, no change needed
    for multi-lattice)

This is a STRUCTURAL toy: the rate constants are illustrative, not the
40+-process DFT-parameterized rate set of the original paper. The goal is
to verify that SPARK's multi-lattice machinery fires the right events in
the right cells, not to reproduce the experimental Pd(100)/PdO reduction
temperature dependence.

Run:
    /home/shidi/ai-chemist/catgo/.venv/bin/python examples/multi_lattice_PdO_reduction.py
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from spark.types import Project, Site, Condition, Action
from spark.engine import KMCEngine


def build_pd_pdo_project():
    """Toy Pd(100) + sqrt5-PdO multi-lattice model."""
    pt = Project()
    pt.set_meta(
        author='SPARK',
        model_name='Pd100_PdO_multilattice_toy',
        model_dimension=2,
    )

    # --- Species ---
    # 'frozen' plays the kmcos convention "disabled-site marker" role
    # (no process accepts 'frozen' as a condition, hence sites stay locked
    # until a destruct/oxidize cross-layer process flips them).
    pt.add_species(name='empty', color='#ffffff')
    pt.add_species(name='CO',    color='#444444')
    pt.add_species(name='O',     color='#cc2222')
    pt.add_species(name='frozen', color='#666666')   # kmcos "Pd" sentinel

    # --- Layer 1: Pd(100) metal surface ---
    # Per cell: 1 hollow + 1 bridge site. Default species 'frozen' meaning
    # "the metal layer is currently masked by the oxide". Reduction process
    # flips these to 'empty', after which CO/O can adsorb.
    pd100 = pt.add_layer(name='Pd100')
    pd100.sites.append(Site(name='hollow', default_species='frozen'))
    pd100.sites.append(Site(name='bridge', default_species='frozen'))

    # --- Layer 2: sqrt5-PdO surface oxide ---
    # Per cell: 1 oxide-bridge + 1 oxide-hollow + 1 oxide-O site. Default
    # 'O' on the O site (the lattice O atom that gets removed by CO_destruct),
    # 'empty' on the bridge/hollow (CO can adsorb here even on intact oxide).
    pdo = pt.add_layer(name='PdO')
    pdo.sites.append(Site(name='bridge', default_species='empty'))
    pdo.sites.append(Site(name='hollow', default_species='empty'))
    pdo.sites.append(Site(name='Olat',   default_species='O'))

    # spuck = 2 + 3 = 5 sites per cell

    # --- Coords ---
    c_pd_h  = pt.lattice.generate_coord('hollow.(0,0,0).Pd100')
    c_pd_b  = pt.lattice.generate_coord('bridge.(0,0,0).Pd100')
    c_pdo_b = pt.lattice.generate_coord('bridge.(0,0,0).PdO')
    c_pdo_h = pt.lattice.generate_coord('hollow.(0,0,0).PdO')
    c_pdo_O = pt.lattice.generate_coord('Olat.(0,0,0).PdO')

    # --- Surface chemistry on the oxide (paper §2.3) ---
    pt.add_process(
        name='CO_ads_PdO_bridge',
        conditions=[Condition(c_pdo_b, 'empty')],
        actions=[Action(c_pdo_b, 'CO')],
        rate_constant='1e8',
    )
    pt.add_process(
        name='CO_des_PdO_bridge',
        conditions=[Condition(c_pdo_b, 'CO')],
        actions=[Action(c_pdo_b, 'empty')],
        rate_constant='1e7',
    )
    pt.add_process(
        name='CO_ads_PdO_hollow',
        conditions=[Condition(c_pdo_h, 'empty')],
        actions=[Action(c_pdo_h, 'CO')],
        rate_constant='5e7',
    )

    # --- Cross-layer destruct (paper §2.2 "lattice-swap process") ---
    # The hallmark multi-lattice operation: CO at PdO bridge consumes the
    # adjacent lattice-O atom, removes itself as CO2(g), and uncovers the
    # underlying Pd(100) hollow + bridge sites (frozen -> empty).
    pt.add_process(
        name='destruct_CO_oxidizes_Olat',
        conditions=[
            Condition(c_pdo_b, 'CO'),
            Condition(c_pdo_O, 'O'),
            Condition(c_pd_h,  'frozen'),
            Condition(c_pd_b,  'frozen'),
        ],
        actions=[
            Action(c_pdo_b, 'empty'),    # CO leaves as CO2
            Action(c_pdo_O, 'empty'),    # lattice O removed
            Action(c_pd_h,  'empty'),    # uncover Pd hollow
            Action(c_pd_b,  'empty'),    # uncover Pd bridge
        ],
        rate_constant='1e9',
    )

    # --- Surface chemistry on uncovered Pd(100) ---
    pt.add_process(
        name='CO_ads_Pd100_bridge',
        conditions=[Condition(c_pd_b, 'empty')],
        actions=[Action(c_pd_b, 'CO')],
        rate_constant='1e8',
    )
    pt.add_process(
        name='CO_des_Pd100_bridge',
        conditions=[Condition(c_pd_b, 'CO')],
        actions=[Action(c_pd_b, 'empty')],
        rate_constant='1e6',
    )

    # --- Reverse oxidize (slow, included for detailed-balance demo) ---
    pt.add_process(
        name='oxidize_reseal',
        conditions=[
            Condition(c_pdo_b, 'empty'),
            Condition(c_pdo_O, 'empty'),
            Condition(c_pd_h,  'empty'),
            Condition(c_pd_b,  'empty'),
        ],
        actions=[
            Action(c_pdo_b, 'empty'),
            Action(c_pdo_O, 'O'),
            Action(c_pd_h,  'frozen'),
            Action(c_pd_b,  'frozen'),
        ],
        rate_constant='1e3',
    )

    return pt


def species_histogram(eng, label):
    counts = {n: int((eng.lattice == i).sum())
              for i, n in enumerate(eng.species_names)}
    total = sum(counts.values())
    print(f'  {label}: {counts}, total={total}')


def main():
    pt = build_pd_pdo_project()
    print(f'Project built: spuck = {pt.lattice.spuck}, '
          f'{len(pt.lattice.layers)} layers, '
          f'{len(pt.process_list)} processes')

    L = 8
    eng = KMCEngine(pt, size=[L, L], banner=True, print_rates=False)
    print(f'Engine: ncells={eng.ncells}, spuck={eng.spuck}, '
          f'nsites={eng.nsites}')

    # Initial state: every cell has 2 'frozen' Pd sites + 2 'empty' PdO ads
    # sites + 1 'O' lattice site = 5 sites per cell
    species_histogram(eng, 'initial')
    n_O_init = int((eng.lattice == eng.species_id['O']).sum())
    n_frozen_init = int((eng.lattice == eng.species_id['frozen']).sum())
    assert n_O_init == eng.ncells, 'every cell must start with 1 lattice O'
    assert n_frozen_init == 2 * eng.ncells, \
        'every cell must start with 2 frozen Pd sites'

    # Snapshot at intervals
    for step_target in [1000, 5000, 20000, 50000]:
        steps_to_run = step_target - eng.kmc_step
        if steps_to_run > 0:
            eng.do_steps(steps_to_run)
        species_histogram(eng, f'step={eng.kmc_step}, t={eng.kmc_time:.3e}s')

    # Reduction extent: fraction of cells whose lattice O is gone
    n_O_final = int((eng.lattice == eng.species_id['O']).sum())
    reduction_fraction = 1 - n_O_final / eng.ncells
    print(f'\nReduction extent: {reduction_fraction*100:.1f}% '
          f'({eng.ncells - n_O_final}/{eng.ncells} cells reduced)')

    # Process firings — confirm cross-layer destruct actually fired
    print('\nProcess firings:')
    for pid, name in enumerate(eng.process_names):
        n_fired = int(eng.procstat[pid])
        if n_fired > 0:
            print(f'  {name}: {n_fired}')

    # Cross-layer invariant: every cell's lattice-O state must match its
    # underlying Pd hollow state — if O is gone, hollow is uncovered
    # (empty/CO); if O is present, hollow is still frozen.
    o_id = eng.species_id['O']
    frozen_id = eng.species_id['frozen']
    spuck = eng.spuck
    pd_hollow_offset = pt.lattice.site_in_cell_id('Pd100', 'hollow')
    pdo_O_offset     = pt.lattice.site_in_cell_id('PdO', 'Olat')
    n_consistent = 0
    n_violations = 0
    for cid in range(eng.ncells):
        o_state  = eng.lattice[cid * spuck + pdo_O_offset]
        pd_state = eng.lattice[cid * spuck + pd_hollow_offset]
        if o_state == o_id:
            if pd_state == frozen_id:
                n_consistent += 1
            else:
                n_violations += 1
        else:  # O removed
            if pd_state != frozen_id:
                n_consistent += 1
            else:
                n_violations += 1
    print(f'\nCross-layer invariant check: '
          f'{n_consistent}/{eng.ncells} consistent, {n_violations} violations')
    assert n_violations == 0, \
        'Multi-lattice invariant broken: lattice O removed without Pd uncovered'

    print('\nMulti-lattice demo: PASS')
    return eng


if __name__ == '__main__':
    main()
