"""Phase B.2 unit test — cross-layer KMC process firing.

Builds a tiny 2-layer model (PaperLayer + MetalLayer) on a 4x4 super-cell,
fires a deterministic cross-layer "destruct" process, and asserts that
species change correctly on BOTH layers at the right global site indices.

Run:
    /home/shidi/ai-chemist/catgo/.venv/bin/python tests/test_cross_layer_process.py
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from spark.types import Project, Site, Condition, Action


def build_two_layer_project():
    """A minimal 2-layer model:

    Layer 0 ('PaperOxide'): 1 site 'h' default 'X'  (oxide marker)
    Layer 1 ('Metal'):      1 site 'm' default 'frozen'   (Pd_atom blocker)

    Process 'destruct': remove X from PaperOxide, set Metal site to 'empty'.
    Process 'CO_ads': adsorb CO on Metal sites that are 'empty'.

    spuck = 2.
    """
    pt = Project()
    pt.set_meta(model_name='two_layer_test', model_dimension=2)
    pt.add_species(name='empty')
    pt.add_species(name='X')          # oxide species
    pt.add_species(name='frozen')     # disabled-metal sentinel
    pt.add_species(name='CO')

    paper = pt.add_layer(name='PaperOxide')
    paper.sites.append(Site(name='h', default_species='X'))
    metal = pt.add_layer(name='Metal')
    metal.sites.append(Site(name='m', default_species='frozen'))

    co_pap = pt.lattice.generate_coord('h.(0,0,0).PaperOxide')
    co_met = pt.lattice.generate_coord('m.(0,0,0).Metal')

    # Cross-layer destruct: PaperOxide.X -> empty AND Metal.frozen -> empty
    pt.add_process(
        name='destruct',
        conditions=[Condition(co_pap, 'X'), Condition(co_met, 'frozen')],
        actions=[Action(co_pap, 'empty'), Action(co_met, 'empty')],
        rate_constant='1e10',
    )

    # CO adsorption on Metal sites that are 'empty' (post-destruct)
    pt.add_process(
        name='CO_ads',
        conditions=[Condition(co_met, 'empty')],
        actions=[Action(co_met, 'CO')],
        rate_constant='1e8',
    )

    return pt


def test_two_layer_initial_state():
    """spuck=2, Metal sites default to 'frozen', PaperOxide default to 'X'."""
    from spark.engine import KMCEngine

    pt = build_two_layer_project()
    eng = KMCEngine(pt, size=[4, 4], banner=False, print_rates=False)

    assert eng.spuck == 2, f'expected spuck=2, got {eng.spuck}'
    assert eng.ncells == 16
    assert eng.nsites == 32

    # Check default tiling: even nr -> PaperOxide site (default 'X'),
    #                       odd  nr -> Metal site (default 'frozen').
    X_id = eng.species_id['X']
    frozen_id = eng.species_id['frozen']
    for nr in range(eng.nsites):
        s_in_cell = eng._site_in_cell(nr)
        if s_in_cell == 0:
            assert eng.lattice[nr] == X_id, \
                f'paper site {nr} should be X, got {eng.lattice[nr]}'
        else:
            assert eng.lattice[nr] == frozen_id, \
                f'metal site {nr} should be frozen, got {eng.lattice[nr]}'
    print('  initial state: PaperOxide=X (16 sites), Metal=frozen (16 sites) - OK')


def test_destruct_process_fires_correctly():
    """Run KMC and verify destruct flips paper X->empty AND metal frozen->empty
    at the same cell, never creating mismatched pairs."""
    from spark.engine import KMCEngine

    pt = build_two_layer_project()
    eng = KMCEngine(pt, size=[4, 4], banner=False, print_rates=False)

    X_id = eng.species_id['X']
    frozen_id = eng.species_id['frozen']
    empty_id = eng.species_id['empty']

    # Initial: 16 X, 16 frozen, 0 empty
    assert int((eng.lattice == X_id).sum()) == 16
    assert int((eng.lattice == frozen_id).sum()) == 16
    assert int((eng.lattice == empty_id).sum()) == 0

    eng.do_steps(50)

    n_X = int((eng.lattice == X_id).sum())
    n_frozen = int((eng.lattice == frozen_id).sum())
    n_empty = int((eng.lattice == empty_id).sum())
    n_CO = int((eng.lattice == eng.species_id['CO']).sum())
    print(f'  after 50 steps: X={n_X}, frozen={n_frozen}, '
          f'empty={n_empty}, CO={n_CO}, sum={n_X+n_frozen+n_empty+n_CO}')

    # Cross-layer invariant: paper X count and metal frozen count must
    # decrement together (every destruct removes one of each).
    n_destructs = 16 - n_X
    assert 16 - n_frozen >= n_destructs, \
        f'metal layer must lose at least {n_destructs} frozen, lost {16-n_frozen}'

    # Each cell's 2 sites should not have mismatched destruct: if paper is
    # empty, metal at the same cell must be empty or CO (post-CO_ads).
    for cell_id in range(eng.ncells):
        paper_nr = cell_id * eng.spuck + 0
        metal_nr = cell_id * eng.spuck + 1
        paper_sp = eng.lattice[paper_nr]
        metal_sp = eng.lattice[metal_nr]
        if paper_sp == empty_id:
            assert metal_sp in (empty_id, eng.species_id['CO']), \
                f'cell {cell_id}: paper destructed but metal={eng.species_names[metal_sp]}'
        elif paper_sp == X_id:
            assert metal_sp == frozen_id, \
                f'cell {cell_id}: paper still X but metal={eng.species_names[metal_sp]}'

    # At least one destruct should have fired (with rate 1e10 over 50 steps).
    assert n_X < 16, 'no destruct fired in 50 steps'
    print(f'  cross-layer invariant: {n_destructs} destructs, '
          f'all cell-pair states consistent - OK')


def test_total_species_conserved():
    """Total number of sites is conserved across multi-layer steps."""
    from spark.engine import KMCEngine

    pt = build_two_layer_project()
    eng = KMCEngine(pt, size=[5, 5], banner=False, print_rates=False)

    initial_total = eng.nsites
    eng.do_steps(200)
    final_total = sum(int((eng.lattice == i).sum()) for i in range(eng.nspecies))
    assert final_total == initial_total, \
        f'site count broke: {final_total} != {initial_total}'
    print(f'  conserved {initial_total} sites across 200 steps - OK')


if __name__ == '__main__':
    for name, fn in list(globals().items()):
        if name.startswith('test_') and callable(fn):
            print(f'--- {name} ---')
            fn()
            print(f'  PASS')
            print()
    print('Phase B.2 cross-layer process tests: ALL PASS')
