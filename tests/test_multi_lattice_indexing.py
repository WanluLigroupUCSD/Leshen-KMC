"""Phase A.2 unit tests — multi-lattice site indexing bijection.

Verifies:
  1. Round-trip lattice_to_nr ∘ nr_to_lattice = identity on every nr in
     [0, Lx*Ly*Lz*spuck) for both single-layer and multi-layer projects.
  2. Layer/site-name → site_in_cell resolution.
  3. offset_nr correctness across PBCs and across layer boundaries.
  4. Single-layer mode produces site indices identical to the natural
     1D embedding (no regression for existing single-layer engine code).

Run:
    /home/shidi/ai-chemist/catgo/.venv/bin/python -m pytest tests/ -v
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from spark.types import Project, Site
from spark.indexing import (
    lattice_to_nr,
    nr_to_lattice,
    lattice_to_nr_array,
    total_sites,
    offset_nr,
)


# ---------- helpers ----------

def make_single_layer_project():
    pt = Project()
    pt.add_species(name='empty')
    pt.add_species(name='CO')
    layer = pt.add_layer(name='surface')
    layer.sites.append(Site(name='hollow'))
    layer.sites.append(Site(name='bridge'))
    return pt


def make_multi_layer_project():
    pt = Project()
    pt.add_species(name='empty')
    pt.add_species(name='CO')
    pt.add_species(name='O')
    pt.add_species(name='Pd_atom')

    # Layer Pd100: 5 sites
    pd100 = pt.add_layer(name='Pd100')
    for n in ['h1', 'h2', 'h3', 'b1', 'b2']:
        pd100.sites.append(Site(name=n))

    # Layer PdO: 3 sites
    pdo = pt.add_layer(name='PdO')
    for n in ['h1', 'b1', 'Pd1']:
        pdo.sites.append(Site(name=n))

    return pt


# ---------- 1. Round-trip bijection ----------

def test_roundtrip_single_layer_3x3x1():
    pt = make_single_layer_project()
    spuck = pt.lattice.spuck
    assert spuck == 2
    L = (3, 3, 1)
    N = total_sites(L, spuck)
    assert N == 18

    for nr in range(N):
        cx, cy, cz, s = nr_to_lattice(nr, L, spuck)
        nr2 = lattice_to_nr(cx, cy, cz, s, L, spuck)
        assert nr == nr2, f"single-layer roundtrip broke at nr={nr}"


def test_roundtrip_multi_layer_4x3x1():
    pt = make_multi_layer_project()
    spuck = pt.lattice.spuck
    assert spuck == 8
    L = (4, 3, 1)
    N = total_sites(L, spuck)
    assert N == 96

    for nr in range(N):
        cx, cy, cz, s = nr_to_lattice(nr, L, spuck)
        nr2 = lattice_to_nr(cx, cy, cz, s, L, spuck)
        assert nr == nr2, f"multi-layer roundtrip broke at nr={nr}"


def test_roundtrip_3d_2x2x2_spuck5():
    """3D super-cell, spuck=5."""
    L = (2, 2, 2)
    spuck = 5
    N = total_sites(L, spuck)
    assert N == 40
    for nr in range(N):
        cx, cy, cz, s = nr_to_lattice(nr, L, spuck)
        assert lattice_to_nr(cx, cy, cz, s, L, spuck) == nr


# ---------- 2. site name resolution ----------

def test_layer_offset_and_site_resolution():
    pt = make_multi_layer_project()
    L_obj = pt.lattice

    # Layer offsets
    assert L_obj.layer_offset('Pd100') == 0
    assert L_obj.layer_offset('PdO') == 5

    # Within Pd100 layer
    assert L_obj.site_in_cell_id('Pd100', 'h1') == 0
    assert L_obj.site_in_cell_id('Pd100', 'h2') == 1
    assert L_obj.site_in_cell_id('Pd100', 'b2') == 4

    # Within PdO layer (offset by 5)
    assert L_obj.site_in_cell_id('PdO', 'h1') == 5
    assert L_obj.site_in_cell_id('PdO', 'b1') == 6
    assert L_obj.site_in_cell_id('PdO', 'Pd1') == 7

    # Inverse mapping site_in_cell_to_layer
    assert L_obj.site_in_cell_to_layer(0) == ('Pd100', 0)
    assert L_obj.site_in_cell_to_layer(4) == ('Pd100', 4)
    assert L_obj.site_in_cell_to_layer(5) == ('PdO', 0)
    assert L_obj.site_in_cell_to_layer(7) == ('PdO', 2)

    # Same site name in different layers resolves separately
    assert L_obj.site_in_cell_id('Pd100', 'h1') != L_obj.site_in_cell_id('PdO', 'h1')


def test_default_and_substrate_layer():
    pt = make_multi_layer_project()
    # First-added layer becomes default + substrate
    assert pt.lattice.default_layer == 'Pd100'
    assert pt.lattice.substrate_layer == 'Pd100'


def test_error_paths():
    pt = make_multi_layer_project()
    L_obj = pt.lattice

    try:
        L_obj.layer_offset('not_a_layer')
        assert False, "should have raised"
    except KeyError:
        pass

    try:
        L_obj.site_in_cell_id('Pd100', 'not_a_site')
        assert False, "should have raised"
    except KeyError:
        pass

    try:
        L_obj.site_in_cell_to_layer(99)
        assert False, "should have raised"
    except IndexError:
        pass


# ---------- 3. offset_nr correctness ----------

def test_offset_nr_within_layer():
    """Within-layer neighbor: same layer, +1 cell in x."""
    L = (4, 4, 1)
    spuck = 8

    # Source: cell (1,1,0), Pd100.h1 (s=0)
    nr = lattice_to_nr(1, 1, 0, 0, L, spuck)

    # Target: cell (2,1,0), Pd100.h1 (s=0)
    nr_target = offset_nr(nr, dcx=1, dcy=0, dcz=0, ds=0,
                          system_size=L, spuck=spuck)

    cx, cy, cz, s = nr_to_lattice(nr_target, L, spuck)
    assert (cx, cy, cz, s) == (2, 1, 0, 0)


def test_offset_nr_cross_layer():
    """Cross-layer process: Pd100.b2 -> PdO.Pd1 same cell.

    This is the lattice-swap operation pattern from the paper.
    """
    pt = make_multi_layer_project()
    L_obj = pt.lattice

    L = (3, 3, 1)
    spuck = L_obj.spuck

    s_src = L_obj.site_in_cell_id('Pd100', 'b2')   # = 4
    s_tgt = L_obj.site_in_cell_id('PdO', 'Pd1')    # = 7

    nr_src = lattice_to_nr(2, 1, 0, s_src, L, spuck)
    nr_tgt = offset_nr(nr_src, 0, 0, 0, ds=s_tgt, system_size=L, spuck=spuck)

    cx, cy, cz, s = nr_to_lattice(nr_tgt, L, spuck)
    # Same cell, different layer
    assert (cx, cy, cz) == (2, 1, 0)
    assert s == 7
    # Confirm layer id
    assert L_obj.site_in_cell_to_layer(s) == ('PdO', 2)


def test_offset_nr_pbc_wrap():
    """Cell offset across PBC boundary."""
    L = (4, 4, 1)
    spuck = 3

    # Source: cell (3,0,0), site 0
    nr = lattice_to_nr(3, 0, 0, 0, L, spuck)

    # +1 in x wraps to cell (0,0,0)
    nr_target = offset_nr(nr, dcx=1, dcy=0, dcz=0, ds=0,
                          system_size=L, spuck=spuck)
    cx, cy, cz, _ = nr_to_lattice(nr_target, L, spuck)
    assert (cx, cy, cz) == (0, 0, 0)


# ---------- 4. Vectorized batch lookup ----------

def test_lattice_to_nr_array_matches_scalar():
    L = (3, 3, 2)
    spuck = 4
    rng = np.random.default_rng(42)
    cells = rng.integers(-2, 5, size=(20, 3))   # include negatives -> wrap
    sites = rng.integers(0, spuck, size=(20,))

    nr_batch = lattice_to_nr_array(cells, sites, L, spuck)
    nr_scalar = np.array([
        lattice_to_nr(c[0], c[1], c[2], s, L, spuck)
        for c, s in zip(cells, sites)
    ])
    assert np.array_equal(nr_batch, nr_scalar)


# ---------- 5. Single-layer back-compat ----------

def test_single_layer_natural_embedding():
    """When spuck = nsites_layer0, indexing must match the trivial 1D layout
    that single-layer engine code currently assumes."""
    L = (4, 5, 1)
    spuck = 2

    expected_total = 4 * 5 * 1 * spuck
    assert total_sites(L, spuck) == expected_total

    # nr layout: [(c=0,s=0), (c=0,s=1), (c=1,s=0), (c=1,s=1), ...] where
    # c=cx + Lx*cy + Lx*Ly*cz. Verify exactly.
    for cx in range(4):
        for cy in range(5):
            for s in range(spuck):
                expected = spuck * (cx + 4 * cy) + s
                assert lattice_to_nr(cx, cy, 0, s, L, spuck) == expected


if __name__ == '__main__':
    # Run all test_* functions in this module
    for name, fn in list(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn()
                print(f'  PASS  {name}')
            except AssertionError as e:
                print(f'  FAIL  {name}: {e}')
                raise
    print()
    print('Phase A.2 indexing tests: ALL PASS')
