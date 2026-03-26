#!/usr/bin/env python3
"""
Generate all POSCAR files for NO₃RR on Cu(100) DFT calculations.

Usage:
    conda activate /scratch/reny0b/iops/sw/envs/gs
    python generate_all_structures.py [--a0 3.615] [--output-dir ../calculations]

This script generates:
  Phase 0:  Bulk Cu + 8 gas-phase molecules
  Phase 1:  Clean p(4×4) Cu(100) slab + 13 single adsorbates
  Phase 2:  Tier 1 co-adsorption (4 pairs × 2 distances = 8)
  Phase 3:  Tier 2 co-adsorption (6 pairs × 2 distances = 12)
  Phase 4:  Extra NEB IS/FS co-adsorption structures (~12)

Author: Auto-generated for SPARK Phase 2
Date: 2026-03-24
"""

import os
import argparse
import numpy as np
from pathlib import Path

from ase import Atoms
from ase.build import bulk, fcc100, add_adsorbate, molecule
from ase.constraints import FixAtoms
from ase.io import write


# ============================================================================
# Configuration
# ============================================================================

# Experimental Cu lattice constant (PBE will optimize to ~3.63 Å)
DEFAULT_A0 = 3.615

# Slab parameters
N_LAYERS = 4
VACUUM = 15.0        # Å
FIX_BOTTOM = 2       # Fix bottom 2 layers

# Gas-phase box size (asymmetric to break symmetry)
GAS_BOX = [15.0, 16.0, 17.0]


# ============================================================================
# Helper functions
# ============================================================================

def ensure_dir(path):
    """Create directory if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def write_poscar(atoms, filepath, sort=True):
    """Write POSCAR with sorted elements (VASP convention)."""
    ensure_dir(os.path.dirname(filepath))
    write(filepath, atoms, format='vasp', sort=sort, vasp5=True)
    print(f"  Written: {filepath}  ({len(atoms)} atoms)")


def make_slab(a0, size=(4, 4, 4)):
    """
    Build Cu(100) slab with selective dynamics.

    Returns slab with bottom FIX_BOTTOM layers fixed.
    """
    slab = fcc100('Cu', size=size, a=a0, vacuum=VACUUM, periodic=True)

    # Identify bottom layers by z-coordinate
    positions = slab.get_positions()
    z_coords = positions[:, 2]
    unique_z = np.sort(np.unique(np.round(z_coords, 3)))

    # Fix bottom FIX_BOTTOM layers
    fix_z_max = unique_z[FIX_BOTTOM - 1] + 0.1
    fix_indices = [i for i, z in enumerate(z_coords) if z <= fix_z_max]
    slab.set_constraint(FixAtoms(indices=fix_indices))

    return slab


def get_hollow_site(slab, ix, iy, a_surf):
    """
    Get Cartesian (x, y) of a hollow site on Cu(100).

    Cu(100) hollow = 4-fold site, located between 4 surface Cu atoms.
    ix, iy: integer indices of the hollow site in the p(4×4) grid.
    """
    # On Cu(100), surface Cu at (i*a_surf, j*a_surf).
    # Hollow at ((i+0.5)*a_surf, (j+0.5)*a_surf).
    x = (ix + 0.5) * a_surf
    y = (iy + 0.5) * a_surf
    return x, y


def get_bridge_site(slab, ix, iy, direction, a_surf):
    """
    Get Cartesian (x, y) of a bridge site on Cu(100).

    direction: 'x' (along [100]) or 'y' (along [010])
    """
    if direction == 'x':
        x = (ix + 0.5) * a_surf
        y = iy * a_surf
    else:
        x = ix * a_surf
        y = (iy + 0.5) * a_surf
    return x, y


def get_atop_site(slab, ix, iy, a_surf):
    """Get Cartesian (x, y) of an atop site on Cu(100)."""
    x = ix * a_surf
    y = iy * a_surf
    return x, y


def get_top_z(slab):
    """Get z-coordinate of the topmost Cu layer."""
    return max(slab.get_positions()[:, 2])


def add_atom_at_site(slab, symbol, site_xy, height):
    """Add a single atom at (x, y, z_top + height)."""
    z = get_top_z(slab) + height
    atom = Atoms(symbol, positions=[(site_xy[0], site_xy[1], z)])
    return slab + atom


def add_molecule_at_site(slab, mol_atoms, site_xy, height):
    """
    Add a molecule fragment at a site.
    mol_atoms: Atoms object (molecule fragment)
    The molecule center-of-mass is placed at (site_xy[0], site_xy[1], z_top + height).
    """
    z_top = get_top_z(slab)
    mol = mol_atoms.copy()
    # Center molecule
    com = mol.get_center_of_mass()
    mol.translate(np.array([site_xy[0], site_xy[1], z_top + height]) - com)
    return slab + mol


# ============================================================================
# Phase 0: Bulk Cu + Gas-phase molecules
# ============================================================================

def generate_phase0(a0, base_dir):
    """Generate bulk Cu and gas-phase molecule POSCARs."""
    print("\n=== Phase 0: Bulk Cu + Gas-phase Molecules ===")
    phase_dir = os.path.join(base_dir, "phase0_bulk_gas")

    # --- Gas-phase molecules ---
    gas_molecules = {
        'H2': {
            'atoms': Atoms('H2', positions=[(0, 0, 0), (0, 0, 0.74)]),
            'ispin_note': 'singlet',
        },
        'H2O': {
            'atoms': molecule('H2O'),
            'ispin_note': 'singlet',
        },
        'NO': {
            'atoms': Atoms('NO', positions=[(0, 0, 0), (0, 0, 1.15)]),
            'ispin_note': 'ISPIN=2, doublet',
        },
        'NH3': {
            'atoms': molecule('NH3'),
            'ispin_note': 'singlet',
        },
        'N2': {
            'atoms': Atoms('N2', positions=[(0, 0, 0), (0, 0, 1.10)]),
            'ispin_note': 'singlet',
        },
        'N2O': {
            'atoms': Atoms('N2O', positions=[
                (0, 0, -1.16),   # N (terminal)
                (0, 0, 0),       # N (central)
                (0, 0, 1.19),    # O
            ]),
            'ispin_note': 'singlet, linear',
        },
        'NO2': {
            'atoms': Atoms('NO2', positions=[
                (0, 0, 0),                             # N
                (-1.20 * np.sin(np.radians(67)), 0,
                  1.20 * np.cos(np.radians(67))),      # O1
                ( 1.20 * np.sin(np.radians(67)), 0,
                  1.20 * np.cos(np.radians(67))),      # O2
            ]),
            'ispin_note': 'ISPIN=2, doublet, bent 134°',
        },
        'OH': {
            'atoms': Atoms('OH', positions=[(0, 0, 0), (0, 0, 0.97)]),
            'ispin_note': 'ISPIN=2, doublet',
        },
    }

    for name, info in gas_molecules.items():
        mol = info['atoms'].copy()
        # Place in asymmetric box to avoid spurious symmetry
        mol.set_cell(GAS_BOX)
        mol.center()
        mol.set_pbc(True)
        mol_dir = os.path.join(phase_dir, f"gas_{name}")
        write_poscar(mol, os.path.join(mol_dir, "POSCAR"))

        # Write a note about spin
        with open(os.path.join(mol_dir, "NOTE"), 'w') as f:
            f.write(f"Molecule: {name}\n")
            f.write(f"Spin: {info['ispin_note']}\n")
            f.write(f"Box: {GAS_BOX[0]}x{GAS_BOX[1]}x{GAS_BOX[2]} Å\n")
            f.write(f"K-points: Gamma only (1x1x1)\n")

    print(f"  Phase 0 total: {len(gas_molecules)} gas-phase POSCARs")


# ============================================================================
# Phase 1: Clean slab + single adsorbates
# ============================================================================

def generate_phase1(a0, base_dir):
    """Generate clean slab and 13 single-adsorbate structures."""
    print("\n=== Phase 1: Clean Slab + 13 Single Adsorbates ===")
    phase_dir = os.path.join(base_dir, "phase1_single_ads")

    a_surf = a0 / np.sqrt(2)
    slab = make_slab(a0)
    z_top = get_top_z(slab)

    # --- Clean slab ---
    write_poscar(slab.copy(), os.path.join(phase_dir, "clean_slab", "POSCAR"))

    # --- Adsorbate definitions ---
    # Each: (name, build_function)
    # Heights are initial guesses; VASP optimization will find the true minimum
    adsorbates = {}

    # Hollow site reference: (1, 1) gives a hollow site away from cell edges
    h_xy = get_hollow_site(slab, 1, 1, a_surf)

    # Bridge site reference
    b_xy = get_bridge_site(slab, 1, 1, 'x', a_surf)

    # Atop site reference
    a_xy = get_atop_site(slab, 1, 1, a_surf)

    # *H @ hollow (height ~0.9 Å above top Cu)
    def make_H():
        s = slab.copy()
        return add_atom_at_site(s, 'H', h_xy, 0.9)
    adsorbates['H_hollow'] = make_H

    # *O @ hollow (height ~0.6 Å)
    def make_O():
        s = slab.copy()
        return add_atom_at_site(s, 'O', h_xy, 0.6)
    adsorbates['O_hollow'] = make_O

    # *N @ hollow (height ~0.3 Å, deeply embedded)
    def make_N():
        s = slab.copy()
        return add_atom_at_site(s, 'N', h_xy, 0.3)
    adsorbates['N_hollow'] = make_N

    # *OH @ bridge/hollow (O at hollow, H pointing up-away)
    def make_OH():
        s = slab.copy()
        oh = Atoms('OH', positions=[(0, 0, 0), (0, 0, 0.97)])
        return add_molecule_at_site(s, oh, h_xy, 1.5)
    adsorbates['OH_hollow'] = make_OH

    # *NO @ hollow (N down, height ~1.3 Å for N)
    def make_NO():
        s = slab.copy()
        no = Atoms('NO', positions=[(0, 0, 0), (0, 0, 1.15)])  # N at bottom
        return add_molecule_at_site(s, no, h_xy, 1.4)
    adsorbates['NO_hollow'] = make_NO

    # *NO₂ @ bridge (O atoms pointing down, bidentate)
    def make_NO2():
        s = slab.copy()
        # N on top, two O atoms pointing down towards bridge Cu atoms
        angle = np.radians(60)  # O-N-O half-angle for bidentate
        d_NO = 1.24
        no2 = Atoms('NO2', positions=[
            (0, 0, d_NO * np.cos(angle)),           # N (top)
            (-d_NO * np.sin(angle), 0, 0),           # O1
            ( d_NO * np.sin(angle), 0, 0),           # O2
        ])
        return add_molecule_at_site(s, no2, b_xy, 2.0)
    adsorbates['NO2_bridge'] = make_NO2

    # *NO₃ @ O-down tridentate/bidentate
    def make_NO3():
        s = slab.copy()
        # Planar NO₃, N in center, 3 O atoms 120° apart, O atoms down
        d_NO = 1.25
        no3 = Atoms('NO3', positions=[
            (0, 0, 0.8),                                        # N (top)
            (d_NO, 0, 0),                                       # O1
            (-d_NO * 0.5, d_NO * np.sqrt(3)/2, 0),             # O2
            (-d_NO * 0.5, -d_NO * np.sqrt(3)/2, 0),            # O3
        ])
        return add_molecule_at_site(s, no3, h_xy, 2.0)
    adsorbates['NO3_hollow'] = make_NO3

    # *NOH @ hollow (N down, O-H pointing up)
    def make_NOH():
        s = slab.copy()
        noh = Atoms('NOH', positions=[
            (0, 0, 0),         # N (bottom, on surface)
            (0, 0, 1.20),      # O
            (0, 0.80, 1.65),   # H
        ])
        return add_molecule_at_site(s, noh, h_xy, 1.3)
    adsorbates['NOH_hollow'] = make_NOH

    # *NHOH @ hollow/bridge
    def make_NHOH():
        s = slab.copy()
        nhoh = Atoms('NHOH', positions=[
            (0, 0, 0),           # N (bottom)
            (-0.50, 0, 0.80),    # H on N
            (0, 0, 1.35),        # O
            (0, 0.80, 1.80),     # H on O
        ])
        return add_molecule_at_site(s, nhoh, h_xy, 1.3)
    adsorbates['NHOH_hollow'] = make_NHOH

    # *NH @ hollow (height ~1.3 Å)
    def make_NH():
        s = slab.copy()
        nh = Atoms('NH', positions=[(0, 0, 0), (0, 0, 1.02)])
        return add_molecule_at_site(s, nh, h_xy, 1.0)
    adsorbates['NH_hollow'] = make_NH

    # *NH₂ @ bridge
    def make_NH2():
        s = slab.copy()
        nh2 = Atoms('NH2', positions=[
            (0, 0, 0),
            (-0.47, 0, 0.82),
            (0.47, 0, 0.82),
        ])
        return add_molecule_at_site(s, nh2, b_xy, 1.6)
    adsorbates['NH2_bridge'] = make_NH2

    # *NH₃ @ atop (weak adsorption, N down)
    def make_NH3():
        s = slab.copy()
        nh3 = molecule('NH3')
        # Flip so N points down
        nh3.positions[:, 2] = -nh3.positions[:, 2]
        return add_molecule_at_site(s, nh3, a_xy, 2.2)
    adsorbates['NH3_atop'] = make_NH3

    # *N₂O @ atop/bridge (weak, linear, N-end down)
    def make_N2O():
        s = slab.copy()
        n2o = Atoms('N2O', positions=[
            (0, 0, 0),       # N (on surface)
            (0, 0, 1.13),    # N
            (0, 0, 2.32),    # O (top)
        ])
        return add_molecule_at_site(s, n2o, a_xy, 2.0)
    adsorbates['N2O_atop'] = make_N2O

    # Generate all
    count = 0
    for name, builder in adsorbates.items():
        atoms = builder()
        write_poscar(atoms, os.path.join(phase_dir, name, "POSCAR"))
        count += 1

    print(f"  Phase 1 total: 1 clean + {count} adsorbates = {1+count} POSCARs")


# ============================================================================
# Phase 2-3: Co-adsorption for lateral interactions
# ============================================================================

def generate_coads(a0, base_dir):
    """Generate co-adsorption structures for lateral interaction calculation."""
    print("\n=== Phase 2-3: Co-adsorption (Lateral Interactions) ===")

    a_surf = a0 / np.sqrt(2)
    slab = make_slab(a0)

    # NN distances on Cu(100):
    # 1NN: along [110], d = a0/sqrt(2) ≈ 2.566 Å → hollow sites separated by 1 unit
    # 2NN: along [100], d = a0 ≈ 3.63 Å → hollow sites separated by sqrt(2) units (diagonal)

    # For hollow-site species on p(4×4):
    # 1NN pair: sites at (0,0) and (1,0) → separated by a_surf
    # 2NN pair: sites at (0,0) and (1,1) → separated by a_surf*sqrt(2)=a0

    def make_two_atoms_hollow(sym1, h1, sym2, h2, ix1, iy1, ix2, iy2):
        """Place two atoms at hollow sites."""
        s = slab.copy()
        z_top = get_top_z(s)
        xy1 = get_hollow_site(s, ix1, iy1, a_surf)
        xy2 = get_hollow_site(s, ix2, iy2, a_surf)
        s += Atoms(sym1, positions=[(xy1[0], xy1[1], z_top + h1)])
        s += Atoms(sym2, positions=[(xy2[0], xy2[1], z_top + h2)])
        return s

    def make_NO_at_hollow(slab_in, ix, iy):
        """Add *NO at hollow site (N down, O up)."""
        s = slab_in
        z_top = get_top_z(s)
        xy = get_hollow_site(s, ix, iy, a_surf)
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 1.0)])
        s += Atoms('O', positions=[(xy[0], xy[1], z_top + 2.15)])
        return s

    def make_OH_at_hollow(slab_in, ix, iy):
        """Add *OH at hollow site."""
        s = slab_in
        z_top = get_top_z(s)
        xy = get_hollow_site(s, ix, iy, a_surf)
        s += Atoms('O', positions=[(xy[0], xy[1], z_top + 1.2)])
        s += Atoms('H', positions=[(xy[0], xy[1], z_top + 2.17)])
        return s

    def make_H_at_hollow(slab_in, ix, iy):
        """Add *H at hollow site."""
        s = slab_in
        z_top = get_top_z(s)
        xy = get_hollow_site(s, ix, iy, a_surf)
        s += Atoms('H', positions=[(xy[0], xy[1], z_top + 0.9)])
        return s

    def make_N_at_hollow(slab_in, ix, iy):
        """Add *N at hollow site."""
        s = slab_in
        z_top = get_top_z(s)
        xy = get_hollow_site(s, ix, iy, a_surf)
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 0.3)])
        return s

    def make_O_at_hollow(slab_in, ix, iy):
        """Add *O at hollow site."""
        s = slab_in
        z_top = get_top_z(s)
        xy = get_hollow_site(s, ix, iy, a_surf)
        s += Atoms('O', positions=[(xy[0], xy[1], z_top + 0.6)])
        return s

    # 1NN pair: hollow sites at (0,0) and (1,0) → separated by a_surf along x
    # 2NN pair: hollow sites at (0,0) and (1,1) → separated by a_surf*sqrt(2) diagonally
    site_1nn = [(0, 0), (1, 0)]
    site_2nn = [(0, 0), (1, 1)]

    # --- Tier 1 (Phase 2) ---
    tier1_dir = os.path.join(base_dir, "phase2_tier1_coads")
    tier1_pairs = [
        ("NO_NO", make_NO_at_hollow, make_NO_at_hollow, "R12 NEB IS reuse"),
        ("NO_H",  make_NO_at_hollow, make_H_at_hollow,  "R6 NEB IS reuse"),
        ("H_H",   make_H_at_hollow,  make_H_at_hollow,  "R16 NEB IS reuse"),
        ("OH_NO", make_OH_at_hollow, make_NO_at_hollow,  "R4 NEB FS reuse"),
    ]

    count = 0
    for pair_name, add_A, add_B, note in tier1_pairs:
        for dist_label, sites in [("1NN", site_1nn), ("2NN", site_2nn)]:
            s = slab.copy()
            s = add_A(s, sites[0][0], sites[0][1])
            s = add_B(s, sites[1][0], sites[1][1])
            dirname = f"{pair_name}_{dist_label}"
            filepath = os.path.join(tier1_dir, dirname, "POSCAR")
            write_poscar(s, filepath)
            # Write note
            with open(os.path.join(tier1_dir, dirname, "NOTE"), 'w') as f:
                f.write(f"Co-adsorption: {pair_name} @ {dist_label}\n")
                f.write(f"Note: {note}\n")
                f.write(f"ε = E(co-ads) − E(slab+A) − E(slab+B) + E(clean)\n")
            count += 1

    print(f"  Tier 1: {count} co-adsorption POSCARs")

    # --- Tier 2 (Phase 3) ---
    tier2_dir = os.path.join(base_dir, "phase3_tier2_coads")
    tier2_pairs = [
        ("OH_OH", make_OH_at_hollow, make_OH_at_hollow, "OH poisoning"),
        ("OH_H",  make_OH_at_hollow, make_H_at_hollow,  "R5 NEB IS reuse"),
        ("N_NO",  make_N_at_hollow,  make_NO_at_hollow,  "alternative N-N coupling"),
        ("N_N",   make_N_at_hollow,  make_N_at_hollow,   "R14 NEB IS reuse"),
        ("O_NO",  make_O_at_hollow,  make_NO_at_hollow,  "deoxidation byproduct"),
        ("N_H",   make_N_at_hollow,  make_H_at_hollow,   "N hydrogenation"),
    ]

    count2 = 0
    for pair_name, add_A, add_B, note in tier2_pairs:
        for dist_label, sites in [("1NN", site_1nn), ("2NN", site_2nn)]:
            s = slab.copy()
            s = add_A(s, sites[0][0], sites[0][1])
            s = add_B(s, sites[1][0], sites[1][1])
            dirname = f"{pair_name}_{dist_label}"
            filepath = os.path.join(tier2_dir, dirname, "POSCAR")
            write_poscar(s, filepath)
            with open(os.path.join(tier2_dir, dirname, "NOTE"), 'w') as f:
                f.write(f"Co-adsorption: {pair_name} @ {dist_label}\n")
                f.write(f"Note: {note}\n")
            count2 += 1

    print(f"  Tier 2: {count2} co-adsorption POSCARs")


# ============================================================================
# Phase 4: Extra NEB IS/FS structures
# ============================================================================

def generate_neb_extra(a0, base_dir):
    """Generate additional co-adsorption structures needed as NEB IS/FS."""
    print("\n=== Phase 4: Extra NEB IS/FS Structures ===")
    phase_dir = os.path.join(base_dir, "phase4_neb_isfs")

    a_surf = a0 / np.sqrt(2)
    slab = make_slab(a0)

    # 1NN sites for IS structures (*A + *H @ adjacent hollows)
    site_A = (0, 0)
    site_H = (1, 0)

    count = 0

    # Helper to add *H at hollow
    def add_H(s, ix, iy):
        z_top = get_top_z(s)
        xy = get_hollow_site(s, ix, iy, a_surf)
        return s + Atoms('H', positions=[(xy[0], xy[1], z_top + 0.9)])

    # --- IS structures: *X + *H @ 1NN ---

    # R2 IS: *NO₃ + *H
    def make_NO3_H_IS():
        s = slab.copy()
        z_top = get_top_z(s)
        xy = get_hollow_site(s, site_A[0], site_A[1], a_surf)
        d_NO = 1.25
        # NO₃ at hollow, O-down
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 2.0)])
        s += Atoms('O', positions=[(xy[0]+d_NO, xy[1], z_top + 1.2)])
        s += Atoms('O', positions=[(xy[0]-d_NO*0.5, xy[1]+d_NO*0.866, z_top + 1.2)])
        s += Atoms('O', positions=[(xy[0]-d_NO*0.5, xy[1]-d_NO*0.866, z_top + 1.2)])
        s = add_H(s, site_H[0], site_H[1])
        return s

    write_poscar(make_NO3_H_IS(),
                 os.path.join(phase_dir, "R2_IS_NO3_H", "POSCAR"))
    count += 1

    # R4 IS: *NO₂ + *H
    def make_NO2_H_IS():
        s = slab.copy()
        z_top = get_top_z(s)
        xy = get_bridge_site(s, site_A[0], site_A[1], 'x', a_surf)
        d_NO = 1.24
        angle = np.radians(60)
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 2.0)])
        s += Atoms('O', positions=[(xy[0]-d_NO*np.sin(angle), xy[1], z_top + 1.2)])
        s += Atoms('O', positions=[(xy[0]+d_NO*np.sin(angle), xy[1], z_top + 1.2)])
        s = add_H(s, site_H[0], site_H[1])
        return s

    write_poscar(make_NO2_H_IS(),
                 os.path.join(phase_dir, "R4_IS_NO2_H", "POSCAR"))
    count += 1

    # R7 IS: *NOH + *H
    def make_NOH_H_IS():
        s = slab.copy()
        z_top = get_top_z(s)
        xy = get_hollow_site(s, site_A[0], site_A[1], a_surf)
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 1.3)])
        s += Atoms('O', positions=[(xy[0], xy[1], z_top + 2.5)])
        s += Atoms('H', positions=[(xy[0], xy[1]+0.8, z_top + 2.95)])
        s = add_H(s, site_H[0], site_H[1])
        return s

    write_poscar(make_NOH_H_IS(),
                 os.path.join(phase_dir, "R7_IS_NOH_H", "POSCAR"))
    count += 1

    # R8 IS: *NHOH + *H
    def make_NHOH_H_IS():
        s = slab.copy()
        z_top = get_top_z(s)
        xy = get_hollow_site(s, site_A[0], site_A[1], a_surf)
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 1.3)])
        s += Atoms('H', positions=[(xy[0]-0.5, xy[1], z_top + 2.1)])
        s += Atoms('O', positions=[(xy[0], xy[1], z_top + 2.65)])
        s += Atoms('H', positions=[(xy[0], xy[1]+0.8, z_top + 3.1)])
        s = add_H(s, site_H[0], site_H[1])
        return s

    write_poscar(make_NHOH_H_IS(),
                 os.path.join(phase_dir, "R8_IS_NHOH_H", "POSCAR"))
    count += 1

    # R9 IS: *NH + *H
    def make_NH_H_IS():
        s = slab.copy()
        z_top = get_top_z(s)
        xy = get_hollow_site(s, site_A[0], site_A[1], a_surf)
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 1.0)])
        s += Atoms('H', positions=[(xy[0], xy[1], z_top + 2.02)])
        s = add_H(s, site_H[0], site_H[1])
        return s

    write_poscar(make_NH_H_IS(),
                 os.path.join(phase_dir, "R9_IS_NH_H", "POSCAR"))
    count += 1

    # R10 IS: *NH₂ + *H
    def make_NH2_H_IS():
        s = slab.copy()
        z_top = get_top_z(s)
        xy_b = get_bridge_site(s, site_A[0], site_A[1], 'x', a_surf)
        s += Atoms('N', positions=[(xy_b[0], xy_b[1], z_top + 1.6)])
        s += Atoms('H', positions=[(xy_b[0]-0.47, xy_b[1], z_top + 2.42)])
        s += Atoms('H', positions=[(xy_b[0]+0.47, xy_b[1], z_top + 2.42)])
        s = add_H(s, site_H[0], site_H[1])
        return s

    write_poscar(make_NH2_H_IS(),
                 os.path.join(phase_dir, "R10_IS_NH2_H", "POSCAR"))
    count += 1

    # R15 IS: H₂O on slab (Volmer IS)
    def make_H2O_slab():
        s = slab.copy()
        z_top = get_top_z(s)
        xy = get_atop_site(s, 1, 1, a_surf)
        # H₂O with O pointing down
        s += Atoms('O', positions=[(xy[0], xy[1], z_top + 2.5)])
        s += Atoms('H', positions=[(xy[0]-0.76, xy[1], z_top + 3.1)])
        s += Atoms('H', positions=[(xy[0]+0.76, xy[1], z_top + 3.1)])
        return s

    write_poscar(make_H2O_slab(),
                 os.path.join(phase_dir, "R15_IS_H2O_slab", "POSCAR"))
    count += 1

    # --- FS structures that need new calculations ---

    # R2 FS: *NO₂ + *OH @ 1NN
    def make_NO2_OH_FS():
        s = slab.copy()
        z_top = get_top_z(s)
        # NO₂ at bridge
        xy_b = get_bridge_site(s, site_A[0], site_A[1], 'x', a_surf)
        d_NO = 1.24
        angle = np.radians(60)
        s += Atoms('N', positions=[(xy_b[0], xy_b[1], z_top + 2.0)])
        s += Atoms('O', positions=[(xy_b[0]-d_NO*np.sin(angle), xy_b[1], z_top + 1.2)])
        s += Atoms('O', positions=[(xy_b[0]+d_NO*np.sin(angle), xy_b[1], z_top + 1.2)])
        # OH at adjacent hollow
        xy_oh = get_hollow_site(s, site_H[0], site_H[1], a_surf)
        s += Atoms('O', positions=[(xy_oh[0], xy_oh[1], z_top + 1.2)])
        s += Atoms('H', positions=[(xy_oh[0], xy_oh[1], z_top + 2.17)])
        return s

    write_poscar(make_NO2_OH_FS(),
                 os.path.join(phase_dir, "R2_FS_NO2_OH", "POSCAR"))
    count += 1

    # R8 FS: *NH + H₂O above
    def make_NH_H2O_FS():
        s = slab.copy()
        z_top = get_top_z(s)
        xy = get_hollow_site(s, site_A[0], site_A[1], a_surf)
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 1.0)])
        s += Atoms('H', positions=[(xy[0], xy[1], z_top + 2.02)])
        # H₂O desorbing above
        xy2 = get_hollow_site(s, site_H[0], site_H[1], a_surf)
        s += Atoms('O', positions=[(xy2[0], xy2[1], z_top + 3.5)])
        s += Atoms('H', positions=[(xy2[0]-0.76, xy2[1], z_top + 4.1)])
        s += Atoms('H', positions=[(xy2[0]+0.76, xy2[1], z_top + 4.1)])
        return s

    write_poscar(make_NH_H2O_FS(),
                 os.path.join(phase_dir, "R8_FS_NH_H2O", "POSCAR"))
    count += 1

    # R14 FS: N₂ above + 2 empty sites
    def make_N2_desorb_FS():
        s = slab.copy()
        z_top = get_top_z(s)
        # N₂ molecule 4 Å above surface (desorbing)
        xy = get_hollow_site(s, 0, 0, a_surf)
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 4.0)])
        s += Atoms('N', positions=[(xy[0], xy[1], z_top + 5.1)])
        return s

    write_poscar(make_N2_desorb_FS(),
                 os.path.join(phase_dir, "R14_FS_N2_desorb", "POSCAR"))
    count += 1

    # R15 FS: *H + OH⁻ above (Volmer FS)
    def make_H_OHminus_FS():
        s = slab.copy()
        z_top = get_top_z(s)
        # *H at hollow
        xy_h = get_hollow_site(s, site_A[0], site_A[1], a_surf)
        s += Atoms('H', positions=[(xy_h[0], xy_h[1], z_top + 0.9)])
        # OH⁻ desorbing above
        xy2 = get_hollow_site(s, site_H[0], site_H[1], a_surf)
        s += Atoms('O', positions=[(xy2[0], xy2[1], z_top + 4.0)])
        s += Atoms('H', positions=[(xy2[0], xy2[1], z_top + 4.97)])
        return s

    write_poscar(make_H_OHminus_FS(),
                 os.path.join(phase_dir, "R15_FS_H_OH", "POSCAR"))
    count += 1

    # R16 FS: H₂ above + 2 empty sites (Tafel FS)
    def make_H2_desorb_FS():
        s = slab.copy()
        z_top = get_top_z(s)
        xy = get_hollow_site(s, 0, 0, a_surf)
        s += Atoms('H', positions=[(xy[0], xy[1], z_top + 4.0)])
        s += Atoms('H', positions=[(xy[0], xy[1], z_top + 4.74)])
        return s

    write_poscar(make_H2_desorb_FS(),
                 os.path.join(phase_dir, "R16_FS_H2_desorb", "POSCAR"))
    count += 1

    print(f"  Phase 4 total: {count} NEB IS/FS POSCARs")


# ============================================================================
# Phase TS: NEB diffusion IS/FS pairs
# ============================================================================

def generate_diffusion(a0, base_dir):
    """Generate diffusion NEB endpoint structures."""
    print("\n=== Diffusion NEB Endpoints ===")
    phase_dir = os.path.join(base_dir, "neb_diffusion")

    a_surf = a0 / np.sqrt(2)
    slab = make_slab(a0)

    # *NO diffusion: hollow → hollow (via bridge)
    for label, species_sym, height in [("NO", "N", 1.0), ("H", "H", 0.9)]:
        # IS: at hollow (0,0)
        s_is = slab.copy()
        z_top = get_top_z(s_is)
        xy_is = get_hollow_site(s_is, 1, 1, a_surf)
        if label == "NO":
            s_is += Atoms('N', positions=[(xy_is[0], xy_is[1], z_top + 1.0)])
            s_is += Atoms('O', positions=[(xy_is[0], xy_is[1], z_top + 2.15)])
        else:
            s_is += Atoms('H', positions=[(xy_is[0], xy_is[1], z_top + height)])

        write_poscar(s_is, os.path.join(phase_dir, f"R_diff_{label}", "IS", "POSCAR"))

        # FS: at adjacent hollow (2,1)
        s_fs = slab.copy()
        z_top = get_top_z(s_fs)
        xy_fs = get_hollow_site(s_fs, 2, 1, a_surf)
        if label == "NO":
            s_fs += Atoms('N', positions=[(xy_fs[0], xy_fs[1], z_top + 1.0)])
            s_fs += Atoms('O', positions=[(xy_fs[0], xy_fs[1], z_top + 2.15)])
        else:
            s_fs += Atoms('H', positions=[(xy_fs[0], xy_fs[1], z_top + height)])

        write_poscar(s_fs, os.path.join(phase_dir, f"R_diff_{label}", "FS", "POSCAR"))

    print(f"  Diffusion: 2 reactions × 2 endpoints = 4 POSCARs")


# ============================================================================
# POTCAR generation script
# ============================================================================

def generate_potcar_script(base_dir):
    """Generate a shell script that creates POTCAR files for all directories."""
    print("\n=== Generating POTCAR setup script ===")

    script_content = r"""#!/bin/bash
# =============================================================================
# setup_potcar.sh — Generate POTCAR for all POSCAR directories
#
# Usage:
#     cd <calculations_root>
#     bash setup_potcar.sh
#
# Requirements:
#     VASP_PP_PATH must be set (default: /scratch/reny0b/VASP/pot64)
# =============================================================================

VASP_PP_PATH="${VASP_PP_PATH:-/scratch/reny0b/VASP/pot64}"
PP_DIR="${VASP_PP_PATH}/potpaw_PBE"

echo "POTCAR library: ${PP_DIR}"
echo ""

# Element POTCAR paths (standard PAW_PBE, matching DFT_methodology.md)
declare -A POTCAR_FILES
POTCAR_FILES[Cu]="${PP_DIR}/Cu/POTCAR"
POTCAR_FILES[N]="${PP_DIR}/N/POTCAR"
POTCAR_FILES[O]="${PP_DIR}/O/POTCAR"
POTCAR_FILES[H]="${PP_DIR}/H/POTCAR"

# Verify all POTCARs exist
for elem in Cu N O H; do
    if [[ ! -f "${POTCAR_FILES[$elem]}" ]]; then
        echo "ERROR: POTCAR not found: ${POTCAR_FILES[$elem]}"
        exit 1
    fi
done
echo "All element POTCARs verified."
echo ""

# Function: read POSCAR, determine element order, concatenate POTCAR
generate_potcar() {
    local dir="$1"
    local poscar="${dir}/POSCAR"

    if [[ ! -f "$poscar" ]]; then
        return
    fi

    # Read element line (line 6 in VASP5 format)
    local elements
    elements=$(sed -n '6p' "$poscar")

    # Build POTCAR by concatenating in order
    local potcar="${dir}/POTCAR"
    > "$potcar"

    for elem in $elements; do
        if [[ -n "${POTCAR_FILES[$elem]}" ]]; then
            cat "${POTCAR_FILES[$elem]}" >> "$potcar"
        else
            echo "  WARNING: Unknown element '$elem' in $poscar"
            return
        fi
    done

    echo "  POTCAR: ${dir} (${elements})"
}

# Find all directories containing POSCAR and generate POTCAR
count=0
while IFS= read -r -d '' poscar_file; do
    dir=$(dirname "$poscar_file")
    generate_potcar "$dir"
    ((count++))
done < <(find . -name "POSCAR" -print0)

echo ""
echo "Done. Generated ${count} POTCAR files."
"""

    script_path = os.path.join(base_dir, "setup_potcar.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    print(f"  Written: {script_path}")


# ============================================================================
# Summary
# ============================================================================

def print_summary(base_dir):
    """Print summary of generated structures."""
    print("\n" + "=" * 70)
    print("SUMMARY OF GENERATED STRUCTURES")
    print("=" * 70)

    total = 0
    for root, dirs, files in os.walk(base_dir):
        if "POSCAR" in files:
            total += 1

    print(f"""
Directory: {base_dir}

  Phase 0: Bulk Cu + 8 gas-phase molecules ........ 9
  Phase 1: Clean slab + 13 adsorbates ............. 14
  Phase 2: Tier 1 co-adsorption (4×2) ............. 8
  Phase 3: Tier 2 co-adsorption (6×2) ............. 12
  Phase 4: Extra NEB IS/FS ....................... ~12
  Diffusion: NEB endpoints ........................ 4
  ─────────────────────────────────────────────────
  Total POSCAR files found: {total}

Next steps:
  1. Run bulk Cu optimization → get a₀(PBE)
  2. If a₀ differs significantly from {DEFAULT_A0} Å, re-run this script
     with --a0 <optimized_value> to regenerate slab structures
  3. Run: bash setup_potcar.sh   (generates all POTCAR files)
  4. Copy INCAR/KPOINTS from /scratch/reny0b/zls/script/vasp/
  5. Submit jobs with cp-vasp.sh
""")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate all POSCAR files for NO₃RR on Cu(100)")
    parser.add_argument('--a0', type=float, default=DEFAULT_A0,
                        help=f'Cu lattice constant in Å (default: {DEFAULT_A0})')
    parser.add_argument('--output-dir', type=str, default='../calculations',
                        help='Output directory (default: ../calculations)')
    args = parser.parse_args()

    base_dir = os.path.abspath(args.output_dir)
    print(f"Cu lattice constant: {args.a0} Å")
    print(f"Output directory: {base_dir}")

    generate_phase0(args.a0, base_dir)
    generate_phase1(args.a0, base_dir)
    generate_coads(args.a0, base_dir)
    generate_neb_extra(args.a0, base_dir)
    generate_diffusion(args.a0, base_dir)
    generate_potcar_script(base_dir)
    print_summary(base_dir)


if __name__ == '__main__':
    main()
