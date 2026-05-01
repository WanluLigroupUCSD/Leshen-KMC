"""
SPARK multi-lattice site-index utilities.

Implements the bijection between a global 1D site index ``nr`` and the tuple
``(cell_x, cell_y, cell_z, site_within_cell)`` for a (multi-) lattice with
``spuck`` (Sites Per Unit Cell, summed across all layers) and a 3D periodic
super-cell of size ``(Lx, Ly, Lz)``.

The encoding follows Hoffmann-Reuter-Scheffler 2015 multi-lattice kMC: layer
information is folded into ``site_within_cell`` ∈ [0, spuck), so the engine
sees a single flat array of length ``Lx*Ly*Lz*spuck``. The mapping reduces
to ordinary single-lattice kMC when only one layer is registered.

All functions are pure (no side effects) and operate on plain ints / tuples
/ numpy arrays. They are designed to be hot-path safe; in particular,
``lattice_to_nr_array`` vectorizes over arrays for batched lookups.

Indexing convention (matches kmcos lattice.mpy and is documented in
docs/multi_lattice_design.md §1):

    nr = spuck * (cx + Lx*cy + Lx*Ly*cz) + s          # forward
    cz = (nr // spuck) // (Lx*Ly)                      # inverse
    cy = ((nr // spuck) % (Lx*Ly)) // Lx
    cx = (nr // spuck) % Lx
    s  = nr % spuck

Site-within-cell ``s`` decomposes into (layer, site_in_layer) by table
lookup on ``Lattice.site_in_cell_to_layer`` (see spark/types.py).
"""

import numpy as np


def lattice_to_nr(cx, cy, cz, s, system_size, spuck):
    """Forward map: (cell_xyz, site_in_cell) -> 1D index.

    Wraps cell coords with periodic boundaries.

    Parameters
    ----------
    cx, cy, cz : int
        Cell indices (will be wrapped modulo system_size).
    s : int
        Site-within-cell index in [0, spuck).
    system_size : tuple of int
        (Lx, Ly, Lz) super-cell dimensions in unit-cells.
    spuck : int
        Sites per unit cell (sum over all layers).

    Returns
    -------
    nr : int
        Global flat site index in [0, Lx*Ly*Lz*spuck).
    """
    Lx, Ly, Lz = system_size
    cx %= Lx
    cy %= Ly
    cz %= Lz
    return spuck * (cx + Lx * cy + Lx * Ly * cz) + s


def nr_to_lattice(nr, system_size, spuck):
    """Inverse map: 1D index -> (cell_xyz, site_in_cell).

    Parameters
    ----------
    nr : int
        Global flat site index in [0, Lx*Ly*Lz*spuck).
    system_size : tuple of int
        (Lx, Ly, Lz) super-cell dimensions in unit-cells.
    spuck : int
        Sites per unit cell (sum over all layers).

    Returns
    -------
    cx, cy, cz, s : int
    """
    Lx, Ly, _ = system_size
    cell_block = nr // spuck
    s = nr % spuck
    cx = cell_block % Lx
    cy = (cell_block // Lx) % Ly
    cz = cell_block // (Lx * Ly)
    return cx, cy, cz, s


def lattice_to_nr_array(cells, sites, system_size, spuck):
    """Vectorized forward map for batched lookups.

    Parameters
    ----------
    cells : (N, 3) int array
        Cell coordinates per site (will be wrapped).
    sites : (N,) int array
        Site-within-cell index per site.
    system_size, spuck : as in ``lattice_to_nr``.

    Returns
    -------
    nrs : (N,) int array
    """
    Lx, Ly, Lz = system_size
    cells = np.asarray(cells)
    sites = np.asarray(sites)
    cx = cells[:, 0] % Lx
    cy = cells[:, 1] % Ly
    cz = cells[:, 2] % Lz
    return spuck * (cx + Lx * cy + Lx * Ly * cz) + sites


def total_sites(system_size, spuck):
    """Convenience: total number of global sites = Lx * Ly * Lz * spuck."""
    Lx, Ly, Lz = system_size
    return Lx * Ly * Lz * spuck


def offset_nr(nr, dcx, dcy, dcz, ds, system_size, spuck):
    """Apply a (cell-offset, site-offset) shift to ``nr``, with PBC.

    Used by the engine to resolve a process condition/action coord whose
    ``Coord.offset`` is relative to the firing cell. ``ds`` is the absolute
    site-within-cell of the *target* (already resolved from layer + site_name
    via ``Lattice.site_in_cell_id``), NOT a delta — passing ``s + ds`` would
    silently wrap into a different layer when ``s + ds >= spuck``.

    Parameters
    ----------
    nr : int
        Source global index.
    dcx, dcy, dcz : int
        Cell-offset to apply (added to source's cell coords, then wrapped).
    ds : int
        Absolute target site-within-cell ∈ [0, spuck). Replaces the source's
        site-within-cell entirely.
    system_size, spuck : as above.

    Returns
    -------
    nr_target : int
    """
    cx, cy, cz, _ = nr_to_lattice(nr, system_size, spuck)
    return lattice_to_nr(cx + dcx, cy + dcy, cz + dcz, ds, system_size, spuck)
