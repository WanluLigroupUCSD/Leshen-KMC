"""
Dynamic surface representation — graph/weak-lattice with mutable site identity.

Each site (node) carries:
  - atom_type: which element occupies this site (mutable — changes on segregation)
  - site_type: geometric site class derived from local topology (mutable)
  - adsorbate: which species is adsorbed here (mutable — changes on adsorption/desorption)
  - frozen: whether this site participates in structural events

The neighbor graph is stored as an adjacency list and can be updated when
topology changes (reconstruction events).
"""

import numpy as np
from collections import defaultdict


# ── Predefined surface presets ────────────────────────────────────────────

METAL_PRESETS = {
    'Pd': {'mass': 106.42, 'z': 46, 'color': '#1565C0'},
    'Au': {'mass': 196.97, 'z': 79, 'color': '#FFD700'},
    'Pt': {'mass': 195.08, 'z': 78, 'color': '#B0B0B0'},
    'Cu': {'mass': 63.546, 'z': 29, 'color': '#B87333'},
    'Ag': {'mass': 107.87, 'z': 47, 'color': '#C0C0C0'},
    'Ni': {'mass': 58.693, 'z': 28, 'color': '#727272'},
    'Ir': {'mass': 192.22, 'z': 77, 'color': '#2E7D32'},
    'Ru': {'mass': 101.07, 'z': 44, 'color': '#4E342E'},
    'Rh': {'mass': 102.91, 'z': 45, 'color': '#546E7A'},
    'Co': {'mass': 58.933, 'z': 27, 'color': '#3949AB'},
}

ADSORBATE_EMPTY = 0


class SiteNode:
    """
    A single surface site with mutable identity.

    Attributes
    ----------
    index : int
        Unique site index in the surface graph.
    atom_type : str
        Element occupying this site (e.g., 'Pd', 'Au'). Mutable.
    site_type : int
        Geometric site class (0=default). Derived from local environment.
    adsorbate : int
        Adsorbed species ID (0=empty). Mutable.
    position : ndarray, shape (2,) or (3,)
        Fractional or Cartesian position on the surface.
    layer : int
        Layer index (0=surface, 1=subsurface, ...).
    frozen : bool
        If True, this site cannot participate in structural events.
    """

    __slots__ = ('index', 'atom_type', 'site_type', 'adsorbate',
                 'position', 'layer', 'frozen', '_env_hash')

    def __init__(self, index, atom_type='Pd', site_type=0,
                 adsorbate=ADSORBATE_EMPTY, position=None,
                 layer=0, frozen=False):
        self.index = index
        self.atom_type = atom_type
        self.site_type = site_type
        self.adsorbate = adsorbate
        self.position = np.zeros(2) if position is None else np.asarray(position)
        self.layer = layer
        self.frozen = frozen
        self._env_hash = None  # cached, invalidated on neighbor change

    @property
    def is_empty(self):
        return self.adsorbate == ADSORBATE_EMPTY

    @property
    def is_occupied(self):
        return self.adsorbate != ADSORBATE_EMPTY

    def invalidate_cache(self):
        """Mark cached environment hash as stale."""
        self._env_hash = None

    def __repr__(self):
        ads = f"+sp{self.adsorbate}" if self.is_occupied else ""
        return f"Site({self.index}, {self.atom_type}{ads}, L{self.layer})"


class DynamicSurface:
    """
    Graph-based dynamic surface with mutable site identities.

    The surface is represented as an undirected graph where nodes are
    SiteNodes and edges represent nearest-neighbor connectivity.
    Both node properties (atom_type, adsorbate) and graph topology
    can evolve during simulation.

    Parameters
    ----------
    n_sites : int
        Total number of sites.
    """

    def __init__(self, n_sites=0):
        self.sites = []         # list[SiteNode]
        self.neighbors = []     # list[list[int]] — adjacency list
        self.species_names = {ADSORBATE_EMPTY: 'empty'}
        self._species_counter = 1
        self._atom_types = set()

        for i in range(n_sites):
            self.sites.append(SiteNode(index=i))
            self.neighbors.append([])

    # ── Properties ────────────────────────────────────────────────────

    @property
    def n_sites(self):
        return len(self.sites)

    @property
    def atom_types(self):
        """Set of all atom types currently on the surface."""
        return set(s.atom_type for s in self.sites)

    def register_species(self, name):
        """Register an adsorbate species and return its ID."""
        for sid, sname in self.species_names.items():
            if sname == name:
                return sid
        sid = self._species_counter
        self.species_names[sid] = name
        self._species_counter += 1
        return sid

    def species_id(self, name):
        """Get species ID by name."""
        for sid, sname in self.species_names.items():
            if sname == name:
                return sid
        raise KeyError(f"Unknown species: {name}")

    def species_name(self, sid):
        """Get species name by ID."""
        return self.species_names.get(sid, f"sp{sid}")

    # ── Graph construction ────────────────────────────────────────────

    def add_site(self, atom_type='Pd', position=None, layer=0, frozen=False):
        """Add a site and return its index."""
        idx = len(self.sites)
        site = SiteNode(index=idx, atom_type=atom_type,
                        position=position, layer=layer, frozen=frozen)
        self.sites.append(site)
        self.neighbors.append([])
        self._atom_types.add(atom_type)
        return idx

    def add_edge(self, i, j):
        """Add an undirected edge (neighbor connection) between sites i and j."""
        if j not in self.neighbors[i]:
            self.neighbors[i].append(j)
        if i not in self.neighbors[j]:
            self.neighbors[j].append(i)

    def remove_edge(self, i, j):
        """Remove edge between sites i and j."""
        if j in self.neighbors[i]:
            self.neighbors[i].remove(j)
        if i in self.neighbors[j]:
            self.neighbors[j].remove(i)

    def get_nn(self, site_idx):
        """Get nearest-neighbor site indices."""
        return self.neighbors[site_idx]

    def get_nn_sites(self, site_idx):
        """Get nearest-neighbor SiteNode objects."""
        return [self.sites[j] for j in self.neighbors[site_idx]]

    # ── Local environment queries ─────────────────────────────────────

    def nn_composition(self, site_idx):
        """
        Get nearest-neighbor composition as a dict {atom_type: count}.

        Example: {'Pd': 4, 'Au': 2} for a site with 4 Pd and 2 Au neighbors.
        """
        comp = defaultdict(int)
        for j in self.neighbors[site_idx]:
            comp[self.sites[j].atom_type] += 1
        return dict(comp)

    def nn_adsorbate_count(self, site_idx):
        """Count occupied nearest neighbors."""
        return sum(1 for j in self.neighbors[site_idx]
                   if self.sites[j].is_occupied)

    def coordination_number(self, site_idx):
        """Number of nearest neighbors."""
        return len(self.neighbors[site_idx])

    # ── Site mutations ────────────────────────────────────────────────

    def set_adsorbate(self, site_idx, species_id):
        """Adsorb a species on a site."""
        self.sites[site_idx].adsorbate = species_id
        self._invalidate_neighborhood(site_idx)

    def clear_adsorbate(self, site_idx):
        """Remove adsorbate from a site."""
        self.sites[site_idx].adsorbate = ADSORBATE_EMPTY
        self._invalidate_neighborhood(site_idx)

    def swap_atoms(self, i, j):
        """
        Swap atom types between sites i and j (segregation event).
        Adsorbates stay on their original sites.
        """
        self.sites[i].atom_type, self.sites[j].atom_type = \
            self.sites[j].atom_type, self.sites[i].atom_type
        self._invalidate_neighborhood(i)
        self._invalidate_neighborhood(j)

    def convert_site(self, site_idx, new_atom_type):
        """
        Convert a site to a new atom type (site-conversion event).
        """
        self.sites[site_idx].atom_type = new_atom_type
        self._invalidate_neighborhood(site_idx)

    def _invalidate_neighborhood(self, site_idx):
        """Invalidate env cache for a site and all its neighbors."""
        self.sites[site_idx].invalidate_cache()
        for j in self.neighbors[site_idx]:
            self.sites[j].invalidate_cache()

    def get_affected_sites(self, site_idx):
        """Return set of site indices affected by a change at site_idx."""
        affected = {site_idx}
        for j in self.neighbors[site_idx]:
            affected.add(j)
            # Second-shell for events that check NN of NN
            for k in self.neighbors[j]:
                affected.add(k)
        return affected

    # ── Bulk observables ──────────────────────────────────────────────

    def get_coverage(self, species_id=None):
        """Fractional coverage of a species (or total if None)."""
        if species_id is not None:
            return sum(1 for s in self.sites if s.adsorbate == species_id) / self.n_sites
        return sum(1 for s in self.sites if s.is_occupied) / self.n_sites

    def get_composition(self, layer=None):
        """Surface composition as {atom_type: fraction}."""
        if layer is not None:
            subset = [s for s in self.sites if s.layer == layer]
        else:
            subset = self.sites
        n = len(subset)
        if n == 0:
            return {}
        comp = defaultdict(int)
        for s in subset:
            comp[s.atom_type] += 1
        return {k: v / n for k, v in comp.items()}

    def get_site_type_distribution(self):
        """Count sites by (atom_type, adsorbate) combination."""
        dist = defaultdict(int)
        for s in self.sites:
            key = (s.atom_type, self.species_name(s.adsorbate))
            dist[key] += 1
        return dict(dist)

    # ── Factory methods ───────────────────────────────────────────────

    @classmethod
    def fcc111(cls, composition, size=(10, 10), lattice_const=None,
               n_layers=1, periodic=True):
        """
        Build an fcc(111) surface slab with given composition.

        Parameters
        ----------
        composition : dict
            e.g. {'Pd': 0.5, 'Au': 0.5} — mole fractions.
        size : tuple of int
            (nx, ny) number of unit cells.
        lattice_const : float, optional
            Lattice constant in Angstrom. If None, uses Vegard's law average.
        n_layers : int
            Number of atomic layers (1=surface only, 2=surface+subsurface).
        periodic : bool
            Periodic boundary conditions.

        Returns
        -------
        DynamicSurface
        """
        nx, ny = size
        elements = list(composition.keys())
        fractions = np.array([composition[e] for e in elements])
        fractions = fractions / fractions.sum()

        # Vegard's law lattice constant
        if lattice_const is None:
            a_vals = {'Pd': 3.89, 'Au': 4.08, 'Pt': 3.92, 'Cu': 3.61,
                      'Ag': 4.09, 'Ni': 3.52, 'Ir': 3.84, 'Ru': 2.71,
                      'Rh': 3.80, 'Co': 3.54}
            lattice_const = sum(fractions[i] * a_vals.get(elements[i], 3.9)
                                for i in range(len(elements)))

        nn_dist = lattice_const / np.sqrt(2)  # NN distance on fcc(111)

        surface = cls()

        # Build sites for each layer
        site_grid = {}  # (layer, ix, iy) -> site_index
        for layer in range(n_layers):
            for ix in range(nx):
                for iy in range(ny):
                    # Random composition assignment
                    atom_type = np.random.choice(elements, p=fractions)
                    # Hex position
                    x = ix * nn_dist + (iy % 2) * nn_dist * 0.5
                    y = iy * nn_dist * np.sqrt(3) / 2
                    pos = np.array([x, y])
                    frozen = (layer > 0)  # freeze subsurface by default
                    idx = surface.add_site(
                        atom_type=atom_type,
                        position=pos,
                        layer=layer,
                        frozen=frozen,
                    )
                    site_grid[(layer, ix, iy)] = idx

        # Build NN connectivity
        # fcc(111) in-plane: 6 neighbors
        nn_offsets_even = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (1, 1)]
        nn_offsets_odd = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (-1, 1)]

        for layer in range(n_layers):
            for ix in range(nx):
                for iy in range(ny):
                    idx = site_grid[(layer, ix, iy)]
                    offsets = nn_offsets_even if iy % 2 == 0 else nn_offsets_odd
                    for dx, dy in offsets:
                        jx = (ix + dx) % nx if periodic else ix + dx
                        jy = (iy + dy) % ny if periodic else iy + dy
                        if not periodic and (jx < 0 or jx >= nx or jy < 0 or jy >= ny):
                            continue
                        key = (layer, jx, jy)
                        if key in site_grid:
                            surface.add_edge(idx, site_grid[key])

                    # Inter-layer connections (3 neighbors to layer below)
                    if layer > 0:
                        below_layer = layer - 1
                        for dx, dy in [(0, 0), (1, 0), (0, 1)]:
                            jx = (ix + dx) % nx if periodic else ix + dx
                            jy = (iy + dy) % ny if periodic else iy + dy
                            if not periodic and (jx < 0 or jx >= nx or jy < 0 or jy >= ny):
                                continue
                            key = (below_layer, jx, jy)
                            if key in site_grid:
                                surface.add_edge(idx, site_grid[key])

        return surface

    @classmethod
    def square(cls, composition, size=(10, 10), lattice_const=3.0,
               periodic=True):
        """
        Build a simple square lattice surface (for testing).

        Parameters
        ----------
        composition : dict
            e.g. {'Pd': 0.5, 'Au': 0.5}
        size : tuple of int
        """
        nx, ny = size
        elements = list(composition.keys())
        fractions = np.array([composition[e] for e in elements])
        fractions = fractions / fractions.sum()

        surface = cls()
        site_grid = {}

        for ix in range(nx):
            for iy in range(ny):
                atom_type = np.random.choice(elements, p=fractions)
                pos = np.array([ix * lattice_const, iy * lattice_const])
                idx = surface.add_site(atom_type=atom_type, position=pos)
                site_grid[(ix, iy)] = idx

        # 4 NN
        for ix in range(nx):
            for iy in range(ny):
                idx = site_grid[(ix, iy)]
                for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                    jx = (ix + dx) % nx if periodic else ix + dx
                    jy = (iy + dy) % ny if periodic else iy + dy
                    if not periodic and (jx < 0 or jx >= nx or jy < 0 or jy >= ny):
                        continue
                    surface.add_edge(idx, site_grid[(jx, jy)])

        return surface

    # ── I/O ───────────────────────────────────────────────────────────

    def summary(self):
        """Print surface summary."""
        print(f"DynamicSurface: {self.n_sites} sites")
        print(f"  Composition: {self.get_composition()}")
        print(f"  Coverage: {self.get_coverage():.4f}")
        print(f"  Adsorbate species: {self.species_names}")
        if self.n_sites > 0:
            cn = [self.coordination_number(i) for i in range(self.n_sites)]
            print(f"  Coordination: mean={np.mean(cn):.1f}, "
                  f"range=[{min(cn)}, {max(cn)}]")

    def __repr__(self):
        comp = self.get_composition()
        comp_str = '/'.join(f"{k}:{v:.0%}" for k, v in comp.items())
        return f"DynamicSurface(n={self.n_sites}, {comp_str})"
