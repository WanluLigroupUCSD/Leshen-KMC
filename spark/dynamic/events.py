"""
Event system — catalytic and structural events with unified interface.

All events share a common base class and can coexist in the same event pool.
The KMC engine selects among both catalytic and structural events at each step.

Event types:
  Catalytic:  Adsorption, Desorption, SurfaceReaction, Diffusion
  Structural: SiteConversion, Segregation
"""

from enum import Enum, auto


class EventType(Enum):
    """Classification of event types."""
    ADSORPTION = auto()
    DESORPTION = auto()
    SURFACE_REACTION = auto()
    DIFFUSION = auto()
    SITE_CONVERSION = auto()
    SEGREGATION = auto()


class Event:
    """
    Base class for all KMC events.

    An event describes what happens (type), where (site indices),
    and how fast (rate). The rate is assigned externally by a
    RateEstimator based on the local environment.

    Parameters
    ----------
    event_type : EventType
    site : int
        Primary site index.
    rate : float
        Rate constant in s^-1 (assigned by RateEstimator).
    name : str, optional
        Human-readable label.
    """

    __slots__ = ('event_type', 'site', 'rate', 'name')

    def __init__(self, event_type, site, rate=0.0, name=''):
        self.event_type = event_type
        self.site = site
        self.rate = rate
        self.name = name

    def execute(self, surface):
        """
        Execute this event on the surface, modifying state in-place.

        Parameters
        ----------
        surface : DynamicSurface

        Returns
        -------
        affected : set of int
            Site indices whose local environment changed.
        """
        raise NotImplementedError

    def is_possible(self, surface):
        """
        Check whether this event can currently occur.

        Parameters
        ----------
        surface : DynamicSurface

        Returns
        -------
        bool
        """
        raise NotImplementedError

    def __repr__(self):
        return f"{self.name or self.event_type.name}@{self.site}(k={self.rate:.2e})"


# ── Catalytic Events ──────────────────────────────────────────────────


class CatalyticEvent(Event):
    """Base for catalytic events."""
    pass


class Adsorption(CatalyticEvent):
    """
    Adsorption: gas species → surface site.

    Conditions: site must be empty.
    Action: set adsorbate on site.

    Parameters
    ----------
    site : int
    species_id : int
        Adsorbate species to place.
    rate : float
    name : str
    """

    __slots__ = ('species_id',)

    def __init__(self, site, species_id, rate=0.0, name=''):
        super().__init__(EventType.ADSORPTION, site, rate,
                         name or f'ads_sp{species_id}')
        self.species_id = species_id

    def is_possible(self, surface):
        return surface.sites[self.site].is_empty

    def execute(self, surface):
        surface.set_adsorbate(self.site, self.species_id)
        return surface.get_affected_sites(self.site)


class Desorption(CatalyticEvent):
    """
    Desorption: surface adsorbate → gas phase.

    Conditions: site must have the specified adsorbate.
    Action: clear adsorbate.
    """

    __slots__ = ('species_id',)

    def __init__(self, site, species_id, rate=0.0, name=''):
        super().__init__(EventType.DESORPTION, site, rate,
                         name or f'des_sp{species_id}')
        self.species_id = species_id

    def is_possible(self, surface):
        return surface.sites[self.site].adsorbate == self.species_id

    def execute(self, surface):
        surface.clear_adsorbate(self.site)
        return surface.get_affected_sites(self.site)


class SurfaceReaction(CatalyticEvent):
    """
    Surface reaction: reactant(s) at site(s) → product(s).

    Parameters
    ----------
    site : int
        Primary site.
    partner_site : int or None
        Second site for bimolecular reactions (e.g., Langmuir-Hinshelwood).
    reactant_species : tuple of int
        Species consumed (at site, partner_site).
    product_species : tuple of int
        Species produced (at site, partner_site). Use 0 for empty.
    rate : float
    tof_count : dict
        TOF counters (e.g., {'CO2_production': 1}).
    """

    __slots__ = ('partner_site', 'reactant_species', 'product_species',
                 'tof_count')

    def __init__(self, site, partner_site=None,
                 reactant_species=(0,), product_species=(0,),
                 rate=0.0, name='', tof_count=None):
        super().__init__(EventType.SURFACE_REACTION, site, rate, name)
        self.partner_site = partner_site
        self.reactant_species = tuple(reactant_species)
        self.product_species = tuple(product_species)
        self.tof_count = tof_count or {}

    def is_possible(self, surface):
        # Check primary site
        if surface.sites[self.site].adsorbate != self.reactant_species[0]:
            return False
        # Check partner site
        if self.partner_site is not None and len(self.reactant_species) > 1:
            if surface.sites[self.partner_site].adsorbate != self.reactant_species[1]:
                return False
        return True

    def execute(self, surface):
        # Apply products
        surface.set_adsorbate(self.site, self.product_species[0])
        affected = surface.get_affected_sites(self.site)
        if self.partner_site is not None and len(self.product_species) > 1:
            surface.set_adsorbate(self.partner_site, self.product_species[1])
            affected |= surface.get_affected_sites(self.partner_site)
        return affected


class Diffusion(CatalyticEvent):
    """
    Surface diffusion: adsorbate hops from site to neighbor.

    Conditions: source occupied, target empty.
    Action: move adsorbate.
    """

    __slots__ = ('target_site',)

    def __init__(self, site, target_site, rate=0.0, name=''):
        super().__init__(EventType.DIFFUSION, site, rate,
                         name or f'diff_{site}->{target_site}')
        self.target_site = target_site

    def is_possible(self, surface):
        return (surface.sites[self.site].is_occupied
                and surface.sites[self.target_site].is_empty)

    def execute(self, surface):
        sp = surface.sites[self.site].adsorbate
        surface.clear_adsorbate(self.site)
        surface.set_adsorbate(self.target_site, sp)
        affected = surface.get_affected_sites(self.site)
        affected |= surface.get_affected_sites(self.target_site)
        return affected


# ── Structural Events ─────────────────────────────────────────────────


class StructuralEvent(Event):
    """Base for structural evolution events."""
    pass


class SiteConversion(StructuralEvent):
    """
    Site conversion: change atom type at a site.

    Models: motif switching, local oxidation-state change, reconstruction.
    Conditions: site must be of the specified original type.
    """

    __slots__ = ('from_type', 'to_type')

    def __init__(self, site, from_type, to_type, rate=0.0, name=''):
        super().__init__(EventType.SITE_CONVERSION, site, rate,
                         name or f'convert_{from_type}->{to_type}')
        self.from_type = from_type
        self.to_type = to_type

    def is_possible(self, surface):
        site = surface.sites[self.site]
        return site.atom_type == self.from_type and not site.frozen

    def execute(self, surface):
        surface.convert_site(self.site, self.to_type)
        return surface.get_affected_sites(self.site)


class Segregation(StructuralEvent):
    """
    Segregation: swap atom types between a surface site and a neighbor.

    Models: thermodynamic or adsorbate-driven segregation.
    Conditions: sites must have different atom types, neither frozen.
    """

    __slots__ = ('partner_site',)

    def __init__(self, site, partner_site, rate=0.0, name=''):
        super().__init__(EventType.SEGREGATION, site, rate,
                         name or f'seg_{site}<->{partner_site}')
        self.partner_site = partner_site

    def is_possible(self, surface):
        s1 = surface.sites[self.site]
        s2 = surface.sites[self.partner_site]
        return (s1.atom_type != s2.atom_type
                and not s1.frozen and not s2.frozen)

    def execute(self, surface):
        surface.swap_atoms(self.site, self.partner_site)
        affected = surface.get_affected_sites(self.site)
        affected |= surface.get_affected_sites(self.partner_site)
        return affected


# ── Event Generator ───────────────────────────────────────────────────


class EventGenerator:
    """
    Generates all possible events for a site based on the surface state.

    This is the bridge between the surface representation and the
    rate estimator. For each site, it enumerates all events that
    are currently possible given the local environment.

    Parameters
    ----------
    species_list : list of int
        Adsorbate species IDs that can adsorb.
    enable_diffusion : bool
        Generate diffusion events.
    enable_segregation : bool
        Generate segregation events between different-type neighbors.
    enable_site_conversion : bool
        Generate site conversion events.
    conversion_rules : list of (str, str)
        Allowed atom type conversions, e.g. [('Pd', 'Au'), ('Au', 'Pd')].
    """

    def __init__(self, species_list=None, enable_diffusion=True,
                 enable_segregation=True, enable_site_conversion=False,
                 conversion_rules=None):
        self.species_list = species_list or []
        self.enable_diffusion = enable_diffusion
        self.enable_segregation = enable_segregation
        self.enable_site_conversion = enable_site_conversion
        self.conversion_rules = conversion_rules or []

    def generate(self, surface, site_idx):
        """
        Generate all possible events at a site.

        Parameters
        ----------
        surface : DynamicSurface
        site_idx : int

        Returns
        -------
        events : list of Event
            Events with rate=0 (rates assigned later by RateEstimator).
        """
        events = []
        site = surface.sites[site_idx]

        # ── Catalytic events ──

        if site.is_empty:
            # Adsorption of each species
            for sp_id in self.species_list:
                events.append(Adsorption(site_idx, sp_id))
        else:
            # Desorption
            events.append(Desorption(site_idx, site.adsorbate))

            # Diffusion to empty neighbors
            if self.enable_diffusion:
                for j in surface.neighbors[site_idx]:
                    if surface.sites[j].is_empty:
                        events.append(Diffusion(site_idx, j))

            # Bimolecular reactions with occupied neighbors
            for j in surface.neighbors[site_idx]:
                nn = surface.sites[j]
                if nn.is_occupied and nn.adsorbate != site.adsorbate:
                    # Generic LH reaction: A* + B* → product
                    events.append(SurfaceReaction(
                        site=site_idx,
                        partner_site=j,
                        reactant_species=(site.adsorbate, nn.adsorbate),
                        product_species=(0, 0),  # both desorb
                        name=f'rxn_sp{site.adsorbate}+sp{nn.adsorbate}',
                    ))

        # ── Structural events ──

        if not site.frozen:
            # Segregation with different-type neighbors
            if self.enable_segregation:
                for j in surface.neighbors[site_idx]:
                    nn = surface.sites[j]
                    if (nn.atom_type != site.atom_type
                            and not nn.frozen):
                        events.append(Segregation(site_idx, j))

            # Site conversion
            if self.enable_site_conversion:
                for from_t, to_t in self.conversion_rules:
                    if site.atom_type == from_t:
                        events.append(SiteConversion(site_idx, from_t, to_t))

        return events

    def generate_all(self, surface):
        """Generate events for all sites. Returns flat list."""
        all_events = []
        for i in range(surface.n_sites):
            all_events.extend(self.generate(surface, i))
        return all_events
