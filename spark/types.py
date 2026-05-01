"""Data model definitions for KMC model construction (kmos-compatible API)."""

import numpy as np
import json


class Species:
    """A chemical species that can occupy a lattice site."""

    def __init__(self, name, color=None, representation=None, tags=None):
        self.name = name
        self.color = color
        self.representation = representation
        self.tags = tags or ''
        self.id = None  # assigned when added to Project


class Site:
    """A site within a lattice layer."""

    def __init__(self, name, pos=(0.5, 0.5, 0.5), default_species='empty',
                 site_type=None):
        self.name = name
        self.pos = tuple(pos) if not isinstance(pos, tuple) else pos
        self.default_species = default_species
        self.site_type = site_type  # int or str label for site type


class Layer:
    """A layer in the lattice (collection of sites)."""

    def __init__(self, name='default'):
        self.name = name
        self.sites = []

    def add_site(self, name, pos=(0.5, 0.5, 0.5), default_species='empty',
                 site_type=None):
        site = Site(name, pos, default_species, site_type=site_type)
        self.sites.append(site)
        return site


class Coord:
    """A coordinate reference for conditions/actions."""

    def __init__(self, name=None, offset=(0, 0, 0), layer=None, site=None):
        self.name = name
        self.offset = tuple(offset)
        self.layer = layer
        self.site = site

    def __repr__(self):
        return f"Coord({self.site}, offset={self.offset}, layer={self.layer})"


class Condition:
    """A condition that must be met for a process to occur."""

    def __init__(self, coord, species):
        self.coord = coord
        self.species = species

    def __repr__(self):
        return f"Condition({self.species}@{self.coord})"


class Action:
    """An action that changes species at a site when a process executes."""

    def __init__(self, coord, species):
        self.coord = coord
        self.species = species

    def __repr__(self):
        return f"Action({self.species}@{self.coord})"


class Parameter:
    """An adjustable model parameter."""

    def __init__(self, name, value=0.0, adjustable=False,
                 min=0.0, max=0.0, scale='linear'):
        self.name = name
        self.value = value
        self.adjustable = adjustable
        self.min = min
        self.max = max
        self.scale = scale


class Process:
    """An elementary process/reaction in the KMC model."""

    def __init__(self, name, conditions=None, actions=None,
                 rate_constant='0', tof_count=None, enabled=True,
                 reverse_of=None, site_type=None):
        self.name = name
        self.conditions = conditions or []
        self.actions = actions or []
        self.rate_constant = rate_constant
        self.tof_count = tof_count or {}
        self.enabled = enabled
        self.reverse_of = reverse_of  # name of the reverse process (for BEP)
        self.site_type = site_type  # required site type (int) or None
        self.id = None  # assigned when added to Project


class LateralInteraction:
    """Pairwise lateral interaction between neighboring adsorbates.

    Parameters
    ----------
    species1, species2 : str
        Names of the two interacting species.
    energy : float
        Interaction energy in eV. Positive = repulsive.
    """

    def __init__(self, species1, species2, energy):
        self.species1 = species1
        self.species2 = species2
        self.energy = energy

    def __repr__(self):
        return (f"LateralInteraction({self.species1}-{self.species2}, "
                f"E={self.energy} eV)")


class BEPRelation:
    """Bronsted-Evans-Polanyi relation for a process.

    Modifies activation energy based on reaction enthalpy change:
        Ea(env) = Ea(0) + alpha * (DeltaH(env) - DeltaH(0))

    Parameters
    ----------
    process_name : str
        Name of the process this BEP applies to.
    alpha : float
        BEP slope (proximity factor), typically 0.0-1.0. Default 0.5.
    """

    def __init__(self, process_name, alpha=0.5):
        self.process_name = process_name
        self.alpha = alpha

    def get_reverse(self):
        """Return BEP relation for the reverse process."""
        return BEPRelation(self.process_name + '_rev', alpha=1.0 - self.alpha)

    def delta_Ea(self, delta_delta_H):
        """Compute change in activation energy from change in reaction enthalpy."""
        return self.alpha * delta_delta_H

    def __repr__(self):
        return f"BEPRelation({self.process_name}, alpha={self.alpha})"


class Lattice:
    """Lattice geometry definition.

    Multi-lattice support (Hoffmann-Reuter-Scheffler 2015): all coexisting
    sub-lattices share one physical ``cell`` and live as separate ``Layer``
    objects in ``layers``. The user supplies a commensurate super-cell large
    enough to embed every primitive cell. ``default_layer`` and
    ``substrate_layer`` follow kmcos convention (former is the layer the
    model starts in by default; latter is the always-present scaffold).
    """

    def __init__(self):
        self.cell = np.diag([1.0, 1.0, 1.0])
        self.layers = []
        self.default_layer = None
        self.substrate_layer = None

    @property
    def spuck(self):
        """Sites Per Unit Cell — total sites across all layers."""
        return sum(len(layer.sites) for layer in self.layers)

    def layer_offset(self, layer_name):
        """Return the cumulative site offset within a unit cell at which
        ``layer_name`` begins. Layer 0's sites occupy [0, len(layer0.sites));
        layer 1's sites occupy [len(layer0.sites), len(layer0.sites)+len(layer1.sites)); etc.
        """
        offset = 0
        for layer in self.layers:
            if layer.name == layer_name:
                return offset
            offset += len(layer.sites)
        raise KeyError(f"Layer '{layer_name}' not found")

    def site_in_cell_id(self, layer_name, site_name):
        """Resolve (layer_name, site_name) to its index in [0, spuck)."""
        offset = self.layer_offset(layer_name)
        for layer in self.layers:
            if layer.name == layer_name:
                for i, site in enumerate(layer.sites):
                    if site.name == site_name:
                        return offset + i
                raise KeyError(
                    f"Site '{site_name}' not found in layer '{layer_name}'")
        raise KeyError(f"Layer '{layer_name}' not found")

    def site_in_cell_to_layer(self, s):
        """Inverse: which layer does site-in-cell index ``s`` belong to?
        Returns (layer_name, site_in_layer_idx).
        """
        offset = 0
        for layer in self.layers:
            n = len(layer.sites)
            if s < offset + n:
                return layer.name, s - offset
            offset += n
        raise IndexError(f"site_in_cell index {s} out of range [0, {offset})")

    def generate_coord(self, description):
        """
        Parse coordinate description like 'hollow.(0,0,0).simple_cubic'
        or just 'hollow' for the default site at offset (0,0,0).
        """
        parts = description.split('.')
        site_name = parts[0]
        if len(parts) > 1:
            offset = eval(parts[1])
        else:
            offset = (0, 0, 0)
        layer_name = parts[2] if len(parts) > 2 else (
            self.layers[0].name if self.layers else 'default')
        return Coord(name=description, offset=offset,
                     layer=layer_name, site=site_name)


class Project:
    """
    Container for a complete KMC model definition.

    Usage follows kmos conventions:
        pt = Project()
        pt.set_meta(model_name='my_model', model_dimension=2)
        pt.add_species(name='empty')
        pt.add_species(name='CO')
        layer = pt.add_layer(name='surface')
        layer.sites.append(Site(name='hollow'))
        coord = pt.lattice.generate_coord('hollow')
        pt.add_process(name='adsorption',
                       conditions=[Condition(coord, 'empty')],
                       actions=[Action(coord, 'CO')],
                       rate_constant='p_CO*bar*A/sqrt(2*pi*m_CO*umass/beta)')

        # Lateral interactions
        pt.add_lateral_interaction('CO', 'CO', energy=0.10)

        # BEP relation
        pt.add_bep_relation('CO_desorption', alpha=0.5)

        # Diffusion
        pt.add_diffusion('CO', rate_constant='1e8*exp(-0.5*eV*beta)')
    """

    def __init__(self):
        self.meta = {
            'author': '',
            'email': '',
            'model_name': '',
            'model_dimension': 2,
        }
        self.species_list = []
        self.species_map = {}
        self.parameter_list = []
        self.parameter_map = {}
        self.process_list = []
        self.lattice = Lattice()
        self.filename = None
        self.lateral_interactions = []
        self.bep_relations = {}  # process_name -> BEPRelation

    def set_meta(self, **kwargs):
        self.meta.update(kwargs)

    def add_species(self, name, **kwargs):
        """Add a species. First species added becomes the default."""
        sp = Species(name=name, **kwargs)
        sp.id = len(self.species_list)
        self.species_list.append(sp)
        self.species_map[name] = sp
        return sp

    def add_layer(self, name=None, layer=None):
        if layer is None:
            layer = Layer(name=name or 'default')
        self.lattice.layers.append(layer)
        if self.lattice.default_layer is None:
            self.lattice.default_layer = layer.name
        if self.lattice.substrate_layer is None:
            self.lattice.substrate_layer = layer.name
        return layer

    def add_parameter(self, name, value=0.0, **kwargs):
        p = Parameter(name=name, value=value, **kwargs)
        self.parameter_list.append(p)
        self.parameter_map[name] = p
        return p

    def add_process(self, name, **kwargs):
        proc = Process(name=name, **kwargs)
        proc.id = len(self.process_list)
        self.process_list.append(proc)
        return proc

    def parse_and_add_process(self, description):
        """
        Parse shorthand: 'name; species@site->species@site; rate_expr'
        Example: 'CO_desorption; CO@hollow->empty@hollow; k_des*exp(-E_des*beta*eV)'
        """
        parts = [p.strip() for p in description.split(';')]
        if len(parts) < 3:
            raise ValueError(
                "Expected format: 'name; conditions->actions; rate_constant'")

        name = parts[0]
        rate_constant = parts[2]

        transitions = parts[1].split('->')
        cond_str = transitions[0].strip()
        act_str = transitions[1].strip() if len(transitions) > 1 else ''

        conditions = []
        actions = []

        for item in cond_str.split('+'):
            item = item.strip()
            if not item:
                continue
            if '@' in item:
                sp, site = item.split('@')
                coord = self.lattice.generate_coord(site.strip())
                conditions.append(Condition(coord, sp.strip()))

        for item in act_str.split('+'):
            item = item.strip()
            if not item:
                continue
            if '@' in item:
                sp, site = item.split('@')
                coord = self.lattice.generate_coord(site.strip())
                actions.append(Action(coord, sp.strip()))

        return self.add_process(name=name, conditions=conditions,
                                actions=actions, rate_constant=rate_constant)

    # ------------------------------------------------------------------
    # Lateral interactions
    # ------------------------------------------------------------------

    def add_lateral_interaction(self, species1, species2, energy):
        """
        Add a pairwise lateral interaction between nearest-neighbor adsorbates.

        Parameters
        ----------
        species1, species2 : str
            Names of interacting species.
        energy : float
            Interaction energy in eV. Positive = repulsive, negative = attractive.
        """
        li = LateralInteraction(species1, species2, energy)
        self.lateral_interactions.append(li)
        return li

    # ------------------------------------------------------------------
    # BEP relations
    # ------------------------------------------------------------------

    def add_bep_relation(self, process_name, alpha=0.5):
        """
        Add a BEP relation for a process.

        The activation energy is modified by:
            Ea(env) = Ea(0) + alpha * (DeltaH(env) - DeltaH(0))

        where DeltaH(env) includes lateral interaction contributions.

        Parameters
        ----------
        process_name : str
            Name of the process.
        alpha : float
            BEP slope (proximity factor). Default 0.5.
        """
        bep = BEPRelation(process_name, alpha=alpha)
        self.bep_relations[process_name] = bep
        return bep

    # ------------------------------------------------------------------
    # Diffusion helper
    # ------------------------------------------------------------------

    def add_diffusion(self, species, rate_constant, empty_species='empty',
                      name_prefix=None, tof_count=None, site_type=None):
        """
        Add diffusion processes for a species in all nearest-neighbor directions.

        Generates one process per NN direction:
            species@center + empty@neighbor -> empty@center + species@neighbor

        Parameters
        ----------
        species : str
            The diffusing species name.
        rate_constant : str or float
            Rate constant expression for the hop.
        empty_species : str
            Name of the empty/vacant species. Default 'empty'.
        name_prefix : str, optional
            Prefix for process names. Default: '{species}_diff'.
        tof_count : dict, optional
            TOF counters for the diffusion process.
        site_type : int, optional
            Required site type.

        Returns
        -------
        list of Process
            The generated diffusion processes.
        """
        ndim = self.meta.get('model_dimension', 2)
        prefix = name_prefix or f'{species}_diff'

        if ndim == 1:
            offsets = [(1, 0, 0), (-1, 0, 0)]
            labels = ['right', 'left']
        elif ndim == 2:
            offsets = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
            labels = ['right', 'left', 'up', 'down']
        else:
            offsets = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                       (0, 0, 1), (0, 0, -1)]
            labels = ['right', 'left', 'up', 'down', 'front', 'back']

        site_name = 'default'
        if self.lattice.layers:
            layer = self.lattice.layers[0]
            if layer.sites:
                site_name = layer.sites[0].name

        processes = []
        for offset, label in zip(offsets, labels):
            center = self.lattice.generate_coord(site_name)
            neighbor = Coord(offset=offset, site=site_name,
                             layer=center.layer)

            proc = self.add_process(
                name=f'{prefix}_{label}',
                conditions=[
                    Condition(center, species),
                    Condition(neighbor, empty_species),
                ],
                actions=[
                    Action(center, empty_species),
                    Action(neighbor, species),
                ],
                rate_constant=rate_constant,
                tof_count=tof_count or {},
                site_type=site_type,
            )
            processes.append(proc)

        return processes

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------

    def save(self, filename=None):
        """Save project to JSON file."""
        fn = filename or self.filename
        if fn is None:
            raise ValueError("No filename specified")
        data = {
            'meta': self.meta,
            'species': [{'name': s.name, 'color': s.color, 'tags': s.tags}
                        for s in self.species_list],
            'parameters': [{'name': p.name, 'value': p.value,
                            'adjustable': p.adjustable,
                            'min': p.min, 'max': p.max}
                           for p in self.parameter_list],
            'processes': [{'name': p.name,
                           'rate_constant': p.rate_constant,
                           'tof_count': p.tof_count,
                           'reverse_of': p.reverse_of,
                           'site_type': p.site_type,
                           'conditions': [
                               {'offset': c.coord.offset, 'species': c.species}
                               for c in p.conditions],
                           'actions': [
                               {'offset': a.coord.offset, 'species': a.species}
                               for a in p.actions]}
                          for p in self.process_list],
            'lattice': {
                'cell': self.lattice.cell.tolist()
                if isinstance(self.lattice.cell, np.ndarray)
                else self.lattice.cell,
                'layers': [
                    {'name': l.name,
                     'sites': [{'name': s.name, 'pos': s.pos,
                                'default_species': s.default_species,
                                'site_type': s.site_type}
                               for s in l.sites]}
                    for l in self.lattice.layers],
            },
            'lateral_interactions': [
                {'species1': li.species1, 'species2': li.species2,
                 'energy': li.energy}
                for li in self.lateral_interactions],
            'bep_relations': [
                {'process_name': b.process_name, 'alpha': b.alpha}
                for b in self.bep_relations.values()],
        }
        with open(fn, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Model saved to {fn}")

    def summary(self):
        """Print model summary."""
        print(f"Model: {self.meta.get('model_name', 'unnamed')}")
        print(f"  Dimension: {self.meta.get('model_dimension', 2)}")
        print(f"  Species ({len(self.species_list)}): "
              f"{', '.join(s.name for s in self.species_list)}")
        print(f"  Parameters ({len(self.parameter_list)}): "
              f"{', '.join(p.name for p in self.parameter_list)}")
        print(f"  Processes ({len(self.process_list)}):")
        for p in self.process_list:
            conds = ', '.join(
                f"{c.species}@{c.coord.offset}" for c in p.conditions)
            acts = ', '.join(
                f"{a.species}@{a.coord.offset}" for a in p.actions)
            extra = ''
            if p.site_type is not None:
                extra += f' [site_type={p.site_type}]'
            if p.reverse_of:
                extra += f' [reverse_of={p.reverse_of}]'
            print(f"    {p.name}: [{conds}] -> [{acts}]{extra}")
            print(f"      rate = {p.rate_constant}")
        if self.lateral_interactions:
            print(f"  Lateral interactions ({len(self.lateral_interactions)}):")
            for li in self.lateral_interactions:
                sign = 'repulsive' if li.energy > 0 else 'attractive'
                print(f"    {li.species1}-{li.species2}: "
                      f"{li.energy:+.4f} eV ({sign})")
        if self.bep_relations:
            print(f"  BEP relations ({len(self.bep_relations)}):")
            for name, bep in self.bep_relations.items():
                print(f"    {name}: alpha={bep.alpha}")
