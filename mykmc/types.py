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

    def __init__(self, name, pos=(0.5, 0.5, 0.5), default_species='empty'):
        self.name = name
        self.pos = tuple(pos) if not isinstance(pos, tuple) else pos
        self.default_species = default_species


class Layer:
    """A layer in the lattice (collection of sites)."""

    def __init__(self, name='default'):
        self.name = name
        self.sites = []

    def add_site(self, name, pos=(0.5, 0.5, 0.5), default_species='empty'):
        site = Site(name, pos, default_species)
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
                 rate_constant='0', tof_count=None, enabled=True):
        self.name = name
        self.conditions = conditions or []
        self.actions = actions or []
        self.rate_constant = rate_constant
        self.tof_count = tof_count or {}
        self.enabled = enabled
        self.id = None  # assigned when added to Project


class Lattice:
    """Lattice geometry definition."""

    def __init__(self):
        self.cell = np.diag([1.0, 1.0, 1.0])
        self.layers = []

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
                                'default_species': s.default_species}
                               for s in l.sites]}
                    for l in self.lattice.layers],
            },
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
            print(f"    {p.name}: [{conds}] -> [{acts}]")
            print(f"      rate = {p.rate_constant}")
