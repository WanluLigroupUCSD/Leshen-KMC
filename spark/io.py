"""
Model I/O for SPARK.

Primary format: VASP-style TAG = VALUE (.sparkin)

    SYSTEM    = HER_Pt111
    SURFACE   = Pt(111)
    TEMP      = 298              (K)
    POTENTIAL = -0.2             (V vs RHE)
    SPECIES   = H
    REACTION  = Volmer_fwd : * -> H : 0.67 : electrochemical : 0.5
    LATERAL   = H-H : 0.10
    BEP       = Volmer_fwd : 0.5

Also supports YAML (.yaml) and JSON (.json).

Usage:
    from spark import load_model, KMCEngine
    pt = load_model('her.sparkin')
    engine = KMCEngine(pt, size=[20, 20])
"""

import os
import re
import numpy as np
from .types import (
    Project, Site, Layer, Condition, Action, Coord, Parameter,
)

# ============================================================================
#  Surface presets: name -> (lattice_constant_A, vacuum_A, site_name)
# ============================================================================

SURFACE_PRESETS = {
    'Pt(111)': (2.775, 15.0, 'top'),
    'Pt(100)': (2.775, 15.0, 'hollow'),
    'Cu(111)': (2.556, 15.0, 'top'),
    'Cu(100)': (2.556, 15.0, 'hollow'),
    'Pd(111)': (2.751, 15.0, 'top'),
    'Pd(100)': (2.751, 15.0, 'hollow'),
    'Au(111)': (2.884, 15.0, 'top'),
    'Au(100)': (2.884, 15.0, 'hollow'),
    'Ag(111)': (2.889, 15.0, 'top'),
    'Ni(111)': (2.492, 15.0, 'top'),
    'Ru(0001)': (2.706, 15.0, 'top'),
    'Ir(111)': (2.715, 15.0, 'top'),
    'Mo(110)': (2.725, 15.0, 'top'),
    'Fe(110)': (2.482, 15.0, 'top'),
}

# Tags that can appear multiple times
_MULTI_TAGS = {'REACTION', 'DIFFUSION', 'LATERAL', 'BEP', 'PARAMETER'}


# ============================================================================
#  Unified loader
# ============================================================================

def load_model(path):
    """
    Load a SPARK model from file. Auto-detects format by extension.

    Supported: .sparkin (VASP-style), .yaml/.yml, .json
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.yaml', '.yml'):
        return load_yaml(path)
    elif ext == '.json':
        return _load_json(path)
    else:
        return load_sparkin(path)


# ============================================================================
#  SPARKIN format: VASP-style TAG = VALUE
# ============================================================================

def load_sparkin(path):
    """
    Load model from VASP-style .sparkin file.

    Tags:
        SYSTEM    = model name
        SURFACE   = preset (Pt(111), Cu(100), ...)
        LATTICE   = a b c  (custom unit cell, Angstrom)
        TEMP      = temperature (K)
        POTENTIAL = applied potential (V vs RHE)
        DIMENSION = 1, 2, or 3 (default 2)
        SPECIES   = sp1 sp2 ...  (empty is implicit)
        PARAMETER = name  value
        REACTION  = name : equation : Ea : [type : beta]
        REACTION  = name : equation : rate=expression
        DIFFUSION = species : Ea
        LATERAL   = sp1-sp2 : energy(eV)
        BEP       = process : alpha
        LSIZE     = nx ny
        NEQUIL    = equilibration steps
        NSTEPS    = production steps
    """
    with open(path, 'r') as f:
        lines = f.readlines()
    return _parse_sparkin(lines, source=path)


def loads_sparkin(text):
    """Load model from a SPARKIN-format string."""
    return _parse_sparkin(text.splitlines(), source='<string>')


def _strip_inline_comment(value):
    """Remove VASP-style inline comments: (description) at end of value."""
    # Remove trailing (...) comment
    result = re.sub(r'\s*\([^)]*\)\s*$', '', value)
    return result.strip()


def _parse_sparkin(lines, source='<unknown>'):
    """Parse VASP-style TAG = VALUE format."""
    tags = {}  # single-value tags
    multi = {}  # multi-value tags (REACTION, LATERAL, etc.)

    for lineno, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line or line.startswith('#') or line.startswith('!'):
            continue

        if '=' not in line:
            continue

        # Split TAG = VALUE
        idx = line.index('=')
        tag = line[:idx].strip().upper()
        value = _strip_inline_comment(line[idx + 1:])

        if tag in _MULTI_TAGS:
            multi.setdefault(tag, []).append((value, lineno))
        else:
            tags[tag] = (value, lineno)

    return _build_from_tags(tags, multi, source)


def _build_from_tags(tags, multi, source):
    """Build Project from parsed tags."""
    pt = Project()

    def get(tag, default=None):
        if tag in tags:
            return tags[tag][0]
        return default

    # --- Metadata ---
    name = get('SYSTEM', 'unnamed')
    dim = int(get('DIMENSION', '2'))
    pt.set_meta(model_name=name, model_dimension=dim)

    # --- Species ---
    species_str = get('SPECIES', '')
    pt.add_species(name='empty')  # always implicit
    for sp in species_str.split():
        if sp and sp != 'empty' and sp != '*':
            pt.add_species(name=sp)

    # --- Lattice ---
    surface = get('SURFACE')
    if surface and surface in SURFACE_PRESETS:
        a, vac, site_name = SURFACE_PRESETS[surface]
        pt.lattice.cell = np.diag([a, a, vac])
        lname = surface.replace('(', '').replace(')', '')
        layer = pt.add_layer(name=lname)
        layer.sites.append(Site(name=site_name, pos=(0.5, 0.5, 0.5),
                                default_species='empty'))
    elif 'LATTICE' in tags:
        parts = tags['LATTICE'][0].split()
        a, b, c = float(parts[0]), float(parts[1]), float(parts[2])
        pt.lattice.cell = np.diag([a, b, c])
        site_name = parts[3] if len(parts) > 3 else 'site'
        layer = pt.add_layer(name='default')
        layer.sites.append(Site(name=site_name, pos=(0.5, 0.5, 0.5),
                                default_species='empty'))
    else:
        layer = pt.add_layer(name='default')
        layer.sites.append(Site(name='site', pos=(0.5, 0.5, 0.5),
                                default_species='empty'))

    # --- Parameters ---
    temp = get('TEMP')
    if temp is not None:
        pt.add_parameter(name='T', value=float(temp))
    potential = get('POTENTIAL')
    if potential is not None:
        pt.add_parameter(name='U', value=float(potential))

    for pval, lineno in multi.get('PARAMETER', []):
        parts = pval.split(None, 1)
        pname = parts[0]
        pvalue = parts[1] if len(parts) > 1 else '0'
        try:
            pvalue = float(pvalue)
        except ValueError:
            pass
        pt.add_parameter(name=pname, value=pvalue)

    # --- Reactions ---
    for rval, lineno in multi.get('REACTION', []):
        try:
            _parse_sparkin_reaction(pt, rval)
        except Exception as e:
            raise ValueError(f"{source}:{lineno}: {e}") from e

    for dval, lineno in multi.get('DIFFUSION', []):
        try:
            _parse_sparkin_diffusion(pt, dval)
        except Exception as e:
            raise ValueError(f"{source}:{lineno}: {e}") from e

    # --- Lateral interactions ---
    for lval, lineno in multi.get('LATERAL', []):
        parts = [p.strip() for p in lval.split(':')]
        pair = parts[0]
        energy = float(parts[1])
        sp1, sp2 = pair.split('-', 1)
        pt.add_lateral_interaction(sp1.strip(), sp2.strip(), energy)

    # --- BEP ---
    for bval, lineno in multi.get('BEP', []):
        parts = [p.strip() for p in bval.split(':')]
        proc = parts[0]
        alpha = float(parts[1])
        pt.add_bep_relation(process_name=proc, alpha=alpha)

    # --- Store simulation settings as metadata ---
    lsize = get('LSIZE')
    if lsize:
        pt.meta['lsize'] = [int(x) for x in lsize.split()]
    nequil = get('NEQUIL')
    if nequil:
        pt.meta['nequil'] = int(nequil)
    nsteps = get('NSTEPS')
    if nsteps:
        pt.meta['nsteps'] = int(nsteps)
    nsample = get('NSAMPLE')
    if nsample:
        pt.meta['nsample'] = int(nsample)

    return pt


def _parse_sparkin_reaction(pt, value):
    """
    Parse REACTION value.

    Formats:
        name : equation : Ea                           (thermal)
        name : equation : Ea : electrochemical : beta  (PCET)
        name : equation : rate=expression              (custom rate)
    """
    fields = [f.strip() for f in value.split(':')]
    name = fields[0]
    equation = fields[1]

    # Check for custom rate
    if len(fields) >= 3 and fields[2].startswith('rate='):
        rate_expr = fields[2][5:]  # strip 'rate='
        conditions, actions = _parse_reaction_string(pt, equation)
        pt.add_process(name=name, conditions=conditions, actions=actions,
                       rate_constant=rate_expr, tof_count={name: 1})
        return

    Ea = float(fields[2]) if len(fields) >= 3 else 0.0
    rxn_type = fields[3].strip().lower() if len(fields) >= 4 else 'thermal'
    beta = float(fields[4]) if len(fields) >= 5 else None

    rxn_dict = {'beta': beta} if beta is not None else {}
    rate_expr = _make_rate_expr(Ea, rxn_type, rxn_dict)
    conditions, actions = _parse_reaction_string(pt, equation)

    pt.add_process(name=name, conditions=conditions, actions=actions,
                   rate_constant=rate_expr, tof_count={name: 1})


def _parse_sparkin_diffusion(pt, value):
    """Parse DIFFUSION value: species : Ea"""
    fields = [f.strip() for f in value.split(':')]
    species = fields[0]
    Ea = float(fields[1]) if len(fields) >= 2 else 0.0
    rate_expr = f'kB*T/h*exp(-{Ea}*eV/(kB*T))'
    pt.add_diffusion(species=species, rate_constant=rate_expr,
                     name_prefix=f'{species}_diff')


# ============================================================================
#  YAML format (alternative for complex models)
# ============================================================================

def load_yaml(path):
    """Load a SPARK model from a YAML file."""
    import yaml
    with open(path, 'r') as f:
        data = yaml.safe_load(f)
    return _build_project_yaml(data, source=path)


def loads_yaml(text):
    """Load a SPARK model from a YAML string."""
    import yaml
    data = yaml.safe_load(text)
    return _build_project_yaml(data, source='<string>')


def _build_project_yaml(data, source='<unknown>'):
    """Build a Project from parsed YAML data."""
    pt = Project()

    # Metadata
    name = data.get('name', data.get('model', {}).get('name', 'unnamed'))
    dim = data.get('dimension', data.get('model', {}).get('dimension', 2))
    pt.set_meta(model_name=name, model_dimension=dim)

    # Species
    species_section = data.get('species', [])
    species_names = [s if isinstance(s, str) else s['name'] for s in species_section]
    if not any(n in ('empty', '*') for n in species_names):
        pt.add_species(name='empty')
    for sp in species_section:
        if isinstance(sp, str):
            pt.add_species(name='empty' if sp == '*' else sp)
        elif isinstance(sp, dict):
            pt.add_species(name=sp['name'], color=sp.get('color'))

    # Lattice
    surface = data.get('surface')
    if surface and surface in SURFACE_PRESETS:
        a, vac, site_name = SURFACE_PRESETS[surface]
        pt.lattice.cell = np.diag([a, a, vac])
        layer = pt.add_layer(name=surface.replace('(', '').replace(')', ''))
        layer.sites.append(Site(name=site_name, pos=(0.5, 0.5, 0.5),
                                default_species='empty'))
    elif 'lattice' in data:
        lat = data['lattice']
        cell = lat.get('cell', [1, 1, 1])
        if isinstance(cell, list) and len(cell) == 3 and not isinstance(cell[0], list):
            pt.lattice.cell = np.diag(cell)
        else:
            pt.lattice.cell = np.array(cell, dtype=float)
        layer = pt.add_layer(name=lat.get('layer_name', 'default'))
        for s in lat.get('sites', [{'name': 'site'}]):
            if isinstance(s, str):
                layer.sites.append(Site(name=s, pos=(0.5, 0.5, 0.5),
                                        default_species='empty'))
            else:
                layer.sites.append(Site(
                    name=s['name'],
                    pos=tuple(s.get('position', [0.5, 0.5, 0.5])),
                    default_species=s.get('default_species', 'empty'),
                    site_type=s.get('site_type')))
    else:
        layer = pt.add_layer(name='default')
        layer.sites.append(Site(name='site', pos=(0.5, 0.5, 0.5),
                                default_species='empty'))

    # Parameters
    if 'temperature' in data:
        pt.add_parameter(name='T', value=float(data['temperature']))
    if 'potential' in data:
        pt.add_parameter(name='U', value=float(data['potential']))
    for pname, val in data.get('parameters', {}).items():
        if pname == 'T' and 'temperature' in data:
            continue
        if pname == 'U' and 'potential' in data:
            continue
        if isinstance(val, dict):
            pt.add_parameter(name=pname, value=val.get('value', 0.0),
                             adjustable=True, min=val.get('min', 0.0),
                             max=val.get('max', 0.0))
        else:
            pt.add_parameter(name=pname, value=val)

    # Reactions
    reactions = data.get('reactions', [])
    if isinstance(reactions, dict):
        for rname, rc in reactions.items():
            rc['name'] = rname
            _add_reaction_yaml(pt, rc, source)
    elif isinstance(reactions, list):
        for rxn in reactions:
            _add_reaction_yaml(pt, rxn, source)

    # Lateral
    lateral = data.get('lateral', data.get('lateral_interactions', []))
    if isinstance(lateral, dict):
        for pair, energy in lateral.items():
            sp1, sp2 = pair.split('-', 1)
            pt.add_lateral_interaction(sp1.strip(), sp2.strip(), energy)
    elif isinstance(lateral, list):
        for lat in lateral:
            sp = lat['species']
            pt.add_lateral_interaction(sp[0], sp[1], lat['energy'])

    # BEP
    beps = data.get('bep', data.get('bep_relations', []))
    if isinstance(beps, dict):
        for proc, alpha in beps.items():
            pt.add_bep_relation(process_name=proc, alpha=alpha)
    elif isinstance(beps, list):
        for b in beps:
            pt.add_bep_relation(process_name=b['process'],
                                alpha=b.get('alpha', 0.5))

    return pt


def _add_reaction_yaml(pt, rxn, source):
    """Parse a YAML reaction entry."""
    name = rxn['name']
    rxn_type = rxn.get('type', 'thermal')

    if rxn_type == 'diffusion':
        Ea = rxn.get('Ea', 0.0)
        rate_expr = f'kB*T/h*exp(-{Ea}*eV/(kB*T))'
        pt.add_diffusion(species=rxn['species'], rate_constant=rate_expr,
                         name_prefix=f'{rxn["species"]}_diff')
        return

    reaction_str = rxn.get('reaction', '')
    if not reaction_str:
        raise ValueError(f"{source}: reaction '{name}' missing 'reaction'")

    Ea = rxn.get('Ea', 0.0)
    rate_expr = rxn.get('rate', _make_rate_expr(Ea, rxn_type, rxn))
    conditions, actions = _parse_reaction_string(pt, reaction_str)
    pt.add_process(name=name, conditions=conditions, actions=actions,
                   rate_constant=str(rate_expr),
                   tof_count=rxn.get('tof_count', {name: 1}),
                   reverse_of=rxn.get('reverse_of'),
                   site_type=rxn.get('site_type'))


# ============================================================================
#  JSON format (Project.save() compatibility)
# ============================================================================

def _load_json(path):
    """Load from JSON."""
    import json
    with open(path, 'r') as f:
        data = json.load(f)

    pt = Project()
    pt.set_meta(**data.get('meta', {}))
    for sp in data.get('species', []):
        pt.add_species(name=sp['name'], color=sp.get('color'),
                       tags=sp.get('tags', ''))
    for p in data.get('parameters', []):
        pt.add_parameter(name=p['name'], value=p['value'],
                         adjustable=p.get('adjustable', False),
                         min=p.get('min', 0.0), max=p.get('max', 0.0))
    lat = data.get('lattice', {})
    cell = lat.get('cell', [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    pt.lattice.cell = np.array(cell, dtype=float)
    for l in lat.get('layers', []):
        layer = pt.add_layer(name=l['name'])
        for s in l.get('sites', []):
            layer.sites.append(Site(
                name=s['name'], pos=tuple(s.get('pos', [0.5, 0.5, 0.5])),
                default_species=s.get('default_species', 'empty'),
                site_type=s.get('site_type')))
    layer_name = pt.lattice.layers[0].name if pt.lattice.layers else 'default'
    site_name = (pt.lattice.layers[0].sites[0].name
                 if pt.lattice.layers and pt.lattice.layers[0].sites else 'site')
    for proc in data.get('processes', []):
        conditions, actions = [], []
        for c in proc.get('conditions', []):
            coord = Coord(offset=tuple(c['offset']), site=site_name, layer=layer_name)
            conditions.append(Condition(coord, c['species']))
        for a in proc.get('actions', []):
            coord = Coord(offset=tuple(a['offset']), site=site_name, layer=layer_name)
            actions.append(Action(coord, a['species']))
        pt.add_process(name=proc['name'], conditions=conditions, actions=actions,
                       rate_constant=proc.get('rate_constant', '0'),
                       tof_count=proc.get('tof_count', {}),
                       reverse_of=proc.get('reverse_of'),
                       site_type=proc.get('site_type'))
    for li in data.get('lateral_interactions', []):
        pt.add_lateral_interaction(li['species1'], li['species2'], li['energy'])
    for b in data.get('bep_relations', []):
        pt.add_bep_relation(b['process_name'], b.get('alpha', 0.5))
    return pt


# ============================================================================
#  Shared helpers
# ============================================================================

def _make_rate_expr(Ea, rxn_type, rxn):
    """Generate rate expression from Ea and type."""
    if rxn_type == 'thermal':
        return f'kB*T/h*exp(-{Ea}*eV/(kB*T))'
    elif rxn_type in ('electrochemical', 'pcet', 'PCET'):
        beta = rxn.get('beta', 'beta_BV')
        if isinstance(beta, (int, float)):
            if beta >= 0:
                return f'kB*T/h*exp(-({Ea}+{beta}*U)*eV/(kB*T))'
            else:
                return f'kB*T/h*exp(-({Ea}{beta}*U)*eV/(kB*T))'
        return f'kB*T/h*exp(-({Ea}+{beta}*U)*eV/(kB*T))'
    elif rxn_type in ('adsorption', 'hertz-knudsen'):
        p = rxn.get('pressure', 'p')
        m = rxn.get('mass', 'm')
        area = rxn.get('area', 'A_site')
        return f'{p}*bar*{area}/sqrt(2*pi*{m}*umass/(kB*T))'
    else:
        raise ValueError(f"Unknown reaction type: {rxn_type}")


def _parse_reaction_string(pt, reaction_str):
    """Parse 'A + B -> C + D'. * = empty."""
    if '->' not in reaction_str:
        raise ValueError(f"Reaction must contain '->': {reaction_str}")
    lhs, rhs = reaction_str.split('->')
    lhs_tokens = [t.strip() for t in lhs.split('+')]
    rhs_tokens = [t.strip() for t in rhs.split('+')]
    if len(lhs_tokens) != len(rhs_tokens):
        raise ValueError(f"Mismatch in '{reaction_str}'")
    conditions, actions = [], []
    for i, (lt, rt) in enumerate(zip(lhs_tokens, rhs_tokens)):
        ls, lo, lsite = _parse_token(pt, lt, i)
        rs, ro, rsite = _parse_token(pt, rt, i)
        coord = Coord(offset=lo, site=lsite or rsite,
                      layer=pt.lattice.layers[0].name)
        conditions.append(Condition(coord, ls))
        actions.append(Action(coord, rs))
    return conditions, actions


def _parse_token(pt, token, position):
    """Parse species token. * = empty. (nn)/(adj) = neighbor."""
    token = token.strip()
    offset = (0, 0, 0)
    site_name = None
    for tag in ('(neighbor)', '(adj)', '(nn)', '(NN)'):
        if token.endswith(tag):
            token = token[:-len(tag)].strip()
            offset = (1, 0, 0)
            break
    else:
        if position > 0:
            offset = (1, 0, 0)
    if '@' in token:
        species, site_spec = token.split('@', 1)
        species = species.strip()
        if '.(' in site_spec:
            site_name, off_str = site_spec.split('.(', 1)
            offset = eval('(' + off_str)
        else:
            site_name = site_spec.strip()
    else:
        species = token
    if species == '*':
        species = 'empty'
    if site_name is None and pt.lattice.layers:
        layer = pt.lattice.layers[0]
        if layer.sites:
            site_name = layer.sites[0].name
    return species, offset, site_name


# ============================================================================
#  Export
# ============================================================================

def project_to_sparkin(project, path=None):
    """Export a Project to .sparkin (VASP-style) format."""
    lines = []
    name = project.meta.get('model_name', 'unnamed')
    lines.append(f'SYSTEM    = {name}')

    cell = project.lattice.cell
    if isinstance(cell, np.ndarray) and np.allclose(cell, np.diag(np.diag(cell))):
        diag = [float(cell[i, i]) for i in range(3)]
        matched = None
        for surf, (a, vac, _) in SURFACE_PRESETS.items():
            if abs(diag[0] - a) < 0.01 and abs(diag[1] - a) < 0.01:
                matched = surf
                break
        if matched:
            lines.append(f'SURFACE   = {matched}')
        else:
            lines.append(f'LATTICE   = {diag[0]:.4f} {diag[1]:.4f} {diag[2]:.4f}')

    for p in project.parameter_list:
        if p.name == 'T':
            lines.append(f'TEMP      = {p.value}')
        elif p.name == 'U':
            lines.append(f'POTENTIAL = {p.value}')

    sp = [s.name for s in project.species_list if s.name != 'empty']
    if sp:
        lines.append(f'SPECIES   = {" ".join(sp)}')

    for p in project.parameter_list:
        if p.name not in ('T', 'U'):
            lines.append(f'PARAMETER = {p.name}  {p.value}')

    lines.append('')
    for proc in project.process_list:
        conds = ['*' if c.species == 'empty' else c.species for c in proc.conditions]
        acts = ['*' if a.species == 'empty' else a.species for a in proc.actions]
        rxn_str = ' + '.join(conds) + ' -> ' + ' + '.join(acts)
        lines.append(f'REACTION  = {proc.name} : {rxn_str} : rate={proc.rate_constant}')

    lines.append('')
    for li in project.lateral_interactions:
        lines.append(f'LATERAL   = {li.species1}-{li.species2} : {li.energy}')
    for pname, bep in project.bep_relations.items():
        lines.append(f'BEP       = {pname} : {bep.alpha}')

    text = '\n'.join(lines) + '\n'
    if path:
        with open(path, 'w') as f:
            f.write(text)
    else:
        return text


def project_to_yaml(project, path=None):
    """Export a Project to YAML format."""
    import yaml
    data = {'name': project.meta.get('model_name', 'unnamed')}
    sp_list = [s.name for s in project.species_list if s.name != 'empty']
    data['species'] = sp_list

    cell = project.lattice.cell
    if isinstance(cell, np.ndarray) and np.allclose(cell, np.diag(np.diag(cell))):
        diag = [float(cell[i, i]) for i in range(3)]
        matched = None
        for surf, (a, vac, _) in SURFACE_PRESETS.items():
            if abs(diag[0] - a) < 0.01 and abs(diag[1] - a) < 0.01:
                matched = surf
                break
        data['surface'] = matched if matched else None
        if not matched:
            data['lattice'] = {'cell': diag}

    for p in project.parameter_list:
        if p.name == 'T':
            data['temperature'] = p.value
        elif p.name == 'U':
            data['potential'] = p.value

    rxns = {}
    for proc in project.process_list:
        conds = ['*' if c.species == 'empty' else c.species for c in proc.conditions]
        acts = ['*' if a.species == 'empty' else a.species for a in proc.actions]
        rxns[proc.name] = {
            'reaction': ' + '.join(conds) + ' -> ' + ' + '.join(acts),
            'rate': proc.rate_constant,
        }
    data['reactions'] = rxns

    if project.lateral_interactions:
        data['lateral'] = {f'{li.species1}-{li.species2}': li.energy
                           for li in project.lateral_interactions}
    if project.bep_relations:
        data['bep'] = {n: b.alpha for n, b in project.bep_relations.items()}

    text = yaml.dump(data, default_flow_style=False, sort_keys=False)
    if path:
        with open(path, 'w') as f:
            f.write(text)
    else:
        return text


# Legacy aliases
load_spark = load_sparkin
loads_spark = loads_sparkin
project_to_spark = project_to_sparkin
