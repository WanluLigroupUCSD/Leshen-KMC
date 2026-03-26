#!/usr/bin/env python3
"""
Benchmark: HER on Pt(111) — SPARK vs kmos vs Zacros

Runs the identical HER model on all three KMC engines and compares:
  - H coverage (θ_H) vs potential
  - Current density (j) vs potential
  - Polarization curves

DFT parameters (from Li et al. 2024):
  Ea_Volmer_fwd  = 0.67 eV
  Ea_Volmer_rev  = 0.62 eV
  Ea_Tafel       = 0.85 eV
  Ea_Heyrovsky   = 0.70 eV
  Ea_diffusion   = 0.10 eV   (not used for speed)
  eps_HH         = +0.10 eV  (repulsive NN)
  alpha_BEP      = 0.5       (Volmer only)
  beta_BV        = 0.5       (Butler-Volmer symmetry factor)

Lattice: square (consistent across all 3 engines), a = 2.775 Å, 20×20
T = 298 K, U = 0 to -0.5 V vs RHE
"""

import os
import sys
import time
import shutil
import subprocess
import numpy as np

# Add SPARK to path
PROJ_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJ_DIR)

# ============================================================================
#  Constants & Parameters
# ============================================================================
kB = 1.380649e-23       # J/K
h_planck = 6.62607015e-34  # J·s
eV_J = 1.602176634e-19  # J/eV
e_charge = eV_J         # C

T = 298.0
A_PT = 2.775  # Angstrom
A_SITE = (A_PT * 1e-10) ** 2  # m² per site

EA_VOLMER_FWD = 0.67
EA_VOLMER_REV = 0.62
EA_TAFEL = 0.85
EA_HEYROVSKY = 0.70
EPS_HH = 0.10
ALPHA_BEP = 0.5
BETA_BV = 0.5

LATTICE_SIZE = 20
EQUIL_STEPS = 50000
PROD_STEPS = 50000

U_POINTS = np.array([0.0, -0.05, -0.10, -0.15, -0.20, -0.25, -0.30,
                      -0.35, -0.40, -0.45, -0.50])

BENCHMARK_DIR = os.path.dirname(os.path.abspath(__file__))
ZACROS_EXE = '/scratch/reny0b/zls/zacros_4.0/build/zacros.x'


def tof_to_j(tof, n_electrons=2):
    """Convert TOF (events/s/site) to current density (mA/cm²)."""
    j = tof * n_electrons * e_charge / A_SITE  # A/m²
    return j * 0.1  # mA/cm²


def compute_rate(Ea_eff, T=298.0):
    """Compute TST rate constant k = kB*T/h * exp(-Ea*eV/(kB*T))."""
    return (kB * T / h_planck) * np.exp(-Ea_eff * eV_J / (kB * T))


# ============================================================================
#  1. SPARK
# ============================================================================
def run_spark():
    """Run HER benchmark using SPARK engine."""
    print("\n" + "=" * 60)
    print("  SPARK")
    print("=" * 60)

    from spark.types import Project as LProject, Site, Condition, Action, Coord
    from spark.engine import KMCEngine

    results = {
        'U': U_POINTS.copy(),
        'theta_H': np.zeros(len(U_POINTS)),
        'j_total': np.zeros(len(U_POINTS)),
    }

    for i, U in enumerate(U_POINTS):
        pt = LProject()
        pt.set_meta(model_name='HER_Pt111', model_dimension=2)
        pt.add_species(name='empty', color='#ffffff')
        pt.add_species(name='H', color='#ff4444')

        layer = pt.add_layer(name='Pt111')
        layer.sites.append(
            Site(name='top', pos=(0.5, 0.5, 0.5), default_species='empty'))
        pt.lattice.cell = np.diag([A_PT, A_PT, 15.0])

        pt.add_parameter(name='T', value=T)
        pt.add_parameter(name='U', value=U)
        pt.add_parameter(name='beta_BV', value=BETA_BV)
        pt.add_parameter(name='Ea_Volmer_fwd', value=EA_VOLMER_FWD)
        pt.add_parameter(name='Ea_Volmer_rev', value=EA_VOLMER_REV)
        pt.add_parameter(name='Ea_Tafel', value=EA_TAFEL)
        pt.add_parameter(name='Ea_Heyrovsky', value=EA_HEYROVSKY)

        coord = pt.lattice.generate_coord('top')
        coord_nn = Coord(offset=(1, 0, 0), site='top', layer='Pt111')

        pt.add_process(
            name='Volmer_fwd',
            conditions=[Condition(coord, 'empty')],
            actions=[Action(coord, 'H')],
            rate_constant='kB*T/h*exp(-(Ea_Volmer_fwd + beta_BV*U)*eV/(kB*T))',
            tof_count={'Volmer_fwd': 1},
        )
        pt.add_process(
            name='Volmer_rev',
            conditions=[Condition(coord, 'H')],
            actions=[Action(coord, 'empty')],
            rate_constant='kB*T/h*exp(-(Ea_Volmer_rev - (1-beta_BV)*U)*eV/(kB*T))',
            tof_count={'Volmer_rev': 1},
        )
        pt.add_process(
            name='Tafel',
            conditions=[Condition(coord, 'H'), Condition(coord_nn, 'H')],
            actions=[Action(coord, 'empty'), Action(coord_nn, 'empty')],
            rate_constant='kB*T/h*exp(-Ea_Tafel*eV/(kB*T))',
            tof_count={'H2_Tafel': 1},
        )
        pt.add_process(
            name='Heyrovsky_fwd',
            conditions=[Condition(coord, 'H')],
            actions=[Action(coord, 'empty')],
            rate_constant='kB*T/h*exp(-(Ea_Heyrovsky + beta_BV*U)*eV/(kB*T))',
            tof_count={'H2_Heyrovsky': 1},
        )

        pt.add_lateral_interaction('H', 'H', energy=EPS_HH)
        pt.add_bep_relation('Volmer_fwd', alpha=ALPHA_BEP)
        pt.add_bep_relation('Volmer_rev', alpha=1.0 - ALPHA_BEP)

        engine = KMCEngine(pt, size=[LATTICE_SIZE, LATTICE_SIZE],
                           print_rates=False, banner=False)
        engine.do_steps(EQUIL_STEPS)
        engine.get_tof()  # reset TOF baseline
        engine.do_steps(PROD_STEPS)

        cov = engine.get_coverage()
        tof = engine.get_tof()
        theta_H = cov.get('H', 0.0)
        tof_total = tof.get('H2_Tafel', 0.0) + tof.get('H2_Heyrovsky', 0.0)

        results['theta_H'][i] = theta_H
        results['j_total'][i] = tof_to_j(tof_total)
        print(f"  U={U:+.3f}V  θ_H={theta_H:.4f}  j={results['j_total'][i]:.4e} mA/cm²")

    return results


# ============================================================================
#  2. kmos
# ============================================================================
def build_kmos_model():
    """Build kmos HER model with pairwise lateral interactions."""
    from kmos.types import Project, Site, Layer, Condition, Action, Coord
    from itertools import product as iterproduct

    pt = Project()
    pt.set_meta(author='Benchmark',
                email='benchmark@test.com',
                model_name='HER_Pt111_kmos',
                model_dimension=2,
                debug=0)

    pt.add_species(name='empty', color='#ffffff')
    pt.add_species(name='H', color='#ff4444',
                   representation="Atoms('H')")
    pt.species_list.default_species = 'empty'

    layer = pt.add_layer(name='sq')
    layer.add_site(name='top', pos=[0.5, 0.5, 0.5])

    # Square unit cell (consistent with SPARK's 4-NN square lattice)
    a = A_PT
    pt.lattice.cell = np.diag([a, a, 15.0])

    # Parameters
    pt.add_parameter(name='T', value=T, adjustable=True, min=200, max=500)
    pt.add_parameter(name='U', value=0.0, adjustable=True, min=-1.0, max=0.1)
    pt.add_parameter(name='Ea_Vf', value=EA_VOLMER_FWD)
    pt.add_parameter(name='Ea_Vr', value=EA_VOLMER_REV)
    pt.add_parameter(name='Ea_T', value=EA_TAFEL)
    pt.add_parameter(name='Ea_H', value=EA_HEYROVSKY)
    pt.add_parameter(name='bBV', value=BETA_BV)
    pt.add_parameter(name='eps', value=EPS_HH)
    pt.add_parameter(name='aBEP', value=ALPHA_BEP)

    center = pt.lattice.generate_coord('top')

    # Get 4 nearest neighbors for square lattice
    nn_coords = [
        pt.lattice.generate_coord('top.(1,0,0).sq'),
        pt.lattice.generate_coord('top.(-1,0,0).sq'),
        pt.lattice.generate_coord('top.(0,1,0).sq'),
        pt.lattice.generate_coord('top.(0,-1,0).sq'),
    ]
    n_nn = len(nn_coords)
    print(f"  kmos: Found {n_nn} nearest neighbors (square lattice)")

    # ---------- Volmer_fwd with lateral+BEP (enumerate 2^6=64 configs) ----------
    for idx, nn_config in enumerate(iterproduct(['empty', 'H'], repeat=n_nn)):
        n_H = nn_config.count('H')

        # BEP: Ea_eff = Ea_Vf + bBV*U + aBEP*n_H*eps
        rate = (f'(beta*h)**(-1)*exp(-beta*'
                f'(Ea_Vf + bBV*U + aBEP*{n_H}*eps)*eV)')

        conditions = [Condition(coord=center, species='empty')]
        conditions += [Condition(coord=c, species=s)
                       for c, s in zip(nn_coords, nn_config)]
        actions = [Action(coord=center, species='H')]

        pt.add_process(name=f'volmer_fwd_{idx}',
                       conditions=conditions,
                       actions=actions,
                       rate_constant=rate,
                       tof_count={'Volmer_fwd': 1})

    # ---------- Volmer_rev with lateral+BEP ----------
    for idx, nn_config in enumerate(iterproduct(['empty', 'H'], repeat=n_nn)):
        n_H = nn_config.count('H')

        # BEP: delta_H = -n_H*eps, Ea_eff = Ea_Vr + (1-aBEP)*(-n_H*eps)
        #    = Ea_Vr - (1-aBEP)*n_H*eps
        # With -(1-bBV)*U for reverse PCET
        rate = (f'(beta*h)**(-1)*exp(-beta*'
                f'(Ea_Vr - (1-bBV)*U - (1-aBEP)*{n_H}*eps)*eV)')

        conditions = [Condition(coord=center, species='H')]
        conditions += [Condition(coord=c, species=s)
                       for c, s in zip(nn_coords, nn_config)]
        actions = [Action(coord=center, species='empty')]

        pt.add_process(name=f'volmer_rev_{idx}',
                       conditions=conditions,
                       actions=actions,
                       rate_constant=rate,
                       tof_count={'Volmer_rev': 1})

    # ---------- Heyrovsky with default lateral ----------
    for idx, nn_config in enumerate(iterproduct(['empty', 'H'], repeat=n_nn)):
        n_H = nn_config.count('H')

        # Default lateral: Ea_eff = Ea_H + bBV*U - n_H*eps
        rate = (f'(beta*h)**(-1)*exp(-beta*'
                f'(Ea_H + bBV*U - {n_H}*eps)*eV)')

        conditions = [Condition(coord=center, species='H')]
        conditions += [Condition(coord=c, species=s)
                       for c, s in zip(nn_coords, nn_config)]
        actions = [Action(coord=center, species='empty')]

        pt.add_process(name=f'heyrovsky_{idx}',
                       conditions=conditions,
                       actions=actions,
                       rate_constant=rate,
                       tof_count={'H2_Heyrovsky': 1})

    # ---------- Tafel (no lateral for simplicity, negligible contribution) ----------
    # Need a neighbor coord for 2-site process
    nn_right = pt.lattice.generate_coord('top.(1,0,0).sq')
    pt.add_process(name='tafel',
                   conditions=[Condition(coord=center, species='H'),
                               Condition(coord=nn_right, species='H')],
                   actions=[Action(coord=center, species='empty'),
                            Action(coord=nn_right, species='empty')],
                   rate_constant='(beta*h)**(-1)*exp(-beta*Ea_T*eV)',
                   tof_count={'H2_Tafel': 1})

    return pt


def run_kmos():
    """Build, export, compile, and run kmos HER model."""
    print("\n" + "=" * 60)
    print("  kmos")
    print("=" * 60)

    kmos_dir = os.path.join(BENCHMARK_DIR, 'kmos_HER')
    xml_file = os.path.join(kmos_dir, 'HER_Pt111_kmos.xml')

    # Clean previous build
    if os.path.exists(kmos_dir):
        shutil.rmtree(kmos_dir)
    os.makedirs(kmos_dir, exist_ok=True)

    # 1. Build and save model
    print("  Building kmos model...")
    pt = build_kmos_model()
    pt.filename = xml_file
    pt.save()
    print(f"  Saved: {xml_file}")
    n_proc = len(pt.process_list)
    print(f"  Total processes: {n_proc}")

    # 2. Export source using kmos API
    print("  Exporting kmos source...")
    import builtins
    _orig_input = builtins.input
    builtins.input = lambda *a: 'n'  # suppress JANAF prompts
    try:
        from kmos.io import export_source
        model_dir = os.path.join(kmos_dir, 'HER_Pt111_kmos_local_smart')
        export_source(pt, export_dir=model_dir, code_generator='local_smart')
    except Exception as e:
        print(f"  kmos export error: {e}")
        return None
    finally:
        builtins.input = _orig_input
    print("  Export done.")

    # 3. Manual compilation (fix f2py/numpy compatibility)
    print("  Compiling kmos Fortran code...")
    src_dir = model_dir  # kmos exports files directly into model_dir

    # Fix kind_values.f90 with hardcoded values
    kv_file = os.path.join(src_dir, 'kind_values.f90')
    with open(kv_file, 'w') as f:
        f.write('module kind_values\nimplicit none\n'
                'integer, parameter :: rsingle = 8\n'
                'integer, parameter :: rdouble = 8\n'
                'integer, parameter :: ibyte = 1\n'
                'integer, parameter :: ishort = 2\n'
                'integer, parameter :: iint = 4\n'
                'integer, parameter :: ilong = 8\n'
                'end module kind_values\n')

    # Fix proclist.f90: query actual seed_size from Fortran runtime
    proclist_file = os.path.join(src_dir, 'proclist.f90')
    with open(proclist_file) as f:
        proclist_src = f.read()
    # Replace hardcoded seed init with runtime query
    proclist_src = proclist_src.replace(
        '    allocate(seed_arr(seed_size))\n'
        '    seed = seed_in\n'
        '    seed_arr = seed\n'
        '    call random_seed(seed_size)\n'
        '    call random_seed(put=seed_arr)\n'
        '    deallocate(seed_arr)',
        '    call random_seed(size=seed_size)\n'
        '    allocate(seed_arr(seed_size))\n'
        '    seed = seed_in\n'
        '    seed_arr = seed\n'
        '    call random_seed(put=seed_arr)\n'
        '    deallocate(seed_arr)',
    )
    with open(proclist_file, 'w') as f:
        f.write(proclist_src)

    # Compile Fortran sources
    f90_files = ['kind_values.f90', 'base.f90', 'lattice.f90', 'proclist.f90']
    for f in f90_files:
        ret = subprocess.run(
            ['gfortran', '-cpp', '-c', '-fPIC', '-O2', f, '-o', f.replace('.f90', '.o')],
            cwd=src_dir, capture_output=True, text=True)
        if ret.returncode != 0:
            print(f"  Fortran compile error ({f}): {ret.stderr[:200]}")
            return None

    # Generate f2py wrappers
    ret = subprocess.run(
        [sys.executable, '-m', 'numpy.f2py'] + f90_files +
        ['-m', 'kmc_model', '--overwrite-signature'],
        cwd=src_dir, capture_output=True, text=True)

    # Get include paths
    python_inc = subprocess.check_output(
        [sys.executable, '-c',
         "import sysconfig; print(sysconfig.get_path('include'))"],
        text=True).strip()
    numpy_inc = subprocess.check_output(
        [sys.executable, '-c',
         "import numpy; print(numpy.get_include())"],
        text=True).strip()
    f2py_inc = subprocess.check_output(
        [sys.executable, '-c',
         "import numpy.f2py, os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), 'src'))"],
        text=True).strip()
    pyext = subprocess.check_output(
        [sys.executable, '-c',
         "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"],
        text=True).strip()

    # Compile fortranobject.c (provides PyFortran_Type needed at runtime)
    fortranobject_c = os.path.join(f2py_inc, 'fortranobject.c')
    ret = subprocess.run(
        ['cc', '-fPIC', '-c', f'-I{python_inc}', f'-I{numpy_inc}', f'-I{f2py_inc}',
         fortranobject_c, '-o', 'fortranobject.o'],
        cwd=src_dir, capture_output=True, text=True)
    if ret.returncode != 0:
        print(f"  fortranobject.c compile error: {ret.stderr[:300]}")
        return None

    # Compile C wrapper
    ret = subprocess.run(
        ['cc', '-fPIC', '-c', f'-I{python_inc}', f'-I{numpy_inc}', f'-I{f2py_inc}',
         'kmc_modelmodule.c', '-o', 'kmc_modelmodule.o'],
        cwd=src_dir, capture_output=True, text=True)
    if ret.returncode != 0:
        print(f"  C compile error: {ret.stderr[:300]}")
        return None

    # Compile f2py Fortran wrappers
    ret = subprocess.run(
        ['gfortran', '-fPIC', '-c', 'kmc_model-f2pywrappers2.f90',
         '-o', 'kmc_model-f2pywrappers2.o'],
        cwd=src_dir, capture_output=True, text=True)

    # Link into shared library (include fortranobject.o for PyFortran_Type)
    so_name = f'kmc_model{pyext}'
    ret = subprocess.run(
        ['gfortran', '-shared', '-o', so_name,
         'kmc_modelmodule.o', 'fortranobject.o', 'kmc_model-f2pywrappers2.o',
         'kind_values.o', 'base.o', 'lattice.o', 'proclist.o'],
        cwd=src_dir, capture_output=True, text=True)
    if ret.returncode != 0:
        print(f"  Link error: {ret.stderr[:300]}")
        return None

    so_path = os.path.join(src_dir, so_name)
    if not os.path.exists(so_path):
        print("  kmc_model.so not created")
        return None
    print(f"  Compilation done: {so_name}")

    # 4. Run simulations
    print("  Running KMC scans...")
    old_path = sys.path.copy()
    old_cwd = os.getcwd()
    # kmos.run needs kmc_settings.py and kmc_model.so in the path
    sys.path.insert(0, src_dir)
    os.chdir(src_dir)

    try:
        from kmos.run import KMC_Model

        results = {
            'U': U_POINTS.copy(),
            'theta_H': np.zeros(len(U_POINTS)),
            'j_total': np.zeros(len(U_POINTS)),
        }

        for i, U in enumerate(U_POINTS):
            model = KMC_Model(size=[LATTICE_SIZE, LATTICE_SIZE],
                              print_rates=False, banner=False)
            model.parameters.T = T
            model.parameters.U = U

            model.do_steps(EQUIL_STEPS)
            atoms = model.get_atoms(geometry=False)  # reset TOF
            model.do_steps(PROD_STEPS)
            atoms = model.get_atoms(geometry=False)

            # Extract coverage — use occupation header to find H index
            try:
                occ_header = model.get_occupation_header().split()
                h_idx = next(i for i, name in enumerate(occ_header)
                             if name.startswith('H_'))
                theta_H = atoms.occupation[h_idx, 0]
            except (IndexError, AttributeError, StopIteration):
                theta_H = 0.0

            # Extract TOF — use tof header to find H2 production channels
            try:
                tof_data = atoms.tof_data
                tof_header = model.get_tof_header().split()
                tof_tafel = 0.0
                tof_hey = 0.0
                for ti, name in enumerate(tof_header):
                    if 'H2_Tafel' in name:
                        tof_tafel = tof_data[ti]
                    if 'H2_Heyrovsky' in name:
                        tof_hey = tof_data[ti]
                tof_total = tof_tafel + tof_hey
            except (IndexError, AttributeError):
                tof_total = 0.0

            results['theta_H'][i] = theta_H
            results['j_total'][i] = tof_to_j(tof_total)

            model.deallocate()

            print(f"  U={U:+.3f}V  θ_H={theta_H:.4f}  "
                  f"j={results['j_total'][i]:.4e} mA/cm²")

        return results

    except Exception as e:
        print(f"  kmos run error: {e}")
        import traceback
        traceback.print_exc()
        return None
    finally:
        sys.path = old_path
        os.chdir(old_cwd)


# ============================================================================
#  3. Zacros
# ============================================================================
def write_zacros_inputs(run_dir, U, max_time=None):
    """Write Zacros input files for HER at potential U."""
    os.makedirs(run_dir, exist_ok=True)

    # Pre-compute rate constants (in 1/s)
    prefactor = kB * T / h_planck  # ~ 6.25e12 1/s

    k_Vf = prefactor * np.exp(-(EA_VOLMER_FWD + BETA_BV * U) * eV_J / (kB * T))
    k_Vr = prefactor * np.exp(-(EA_VOLMER_REV - (1 - BETA_BV) * U) * eV_J / (kB * T))
    k_Tafel = prefactor * np.exp(-EA_TAFEL * eV_J / (kB * T))
    k_Hey = prefactor * np.exp(-(EA_HEYROVSKY + BETA_BV * U) * eV_J / (kB * T))

    # For Zacros: rate = pre_expon * exp(-activ_eng / (kB*T))
    # We set activ_eng = 0 and use the full rate as pre_expon
    # BUT: we want lateral interactions via cluster expansion + prox_factor
    # So we keep the structure: pre_expon = kB*T/h, activ_eng = Ea_eff(U)

    Ea_Vf_eff = EA_VOLMER_FWD + BETA_BV * U
    Ea_Vr_eff = EA_VOLMER_REV - (1 - BETA_BV) * U
    Ea_Hey_eff = EA_HEYROVSKY + BETA_BV * U

    # Estimate max_time from rates
    if max_time is None:
        rate_max = max(k_Vf, k_Vr, k_Hey, k_Tafel)
        n_sites = LATTICE_SIZE * LATTICE_SIZE
        total_steps = EQUIL_STEPS + PROD_STEPS
        max_time = total_steps / (rate_max * n_sites * 0.5)
        max_time = max(max_time * 5, 1e-6)  # generous margin

    # --- simulation_input.dat ---
    sim_input = f"""random_seed               12345

temperature               {T}
pressure                  1.0

n_gas_species             1
gas_specs_names           H2
gas_energies              0.000
gas_molec_weights         2.016
gas_molar_fracs           0.000

n_surf_species            1
surf_specs_names          H*
surf_specs_dent           1

snapshots                 on event {max(PROD_STEPS // 10, 1000)}
process_statistics        on event {max(PROD_STEPS // 10, 1000)}
species_numbers           on event {max(PROD_STEPS // 100, 100)}

event_report              off

max_steps                 {EQUIL_STEPS + PROD_STEPS}
max_time                  infinity

wall_time                 600

no_restart

finish
"""
    with open(os.path.join(run_dir, 'simulation_input.dat'), 'w') as f:
        f.write(sim_input)

    # --- lattice_input.dat ---
    # Square/rectangular lattice (consistent with SPARK's 4-NN)
    lattice_input = f"""lattice default_choice

rectangular_periodic {A_PT:.6f} {LATTICE_SIZE} {LATTICE_SIZE}

end_lattice
"""
    with open(os.path.join(run_dir, 'lattice_input.dat'), 'w') as f:
        f.write(lattice_input)

    # --- energetics_input.dat ---
    # H* point energy = 0 (absorbed into Ea)
    # H*-H* 1NN pair = +0.10 eV (repulsive)
    energetics_input = f"""energetics

cluster H_point
  sites 1
  lattice_state
    1 H*  1
  graph_multiplicity 1
  cluster_eng  0.000
end_cluster

cluster HH_pair_1NN
  sites 2
  neighboring 1-2
  lattice_state
    1 H*  1
    2 H*  1
  graph_multiplicity 2
  cluster_eng  {EPS_HH:.4f}
end_cluster

end_energetics
"""
    with open(os.path.join(run_dir, 'energetics_input.dat'), 'w') as f:
        f.write(energetics_input)

    # --- mechanism_input.dat ---
    # All steps are irreversible (fwd and rev defined separately)
    # because forward/reverse have different potential dependence (not detailed balance)
    #
    # prox_factor for BEP:
    #   Volmer_fwd: 0.5 (BEP α = 0.5)
    #   Volmer_rev: 0.5 (BEP α = 0.5, but ΔE is negative → barrier decreases)
    #   Heyrovsky:  1.0 (default lateral = full reactant destabilization)
    #   Tafel:      0.0 (no correction, negligible anyway)

    mechanism_input = f"""mechanism

step Volmer_fwd
  sites 1
  initial
    1 *     1
  final
    1 H*    1
  pre_expon  {prefactor:.6e}
  activ_eng  {max(Ea_Vf_eff, 0.0):.6f}
  prox_factor  {ALPHA_BEP:.4f}
end_step

step Volmer_rev
  sites 1
  initial
    1 H*    1
  final
    1 *     1
  pre_expon  {prefactor:.6e}
  activ_eng  {max(Ea_Vr_eff, 0.0):.6f}
  prox_factor  {1.0 - ALPHA_BEP:.4f}
end_step

step Tafel
  gas_reacs_prods  H2 1
  sites 2
  neighboring 1-2
  initial
    1 H*    1
    2 H*    1
  final
    1 *     1
    2 *     1
  pre_expon  {prefactor:.6e}
  activ_eng  {EA_TAFEL:.6f}
  prox_factor  0.0000
end_step

step Heyrovsky
  gas_reacs_prods  H2 1
  sites 1
  initial
    1 H*    1
  final
    1 *     1
  pre_expon  {prefactor:.6e}
  activ_eng  {max(Ea_Hey_eff, 0.0):.6f}
  prox_factor  1.0000
end_step

end_mechanism
"""
    with open(os.path.join(run_dir, 'mechanism_input.dat'), 'w') as f:
        f.write(mechanism_input)


def parse_zacros_output(run_dir):
    """Parse Zacros output for coverage and event counts.

    Zacros output files use .txt extension:
      specnum_output.txt  — species numbers over time
      procstat_output.txt — process statistics over time
      general_output.txt  — general simulation info
    """
    spec_file = os.path.join(run_dir, 'specnum_output.txt')
    proc_file = os.path.join(run_dir, 'procstat_output.txt')
    gen_file = os.path.join(run_dir, 'general_output.txt')

    theta_H = 0.0
    tof_total = 0.0
    n_sites = LATTICE_SIZE * LATTICE_SIZE

    # Parse specnum_output.txt for coverage
    # Format: Entry  Nevents  Time  Temperature  Energy  H*  H2
    if os.path.exists(spec_file):
        with open(spec_file) as f:
            lines = f.readlines()

        # Find data lines (skip header line with column names)
        data_lines = []
        for line in lines:
            parts = line.split()
            if len(parts) >= 7:
                try:
                    int(parts[0])  # Entry is integer
                    data_lines.append(parts)
                except ValueError:
                    pass

        if data_lines:
            # Use second half for production-phase average
            half = len(data_lines) // 2
            prod_lines = data_lines[half:] if half > 0 else data_lines
            h_counts = []
            for parts in prod_lines:
                try:
                    n_H = int(parts[5])  # H* column (6th, 0-indexed=5)
                    h_counts.append(n_H)
                except (ValueError, IndexError):
                    pass
            if h_counts:
                theta_H = np.mean(h_counts) / n_sites

    # Parse procstat_output.txt for TOF
    # Format is multi-line per configuration entry:
    #   Line 1: "configuration  N  Nevents  Time"
    #   Line 2: per-process rates
    #   Line 3: per-process cumulative counts (Volmer_fwd  Volmer_rev  Tafel  Heyrovsky)
    if os.path.exists(proc_file):
        with open(proc_file) as f:
            content = f.read()

        # Extract all "configuration" blocks
        blocks = content.split('configuration')
        config_data = []
        for block in blocks[1:]:  # skip everything before first "configuration"
            lines = block.strip().split('\n')
            if len(lines) >= 3:
                # First line: entry_num  Nevents  Time
                header_parts = lines[0].split()
                # Third line: cumulative event counts per process
                count_line = lines[2].strip()
                count_parts = count_line.split()
                try:
                    nevents = int(header_parts[1])
                    sim_time = float(header_parts[2])
                    # Counts: Volmer_fwd, Volmer_rev, Tafel, Heyrovsky
                    counts = [int(x) for x in count_parts]
                    config_data.append({
                        'nevents': nevents,
                        'time': sim_time,
                        'counts': counts
                    })
                except (ValueError, IndexError):
                    pass

        if len(config_data) >= 2:
            # Use second half for production TOF
            mid = len(config_data) // 2
            d_start = config_data[mid]
            d_end = config_data[-1]
            dt = d_end['time'] - d_start['time']
            if dt > 0 and len(d_end['counts']) >= 5:
                # Counts order: Overall(0), Volmer_fwd(1), Volmer_rev(2), Tafel(3), Heyrovsky(4)
                n_tafel = d_end['counts'][3] - d_start['counts'][3]
                n_hey = d_end['counts'][4] - d_start['counts'][4]
                tof_total = (n_tafel + n_hey) / (dt * n_sites)

    # Check for errors
    if os.path.exists(gen_file):
        with open(gen_file) as f:
            for line in f:
                if 'error' in line.lower() and 'warning' not in line.lower():
                    print(f"    Zacros error: {line.strip()}")

    return theta_H, tof_total


def run_zacros():
    """Run HER benchmark using Zacros."""
    print("\n" + "=" * 60)
    print("  Zacros")
    print("=" * 60)

    if not os.path.exists(ZACROS_EXE):
        print(f"  Zacros executable not found: {ZACROS_EXE}")
        return None

    zacros_dir = os.path.join(BENCHMARK_DIR, 'zacros_HER')
    if os.path.exists(zacros_dir):
        shutil.rmtree(zacros_dir)
    os.makedirs(zacros_dir, exist_ok=True)

    results = {
        'U': U_POINTS.copy(),
        'theta_H': np.zeros(len(U_POINTS)),
        'j_total': np.zeros(len(U_POINTS)),
    }

    for i, U in enumerate(U_POINTS):
        run_dir = os.path.join(zacros_dir, f'U_{abs(U):.3f}')
        write_zacros_inputs(run_dir, U)

        # Run Zacros
        t0 = time.time()
        ret = subprocess.run(
            [ZACROS_EXE],
            cwd=run_dir, capture_output=True, text=True,
            timeout=300
        )
        dt = time.time() - t0

        if ret.returncode != 0:
            print(f"  U={U:+.3f}V  FAILED ({dt:.1f}s): {ret.stderr[:200]}")
            continue

        theta_H, tof_total = parse_zacros_output(run_dir)
        results['theta_H'][i] = theta_H
        results['j_total'][i] = tof_to_j(tof_total)
        print(f"  U={U:+.3f}V  θ_H={theta_H:.4f}  "
              f"j={results['j_total'][i]:.4e} mA/cm²  ({dt:.1f}s)")

    return results


# ============================================================================
#  4. Plotting
# ============================================================================
def make_comparison_plots(results_dict, outdir=None):
    """Generate comparison plots."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if outdir is None:
        outdir = BENCHMARK_DIR
    os.makedirs(outdir, exist_ok=True)

    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 14,
        'legend.fontsize': 11,
        'figure.dpi': 150,
        'savefig.dpi': 200,
        'savefig.bbox': 'tight',
    })

    colors = {'SPARK': 'C0', 'kmos': 'C1', 'Zacros': 'C2'}
    markers = {'SPARK': 'o', 'kmos': 's', 'Zacros': '^'}

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # (a) Coverage vs U
    ax = axes[0]
    for name, res in results_dict.items():
        if res is not None:
            ax.plot(res['U'], res['theta_H'],
                    f'{markers[name]}-', color=colors[name],
                    ms=7, lw=1.5, label=name)
    ax.set_xlabel('U [V vs RHE]')
    ax.set_ylabel('θ_H')
    ax.set_title('(a) H Coverage')
    ax.set_ylim([-0.02, 0.6])
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (b) Current density vs U (linear scale)
    ax = axes[1]
    for name, res in results_dict.items():
        if res is not None:
            ax.plot(res['U'], res['j_total'],
                    f'{markers[name]}-', color=colors[name],
                    ms=7, lw=1.5, label=name)
    ax.set_xlabel('U [V vs RHE]')
    ax.set_ylabel('j [mA/cm²]')
    ax.set_title('(b) Polarization Curve')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (c) Current density vs U (log scale)
    ax = axes[2]
    for name, res in results_dict.items():
        if res is not None:
            j = np.abs(res['j_total']) + 1e-30
            ax.semilogy(res['U'], j,
                        f'{markers[name]}-', color=colors[name],
                        ms=7, lw=1.5, label=name)
    ax.set_xlabel('U [V vs RHE]')
    ax.set_ylabel('|j| [mA/cm²]')
    ax.set_title('(c) Polarization (log scale)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.suptitle('HER on Pt(111) — KMC Engine Benchmark\n'
                 f'Lattice: {LATTICE_SIZE}×{LATTICE_SIZE} hex, T={T}K, '
                 f'ε(H-H)={EPS_HH}eV, BEP α={ALPHA_BEP}',
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    outfile = os.path.join(outdir, 'fig_benchmark_comparison.png')
    fig.savefig(outfile)
    plt.close(fig)
    print(f"\n  Saved: {outfile}")

    # Print comparison table
    print("\n" + "=" * 80)
    print("  Benchmark Comparison Table")
    print("=" * 80)
    header = f"  {'U [V]':>8s}"
    for name in results_dict:
        if results_dict[name] is not None:
            header += f"  {'θ_'+name:>12s}  {'j_'+name:>12s}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    for i, U in enumerate(U_POINTS):
        row = f"  {U:>+8.3f}"
        for name, res in results_dict.items():
            if res is not None:
                row += f"  {res['theta_H'][i]:>12.4f}  {res['j_total'][i]:>12.4e}"
        print(row)


# ============================================================================
#  Main
# ============================================================================
def main():
    print("=" * 60)
    print("  HER on Pt(111) — KMC Engine Benchmark")
    print(f"  Lattice: {LATTICE_SIZE}×{LATTICE_SIZE}, T={T}K")
    print(f"  Steps: {EQUIL_STEPS} equil + {PROD_STEPS} prod")
    print("=" * 60)

    results = {}

    # 1. SPARK
    t0 = time.time()
    results['SPARK'] = run_spark()
    print(f"  SPARK total time: {time.time() - t0:.1f}s")

    # 2. kmos
    t0 = time.time()
    try:
        results['kmos'] = run_kmos()
    except Exception as e:
        print(f"  kmos failed: {e}")
        import traceback
        traceback.print_exc()
        results['kmos'] = None
    print(f"  kmos total time: {time.time() - t0:.1f}s")

    # 3. Zacros
    t0 = time.time()
    try:
        results['Zacros'] = run_zacros()
    except Exception as e:
        print(f"  Zacros failed: {e}")
        import traceback
        traceback.print_exc()
        results['Zacros'] = None
    print(f"  Zacros total time: {time.time() - t0:.1f}s")

    # 4. Plot comparison
    make_comparison_plots(results)


if __name__ == '__main__':
    main()
