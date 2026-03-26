# SPARK

**SPatial Atomistic Reaction Kinetics**

A Python + Rust toolkit for lattice kinetic Monte Carlo (KMC) and mean-field microkinetic modeling of electrochemical and heterogeneous catalysis.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Rust](https://img.shields.io/badge/rust-1.70+-orange.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## Overview

SPARK provides two complementary engines for simulating surface catalytic reactions:

- **`spark`** (Python) -- research-friendly API with SciPy solvers, DFT energy interpolation, and polarization curve support
- **`spark-rs`** (Rust) -- high-performance binary with Fenwick tree O(log N) site selection and Newton-Raphson steady-state solver

Both engines implement the BKL rejection-free KMC algorithm with pairwise lateral interactions, BEP relations, Butler-Volmer electrochemical rates, and kmos-compatible rate expressions.

### Validated

SPARK has been quantitatively validated against two established KMC codes on an identical HER/Pt(111) test case:

| Engine | Language | Coverage (theta_H) | Current density (j) | Wall time |
|--------|----------|:------------------:|:-------------------:|:---------:|
| **SPARK** | Python | 0.338 | reference | 281 s |
| **kmos** | Python/Fortran | 0.332 | < 2% deviation | 20 s |
| **Zacros 4.0** | Fortran | 0.341 | < 2% deviation | 12 s |

> All three engines agree within stochastic noise across 11 potential points (0 to -0.5 V vs RHE). See [`benchmark/BENCHMARK_SUMMARY.md`](benchmark/BENCHMARK_SUMMARY.md) for full results.

---

## Features

| | Feature | Description |
|---|---------|-------------|
| **KMC** | BKL rejection-free algorithm | O(1) bookkeeping (Python), Fenwick tree O(log N) site selection (Rust) |
| **Spatial** | Neighbor list + lateral interactions | 4-NN (2D) / 6-NN (3D) with PBC, pairwise repulsion/attraction |
| **Electrochemistry** | Butler-Volmer PCET rates | Potential-dependent barriers, BEP relations |
| **Mean-field** | ODE steady-state solver | SciPy fsolve (Python), Newton-Raphson with damped line search (Rust) |
| **Polarization** | j-U curves from DFT data | Cubic spline interpolation of constant-potential DFT barriers |
| **Rate parser** | kmos-compatible expressions | `kB*T/h*exp(-(Ea + beta_BV*U)*eV/(kB*T))` |
| **I/O** | SPARKIN + YAML + JSON + Python API | VASP-style input files, YAML, JSON, or programmatic |

---

## Installation

### Python

Requires Python >= 3.8, NumPy, SciPy, PyYAML (`pip install pyyaml`).

```bash
git clone https://github.com/WanluLigroupUCSD/SPARK.git
cd SPARK
pip install numpy scipy matplotlib pyyaml
```

No `pip install` step -- import directly:

```python
from spark import load_model, KMCEngine, MicroKineticModel
```

### Rust (optional)

```bash
cd spark-rs
cargo build --release
# Binary: spark-rs/target/release/spark
```

---

## Quick Start

### Option A: Input file (recommended)

VASP-style `TAG = VALUE` format. One file, self-documenting:

```
# her_model.sparkin

SYSTEM    = HER_Pt111
SURFACE   = Pt(111)             (auto-sets lattice constant)
TEMP      = 298                 (K)
POTENTIAL = -0.2                (V vs RHE)
SPECIES   = H                   (* is always implicit)

REACTION  = Volmer_fwd : * -> H : 0.67 : electrochemical : 0.5
REACTION  = Volmer_rev : H -> * : 0.62 : electrochemical : -0.5
REACTION  = Tafel : H + H -> * + * : 0.85
REACTION  = Heyrovsky : H -> * : 0.70 : electrochemical : 0.5
DIFFUSION = H : 0.10

LATERAL   = H-H : 0.10         (eV, positive = repulsive)
BEP       = Volmer_fwd : 0.5

LSIZE     = 20 20
NEQUIL    = 100000
NSTEPS    = 100000
```

Load and run:

```python
from spark import load_model, KMCEngine

pt = load_model('her_model.sparkin')
size = pt.meta.get('lsize', [20, 20])
engine = KMCEngine(pt, size=size, print_rates=True)
engine.do_steps(100000)
print(engine.get_coverage())
```

### Option B: Python API

```python
from spark import Project, Site, Condition, Action
import numpy as np

pt = Project()
pt.set_meta(model_name='CO_on_Pd100', model_dimension=2)
pt.add_species(name='empty', color='#ffffff')
pt.add_species(name='CO', color='#ff0000')

layer = pt.add_layer(name='surface')
layer.sites.append(Site(name='hollow', pos=(0.5, 0.5, 0.5), default_species='empty'))
pt.lattice.cell = np.diag([3.5, 3.5, 10.0])

pt.add_parameter(name='T', value=600.0)
pt.add_parameter(name='p_CO', value=1.0)
pt.add_parameter(name='A', value='(3.5*angstrom)**2')
pt.add_parameter(name='deltaG', value=-0.5)

coord = pt.lattice.generate_coord('hollow')
pt.add_process(
    name='CO_adsorption',
    rate_constant='p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)',
    conditions=[Condition(coord=coord, species='empty')],
    actions=[Action(coord=coord, species='CO')],
    tof_count={'CO_adsorption': 1},
)
pt.add_process(
    name='CO_desorption',
    rate_constant='kB*T/h*exp(-(-deltaG)*eV/(kB*T))',
    conditions=[Condition(coord=coord, species='CO')],
    actions=[Action(coord=coord, species='empty')],
    tof_count={'CO_desorption': 1},
)
```

### Running KMC

```python
from spark import KMCEngine

engine = KMCEngine(pt, size=[20, 20], print_rates=True)
engine.do_steps(100000)           # equilibration
engine.get_tof()                  # reset TOF baseline
engine.do_steps(100000)           # production

print(engine.get_coverage())      # {'CO': 0.498, 'empty': 0.502}
print(engine.get_tof())           # {'CO_adsorption': 1.23e6, ...}
```

### Mean-field microkinetics

```python
from spark import MicroKineticModel
from spark.rates import tst_rate

mkm = MicroKineticModel()
mkm.add_species('N2')
mkm.add_species('NNH')
mkm.add_reaction(
    name='N2_hydrogenation',
    reactants={'N2': 1}, products={'NNH': 1},
    rate_fwd=lambda p: tst_rate(2.82, p['T']),
    rate_rev=lambda p: tst_rate(1.50, p['T']),
)
mkm.parameters = {'T': 300, 'U': -0.5}
ss = mkm.solve_steady_state()
mkm.print_summary(ss)
```

### Electrochemical features

```python
# Lateral interactions
pt.add_lateral_interaction('H', 'H', energy=0.10)  # repulsive NN

# BEP relation
pt.add_bep_relation('Volmer_fwd', alpha=0.5)

# Polarization curve
from spark.polarization import PolarizationCurve
pc = PolarizationCurve(mkm, n_electrons={'H2_production': 2}, A_site=(2.775e-10)**2)
results = pc.compute(U_range=np.linspace(-0.5, 0, 11), T=298)
```

---

## Input Format Reference

### SPARKIN format (VASP-style)

`TAG = VALUE` format with `(...)` inline comments. `#` for line comments.

**Tags:**

| Tag | Description | Example |
|-----|-------------|---------|
| `SYSTEM` | Model name | `SYSTEM = HER_Pt111` |
| `SURFACE` | Surface preset | `SURFACE = Pt(111)` |
| `LATTICE` | Custom cell (A) | `LATTICE = 3.5 3.5 10.0` |
| `TEMP` | Temperature (K) | `TEMP = 298` |
| `POTENTIAL` | Potential (V vs RHE) | `POTENTIAL = -0.2` |
| `SPECIES` | Surface species | `SPECIES = H CO` |
| `PARAMETER` | Custom parameter | `PARAMETER = p_CO 1.0` |
| `REACTION` | Elementary step | see below |
| `DIFFUSION` | NN hop (auto-expand) | `DIFFUSION = H : 0.10` |
| `LATERAL` | Pair interaction (eV) | `LATERAL = H-H : 0.10` |
| `BEP` | BEP relation | `BEP = Volmer_fwd : 0.5` |
| `LSIZE` | Lattice size | `LSIZE = 20 20` |
| `NEQUIL` | Equilibration steps | `NEQUIL = 100000` |
| `NSTEPS` | Production steps | `NSTEPS = 100000` |

**REACTION format** -- colon-separated fields:

```
REACTION = name : equation : Ea(eV) [: type : beta]
REACTION = name : equation : rate=expression
```

- Thermal (default): `name : * -> H : 0.85`  -->  `kBT/h * exp(-Ea/kBT)`
- Electrochemical: `name : * -> H : 0.67 : electrochemical : 0.5`  -->  `kBT/h * exp(-(Ea+beta*U)/kBT)`
- Custom rate: `name : CO -> * : rate=kB*T/h*exp(-0.5*eV/(kB*T))`

**Reaction string**: `*` = empty/vacant site. `+` separates multi-site participants (2nd site is neighbor by default).

**Surface presets:**

| Preset | a (A) | Site | Preset | a (A) | Site |
|--------|-------|------|--------|-------|------|
| Pt(111) | 2.775 | top | Cu(111) | 2.556 | top |
| Pt(100) | 2.775 | hollow | Cu(100) | 2.556 | hollow |
| Pd(111) | 2.751 | top | Au(111) | 2.884 | top |
| Ni(111) | 2.492 | top | Ru(0001) | 2.706 | top |

### YAML format (alternative)

Also supported for complex models. See [`examples/her_Pt111.yaml`](examples/her_Pt111.yaml).

---

## Benchmark

HER on Pt(111) with identical parameters across all three engines:

- 4 elementary steps (Volmer fwd/rev, Tafel, Heyrovsky)
- Square lattice 20x20, T = 298 K, U = 0 to -0.5 V
- Lateral interactions: H-H repulsion (+0.10 eV)
- BEP on Volmer (alpha = 0.5)

<p align="center">
  <img src="benchmark/fig_benchmark_comparison.png" width="800">
</p>

**Key results:**
- Current density agreement: average relative error ~1.4% (stochastic noise)
- Coverage plateau at theta_H ~ 0.34, governed by H-H repulsion
- Heyrovsky pathway dominates; Tafel slope ~120 mV/dec
- Speed: Zacros (12s) > kmos (20s) > SPARK-Python (281s)

Full data tables: [`benchmark/BENCHMARK_SUMMARY.md`](benchmark/BENCHMARK_SUMMARY.md)

---

## Project Structure

```
SPARK/
├── spark/                    # Python package
│   ├── io.py                 #   YAML/JSON model loader
│   ├── engine.py             #   BKL KMC engine
│   ├── microkinetic.py       #   Mean-field ODE solver
│   ├── polarization.py       #   Polarization curves & DFT interpolation
│   ├── rates.py              #   Rate expression parser
│   ├── types.py              #   Data model (Project, Species, Process, ...)
│   ├── analysis.py           #   Trajectory recording, TOF, steady-state
│   └── units.py              #   Physical constants
├── spark-rs/                 # Rust implementation
│   └── src/
│       ├── engine.rs         #   BKL + Fenwick tree
│       ├── microkinetic.rs   #   Newton-Raphson solver
│       ├── polarization.rs   #   Cubic spline interpolation
│       └── ...
├── examples/                 # Example input files
│   ├── her_Pt111.sparkin     #   HER on Pt(111) -- VASP-style format
│   ├── co_adsorption.sparkin #   CO/Pd(100) -- VASP-style format
│   ├── her_Pt111.yaml        #   HER on Pt(111) -- YAML alternative
│   └── co_adsorption.yaml    #   CO/Pd(100) -- YAML alternative
├── models/                   # Pre-built reaction models (Python)
│   ├── her_Pt111.py          #   HER on Pt(111)
│   ├── n2_reduction_Mo.py    #   N2 reduction on Mo
│   └── co_adsorption_test.py #   CO/Pd(100) test
├── benchmark/                # 3-engine benchmark (SPARK vs kmos vs Zacros)
├── tutorial.py               # 7-part interactive tutorial
├── validate_her.py           # HER validation suite
└── run_simulation.py         # CLI simulation driver
```

---

## Rust CLI Reference

```bash
spark kmc   --t 300 --u -1.0 --size 20 --steps 1000000   # Lattice KMC
spark mkm   --t 300 --u -1.0                               # Mean-field steady state
spark scan-u --t 300 --u-min -0.5 --u-max -3.0             # Potential scan
spark scan-t --u -1.0 --t-min 250 --t-max 500              # Temperature scan
spark polarization --dft-data dft.json --t 300              # Polarization curve
spark validate                                              # Langmuir isotherm test
```

All commands support `--cycle` flag for pure catalytic cycle mode (no adsorption/desorption).

---

## Algorithm

### KMC: BKL Rejection-Free (VSSM / n-fold way)

Each step:
1. Compute cumulative rates R_i
2. Draw random numbers r1, r2, r3
3. Advance time: dt = -ln(r1) / R_total
4. Select process (binary search on r2 * R_total)
5. Select site (uniform from available sites)
6. Execute and update bookkeeping

Rust optimizations: Fenwick tree site selection, flat lateral energy arrays, zero-allocation coordinates, periodic rebuild to prevent float drift.

### Mean-Field Solver

Solves dtheta_i/dt = sum_j(nu_ij * r_j) to steady state.

- **Python**: `scipy.optimize.fsolve`
- **Rust**: Euler pre-equilibration (dt: 1e-15 -> 1s) + Newton-Raphson with damped line search; typically 5-10 iterations

---

## Performance

| Benchmark | Rust | Python |
|-----------|------|--------|
| CO Langmuir validation (5 x 700k steps) | 3.7 s | ~80 s |
| HER benchmark (11 potentials, 100k steps each) | -- | 281 s |
| Potential scan (31 points, MKM) | 0.6 s | ~2 s |
| Binary size | 1.6 MB | -- |

---

## Documentation

- **[`docs/tutorial.md`](docs/tutorial.md)** -- Complete tutorial: SPARKIN format, KMC, MKM, polarization, advanced features
- [`examples/`](examples/) -- Example input files (.sparkin, .yaml)
- [`tutorial.py`](tutorial.py) -- Interactive Python tutorial (7 parts)
- [`benchmark/BENCHMARK_SUMMARY.md`](benchmark/BENCHMARK_SUMMARY.md) -- Full 3-engine benchmark data
- [`docs/`](docs/) -- DFT methodology, research plans, software comparison

---

## Citation

```bibtex
@software{spark_kmc,
  title  = {SPARK: SPatial Atomistic Reaction Kinetics},
  author = {Wanlu Li Group},
  url    = {https://github.com/WanluLigroupUCSD/SPARK},
  year   = {2026}
}
```

---

## License

MIT License

## Acknowledgments

- KMC algorithm inspired by [kmos](https://github.com/mhoffman/kmos) (M. Hoffmann et al.)
- Benchmark validated against [Zacros](https://zacros.org/) (M. Stamatakis)
- Developed at the [Wanlu Li Group](https://wanluligroup.ucsd.edu/), UC San Diego
