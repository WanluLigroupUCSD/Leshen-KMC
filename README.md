# SPARK — SPatial Atomistic Reaction Kinetics

A unified Python + Rust toolkit for **kinetic Monte Carlo (KMC)** and **mean-field microkinetic modeling (MKM)** of electrochemical and heterogeneous catalysis.

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Rust 1.70+](https://img.shields.io/badge/rust-1.70+-orange.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## Three KMC Engines, One Framework

SPARK provides **three complementary KMC kernels** that share the same Python API and rate-expression parser, plus a Rust-accelerated hot loop:

| Engine | Module | When to use |
|---|---|---|
| **Lattice KMC** | `spark.engine` | Known reaction network on a periodic surface (e.g. HER on Pt(111), CO₂RR on Cu(100)) |
| **Multi-Lattice KMC** | `spark.engine` (multi-layer) | Morphological transitions between coexisting commensurate phases (e.g. PdO ↔ Pd(100), surface oxide reduction) |
| **Off-Lattice OTF KMC** | `spark.offlattice` | Continuous-coordinate systems with on-the-fly mechanism discovery (e.g. amorphous oxides, surface reconstruction) |
| **Dynamic Catalytic KMC** | `spark.dynamic` | Environment-dependent rates with site identity evolving (e.g. PdAu segregation, alloy ordering) |

All four share the same rate-expression parser (kmos-compatible), microkinetic ODE solver, and polarization-curve infrastructure.

---

## Headline Features

### Lattice KMC (spark.engine)

- BKL rejection-free algorithm with O(1) bookkeeping
- Pairwise lateral interactions, BEP relations, Butler-Volmer electrochemical rates
- **Multi-lattice support** (Hoffmann-Reuter-Scheffler 2015): super-lattice with N coexisting `Layer`s, cross-layer "lattice-swap" elementary processes for morphological transitions. See [`docs/multi_lattice_design.md`](docs/multi_lattice_design.md).
- Quantitatively validated against `kmos` and `Zacros 4.0` on HER/Pt(111) (< 2% deviation across 11 potential points)

### Off-Lattice On-the-Fly KMC (spark.offlattice)

- Ported from openFLY (C++): dimer saddle search, 3-step environment matching, mechanism catalogue with symmetry exploitation, Basin / SuperBasin / SuperCache (bac-MRM) acceleration
- **Rust acceleration** via `spark_rs` PyO3 wheel — 3.08× wall-time speedup on Cu(100) slab + 20× saddle-search convergence-rate gain (rotor-sign bug discovered + fixed during Rust port; mirrored to Python)

### Dynamic Catalytic KMC (spark.dynamic)

- DynamicSurface with mutable site identity (graph-based, weak-lattice)
- Environment-dependent rates via `RateEstimator` protocol; swappable backends (lookup table / ML surrogate / GNN)
- Unified catalytic + structural event system with EventCache (100 % hit rate on PdAu CO oxidation benchmark, 50 k steps)

### Microkinetic Modeling (spark.microkinetic)

- Mean-field ODE steady-state solver with continuation-fsolve + LSODA-relax + fsolve polish (handles ill-conditioned Jacobians at near-onset polarization, validated against kMCOS)
- Cubic-spline interpolation of constant-potential DFT barriers
- Tafel-slope analysis, polarization curves, j(U) plotting

---

## Installation

### Quick: Python only (no Rust)

```bash
pip install numpy scipy ase pyyaml
git clone https://github.com/WanluLigroupUCSD/SPARK.git
cd SPARK
pip install -e .
```

`spark` is now importable. All Python-only features work, including off-lattice KMC (slower path).

### Full: Python + Rust acceleration (recommended for off-lattice)

You need `rustc 1.70+` and `cargo`:

```bash
# Install Rust toolchain (Linux/macOS)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Install maturin + clone SPARK
pip install maturin
git clone https://github.com/WanluLigroupUCSD/SPARK.git
cd SPARK

# Build the spark_rs wheel
cd spark-rs
maturin build --release         # produces target/wheels/spark_rs-*.whl
pip install --force-reinstall target/wheels/spark_rs-*.whl
cd ..

# Install the Python package
pip install -e .
```

Verify:

```python
>>> import spark
>>> import spark_rs
>>> spark_rs.hello()
'spark-rs is alive'
>>> spark_rs.version()
'0.4.0'
```

### HPC: cross-account toolchain (KAUST Shaheen)

If you don't have Rust installed but a colleague does, you can use their toolchain:

```bash
# Owner of the rust install:
chmod -R o+rX ~/.cargo ~/.rustup
chmod o+x ~

# Builder:
export RUSTUP_HOME=/home/<rust-owner>/.rustup
export PATH=/home/<rust-owner>/.cargo/bin:$PATH
export CARGO_HOME=/scratch/<your-id>/.cargo    # private build cache
cd /path/to/SPARK/spark-rs
python -m maturin build --release
python -m pip install --force-reinstall target/wheels/spark_rs-*.whl
```

This is what we did for ai-kmc on Shaheen (reny0b's rustc 1.94.0 → wangc0i's gs env wheel).

---

## Quick Start

### Lattice KMC: HER on Pt(111)

```python
from spark.types import Project, Site, Condition, Action
from spark.engine import KMCEngine

pt = Project()
pt.set_meta(model_name="HER_Pt111", model_dimension=2)
pt.add_species(name="empty")
pt.add_species(name="H")

layer = pt.add_layer(name="Pt111")
layer.sites.append(Site(name="top", default_species="empty"))
top = pt.lattice.generate_coord("top")

pt.add_process(
    name="Volmer",
    conditions=[Condition(top, "empty")],
    actions=[Action(top, "H")],
    rate_constant="1e8 * exp(-0.5*beta*eV * (E_Volmer + beta_BV*U))",
)

eng = KMCEngine(pt, size=[20, 20])
eng.do_steps(10000)
print(eng.get_coverage("H"))
```

See `models/her_Pt111.py` and `validate_her.py` for a complete benchmark.

### Multi-Lattice KMC: Pd(100) ↔ √5-PdO oxide reduction

```python
from spark.types import Project, Site, Condition, Action
from spark.engine import KMCEngine

pt = Project()
pt.set_meta(model_name="Pd100_PdO", model_dimension=2)
pt.add_species(name="empty")
pt.add_species(name="CO")
pt.add_species(name="O")
pt.add_species(name="frozen")  # disabled-site marker (see docs/multi_lattice_design.md)

# Two coexisting sub-lattices on a shared super-cell
pd100 = pt.add_layer(name="Pd100")
pd100.sites.append(Site(name="hollow", default_species="frozen"))
pd100.sites.append(Site(name="bridge", default_species="frozen"))

pdo = pt.add_layer(name="PdO")
pdo.sites.append(Site(name="bridge", default_species="empty"))
pdo.sites.append(Site(name="Olat",   default_species="O"))

# Cross-layer "destruct" process: CO + lattice O → CO2 (g), uncovers Pd hollow
pt.add_process(
    name="destruct",
    conditions=[
        Condition(pt.lattice.generate_coord("bridge.(0,0,0).PdO"), "CO"),
        Condition(pt.lattice.generate_coord("Olat.(0,0,0).PdO"),   "O"),
        Condition(pt.lattice.generate_coord("hollow.(0,0,0).Pd100"), "frozen"),
    ],
    actions=[
        Action(pt.lattice.generate_coord("bridge.(0,0,0).PdO"), "empty"),
        Action(pt.lattice.generate_coord("Olat.(0,0,0).PdO"),   "empty"),
        Action(pt.lattice.generate_coord("hollow.(0,0,0).Pd100"), "empty"),
    ],
    rate_constant="1e9",
)
# ... add ads/des processes on each layer ...

eng = KMCEngine(pt, size=[8, 8])
eng.do_steps(50000)
```

See `examples/multi_lattice_PdO_reduction.py` for the full toy model.

### Off-Lattice OTF KMC (Rust-accelerated)

```python
import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from spark_rs import dimer_find_saddle, fire_minimize

atoms = bulk("Cu", "fcc", a=3.6) * (3, 3, 3)
del atoms[0]                     # vacancy
atoms.calc = EMT()

def force_callback(positions):
    atoms.set_positions(positions)
    return atoms.get_potential_energy(), atoms.get_forces()

# Saddle search via Rust dimer
result = dimer_find_saddle(
    atoms.get_positions(),
    np.random.randn(len(atoms), 3) / np.sqrt(3*len(atoms)),
    force_callback,
    f_tol=0.05, max_steps=100,
)
print("status:", result["status"], "  energy:", result["energy"])
```

See `tests/` for cargo unit tests + `phase_d_bench2.py` for a Cu(100) benchmark comparing Rust vs Python.

### Dynamic Catalytic KMC: PdAu(111) CO oxidation with segregation

```python
from spark.dynamic import DynamicSurface, DynamicKMCEngine, EventGenerator, LookupTableEstimator

surf = DynamicSurface(n_sites=225)
# ... populate sites with Pd/Au atom_type, build NN graph ...

engine = DynamicKMCEngine(
    surface=surf,
    event_generator=EventGenerator(species_list=["empty", "CO", "O"]),
    rate_estimator=LookupTableEstimator(),
    temperature=500.0,
)
engine.run(max_steps=50000)
```

See `examples/dynamic_PdAu_CO_oxidation.py` for the full setup.

---

## Repository Layout

```
SPARK/
├── spark/                         # Python package
│   ├── types.py                   # Project / Layer / Site / Process / LateralInteraction
│   ├── engine.py                  # Lattice + multi-lattice KMC engine
│   ├── microkinetic.py            # Mean-field ODE solver (with v3 SS-robustness patch)
│   ├── polarization.py            # j-U polarization curves from DFT data
│   ├── rates.py                   # kmos-compatible rate expression parser
│   ├── indexing.py                # Multi-lattice spuck-based 1D site indexing
│   ├── offlattice/                # On-the-fly off-lattice KMC (Python reference)
│   │   ├── engine.py              # SKMCEngine main loop (10 steps)
│   │   ├── catalogue.py           # 3-step environment matching + symmetry
│   │   ├── saddle.py              # Dimer search + SaddleMaster (use spark_rs for hot loop)
│   │   ├── minimize.py            # scipy L-BFGS-B minimizer (use spark_rs for hot loop)
│   │   ├── basin.py               # Basin acceleration
│   │   ├── superbasin.py          # SuperBasin acceleration (bac-MRM)
│   │   ├── cache.py               # SuperCache event caching
│   │   └── ...
│   └── dynamic/                   # Dynamic catalytic KMC
│       ├── surface.py             # Mutable graph-based surface
│       ├── events.py              # Catalytic + structural event generators
│       └── ...
│
├── spark-rs/                      # Rust acceleration crate (PyO3 wheel)
│   ├── src/
│   │   ├── lib.rs                 # Crate root + #[pymodule] _native registration
│   │   ├── main.rs                # Standalone CLI binary `spark`
│   │   ├── engine.rs              # Lattice KMC core (Fenwick tree)
│   │   ├── microkinetic.rs        # MKM (Newton-Raphson with damped LS)
│   │   ├── polarization.rs        # j-U interpolation
│   │   ├── python_bindings.rs     # PyO3 dimer_find_saddle + fire_minimize + PyCalculator
│   │   └── offlattice/
│   │       ├── calc.rs            # Calculator trait + MullerBrown + QuadraticSaddle
│   │       ├── minimize.rs        # FIRE pure Rust
│   │       ├── saddle.rs          # DimerSearch (rotate + effective_gradient + collision)
│   │       ├── catalogue.rs       # Environment matching + mechanism storage
│   │       ├── basin.rs / superbasin.rs / cache.rs
│   │       └── engine.rs          # SKMCEngine framework
│   ├── python/spark_rs/           # Python package (re-exports from _native.so)
│   ├── Cargo.toml                 # `python` feature gates pyo3/numpy
│   └── pyproject.toml             # maturin build-backend
│
├── tests/                         # Python unit tests
│   ├── test_multi_lattice_indexing.py    # 11 tests, multi-lattice site indexing
│   └── test_cross_layer_process.py       # 3 tests, cross-layer process firing
│
├── examples/
│   ├── dynamic_PdAu_CO_oxidation.py
│   ├── offlattice_fe_vacancy.py
│   └── multi_lattice_PdO_reduction.py    # 64-cell Pd100+PdO toy, full reduction in 50k steps
│
├── models/
│   ├── her_Pt111.py
│   ├── n2_reduction_Mo.py
│   └── co_adsorption_test.py
│
├── benchmark/                     # Validation against kmos / Zacros 4.0
├── docs/
│   ├── multi_lattice_design.md   # Hoffmann-Reuter 2015 algorithm + SPARK implementation
│   ├── DFT_methodology.md
│   ├── tutorial.md
│   └── ...
│
├── validate_her.py                # 6-test HER/Pt(111) regression suite (~210 s)
├── tutorial.py
└── README.md
```

---

## Validation & Tests

| Test | Coverage | Wall time | Last status |
|---|---|---|---|
| `validate_her.py` | 6 HER tests (Langmuir / Tafel / lateral / diffusion / KMC vs MKM / size convergence) | ~210 s | **6/6 PASS** (2026-05-01) |
| `tests/test_multi_lattice_indexing.py` | 11 tests on multi-lattice spuck-based 1D indexing | < 1 s | **11/11 PASS** |
| `tests/test_cross_layer_process.py` | 3 tests on cross-layer process firing | < 1 s | **3/3 PASS** |
| `examples/multi_lattice_PdO_reduction.py` | 8×8 super-cell, 50k steps, full Pd100/PdO reduction cycle | ~1 s | **64/64 cells reduced, 0 invariant violations** |
| `cargo test --release --lib` (in `spark-rs/`) | 33 Rust unit tests covering Calculator / FIRE / DimerSearch / catalogue / basin etc. | < 5 s | **33/33 PASS** |
| Cu(100) 47-atom slab benchmark (N=20 dimer searches) | Rust vs Python equivalence + speedup | ~7 s | Rust 20/20 success, 3.08× faster than Python; Python 1/20 success before sign fix → 20/20 after |

Run them all:

```bash
# Python regressions
python validate_her.py
python -m pytest tests/                    # if pytest installed
python examples/multi_lattice_PdO_reduction.py

# Rust unit tests (no Python needed)
cd spark-rs && cargo test --release --lib
```

---

## Performance Notes

### Off-lattice OTF saddle search

| | Per dimer search (47-atom Cu(100) slab) | Convergence rate (N=20) |
|---|---|---|
| Python `spark.offlattice.DimerSearch` | 260.6 ms | 1/20 (before 2026-05-01 rotor-sign fix) → ~20/20 (after) |
| Rust `spark_rs.dimer_find_saddle` | 84.5 ms | 20/20 |

**3.08× wall-time speedup** + **20× convergence-rate gain pre-fix** (Python now matches Rust after `8fd2166`). Speedup attributable to:
- Rust dimer rotation arithmetic (no numpy alloc per inner step)
- Trust-radius adapt + collision check inside Rust
- Single Python ↔ Rust crossing per force eval (instead of per-arithmetic)

Force eval (ASE EMT or any ASE calculator) stays Python — that's why the speedup is 3× rather than 30×. For systems where force eval dominates (DFT, MACE), the proportional gain is smaller; for systems where Python overhead dominates (small EMT), it's larger.

### Lattice KMC cross-validation (HER/Pt(111))

| Engine | Language | θ_H | j (current density) | Wall time |
|---|---|---|---|---|
| **SPARK** | Python | 0.338 | reference | 281 s |
| `kmos` | Python + Fortran | 0.332 | < 2% deviation | 20 s |
| `Zacros 4.0` | Fortran | 0.341 | < 2% deviation | 12 s |

All three agree within stochastic noise across 11 potential points (0 to -0.5 V vs RHE). See [`benchmark/BENCHMARK_SUMMARY.md`](benchmark/BENCHMARK_SUMMARY.md).

---

## Algorithm References

| Module | Primary reference |
|---|---|
| Lattice BKL KMC | Bortz, Kalos, Lebowitz, *J. Comput. Phys.* **17**, 10 (1975) |
| **Multi-Lattice kMC** | **Hoffmann, Scheffler, Reuter, *ACS Catal.* **5**, 1199 (2015), DOI 10.1021/cs501352t** |
| Dimer saddle search | Henkelman & Jónsson, *J. Chem. Phys.* **111**, 7010 (1999) |
| FIRE minimizer | Bitzek, Koskinen, Gähler, Moseler, Gumbsch, *PRL* **97**, 170201 (2006) |
| Off-lattice on-the-fly KMC framework | openFLY (C++); ported into SPARK as `spark.offlattice` and `spark-rs/src/offlattice/` |
| bac-MRM SuperBasin | Mason, Hudson, Pellet-Mary, *J. Chem. Theory Comput.* **17**, 5779 (2021) |
| Butler-Volmer PCET barriers | Nørskov et al., *J. Phys. Chem. B* **108**, 17886 (2004) |
| BEP / kinetic scaling | Brønsted, *Chem. Rev.* **5**, 231 (1928); Evans & Polanyi, *Trans. Faraday Soc.* **32**, 1340 (1936) |

---

## Citing SPARK

If SPARK contributes to a publication, please cite the relevant algorithm references above and mention `SPARK <commit-hash>` from this repository.

---

## License

MIT — see [`LICENSE`](LICENSE).

The Rust crate `spark-rs` reads `kmcos` (GPL v3) source for **algorithmic ideas only** in implementing multi-lattice support; no code is copied. SPARK remains MIT-licensed.

---

## Status (2026-05-01)

- ✅ Lattice KMC + multi-lattice (Hoffmann-Reuter 2015) end-to-end
- ✅ Off-lattice OTF KMC algorithm complete (Python reference + Rust hot-loop acceleration)
- ✅ Dynamic catalytic KMC V1 (DynamicSurface + EventCache) with PdAu CO oxidation benchmark
- ✅ Microkinetic v3 SS-robustness patch (handles ill-conditioned Jacobians at near-onset polarization)
- ✅ Rust acceleration via PyO3 wheel — 3.08× speedup, abi3 forward-compatible (one wheel for cpython ≥ 3.10)
- ⏳ Off-lattice catalytic demo (Cu(100) + adsorbates, S2.10 in roadmap) — pending; current off-lattice example is BCC Fe vacancy
- ⏳ `spark.offlattice` Mechanism doesn't yet support adsorption/desorption (`delta_fwd.shape=(n_local, 3)` locks atom count) — workaround on roadmap

See [`docs/multi_lattice_design.md`](docs/multi_lattice_design.md) for design notes on the multi-lattice port.
