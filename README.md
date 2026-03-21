# KMC

**Kinetic Monte Carlo & Microkinetic Modeling for Heterogeneous Catalysis**

动力学蒙特卡洛与微观动力学建模软件 —— 面向多相催化

---

A dual-language (Python + Rust) toolkit for lattice Kinetic Monte Carlo (KMC) simulation and mean-field microkinetic modeling, designed for electrochemical N₂ reduction and general heterogeneous catalysis research. The Python package (`mykmc`) provides a flexible, research-friendly API with SciPy ODE solvers and polarization curve support, while the Rust binary (`mykmc-rs`) delivers high performance with Fenwick tree O(log N) site selection, Newton-Raphson steady-state solver, and zero-allocation coordinate handling.

双语言（Python + Rust）工具包，用于晶格动力学蒙特卡洛（KMC）模拟与平均场微观动力学建模，适用于电化学氮还原及通用多相催化研究。Python 包（`mykmc`）提供灵活的研究接口，支持 SciPy ODE 求解器与极化曲线计算；Rust 二进制（`mykmc-rs`）通过 Fenwick 树 O(log N) 位点选择、Newton-Raphson 稳态求解器和零分配坐标处理实现高性能计算。

## Table of Contents / 目录

- [Features / 功能特性](#features--功能特性)
- [Project Structure / 项目结构](#project-structure--项目结构)
- [Installation / 安装](#installation--安装)
- [Quick Start / 快速开始](#quick-start--快速开始)
- [Two Modes / 两种模式](#two-modes--两种模式)
- [Python Package / Python 包](#python-package--python-包)
- [Rust Command Reference / Rust 命令参考](#rust-command-reference--rust-命令参考)
- [Polarization Curve / 极化曲线](#polarization-curve--极化曲线)
- [Algorithm / 算法](#algorithm--算法)
- [Built-in Model / 内置模型](#built-in-model--内置模型)
- [Architecture / 代码架构](#architecture--代码架构)
- [Tutorial / 教程](#tutorial--教程)
- [Validation / 验证](#validation--验证)
- [Performance / 性能](#performance--性能)
- [Citation / 引用](#citation--引用)

---

## Features / 功能特性

| Feature | Description |
|---------|-------------|
| **BKL KMC Engine** | Rejection-free algorithm with Fenwick tree O(log N) site selection (Rust), O(1) site bookkeeping (Python & Rust) / BKL 无拒绝算法，Fenwick 树 O(log N) 位点选择（Rust），O(1) 位点簿记 |
| **Spatial KMC** | Neighbor list (4-NN 2D / 6-NN 3D with PBC), pairwise lateral interactions, BEP relations, surface diffusion, site types / 邻居列表、横向相互作用、BEP 关系、表面扩散、位点类型 |
| **Mean-Field ODE Solver** | Newton-Raphson with damped line search (Rust), SciPy fsolve (Python) / Newton-Raphson 阻尼线搜索（Rust），SciPy fsolve（Python） |
| **Polarization Curve** | Current density vs. potential from Butler-Volmer or constant-potential DFT data (Python + Rust) / 极化曲线计算，支持 Butler-Volmer 和恒电位 DFT 数据（Python + Rust） |
| **DFT Energy Interpolation** | Potential-dependent barriers from interpolated DFT energy surfaces (Python) / DFT 能量面插值，电位依赖活化能（Python） |
| **Two Modes** | Full model (ads/des/reactions) and cycle mode (pure chemistry) / 完整模型与纯化学循环模式 |
| **Rate Expression Parser** | kmos-compatible string expressions with physical constants (Python & Rust) / 兼容 kmos 的速率表达式解析器 |
| **Electrochemistry** | Butler-Volmer PCET rates with potential dependence / Butler-Volmer 电化学速率 |
| **Parameter Scanning** | Built-in potential and temperature scanning / 内置电位与温度扫描 |
| **JSON Model I/O** | Import/export models in JSON format / JSON 格式模型导入导出 |
| **Validated** | Matches analytical Langmuir isotherm within <1% / Langmuir 等温线验证误差 <1% |
| **High Performance** | Rust binary: ~1.6 MB standalone, Fenwick tree + flat array lookups + zero-alloc coordinates / Rust 二进制约 1.6 MB，Fenwick 树 + 平铺数组查找 + 零堆分配坐标 |

---

## Project Structure / 项目结构

```
KMC/
├── mykmc/                         # Python package / Python 包
│   ├── __init__.py                # Public API exports / 公共 API 导出
│   ├── engine.py                  # BKL KMC engine / KMC 引擎
│   ├── microkinetic.py            # Mean-field ODE solver (scipy) / 平均场 ODE 求解器
│   ├── polarization.py            # Polarization curve & DFT interpolation / 极化曲线
│   ├── rates.py                   # Rate expression parser & rate functions / 速率表达式
│   ├── analysis.py                # Trajectory recording, TOF, steady-state / 分析工具
│   ├── types.py                   # Data model: Project, Species, Process / 数据模型
│   └── units.py                   # Physical constants (kB, h, eV, bar, ...) / 物理常数
├── mykmc-rs/                      # Rust implementation / Rust 高性能实现
│   ├── src/
│   │   ├── main.rs                # CLI entry point with clap / CLI 入口
│   │   ├── engine.rs              # BKL KMC engine with O(1) bookkeeping / KMC 引擎
│   │   ├── microkinetic.rs        # Mean-field ODE solver (RK4 + adaptive Euler) / 求解器
│   │   ├── rates.rs               # Rate expression parser / 速率表达式解析器
│   │   ├── model.rs               # Data model / 数据模型
│   │   ├── models.rs              # Built-in N₂ reduction & CO test / 内置模型
│   │   ├── analysis.rs            # Trajectory recording, steady-state / 分析工具
│   │   └── units.rs               # Physical constants / 物理常数
│   ├── Cargo.toml
│   └── README.md                  # Rust-specific documentation / Rust 专项文档
├── models/                        # Pre-built reaction models / 预置反应模型
│   ├── n2_reduction_Mo.py         # N₂ reduction on Mo surface / Mo 表面 N₂ 还原
│   ├── co_adsorption_test.py      # CO/Pd(100) test model / CO 吸附测试模型
│   └── example_dft_energies.json  # Example DFT energy data / DFT 能量数据示例
├── run_simulation.py              # Main simulation script (CLI) / 主模拟脚本
├── tutorial.py                    # Comprehensive 7-part tutorial / 七部分完整教程
└── co_model.json                  # JSON model: CO on Pd(100) / JSON 模型文件
```

---

## Installation / 安装

### Python / Python 安装

**Prerequisites / 前提条件：** Python ≥ 3.8, NumPy, SciPy

```bash
git clone https://github.com/leshenzhang/KMC.git
cd KMC
pip install numpy scipy
```

No installation step needed — import `mykmc` directly from the project root:

无需安装步骤 —— 直接从项目根目录导入 `mykmc`：

```python
import sys; sys.path.insert(0, '/path/to/KMC')
from mykmc import Project, KMCEngine, MicroKineticModel
```

### Rust / Rust 安装

**Prerequisites / 前提条件：** Rust toolchain (≥ 1.70): https://rustup.rs

如果尚未安装 Rust：

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

### Build / 编译

```bash
cd KMC/mykmc-rs
cargo build --release
```

The binary is at `mykmc-rs/target/release/mykmc`.

编译后的可执行文件位于 `mykmc-rs/target/release/mykmc`。

### HPC Note / 超算注意事项

On Cray HPC systems (e.g., Shaheen), `cc` is the Cray compiler wrapper which conflicts with the Rust linker. A `.cargo/config.toml` is included to use `/usr/bin/gcc` instead:

在 Cray 超算系统上，`cc` 是 Cray 编译器包装器，与 Rust 链接器冲突。项目已包含 `.cargo/config.toml` 配置使用 `/usr/bin/gcc`：

```toml
[target.x86_64-unknown-linux-gnu]
linker = "/usr/bin/gcc"
```

If Rust is not installed system-wide, install to a custom path:

如果无法全局安装 Rust，可安装到自定义路径：

```bash
export RUSTUP_HOME=/scratch/your_user/rust_env/rustup
export CARGO_HOME=/scratch/your_user/rust_env/cargo
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --no-modify-path
export PATH="$CARGO_HOME/bin:$PATH"
```

---

## Quick Start / 快速开始

### Python

```bash
# Run both mean-field and KMC models
# 同时运行平均场和 KMC 模型
python run_simulation.py

# Mean-field microkinetic model only
# 仅运行平均场微观动力学模型
python run_simulation.py --mkm-only --T 300 --U -0.5

# KMC simulation only
# 仅运行 KMC 模拟
python run_simulation.py --kmc-only --T 300 --U -0.5 --lattice-size 20

# Scan applied potential
# 扫描施加电位
python run_simulation.py --scan-potential

# Scan temperature
# 扫描温度
python run_simulation.py --scan-temp

# Compute polarization curve (Butler-Volmer)
# 计算极化曲线（Butler-Volmer 模型）
python run_simulation.py --polarization

# Polarization curve with DFT data
# 使用 DFT 数据计算极化曲线
python run_simulation.py --polarization --dft-data models/example_dft_energies.json
```

### Rust

```bash
cd mykmc-rs

# Validate the installation (CO/Pd100 Langmuir isotherm test)
# 验证安装（CO/Pd100 Langmuir 等温线测试）
./target/release/mykmc validate

# Run mean-field microkinetic model for N₂ reduction
# 运行 N₂ 还原平均场微观动力学模型
./target/release/mykmc mkm --cycle --t 300 --u -1.0

# Run lattice KMC simulation
# 运行晶格 KMC 模拟
./target/release/mykmc kmc --cycle --t 300 --u -5.0 --size 20 --steps 500000

# Scan applied potential
# 扫描施加电位
./target/release/mykmc scan-u --cycle --t 300 --u-min -0.5 --u-max -3.0

# Scan temperature
# 扫描温度
./target/release/mykmc scan-t --cycle --u -1.0 --t-min 250 --t-max 500
```

---

## Two Modes / 两种模式

The software supports two complementary simulation modes (both Python and Rust):

软件支持两种互补的模拟模式（Python 和 Rust 均支持）：

### Full Mode (default) / 完整模式（默认）

Includes adsorption, desorption, and all reaction steps. Supports empty sites on the lattice.

包含吸附、脱附和所有反应步骤。格点上允许空位。

```bash
# Python
python run_simulation.py --T 300 --U -1.0

# Rust
./target/release/mykmc mkm --t 300 --u -1.0
./target/release/mykmc kmc --t 300 --u -1.0 --size 20
```

- 11 species: empty, N₂, NNH, HNNH, NNH₂, HNNH₂, H₂NNH₂, N, NH, NH₂, NH₃
- N₂ adsorption/desorption (Hertz-Knudsen kinetics)
- NH₃ desorption
- Suitable for studying coverage effects and competitive adsorption

### Cycle Mode (`--cycle`) / 催化循环模式

Pure chemical reaction steps only. No adsorption, desorption, migration, or surface defects. All sites are always occupied (Σθᵢ = 1). NH₃ is released instantly and the site regenerates to N₂.

仅包含纯化学反应步骤。不考虑吸附、脱附、迁移和表面缺陷。所有位点始终被占据（Σθᵢ = 1）。NH₃ 即时释放，位点回到 N₂。

```bash
# Rust
./target/release/mykmc mkm --cycle --t 300 --u -1.0
./target/release/mykmc kmc --cycle --t 300 --u -5.0 --size 20
```

- 9 species: N₂, NNH, HNNH, NNH₂, HNNH₂, H₂NNH₂, N, NH, NH₂
- No empty sites
- Catalytic cycle: NH₂ + H⁺/e⁻ → N₂ + NH₃↑ (cycle closure)
- Suitable for fast microkinetic screening and mechanism analysis

### Comparison / 对比

| | Full Mode / 完整模式 | Cycle Mode / 循环模式 |
|---|---|---|
| Species / 物种数 | 11 (with empty, NH₃) | 9 (no empty, no NH₃) |
| Coverage / 覆盖度 | Σθᵢ ≤ 1, empty sites | Σθᵢ = 1, always |
| Adsorption / 吸附 | Hertz-Knudsen | None / 无 |
| Desorption / 脱附 | N₂, NH₃ desorption | None / 无 |
| Use case / 用途 | Full KMC with spatial effects | Fast microkinetic screening |

---

## Python Package / Python 包

### Core API / 核心 API

#### 1. Define a model / 定义模型

```python
from mykmc import Project, Site, Layer, Condition, Action
import numpy as np

pt = Project()
pt.set_meta(author='User', model_name='CO_on_Pd100', model_dimension=2)

# Species (first = default) / 物种（第一个为默认）
pt.add_species(name='empty', color='#ffffff')
pt.add_species(name='CO', color='#ff0000')

# Lattice / 晶格
layer = pt.add_layer(name='surface')
layer.sites.append(Site(name='hollow', pos=(0.5, 0.5, 0.5), default_species='empty'))
pt.lattice.cell = np.diag([3.5, 3.5, 10.0])

# Parameters / 参数
pt.add_parameter(name='T', value=600.0, adjustable=True, min=400, max=800)
pt.add_parameter(name='p_CO', value=1.0, adjustable=True, min=1e-10, max=100)
pt.add_parameter(name='A', value='(3.5*angstrom)**2')
pt.add_parameter(name='deltaG', value=-0.5, adjustable=True, min=-1.3, max=0.3)

# Processes / 基元过程
coord = pt.lattice.generate_coord('hollow')
pt.add_process(
    name='CO_adsorption',
    rate_constant='p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)',
    conditions=[Condition(coord=coord, species='empty')],
    actions=[Action(coord=coord, species='CO')],
    tof_count={'CO_adsorption': 1},
)
```

#### 2. Run KMC simulation / 运行 KMC 模拟

```python
from mykmc import KMCEngine, TrajectoryRecorder

engine = KMCEngine(pt, size=[20, 20], print_rates=True, banner=True)
engine.parameters.T = 600.0

recorder = TrajectoryRecorder(engine)
for _ in range(200):
    engine.do_steps(5000)
    recorder.record()

engine.print_coverages()
print(engine.get_tof())
```

#### 3. Mean-field microkinetic model / 平均场微观动力学模型

```python
from mykmc import MicroKineticModel
from mykmc.rates import tst_rate, electrochemical_rate

mkm = MicroKineticModel()
mkm.add_species('N2')
mkm.add_species('NNH')
mkm.add_species('NH3')

mkm.add_reaction(
    name='N2_hydrogenation',
    reactants={'N2': 1},
    products={'NNH': 1},
    rate_fwd=lambda p: tst_rate(2.82, p['T']),
    rate_rev=lambda p: tst_rate(1.50, p['T']),
)

mkm.parameters = {'T': 300, 'U': -0.5}

# Steady state / 稳态求解
ss = mkm.solve_steady_state()
mkm.print_summary(ss)

# Time-dependent ODE / 瞬态 ODE 求解
sol = mkm.solve_ode(t_span=(0, 1e3), t_eval=np.logspace(-6, 3, 500))

# Parameter scanning / 参数扫描
results = mkm.scan_parameter('U', np.linspace(-2, 0, 41), observable='NH3_production')
```

#### 4. Polarization curve / 极化曲线

```python
from mykmc.polarization import PolarizationCurve, EnergyLandscape, load_energy_data

# Mode 1: Butler-Volmer model / 模式1：Butler-Volmer 模型
pc = PolarizationCurve(mkm, n_electrons={'NH3_production': 3}, A_site=(3.2e-10)**2)
results = pc.compute(U_range=np.linspace(-2, 0, 41), T=300, verbose=True)
pc.print_results()
pc.save('polarization_curve.dat')

# Mode 2: DFT constant-potential data / 模式2：恒电位 DFT 数据
landscape = load_energy_data('models/example_dft_energies.json')
```

### Rate Expression Parser / 速率表达式解析器

Supports kmos-compatible string expressions with:

支持兼容 kmos 的字符串速率表达式：

- Physical constants / 物理常数: `kB`, `h`, `eV`, `bar`, `angstrom`, `umass`, `pi`
- Derived variables / 导出变量: `beta` = 1/(kB×T)
- Molecular masses / 分子质量: `m_N2`, `m_CO`, `m_NH3`, etc.
- Math functions / 数学函数: `exp()`, `log()`, `sqrt()`, `sin()`, `cos()`, `pow()`, `abs()`
- All user-defined parameters / 所有用户定义参数
- Operators / 运算符: `+`, `-`, `*`, `/`, `**` (power)

**Example / 示例：**

```
kB*T/h*exp(-(Ea_N2_to_NNH + beta_BV*U)*eV/(kB*T))
```

---

## Rust Command Reference / Rust 命令参考

### `mykmc kmc` — Lattice KMC Simulation / 晶格 KMC 模拟

```
mykmc kmc [OPTIONS]

Options:
  --t <T>          Temperature [K] (default: 300) / 温度
  --p-n2 <P_N2>   N₂ partial pressure [bar] (default: 1, full mode only) / N₂ 分压
  --u <U>          Applied potential [V vs RHE] (default: -1) / 施加电位
  --size <SIZE>    Lattice side length (default: 20, creates SIZE×SIZE grid) / 格点边长
  --steps <STEPS>  Number of KMC steps (default: 1000000) / KMC 步数
  --cycle          Use pure catalytic cycle mode / 使用纯化学循环模式
```

**Example / 示例：**

```bash
# Full model, 30×30 lattice, 2 million steps
./target/release/mykmc kmc --t 300 --u -1.0 --size 30 --steps 2000000

# Cycle mode, compact lattice
./target/release/mykmc kmc --cycle --t 300 --u -5.0 --size 15 --steps 500000
```

### `mykmc mkm` — Mean-Field Microkinetic Model / 平均场微观动力学

```
mykmc mkm [OPTIONS]

Options:
  --t <T>        Temperature [K] (default: 300)
  --p-n2 <P_N2>  N₂ partial pressure [bar] (default: 1, full mode only)
  --u <U>        Applied potential [V vs RHE] (default: -1)
  --cycle        Use pure catalytic cycle mode
```

**Example / 示例：**

```bash
./target/release/mykmc mkm --cycle --t 300 --u -1.5
```

**Output includes / 输出包含：**
- Steady-state surface coverages / 稳态表面覆盖度
- Turn-over frequencies (TOF) / 转化频率
- Parameter summary / 参数摘要

### `mykmc scan-u` — Potential Scan / 电位扫描

```
mykmc scan-u [OPTIONS]

Options:
  --t <T>              Temperature [K] (default: 300)
  --u-min <U_MIN>      Start potential [V] (default: -0.5)
  --u-max <U_MAX>      End potential [V] (default: -2.0)
  --u-steps <U_STEPS>  Number of scan points (default: 31)
  --cycle              Use cycle mode
```

**Example / 示例：**

```bash
./target/release/mykmc scan-u --cycle --t 300 --u-min -0.5 --u-max -3.0 --u-steps 51
```

### `mykmc scan-t` — Temperature Scan / 温度扫描

```
mykmc scan-t [OPTIONS]

Options:
  --u <U>              Applied potential [V] (default: -1)
  --t-min <T_MIN>      Start temperature [K] (default: 250)
  --t-max <T_MAX>      End temperature [K] (default: 500)
  --t-steps <T_STEPS>  Number of scan points (default: 11)
  --cycle              Use cycle mode
```

### `mykmc validate` — Validation Test / 验证测试

Runs CO adsorption/desorption on Pd(100) and compares KMC results with the analytical Langmuir isotherm at 5 different free energies.

运行 CO/Pd(100) 吸附脱附模拟，在 5 个不同自由能下与解析 Langmuir 等温线对比。

```bash
./target/release/mykmc validate
```

### `mykmc polarization` — Polarization Curve from DFT / DFT 极化曲线

Compute electrochemical polarization curve (current density vs. potential) directly from constant-potential DFT data. Uses cubic spline interpolation for potential-dependent barriers and steady-state continuation for robust convergence.

从恒电位 DFT 数据直接计算电化学极化曲线（电流密度 vs. 电位）。使用三次样条插值处理电位依赖的能垒，稳态延续法确保收敛稳定性。

```
mykmc polarization [OPTIONS]

Options:
  --dft-data <PATH>    Path to DFT energy JSON file (required) / DFT 能量 JSON 文件路径
  --t <T>              Temperature [K] (default: 300) / 温度
  --p-n2 <P_N2>        N₂ partial pressure [bar] (default: 1) / N₂ 分压
  --u-min <U_MIN>      Start potential [V vs RHE] (default: -2.0) / 起始电位
  --u-max <U_MAX>      End potential [V vs RHE] (default: 0.0) / 终止电位
  --u-steps <U_STEPS>  Number of scan points (default: 50) / 扫描点数
  --output <FILE>      Output filename (default: polarization_curve.dat) / 输出文件名
```

**Example / 示例：**

```bash
# Compute polarization curve from DFT data
./target/release/mykmc polarization --dft-data models/example_dft_energies.json --t 300 --u-min -2.0 --u-max 0.0 --u-steps 100

# Custom output file
./target/release/mykmc polarization --dft-data my_dft.json --output my_curve.dat
```

**DFT JSON format / DFT JSON 格式：**

```json
{
  "potentials": [-2.0, -1.5, -1.0, -0.5, 0.0],
  "state_energies": {
    "clean": [0.0, 0.0, 0.0, 0.0, 0.0],
    "N2*": [-0.80, -0.72, -0.60, -0.50, -0.38],
    "NNH*": [1.50, 1.65, 1.82, 2.00, 2.20]
  },
  "ts_energies": {
    "N2*->NNH*": [1.80, 2.00, 2.30, 2.60, 3.00]
  }
}
```

---

## Polarization Curve / 极化曲线

Both Python and Rust support computing polarization curves. The Rust implementation uses a custom cubic spline interpolator (no external dependency), while Python uses scipy.interpolate.

Python 和 Rust 均支持极化曲线计算。Rust 实现使用自研三次样条插值器（无外部依赖），Python 使用 scipy.interpolate。

### Mode 1: Butler-Volmer model / 模式1：Butler-Volmer 模型

Uses the `MicroKineticModel` with electrochemical rate functions. The TOF is converted to current density via:

使用 `MicroKineticModel` 和电化学速率函数，TOF 通过以下公式转换为电流密度：

**j = TOF × n_e × e × N_site / A_geo**

where n_e is the number of electrons per turnover.

其中 n_e 为每次转化的电子数。

### Mode 2: DFT constant-potential data / 模式2：恒电位 DFT 数据

Provide DFT-computed energies of adsorbates and transition states at multiple potentials. The `EnergyLandscape` class interpolates barriers as functions of potential:

提供多个电位下 DFT 计算的吸附物和过渡态能量。`EnergyLandscape` 类将活化能插值为电位的函数：

```json
{
  "potentials": [-2.0, -1.5, -1.0, -0.5, 0.0],
  "state_energies": {
    "clean":   [0.00, 0.00, 0.00, 0.00, 0.00],
    "N2*":     [-0.80, -0.72, -0.60, -0.50, -0.38],
    "NNH*":    [1.50, 1.65, 1.82, 2.00, 2.20]
  },
  "ts_energies": {
    "N2*->NNH*": [1.80, 2.00, 2.30, 2.60, 3.00]
  }
}
```

Forward and reverse barriers are automatically computed: Ea_fwd = E_TS − E_IS, Ea_rev = E_TS − E_FS.

正向和逆向活化能自动计算：Ea_fwd = E_TS − E_IS, Ea_rev = E_TS − E_FS。

---

## Algorithm / 算法

### KMC: BKL Rejection-Free Algorithm / BKL 无拒绝算法

The KMC engine implements the Variable Step-Size Method (VSSM), also known as the BKL algorithm or n-fold way. This is a rejection-free algorithm where every step advances the simulation clock and changes the configuration.

KMC 引擎实现了可变步长法（VSSM），也称 BKL 算法或 n-fold way。这是一种无拒绝算法——每一步都推进模拟时钟并改变构型。

**Each KMC step / 每步 KMC：**

1. Compute cumulative rates: R_i = Σ_{j≤i} k_j × n_j
2. Draw 3 random numbers: r₁, r₂, r₃ ~ U(0,1)
3. Advance time: Δt = −ln(r₁) / R_total (Poisson process)
4. Select process by binary search on cumulative rates using r₂
5. Select site uniformly from available sites using r₃
6. Execute process and update bookkeeping

**Key optimizations (Rust) / 关键优化（Rust）：**

| Optimization | Complexity | Description |
|---|---|---|
| **Fenwick tree site selection** | O(log N) | Binary Indexed Tree for rate-weighted site selection under lateral interactions (was O(N) linear scan) / Fenwick 树实现横向相互作用下的加权位点选择 |
| **Flat lateral energy array** | O(1) | `lateral_energy[sp1 * N + sp2]` replaces HashMap lookup / 平铺数组替代 HashMap |
| **Direct-indexed availability** | O(1) | `site_in_avail[proc][site]` replaces HashMap / 直接索引数组替代 HashMap |
| **Fixed [i32;3] coordinates** | 0 alloc | Stack-allocated coordinates, no heap allocation in hot path / 栈上固定大小坐标，热路径零堆分配 |
| **Neighbor-list fast path** | 5 vs 9 sites | Single-site processes use neighbor list directly (5 sites in 2D vs 9 from grid) / 单位点过程直接用邻居列表 |
| **Cached beta_thermal** | O(1) | `eV/(kB*T)` cached, updated only on parameter change / 缓存热力学因子 |
| **Periodic Fenwick rebuild** | correctness | Rebuild every 100K steps to prevent floating-point drift / 每10万步重建防浮点漂移 |
| **OnceLock rate constants** | O(1) | Physical constants cached via `std::sync::OnceLock`, only user params merged per call / 物理常数一次性缓存 |

### Mean-Field Microkinetic Solver / 平均场微观动力学求解器

Solves the coupled ODE system:

求解耦合 ODE 方程组：

dθᵢ/dt = Σⱼ νᵢⱼ · rⱼ

where νᵢⱼ is the stoichiometric coefficient and rⱼ is the net rate of reaction j.

其中 νᵢⱼ 为化学计量系数，rⱼ 为反应 j 的净速率。

**Python solver / Python 求解器：**
- `scipy.integrate.solve_ivp` — adaptive time-stepping ODE integration / 自适应步长 ODE 积分
- `scipy.optimize.fsolve` — steady-state root finding / 稳态根求解

**Rust solver / Rust 求解器：**
- **Phase 1**: Euler pre-equilibration with progressive dt (10⁻¹⁵ → 10⁰ s) to reach vicinity of steady state / 阶段1：自适应 Euler 预平衡
- **Phase 2**: Newton-Raphson with finite-difference Jacobian, damped line search with backtracking, Gaussian elimination with partial pivoting / 阶段2：Newton-Raphson + 有限差分 Jacobian + 阻尼线搜索 + 部分选主元高斯消元
- Typically converges in 5-10 Newton iterations (vs 8000 Euler steps in v0.2) / 通常 5-10 步 Newton 迭代收敛（v0.2 需 8000 步 Euler）
- Fallback to Euler if Newton stalls (singular Jacobian or no descent) / Newton 停滞时自动回退 Euler
- Convergence check: max|dθ/dt| < 10⁻¹⁵ / 收敛判据

---

## Built-in Model / 内置模型

### Electrochemical N₂ Reduction on Mo Surface / Mo 表面电化学 N₂ 还原

The built-in model implements the nitrogen reduction reaction (NRR) with two competing pathways:

内置模型实现了氮还原反应（NRR），包含两条竞争路径：

**Distal pathway / 远端路径：**

```
N₂* → NNH* → NNH₂* → N* + NH₃↑ → NH* → NH₂* → NH₃↑
```

**Alternating pathway / 交替路径：**

```
N₂* → NNH* → HNNH* → HNNH₂* → H₂NNH₂* → NH₂* → NH₃↑
```

**DFT activation barriers / DFT 活化能 [eV]：**

| Step / 步骤 | Barrier / 活化能 | Type / 类型 |
|---|---|---|
| N₂* → NNH* | 2.82 | PCET |
| NNH* → HNNH* | 1.14 | PCET (alternating) |
| NNH* → NNH₂* | 3.22 | PCET (distal) |
| HNNH* → HNNH₂* | 3.32 | PCET |
| NNH₂* → N* + NH₃ | 2.94 | Thermal (N-N cleavage) |
| HNNH₂* → H₂NNH₂* | 3.36 | PCET |
| H₂NNH₂* → NH₂* | 5.14 | Thermal (N-N cleavage) |
| N* → NH* | 2.01 | PCET |
| NH* → NH₂* | 2.68 | PCET |
| NH₂* → NH₃* | 4.35 | PCET |

PCET steps follow Butler-Volmer kinetics: k = (kB·T/h) · exp(−(Ea + β_BV·U)·eV/(kB·T))

Thermal steps (N-N bond cleavage) are potential-independent: k = (kB·T/h) · exp(−Ea·eV/(kB·T))

PCET 步骤遵循 Butler-Volmer 动力学，热反应步骤（N-N 断键）不受电位影响。

---

## Architecture / 代码架构

### Python (`mykmc/`)

```
mykmc/
├── __init__.py          Public API: Project, KMCEngine, MicroKineticModel, ...
├── types.py             Data model: Project, Species, Site, Layer, Process, Condition, Action
├── engine.py            BKL KMC engine with rejection-free algorithm
├── microkinetic.py      Mean-field ODE solver (scipy solve_ivp + fsolve)
├── polarization.py      InterpolatedBarrier, EnergyLandscape, PolarizationCurve
├── rates.py             Rate expression parser + TST/Hertz-Knudsen/Butler-Volmer functions
├── analysis.py          TrajectoryRecorder, run_to_steady_state
└── units.py             Physical constants: kB, h, eV, bar, angstrom, umass, ...
```

### Rust (`mykmc-rs/src/`)

```
src/
├── main.rs              CLI entry point with clap / CLI 入口
├── lib.rs               Module declarations / 模块声明
├── units.rs             Physical constants (kB, h, eV, bar, ...) / 物理常数
├── model.rs             Data model: Project, Species, Process, LateralInteraction, BEPRelation / 数据模型
├── rates.rs             Rate expression parser (OnceLock cached) + TST/HK/BV functions / 速率计算
├── engine.rs            BKL KMC engine: Fenwick tree, flat arrays, neighbor list, lateral/BEP / KMC 引擎
├── microkinetic.rs      Mean-field solver: Euler pre-equilibration + Newton-Raphson + Gauss elimination / 求解器
├── analysis.rs          Trajectory recording, steady-state detection / 分析工具
├── models.rs            Built-in N₂ reduction (full + cycle) and CO test / 内置模型
└── polarization.rs      Cubic spline interpolation, DFT energy landscape, polarization curve / 极化曲线
```

### Key Design / 设计要点

- **Dual-language** — Python for flexibility & analysis, Rust for performance / 双语言：Python 灵活分析，Rust 高性能计算
- **kmos-compatible API** — same concepts: Project, Species, Process, Condition, Action / 兼容 kmos 的 API 概念
- **String-based rate expressions** — parsed at runtime, no recompilation needed / 字符串速率表达式，运行时解析
- **Generic engine** — model-agnostic KMC core, easy to add new models / 通用引擎，易于添加新模型
- **JSON model I/O** — share models between Python and Rust / JSON 模型在 Python 和 Rust 间共享

---

## Tutorial / 教程

A comprehensive 7-part tutorial (`tutorial.py`) is included, covering the full workflow from model building to analysis:

项目包含完整的七部分教程（`tutorial.py`），覆盖从模型构建到分析的全流程：

| Part / 部分 | Topic / 主题 |
|---|---|
| 1 | Build a simple model — CO adsorption/desorption on Pd(100) / 构建简单模型 |
| 2 | Run KMC simulation / 运行 KMC 模拟 |
| 3 | Dynamic parameter adjustment & coverage analysis / 参数动态调节与覆盖度分析 |
| 4 | Build complex reaction networks — N₂ reduction / 构建复杂反应网络 |
| 5 | Mean-field microkinetic ODE solving / 平均场微观动力学 ODE 求解 |
| 6 | Parameter scanning (potential, temperature) & rate control analysis / 参数扫描与速率控制度分析 |
| 7 | Results visualization / 结果可视化 |

```bash
python tutorial.py           # Run all parts / 运行全部教程
python tutorial.py --part 1  # Run part 1 only / 只运行第1部分
python tutorial.py --part 5  # Run part 5 only / 只运行第5部分
```

---

## Validation / 验证

The KMC engine is validated against the analytical Langmuir isotherm for CO adsorption/desorption on Pd(100):

KMC 引擎通过 CO/Pd(100) 吸附脱附与解析 Langmuir 等温线对比验证：

```
$ ./target/release/mykmc validate

T=600 K, p_CO=1 bar, 30x30 lattice

    dG[eV]        KMC        Langmuir   Status
------------------------------------------------
      -0.5     0.9998          0.9999     PASS
      -0.3     0.9964          0.9970     PASS
      -0.1     0.8729          0.8737     PASS
       0.0     0.4986          0.5000     PASS
       0.1     0.1262          0.1263     PASS

All PASSED!
```

Note: KMC is stochastic — exact values vary slightly between runs, but all stay within <1% of the analytical Langmuir isotherm.

注意：KMC 是随机模拟，每次运行数值会有微小波动，但均在解析 Langmuir 等温线的 1% 误差范围内。

All tests pass with relative error < 1%.

所有测试通过，相对误差 < 1%。

---

## Performance / 性能

| Benchmark | Result |
|---|---|
| CO validation (5×700k steps) / CO 验证 | **3.7 s** |
| Potential scan (31 points, Newton-Raphson) / 电位扫描 | **0.6 s** |
| Binary size / 可执行文件大小 | **1.6 MB** |
| Compilation / 编译时间 | ~7 s (release) |
| Unit tests / 单元测试 | **19/19 pass** |

### Algorithmic Complexity / 算法复杂度

| Operation | Without Lateral | With Lateral (Fenwick) |
|---|---|---|
| Process selection / 过程选择 | O(log N_proc) | O(log N_proc) |
| Site selection / 位点选择 | O(1) | **O(log N_sites)** |
| Site add/remove / 位点增删 | O(1) | **O(log N_sites)** |
| MKM steady state / MKM 稳态 | — | **5-10 Newton iterations** |

与 Python 版本相比：KMC 模拟通常快 20-100 倍，含横向相互作用的大格点可达 100-500 倍。

---

## Citation / 引用

If you use this software in your research, please cite:

如果在研究中使用本软件，请引用：

```
KMC: Kinetic Monte Carlo & Microkinetic Modeling for Heterogeneous Catalysis
https://github.com/leshenzhang/KMC
```

---

## License / 许可证

MIT License

---

## Acknowledgments / 致谢

- Algorithm inspired by [kmos](https://github.com/mhoffman/kmos) (M. Hoffmann et al.)
- DFT activation barriers from VASP calculations on Mo surface
- Developed at Wanlu Li Group, UC San Diego

算法参考 kmos 软件（M. Hoffmann 等），DFT 活化能来自 VASP 计算（Mo 表面），开发于加州大学圣地亚哥分校 Wanlu Li 课题组。
