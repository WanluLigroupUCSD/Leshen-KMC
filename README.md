# Leshen-KMC

**Kinetic Monte Carlo & Microkinetic Modeling for Heterogeneous Catalysis**

动力学蒙特卡洛与微观动力学建模软件 —— 面向多相催化

---

A high-performance Rust implementation of lattice Kinetic Monte Carlo (KMC) and mean-field microkinetic modeling, designed for electrochemical N₂ reduction and general heterogeneous catalysis research.

高性能 Rust 实现的晶格动力学蒙特卡洛（KMC）与平均场微观动力学建模软件，适用于电化学氮还原及通用多相催化研究。

## Table of Contents / 目录

- [Features / 功能特性](#features--功能特性)
- [Installation / 安装](#installation--安装)
- [Quick Start / 快速开始](#quick-start--快速开始)
- [Two Modes / 两种模式](#two-modes--两种模式)
- [Command Reference / 命令参考](#command-reference--命令参考)
- [Algorithm / 算法](#algorithm--算法)
- [Built-in Model / 内置模型](#built-in-model--内置模型)
- [Architecture / 代码架构](#architecture--代码架构)
- [Validation / 验证](#validation--验证)
- [Citation / 引用](#citation--引用)

---

## Features / 功能特性

| Feature | Description |
|---------|-------------|
| **BKL KMC Engine** | Rejection-free variable step-size algorithm with O(1) site bookkeeping / BKL 无拒绝算法，O(1) 位点簿记 |
| **Mean-Field ODE Solver** | Adaptive time-stepping steady-state solver (Euler + RK4) / 自适应步长稳态求解器 |
| **Two Modes** | Full model (ads/des/reactions) and cycle mode (pure chemistry) / 完整模型与纯化学循环模式 |
| **Rate Expression Parser** | kmos-compatible string expressions with physical constants / 兼容 kmos 的速率表达式解析器 |
| **Electrochemistry** | Butler-Volmer PCET rates with potential dependence / Butler-Volmer 电化学速率 |
| **Parameter Scanning** | Built-in potential and temperature scanning / 内置电位与温度扫描 |
| **Validated** | Matches analytical Langmuir isotherm within <1% / Langmuir 等温线验证误差 <1% |
| **Single Binary** | ~1.3 MB standalone executable, no dependencies / 约 1.3 MB 独立可执行文件，无外部依赖 |

---

## Installation / 安装

### Prerequisites / 前提条件

- **Rust toolchain** (≥ 1.70): https://rustup.rs
- **GCC** (for linking on HPC systems where `cc` is the Cray compiler)

如果尚未安装 Rust：

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

### Build / 编译

```bash
git clone https://github.com/WanluLigroupUCSD/Leshen-KMC.git
cd Leshen-KMC
cargo build --release
```

The binary is at `target/release/mykmc`.

编译后的可执行文件位于 `target/release/mykmc`。

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

```bash
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

The software supports two complementary simulation modes:

软件支持两种互补的模拟模式：

### Full Mode (default) / 完整模式（默认）

Includes adsorption, desorption, and all reaction steps. Supports empty sites on the lattice.

包含吸附、脱附和所有反应步骤。格点上允许空位。

```bash
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

## Command Reference / 命令参考

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

**Key optimizations / 关键优化：**
- O(1) available-site add/remove via swap-with-last trick
- Binary search for process selection: O(log N_proc)
- Periodic boundary conditions

### Mean-Field Microkinetic Solver / 平均场微观动力学求解器

Solves the coupled ODE system:

dθᵢ/dt = Σⱼ νᵢⱼ · rⱼ

where νᵢⱼ is the stoichiometric coefficient and rⱼ is the net rate of reaction j.

求解耦合 ODE 方程组，其中 νᵢⱼ 为化学计量系数，rⱼ 为反应 j 的净速率。

**Steady-state solver / 稳态求解器：**
- Adaptive time-stepping Euler integration
- Progressive dt from 10⁻¹⁵ to 10⁶ s
- Convergence check: max|dθ/dt| < 10⁻¹⁵

### Rate Expression Parser / 速率表达式解析器

Supports kmos-compatible string expressions with:

- Physical constants: `kB`, `h`, `eV`, `bar`, `angstrom`, `umass`, `pi`
- Derived variables: `beta` = 1/(kB×T)
- Molecular masses: `m_N2`, `m_CO`, `m_NH3`, etc.
- Math functions: `exp()`, `log()`, `sqrt()`, `sin()`, `cos()`, `pow()`, `abs()`
- All user-defined parameters
- Operators: `+`, `-`, `*`, `/`, `**` (power)

**Example / 示例：**

```
kB*T/h*exp(-(Ea_N2_to_NNH + beta_BV*U)*eV/(kB*T))
```

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

```
src/
├── main.rs           CLI entry point with clap / CLI 入口
├── lib.rs            Module declarations / 模块声明
├── units.rs          Physical constants (kB, h, eV, bar, ...) / 物理常数
├── model.rs          Data model: Project, Species, Process, Condition, Action / 数据模型
├── rates.rs          Rate expression parser + Arrhenius/TST/HK/BV functions / 速率计算
├── engine.rs         BKL KMC engine with O(1) bookkeeping / KMC 引擎
├── microkinetic.rs   Mean-field ODE solver (RK4 + adaptive Euler) / 平均场求解器
├── analysis.rs       Trajectory recording, steady-state detection / 分析工具
└── models.rs         Built-in N₂ reduction (full + cycle) and CO test / 内置模型
```

### Key Design / 设计要点

- **Pure Rust, zero external runtime dependencies** — single static binary / 纯 Rust，无运行时依赖，单一静态二进制
- **kmos-compatible API** — same concepts: Project, Species, Process, Condition, Action / 兼容 kmos 的 API 概念
- **String-based rate expressions** — parsed at runtime, no recompilation needed / 字符串速率表达式，运行时解析
- **Generic engine** — model-agnostic KMC core, easy to add new models / 通用引擎，易于添加新模型

---

## Validation / 验证

The KMC engine is validated against the analytical Langmuir isotherm for CO adsorption/desorption on Pd(100):

KMC 引擎通过 CO/Pd(100) 吸附脱附与解析 Langmuir 等温线对比验证：

```
T = 600 K, p_CO = 1 bar, 30×30 lattice

  ΔG [eV]    θ_KMC    θ_Langmuir    Status
    -0.5     0.9999       0.9999      PASS
    -0.3     0.9968       0.9970      PASS
    -0.1     0.8742       0.8737      PASS
     0.0     0.5002       0.5000      PASS
    +0.1     0.1261       0.1263      PASS
```

All tests pass with relative error < 1%.

所有测试通过，相对误差 < 1%。

---

## Performance / 性能

| Benchmark | Result |
|---|---|
| CO validation (5 models × 700k steps) / CO 验证 | **3.7 s** |
| Binary size / 可执行文件大小 | **1.3 MB** |
| Compilation / 编译时间 | ~8 s (release) |

Compared to Python: ~24× speedup for KMC simulations.

与 Python 版本相比：KMC 模拟约 24 倍加速。

---

## Citation / 引用

If you use this software in your research, please cite:

如果在研究中使用本软件，请引用：

```
Leshen-KMC: Kinetic Monte Carlo & Microkinetic Modeling for Heterogeneous Catalysis
https://github.com/WanluLigroupUCSD/Leshen-KMC
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
