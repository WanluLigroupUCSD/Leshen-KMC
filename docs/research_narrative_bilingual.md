# Research Narrative: Spatial KMC Reveals Selectivity Control in Electrochemical NO₃RR on Cu
# 研究思路：空间KMC揭示Cu上电化学NO₃RR选择性调控机制

---

## 1. Research Background and Motivation / 研究背景与动机

### 1.1 The selectivity puzzle in NO₃RR / NO₃RR中的选择性难题

Electrochemical nitrate reduction (NO₃RR) on Cu is a promising route for sustainable ammonia synthesis and wastewater remediation. However, NO₃RR produces multiple products — NH₃ (desired, 8e⁻), N₂ (10e⁻), N₂O (8e⁻), and the parasitic H₂ (HER, 2e⁻) — and the selectivity between them remains poorly understood at the mechanistic level.

电化学硝酸盐还原(NO₃RR)在Cu上是实现可持续合成氨和废水治理的有前景路径。然而，NO₃RR会产生多种产物——NH₃（目标产物，8e⁻）、N₂（10e⁻）、N₂O（8e⁻）以及寄生的H₂（HER，2e⁻）——产物间的选择性机制尚不清楚。

### 1.2 The critical branching point: *NO / 关键分支点：*NO

The surface intermediate *NO sits at the crossroads of the entire reaction network:

表面中间体*NO处于整个反应网络的十字路口：

- **Hydrogenation pathway / 加氢路径**: *NO + H⁺ + e⁻ → *NOH → ... → NH₃ (Eₐ ≈ 0.08 eV)
- **Coupling pathway / 偶联路径**: *NO + *NO(adj) → *N₂O → N₂/N₂O (Eₐ ≈ 0.05–0.10 eV)

The activation barriers for these two pathways are nearly identical. Therefore, the selectivity is **not** determined by kinetic barriers alone, but by the **spatial distribution** of *NO on the surface — whether *NO molecules are dispersed (favoring hydrogenation → NH₃) or clustered (favoring N-N coupling → N₂).

两条路径的活化能几乎相同。因此，选择性**并非**由动力学势垒单独决定，而是由*NO在表面上的**空间分布**决定——*NO分子是分散分布（有利于加氢→NH₃）还是岛状聚集（有利于N-N偶联→N₂）。

### 1.3 Why mean-field kinetics fails here / 为什么平均场动力学在此失效

Conventional mean-field microkinetic models (MKM) assume spatially uniform adsorbate distributions: r = k·θ_A·θ_B. This assumption breaks down when:

传统的平均场微动力学模型(MKM)假设吸附物空间均匀分布：r = k·θ_A·θ_B。这一假设在以下条件下失效：

1. **Strong lateral interactions** drive spatial ordering (islands, depletion zones) / **强横向相互作用**驱动空间有序化（岛状、空缺带）
2. **Bimolecular surface reactions** (e.g., *NO + *NO) require adjacent sites; r ∝ θ² overestimates the rate when adsorbates are clustered / **双分子表面反应**（如*NO + *NO）需要相邻位点；当吸附物聚集时，r ∝ θ²高估反应速率
3. **Finite diffusion** prevents equilibrium mixing / **有限扩散**阻止平衡混合

Prior KMC studies in other systems have shown MKM errors of **orders of magnitude** (e.g., CO oxidation on RuO₂) and **2.3× coverage deviation** in HER on Pt (our Phase 1 result).

前人对其他体系的KMC研究表明MKM误差可达**数个数量级**（如RuO₂上CO氧化），以及HER on Pt的**2.3倍覆盖度偏差**（我们Phase 1的结果）。

### 1.4 Literature gap / 文献空白

As of March 2026, **no kinetic Monte Carlo study of NO₃RR exists**. All existing NO₃RR selectivity models rely on mean-field MKM [refs: ACS Catal. 2019; Faraday Discuss. 2023]. This is the first study to apply spatial KMC with ab initio parameters to NO₃RR.

截至2026年3月，**NO₃RR领域不存在任何动力学蒙特卡洛研究**。现有NO₃RR选择性模型全部依赖平均场MKM。本研究是首个将基于从头算参数的空间KMC应用于NO₃RR的工作。

---

## 2. Core Scientific Question / 核心科学问题

> **Does the spatial distribution of *NO on Cu(100), governed by lateral interactions, surface diffusion, and *H competition, determine the NH₃ vs N₂ selectivity in electrochemical NO₃RR?**

> **Cu(100)上*NO的空间分布（由横向相互作用、表面扩散和*H竞争共同控制）是否决定了电化学NO₃RR中NH₃与N₂的选择性？**

This question is:
- **Experimentally accessible but mechanistically unresolved** — high *NO coverage correlates with N₂ selectivity, but the causal mechanism is unknown / **实验可观察但机理不明**——高*NO覆盖度与N₂选择性正相关，但因果机制未知
- **Impossible to answer with mean-field models** — MKM assumes θ_NO is spatially uniform / **平均场模型无法回答**——MKM假设θ_NO空间均匀
- **A unique capability of spatial KMC** — the only method that tracks individual adsorbate positions and their evolution / **空间KMC的独有能力**——唯一追踪单个吸附物位置及其演化的方法

---

## 3. Research Approach / 研究方法

### 3.1 Overall strategy / 总体策略

We develop **SPARK**, a general-purpose electrochemical KMC framework (Python + Rust), and apply it to NO₃RR on Cu(100) with a fully ab initio parameter set from our own GC-DFT calculations.

我们开发**SPARK**——通用电催化KMC框架（Python + Rust），将其应用于Cu(100)上的NO₃RR，所有参数均来自我们自己的GC-DFT从头算计算。

### 3.2 Ab initio DFT parameters / 从头算DFT参数

All kinetic parameters are obtained from our own DFT calculations — no literature values are borrowed.

所有动力学参数均来自我们自己的DFT计算——不借用文献数值。

**DFT methodology / DFT方法学**:

| Item / 项目 | Choice / 选择 |
|---|---|
| Surface / 表面 | Cu(100), p(4×4) supercell, 4-layer slab (64 Cu atoms) |
| Functional / 泛函 | PBE + DFT-D3(BJ) dispersion correction |
| Smearing / 展宽 | Gaussian (ISMEAR=0), σ = 0.05 eV |
| Transition states / 过渡态 | CI-NEB for all steps (thermal + PCET thermodynamic analog) |
| Barrier calculation / 能垒计算 | **GC-DFT** (VASPsol implicit solvation + variable charge + grand potential Ω(U) quadratic fit) |
| PCET fallback / PCET备选 | BEP relation (Eₐ = α·ΔG + Eₐ₀, α ≈ 0.5) if CI-NEB fails to converge |

**Calculation scope / 计算规模** (~390–400 VASP jobs):

| Stage / 阶段 | Content / 内容 | Count / 数量 |
|---|---|---|
| Structure optimization / 结构优化 | Bulk, gas phase, adsorbates, co-adsorption | ~57 |
| CI-NEB transition states / CI-NEB过渡态 | 14 elementary steps × 5 images | ~14 NEB jobs |
| Frequency calculations / 频率计算 | IS + TS for all reactions | ~30 |
| GC-DFT constant-potential / GC-DFT恒电势 | 15 steps × 7 charges × 2 (IS+TS) | ~210 |

**Lateral interactions / 横向相互作用**: Pairwise interaction energies ε(A,B,d) are computed for 10 adsorbate pairs at 1NN and 2NN distances:

成对相互作用能ε(A,B,d)在1NN和2NN距离处对10对吸附物进行计算：

```
ε(A,B,d) = E(slab+A+B@d) − E(slab+A) − E(slab+B) + E(clean slab)
```

Tier 1 (critical / 关键): *NO–*NO, *NO–*H, *H–*H, *OH–*NO
Tier 2 (recommended / 建议): *OH–*OH, *OH–*H, *N–*NO, *N–*N, *O–*NO, *N–*H

### 3.3 KMC model / KMC模型

**Reaction network / 反应网络** (~20 elementary steps):

```
NO₃⁻(aq) ─→ *NO₃ ─(PCET)→ *NO₂ ─(PCET)→ *NO ─┬─(PCET)→ *NOH → *NHOH → *NH → *NH₂ → *NH₃(aq)  [NH₃产物]
                                                  └─(thermal)→ *N₂O → N₂O(g) / N₂(g)              [N₂/N₂O产物]
HER competing / HER竞争: H₂O + e⁻ → *H → H₂                                                       [H₂副产物]
```

**Rate expressions / 速率表达式**:

- PCET steps / PCET步骤: k(U) = (k_BT/h)·exp(−Eₐ(U)/k_BT), where Eₐ(U) from GC-DFT quadratic fit / 其中Eₐ(U)由GC-DFT二次拟合得到
- Thermal steps / 热力学步骤: k = (k_BT/h)·exp(−Eₐ/k_BT)
- Lateral interaction correction / 横向相互作用修正: k(site) = k_base × exp(+E_lat/k_BT)

**Simulation parameters / 模拟参数**: 50×50 square lattice (Cu(100)), T = 298 K, 10⁶ equilibration + 5×10⁶ production steps, 5–10 replicas with different random seeds.

### 3.4 Four systematic studies / 四项系统研究

#### Study A: *NO spatial distribution vs lateral interaction strength / *NO空间分布 vs 横向相互作用强度

- **Question / 问题**: How does ε(*NO–*NO) control the spatial pattern of *NO?
  ε(*NO–*NO)如何控制*NO的空间分布模式？
- **Method / 方法**: Fix U = −0.6 V, sweep ε(*NO–*NO) from −0.10 to +0.20 eV
  固定U = −0.6 V，扫描ε(*NO–*NO)从−0.10到+0.20 eV
- **Observables / 观测量**: Coverage θ_NO, radial distribution function g(r), cluster size distribution, FE(N₂)/FE(NH₃) vs ε
  覆盖度θ_NO、径向分布函数g(r)、团簇尺寸分布、FE(N₂)/FE(NH₃) vs ε
- **Expected result / 预期结果**: Attractive ε → island clustering → enhanced N-N coupling → more N₂; Repulsive ε → dispersed *NO → more NH₃
  吸引ε → 岛状聚集 → 增强N-N偶联 → 更多N₂；排斥ε → *NO分散 → 更多NH₃

#### Study B: NH₃/N₂ selectivity vs potential (KMC vs MKM) / NH₃/N₂选择性 vs 电位（KMC对比MKM）

- **Question / 问题**: At which potentials does mean-field fail? How large is the selectivity deviation?
  哪些电位下平均场失效？选择性偏差有多大？
- **Method / 方法**: Sweep U from −0.2 to −1.0 V, run both KMC and MKM with identical parameters
  扫描U从−0.2到−1.0 V，用相同参数同时运行KMC和MKM
- **Expected result / 预期结果**: KMC ≈ MKM at low overpotential; KMC ≠ MKM at high overpotential due to spatial phase separation
  低过电位KMC ≈ MKM；高过电位因空间相分离导致KMC ≠ MKM

#### Study C: *NO diffusion rate vs selectivity / *NO扩散速率对选择性的影响

- **Question / 问题**: Does the rate of *NO diffusion determine when KMC deviates from MKM?
  *NO的扩散速率是否决定KMC何时偏离MKM？
- **Method / 方法**: Fix U = −0.6 V, sweep diffusion barrier E_diff(NO) from 0.05 to 0.50 eV
  固定U = −0.6 V，扫描扩散势垒E_diff(NO)从0.05到0.50 eV
- **Expected result / 预期结果**: Fast diffusion → well-mixed → KMC ≈ MKM; Slow diffusion → spatial inhomogeneity → KMC ≠ MKM. Quantify via Damköhler number Da = k_reaction/k_diffusion.
  快扩散 → 充分混合 → KMC ≈ MKM；慢扩散 → 空间不均匀 → KMC ≠ MKM。用Damköhler数Da = k_reaction/k_diffusion定量化。

#### Study D: Polarization curves + Faradaic efficiency vs potential / 极化曲线 + 法拉第效率 vs 电位

- **Question / 问题**: Can SPARK quantitatively predict experimental observables?
  SPARK能否定量预测实验可观测量？
- **Method / 方法**: High-resolution polarization curves (17 potentials), compare with experimental FE(U) data
  高分辨极化曲线（17个电位点），与实验FE(U)数据对比
- **Observables / 观测量**: j_total(U), j_NH₃(U), j_N₂(U), j_H₂(U), Tafel slope, degree of rate control
  总电流密度j_total(U)、各产物偏电流密度、Tafel斜率、速率控制度

---

## 4. Expected Scientific Contributions / 预期科学贡献

| # | Contribution / 贡献 | Novelty / 新颖性 |
|---|---|---|
| 1 | **First KMC study of NO₃RR** / **首个NO₃RR的KMC研究** | Complete literature gap / 完全文献空白 |
| 2 | **Spatial effects quantitatively determine NH₃/N₂ selectivity** / **空间效应定量决定NH₃/N₂选择性** | Mean-field cannot address this / MKM无法回答 |
| 3 | **ε(*NO–*NO) as a selectivity tuning knob** / **ε(*NO–*NO)作为选择性调控旋钮** | Links DFT lateral interactions to macroscopic selectivity / 连接DFT横向相互作用与宏观选择性 |
| 4 | **Diffusion rate threshold for mean-field breakdown** / **平均场失效的扩散速率阈值** | Practical guide for when MKM is insufficient / 判断MKM何时不足的实用指南 |
| 5 | **First general-purpose electrochemical KMC framework** / **首个通用电催化KMC框架** | Software contribution / 软件贡献 |
| 6 | **Fully ab initio parameter set via GC-DFT** / **基于GC-DFT的完全从头算参数集** | Self-consistent, not literature patchwork / 自洽参数，非文献拼凑 |

---

## 5. Technical Innovation / 技术创新

### 5.1 SPARK framework / SPARK框架

Existing KMC software (Zacros, kmos, SPPARKS) are designed for thermal catalysis and lack native electrochemical support. SPARK provides:

现有KMC软件（Zacros、kmos、SPPARKS）均为热催化设计，不原生支持电化学。SPARK提供：

- **Butler-Volmer electrochemistry** — potential-dependent PCET rates / **Butler-Volmer电化学**——电位依赖的PCET速率
- **GC-DFT constant-potential barriers** — Eₐ(U) = aU² + bU + c from grand potential fitting / **GC-DFT恒电势能垒**——从巨势拟合得到Eₐ(U)
- **Native polarization curves** — j(U), FE(U) as direct simulation output / **原生极化曲线**——j(U)、FE(U)作为直接模拟输出
- **Pairwise lateral interactions** — site-dependent rate modification / **成对横向相互作用**——位点依赖的速率修正
- **Multi-product tracking** — TOF, partial current density, and Faradaic efficiency for each product / **多产物追踪**——每种产物的TOF、偏电流密度和法拉第效率
- **Dual engine** — Python prototype (~1,500 steps/s) + Rust high-performance engine (Fenwick tree O(log N) site selection) / **双引擎**——Python原型 + Rust高性能引擎（Fenwick树O(log N)位点选择）

### 5.2 GC-DFT parameterization workflow / GC-DFT参数化工作流

Unlike prior KMC studies that use CHE (computational hydrogen electrode) for PCET barriers, we employ the more rigorous **grand-canonical DFT (GC-DFT)** approach:

与此前使用CHE（计算氢电极）处理PCET能垒的KMC研究不同，我们采用更严格的**巨正则DFT (GC-DFT)** 方法：

1. Optimize IS and TS structures at neutral charge / 在中性电荷下优化IS和TS结构
2. Perform 7 variable-charge single-point calculations for each (NELECT ± 0.5/1.0/1.5 + neutral) with VASPsol implicit solvation / 对每个结构进行7个变电荷单点能计算（NELECT ± 0.5/1.0/1.5 + 中性），使用VASPsol隐式溶剂化
3. Extract work function → convert to electrode potential U vs SHE / 提取功函数 → 转换为电极电位U vs SHE
4. Fit grand potential Ω(U) = aU² + bU + c for both IS and TS / 对IS和TS分别拟合巨势Ω(U) = aU² + bU + c
5. Eₐ(U) = Ω_TS(U) − Ω_IS(U) + ΔZPE / Eₐ(U) = Ω_TS(U) − Ω_IS(U) + ΔZPE

This captures the potential-dependent charge reorganization that CHE ignores.

这捕捉了CHE忽略的电位依赖的电荷重组效应。

---

## 6. Research Significance / 研究意义

This work advances the field in three dimensions:

本研究从三个维度推进领域发展：

**Fundamental understanding / 基础理解**: Reveals, for the first time, how *NO spatial distribution — not just energetics — controls NH₃ vs N₂ selectivity in NO₃RR. Establishes ε(*NO–*NO) as the microscopic origin of macroscopic selectivity.

首次揭示*NO的空间分布（而非仅仅是能量学）如何控制NO₃RR中NH₃与N₂的选择性。将ε(*NO–*NO)确立为宏观选择性的微观起源。

**Methodology / 方法学**: Demonstrates a complete workflow from ab initio GC-DFT parameters → spatial KMC → experimentally comparable observables (FE, j, Tafel slope). Provides a quantitative criterion for when mean-field MKM is insufficient.

展示了从从头算GC-DFT参数 → 空间KMC → 实验可比较观测量（FE、j、Tafel斜率）的完整工作流。提供了判断平均场MKM何时不足的定量标准。

**Tool development / 工具开发**: SPARK fills the gap of a general-purpose electrochemical KMC framework. Its Python + Rust dual-engine architecture balances usability and performance, enabling community adoption.

SPARK填补了通用电催化KMC框架的空白。其Python + Rust双引擎架构兼顾易用性和性能，便于社区推广。

**Target journal / 目标期刊**: ACS Catalysis

**Extension / 扩展方向**: NO₃RR → Urea electrosynthesis (incremental development: +CO₂RR sub-pathway + C-N coupling → Nature Catalysis level)

NO₃RR → 尿素电合成（增量开发：+CO₂RR子路径 + C-N偶联 → Nature Catalysis级别）
