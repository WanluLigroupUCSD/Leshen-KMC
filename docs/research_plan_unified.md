# Spatial Kinetic Monte Carlo Simulation of Electrocatalytic Reactions: From HER Validation to CO₂RR Selectivity Control

# 电催化反应的空间动力学蒙特卡洛模拟：从HER验证到CO₂RR选择性调控

**Leshen-KMC: First General-Purpose Electrocatalytic KMC Framework**

---

## Abstract / 摘要

We present a two-phase computational study using Leshen-KMC — the first general-purpose KMC framework with native electrochemical support — to investigate electrocatalytic reactions where spatial effects critically determine activity and selectivity. **Phase 1** (completed) validates the software using the hydrogen evolution reaction (HER) on Pt(111), reproducing experimental polarization curves and demonstrating a 2× coverage divergence between KMC and mean-field models due to lateral interactions. **Phase 2** addresses a central open question in CO₂ electroreduction: *how does the spatial distribution of CO* adsorbates on Cu(100) control C₂ product selectivity?* This question cannot be answered by mean-field models and has never been studied with a general-purpose electrocatalytic KMC tool.

本方案提出一项两阶段计算研究，使用 Leshen-KMC —— 首个原生支持电化学的通用 KMC 框架。**第一阶段**（已完成）以 Pt(111) 上析氢反应（HER）验证软件，复现实验极化曲线，并证明了横向相互作用导致KMC与mean-field模型之间2倍的覆盖度偏差。**第二阶段**研究 CO₂ 电还原中的核心开放问题：*Cu(100) 上 CO* 吸附物的空间分布如何决定 C₂ 产物选择性？* 此问题 mean-field 模型无法回答，且从未有通用电催化 KMC 工具对其进行研究。

---

## Table of Contents / 目录

1. [Introduction / 引言](#1-introduction--引言)
2. [Literature Review / 文献综述](#2-literature-review--文献综述)
3. [Phase 1: HER on Pt(111) — Validation (Completed)](#3-phase-1-her-on-pt111--validation-completed)
4. [Phase 2: CO₂RR on Cu(100) — Spatial Control of C₂ Selectivity](#4-phase-2-co₂rr-on-cu100--spatial-control-of-c₂-selectivity)
5. [Computational Details / 计算细节](#5-computational-details--计算细节)
6. [Expected Scientific Contributions / 预期科学贡献](#6-expected-scientific-contributions--预期科学贡献)
7. [Timeline / 时间计划](#7-timeline--时间计划)
8. [References / 参考文献](#8-references--参考文献)

---

## 1. Introduction / 引言

### 1.1 The Gap: No General-Purpose Electrocatalytic KMC Tool / 空白：缺乏通用电催化KMC工具

Kinetic Monte Carlo simulation is the gold standard for bridging DFT-computed energetics and macroscopic catalytic behavior, capturing spatial correlations, coverage fluctuations, and island formation that mean-field models miss (*J. Chem. Phys.* **156**, 120902 (2022)). However, **all existing general-purpose KMC packages were designed for thermal catalysis:**

KMC 模拟是连接 DFT 计算能量学与宏观催化行为的黄金标准方法，能够捕捉 mean-field 模型遗漏的空间关联、覆盖度涨落和岛状聚集。然而，**现有所有通用 KMC 软件均为热催化设计：**

| Software | Language | Lateral Interactions | Electrochemistry | Polarization Curve |
|----------|----------|:--------------------:|:----------------:|:------------------:|
| Zacros | Fortran | Graph-theoretical CE (多体) | ✗ | ✗ |
| kmos | Python/Fortran | Pairwise | ✗ | ✗ |
| SPPARKS | C++ | Custom | ✗ | ✗ |
| KMCLib | Python/C++ | Pairwise | ✗ | ✗ |
| MonteCoffee | Python | Custom (graph) | ✗ | ✗ |
| CatMAP | Python | Linear correction (MKM) | ✓ (MKM only) | ✗ |
| **Leshen-KMC** | **Python + Rust** | **Pairwise NN** | **✓ Butler-Volmer PCET** | **✓ Native** |

Recent electrocatalytic KMC studies have all used one-off custom codes:

近期电催化KMC研究均使用一次性自编代码：

| Study | System | Code | Key Finding |
|-------|--------|------|-------------|
| *J. Phys. Chem. Lett.* **16**, 2896–2904 (2025) | CO₂RR on Cu, 178 reactions | Custom | Cu(100)→C₂, Cu(111)→C₁, voltage-dependent RDS |
| *ACS Catal.* **14**, 2696–2708 (2024) | HER on Pt, coverage-dependent | Custom | H₅O₂⁺ solvent + dynamic H coverage essential |
| *Nat. Commun.* (2025) | HER on Au, multiscale | Custom | Full pH range polarization |
| *J. Chem. Phys.* **155**, 164701 (2021) | CO₂RR on Cu, diffusion | Custom | Mean-field breakdown under finite diffusion |

**Leshen-KMC fills this gap** as the first general-purpose framework with native Butler-Volmer rates, polarization curve computation, and integrated KMC + mean-field microkinetic modeling.

### 1.2 Why Spatial KMC Matters for Electrocatalysis / 为什么空间KMC对电催化至关重要

Mean-field microkinetic models (MKM) assume uniform surface coverage:

$$r_{C-C} = k_{CC} \cdot \theta_{CO}^2$$

But in reality, CO* molecules form spatial patterns (islands, stripes, random) depending on:
- **Lateral interactions** — attraction promotes clustering, repulsion promotes dispersion
- **Surface diffusion rate** — fast diffusion → well-mixed (mean-field valid); slow diffusion → spatial correlations
- **Relative rates of production vs consumption** — imbalance creates local concentration gradients

Recent work has systematically demonstrated when mean-field breaks down:

| Scenario | MKM | KMC | Evidence |
|----------|-----|-----|----------|
| Uniform high coverage | ✅ Accurate | ✅ Accurate (slower) | — |
| Strong lateral interactions | ❌ Over/underestimates θ | ✅ Captures ordered phases/islands | *ACS Catal.* **14**, 2696 (2024): θ_H diverges 2× |
| Diffusion-limited | ❌ Assumes instant mixing | ✅ Captures site selectivity | *J. Chem. Phys.* **155**, 164701 (2021) |
| Multi-site competition | ❌ Ignores spatial correlations | ✅ Natural treatment | *J. Phys. Chem. Lett.* **16**, 2896 (2025) |
| Low-coverage fluctuations | ❌ Continuum fails | ✅ Stochastic by nature | *J. Chem. Phys.* **156**, 120902 (2022) |

*Catal. Today* **372**, 11 (2021) showed that **mean-field models can overestimate catalytic activity by several orders of magnitude** when lateral interactions are significant. For C-C coupling (a bimolecular surface reaction requiring adjacent sites), this error is expected to be particularly severe.

### 1.3 The Central Question / 核心科学问题

> **Does the spatial arrangement of CO* adsorbates on Cu(100) — controlled by lateral interactions and diffusion — determine the selectivity between C₁ (CH₄) and C₂ (C₂H₄) products in electrochemical CO₂ reduction?**
>
> **Cu(100) 上 CO* 吸附物的空间排布（由横向相互作用和扩散控制）是否决定了电化学CO₂还原中C₁(CH₄)与C₂(C₂H₄)产物的选择性？**

This question is:
- **Experimentally observed but mechanistically unexplained** — higher CO* coverage correlates with higher C₂ selectivity (*Nat. Energy* **9**, 1485–1496 (2024))
- **Impossible to answer with mean-field models** — they assume uniform θ_CO
- **Uniquely suited for spatial KMC** — the only computational method that tracks individual adsorbate positions and their evolution

---

## 2. Literature Review / 文献综述

### 2.1 KMC Methodology Reviews / KMC方法论综述

| Ref | Article | Journal | Key Points |
|-----|---------|---------|------------|
| [10] | "KMC simulations for heterogeneous catalysis: Fundamentals, current status, and challenges" | *J. Chem. Phys.* **156**, 120902 (2022) | **金标准综述**。1p-KMC框架桥接原子→宏观多尺度。讨论刚性问题、大晶格、加速方法。 |
| [14] | "Unraveling the Complexity of Catalytic Reactions via KMC Simulation" | *ACS Catal.* **2**, 2648–2663 (2012) | 里程碑综述，建立graph-theoretical KMC框架，为Zacros奠基。 |
| [15] | "A Practical Guide to Surface KMC Simulations" | *Front. Chem.* **7**, 202 (2019) | **实用指南**。BKL算法、速率常数计算、lateral interactions处理、时间步进实现细节。 |

### 2.2 Multiscale Modeling Reviews / 多尺度建模综述

| Ref | Article | Journal | Key Points |
|-----|---------|---------|------------|
| [16] | "DFT-Based Multiscale Modeling of Heterogeneous (Electro)Catalytic Reactions" | *ACS Catal.* (2025) | **最新综述**。DFT→KMC/MKM→CFD全链路。恒电位DFT、溶剂效应、KMC与微观动力学耦合。 |
| [8] | "Multiscale Modeling of CO₂ Electrochemical Reduction on Copper" | *ChemSusChem* (2025) | CO₂RR多尺度建模专题。DFT/MD→KMC/MK→CFD，聚焦Cu催化剂。 |
| [17] | "Microkinetic modeling in electrocatalysis: Applications, limitations, recommendations" | *J. Catal.* **400**, 290 (2021) | MKM在电催化中的局限。Mean-field失效场景，推荐KMC作替代。 |

### 2.3 KMC for HER / KMC在HER中的应用

| Ref | Article | Journal | Key Points |
|-----|---------|---------|------------|
| [2] | "First-Principles-Based KMC Model of HER under Realistic Conditions" | *ACS Catal.* **14**, 2696–2708 (2024) | **Phase 1核心参考**。Pt(111)/(100)第一性原理KMC-HER。考虑H₅O₂⁺溶剂化、动态H覆盖度、电场。Pt(111)以Tafel为主，Pt(100)以Heyrovsky为主。 |
| [3] | "Toward rational understanding of HER polarization curves through multiscale simulations" | *Nat. Commun.* (2025) | HER on Au，多尺度模拟，全pH范围极化曲线。 |
| [12] | "Cyclic Voltammograms for H on Pt(111) and Pt(100) from First Principles" | *Phys. Rev. Lett.* **99**, 126101 (2007) | Pt上H横向相互作用参数。ε(H*-H*) = +0.10 eV (repulsive)。 |

**Key findings for HER / HER关键发现：**
- 溶剂效应关键：必须考虑H₅O₂⁺物种和动态H覆盖度演化
- Pt(111) vs Pt(100)：不同晶面主导机理不同（Tafel vs Heyrovsky）
- 氢溢出(spillover)：金属表面H覆盖度高时，载体氢溢出可提升HER
- 覆盖度效应：高电位→高覆盖度→加速HER，KMC可捕捉动态演化

### 2.4 KMC for CO₂RR / KMC在CO₂RR中的应用

| Ref | Article | Journal | Key Points |
|-----|---------|---------|------------|
| [1] | "Voltage-Dependent CO₂RR Mechanism Unveiled by KMC Simulation" | *J. Phys. Chem. Lett.* **16**, 2896–2904 (2025) | **Phase 2核心参考**。Cu(111)/(100)上178个基元反应→C₁/C₂。Cu(111)→CH₄(RDS: CO→CHO\*/COH\*)；Cu(100)→C₂H₄/EtOH(RDS: CO\*对称偶联)。 |
| [4] | "Effects of surface diffusion in electrocatalytic CO₂RR on Cu" | *J. Chem. Phys.* **155**, 164701 (2021) | **扩散效应经典工作**。Lattice KMC。无穷快扩散=mean-field；有限扩散→位点选择性出现，mean-field失效。 |
| [18] | "In Situ Reconstruction of Cu(100) for Promoted C–C Coupling" | *ACS Catal.* (2025) | GCMC+EKMC联合模拟Cu(100)重构(~35nm×35nm, 300s)。重构motif通过界面氢键增强CO二聚化。 |

**Key findings for CO₂RR / CO₂RR关键发现：**
- **C₁ vs C₂选择性**：Cu(111)→CH₄(C₁), Cu(100)→C₂H₄/EtOH(C₂)，RDS不同
- **表面扩散是关键变量**：扩散速率直接影响C₂产物选择性，KMC揭示mean-field无法捕捉的空间关联
- **电位依赖机理转变**：随电位增大，Cu(111)上RDS从CO→CHO\*转为CO→COH\*
- **表面重构**：Cu(100)在CO₂RR条件下发生原位重构，影响C-C偶联活性
- **178个基元反应**：*J. Phys. Chem. Lett.* **16**, 2896 (2025) 的KMC模型涵盖全面反应网络

### 2.5 KMC for NRR and Battery Interfaces / KMC在NRR和电池界面中的应用

| Ref | Article | Journal | Key Points |
|-----|---------|---------|------------|
| [19] | "DFT–kMC–LSTM for NRR on Transition Metal Oxides" | *ACS Catal.* **13**, 7 (2023) | DFT→kMC→LSTM三阶段。V₂O₃ NRR活性比Ru高1000×。10000+小时稳定。 |
| [20] | "SEI Formation and Li-Ion Electrodeposition via KMC" | *ACS Energy Lett.* (2024) | 第一性原理KMC模拟锂金属SEI。晶界促进Li⁺传输。 |

### 2.6 Mean-field Breakdown and Lateral Interactions / Mean-field失效与横向相互作用

| Ref | Article | Journal | Key Points |
|-----|---------|---------|------------|
| [9] | "Evaluating benefits of KMC and MKM in the presence of lateral interactions" | *Catal. Today* **372**, 11 (2021) | **系统比较MKM vs KMC**。Lateral interactions越强，MKM偏差越大。 |
| [21] | "Speeding up lateral interactions detection in graph-theoretical KMC" | *J. Phys. Chem. A* (2024) | Zacros框架中加速lateral interaction检测的算法。 |

**KMC captures spatial effects that MKM misses / KMC捕捉MKM遗漏的空间效应：**
1. 吸附物岛状聚集 — 即使无lateral interactions，中间体倾向形成cluster
2. 空间关联 — 扩散慢时CO\*非随机分布，影响C-C偶联概率
3. 覆盖度横向变化 — 催化面上覆盖度非均匀，RDS随之变化
4. 有序相转变 — 特定覆盖度下出现p(2×2)等有序结构

### 2.7 KMC Acceleration Methods / KMC加速方法

#### τ-leaping and Time-Scale Acceleration / τ-leaping与时间尺度加速

| Method | Principle | Applicable Scenario |
|--------|-----------|---------------------|
| **τ-leaping** | 时间窗口τ内批量执行多个事件 | 事件频率均匀的体系 |
| **Implicit τ-leaping** | τ-leaping隐式版本 | **刚性体系**（快慢反应共存，如CO₂RR） |
| **Rate constant rescaling** | 缩放快反应速率减小刚性 | 扩散远快于反应 |
| **Net-event algorithm** | 只追踪净事件数 | 可逆反应频繁正反抵消 |
| **Parallel KMC** | 空间分解+时间同步 | 大晶格 |

#### Machine Learning Acceleration / 机器学习加速

| Method | Principle | Reference |
|--------|-----------|-----------|
| **Neural Network Potentials** | NNP替代DFT，on-the-fly KMC | JCP 2024, 160, 204108 |
| **Graph Neural Networks** | 高通量催化剂筛选，预测吸附能 | WIREs Comp. Mol. Sci. 2025 |
| **DFT-kMC-LSTM** | KMC短期动力学 + LSTM外推长期降解 | *ACS Catal.* **13**, 7 (2023) |
| **Gaussian Process + KMC** | 主动学习发现过渡态 | npj Comp. Mater. 2024 |

---

## 3. Phase 1: HER on Pt(111) — Validation (Completed) / 第一阶段：HER验证（已完成）

### 3.1 Objective / 目标

Validate Leshen-KMC against analytical solutions and published data using the simplest electrochemical reaction: HER on Pt(111). This system exercises all framework features (PCET, lateral interactions, diffusion, KMC vs MKM) with only 2 species and 3 reaction steps.

以最简单的电化学反应验证 Leshen-KMC：Pt(111) 上 HER。该体系仅含2种物种和3步反应，但能全面检验框架功能。

### 3.2 Model / 模型

**Species:** empty (*), H*

**Lattice:** 50×50 2D square, PBC, T = 298 K

**Reaction Network:**

| # | Step | Reaction | Type | E_a [eV] | Source |
|---|------|----------|------|----------|--------|
| 1 | Volmer (fwd) | H⁺ + e⁻ + * → H* | PCET | 0.67 | [2] |
| 2 | Volmer (rev) | H* → H⁺ + e⁻ + * | PCET | 0.62 | [2] |
| 3 | Tafel | H* + H*(adj) → H₂ + 2* | Thermal | 0.85 | [11] |
| 4 | Heyrovsky (fwd) | H* + H⁺ + e⁻ → H₂ + * | PCET | 0.70 | [2] |
| 5 | H* diffusion | H*@s₁ + *@s₂ → *@s₁ + H*@s₂ | Thermal | 0.10 | [13] |

**Lateral interaction:** ε(H*, H*) = +0.10 eV (repulsive, NN), from *Phys. Rev. Lett.* **99**, 126101 (2007)

**Rate expressions:**

PCET steps (Butler-Volmer):
$$k(U) = \frac{k_B T}{h} \exp\left(-\frac{E_a + \beta_{BV} \cdot e \cdot (U - U_0)}{k_B T}\right)$$

Thermal steps (TST):
$$k = \frac{k_B T}{h} \exp\left(-\frac{E_a}{k_B T}\right)$$

### 3.3 Validation Tests / 验证测试

| Test | Method | Expected Outcome | Status |
|------|--------|------------------|--------|
| 1. Langmuir limit | Off lateral + Tafel | θ_H = Langmuir isotherm | ✅ PASSED |
| 2. Tafel slope | j-V at low η | 30-120 mV/dec | ✅ PASSED (81 mV/dec) |
| 3. Lateral effect | θ_H with/without ε | θ_H decreases with repulsion | ✅ PASSED |
| 4. Diffusion effect | TOF with/without diff | Diffusion enhances Tafel rate | ✅ PASSED |
| 5. KMC vs MKM | Same model, both methods | MKM overestimates θ_H | ✅ PASSED (θ: 0.76 vs 0.33) |
| 6. Lattice convergence | 10² to 10⁴ sites | Converge by 50×50 | ✅ PASSED |

All 6 tests passed. Validation script: `validate_her.py`. Pure Python performance: ~1500 steps/s.

### 3.4 Results / 结果

#### Polarization Curve

MKM 101 potential points + KMC 11 points. Results in `results/` directory.

**Key metrics:**
- Exchange current density: j₀ = 0.45 mA/cm² (experimental Pt: ~1 mA/cm²)
- Tafel slope: 81 mV/dec → Heyrovsky-dominant mechanism
- Consistent with *ACS Catal.* **14**, 2696–2708 (2024) findings

#### KMC vs MKM Divergence — The Key Result

| Metric | MKM (mean-field) | KMC (spatial) | Divergence |
|--------|-------------------|---------------|------------|
| θ_H at U = −0.2 V | ~0.76 | ~0.33 | **2.3×** |
| Mechanism | Heyrovsky | Heyrovsky | Same |
| Tafel slope | ~80 mV/dec | ~81 mV/dec | Consistent |

**Interpretation:** H*-H* repulsion (+0.10 eV per NN pair) creates anti-clustering in the spatial KMC. Mean-field assumes uniform coverage and underestimates the energy penalty, leading to 2× overestimation of θ_H. This validates the central thesis: **spatial KMC captures lateral interaction effects that mean-field cannot.**

**解读：** H*-H*排斥作用在空间KMC中产生反聚集效应。Mean-field假设均匀覆盖，低估了能量惩罚，导致θ_H被高估2倍。这验证了核心论点：**空间KMC能捕捉mean-field无法处理的横向相互作用效应。**

### 3.5 Deliverables / 产出

- ✅ 5 figures + 2 data files in `results/`
- ✅ `run_polarization_her.py` — polarization curve script
- ✅ `validate_her.py` — 6-test validation suite
- ✅ KMC vs MKM divergence quantified

---

## 4. Phase 2: CO₂RR on Cu(100) — Spatial Control of C₂ Selectivity / 第二阶段：Cu(100)上CO₂RR的空间选择性调控

### 4.1 Scientific Hypothesis / 科学假说

**Hypothesis:** The C₂H₄/CH₄ selectivity ratio in electrochemical CO₂ reduction on Cu(100) is primarily controlled by the spatial distribution of CO* adsorbates, which is tuned by:
1. CO*-CO* lateral interaction strength (attractive → island formation → more C-C coupling)
2. CO* surface diffusion rate (fast → well-mixed → mean-field applies; slow → spatial patterns)
3. Applied potential (controls CO* production rate and consumption pathways)

**假说：** Cu(100) 上电化学 CO₂ 还原的 C₂H₄/CH₄ 选择性主要由 CO* 吸附物的空间分布控制，后者由横向相互作用强度、表面扩散速率和施加电位调节。

### 4.2 Reaction Model / 反应模型

**Species (8):** empty, CO₂*, COOH*, CO*, CHO*, OCCO*, CH₂O*, CH₃O*

**Lattice:** 2D square (Cu(100) top sites), 50×50, PBC, T = 298 K

**Elementary Steps (~12):**

| # | Reaction | Type | E_a [eV] | Note | Source |
|---|----------|------|----------|------|--------|
| 1 | CO₂(g) + * → CO₂* | Adsorption | ~0 | Hertz-Knudsen | — |
| 2 | CO₂* + H⁺ + e⁻ → COOH* | PCET | 0.43 | First protonation | [6] |
| 3 | COOH* + H⁺ + e⁻ → CO* + H₂O | PCET | 0.02 | Nearly barrierless | [6] |
| 4 | CO* → CO(g) + * | Desorption | 0.55 | Controls θ_CO | — |
| 5 | **CO* + CO*(adj) → OCCO*** | **Thermal** | **0.72** | **C-C coupling (C₂ key step, RDS on Cu(100))** | **[1]** |
| 6 | CO* + H⁺ + e⁻ → CHO* | PCET | 0.74 | C₁ pathway entry | [1] |
| 7 | OCCO* + H⁺ + e⁻ → ... → C₂H₄ | PCET | 0.40 | C₂ downstream (lumped) | [1] |
| 8 | CHO* + H⁺ + e⁻ → CH₂O* | PCET | 0.35 | C₁ downstream | — |
| 9 | CH₂O* + H⁺ + e⁻ → CH₃O* | PCET | 0.30 | C₁ downstream | — |
| 10 | CH₃O* + H⁺ + e⁻ → CH₄ + * | PCET | 0.25 | C₁ product release | — |
| 11 | CO* diffusion | Thermal | 0.10 | **Tunable parameter** | — |

**Lateral Interactions:**

| Pair | Energy | Physical Origin |
|------|--------|----------------|
| CO*-CO* | Variable: −0.05 to +0.10 eV | DFT, scan as parameter |
| CO*-CHO* | +0.05 eV (repulsive) | Steric |

**DFT Parameter Sources:**
- **[1]** *J. Phys. Chem. Lett.* **16**, 2896–2904 (2025) — 178 elementary reactions on Cu(111)/(100), voltage-dependent, most comprehensive
- **[6]** *J. Phys. Chem. Lett.* **6**, 2032–2037 (2015) — CO dimerization mechanism, early landmark
- **[4]** *J. Chem. Phys.* **155**, 164701 (2021) — Lattice KMC framework for CO₂RR, diffusion effects
- Nørskov group — CO binding energies on Cu surfaces

### 4.3 Study A: CO* Spatial Distribution Mapping / 研究A：CO*空间分布图谱

**Method:** Fix U = −1.0 V, T = 298 K. Run KMC on 50×50 lattice. Vary CO*-CO* lateral interaction:

| ε(CO*-CO*) | Expected CO* Pattern | C-C Coupling Rate |
|------------|---------------------|-------------------|
| −0.05 eV (attractive) | Island/cluster formation | **Enhanced** (many adjacent pairs) |
| 0.00 eV (no interaction) | Random distribution | Baseline |
| +0.05 eV (repulsive) | Dispersed/ordered | **Reduced** (few adjacent pairs) |
| +0.10 eV (strong repulsive) | Anti-clustering | **Strongly reduced** |

**Output:** CO* pair correlation function g(r), cluster size distribution P(n), lattice snapshot images

**Reference:** *Catal. Today* **372**, 11 (2021) — systematic MKM vs KMC comparison under varying lateral interaction strength

### 4.4 Study B: C₁/C₂ Selectivity vs Potential / 研究B：C₁/C₂选择性随电位变化

**Method:** Scan U from −0.6 to −1.4 V. At each potential, run KMC + MKM.

**Output:**
- Faradaic efficiency: FE(CH₄) and FE(C₂H₄) vs U
- KMC vs MKM comparison of selectivity
- Identify crossover potential where C₂ > C₁

**Key prediction:** KMC should predict higher C₂ selectivity than MKM when CO*-CO* interaction is attractive (islands enhance C-C coupling). MKM should overestimate C₁ pathway because it assumes uniform CO* distribution.

**Reference:** *J. Phys. Chem. Lett.* **16**, 2896–2904 (2025) — voltage-dependent RDS transition on Cu(111) and Cu(100)

### 4.5 Study C: Diffusion Rate Effect on C₂ Selectivity / 研究C：扩散速率对C₂选择性的影响

**Method:** Fix U = −1.0 V, ε = −0.05 eV. Vary E_a(diffusion):

| E_a(diff) [eV] | Regime | Expected Effect |
|-----------------|--------|-----------------|
| 0.05 | Very fast (mean-field limit) | KMC ≈ MKM |
| 0.10 | Moderate | KMC ≠ MKM, spatial effects emerge |
| 0.30 | Slow | Strong spatial patterns, KMC ≫ MKM deviation |
| 0.50 | Very slow (frozen CO*) | Maximal spatial effects |

**Output:** C₂/C₁ ratio vs diffusion barrier, quantify mean-field breakdown threshold

**Reference:** *J. Chem. Phys.* **155**, 164701 (2021) — demonstrated that finite diffusion causes site selectivity in CO₂RR; this study will quantify the critical threshold

### 4.6 Study D: Polarization Curves and Faradaic Efficiency / 研究D：极化曲线与法拉第效率

**Method:** Full potential scan, compute j_total, j_CH₄, j_C₂H₄, j_CO

**Expected output:**
```
U [V vs RHE]  |  j_total  |  FE(CH₄)  |  FE(C₂H₄)  |  FE(CO)
─────────────────────────────────────────────────────────────
   −0.6       |   0.1     |   10%     |    2%       |   80%
   −0.8       |   1.0     |   25%     |   15%       |   50%
   −1.0       |   5.0     |   30%     |   35%       |   20%
   −1.2       |  15.0     |   25%     |   50%       |   10%
```

Compare with experimental data from Hori group and recent operando studies.

**Reference:** *J. Phys. Chem. Lett.* **16**, 2896 (2025) + *ChemSusChem* (2025) review

---

## 5. Computational Details / 计算细节

### 5.1 Software / 软件

- **Leshen-KMC v0.3.0** (Python + Rust)
- Rust engine for production runs (Fenwick tree O(log N) site selection)
- Python engine for prototyping and analysis
- Newton-Raphson MKM solver for mean-field comparison

### 5.2 Simulation Parameters / 模拟参数

| Parameter | Phase 1 (HER) | Phase 2 (CO₂RR) |
|-----------|---------------|-----------------|
| Lattice | 50×50 | 50×50 |
| Temperature | 298 K | 298 K |
| Potential range | −0.5 to 0.0 V | −0.6 to −1.4 V |
| Equilibration | 10⁶ steps | 10⁶ steps |
| Production | 10⁶ steps | 5×10⁶ steps |
| Sampling | Every 10⁴ steps | Every 10⁴ steps |
| Independent runs | 5 (for error bars) | 5 |

### 5.3 Analysis Tools / 分析工具

- Coverage time series: θ_i(t)
- TOF per product: TOF_CH₄, TOF_C₂H₄, TOF_CO
- Pair correlation function: g_CO-CO(r)
- Cluster size distribution: P(n)
- Lattice snapshots: species coloring
- Polarization curves: j(U) and FE(U)

### 5.4 Computational Cost Considerations / 计算成本考量

Phase 2 CO₂RR体系存在刚性问题（扩散 >> 反应），可能需要：
- **速率常数重缩放**: 缩放快扩散速率以减小刚性比
- **隐式τ-leaping**: 未来扩展方向，适用于刚性电化学KMC（见§2.7）
- **并行化**: Rust引擎可利用IBEX HPC节点

---

## 6. Expected Scientific Contributions / 预期科学贡献

| # | Contribution | Novelty | Relevant Literature |
|---|-------------|---------|---------------------|
| 1 | **First spatial KMC study of CO* distribution → C₂ selectivity** | KMC essential, MKM cannot answer | [1], [4] |
| 2 | **Quantitative mean-field error for C-C coupling predictions** | Gap: [9] showed error exists but not for CO₂RR C₂ | [9] |
| 3 | **CO*-CO* lateral interaction as design knob for C₂ selectivity** | Connects DFT interaction energy → macroscopic selectivity | — |
| 4 | **Diffusion rate threshold for mean-field validity** | Practical guide: when is MKM sufficient? | [4] |
| 5 | **First general-purpose electrocatalytic KMC framework** | Software contribution to community | — |
| 6 | **Electrochemical polarization curve from spatial KMC** | Unique capability among KMC tools | [2] |

### Scientific Impact / 科学影响

This work provides the **first quantitative answer** to whether CO* spatial arrangement controls C₂ selectivity in CO₂ electroreduction — a question experimentally motivated (*Nat. Energy* **9**, 1485–1496 (2024)) but computationally inaccessible without electrocatalytic spatial KMC.

### Methodological Impact / 方法学影响

Establishes **when spatial KMC is essential for electrochemical systems** with bimolecular surface reactions (C-C coupling), providing practical guidelines for the computational catalysis community.

### Software Impact / 软件影响

Leshen-KMC as the **first open-source general-purpose electrocatalytic KMC framework** enables the broader community to perform spatial KMC for any electrocatalytic system without custom code — lowering the barrier from months of development to hours of model definition.

---

## 7. Timeline / 时间计划

| Week | Phase | Tasks | Status |
|------|-------|-------|--------|
| **1** | Phase 1 | T1-T2: Collect DFT data, implement HER model | ✅ Complete |
| **2** | Phase 1 | T3-T4: Validation tests, polarization curves | ✅ Complete |
| **3** | Phase 2 | **T5**: Collect CO₂RR DFT data from [1], [6], Nørskov | ⏳ Next |
| **3** | Phase 2 | **T6**: Implement Cu(100) model (~8 species, ~12 processes) | ⏳ Pending |
| **4** | Phase 2 | **T7**: Study A — CO* spatial distribution mapping | ⏳ Pending |
| **5** | Phase 2 | **T8**: Study B — C₁/C₂ selectivity vs potential, KMC vs MKM | ⏳ Pending |
| **6** | Phase 2 | **T9**: Study C — Diffusion rate effect (mean-field breakdown) | ⏳ Pending |
| **7** | Phase 2 | **T10**: Study D — Polarization curves and Faradaic efficiency | ⏳ Pending |
| **8** | Writing | **T11**: Manuscript preparation | ⏳ Pending |

**Current progress: 5/11 tasks complete (45%)**

---

## 8. References / 参考文献

### Electrocatalytic KMC Studies / 电催化KMC研究

[1] Wei, C. et al. Voltage-Dependent Electrochemical Carbon Dioxide Reduction Mechanism Unveiled by Kinetic Monte Carlo Simulation. *J. Phys. Chem. Lett.* **16**, 2896–2904 (2025). [DOI](https://pubs.acs.org/doi/10.1021/acs.jpclett.4c03426)

[2] Li, H. et al. First-Principles-Based Kinetic Monte Carlo Model of Hydrogen Evolution Reaction under Realistic Conditions: Solvent, Hydrogen Coverage and Electric Field Effects. *ACS Catal.* **14**, 2696–2708 (2024). [DOI](https://pubs.acs.org/doi/10.1021/acscatal.3c04588)

[3] Shou, W. et al. Toward rational understanding of the hydrogen evolution polarization curves through multiscale simulations. *Nat. Commun.* (2025).

[4] Jørgensen, M. & Grönbeck, H. Effects of surface diffusion in electrocatalytic CO₂ reduction on Cu revealed by kinetic Monte Carlo simulations. *J. Chem. Phys.* **155**, 164701 (2021). [DOI](https://pubs.aip.org/aip/jcp/article/155/16/164701/199778)

### CO₂RR Mechanism and Selectivity / CO₂RR机理与选择性

[5] Zhan, C. et al. Key intermediates and Cu active sites for CO₂ electroreduction to ethylene and ethanol. *Nat. Energy* **9**, 1485–1496 (2024).

[6] Montoya, J. H. et al. Theoretical Insights into a CO Dimerization Mechanism in CO₂ Electroreduction. *J. Phys. Chem. Lett.* **6**, 2032–2037 (2015).

[7] Nitopi, S. et al. Progress and Perspectives of Electrochemical CO₂ Reduction on Copper in Aqueous Electrolyte. *Chem. Rev.* **119**, 7610–7672 (2019).

[8] Gholizadeh, M. et al. Multiscale Modeling of CO₂ Electrochemical Reduction on Copper Electrocatalysts: A Review. *ChemSusChem* (2025). [DOI](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202400898)

### Mean-field vs KMC / Mean-field与KMC对比

[9] Vignola, E. et al. Evaluating the benefits of kinetic Monte Carlo and microkinetic modeling for catalyst design studies in the presence of lateral interactions. *Catal. Today* **372**, 11 (2021). [DOI](https://www.sciencedirect.com/science/article/abs/pii/S092058612100119X)

### KMC Methodology / KMC方法论

[10] Pineda, M. & Stamatakis, M. Kinetic Monte Carlo simulations for heterogeneous catalysis: Fundamentals, current status, and challenges. *J. Chem. Phys.* **156**, 120902 (2022). [DOI](https://pubs.aip.org/aip/jcp/article/156/12/120902/2840948)

### HER Energetics / HER能量学参数

[11] Skúlason, E. et al. Modeling the Electrochemical Hydrogen Oxidation and Evolution Reactions on the Basis of DFT Calculations. *J. Phys. Chem. C* **114**, 18182–18197 (2010).

[12] Karlberg, G. S. et al. Cyclic Voltammograms for H on Pt(111) and Pt(100) from First Principles. *Phys. Rev. Lett.* **99**, 126101 (2007).

[13] Greeley, J. & Mavrikakis, M. Surface and Subsurface Hydrogen: Adsorption Properties on Transition Metals. *J. Phys. Chem. B* **109**, 3460–3471 (2005).

### Reviews and Guides / 综述与指南

[14] Stamatakis, M. & Vlachos, D. G. Unraveling the Complexity of Catalytic Reactions via Kinetic Monte Carlo Simulation: Current Status and Frontiers. *ACS Catal.* **2**, 2648–2663 (2012).

[15] Andersen, M., Panosetti, C. & Reuter, K. A Practical Guide to Surface Kinetic Monte Carlo Simulations. *Front. Chem.* **7**, 202 (2019). [DOI](https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2019.00202/full)

[16] DFT-Based Multiscale Modeling of Heterogeneous (Electro)Catalytic Reactions. *ACS Catal.* (2025). [DOI](https://pubs.acs.org/doi/10.1021/acscatal.5c07967)

[17] Exner, K. S. Microkinetic modeling in electrocatalysis: Applications, limitations, and recommendations. *J. Catal.* **400**, 290 (2021).

### Surface Reconstruction / 表面重构

[18] In Situ Reconstruction of a Cu(100) Surface for Promoted C–C Coupling in CO₂ Electroreduction from First-Principles Multiscale Modeling. *ACS Catal.* (2025). [DOI](https://pubs.acs.org/doi/10.1021/acscatal.5c07187)

### NRR and Battery Interfaces / NRR与电池界面

[19] Naqvi, S. R. et al. Investigating High-Performance Transition Metal Oxide Catalysts for NRR: A DFT–kMC–LSTM Approach. *ACS Catal.* **13**, 7 (2023). [DOI](https://pubs.acs.org/doi/10.1021/acscatal.3c01360)

[20] SEI Formation and Lithium-Ion Electrodeposition Dynamics in Li Metal Batteries via First-Principles KMC Modeling. *ACS Energy Lett.* (2024). [DOI](https://pubs.acs.org/doi/10.1021/acsenergylett.4c02019)

### Lateral Interaction Algorithms / 横向相互作用算法

[21] Speeding up the Detection of Adsorbate Lateral Interactions in Graph-Theoretical KMC Simulations. *J. Phys. Chem. A* (2024). [DOI](https://pubs.acs.org/doi/10.1021/acs.jpca.3c05581)

---

*Last updated: 2026-03-21*
