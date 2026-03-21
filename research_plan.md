# Research Plan: Spatial Kinetic Monte Carlo Simulation of Electrocatalytic Reactions

# 研究方案：电催化反应的空间动力学蒙特卡洛模拟

**Leshen-KMC: First General-Purpose Electrocatalytic KMC Framework**

---

## Abstract / 摘要

We propose a two-phase computational study using Leshen-KMC — the first general-purpose KMC framework with native electrochemical support — to investigate electrocatalytic reactions where spatial effects critically determine activity and selectivity. **Phase 1** validates the software using the hydrogen evolution reaction (HER) on Pt(111), reproducing experimental polarization curves. **Phase 2** addresses a central open question in CO₂ electroreduction: *how does the spatial distribution of CO* adsorbates on Cu(100) control C₂ product selectivity?* This question cannot be answered by mean-field models and has never been studied with a general-purpose electrocatalytic KMC tool.

本方案提出一项两阶段计算研究，使用 Leshen-KMC —— 首个原生支持电化学的通用 KMC 框架。**第一阶段**以 Pt(111) 上析氢反应（HER）验证软件，复现实验极化曲线。**第二阶段**研究 CO₂ 电还原中的核心开放问题：*Cu(100) 上 CO* 吸附物的空间分布如何决定 C₂ 产物选择性？* 此问题 mean-field 模型无法回答，且从未有通用电催化 KMC 工具对其进行研究。

---

## 1. Introduction / 引言

### 1.1 The Gap: No General-Purpose Electrocatalytic KMC Tool

Kinetic Monte Carlo simulation is the gold standard for bridging DFT-computed energetics and macroscopic catalytic behavior, capturing spatial correlations, coverage fluctuations, and island formation that mean-field models miss. However, **all existing general-purpose KMC packages were designed for thermal catalysis:**

KMC 模拟是连接 DFT 计算能量学与宏观催化行为的黄金标准方法，能够捕捉 mean-field 模型遗漏的空间关联、覆盖度涨落和岛状聚集。然而，**现有所有通用 KMC 软件均为热催化设计：**

| Software | Language | Electrochemistry | Polarization Curve |
|----------|----------|:-----------------:|:------------------:|
| Zacros | Fortran | ✗ | ✗ |
| kmos | Python/Fortran | ✗ | ✗ |
| SPPARKS | C++ | ✗ | ✗ |
| KMCLib | Python/C++ | ✗ | ✗ |
| MonteCoffee | Python | ✗ | ✗ |
| **Leshen-KMC** | **Python + Rust** | **✓ Butler-Volmer PCET** | **✓ Native** |

Recent electrocatalytic KMC studies have all used one-off custom codes:

| Study | System | Code |
|-------|--------|------|
| Wei et al., JPCL 2025 | CO₂RR on Cu, 178 reactions, LSV prediction | Custom |
| Li et al., ACS Catal. 2024 | HER on Pt, coverage-dependent, polarization | Custom |
| Shou et al., Nat. Comm. 2025 | HER on Au, multiscale, full pH range | Custom |
| Deshpande et al., JCP 2021 | CO₂RR on Cu, surface diffusion effects | Custom |

**Leshen-KMC fills this gap** as the first general-purpose framework with native Butler-Volmer rates, polarization curve computation, and integrated KMC + mean-field microkinetic modeling.

### 1.2 Why Spatial KMC Matters for Electrocatalysis

Mean-field microkinetic models assume uniform surface coverage:

$$r_{C-C} = k_{CC} \cdot \theta_{CO}^2$$

But in reality, CO* molecules form spatial patterns (islands, stripes, random) depending on:
- **Lateral interactions** — attraction promotes clustering, repulsion promotes dispersion
- **Surface diffusion rate** — fast diffusion → well-mixed (mean-field valid); slow diffusion → spatial correlations
- **Relative rates of production vs consumption** — imbalance creates local concentration gradients

Deshpande et al. (Catal. Sci. Tech. 2021) showed that **mean-field models can overestimate catalytic activity by several orders of magnitude** when lateral interactions are significant. For C-C coupling (a bimolecular surface reaction requiring adjacent sites), this error is expected to be particularly severe.

### 1.3 The Central Question

> **Does the spatial arrangement of CO* adsorbates on Cu(100) — controlled by lateral interactions and diffusion — determine the selectivity between C₁ (CH₄) and C₂ (C₂H₄) products in electrochemical CO₂ reduction?**

This question is:
- **Experimentally observed but mechanistically unexplained** — higher CO* coverage correlates with higher C₂ selectivity (Nature Energy 2024)
- **Impossible to answer with mean-field models** — they assume uniform θ_CO
- **Uniquely suited for spatial KMC** — the only computational method that tracks individual adsorbate positions and their evolution

---

## 2. Phase 1: Software Validation — HER on Pt(111)

### 2.1 Objective / 目标

Validate Leshen-KMC against analytical solutions and published computational/experimental data using the simplest electrochemical reaction: hydrogen evolution on Pt(111).

### 2.2 Model

**Species:** empty (*), H*

**Reactions:**

| # | Step | Reaction | Type | E_a [eV] | Source |
|---|------|----------|------|----------|--------|
| 1 | Volmer (fwd) | H⁺ + e⁻ + * → H* | PCET | 0.67 | Li 2024 |
| 2 | Volmer (rev) | H* → H⁺ + e⁻ + * | PCET | 0.62 | Li 2024 |
| 3 | Tafel | H* + H*(adj) → H₂ + 2* | Thermal | 0.85 | Skúlason 2010 |
| 4 | Heyrovsky (fwd) | H* + H⁺ + e⁻ → H₂ + * | PCET | 0.70 | Li 2024 |
| 5 | H* diffusion | H*@s₁ + *@s₂ → *@s₁ + H*@s₂ | Thermal | 0.10 | Greeley 2004 |

**Lateral interaction:** ε(H*, H*) = +0.10 eV (repulsive, NN), from Karlberg et al. 2007

**BEP:** α = 0.5 for Volmer step

**Lattice:** 50×50 2D square, PBC, T = 298 K

### 2.3 Validation Tests

| Test | Method | Expected Outcome | Validates |
|------|--------|------------------|-----------|
| 1. Langmuir limit | Turn off lateral + Tafel | θ_H = Langmuir isotherm | Basic engine + PCET rates |
| 2. Tafel slope | j-V at low η | 30-120 mV/dec depending on RDS | Butler-Volmer implementation |
| 3. Lateral effect | θ_H with/without ε | θ_H decreases with repulsion | Lateral interaction engine |
| 4. Diffusion effect | TOF_Tafel with/without diff | Diffusion enhances Tafel rate | Diffusion + spatial events |
| 5. KMC vs MKM | Same model, both methods | MKM overestimates θ_H at high coverage | Spatial correlation effects |
| 6. Lattice convergence | 10² to 10⁴ sites | Converge by 50×50 | Finite-size effects |

### 2.4 Deliverables

- Validated polarization curve j(U) comparable to experimental Pt HER
- Quantified KMC vs MKM divergence at high H* coverage
- Demonstrated Tafel slope analysis from KMC data

---

## 3. Phase 2: CO₂RR on Cu(100) — Spatial Control of C₂ Selectivity

### 3.1 Scientific Hypothesis / 科学假说

**Hypothesis:** The C₂H₄/CH₄ selectivity ratio in electrochemical CO₂ reduction on Cu(100) is primarily controlled by the spatial distribution of CO* adsorbates, which is tuned by:
1. CO*-CO* lateral interaction strength (attractive → island formation → more C-C coupling)
2. CO* surface diffusion rate (fast → well-mixed → mean-field applies; slow → spatial patterns)
3. Applied potential (controls CO* production rate and consumption pathways)

**假说：** Cu(100) 上电化学 CO₂ 还原的 C₂H₄/CH₄ 选择性主要由 CO* 吸附物的空间分布控制，后者由横向相互作用强度、表面扩散速率和施加电位调节。

### 3.2 Reaction Model / 反应模型

**Species (8):** empty, CO₂*, COOH*, CO*, CHO*, OCCO*, CH₂O*, CH₃O*

**Lattice:** 2D square (Cu(100) top sites), 50×50, PBC

**Elementary Steps (~12):**

| # | Reaction | Type | E_a [eV] | Note |
|---|----------|------|----------|------|
| 1 | CO₂(g) + * → CO₂* | Adsorption | ~0 | Hertz-Knudsen |
| 2 | CO₂* + H⁺ + e⁻ → COOH* | PCET | 0.43 | Montoya 2015 |
| 3 | COOH* + H⁺ + e⁻ → CO* + H₂O | PCET | 0.02 | Nearly barrierless |
| 4 | CO* → CO(g) + * | Desorption | 0.55 | Controls CO* coverage |
| 5 | **CO* + CO*(adj) → OCCO*** | **Thermal** | **0.72** | **C-C coupling (C₂ key step)** |
| 6 | CO* + H⁺ + e⁻ → CHO* | PCET | 0.74 | C₁ pathway entry |
| 7 | OCCO* + H⁺ + e⁻ → ... → C₂H₄ | PCET | 0.40 | C₂ downstream (lumped) |
| 8 | CHO* + H⁺ + e⁻ → CH₂O* | PCET | 0.35 | C₁ downstream |
| 9 | CH₂O* + H⁺ + e⁻ → CH₃O* | PCET | 0.30 | C₁ downstream |
| 10 | CH₃O* + H⁺ + e⁻ → CH₄ + * | PCET | 0.25 | C₁ product release |
| 11 | CO* diffusion | Thermal | 0.10 | Tunable parameter |

**Lateral Interactions:**

| Pair | Energy | Physical Origin |
|------|--------|----------------|
| CO*-CO* | Variable: −0.05 to +0.10 eV | DFT, scan as parameter |
| CO*-CHO* | +0.05 eV (repulsive) | Steric |

**BEP:** α = 0.5 for Volmer step on PCET reactions

### 3.3 Key Simulations / 关键模拟

#### Study A: CO* Spatial Distribution Mapping

**Method:** Fix U = −1.0 V, T = 298 K. Run KMC on 50×50 lattice. Vary CO*-CO* lateral interaction:

| ε(CO*-CO*) | Expected CO* Pattern | C-C Coupling Rate |
|------------|---------------------|-------------------|
| −0.05 eV (attractive) | Island/cluster formation | **Enhanced** (many adjacent pairs) |
| 0.00 eV (no interaction) | Random distribution | Baseline |
| +0.05 eV (repulsive) | Dispersed/ordered | **Reduced** (few adjacent pairs) |
| +0.10 eV (strong repulsive) | Anti-clustering | **Strongly reduced** |

**Output:** CO* pair correlation function g(r), cluster size distribution, snapshot images

#### Study B: C₁/C₂ Selectivity vs Potential

**Method:** Scan U from −0.6 to −1.4 V. At each potential, run KMC + MKM.

**Output:**
- Faradaic efficiency: FE(CH₄) and FE(C₂H₄) vs U
- KMC vs MKM comparison of selectivity
- Identify crossover potential where C₂ > C₁

**Key prediction:** KMC should predict higher C₂ selectivity than MKM when CO*-CO* interaction is attractive (islands enhance C-C coupling). MKM should overestimate C₁ pathway because it assumes uniform CO* distribution.

#### Study C: Diffusion Rate Effect

**Method:** Fix U = −1.0 V, ε = −0.05 eV. Vary E_a(diffusion):

| E_a(diff) [eV] | Regime | Expected Effect |
|-----------------|--------|-----------------|
| 0.05 | Very fast (mean-field limit) | KMC ≈ MKM |
| 0.10 | Moderate | KMC ≠ MKM, spatial effects emerge |
| 0.30 | Slow | Strong spatial patterns, KMC ≫ MKM deviation |
| 0.50 | Very slow | Frozen CO*, maximal spatial effects |

**Output:** C₂/C₁ ratio vs diffusion barrier, quantify mean-field breakdown point

#### Study D: Polarization Curve and Faradaic Efficiency

**Method:** Full potential scan, compute j_total, j_CH₄, j_C₂H₄

**Output:**
```
U [V vs RHE]  |  j_total  |  FE(CH₄)  |  FE(C₂H₄)  |  FE(CO)
─────────────────────────────────────────────────────────────
   −0.6       |   0.1     |   10%     |    2%       |   80%
   −0.8       |   1.0     |   25%     |   15%       |   50%
   −1.0       |   5.0     |   30%     |   35%       |   20%
   −1.2       |  15.0     |   25%     |   50%       |   10%
```

Compare with experimental data from Hori et al. and recent operando studies.

### 3.4 Expected Scientific Contributions / 预期科学贡献

| # | Contribution | Novelty |
|---|-------------|---------|
| 1 | **First spatial KMC study of CO* distribution → C₂ selectivity** | KMC essential, MKM cannot answer |
| 2 | **Quantitative mean-field error for C-C coupling predictions** | Literature gap: Deshpande showed error exists but not for CO₂RR C₂ |
| 3 | **CO*-CO* lateral interaction as design knob for C₂ selectivity** | Connects DFT interaction energy to macroscopic selectivity |
| 4 | **Diffusion rate threshold for mean-field validity** | Practical guide: when is MKM sufficient? |
| 5 | **First general-purpose electrocatalytic KMC framework application** | Software contribution to community |
| 6 | **Electrochemical polarization curve from spatial KMC** | Unique capability, no other KMC software can do this |

---

## 4. Computational Details / 计算细节

### 4.1 Software

- **Leshen-KMC v0.3.0** (Python + Rust)
- Rust engine for production runs (Fenwick tree O(log N) site selection)
- Python engine for prototyping and analysis
- Newton-Raphson MKM solver for mean-field comparison

### 4.2 Simulation Parameters

| Parameter | Phase 1 (HER) | Phase 2 (CO₂RR) |
|-----------|---------------|-----------------|
| Lattice | 50×50 | 50×50 |
| Temperature | 298 K | 298 K |
| Potential range | −0.5 to 0.0 V | −0.6 to −1.4 V |
| Equilibration | 10⁶ steps | 10⁶ steps |
| Production | 10⁶ steps | 5×10⁶ steps |
| Sampling | Every 10⁴ steps | Every 10⁴ steps |
| Independent runs | 5 (for error bars) | 5 |

### 4.3 Analysis Tools

- Coverage time series: θ_i(t)
- TOF per product: TOF_CH₄, TOF_C₂H₄
- Pair correlation function: g_CO-CO(r)
- Cluster size distribution: P(n)
- Lattice snapshots: species coloring
- Polarization curves: j(U) and FE(U)

---

## 5. Timeline / 时间计划

| Week | Phase | Tasks |
|------|-------|-------|
| **1** | Phase 1 | T1-T2: Collect DFT data, implement HER model (Python + Rust) |
| **2** | Phase 1 | T3-T4: Run validation tests, generate polarization curves |
| **3** | Phase 2 | T5-T6: Collect CO₂RR DFT data, implement Cu(100) model |
| **4** | Phase 2 | T7: CO* spatial distribution analysis (Study A) |
| **5** | Phase 2 | T8: C₁/C₂ selectivity vs potential (Study B), KMC vs MKM |
| **6** | Phase 2 | T9: Diffusion rate study (Study C) |
| **7** | Phase 2 | T10: Polarization curves and Faradaic efficiency (Study D) |
| **8** | Writing | Manuscript preparation |

---

## 6. Significance / 研究意义

### 6.1 Scientific Impact

This work will provide the **first quantitative answer** to whether CO* spatial arrangement controls C₂ selectivity in CO₂ electroreduction — a question that has been experimentally motivated (Nature Energy 2024) but computationally inaccessible due to the lack of electrocatalytic spatial KMC tools.

### 6.2 Methodological Impact

The study demonstrates that **spatial KMC is essential for electrochemical systems with bimolecular surface reactions** (C-C coupling), establishing when mean-field models fail and KMC is required. This provides practical guidelines for the computational catalysis community.

### 6.3 Software Impact

Leshen-KMC, as the **first open-source general-purpose electrocatalytic KMC framework**, enables the broader community to perform spatial KMC simulations for any electrocatalytic system without writing custom code — lowering the barrier from months of development to hours of model definition.

---

## 7. References / 参考文献

### Electrocatalytic KMC Studies

1. Wei, C. et al. Voltage-Dependent Electrochemical Carbon Dioxide Reduction Mechanism Unveiled by Kinetic Monte Carlo Simulation. *J. Phys. Chem. Lett.* **16**, 2896–2904 (2025).

2. Li, H. et al. First-Principles-Based Kinetic Monte Carlo Model of Hydrogen Evolution Reaction under Realistic Conditions. *ACS Catal.* **14**, 2696–2708 (2024).

3. Shou, W. et al. Toward rational understanding of the hydrogen evolution polarization curves through multiscale simulations. *Nat. Commun.* (2025).

4. Deshpande, S. et al. Effects of surface diffusion in electrocatalytic CO₂ reduction on Cu revealed by kinetic Monte Carlo simulations. *J. Chem. Phys.* **155**, 164701 (2021).

### CO₂RR Mechanism and Selectivity

5. Zhan, C. et al. Key intermediates and Cu active sites for CO₂ electroreduction to ethylene and ethanol. *Nat. Energy* **9**, 1485–1496 (2024).

6. Montoya, J. H. et al. Theoretical Insights into a CO Dimerization Mechanism in CO₂ Electroreduction. *J. Phys. Chem. Lett.* **6**, 2032–2037 (2015).

7. Nitopi, S. et al. Progress and Perspectives of Electrochemical CO₂ Reduction on Copper in Aqueous Electrolyte. *Chem. Rev.* **119**, 7610–7672 (2019).

8. Gholizadeh, M. et al. Multiscale Modeling of CO₂ Electrochemical Reduction on Copper Electrocatalysts. *ChemSusChem* (2025).

### KMC Methodology

9. Deshpande, S. et al. Evaluating the benefits of kinetic Monte Carlo and microkinetic modeling for catalyst design studies in the presence of lateral interactions. *J. Catal.* **401**, 113–119 (2021).

10. Pineda, M. & Stamatakis, M. Kinetic Monte Carlo simulations for heterogeneous catalysis: Fundamentals, current status, and challenges. *J. Chem. Phys.* **156**, 120902 (2022).

### HER Energetics

11. Skúlason, E. et al. Modeling the Electrochemical Hydrogen Oxidation and Evolution Reactions on the Basis of Density Functional Theory Calculations. *J. Phys. Chem. C* **114**, 18182–18197 (2010).

12. Karlberg, G. S. et al. Cyclic Voltammograms for H on Pt(111) and Pt(100) from First Principles. *Phys. Rev. Lett.* **99**, 126101 (2007).

13. Greeley, J. & Mavrikakis, M. Surface and Subsurface Hydrogen: Adsorption Properties on Transition Metals. *J. Phys. Chem. B* **109**, 3460–3471 (2005).
