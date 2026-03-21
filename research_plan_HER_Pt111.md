# Research Plan: Validating Leshen-KMC with Electrochemical HER on Pt(111)

# 研究方案：基于 Pt(111) 析氢反应验证 Leshen-KMC 电催化模拟能力

---

## 1. Motivation / 研究动机

### 1.1 Background / 背景

Kinetic Monte Carlo (KMC) simulation has become an essential tool for bridging the gap between DFT-calculated energetics and experimentally observed catalytic activity. However, **all existing general-purpose KMC software packages (Zacros, kmos, SPPARKS, KMCLib, MonteCoffee) are designed exclusively for thermal catalysis** — none natively supports electrochemical potential-dependent rates, Butler-Volmer kinetics, or polarization curve computation.

动力学蒙特卡洛（KMC）模拟已成为连接 DFT 计算与实验催化活性的关键工具。然而，**现有所有通用 KMC 软件（Zacros、kmos、SPPARKS、KMCLib、MonteCoffee）均仅支持热催化**，不原生支持电化学电位依赖速率、Butler-Volmer 动力学或极化曲线计算。

Recent electrocatalytic KMC studies (Wei et al., JPCL 2025; Li et al., ACS Catal. 2024; Shou et al., Nat. Comm. 2025) have all used **one-off custom codes**, highlighting the need for a general-purpose electrocatalytic KMC framework.

近期的电催化 KMC 研究均使用**一次性自编代码**，凸显了通用电催化 KMC 框架的迫切需求。

### 1.2 Objective / 目标

**Validate Leshen-KMC** — the first general-purpose KMC framework for electrocatalysis — using the hydrogen evolution reaction (HER) on Pt(111) as a benchmark. This system is ideal because:

**验证 Leshen-KMC** —— 首个面向电催化的通用 KMC 框架 —— 以 Pt(111) 上析氢反应（HER）为基准。选择此体系因为：

1. **Simplest electrochemical reaction** — only 2 species (H*, empty), 3 reaction steps
2. **Exercises all software features** — PCET rates, lateral interactions, surface diffusion, polarization curves, KMC vs MKM comparison
3. **Most well-characterized system** — extensive experimental polarization data available
4. **Coverage-dependent kinetics** — H*-H* repulsion at high coverage causes mean-field breakdown, where spatial KMC is essential
5. **DFT data fully available** — Li et al. (ACS Catal. 2024) provides complete energetics

---

## 2. Model Description / 模型描述

### 2.1 Reaction Network / 反应网络

The HER on Pt(111) in acidic media proceeds via three elementary steps:

酸性介质中 Pt(111) 上 HER 包含三个基元步骤：

| Step | Reaction | Type | Description |
|------|----------|------|-------------|
| **Volmer** | H⁺(aq) + e⁻ + * → H* | PCET | Electrochemical adsorption / 电化学吸附 |
| **Tafel** | H* + H* → H₂(g) + 2* | Thermal | Recombinative desorption (requires adjacent sites) / 复合脱附（需相邻位点） |
| **Heyrovsky** | H* + H⁺(aq) + e⁻ → H₂(g) + * | PCET | Electrochemical desorption / 电化学脱附 |

Additional process:

| Step | Reaction | Type |
|------|----------|------|
| **H* diffusion** | H*@site₁ + *@site₂ → *@site₁ + H*@site₂ | Thermal |

### 2.2 Rate Expressions / 速率表达式

**PCET steps (Volmer, Heyrovsky) — Butler-Volmer kinetics:**

$$k(U) = \frac{k_B T}{h} \exp\left(-\frac{E_a + \beta_{BV} \cdot e \cdot (U - U_0)}{k_B T}\right)$$

where β_BV is the symmetry factor (typically 0.5), U is the applied potential vs RHE, and U₀ is the equilibrium potential.

**Thermal steps (Tafel, diffusion) — Transition State Theory:**

$$k = \frac{k_B T}{h} \exp\left(-\frac{E_a}{k_B T}\right)$$

### 2.3 DFT Parameters / DFT 参数

From Li et al. (ACS Catal. 2024) and related literature:

| Parameter | Value | Source |
|-----------|-------|--------|
| E_a(Volmer, fwd) at θ→0 | 0.67 eV | Li et al. 2024 |
| E_a(Volmer, rev) at θ→0 | 0.62 eV | Li et al. 2024 |
| E_a(Tafel) at θ→0 | 0.85 eV | Skúlason et al. 2010 |
| E_a(Heyrovsky, fwd) | 0.70 eV | Li et al. 2024 |
| E_a(diffusion) | 0.10 eV | Greeley & Mavrikakis 2004 |
| ε(H*-H*) lateral | +0.10 eV (repulsive) | Karlberg et al. 2007 |
| β_BV (symmetry factor) | 0.5 | Standard |
| T | 298 K | Room temperature |

### 2.4 Lattice / 晶格

- 2D square lattice (approximating Pt(111) top sites)
- Periodic boundary conditions
- Size: 50×50 (2500 sites) for production runs
- 4 nearest neighbors per site

### 2.5 Lateral Interactions / 横向相互作用

Pairwise nearest-neighbor H*-H* repulsion:
- ε(H*, H*) = +0.10 eV (repulsive)
- At θ_H = 1 ML, each H* has 4 repulsive neighbors → effective binding weakened by 0.40 eV

BEP relation for Volmer step:
- α = 0.5 (linear scaling between ΔE_ads and ΔE_a)

---

## 3. Computational Methods / 计算方法

### 3.1 Leshen-KMC Features Used / 使用的 Leshen-KMC 功能

| Feature | Role in HER Model |
|---------|-------------------|
| **Butler-Volmer PCET** | Potential-dependent Volmer/Heyrovsky rates |
| **Pairwise lateral interactions** | H*-H* repulsion at NN sites |
| **BEP relations** | Coverage-dependent Volmer barrier |
| **Surface diffusion** | H* hop between NN sites |
| **Neighbor list** | Tafel step requires adjacent H* pair |
| **KMC engine** | Stochastic dynamics with spatial resolution |
| **MKM solver** | Mean-field comparison (Newton-Raphson) |
| **Polarization curve** | j-V output for experimental comparison |

### 3.2 Simulation Protocol / 模拟流程

```
For each potential U in [-0.5, -0.4, ..., 0.0] V vs RHE:
    1. Initialize 50×50 lattice, all sites empty
    2. Set T = 298 K, U = U_i
    3. Equilibrate: 10⁶ KMC steps (discard)
    4. Production: 10⁶ KMC steps, record every 10⁴ steps
    5. Compute average θ_H, TOF(H₂)
    6. Convert: j = TOF × n_e × e / A_site  (n_e = 2 per H₂)
```

### 3.3 KMC vs MKM Comparison / KMC 与 MKM 对比

Run the **identical reaction model** with:
- **KMC**: Full spatial simulation on 50×50 lattice
- **MKM**: Mean-field ODE steady-state (Newton-Raphson)
- Compare: θ_H(U), TOF(U), j(U) from both methods
- Quantify mean-field error at high coverage

---

## 4. Validation Tests / 验证测试

### Test 1: Langmuir Limit / Langmuir 极限

**Turn off lateral interactions and Tafel step.** Only Volmer adsorption/desorption remains.

Expected: θ_H should follow the Langmuir isotherm:

$$\theta_H = \frac{K \cdot a_{H^+}}{1 + K \cdot a_{H^+}}$$

This validates the basic KMC engine and Butler-Volmer rate implementation.

### Test 2: Tafel Slope / Tafel 斜率

At low overpotential (η < 50 mV), the polarization curve should follow:

$$\eta = -\frac{2.303 \cdot k_B T}{\alpha \cdot e} \cdot \log_{10}\left(\frac{j}{j_0}\right)$$

- Volmer-limited: ~120 mV/decade (α = 0.5)
- Tafel-limited: ~30 mV/decade
- Heyrovsky-limited: ~40 mV/decade

The simulated Tafel slope reveals the rate-determining step.

### Test 3: Coverage Effect / 覆盖度效应

Compare θ_H with and without lateral interactions at fixed U = -0.2 V:

- Without lateral: θ_H ≈ 0.95 (Langmuir at strong binding)
- With lateral (ε = +0.10 eV): θ_H should decrease significantly

This validates the lateral interaction implementation.

### Test 4: Diffusion Effect / 扩散效应

Compare Tafel pathway activity with and without H* diffusion:

- Without diffusion: Tafel rate limited by random adjacency of H*
- With diffusion: H* can migrate to find partners, Tafel rate increases

This validates the diffusion mechanism and its coupling with spatial events.

### Test 5: KMC vs MKM Divergence / KMC 与 MKM 偏差

At high θ_H (U < -0.3 V), mean-field assumes uniform coverage but KMC may show:
- Non-uniform H* distribution
- Island/domain formation at certain conditions
- Quantitative TOF differences

This demonstrates the unique value of spatial KMC over mean-field models.

### Test 6: Lattice Size Convergence / 格点尺寸收敛

Run at U = -0.2 V with lattice sizes: 10×10, 20×20, 50×50, 100×100.

Expected: Results converge by 50×50 (2500 sites). This establishes the minimum lattice size for production simulations.

---

## 5. Expected Results / 预期结果

### 5.1 Polarization Curve / 极化曲线

```
          │
    log|j|│          ○ ○ ○ ○    ← KMC
          │        ○         ○
          │      ○    Tafel slope
          │    ○     ~30-40 mV/dec
          │  ○
          │○
          └──────────────────────
          0    -0.1   -0.2  -0.3  U [V vs RHE]
```

Expected exchange current density: j₀ ≈ 1-10 mA/cm² (consistent with experimental Pt HER).

### 5.2 KMC vs MKM Comparison / KMC 与 MKM 对比

| U [V vs RHE] | θ_H (MKM) | θ_H (KMC) | Relative error |
|--------------|-----------|-----------|----------------|
| 0.0 | ~0.7 | ~0.65 | ~7% |
| -0.1 | ~0.85 | ~0.80 | ~6% |
| -0.2 | ~0.95 | ~0.88 | ~7% |
| -0.3 | ~0.99 | ~0.92 | ~7% |

Mean-field overestimates θ_H because it ignores spatial correlations from H*-H* repulsion.

### 5.3 Demonstration of Software Advantages / 软件优势展示

| Advantage | Demonstrated by |
|-----------|-----------------|
| **First electrocatalytic KMC framework** | HER model with Butler-Volmer PCET — impossible in Zacros/kmos |
| **KMC + MKM same framework** | Same model definition, two solvers, direct comparison |
| **Polarization curve output** | j-V curve from KMC — unique capability |
| **Spatial effects matter** | KMC vs MKM divergence at high coverage |
| **Coverage-dependent kinetics** | Lateral interaction + BEP correctly modifies rates |
| **High performance** | Fenwick tree enables 50×50 lattice with lateral interactions |

---

## 6. Timeline / 时间计划

| Phase | Task | Duration |
|-------|------|----------|
| **Phase 1** | Collect DFT parameters from literature | 1 day |
| **Phase 2** | Implement HER model (Python + Rust) | 1-2 days |
| **Phase 3** | Run validation tests (Tests 1-6) | 1-2 days |
| **Phase 4** | Generate polarization curves + benchmarks | 1 day |
| **Phase 5** | Analysis and documentation | 1 day |
| **Total** | | **5-7 days** |

---

## 7. Future Extension / 后续扩展

After HER validation:

1. **CO₂RR on Cu(100)** — Multi-product selectivity (CH₄, C₂H₄, CO), CO* dimerization, ~15 processes
2. **NRR on Mo** — Extend existing built-in model with lateral interactions
3. **OER on IrO₂** — Oxygen evolution with multiple oxidation states
4. **Multi-facet catalysts** — Use site types for step/terrace differentiation

---

## 8. References / 参考文献

1. Li, H. et al. First-Principles-Based Kinetic Monte Carlo Model of Hydrogen Evolution Reaction under Realistic Conditions. *ACS Catal.* **14**, 2696–2708 (2024). DOI: 10.1021/acscatal.3c04588

2. Wei, C. et al. Voltage-Dependent Electrochemical Carbon Dioxide Reduction Mechanism Unveiled by Kinetic Monte Carlo Simulation. *J. Phys. Chem. Lett.* **16**, 2896–2904 (2025). DOI: 10.1021/acs.jpclett.4c03426

3. Shou, W. et al. Toward rational understanding of the hydrogen evolution polarization curves through multiscale simulations. *Nat. Commun.* (2025). DOI: 10.1038/s41467-025-67493-y

4. Skúlason, E. et al. Modeling the Electrochemical Hydrogen Oxidation and Evolution Reactions on the Basis of Density Functional Theory Calculations. *J. Phys. Chem. C* **114**, 18182–18197 (2010).

5. Karlberg, G. S. et al. Cyclic Voltammograms for H on Pt(111) and Pt(100) from First Principles. *Phys. Rev. Lett.* **99**, 126101 (2007).

6. Greeley, J. & Mavrikakis, M. Surface and Subsurface Hydrogen: Adsorption Properties on Transition Metals and Near-Surface Alloys. *J. Phys. Chem. B* **109**, 3460–3471 (2005).

7. Pineda, M. & Stamatakis, M. Kinetic Monte Carlo simulations for heterogeneous catalysis: Fundamentals, current status, and challenges. *J. Chem. Phys.* **156**, 120902 (2022).

8. Deshpande, S. et al. Evaluating the benefits of kinetic Monte Carlo and microkinetic modeling for catalyst design studies in the presence of lateral interactions. *Catal. Sci. Technol.* **11**, 4946 (2021).
