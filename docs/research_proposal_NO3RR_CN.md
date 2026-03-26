# 基于第一性原理的空间KMC模拟揭示Cu(100)上NO₃RR的NH₃/N₂选择性调控机理

# First-Principles Spatial KMC Simulation Reveals NH₃/N₂ Selectivity Control in Electrochemical NO₃RR on Cu(100)

## 课题思路 / Research Proposal

> SPARK Phase 2 | 2026-03-24
> Author: Ren (UCSD)

---

## 一、研究背景与意义 / 1. Background and Significance

### 1.1 电化学硝酸盐还原（NO₃RR）的重要性 / Importance of Electrochemical NO₃RR

电化学硝酸盐还原反应（NO₃RR）是同时解决环境污染和绿色化学品合成两大挑战的关键技术。地下水和工业废水中的硝酸盐（NO₃⁻）是全球性的水污染问题，而电化学方法可将其选择性还原为高附加值的氨（NH₃），实现"变废为宝"。

Electrochemical nitrate reduction (NO₃RR) addresses both environmental pollution and green chemical synthesis. Nitrate (NO₃⁻) is a major water contaminant, and electrochemistry can selectively convert it into valuable ammonia (NH₃).

Cu基催化剂是NO₃RR最有前景的非贵金属催化剂，但面临核心挑战：**产物选择性的精确调控**。NO₃RR可生成NH₃（8e⁻，目标产物）、N₂（10e⁻）、N₂O（8e⁻）和NO₂⁻（2e⁻），其中NH₃/N₂选择性问题尤为突出——两条路径在关键中间体\*NO处分叉，且活化能几乎相同（Ea ≈ 0.05–0.10 eV），使传统动力学分析失效。

Cu-based catalysts are the most promising non-noble-metal catalysts for NO₃RR, but the core challenge is **precise product selectivity control**. NO₃RR produces NH₃ (8e⁻, target), N₂ (10e⁻), N₂O (8e⁻), and NO₂⁻ (2e⁻). The NH₃/N₂ selectivity problem is particularly critical — the two pathways branch at the key intermediate \*NO with nearly identical activation barriers (Ea ≈ 0.05–0.10 eV), rendering traditional kinetic analysis ineffective.

### 1.2 Mean-field微动力学的局限性 / Limitations of Mean-Field Microkinetics

现有NO₃RR理论研究主要依赖mean-field微动力学模型（MKM），其核心假设是吸附物在催化剂表面**均匀分布**：

Existing NO₃RR theoretical studies rely on mean-field microkinetic models (MKM), which assume adsorbates are **uniformly distributed** on the catalyst surface:

$$r_{bimol} = k \cdot \theta_A \cdot \theta_B$$

然而，多个体系的研究已表明这一假设在强横向相互作用和扩散受限条件下严重失效：

However, studies on multiple systems have shown this assumption breaks down severely under strong lateral interactions and diffusion-limited conditions:

| 体系 / System | MKM预测 / MKM Prediction | KMC/实验 / KMC/Experiment | 偏差原因 / Cause |
|------|--------|---------|---------|
| CO氧化/RuO₂ | 活性高 / High activity | 低数个数量级 / Orders of magnitude lower | 空间关联 / Spatial correlations |
| NH₃氧化/RuO₂ | NO选择性低 / Low NO selectivity | 93% NO | 空间效应主导 / Spatial effects dominate |
| HER/Pt(111) | θ_H = 0.76 | θ_H = 0.33 | H\*-H\*排斥→反聚集 / H\*-H\* repulsion → anti-clustering |

### 1.3 空间KMC的不可替代性 / Irreplaceability of Spatial KMC

Kinetic Monte Carlo（KMC）是唯一能同时追踪每个表面位点占据状态和逐事件时间演化的模拟方法：

KMC is the only simulation method that simultaneously tracks the occupancy of every surface site and time evolution event-by-event:

- **空间分辨 / Spatial resolution**：追踪每个吸附物位置，捕捉岛状聚集、有序相 / Tracks each adsorbate position, captures island formation and ordered phases
- **真实双分子动力学 / True bimolecular dynamics**：N-N偶联速率由实际相邻\*NO对数决定 / N-N coupling rate determined by actual number of adjacent \*NO pairs
- **横向相互作用的非平均场效应 / Non-mean-field lateral interaction effects**：排斥/吸引导致空间重分布 / Repulsion/attraction causes spatial redistribution
- **显式扩散 / Explicit diffusion**：有限扩散打破局部平衡 / Finite diffusion breaks local equilibrium

### 1.4 文献空白 / Literature Gap

截至2026年3月 / As of March 2026:
- 电催化KMC已有少量研究（CO₂RR [1]、HER [2,3]），但**NO₃RR的KMC研究完全空白** / A few electrocatalytic KMC studies exist (CO₂RR [1], HER [2,3]), but **NO₃RR KMC is entirely absent**
- NO₃RR的MKM有若干报道 [20,24]，但均无法回答空间选择性问题 / NO₃RR MKM has several reports [20,24], but none can address spatial selectivity
- 不存在通用的电催化KMC软件框架 / No general-purpose electrocatalytic KMC software framework exists

---

## 二、核心科学问题 / 2. Central Scientific Question

> **Cu(100)上\*NO吸附物的空间分布——由横向相互作用（ε）、扩散势垒和\*H竞争共同控制——是否是决定电化学NO₃RR中NH₃与N₂选择性的关键因素？**
>
> **Does the spatial distribution of \*NO adsorbates on Cu(100) — controlled by lateral interactions (ε), diffusion barriers, and \*H competition — determine the NH₃ vs. N₂ selectivity in electrochemical NO₃RR?**

子问题 / Sub-questions:

1. **ε(\*NO–\*NO)如何决定\*NO空间分布？** / How does ε(\*NO–\*NO) determine \*NO spatial distribution?
   - ε < 0 (吸引/attractive) → 岛状聚集/island clustering → N₂↑?
   - ε > 0 (排斥/repulsive) → 分散/dispersed → NH₃↑?

2. **MKM在什么条件下失效？** / Under what conditions does MKM fail?
   - 低过电位 vs 高过电位? / Low vs. high overpotential?

3. **\*NO扩散速率是否是KMC-MKM偏差的控制参数？** / Is \*NO diffusion rate the control parameter for KMC-MKM deviation?
   - 临界Damköhler数 Da = k_reaction / k_diffusion? / Critical Damköhler number?

4. **KMC能否定量预测实验可观测量？** / Can KMC quantitatively predict experimental observables?
   - 极化曲线 j(U), 法拉第效率 FE(U) / Polarization curves, Faradaic efficiency

---

## 三、研究体系 / 3. System Selection

### 3.1 体系选择：NO₃RR on Cu(100) / System: NO₃RR on Cu(100)

经过系统对比评估（40+文献调研），最终确认NO₃RR on Cu(100)：

After systematic evaluation (40+ literature review), NO₃RR on Cu(100) was selected:

| 评估维度 / Criterion | NO₃RR优势 / NO₃RR Advantage |
|---------|----------|
| **KMC核心价值 / KMC core value** | \*NO+\*NO偶联是双分子表面反应 / \*NO+\*NO coupling is bimolecular surface reaction |
| **选择性问题 / Selectivity** | NH₃ vs N₂势垒几乎相同，由空间分布决定 / Nearly equal barriers; spatial distribution determines selectivity |
| **复杂度 / Complexity** | 12种物种，~20步（适中） / 12 species, ~20 steps (moderate) |
| **DFT可行性 / DFT feasibility** | Cu(100)计算成熟 / Cu(100) well-established |
| **实验验证 / Experimental validation** | FE vs U完整数据 / Complete FE vs. U data available |
| **文献空白 / Novelty** | 完全空白 / Completely unexplored |

### 3.2 方法框架 / Methodological Framework

```
         ┌─────────────────────────────────────────────────────┐
         │     第一性原理参数 / First-Principles Parameters      │
         │                                                      │
         │  结构优化    →  CI-NEB过渡态  →  GC-DFT恒电势Ea        │
         │  Structure opt   TS search      Constant-potential Ea │
         │  吸附能/ΔE       所有步骤TS      Ea(U) quadratic fit   │
         │  横向相互作用ε    + 频率(ZPE)                          │
         └─────────────────────┬───────────────────────────────┘
                               │
                               ▼
         ┌─────────────────────────────────────────────────────┐
         │          KMC + MKM 模拟 / KMC + MKM Simulation       │
         │                                                      │
         │  KMC: 50×50正方晶格, BKL算法, 显式空间追踪             │
         │  KMC: 50×50 square lattice, BKL algorithm, explicit   │
         │  MKM: Mean-field ODE, 相同参数集 / same parameter set  │
         └─────────────────────┬───────────────────────────────┘
                               │
                               ▼
         ┌─────────────────────────────────────────────────────┐
         │   四组科学研究 / Four Systematic Studies (A–D)         │
         │                                                      │
         │  A: ε → *NO空间分布 → 选择性 / ε → spatial → FE       │
         │  B: FE vs U (KMC vs MKM) / MKM failure conditions    │
         │  C: 扩散 → Da数 / Diffusion → critical Da number      │
         │  D: 极化曲线 vs 实验 / Polarization vs. experiment     │
         └─────────────────────────────────────────────────────┘
```

---

## 四、技术路线 / 4. Technical Approach

### 4.1 Phase 1: HER on Pt(111) 验证（已完成）/ HER Validation (Completed)

6项验证测试全部通过 / All 6 validation tests passed:
- Langmuir极限 / Langmuir limit ✅
- Tafel斜率 81 mV/dec ✅
- 横向相互作用效应 / Lateral interaction effect ✅
- 扩散效应 / Diffusion effect ✅
- KMC vs MKM: θ_H偏差2.3× / θ_H deviation 2.3× ✅
- 晶格收敛 / Lattice convergence ✅

**关键结论 / Key finding**: H\*-H\*排斥在KMC中产生反聚集，MKM高估覆盖度2.3倍。

H\*-H\* repulsion produces anti-clustering in KMC; MKM overestimates coverage by 2.3×.

### 4.2 Phase 2: NO₃RR on Cu(100) 第一性原理KMC / First-Principles KMC

#### 4.2.1 反应网络 / Reaction Network

**总反应 / Overall reactions:**
```
NO₃⁻ + 9H⁺ + 8e⁻  → NH₃ + 3H₂O       (target, 8e⁻)
2NO₃⁻ + 12H⁺ + 10e⁻ → N₂ + 6H₂O      (byproduct, 10e⁻)
2NO₃⁻ + 10H⁺ + 8e⁻  → N₂O + 5H₂O     (byproduct, 8e⁻)
2H⁺ + 2e⁻           → H₂               (HER competition, 2e⁻)
```

**12种表面物种 / 12 surface species**: \*, \*NO₃, \*NO₂, **\*NO** (key branch point), \*NOH, \*NHOH, \*NH, \*NH₂, \*NH₃, \*N₂O, \*OH, \*H

**关键分支点 / Key branching point — \*NO:**
```
                    *NO
                   ╱    ╲
     [加氢/Hydrogenation]   [N-N偶联/Coupling]
      *NO + H⁺+e⁻           *NO + *NO(adj)
      Ea ≈ 0.08 eV           Ea ≈ 0.05-0.10 eV
           │                      │
        *NOH                   *N₂O
           │                      │
        ...→ NH₃              N₂O / N₂
```

**两条路径势垒几乎相同！选择性由\*NO的空间分布决定——这是KMC的核心价值。**

Both pathways have nearly identical barriers! Selectivity is determined by the spatial distribution of \*NO — this is the core value of KMC.

#### 4.2.2 DFT计算方案 / DFT Calculation Plan

所有KMC参数均从第一性原理计算获取 / All KMC parameters obtained from first principles (ab initio).

**表面模型 / Surface model**:
- Cu(100) fcc, a₀ = 3.615 Å（实验值，直接用于切面，不做bulk优化）
- Cu(100) fcc, a₀ = 3.615 Å (experimental, used directly for slab construction without bulk optimization)
- p(4×4) slab, 4层64Cu, 底2层固定/顶2层弛豫, 15Å真空
- p(4×4) slab, 4 layers × 64 Cu, bottom 2 fixed / top 2 relaxed, 15 Å vacuum

**计算方法 / Methods**: PBE + DFT-D3(BJ), ISMEAR=0/σ=0.05, ENCUT=450 eV, 2×2×1 k-points

**活化能: GC-DFT / Activation energies: GC-DFT**:
- IS和TS各7个变电荷单点（NELECT ± 0.5/1.0/1.5 + 中性）+ VASPsol
- 7 variable-charge single points for each IS and TS (NELECT ± 0.5/1.0/1.5 + neutral) + VASPsol
- Ω(U) = aU² + bU + c → Ea(U) = Ω_TS(U) − Ω_IS(U) + ΔZPE

**过渡态: CI-NEB / Transition states: CI-NEB**:
- 所有步骤（热力学 + PCET热力学类比）均用CI-NEB / All steps via CI-NEB
- PCET建模为 \*A + \*H → \*B（热力学类比）/ PCET modeled as \*A + \*H → \*B (thermodynamic analog)
- 备选: BEP关系 / Fallback: BEP relation if CI-NEB fails

**计算规模 / Computational scale:**

| 类别 / Category | 计算数 / Count | 方法 / Method |
|------|-------|------|
| 结构优化 / Structure optimization | ~55 | IBRION=2 |
| CI-NEB过渡态 / TS search | ~14 | VTST-VASP |
| 频率 / Frequency (ZPE) | ~30 | IBRION=5 |
| GC-DFT变电荷 / Variable charge | ~210 | VASPsol |
| **总计 / Total** | **~310** | |

#### 4.2.3 横向相互作用 / Lateral Interactions

$$\varepsilon(A,B,d) = E(slab{+}A{+}B@d) - E(slab{+}A) - E(slab{+}B) + E(clean)$$

KMC速率修正 / KMC rate correction:
$$k(\text{site}) = k_{\text{base}} \times \exp\left(+\frac{\sum_j \varepsilon(i,j,d_{ij})}{k_BT}\right)$$

**Tier 1（必须 / Essential）**: \*NO–\*NO, \*NO–\*H, \*H–\*H, \*OH–\*NO
**Tier 2（建议 / Recommended）**: \*OH–\*OH, \*OH–\*H, \*N–\*NO, \*N–\*N, \*O–\*NO, \*N–\*H

### 4.3 四组系统性科学研究 / Four Systematic Studies

#### Study A: \*NO空间分布 vs ε / \*NO Spatial Distribution vs. ε

**问题 / Question**: ε(\*NO–\*NO)如何决定\*NO空间分布？/ How does ε determine \*NO patterns?
**方法 / Method**: 固定U=−0.6V, 扫描ε从−0.10到+0.20 eV / Fix U=−0.6V, scan ε from −0.10 to +0.20 eV
**观测量 / Observables**: θ_NO, g(r), cluster size, FE(N₂)/FE(NH₃) vs ε
**预期 / Expected**: ε<0→岛状→N₂↑; ε>0→分散→NH₃↑ / ε<0→islands→N₂↑; ε>0→dispersed→NH₃↑

#### Study B: 选择性 vs 电位 (KMC vs MKM) / Selectivity vs. Potential

**问题 / Question**: 哪些电位下MKM失效？/ At which potentials does MKM fail?
**方法 / Method**: 扫描U: −0.2~−1.0V, KMC与MKM同时运行 / Scan U, run KMC and MKM in parallel
**观测量 / Observables**: FE(NH₃/N₂/N₂O/H₂) vs U, |FE_KMC − FE_MKM| vs U
**预期 / Expected**: 高覆盖度时KMC≠MKM / Deviation at high coverage

#### Study C: \*NO扩散 vs 选择性 / \*NO Diffusion vs. Selectivity

**问题 / Question**: 扩散速率是否控制KMC-MKM偏差？/ Does diffusion control KMC-MKM deviation?
**方法 / Method**: 固定U=−0.6V, 扫描E_diff: 0.05~0.50 eV / Fix U, scan diffusion barrier
**观测量 / Observables**: FE ratio vs E_diff, Da = k_rxn/k_diff
**预期 / Expected**: 快扩散→KMC≈MKM; 慢扩散→KMC≠MKM / Fast→agreement; slow→deviation

#### Study D: 极化曲线 + FE vs 实验 / Polarization + FE vs. Experiment

**问题 / Question**: KMC能否定量预测实验？/ Can KMC quantitatively predict experiments?
**方法 / Method**: U: −0.2~−1.0V, 17点 / 17-point potential scan
**观测量 / Observables**: j(U), FE(U) vs experimental data, Tafel slope

---

## 五、创新点 / 5. Innovation Points

| # | 创新点 / Innovation | 意义 / Significance |
|---|--------|------|
| 1 | **首个NO₃RR KMC研究** / First NO₃RR KMC study | 填补完全空白 / Fills complete gap |
| 2 | **空间效应定量决定选择性** / Spatial effects quantitatively determine selectivity | MKM无法回答 / MKM cannot address |
| 3 | **ε作为选择性调控旋钮** / ε as selectivity tuning knob | DFT→KMC→宏观选择性桥梁 / DFT→KMC→macroscopic bridge |
| 4 | **扩散→MKM失效阈值** / Diffusion → MKM failure threshold | 方法选择实用指南 / Practical guide for method selection |
| 5 | **全参数ab initio + GC-DFT** / All-ab-initio + GC-DFT | 最高精度 / Highest accuracy |
| 6 | **首个通用电催化KMC框架** / First general electrocatalytic KMC framework | Python+Rust开源 / Open-source |

---

## 六、预期成果 / 6. Expected Outcomes

### 论文 / Paper

**暂定标题 / Tentative title**: *Spatial Kinetic Monte Carlo Simulation Reveals Selectivity Control Mechanism of Electrochemical NO₃RR on Cu(100): The Decisive Role of \*NO Spatial Distribution*

**目标期刊 / Target journal**: ACS Catalysis

**核心Figure / Key figures**:
1. 反应网络 + \*NO分支点 / Reaction network + \*NO branching
2. Study A: \*NO快照(吸引/排斥) + g(r) + FE vs ε / Snapshots + g(r) + FE vs ε
3. Study B: FE vs U (KMC vs MKM对比) / FE vs U comparison
4. Study C: Da数 vs mean-field breakdown / Da vs. MKM failure
5. Study D: 极化曲线 + FE vs 实验 / Polarization + FE vs. experiment

### 开源软件 / Open-Source Software

**SPARK v1.0**: Python (~4000行/lines) + Rust (~7000行/lines)

### 扩展方向 / Future Directions

NO₃RR → 尿素电合成 → Nature Catalysis级 / NO₃RR → urea electrosynthesis → Nature Catalysis level

---

## 参考文献 / References

[1] Wei, C. et al. *J. Phys. Chem. Lett.* **16**, 2896 (2025). — CO₂RR KMC
[2] Li, H. et al. *ACS Catal.* **14**, 2696 (2024). — HER KMC
[3] Shou, W. et al. *Nat. Commun.* (2025). — HER multiscale
[4] Jørgensen, M. & Grönbeck, H. *J. Chem. Phys.* **155**, 164701 (2021). — CO₂RR diffusion
[9] Vignola, E. et al. *Catal. Today* **372**, 11 (2021). — KMC vs MKM
[10] Pineda, M. & Stamatakis, M. *J. Chem. Phys.* **156**, 120902 (2022). — KMC review
[20] Liu, J.-X. et al. *ACS Catal.* **9**, 7052 (2019). — NO₃RR MKM
[21] Calle-Vallejo, F. *Faraday Discuss.* (2023). — Cu NO₃RR DFT
[23] *J. Phys. Chem. C* (2025). — NO₃RR on Cu(111) GC-DFT
[24] *ACS Catal.* (2024). — NO₃RR acidic MKM
[25] *Angew. Chem.* (2025). — CuPd N-N vs N-H selectivity
