# Spatial KMC Simulation of Electrochemical NO₃RR: From HER Validation to Selectivity Control

# 电化学NO₃RR的空间KMC模拟：从HER验证到选择性调控

**SPARK Phase 2 Complete Research Plan**

> Last updated: 2026-03-23
> 体系: NO₃RR on Cu(111) | 所有DFT参数从头算(ab initio)
> 决策历程: CO₂RR → 尿素/NO₃RR调研(40+文献) → 最终确认NO₃RR

---

## 目录

1. [引言与核心科学问题](#1-引言与核心科学问题)
2. [文献综述](#2-文献综述)
3. [Phase 1: HER on Pt(111) — 验证（已完成）](#3-phase-1-her-on-pt111--验证已完成)
4. [Phase 2: 体系选择与决策过程](#4-phase-2-体系选择与决策过程)
5. [Phase 2: NO₃RR反应网络与KMC模型](#5-phase-2-no₃rr反应网络与kmc模型)
6. [Phase 2: DFT计算方案](#6-phase-2-dft计算方案)
7. [Phase 2: KMC方法论要点](#7-phase-2-kmc方法论要点)
8. [Phase 2: 研究方案 (Study A–D)](#8-phase-2-研究方案-study-ad)
9. [计算细节与参数](#9-计算细节与参数)
10. [预期科学贡献](#10-预期科学贡献)
11. [任务列表与时间计划](#11-任务列表与时间计划)
12. [参考文献](#12-参考文献)

---

## 1. 引言与核心科学问题

### 1.1 空白：缺乏通用电催化KMC工具

现有通用KMC软件均为热催化设计，不原生支持电化学：

| Software | 语言 | 横向相互作用 | 电化学 | 极化曲线 |
|----------|------|:----------:|:-----:|:------:|
| Zacros | Fortran | CE (多体) | ✗ | ✗ |
| kmos | Python/Fortran | Pairwise | ✗ | ✗ |
| SPPARKS | C++ | Custom | ✗ | ✗ |
| **SPARK** | **Python + Rust** | **Pairwise NN** | **✓ Butler-Volmer** | **✓ Native** |

近年电催化KMC研究均使用一次性代码 [1–4]。**SPARK 是首个通用电催化KMC框架。**

### 1.2 为什么空间KMC对电催化至关重要

Mean-field微动力学(MKM)假设吸附物均匀分布：r = k·θ_A·θ_B。但实际上吸附物形成空间结构（岛状、有序相），导致MKM预测偏差可达**数个数量级** [9]。

| 场景 | MKM | KMC | 已知证据 |
|------|-----|-----|---------|
| 强横向相互作用 | 覆盖度偏差 | 捕捉有序相/岛状 | θ_H偏差2× [2] |
| 扩散受限 | 假设瞬时混合 | 捕捉空间关联 | CO₂RR选择性变化 [4] |
| 双分子选择性步骤 | r∝θ²假设 | 真实pair相关 | NH₃氧化93% NO选择性 [ref] |
| 低覆盖度涨落 | 连续失效 | 天然随机性 | [10] |

### 1.3 核心科学问题

> **Cu(111)上*NO吸附物的空间分布（由横向相互作用、扩散和*H竞争控制）是否决定了电化学NO₃RR中NH₃与N₂的选择性？**

这个问题：
- **实验可观察但机理不明** — 高*NO覆盖度与N₂选择性正相关
- **mean-field无法回答** — 它假设θ_NO均匀
- **KMC独有能力** — 唯一追踪单个吸附物位置及其演化的方法
- **完全文献空白** — 截至2026年3月，不存在任何NO₃RR KMC研究

---

## 2. 文献综述

### 2.1 KMC方法论

| 文献 | 要点 |
|------|------|
| *J. Chem. Phys.* **156**, 120902 (2022) [10] | KMC金标准综述，1p-KMC，刚性问题，加速方法 |
| *Front. Chem.* **7**, 202 (2019) [15] | 实用指南：BKL算法、速率常数、lateral interactions实现 |
| *ACS Catal.* **2**, 2648 (2012) [14] | Graph-theoretical KMC，Zacros奠基 |
| *J. Catal.* **400**, 290 (2021) [17] | MKM在电催化中的局限，推荐KMC替代 |

### 2.2 电催化KMC应用

| 文献 | 体系 | 发现 |
|------|------|------|
| *J. Phys. Chem. Lett.* **16**, 2896 (2025) [1] | CO₂RR on Cu, 178反应 | Cu(100)→C₂, Cu(111)→C₁, 电位依赖RDS |
| *ACS Catal.* **14**, 2696 (2024) [2] | HER on Pt | H₅O₂⁺溶剂 + 动态H覆盖度 |
| *Nat. Commun.* (2025) [3] | HER on Au | 多尺度，全pH极化曲线 |
| *J. Chem. Phys.* **155**, 164701 (2021) [4] | CO₂RR on Cu, 扩散 | 有限扩散破坏mean-field |

### 2.3 NO₃RR文献

| 文献 | 催化剂 | 数据 |
|------|--------|------|
| *ACS Catal.* **9**, 7052 (2019) [20] | 多种TM | MKM参数, volcano plot, 标度关系 |
| *Faraday Discuss.* (2023) [21] | Cu | Grand-canonical DFT, U_L=−0.23V |
| *EES Catalysis* (2024) [22] | Cu基 | NRA3路径Ea，扩散势垒 |
| *J. Phys. Chem. C* (2025) [23] | Cu(111) | GC-DFT完整路径，HER交叉电位 |
| *ACS Catal.* (2024) [24] | 酸性TM | 更新MKM |
| *Angew. Chem.* (2025) [25] | CuPd | N-N vs N-H选择性调控 |
| *Commun. Chem.* (2025) [26] | 综述 | NO₃RR to NH₃最新进展 |

**关键空白：NO₃RR MKM有 [20]，但KMC完全空白。**

### 2.4 横向相互作用

| 文献 | 方法 | 要点 |
|------|------|------|
| Curulla Ferré et al., *ChemPhysChem* **6**, 1009 (2005) [30] | DFT | 成对可加性近似验证 |
| Schneider et al., *J. Catal.* **286**, 88 (2012) [31] | DFT+KMC | O/Pt(111) 覆盖度依赖模型 |
| Hess et al., *J. Comput. Chem.* **40**, 2664 (2019) [32] | 算法 | Cluster expansion高效实现 |

---

## 3. Phase 1: HER on Pt(111) — 验证（已完成）

### 3.1 模型

- **物种**: empty (*), H*
- **晶格**: 50×50 2D, PBC, T = 298K
- **反应**: Volmer fwd/rev (PCET), Tafel (thermal), Heyrovsky (PCET), H*扩散
- **横向相互作用**: ε(H*–H*) = +0.10 eV (排斥) [12]
- **参数来源**: Li 2024 [2], Skúlason 2010 [11], Karlberg 2007 [12]

### 3.2 验证结果

| 测试 | 方法 | 结果 |
|------|------|------|
| 1. Langmuir极限 | 关闭lateral+Tafel | ✅ θ_H = Langmuir解析解 |
| 2. Tafel斜率 | MKM极化曲线 | ✅ 81 mV/dec |
| 3. 横向相互作用 | 有/无ε对比 | ✅ 排斥降低θ_H |
| 4. 扩散效应 | 有/无diffusion | ✅ 扩散增强Tafel速率 |
| 5. KMC vs MKM | 同模型双方法 | ✅ θ_H: MKM 0.76 vs KMC 0.33 (**2.3×偏差**) |
| 6. 晶格收敛 | 10²→10⁴位点 | ✅ 50×50收敛 |

**关键结果**: j₀ = 0.45 mA/cm², Heyrovsky主导。H*–H*排斥在KMC中产生反聚集效应，MKM高估覆盖度2.3倍。

**产出**: `validate_her.py`, `run_polarization_her.py`, `results/`下5张图

---

## 4. Phase 2: 体系选择与决策过程

### 4.1 决策历程

| 日期 | 决策 | 原因 |
|------|------|------|
| 2026-03-21 | 原计划CO₂RR on Cu(100) | C₁/C₂选择性，CO*偶联 |
| 2026-03-22 | 改为多产物体系调研 | 用户不再关注CO₂RR，想体现KMC空间优势 |
| 2026-03-22 | 6项并行调研(40+文献) | 综合评估尿素电合成 vs NO₃RR |
| 2026-03-22 | **最终确认NO₃RR on Cu(111)** | 评分37.5 vs 33.5, 参数最充足, ab initio要求 |

### 4.2 NO₃RR vs 尿素电合成 vs CO₂RR 对比

| 维度 | CO₂RR (原) | 尿素电合成 | **NO₃RR (最终)** |
|------|-----------|----------|----------------|
| 电子转移 | 2–12 | 16 | **8 (NH₃), 10 (N₂)** |
| 物种数 | ~8 | ~25–28 | **~12–15** |
| 步骤数 | ~12 | ~25–35 | **~20–25** |
| 双分子关键步 | CO*+CO* | CO*+NHx | **NO*+NO*** |
| DFT参数 | 充分 | 分散 | **集中且完整** |
| KMC文献 | 有 [1,4] | 空白 | **完全空白** |
| 实现难度 | 中 | 高 | **适中** |
| 评分 | — | 33.5 | **37.5** |

### 4.3 NO₃RR为什么特别适合KMC

1. **N-N偶联是双分子表面反应** — *NO+*NO→*N₂O，需相邻位点，mean-field假设r∝θ²失效
2. **两条路径势垒几乎相同** — *NO加氢Ea≈0.08 eV vs *NO偶联Ea≈0.05–0.10 eV → 选择性由空间分布决定
3. **三方覆盖度竞争** — *NO vs *H vs *OH，空间相分离效应
4. **丰富的实验数据** — FE vs U完整数据可验证KMC预测

---

## 5. Phase 2: NO₃RR反应网络与KMC模型

### 5.1 总反应

```
NO₃⁻ + 9H⁺ + 8e⁻  → NH₃ + 3H₂O       (主产物, 8e⁻)
2NO₃⁻ + 12H⁺ + 10e⁻ → N₂ + 6H₂O      (副产物, 10e⁻)
2NO₃⁻ + 10H⁺ + 8e⁻  → N₂O + 5H₂O     (副产物, 8e⁻)
2H⁺ + 2e⁻           → H₂               (HER竞争, 2e⁻)
```

### 5.2 物种清单 (12种)

| # | 物种 | 描述 | 预期吸附位点 |
|---|------|------|------------|
| 1 | * | 空位 | — |
| 2 | *NO₃ | 吸附硝酸根 | bidentate (O,O) |
| 3 | *NO₂ | 吸附亚硝酸根 | bidentate |
| 4 | **\*NO** | **关键分支中间体** | **fcc hollow** |
| 5 | *NOH | 氮氧氢化物 | atop/bridge |
| 6 | *NHOH | 羟胺中间体 | bridge/hollow |
| 7 | *NH | 亚胺 | fcc hollow |
| 8 | *NH₂ | 氨基 | bridge |
| 9 | *NH₃ | 吸附氨 | atop |
| 10 | *N₂O | N-N偶联产物 | atop/bridge |
| 11 | *OH | 吸附羟基 | fcc hollow/bridge |
| 12 | *H | 吸附氢 (HER竞争) | fcc hollow |

### 5.3 基元步骤 (~20–25步)

**Phase I: 脱氧 (NO₃⁻ → *NO)**

| Step | 反应 | 类型 | Ea (eV) |
|------|------|------|---------|
| R1 | NO₃⁻(aq) + * → *NO₃ | 吸附 | ~0 |
| R2 | *NO₃ + H⁺+e⁻ → *NO₂ + *OH | PCET | ~0.3–0.5 |
| R3 | *NO₂ → NO₂⁻(aq) | 脱附 | — (2e⁻副产物) |
| R4 | *NO₂ + H⁺+e⁻ → *NO + *OH | PCET | ~0.3–0.5 |
| R5 | *OH + H⁺+e⁻ → H₂O + * | PCET | ~0.1 |

**Phase II: *NO加氢 → NH₃ (NRA3路径)**

| Step | 反应 | 类型 | Ea (eV) |
|------|------|------|---------|
| R6 | **\*NO + H⁺+e⁻ → \*NOH** | **PCET** | **~0.08** (PLS入口) |
| R7 | *NOH + H⁺+e⁻ → *NHOH | PCET | ~0.3 |
| R8 | *NHOH + H⁺+e⁻ → *NH + H₂O | PCET | ~0.23 |
| R9 | *NH + H⁺+e⁻ → *NH₂ | PCET | ~0.3 |
| R10 | *NH₂ + H⁺+e⁻ → *NH₃ | PCET | ~0.2 |
| R11 | *NH₃ → NH₃(aq) + * | 脱附 | ~0.37 |

**Phase III: N-N偶联 → N₂/N₂O (KMC核心)**

| Step | 反应 | 类型 | Ea (eV) |
|------|------|------|---------|
| R12 | ***NO + \*NO(adj) → \*N₂O + \*** | **Thermal, 双分子** | **~0.05–0.10** |
| R13 | *N₂O → N₂O(g) + * | 脱附 | ~0.3 |
| R14 | *N + *N → *N₂ (或 N₂(g)) | Thermal | ~0.45 |

**Phase IV: HER竞争**

| Step | 反应 | 类型 | Ea (eV) |
|------|------|------|---------|
| R15 | H₂O + * + e⁻ → *H + OH⁻ | PCET (Volmer) | ~0.7 |
| R16 | *H + *H → H₂ + 2* | Thermal (Tafel) | ~0.25 |
| R17 | *H + H⁺+e⁻ → H₂ + * | PCET (Heyrovsky) | ~0.7 |

**扩散**

| Step | 反应 | 类型 | Ea (eV) |
|------|------|------|---------|
| R18 | *NO扩散 (fcc→fcc) | Thermal | ~0.1–0.3 |
| R19 | *H扩散 | Thermal | ~0.1–0.2 |

### 5.4 关键分支点：*NO

```
                    *NO
                   ╱    ╲
          [加氢路径]      [偶联路径]
         *NO + H⁺+e⁻     *NO + *NO(adj)
         Ea ≈ 0.08 eV     Ea ≈ 0.05-0.10 eV
              │                 │
           *NOH              *N₂O
              │                 │
           ...→ NH₃          N₂O / N₂
```

**两条路径势垒几乎相同！** → 选择性不由动力学势垒决定，而由*NO的空间分布决定 → **这是KMC的核心价值。**

### 5.5 横向相互作用

**Tier 1 (必须)**

| 对 | 预期ε (eV) | 原因 |
|----|-----------|------|
| *NO–*NO | +0.10~0.30 | N-N偶联选择性核心 |
| *NO–*H | +0.05~0.10 | 加氢vs偶联空间竞争 |
| *H–*H | +0.05~0.10 | HER竞争, H有序化 |
| *OH–*NO | +0.10~0.20 | OH累积对NO的影响 |

**Tier 2 (建议)**

| 对 | 原因 |
|----|------|
| *OH–*OH | OH累积毒化 |
| *OH–*H | OH清除动力学 |
| *N–*NO | 替代N-N偶联路径 |
| *N–*N | 直接N₂形成 |
| *O–*NO | 脱氧副产物影响 |
| *N–*H | N加氢动力学 |

**KMC中的使用**: k(site) = k_base × exp(+E_lat/kBT)，E_lat = Σε(邻居)

### 5.6 实验验证数据

| U (V vs RHE) | FE(NO₂⁻) | FE(NH₃) | 主产物 |
|-------------|-----------|---------|-------|
| +0.10 | ~75% | <10% | NO₂⁻ |
| −0.10 | ~50% | ~30% | 混合 |
| −0.20 | ~20% | ~70% | NH₃ |
| −0.30 | <10% | ~88% | NH₃(最大) |
| −0.40 | <5% | ~85% | NH₃(HER开始) |

Tafel斜率 ~120 mV/dec。

---

## 6. Phase 2: DFT计算方案

> **理论公式和VASP详细操作见**: `docs/DFT_workflow_operational.md` (操作手册)

### 6.1 核心公式

```
ε(A,B,d) = E(slab+A+B@d) − E(slab+A) − E(slab+B) + E(clean_slab)
```

- ε > 0 → 排斥 → 加速反应/脱附
- ε < 0 → 吸引 → 减慢脱附
- 需4个DFT计算/对/距离（单吸附和clean可复用）

### 6.2 Cu(111) slab模型

| 参数 | 值 |
|------|---|
| 表面 | Cu(111), fcc |
| 层数 | 4层 (底2固定, 顶2弛豫) |
| 真空层 | 15 Å |
| 晶格常数 | PBE优化 (~3.63 Å) |
| 最近邻距离 | a₀/√2 ≈ 2.556 Å |
| 偶极修正 | LDIPOL=.TRUE., IDIPOL=3 |

| 超胞 | 原子数 | 可测距离 | k-points | 用途 |
|------|-------|---------|----------|------|
| p(3×3) | 36 Cu | 1NN+2NN | 3×3×1 | 位点测试 |
| **p(4×4)** | **64 Cu** | **1NN+2NN+3NN** | **2×2×1** | **精确值** |

### 6.3 完整DFT计算清单 (~78–86个)

| Phase | 内容 | 计算数 | 超胞 | 预估 |
|-------|------|-------|------|------|
| **Phase 0** | Bulk Cu + 8种气相分子 | 9 | 小 | ~4 CPU·h |
| **Phase 1** | 位点偏好: 10–12物种×4位点 | ~41–49 | p(3×3) | ~60–90 CPU·h |
| **Phase 2** | 精确单吸附: 7关键物种 | 8 | p(4×4) | ~30 CPU·h |
| **Phase 3** | Tier 1成对: 4对×2距离 | 8 | p(4×4) | ~30 CPU·h |
| **Phase 4** | Tier 2成对: 6对×2距离 | 12 | p(4×4) | ~45 CPU·h |
| **总计** | | **~78–86** | | **~170–200 CPU·h** |

64核节点预估 **~3–5天**。

### 6.4 VASP关键参数

```fortran
! 通用slab计算
PREC=Accurate  ENCUT=450  EDIFF=1E-5  EDIFFG=-0.02
IBRION=2  NSW=300  ISIF=2  ISMEAR=1  SIGMA=0.15
ISPIN=2  LDIPOL=.TRUE.  IDIPOL=3  LREAL=Auto

! 气相分子
ISMEAR=0  SIGMA=0.01  K-points: Gamma only
! Box: 15×16×17 Å (三方向不同避免对称性)

! POTCAR: Cu PAW_PBE + N + O + H
! ENCUT=450 ≥ max(ENMAX)×1.13
```

### 6.5 过渡态搜索

**热力学步骤** — CI-NEB或Dimer方法:
- *NO+*NO→*N₂O (最关键! 决定N₂选择性)
- *H+*H→H₂ (Tafel)
- *NO扩散, *H扩散

**PCET步骤** — CHE方法:
```
Ea(U) = Ea(U=0) + β·e·U,  β ≈ 0.5
```
可先用文献Ea₀值 [20,23]，后续替换为自己的DFT。

### 6.6 数据提取

```bash
# 从OUTCAR提取能量
grep "energy(sigma->0)" OUTCAR | tail -1 | awk '{print $NF}'

# 计算ε (同种, 如*NO+*NO)
ε = E(co-ads) − 2×E(single) + E(clean)
```

### 6.7 自洽性检查

| 检查 | 标准 |
|------|------|
| a₀ vs 实验3.615 Å | < 2% |
| p(3×3) vs p(4×4) 单吸附 | ΔE < 30 meV |
| k-points 3×3×1 vs 5×5×1 | ΔE < 10 meV |
| 等价位点一致性 | < 5 meV |
| entropy T*S | < 1 meV/atom |

---

## 7. Phase 2: KMC方法论要点

### 7.1 速率表达式

**PCET步骤 (Butler-Volmer)**:
$$k(U) = \frac{k_BT}{h}\exp\left(-\frac{E_a + \beta \cdot e \cdot (U-U_0)}{k_BT}\right)$$

**热力学步骤 (TST)**:
$$k = \frac{k_BT}{h}\exp\left(-\frac{E_a}{k_BT}\right)$$

**横向相互作用修正**:
$$k(\text{site}) = k_\text{base} \times \exp\left(+\frac{E_\text{lat}}{k_BT}\right)$$

### 7.2 吸附与脱附

**溶液相吸附**: k_ads = k⁰·(c/c°)，KMC自动处理(1−θ)因子（仅在空位发生）

**电化学脱附**: 本质是PCET步骤，用Butler-Volmer处理

**热脱附**: k_des = ν·exp(−|BE|/kBT)，ν ~ 10¹²–10¹³ s⁻¹

### 7.3 多产物追踪

**TOF**: TOF_i = N_desorption_i / (N_sites × t_simulation)

**偏电流密度**: j_i = n_i·e·ρ_sites·TOF_i

**法拉第效率**: FE_i = n_i·TOF_i / Σ(n_j·TOF_j)

| 产物 | n_i (电子数) |
|------|------------|
| NH₃ | 8 |
| N₂ | 10 |
| N₂O | 8 |
| H₂ | 2 |

### 7.4 统计要求

- 每种产物至少 **100–1000个脱附事件** (10%相对误差)
- 少量产物(5% selectivity)需约10,000总事件
- 建议5–10个独立副本(不同随机种子)平均
- 稳态判定：迭代batch-means + 覆盖度收敛

### 7.5 KMC vs MKM已知偏差

| 体系 | MKM预测 | KMC结果 | 偏差原因 |
|------|--------|--------|---------|
| CO氧化/RuO₂ | 活性高 | **低数个数量级** | 空间关联 |
| NH₃氧化/RuO₂ | NO选择性低 | **93% NO** (实验95%) | 空间效应主导 |
| HER/Pt(111) | θ_H=0.76 | **θ_H=0.33** | 横向排斥→反聚集 |
| CO₂RR/Cu | — | 扩散改变C₂比例 | 有限扩散破坏mean-field |

---

## 8. Phase 2: 研究方案 (Study A–D)

### Study A: *NO空间分布 vs 横向相互作用强度

**问题**: ε(*NO–*NO)的符号和大小如何决定*NO空间分布模式？

**方法**: 固定U=−0.6V, 扫描ε(*NO–*NO)从−0.10到+0.20 eV (步长0.02 eV)

**观测量**:
- θ_NO 覆盖度
- g_NO-NO(r) 径向分布函数
- *NO cluster size分布
- N-N偶联速率 vs ε
- FE(N₂)/FE(NH₃) vs ε

**预期**: ε<0(吸引)→岛状聚集→N₂增多; ε>0(排斥)→分散→NH₃增多

### Study B: NH₃/N₂选择性 vs 电位 (KMC vs MKM)

**问题**: 哪些电位下mean-field失效？KMC与MKM的选择性偏差有多大？

**方法**: 扫描U从−0.2到−1.0V (步长50mV), 同时运行KMC和MKM

**观测量**:
- FE(NH₃), FE(N₂), FE(N₂O), FE(H₂) vs U — KMC vs MKM
- θ_NO, θ_H, θ_OH vs U
- j_total vs U — 极化曲线
- |FE_KMC − FE_MKM| vs U

**预期**: 低过电位KMC≈MKM; 高过电位KMC≠MKM(空间相分离)

### Study C: *NO扩散速率对选择性的影响

**问题**: NO*扩散的快慢是否决定KMC与MKM是否一致？

**方法**: 固定U=−0.6V, 扫描E_diff(NO)从0.05到0.50 eV (步长0.05 eV)

**观测量**:
- FE(N₂)/FE(NH₃) vs E_diff
- g_NO-NO(r) vs E_diff
- |KMC−MKM|/MKM vs E_diff → mean-field breakdown阈值
- Damköhler数 Da = k_reaction/k_diffusion

**预期**: 快扩散→KMC≈MKM; 慢扩散→空间不均匀→KMC≠MKM

### Study D: 极化曲线, FE vs U, 偏电流密度

**问题**: SPARK能否定量预测实验可观测量？

**方法**: 高精度极化曲线(U: −0.2~−1.0V, 17点), 与实验数据对比

**观测量**:
- j_total vs U — KMC极化曲线
- j_NH₃, j_N₂, j_N₂O, j_H₂ — 各产物偏电流密度
- FE vs U — 与实验对比 [22,23]
- Tafel斜率, DRC分析

---

## 9. 计算细节与参数

### 9.1 软件

- **SPARK v0.3.0** (Python 3,943行 + Rust 7,020行)
- Rust引擎: Fenwick tree O(log N) 位点选择
- Python引擎: 原型验证, ~1500 steps/s
- MKM: scipy.fsolve (Python) / RK4 (Rust)

### 9.2 模拟参数

| 参数 | Phase 1 (HER) | Phase 2 (NO₃RR) |
|------|--------------|-----------------|
| 晶格 | 50×50 | 50×50 (或100×100 Rust) |
| 温度 | 298 K | 298 K |
| 电位范围 | −0.5~0.0 V | −0.2~−1.0 V vs RHE |
| 平衡化 | 10⁶步 | 10⁶步 |
| 生产 | 10⁶步 | 5×10⁶步 |
| 采样 | 每10⁴步 | 每10⁴步 |
| 副本数 | 5 | 5–10 |

### 9.3 分析工具

- θ_i(t) 覆盖度时间序列
- TOF_i 各产物转换频率
- g(r) pair相关函数
- P(n) cluster大小分布
- 晶格快照可视化
- j(U), FE(U) 极化曲线/法拉第效率

---

## 10. 预期科学贡献

| # | 贡献 | 新颖性 |
|---|------|--------|
| 1 | **首个NO₃RR KMC研究** | 完全文献空白 |
| 2 | **空间效应定量决定NH₃/N₂选择性** | MKM无法回答 |
| 3 | **ε(*NO–*NO)作为选择性调控旋钮** | DFT横向相互作用→宏观选择性 |
| 4 | **扩散速率阈值: mean-field何时失效** | 实用指南 |
| 5 | **首个通用电催化KMC框架** | 软件贡献 |
| 6 | **所有参数ab initio, 非文献拼凑** | 自洽DFT参数集 |

**目标期刊**: ACS Catalysis (首篇NO₃RR KMC)

**扩展方向**: NO₃RR → 尿素电合成 (增量开发: +CO₂RR子路径 + C-N偶联 → Nature Catalysis级)

---

## 11. 任务列表与时间计划

> 详细任务跟踪见 `.autopilot/TODO_TRACKER.json`

| ID | 任务 | 状态 | 依赖 |
|----|------|------|------|
| T1–T4 | Phase 1: HER参数、实现、验证、极化曲线 | ✅ 完成 | — |
| T5 | Phase 2: 体系调研(40+文献) | ✅ 完成 | — |
| T12 | Phase 2: 研究计划文档 | ✅ 完成 | — |
| **T6** | **DFT: 吸附位点+吸附能+横向相互作用ε** | **⏳ 当前** | T5 |
| T6b | DFT: 过渡态搜索(活化能Ea) | ⏳ 待开始 | T6 |
| T7 | Code: 实现NO₃RR KMC模型 | ⏳ | T6,T6b |
| T8 | Study A: ε vs 空间分布 | ⏳ | T7 |
| T9 | Study B: 选择性 vs 电位 | ⏳ | T7 |
| T10 | Study C: 扩散 vs mean-field breakdown | ⏳ | T7 |
| T11 | Study D: 极化曲线 + FE vs U | ⏳ | T8–T10 |
| T14 | 论文撰写 | ⏳ | T11 |

**当前阻塞**: T6 (用户DFT计算, ~78个VASP作业)

**依赖链**: T6→T6b→T7→{T8,T9,T10}→T11→T14

---

## 12. 参考文献

### 电催化KMC研究

[1] Wei, C. et al. Voltage-Dependent CO₂RR Mechanism Unveiled by KMC. *J. Phys. Chem. Lett.* **16**, 2896–2904 (2025). DOI: 10.1021/acs.jpclett.4c03426

[2] Li, H. et al. First-Principles-Based KMC Model of HER. *ACS Catal.* **14**, 2696–2708 (2024). DOI: 10.1021/acscatal.3c04588

[3] Shou, W. et al. HER Polarization Curves through Multiscale Simulations. *Nat. Commun.* (2025).

[4] Jørgensen, M. & Grönbeck, H. Effects of Surface Diffusion in CO₂RR on Cu. *J. Chem. Phys.* **155**, 164701 (2021). DOI: 10.1063/5.0065888

### Mean-field vs KMC

[9] Vignola, E. et al. Evaluating KMC and MKM with Lateral Interactions. *Catal. Today* **372**, 11 (2021). DOI: 10.1016/j.cattod.2020.10.018

### KMC方法论

[10] Pineda, M. & Stamatakis, M. KMC for Heterogeneous Catalysis: Fundamentals. *J. Chem. Phys.* **156**, 120902 (2022). DOI: 10.1063/5.0083251

[11] Skúlason, E. et al. Modeling HER on the Basis of DFT. *J. Phys. Chem. C* **114**, 18182 (2010).

[12] Karlberg, G. S. et al. Cyclic Voltammograms for H on Pt from First Principles. *Phys. Rev. Lett.* **99**, 126101 (2007).

[14] Stamatakis, M. & Vlachos, D. G. Unraveling Complexity via KMC. *ACS Catal.* **2**, 2648 (2012).

[15] Andersen, M. et al. A Practical Guide to Surface KMC. *Front. Chem.* **7**, 202 (2019). DOI: 10.3389/fchem.2019.00202

[17] Exner, K. S. Microkinetic Modeling in Electrocatalysis: Limitations. *J. Catal.* **400**, 290 (2021).

### NO₃RR文献

[20] Liu, J.-X. et al. Activity and Selectivity Trends in NO₃RR. *ACS Catal.* **9**, 7052–7064 (2019). DOI: 10.1021/acscatal.9b02179

[21] Calle-Vallejo, F. Why Copper Catalyzes NO₃RR to NH₃. *Faraday Discuss.* (2023). DOI: 10.1039/D2FD00145D

[22] Cu-based NO₃RR Review. *EES Catalysis* (2024).

[23] New Insights into NO₃RR on Cu(111). *J. Phys. Chem. C* (2025). DOI: 10.1021/acs.jpcc.5c02461

[24] Nitrate Reduction Modeling under Acidic Conditions. *ACS Catal.* (2024). DOI: 10.1021/acscatal.4c06394

[25] CuPd N-N vs N-H Selectivity. *Angew. Chem.* (2025). DOI: 10.1002/anie.202524218

[26] Recent Advances in NO₃RR Mechanisms. *Commun. Chem.* (2025). DOI: 10.1038/s42004-025-01864-w

### 尿素电合成 (备选扩展方向)

[27] Cu Surface EDL Study. *Nat. Commun.* (2024). DOI: 10.1038/s41467-024-45522-6

[28] CuWO₄ C-N Coupling. *Nat. Commun.* (2023). DOI: 10.1038/s41467-023-40273-2

### 横向相互作用方法论

[30] Curulla Ferré et al. Testing Pairwise Additive Potential Approximation. *ChemPhysChem* **6**, 1009 (2005). DOI: 10.1002/cphc.200400399

[31] Wu, Schmidt, et al. Coverage-Dependent Model of Pt-Catalyzed NO Oxidation. *J. Catal.* **286**, 88 (2012). DOI: 10.1016/j.jcat.2011.10.020

[32] Hess, Sander, Stamatakis. Cluster Expansion in Surface KMC. *J. Comput. Chem.* **40**, 2664 (2019). DOI: 10.1002/jcc.26041

---

*本文档合并自: research_plan_unified.md, analysis_NO3RR_for_KMC.md, phase2_deep_research_report.md, DFT_lateral_interactions_guide.md, DFT_workflow_operational.md, multi_product_KMC_analysis.md, adsorption_desorption_KMC_review.md*

*DFT操作手册详见: docs/DFT_workflow_operational.md*

*Last updated: 2026-03-23*
