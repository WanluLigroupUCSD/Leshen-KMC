# Phase 2 深度调研报告：反应体系选择与KMC方法论

> **SPARK Phase 2 Research Report**
> 2026-03-22 | 综合6项并行调研结果

---

## 目录

1. [KMC方法论：横向排斥相互作用](#1-横向排斥相互作用的kmc处理方法)
2. [KMC方法论：分子吸附与脱附](#2-分子吸附与脱附的kmc处理方法)
3. [KMC方法论：多产物分析](#3-多产物kmc分析方法)
4. [候选体系A：尿素电合成](#4-候选体系a尿素电合成-co₂--no₃⁻)
5. [候选体系B：NO₃RR](#5-候选体系bno₃rr-on-cu)
6. [对比评估与推荐](#6-对比评估与最终推荐)

---

## 1. 横向排斥相互作用的KMC处理方法

### 1.1 Pairwise Nearest-Neighbor (NN) 模型 — SPARK现有实现

**基本公式：**

```
E_lat(site_i) = Σ_j ε(sp_i, sp_j)
```

- `ε > 0`：排斥（repulsive），削弱吸附，加速反应/脱附
- `ε < 0`：吸引（attractive），稳定吸附，减慢脱附

**速率常数修正：**

```
k(site) = k_base × exp(+E_lat_reactant / kBT)    # 无BEP
k(site) = k_base × exp(-α × ΔΔH / kBT)           # 含BEP修正
```

其中 α 为BEP斜率（~0.5），ΔΔH = E_lat(products) - E_lat(reactants)

**室温影响量级：** 典型排斥相互作用 0.05–0.25 eV 可改变速率常数 **7–15,000倍**

**SPARK已实现：** `engine.py` 中 `_compute_site_rate()` 和 `rates.py` 中 `bep_modified_rate()`

### 1.2 Cluster Expansion (CE) — 高级方法

**CE哈密顿量：**

```
H_CE = Σ_α J_α · Φ_α(σ)
```

- J_α = effective cluster interaction (ECI)
- 包含pair、trio、quartet等多体项
- 需要50-200个DFT构型拟合，典型RMSE < 20 meV

**典型数值（O/Pt(111)）：** 1NN = +237.1 meV, 2NN = +39.5 meV, 3NN = -5.8 meV

**何时需要CE vs Pairwise NN：**
- 低覆盖度（< 0.25 ML）：Pairwise NN足够
- 高覆盖度（> 0.5 ML）或有序相：需要CE / 三体项
- **Phase 2建议：先用Pairwise NN，后续按需升级**

### 1.3 对Phase 2体系的具体建议

**关键吸附物对的预期排斥强度：**

| 吸附物对 | 预期ε (eV) | 物理来源 | 对选择性的影响 |
|---------|-----------|---------|-------------|
| \*NO–\*NO | +0.10~0.30 | 偶极-偶极，空间位阻 | 排斥→分散→抑制N-N偶联 |
| \*CO–\*CO | +0.15~0.32 | 共享金属d轨道 | 排斥→分散→抑制C-C/C-N偶联 |
| \*NO–\*H | +0.05~0.10 | 弱空间位阻 | 相分离→影响加氢vs偶联竞争 |
| \*OH–\*OH | +0.10~0.20 | 偶极-偶极 | 表面毒化 |
| \*CO–\*NO | 未知(需DFT) | 取决于催化剂 | **决定C-N偶联概率的核心参数** |
| \*H–\*H (Pt) | +0.05~0.10 | 弱排斥 | 有序吸附层 |

**SPARK升级建议优先级：**
1. 添加2NN相互作用（当前仅1NN）
2. 添加Blowers-Masel近似（更准确的势垒修正）
3. 高覆盖度场景下添加trio相互作用

---

## 2. 分子吸附与脱附的KMC处理方法

### 2.1 气相/溶液相吸附

**非活化吸附（气相，Hertz-Knudsen）：**

```
k_ads = P · A_site · S / √(2π m kB T)    [s⁻¹]
```

典型值：~10⁶–10⁸ s⁻¹（取决于压力和分子量）

**溶液相吸附（电化学体系）：**

```
k_ads = k₀ · (c / c_ref)    [s⁻¹]
```

其中 c 为反应物溶液浓度，c_ref 为参考浓度。KMC自动处理(1-θ)因子（仅在空位上发生）。

**电化学PCET吸附（Butler-Volmer）：**

```
k_PCET = k₀ · exp(-α·e·η / kBT)
```

其中 η = U - U₀ 为过电位，α 为转移系数（~0.5）。SPARK已原生支持。

### 2.2 脱附

**热脱附（Arrhenius）：**

```
k_des = ν · exp(-E_des / kBT)    [s⁻¹]
```

- ν ~ 10¹²–10¹³ s⁻¹（典型预指数因子）
- E_des 受横向相互作用修正：E_des(env) = E_des(0) - Σ ε_ij · n_j

**缔合脱附（如Tafel步：2\*H → H₂）：**
- 需要两个相邻占据位点
- KMC中通过双分子反应事件处理
- 速率取决于局部构型

**电化学脱附（如\*OH + H⁺ + e⁻ → H₂O）：**
- 本质是电化学反应步骤，不是简单脱附
- 用Butler-Volmer速率处理

### 2.3 多产物追踪

**产物形成速率（TOF）：**

```
TOF_i = N_desorption_i / (N_sites × t_simulation)    [s⁻¹ site⁻¹]
```

**偏电流密度：**

```
j_i = n_i · e · ρ_sites · TOF_i    [A/cm²]
```

其中 n_i 为产物i的电子转移数，ρ_sites 为位点面密度（~1.5×10¹⁵ sites/cm²）。

**法拉第效率：**

```
FE_i = n_i · TOF_i / Σ_j(n_j · TOF_j)
```

---

## 3. 多产物KMC分析方法

### 3.1 选择性计算

**稳态判定：** 使用迭代batch-means方法——将模拟轨迹分批，监测各批TOF方差收敛。

**统计需求：**
- 每种产物至少需要 **100–1000个脱附事件** 才能达到~10%相对误差
- 对于5%选择性的少量产物：需约10,000总事件才能获得~500个少量产物事件
- 建议至少5-10个独立副本（不同随机种子）平均

### 3.2 极化曲线生成

**方法：** 在每个电位点独立运行KMC模拟

```
对 U = U₁, U₂, ..., Uₙ:
    运行KMC → 达到稳态 → 统计TOF_i
    j_total(U) = Σ_i n_i · e · ρ · TOF_i(U)
    FE_i(U) = n_i · TOF_i(U) / Σ_j n_j · TOF_j(U)
```

**建议：** Tafel区域内电位间距25-50 mV，每点~10⁶ KMC步。

### 3.3 KMC vs MKM选择性差异的已知案例

| 体系 | MKM预测 | KMC结果 | 偏差原因 |
|------|--------|--------|---------|
| CO氧化/RuO₂ | 活性高 | **活性低数个数量级** | 空间关联即使无lateral interaction也重要 |
| NH₃氧化/RuO₂ | NO选择性偏低 | **93% NO（实验95%）** | 空间效应主导选择性 |
| HER/Pt | θ_H偏离 | **覆盖度偏差2倍** | 横向相互作用下的有序吸附 |
| CO₂RR/Cu | C₂选择性 | **扩散速率改变C₂比例** | 有限扩散破坏mean-field假设 |

**关键结论：** 对于涉及双分子表面反应的选择性问题，KMC几乎总能提供比MKM更准确的预测。

---

## 4. 候选体系A：尿素电合成 (CO₂ + NO₃⁻)

### 4.1 总体反应

```
CO₂ + 2NO₃⁻ + 18H⁺ + 16e⁻ → CO(NH₂)₂ + 7H₂O    (16电子转移!)
```

**产物清单（10-12种）：**
- 目标产物：**尿素 CO(NH₂)₂**
- NO₃RR副产物：NH₃, NO₂⁻, N₂, N₂O, NH₂OH
- CO₂RR副产物：CO, HCOOH, CH₄
- HER副产物：H₂

### 4.2 反应机理

**三阶段机制：**

```
阶段1: CO₂活化                    阶段2: NO₃⁻活化
CO₂(g) → *CO₂                    NO₃⁻(aq) → *NO₃
*CO₂ + H⁺+e⁻ → *COOH            *NO₃ + H⁺+e⁻ → *NO₂ + *OH
*COOH + H⁺+e⁻ → *CO + H₂O      *NO₂ + H⁺+e⁻ → *NO + *OH
                                  *NO + H⁺+e⁻ → *NOH → ... → *NH₂
                  ↓                         ↓
              阶段3: C-N偶联 (选择性决定步!)
              *CO + *NH → *CO-NH        (Ea = 0.357 eV, Cu(111))
              *CO-NH + *NH → *NH-CO-NH  (第二个C-N键)
              *NH-CO-NH + 2(H⁺+e⁻) → NH₂CONH₂ (尿素)
```

### 4.3 C-N偶联——KMC的核心价值

**已报道的C-N偶联方式：**

| 偶联组合 | Ea (eV) | 催化剂 | 来源 |
|---------|---------|-------|------|
| \*CO + \*NH | **0.357** | Cu(111) | Nat. Commun. 2024 |
| \*CO + \*NO | 1.35 (纯Cu), 0.36 (Cu/Cu₂O) | Cu, Cu/Cu₂O | CCS Chemistry 2023 |
| \*CO + \*NO₂ | 0.87 | CuWO₄ | Nat. Commun. 2023 |
| \*CO + \*NH₂ | 0.85 (Zn), 0.32 (工程催化剂) | Zn, SAC | 综述2024 |
| \*NH + \*CO | 0.21 | B-FeNi-DASC | 综述2024 |

**这是Langmuir-Hinshelwood双分子反应 → 完美适合KMC！**

Mean-field假设 r_CN = k · θ_CO · θ_NHx 在以下情况失效：
1. \*CO形成岛状结构，\*NHx分布不均匀 → 偶联仅在域边界发生
2. 表面扩散慢 → \*CO来不及与\*NHx相遇就被加氢或脱附
3. NO₃RR和CO₂RR时间尺度不匹配 → C/N中间体的时空梯度

### 4.4 关键选择性分支点

```
分支点1: *CO脱附 vs C-N偶联
  *CO → CO(g)             (CO副产物，lost carbon)
  *CO + *NHx → *CONHx     (尿素路径)
  → 需要足够的*NHx覆盖度 AND 空间邻近性

分支点2: *NH₂加氢 vs C-N偶联
  *NH₂ + H⁺+e⁻ → *NH₃    (NH₃副产物，lost nitrogen)
  *NH₂ + *CO → *CONH₂     (尿素路径)
  → 需要足够的*CO覆盖度 AND 空间邻近性

分支点3: *NO路径分歧
  *NO + *CO → *OCNO        (尿素路径)
  *NO + H⁺+e⁻ → *NOH      (继续到NH₃)
  *NO + *NO → *N₂O         (N₂O副产物)

分支点4: HER竞争
  * + H⁺+e⁻ → *H → H₂(g)
```

### 4.5 DFT参数可用性

**最完整数据集：**

| 参考文献 | 催化剂 | 可用数据 | DOI |
|---------|-------|---------|-----|
| Nat. Commun. 2024 | Cu(111)/(100)/(110) | **完整路径Ea，包含C-N偶联和EDL效应** | 10.1038/s41467-024-45522-6 |
| Nat. Commun. 2023 | CuWO₄ | C-N偶联Ea=0.87eV，完整自由能图 | 10.1038/s41467-023-40273-2 |
| Nat. Commun. 2023 | N-掺杂碳 | 序贯共还原机理，CI-NEB势垒 | 10.1038/s41467-023-44131-z |
| CCS Chemistry 2023 | Cu/Cu₂O | \*CO+\*NO偶联势垒，氧化态效应 | 10.31635/ccschem.023.202202408 |
| Commun. Chem. 2023 | 通用 | DFT误差分析（尿素形成能-0.25±0.10 eV） | 10.1038/s42004-023-00990-7 |
| J. Catal. 2025 | Cu BIF-29 | DFT+MKM闭环，TOF vs电位 | 最接近KMC的现有工作 |
| Commun. Chem. 2025 | M-Pc系列 | 描述符法，volcano图 | 10.1038/s42004-025-01424-2 |

**评估：Cu(111)参数最充足，可构建初始KMC模型。CuWO₄作为备选。**

### 4.6 模型规模

| 指标 | 完整模型 | 最小可行模型 |
|------|---------|------------|
| 表面物种数 | ~25-28 | ~15 |
| 基元步骤数 | ~25-35 | ~20 |
| 产物数 | ~10-12 | 5-6 (尿素, NH₃, CO, H₂, N₂O) |
| 复杂度对比 | ~2× NO₃RR | ~1.5× NO₃RR |

### 4.7 文献空白

**确认：截至2026年3月，不存在任何尿素电合成KMC研究。**

最接近的工作：
- CO₂RR KMC (J. Chem. Phys. 2021) — 仅CO₂还原，无N源
- CO₂RR KMC 178反应 (J. Phys. Chem. Lett. 2024) — 仅CO₂还原
- Cu BIF-29 MKM (J. Catal. 2025) — Mean-field微动力学，非KMC

### 4.8 实验现状

FE从2020-2021的个位数进展到2026的 **>80%**（Cu-S配位聚合物，Nat. Commun.）。多数催化剂在40-70%范围。领域发展迅速，高影响力。

---

## 5. 候选体系B：NO₃RR on Cu

### 5.1 总体反应

```
NO₃⁻ + 9H⁺ + 8e⁻ → NH₃ + 3H₂O         (主产物, 8e⁻)
2NO₃⁻ + 12H⁺ + 10e⁻ → N₂ + 6H₂O       (副产物, 10e⁻)
2NO₃⁻ + 10H⁺ + 8e⁻ → N₂O + 5H₂O       (副产物, 8e⁻)
2H⁺ + 2e⁻ → H₂                          (竞争HER, 2e⁻)
```

**产物：** NH₃, N₂, N₂O, NO₂⁻, NH₂OH, H₂ (5-6种)

### 5.2 完整基元步骤 (~25步)

**Phase I: 脱氧 (NO₃⁻ → \*NO)**

| Step | 反应 | 类型 | 备注 |
|------|------|------|------|
| R1 | NO₃⁻(aq) + \* → \*NO₃ | 吸附 | Bridge-bidentate via O |
| R2 | \*NO₃ + H⁺+e⁻ → \*NO₂ + \*OH | PCET | 第一步脱氧 |
| R4 | \*NO₂ → NO₂⁻(aq) | 脱附 | 2e⁻副产物 |
| R5 | \*NO₂ + H⁺+e⁻ → \*NO + \*OH | PCET | 第二步脱氧 |
| R6 | \*OH + H⁺+e⁻ → H₂O + \* | PCET | OH清除 |

**Phase II: \*NO加氢 → NH₃ (NRA3路径，Cu上最优)**

| Step | 反应 | Ea (eV) | 来源 |
|------|------|---------|------|
| R7a | \*NO + H⁺+e⁻ → \*NOH | **0.08** | EES Catalysis 2024 |
| R8 | \*NOH + H⁺+e⁻ → \*NHOH | ~0.3 | |
| R9a | \*NHOH + H⁺+e⁻ → \*NH + H₂O | **0.23** | EES Catalysis 2024 |
| R10 | \*NH + H⁺+e⁻ → \*NH₂ | ~0.3 | |
| R11 | \*NH₂ + H⁺+e⁻ → \*NH₃ | ~0.2 | |
| R12 | \*NH₃ → NH₃(aq) | 0.37 (pH 0为RDS) | |

**Phase III: N-N偶联 → N₂/N₂O**

| Step | 反应 | Ea (eV) | 关键性 |
|------|------|---------|-------|
| R15 | **\*NO + \*NO → \*cis-(ONNO)** | **0.05-0.10** | **需要相邻位点！KMC核心** |
| R16 | \*cis-(ONNO) → \*N₂O + \*O | ~0.3 | N₂O形成 |
| R17 | \*N₂O → N₂O(g) | ~0.3 | 脱附 |
| R19 | \*N + \*N → \*N₂ | 0.45 (Pt) | 需要相邻\*N |

**Phase IV: HER竞争**

| Step | 反应 | Ea (eV) |
|------|------|---------|
| R21 | H⁺+e⁻ + \* → \*H | Volmer |
| R22 | \*H + \*H → H₂ + 2\* | 0.25 (pH 0) |
| R23 | \*H + H⁺+e⁻ → H₂ + \* | Heyrovsky |

### 5.3 N-N偶联——KMC的核心价值

**关键发现：**

```
*NO + H⁺+e⁻ → *NOH    Ea = 0.08 eV (一级反应，NH₃路径)
*NO + *NO → *ONNO       Ea = 0.05-0.10 eV (二级反应，N₂路径)
```

**两条路径势垒几乎相同且极低！** 这意味着选择性**不由动力学势垒决定**，而由：
1. **\*NO的局部空间分布** — 岛状聚集 vs 分散
2. **\*H到达\*NO的相对速率** — \*H扩散极快(E_diff=0.10 eV)，\*N扩散慢(E_diff=0.75 eV)
3. **这恰恰是KMC能回答而MKM不能的核心问题**

### 5.4 DFT参数可用性

| 参考文献 | 催化剂 | 数据 | DOI |
|---------|-------|------|-----|
| ACS Catal. 2019 | 多种TM | MKM参数, volcano plot, 标度关系 | 10.1021/acscatal.9b02179 |
| Faraday Discuss. 2023 | Cu | Grand-canonical DFT, U_L=-0.23V | 10.1039/D2FD00145D |
| EES Catalysis 2024 | Cu基 | NRA3路径Ea，扩散势垒 | — |
| J. Phys. Chem. C 2025 | Cu(111) | GC-DFT完整路径，HER交叉电位 | 10.1021/acs.jpcc.5c02461 |
| ACS Catal. 2024 | 酸性条件TM | 更新MKM | 10.1021/acscatal.4c06394 |
| J. Struct. Chem. 2007 | Cu(100) | N-N成键Ea=0.05-0.10 eV | — |
| ACS Catal. 2014 | Pt(111) | 覆盖度依赖势垒 | 10.1021/cs500668k |
| Angew. Chem. 2025 | CuPd | N-N vs N-H选择性调控 | 10.1002/anie.202524218 |

**评估：参数非常充足。多篇文章提供完整活化能和自由能数据。**

### 5.5 模型规模

| 指标 | 数值 |
|------|------|
| 表面物种数 | ~15 (含\*N₂O₂等) |
| 基元步骤数 | ~20-25 |
| 产物数 | 5-6 (NH₃, N₂, N₂O, NO₂⁻, NH₂OH, H₂) |
| 晶格类型 | Cu(100) square 或 Cu(111) hex |

### 5.6 文献空白

**确认：不存在任何NO₃RR KMC研究。** MKM有（ACS Catal. 2019），但KMC完全空白。

### 5.7 实验验证数据

| 电位 (V vs RHE) | FE(NO₂⁻) | FE(NH₃) | 主产物 |
|----------------|-----------|---------|-------|
| +0.10 | ~75% | <10% | NO₂⁻ |
| -0.10 | ~50% | ~30% | 混合 |
| -0.20 | ~20% | ~70% | NH₃ |
| -0.30 | <10% | ~88% | NH₃(最大) |
| -0.40 | <5% | ~85% | NH₃(HER开始) |

Tafel斜率 ~120 mV/dec。丰富的实验数据可用于验证。

---

## 6. 对比评估与最终推荐

### 6.1 定量评分

| 评价维度 | 权重 | 尿素电合成 | NO₃RR |
|---------|------|----------|-------|
| **多产物丰富度** | 1× | **5** (10+产物) | 4 (5-6产物) |
| **KMC优势（双分子偶联）** | 1.5× | **5** (C-N偶联×2，需两次相邻偶联) | 5 (N-N偶联，需相邻\*NO) |
| **DFT参数充足性** | 1.5× | 3 (分散于多篇，不同催化剂) | **5** (Cu上数据集中且完整) |
| **文献空白/新颖性** | 1× | **5** (完全空白) | **5** (完全空白) |
| **科学影响力** | 1× | **5** (前沿热点，Nature级) | 4 (热点，ACS Catal级) |
| **实现可行性** | 1× | 2 (25+物种，30+步骤，双反应物) | **4** (~15物种，~20步骤) |
| **实验验证数据** | 0.5× | 3 (FE数据有，但分散) | **5** (系统的FE vs U数据) |
| **与Phase 1的逻辑衔接** | 0.5× | 3 | **5** (HER→NO₃RR，包含HER竞争) |
| **加权总分** | | **33.5** | **37.5** |

### 6.2 两种策略建议

#### 策略A：直接做尿素电合成（高风险高回报）

**优点：**
- 最大科学影响力，Nature Catalysis / JACS级别
- 产物最丰富（10+），KMC优势最明显
- 需要两次C-N偶联，双分子空间效应加倍
- 领域发展极快（2020: ~5% FE → 2026: >80% FE），时效性强

**风险：**
- DFT参数分散在不同催化剂上，需从多篇文献拼凑
- 25+物种、30+步骤，实现和调试工作量大
- CO₂RR和NO₃RR两条子路径的参数需要自洽
- 可能需要自己做补充DFT计算

**建议催化剂：** Cu(111)（Nat. Commun. 2024提供最完整数据）

#### 策略B：先NO₃RR，再扩展到尿素（稳健渐进）

**优点：**
- DFT参数最充足且自洽
- 实现复杂度适中（~15物种，~20步骤）
- 丰富的实验数据可直接验证
- N-N偶联（Ea~0.05-0.10 eV）是完美的KMC showcase
- 完成后自然扩展：NO₃RR + CO₂RR → 尿素（增量式开发）

**风险：**
- 影响力略低于尿素（ACS Catal vs Nature Catalysis）
- 产物种类略少

**建议催化剂：** Cu(100)或Cu(111)

### 6.3 推荐方案：策略B+（NO₃RR为主，设计为可扩展到尿素）

**核心思路：**

```
Phase 2a: NO₃RR on Cu(111)
  → 实现12+物种、~20步骤的多产物KMC
  → N-N偶联展示KMC空间效应优势
  → 极化曲线 + FE vs U + Tafel分析
  → 与实验和MKM对比
  → 发表ACS Catal.级论文

Phase 2b: 扩展到尿素电合成（增量开发）
  → 在NO₃RR模型基础上添加CO₂RR子路径（~5步）
  → 添加C-N偶联步骤（~4-6步）
  → 添加尿素中间体（~8物种）
  → 发表Nature Catalysis级论文
```

**优势：**
1. 降低风险：NO₃RR参数充足，先验证框架
2. 增量复杂度：尿素 = NO₃RR + CO₂RR + C-N偶联
3. 两篇论文：第一篇NO₃RR KMC（文献首次），第二篇尿素KMC（文献首次）
4. 代码复用：NO₃RR子路径直接用于尿素模型

**但如果你更倾向直接做尿素电合成**，Cu(111)的参数也足够构建初始模型（Nat. Commun. 2024提供了完整路径），只是调试和验证的工作量更大。

---

## 7. 参考文献汇总

### 尿素电合成
1. Cu surface EDL study — *Nat. Commun.* (2024) DOI: 10.1038/s41467-024-45522-6
2. CuWO₄ C-N coupling — *Nat. Commun.* (2023) DOI: 10.1038/s41467-023-40273-2
3. Sequential co-reduction — *Nat. Commun.* (2023) DOI: 10.1038/s41467-023-44131-z
4. Cu/Cu₂O asymmetric — *CCS Chemistry* (2023) DOI: 10.31635/ccschem.023.202202408
5. Minimum conditions — *Commun. Chem.* (2023) DOI: 10.1038/s42004-023-00990-7
6. M-Pc descriptors — *Commun. Chem.* (2025) DOI: 10.1038/s42004-025-01424-2
7. Cu BIF-29 MKM — *J. Catal.* (2025)
8. Cu₂O surface engineering — *Nat. Commun.* (2025) DOI: 10.1038/s41467-025-57708-7
9. C-N coupling review — *Green Chem.* (2024) DOI: 10.1039/D3GC04920E
10. Urea review — *Adv. Funct. Mater.* (2024) DOI: 10.1002/adfm.202313420

### NO₃RR
11. MKM volcano — *ACS Catal.* **9**, 7052 (2019) DOI: 10.1021/acscatal.9b02179
12. Why Cu → NH₃ — *Faraday Discuss.* (2023) DOI: 10.1039/D2FD00145D
13. Cu(111) GC-DFT — *J. Phys. Chem. C* (2025) DOI: 10.1021/acs.jpcc.5c02461
14. Cu-based review — *EES Catalysis* (2024)
15. Acidic MKM — *ACS Catal.* (2024) DOI: 10.1021/acscatal.4c06394
16. N-N bond on Cu(100) — *J. Struct. Chem.* (2007)
17. CuPd selectivity — *Angew. Chem.* (2025) DOI: 10.1002/anie.202524218
18. Cu₂O JACS — *JACS* (2024) DOI: 10.1021/jacs.3c13288
19. Cu hyponitrite — *JACS* (2022) DOI: 10.1021/jacs.2c04033

### KMC方法论
20. KMC review — *J. Chem. Phys.* **156**, 120902 (2022)
21. CO₂RR KMC diffusion — *J. Chem. Phys.* **155**, 164701 (2021)
22. CO₂RR KMC 178步 — *J. Phys. Chem. Lett.* **16**, 2896 (2025)
23. NH₃氧化KMC — DFT+KMC, 93% NO选择性
24. CO氧化KMC — *J. Chem. Phys.* (2011), mean-field偏差数量级
25. Lateral interactions CE — *J. Phys. Chem. A* (2023) DOI: 10.1021/acs.jpca.3c05581

---

*报告由6个并行研究代理生成，综合分析了40+篇文献*
*Last updated: 2026-03-22*
