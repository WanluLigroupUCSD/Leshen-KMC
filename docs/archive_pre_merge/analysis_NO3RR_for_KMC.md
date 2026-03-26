# NO₃RR 作为 SPARK Phase 2 研究体系的可行性分析

> 2026-03-21

---

## 1. 为什么 NO₃RR 适合空间 KMC 研究

### 1.1 反应复杂度

NO₃RR 是一个 **8电子、9质子** 转移的多步反应，远比 HER (2e⁻) 和 CO₂RR-to-CO (2e⁻) 复杂：

```
NO₃⁻ + 9H⁺ + 8e⁻ → NH₃ + 3H₂O     (主产物)
2NO₃⁻ + 12H⁺ + 10e⁻ → N₂ + 6H₂O   (副产物)
2NO₃⁻ + 10H⁺ + 8e⁻ → N₂O + 5H₂O   (副产物)
```

### 1.2 KMC 的独特价值 — N-N 偶联是双分子反应

**核心论点**：NO₃RR 中 N₂/N₂O 的形成需要 \*NO + \*NO 偶联（类似 CO₂RR 中 \*CO + \*CO → C₂），这是一个**需要相邻位点的双分子表面反应**。Mean-field 模型假设 \*NO 均匀分布：

$$r_{N-N} = k_{NN} \cdot \theta_{NO}^2 \quad (\text{mean-field})$$

但实际上 \*NO 的空间分布取决于：
- **横向相互作用**（吸引→岛状聚集→增强N-N偶联→更多N₂）
- **表面扩散速率**（慢扩散→空间不均匀→selectivity变化）
- **\*H 竞争覆盖**（H\* 占据位点→稀释 \*NO→抑制N-N偶联→更多NH₃）

**这与 CO₂RR 中 C₂ 选择性问题完全类似，但 NO₃RR 更复杂**：
- CO₂RR: \*CO + \*CO → C₂ vs \*CO + H → C₁
- NO₃RR: \*NO + \*NO → N₂/N₂O vs \*NO + H → NH₃ pathway

### 1.3 与 CO₂RR 的对比优势

| 特性 | CO₂RR (Cu) | NO₃RR (Cu) | KMC 优势 |
|------|-----------|-----------|----------|
| 电子转移数 | 2-12 | 8 (NH₃), 10 (N₂) | NO₃RR 路径更多 |
| 中间物种数 | ~8 | **~10+** | 更丰富的表面竞争 |
| 双分子关键步 | CO\*-CO\* 偶联 | **NO\*-NO\* 偶联** | ✅ 空间效应关键 |
| H\* 竞争 | 弱 (CO₂RR电位下HER弱) | **强** (NO₃RR与HER电位重叠) | ✅ 覆盖度竞争 |
| 产物选择性问题 | C₁ vs C₂ | **NH₃ vs N₂ vs N₂O** | ✅ 三路竞争 |
| 晶面效应 | Cu(111) vs Cu(100) | **Cu(100) vs Cu(111)** | 类似 |
| KMC 文献 | 有 (Wei 2025, Jørgensen 2021) | **空白** | ✅ 明确文献空白 |
| MKM 文献 | 有 | 有 (*ACS Catal.* **9**, 7052 (2019)) | 可做KMC vs MKM对比 |

---

## 2. NO₃RR 反应网络详解

### 2.1 主要路径（NH₃ 路径）

文献共识的最优路径 (*ACS Catal.* (2021); *Faraday Discuss.* (2023); *J. Phys. Chem. C* (2025)):

```
NO₃⁻(aq)                           [溶液相]
    │ adsorption
    ▼
*NO₃ + H⁺ + e⁻ → *NO₂ + *OH       [PCET, 脱氧]  Step 1
    │
*NO₂ + H⁺ + e⁻ → *NO + *OH        [PCET, 脱氧]  Step 2
    │
    ▼ ──── 关键分支点 (*NO) ────
    │                              │
    │ [NH₃路径]                     │ [N₂路径]
    │                              │
*NO + H⁺ + e⁻ → *NOH   [PCET]    *NO + *NO → *N₂O  [Thermal, 双分子!]
    │                              │
*NOH + H⁺ + e⁻ → *NHOH [PCET]    *N₂O → N₂O(g)  或  *N₂O + H⁺ + e⁻ → N₂
    │
*NHOH + H⁺ + e⁻ → *NH + H₂O [PCET, 脱氧]
    │
*NH + H⁺ + e⁻ → *NH₂   [PCET]
    │
*NH₂ + H⁺ + e⁻ → *NH₃  [PCET]
    │
*NH₃ → NH₃(aq)          [Desorption]
```

### 2.2 物种清单

| # | Species | Description |
|---|---------|-------------|
| 1 | \* (empty) | 空位 |
| 2 | \*NO₃ | 吸附硝酸根 |
| 3 | \*NO₂ | 吸附亚硝酸根 |
| 4 | \*NO | **关键分支中间体** |
| 5 | \*NOH | 氮氧氢化物 |
| 6 | \*NHOH | 羟胺中间体 |
| 7 | \*NH | 亚胺 |
| 8 | \*NH₂ | 氨基 |
| 9 | \*NH₃ | 吸附氨 |
| 10 | \*N₂O | 一氧化二氮（N-N偶联产物） |
| 11 | \*OH | 吸附羟基（副产物，需脱附） |
| 12 | \*H | 吸附氢（HER竞争） |

**共12种表面物种**，比 CO₂RR 的8种更多。

### 2.3 基元步骤（~15步）

| # | Reaction | Type | E_a [eV] | Note |
|---|----------|------|----------|------|
| 1 | NO₃⁻(aq) + \* → \*NO₃ | Adsorption | ~0 | 硝酸根吸附 |
| 2 | \*NO₃ + H⁺ + e⁻ → \*NO₂ + \*OH | PCET | ~0.3–0.5 | 第一步脱氧 |
| 3 | \*NO₂ + H⁺ + e⁻ → \*NO + \*OH | PCET | ~0.3–0.5 | 第二步脱氧 |
| 4 | **\*NO + H⁺ + e⁻ → \*NOH** | **PCET** | **~0.6–0.8** | **PLS (potential limiting step)**, NH₃路径入口 |
| 5 | \*NOH + H⁺ + e⁻ → \*NHOH | PCET | ~0.3 | 加氢 |
| 6 | \*NHOH + H⁺ + e⁻ → \*NH + H₂O | PCET | ~0.2–0.4 | 脱氧+加氢 |
| 7 | \*NH + H⁺ + e⁻ → \*NH₂ | PCET | ~0.3 | 加氢 |
| 8 | \*NH₂ + H⁺ + e⁻ → \*NH₃ | PCET | ~0.2 | 加氢 |
| 9 | \*NH₃ → NH₃(aq) + \* | Desorption | ~0.37 | 氨脱附 |
| 10 | **\*NO + \*NO(adj) → \*N₂O + \*** | **Thermal** | **~0.5–0.8** | **N-N偶联！双分子！** |
| 11 | \*N₂O → N₂O(g) + \* | Desorption | ~0.3 | N₂O脱附 |
| 12 | \*OH + H⁺ + e⁻ → H₂O + \* | PCET | ~0.1 | OH清除 |
| 13 | H₂O + \* + e⁻ → \*H + OH⁻ | PCET | ~0.7 | Volmer (HER竞争) |
| 14 | \*H + \*H → H₂ + 2\* | Thermal | ~0.8 | Tafel (HER) |
| 15 | \*NO diffusion | Thermal | ~0.1 | 表面扩散 (可调) |

**注意：** 以上活化能为文献估计值范围，具体数值需从DFT文献提取。关键参考：
- *ACS Catal.* (2021) — Cu上NO₃RR完整路径DFT
- *ACS Catal.* **9**, 7052–7064 (2019) — 过渡金属NO₃RR MKM
- *ACS Catal.* (2024) — 酸性条件NO₃RR MKM
- *J. Phys. Chem. C* (2025) — Cu(111)上NO₃RR新见解
- *Faraday Discuss.* (2023) — 为什么Cu催化NO₃RR到NH₃

### 2.4 关键分支点：\*NO

```
                    *NO
                   ╱    ╲
          [加氢路径]      [偶联路径]
         *NO + H⁺+e⁻     *NO + *NO(adj)
              │                 │
           *NOH              *N₂O
              │                 │
           ...→ NH₃          N₂O / N₂
```

**\*NO 是 NH₃ vs N₂ 选择性的关键决定者**：
- 高 \*NO 覆盖度 + 聚集分布 → 更多 N-N 偶联 → N₂/N₂O
- 低 \*NO 覆盖度 / 分散分布 → 加氢为主 → NH₃
- 高 \*H 覆盖度 → 稀释 \*NO + 加速加氢 → NH₃

**这正是空间 KMC 能回答而 mean-field 不能的核心问题。**

---

## 3. 横向相互作用与空间效应分析

### 3.1 预期横向相互作用

| Pair | Expected Interaction | Physical Origin | Effect on Selectivity |
|------|---------------------|-----------------|----------------------|
| \*NO–\*NO | Variable (need DFT) | Dipole-dipole | 吸引→聚集→促进N₂ |
| \*NO–\*H | Repulsive (+) | Steric/electronic | 排斥→相分离→影响选择性 |
| \*OH–\*NO | Repulsive (+) | Steric | OH累积→毒化表面 |
| \*NO–\*NH | Weak | — | 次要 |

### 3.2 空间 KMC 独有的可观测量

1. **\*NO pair correlation function g_NO-NO(r)**：反映 NO 聚集/分散程度
2. **\*NO cluster size distribution**：岛状结构直接决定 N-N 偶联概率
3. **\*NO–\*H 空间分离**：H-rich区域 vs NO-rich区域的空间分相
4. **N₂O/NH₃ selectivity vs \*NO 空间分布**：空间效应对选择性的定量影响

### 3.3 Mean-field 失效的预期场景

| 场景 | MKM 预测 | KMC 预测 | 原因 |
|------|---------|---------|------|
| \*NO-\*NO 吸引 | r_NN ∝ θ_NO² | r_NN > MKM预测 | 岛状聚集增强相邻概率 |
| \*NO-\*H 排斥 | 均匀混合 | 空间分相 | \*H和\*NO形成不同域 |
| 低 \*NO 覆盖度 | θ_NO²≈0 → N₂≈0 | θ_NO涨落→仍有局部高密度 | 随机涨落 |
| 高 \*OH 覆盖度 | 均匀毒化 | 局部毒化+活性区 | 空间不均匀 |

---

## 4. 可行的 KMC 研究设计

### 4.1 核心科学问题

> **\*NO 吸附物在 Cu 表面上的空间分布（由横向相互作用、扩散和 \*H 竞争控制）是否决定了电化学 NO₃RR 中 NH₃ 与 N₂ 的选择性？**

### 4.2 建议的研究方案

#### Study A: \*NO 空间分布 vs 横向相互作用

扫描 ε(\*NO–\*NO) 从 −0.05 到 +0.10 eV，观察：
- \*NO 聚集程度 (cluster size, g(r))
- N-N 偶联速率变化
- NH₃/N₂ selectivity 变化

#### Study B: \*H 竞争对选择性的影响

扫描 Volmer 势垒（控制 \*H 产生速率）或电位，观察：
- \*H 和 \*NO 的空间共存/分离模式
- KMC vs MKM 的 NH₃/N₂ selectivity 偏差

#### Study C: 扩散速率效应

扫描 \*NO 扩散势垒 (0.05–0.50 eV)：
- 快扩散：\*NO均匀分布，KMC≈MKM
- 慢扩散：空间不均匀，KMC≠MKM
- 确定 mean-field breakdown 阈值

#### Study D: 极化曲线与法拉第效率

全电位扫描，计算：
- j_total, FE(NH₃), FE(N₂), FE(N₂O), FE(H₂)
- KMC vs MKM 对比
- 与实验比较

### 4.3 模拟参数

| Parameter | Value |
|-----------|-------|
| Lattice | 50×50 2D square (Cu(100) top sites) |
| Species | 12 types |
| Elementary steps | ~15 |
| Temperature | 298 K |
| Potential range | −0.2 to −1.0 V vs RHE |
| Equilibration | 10⁶ steps |
| Production | 5×10⁶ steps |

---

## 5. 文献基础与空白分析

### 5.1 已有 DFT 参数来源

| Reference | System | Data Available |
|-----------|--------|----------------|
| *ACS Catal.* (2021) "Superior Nitrate Reduction to Ammonia Performance of Cu Catalysts" | Cu(111)/(100)/(110) | 完整自由能图，所有中间体 |
| *ACS Catal.* **9**, 7052–7064 (2019) "Activity and Selectivity Trends in NO₃RR on Transition Metals" | 多种过渡金属 | Volcano plot, MKM参数, DRC分析 |
| *ACS Catal.* (2024) "Nitrate Reduction Modeling under Acidic Conditions with Late Transition Metals" | Pt, Cu, etc. | 酸性MKM, 电位依赖 |
| *Faraday Discuss.* (2023) "Why copper catalyzes electrochemical reduction of nitrate to ammonia" | Cu | Grand-canonical DFT, 自由能 |
| *JACS* (2024) "Electrocatalytic NO₃RR and NO₂RR Using Cu₂O" | Cu₂O(100)/(111) | Cu(I)/Cu(0)活性物种，机理 |
| *J. Phys. Chem. C* (2025) "New Insights into NO₃RR on Cu(111)" | Cu(111) | 最新DFT自由能 |

### 5.2 文献空白

| 方向 | 现有工作 | 空白 |
|------|---------|------|
| NO₃RR MKM | ✅ *ACS Catal.* **9**, 7052 (2019) | MKM做了，但mean-field |
| NO₃RR KMC | ❌ **完全空白** | **没有任何KMC研究** |
| N-N偶联空间效应 | ❌ 未研究 | **mean-field无法处理** |
| \*NO-\*H空间竞争 | ❌ 未研究 | **KMC独有能力** |
| 覆盖度对NH₃/N₂选择性 | 实验观察 → 机理不明 | **KMC可揭示机理** |

**结论：NO₃RR的KMC研究是一个完全的文献空白，且存在明确的科学问题。**

---

## 6. 与 CO₂RR 方案的对比评估

| 评价维度 | CO₂RR (原Phase 2) | NO₃RR (新方案) | 评估 |
|---------|-------------------|---------------|------|
| 反应复杂度 | 中等 (~12步, 8物种) | **高** (~15步, 12物种) | NO₃RR更符合"复杂反应" |
| 双分子关键步 | CO\*+CO\* (C₂) | \*NO+\*NO (N₂) | 同等适合KMC |
| 覆盖度竞争 | CO\*为主 | **\*NO vs \*H vs \*OH 三方竞争** | NO₃RR空间效应更丰富 |
| 选择性维度 | C₁ vs C₂ (2路) | **NH₃ vs N₂ vs N₂O (3路)** | NO₃RR更复杂 |
| DFT参数可得性 | 充分 (Wei 2025等) | 充分 (多篇DFT文献) | 相当 |
| KMC文献空白 | 部分空白 | **完全空白** | NO₃RR新颖性更强 |
| 应用价值 | CO₂利用 (高影响力) | 废水处理+绿色NH₃合成 | 相当 |
| 实现难度 | 中等 | **较高** (物种多) | NO₃RR需更多编程 |

---

## 7. 推荐方案

### 7.1 总体判断

**NO₃RR 是一个比 CO₂RR 更具挑战性但也更有新颖性的 KMC 研究体系**：
- **文献完全空白**：没有任何NO₃RR-KMC研究
- **科学问题清晰**：\*NO空间分布如何决定NH₃ vs N₂选择性
- **与Phase 1逻辑一致**：HER验证→复杂电催化体系（NO₃RR包含HER竞争）

### 7.2 下一步建议

1. **T5**: 从文献提取Cu(100)上NO₃RR完整DFT参数（活化能、吸附能）
2. **T6**: 实现NO₃RR模型（12物种、~15步）
3. **Study A-D**: 按§4.2执行四项研究

### 7.3 关键参考文献

1. *ACS Catal.* **9**, 7052–7064 (2019) — NO₃RR过渡金属MKM+volcano [DOI](https://pubs.acs.org/doi/abs/10.1021/acscatal.9b02179)
2. *ACS Catal.* (2021) "Superior Nitrate Reduction to NH₃ Performance of Cu" [DOI](https://pubs.acs.org/doi/abs/10.1021/acscatal.1c03666)
3. *ACS Catal.* (2024) "Nitrate Reduction Modeling under Acidic Conditions" [DOI](https://pubs.acs.org/doi/10.1021/acscatal.4c06394)
4. *Faraday Discuss.* (2023) "Why copper catalyzes NO₃RR to NH₃"
5. *JACS* (2024) "Electrocatalytic NO₃RR and NO₂RR Using Cu₂O" [DOI](https://pubs.acs.org/doi/10.1021/jacs.3c13288)
6. *J. Phys. Chem. C* (2025) "New Insights into NO₃RR on Cu(111)" [DOI](https://pubs.acs.org/doi/10.1021/acs.jpcc.5c02461)
7. *Chem Catal.* (2023) "NO₃RR Selectivity at the crossroads between NH₃ and N₂"
8. *Commun. Chem.* (2025) "Recent advances in mechanistic studies for NO₃RR to NH₃" [DOI](https://www.nature.com/articles/s42004-025-01864-w)

---

*Last updated: 2026-03-21*
