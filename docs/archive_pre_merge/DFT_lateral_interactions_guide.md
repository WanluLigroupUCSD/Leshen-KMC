# NO₃RR on Cu：DFT计算横向相互作用完整指南

> SPARK Phase 2 | 2026-03-22
> 所有能量参数从真空DFT计算获取，不使用文献数值

---

## 目录

1. [核心公式](#1-核心公式)
2. [Cu表面slab模型设置](#2-cu表面slab模型设置)
3. [各物种吸附位点](#3-各物种吸附位点确定)
4. [横向相互作用计算方案](#4-横向相互作用计算方案)
5. [完整DFT计算清单](#5-完整dft计算清单)
6. [VASP计算参数](#6-vasp计算参数)
7. [数据提取与KMC参数映射](#7-数据提取与kmc参数映射)
8. [预期结果与验证](#8-预期结果与验证)
9. [参考文献](#9-参考文献)

---

## 1. 核心公式

### 1.1 成对相互作用能（pairwise interaction energy）

对于吸附物A和B在表面上相距d时的相互作用能：

```
ε(A,B,d) = E(slab+A+B@d) − E(slab+A) − E(slab+B) + E(clean_slab)
```

- `ε > 0`：**排斥**（repulsive）— 吸附物不稳定化，加速反应/脱附
- `ε < 0`：**吸引**（attractive）— 吸附物稳定化，减慢脱附

需要 **4个DFT计算** 确定一对吸附物在一个距离的ε值。

**物理验证**：Curulla Ferré et al. (ChemPhysChem 2005, DOI: 10.1002/cphc.200400399) 证明成对可加性近似(pairwise additive approximation)对于表面吸附物-吸附物相互作用是很好的近似，三体项通常可忽略。

### 1.2 等价表述——用结合能差

```
ε(A,B,d) = BE(A, 有B在距离d处) − BE(A, 孤立)
```

其中：
- `BE(A, 有B) = E(slab+A+B) − E(slab+B) − E(A_gas)`
- `BE(A, 孤立) = E(slab+A) − E(slab) − E(A_gas)`

气相能量 E(A_gas) 相消，结果与公式1.1相同。

### 1.3 如果A和B是同种物种（对称情况）

当两个位点等价时 E(slab+A@site1) = E(slab+A@site2)，简化为：

```
ε(A,A,d) = E(co-ads) − 2 × E(single) + E(clean)
```

### 1.4 在KMC中的使用方式

```
k(site) = k_base × exp(+E_lat / kBT)
```

其中 `E_lat = Σ_j ε(sp_i, sp_j)` 对所有被占据的最近邻j求和。

**量级参考**：ε = 0.1 eV → 室温(298K)下速率改变 ~50倍；ε = 0.2 eV → ~2400倍。

---

## 2. Cu表面slab模型设置

### 2.1 表面选择

推荐 **Cu(111)**，理由：
- NO₃RR DFT文献最多的面，便于交叉验证
- fcc hollow / hcp hollow / bridge / atop 四种高对称位点
- 最近邻距离 = a_Cu/√2 = 3.615/√2 = **2.556 Å**

备选 Cu(100)，有4-fold hollow位点，square lattice与KMC更自然对应。

### 2.2 slab参数

| 参数 | 推荐值 | 说明 |
|------|-------|------|
| 层数 | **4层** (底2固定，顶2弛豫) | 最低要求；5层(底3固定，顶2弛豫)更好 |
| 真空层 | **15 Å** | 防止周期性镜像相互作用 |
| 晶格常数 | 先优化bulk Cu (PBE: ~3.63 Å, 实验: 3.615 Å) | 用优化值构建slab |
| 偶极修正 | **必须开启** (LDIPOL=.TRUE., IDIPOL=3) | 非对称slab必须 |

### 2.3 超胞大小

| 超胞 | 原子数(4层) | 单吸附物覆盖度 | 可计算距离 | k-points | 推荐 |
|------|-----------|-------------|---------|----------|------|
| p(2×2) | 16 Cu | 0.25 ML | 仅1NN(镜像干扰大) | 5×5×1 | ❌ 不推荐 |
| **p(3×3)** | **36 Cu** | **0.11 ML** | **1NN + 2NN** | **3×3×1** | **✅ 最低推荐** |
| **p(4×4)** | **64 Cu** | **0.0625 ML** | **1NN + 2NN + 3NN** | **2×2×1** | **✅ 最佳平衡** |
| p(5×5) | 100 Cu | 0.04 ML | 到4NN-5NN | 2×2×1 | 用于精确长程 |

**推荐：p(3×3) 用于初步计算（快），p(4×4) 用于最终精确值。**

**镜像距离检查**：p(3×3)晶胞边长 = 3 × 2.556 = 7.67 Å。1NN对(2.556 Å)到最近镜像的距离 = 7.67 − 2.556 = 5.11 Å，可接受。

### 2.4 Cu(111)上的邻居距离

| 邻居壳层 | 距离 (Å) | 关系 | 重要性 |
|---------|---------|------|-------|
| **1NN** | **2.556** | a/√2 | **必须计算** |
| **2NN** | **4.427** | a√(3/2) | **建议计算** |
| 3NN | 5.112 | a√2 | 可选 |
| 4NN+ | >5.9 | — | 通常可忽略 |

**衰减规律**：横向相互作用大致按 ~1/r³ (偶极-偶极) 或 ~1/r⁵ (基底介导) 衰减。2NN通常为1NN的20-40%，3NN为5-15%。

---

## 3. 各物种吸附位点确定

### 3.1 第一步：位点偏好测试

对每种物种，需在 **4个高对称位点**（fcc hollow, hcp hollow, bridge, atop）分别优化，找到能量最低的位点。

**Cu(111)上各物种的文献预期位点：**

| 物种 | 预期最优位点 | 结合能参考(eV) | 吸附构型 | 自旋 |
|------|-----------|-------------|---------|------|
| \*H | fcc hollow | −2.45 | H朝下 | 无 |
| \*O | fcc hollow | 强(< −4) | O朝下 | 无 |
| \*N | fcc hollow | 强(< −3) | N朝下 | 无 |
| \*OH | fcc hollow 或 bridge | −3.0 ~ −3.5 | O朝下，H倾斜 | 无 |
| \*NO | **fcc hollow** | −1.0 ~ −1.1 | **N朝下，分子垂直** | **有！** |
| \*NOH | atop 或 bridge | — | N朝下 | 需测试 |
| \*NH | fcc hollow | 强 | N朝下 | 需测试 |
| \*NH₂ | bridge | 中等 | N在bridge | 无 |
| \*NH₃ | atop | −0.37 (弱) | N在atop | 无 |
| \*NO₂ | bidentate (O,O) | — | 两个O桥接两个Cu | 需测试 |
| \*NO₃ | bidentate/tridentate | — | O原子配位Cu | 需测试 |
| \*NHOH | bridge 或 hollow | — | N朝下 | 需测试 |
| \*N₂O | atop 或 bridge | 弱 | 可能立即脱附 | 无 |

**⚠️ 重要**：\*NO是自由基（有未配对电子），DFT计算 **必须开自旋极化 (ISPIN=2)**！

### 3.2 位点测试的DFT计算量

| 物种数 | 位点数/物种 | 总计算数 |
|-------|-----------|---------|
| 10-12种 | 4 (fcc/hcp/bridge/atop) | **~40-48个结构优化** |

可以用较小超胞p(2×2)做初步位点筛选，然后用p(3×3)或p(4×4)做精确计算。

---

## 4. 横向相互作用计算方案

### 4.1 需要计算的吸附物对——按优先级

#### Tier 1：必须计算（直接决定NO₃RR选择性）

| 编号 | 吸附物对 | 重要性 | 原因 |
|------|---------|-------|------|
| **P1** | **\*NO + \*NO** | ⭐⭐⭐⭐⭐ | **N-N偶联选择性的核心参数！排斥→分散→抑制N₂；吸引→聚集→促进N₂** |
| **P2** | **\*NO + \*H** | ⭐⭐⭐⭐⭐ | **加氢路径(NH₃) vs 偶联路径(N₂)的空间竞争** |
| **P3** | **\*H + \*H** | ⭐⭐⭐⭐ | HER竞争，\*H有序化 |
| **P4** | **\*OH + \*NO** | ⭐⭐⭐⭐ | OH累积对\*NO的影响 |

#### Tier 2：强烈建议（影响覆盖度分布和动力学）

| 编号 | 吸附物对 | 重要性 | 原因 |
|------|---------|-------|------|
| **P5** | \*OH + \*OH | ⭐⭐⭐ | OH累积毒化效应 |
| **P6** | \*OH + \*H | ⭐⭐⭐ | OH清除(Volmer-Heyrovsky)动力学 |
| **P7** | \*N + \*NO | ⭐⭐⭐ | 替代N-N偶联路径(\*N+\*NO→\*N₂O) |
| **P8** | \*N + \*N | ⭐⭐⭐ | 直接N₂形成 |
| **P9** | \*O + \*NO | ⭐⭐ | 脱氧步骤副产物\*O的影响 |
| **P10** | \*N + \*H | ⭐⭐ | \*N加氢动力学 |

#### Tier 3：精细化（如有计算资源）

| 编号 | 吸附物对 | 重要性 |
|------|---------|-------|
| P11 | \*NH + \*NO | ⭐ |
| P12 | \*NH₂ + \*H | ⭐ |
| P13 | \*O + \*O | ⭐ |
| P14 | \*O + \*H | ⭐ |
| P15 | \*OH + \*O | ⭐ |

### 4.2 每对的DFT构型设置

以 **\*NO + \*NO 在 Cu(111) p(3×3)** 为例：

**fcc hollow位点编号 (俯视图，×标记fcc位点)：**

```
第3行:  ×7   ×8   ×9
第2行:  ×4   ×5   ×6
第1行:  ×1   ×2   ×3
```

| 构型 | 位点 | 距离 | 目的 |
|------|------|------|------|
| 共吸附-1NN | NO@1 + NO@2 | 2.556 Å | ε(NO-NO, 1NN) |
| 共吸附-2NN | NO@1 + NO@5 | 4.427 Å | ε(NO-NO, 2NN) |
| 单吸附 | NO@1 only | 孤立 | E(slab+NO) |
| 清洁表面 | 无吸附物 | — | E(clean) |

**⚠️ 重要注意事项：**
- 确保两个\*NO的初始构型合理（都是N朝下，垂直于表面）
- 弛豫后检查结构是否稳定（没有NO脱附或解离）
- 检查磁矩变化（共吸附时自旋态可能变化）

### 4.3 异种吸附物对的设置

以 **\*NO + \*H** 为例：

| 构型 | 位点 | 目的 |
|------|------|------|
| 共吸附-1NN | NO@fcc1 + H@fcc2 | ε(NO-H, 1NN) |
| 共吸附-2NN | NO@fcc1 + H@fcc5 | ε(NO-H, 2NN) |
| 单\*NO | NO@fcc1 | E(slab+NO) |
| 单\*H | H@fcc1 | E(slab+H) |
| 清洁 | — | E(clean) |

注意：E(slab+NO) 和 E(clean) 已在前面计算过，**可复用**。

---

## 5. 完整DFT计算清单

### 5.1 Phase 0：基础设置

| 编号 | 计算 | 超胞 | 目的 | 数量 |
|------|------|------|------|------|
| B1 | bulk Cu优化 | 常规 | 获取晶格常数 | 1 |
| B2 | 气相分子优化 | 大box(15×16×17 Å) | NO, H₂, H₂O, NH₃, N₂, N₂O, NO₂, OH | ~8 |
| | | | **小计** | **9** |

### 5.2 Phase 1：位点偏好测试

| 编号 | 计算 | 超胞 | 目的 | 数量 |
|------|------|------|------|------|
| S0 | 清洁slab | p(3×3) | E(clean) | 1 |
| S1-S10 | 各物种×4位点 | p(3×3) | 找最优吸附位点 | ~40 |
| | | | **小计** | **~41** |

### 5.3 Phase 2：精确单吸附物（在最优位点）

| 编号 | 计算 | 超胞 | 目的 | 数量 |
|------|------|------|------|------|
| R0 | 清洁slab | p(4×4) | E(clean) 精确 | 1 |
| R1-R7 | 各关键物种×1 | p(4×4) | E(slab+A) 精确 | 7 |
| | 物种：\*NO, \*H, \*OH, \*O, \*N, \*NH, \*NH₂ | | | |
| | | | **小计** | **8** |

### 5.4 Phase 3：Tier 1 成对相互作用

| 编号 | 吸附物对 | 距离 | 超胞 | 数量 |
|------|---------|------|------|------|
| T1a | \*NO + \*NO | 1NN | p(4×4) | 1 |
| T1b | \*NO + \*NO | 2NN | p(4×4) | 1 |
| T2a | \*NO + \*H | 1NN | p(4×4) | 1 |
| T2b | \*NO + \*H | 2NN | p(4×4) | 1 |
| T3a | \*H + \*H | 1NN | p(4×4) | 1 |
| T3b | \*H + \*H | 2NN | p(4×4) | 1 |
| T4a | \*OH + \*NO | 1NN | p(4×4) | 1 |
| T4b | \*OH + \*NO | 2NN | p(4×4) | 1 |
| | | | **小计** | **8** |

### 5.5 Phase 4：Tier 2 成对相互作用

| 编号 | 吸附物对 | 距离 | 超胞 | 数量 |
|------|---------|------|------|------|
| T5a-b | \*OH + \*OH | 1NN, 2NN | p(4×4) | 2 |
| T6a-b | \*OH + \*H | 1NN, 2NN | p(4×4) | 2 |
| T7a-b | \*N + \*NO | 1NN, 2NN | p(4×4) | 2 |
| T8a-b | \*N + \*N | 1NN, 2NN | p(4×4) | 2 |
| T9a-b | \*O + \*NO | 1NN, 2NN | p(4×4) | 2 |
| T10a-b | \*N + \*H | 1NN, 2NN | p(4×4) | 2 |
| | | | **小计** | **12** |

### 5.6 Phase 5（可选）：Tier 3 + 三体项

| 编号 | 计算 | 数量 |
|------|------|------|
| T11-T15 | Tier 3 各对 1NN | 5 |
| 3B1 | \*NO+\*NO+\*NO 紧凑三角 | 1 |
| 3B2 | \*NO+\*NO+\*H 紧凑三角 | 1 |
| | **小计** | **7** |

### 5.7 总计算量汇总

| Phase | 内容 | 计算数 | 超胞 | 预估时间/个 |
|-------|------|-------|------|-----------|
| Phase 0 | 基础设置 | 9 | 小 | 0.5 h |
| Phase 1 | 位点测试 | 41 | p(3×3) 36 Cu | 1-2 h |
| Phase 2 | 精确单吸附 | 8 | p(4×4) 64 Cu | 3-5 h |
| Phase 3 | Tier 1 成对 | 8 | p(4×4) 64 Cu | 3-5 h |
| Phase 4 | Tier 2 成对 | 12 | p(4×4) 64 Cu | 3-5 h |
| Phase 5 | Tier 3 + 三体 | 7 | p(4×4) 64 Cu | 3-5 h |
| **总计** | | **~85** | | |

**预估总机时**：
- Phase 0-1（小超胞）：~41×1.5h = **~62 CPU小时**
- Phase 2-5（大超胞）：~35×4h = **~140 CPU小时**
- **总计约200 CPU小时**（单核计算时间；并行计算可大幅缩短）
- 如用64核节点：**~3-5天可完成全部计算**

---

## 6. VASP计算参数

### 6.1 INCAR（结构优化）

```fortran
! === 基本设置 ===
SYSTEM = Cu111_NO3RR_lateral

! === 泛函 ===
GGA    = PE          ! PBE泛函（推荐：相对能量准确）
                     ! 或 GGA = RP 用RPBE（吸附能绝对值更准）
                     ! 注意：横向相互作用是能量差，PBE和RPBE结果接近

! === 精度 ===
PREC   = Accurate
ENCUT  = 450         ! 截断能 (eV)，需 ≥ POTCAR中ENMAX×1.3
EDIFF  = 1E-5        ! 电子收敛 (eV)
EDIFFG = -0.02       ! 力收敛 (eV/Å)，负值=力判据
LREAL  = Auto        ! 大超胞用Auto加速

! === 离子弛豫 ===
NSW    = 300         ! 最大离子步数
IBRION = 2           ! 共轭梯度法
ISIF   = 2           ! 仅弛豫原子，固定晶胞

! === 电子设置 ===
ISMEAR = 1           ! Methfessel-Paxton展宽（金属）
SIGMA  = 0.15        ! 展宽宽度 (eV)
ALGO   = Fast        ! 电子优化算法
NELM   = 200         ! 最大电子步数

! === 自旋 ===
ISPIN  = 2           ! ⚠️ 必须！NO是自由基
                     ! 即使Cu无磁性，含NO时必须开自旋

! === 偶极修正 ===
LDIPOL = .TRUE.      ! ⚠️ 必须！非对称slab
IDIPOL = 3           ! 沿z方向(表面法线)修正

! === 输出 ===
LORBIT = 11          ! 输出PDOS（可选）
LWAVE  = .FALSE.     ! 不保存波函数（节省磁盘）
LCHARG = .FALSE.     ! 不保存电荷密度（节省磁盘）
```

### 6.2 KPOINTS

**p(3×3) 超胞：**
```
K-Points
0
Gamma
3 3 1
0 0 0
```

**p(4×4) 超胞：**
```
K-Points
0
Gamma
2 2 1
0 0 0
```

### 6.3 POSCAR设置要点

```
Cu(111) p(3x3) slab, 4 layers
1.0
  7.6690  0.0000  0.0000     ! 3 × a_Cu/√2
  3.8345  6.6407  0.0000     ! 120°六角
  0.0000  0.0000  25.0000    ! slab + 15 Å真空
Cu  N  O  H                  ! 元素（按需添加）
36  ...                      ! 原子数
Selective Dynamics           ! ⚠️ 必须！固定底层
Cartesian                    ! 或 Direct
... (坐标)
```

**Selective Dynamics标记：**
- 底2层Cu原子：`F F F`（固定）
- 顶2层Cu原子：`T T T`（弛豫）
- 所有吸附物原子：`T T T`（弛豫）

### 6.4 气相分子计算设置

```fortran
! 气相分子用大box
ISPIN  = 2           ! NO有自旋
ISMEAR = 0           ! Gaussian展宽（分子/绝缘体）
SIGMA  = 0.01
IBRION = 2
NSW    = 100
EDIFFG = -0.01
! box大小：15 × 16 × 17 Å (三个方向不同，避免对称性问题)
! k-points: Gamma only (1×1×1)
```

### 6.5 推荐POTCAR

```
Cu: PAW_PBE Cu 22Jun2005     (ENMAX = 295 eV)
N:  PAW_PBE N 08Apr2002      (ENMAX = 400 eV)
O:  PAW_PBE O 08Apr2002      (ENMAX = 400 eV)
H:  PAW_PBE H 15Jun2001      (ENMAX = 250 eV)
```

ENCUT = 450 eV ≥ max(ENMAX) × 1.13，满足精度要求。

---

## 7. 数据提取与KMC参数映射

### 7.1 从VASP输出提取能量

从每个OUTCAR文件提取最终总能量（TOTEN或E0）：

```bash
grep "energy  without entropy" OUTCAR | tail -1
# 输出: energy  without entropy=  -XXX.XXXXXXXX  energy(sigma->0) =  -XXX.XXXXXXXX
# 使用 "energy(sigma->0)" 的值 (即 E0 = F + TS → S→0 外推)
```

### 7.2 计算示例

假设得到以下DFT总能量（虚拟数值，单位eV）：

```
E(clean_slab)         = -180.000
E(slab + NO@fcc1)     = -186.500
E(slab + NO@fcc2)     = -186.500   (等价位点，应相同)
E(slab + NO@1 + NO@2) = -192.850   (1NN共吸附)
E(slab + NO@1 + NO@5) = -192.950   (2NN共吸附)
```

**计算ε：**

```
ε(NO-NO, 1NN) = -192.850 - (-186.500) - (-186.500) + (-180.000)
              = -192.850 + 186.500 + 186.500 - 180.000
              = +0.150 eV (排斥)

ε(NO-NO, 2NN) = -192.950 + 186.500 + 186.500 - 180.000
              = +0.050 eV (弱排斥)
```

### 7.3 结果整理表格模板

| 吸附物对 | ε_1NN (eV) | ε_2NN (eV) | 性质 | 影响 |
|---------|-----------|-----------|------|------|
| \*NO–\*NO | +0.150 | +0.050 | 排斥 | 分散\*NO，抑制N-N偶联 |
| \*NO–\*H | +0.060 | +0.015 | 弱排斥 | \*NO和\*H相分离 |
| \*H–\*H | +0.045 | +0.010 | 弱排斥 | \*H有序分布 |
| \*OH–\*NO | +0.120 | +0.035 | 排斥 | OH毒化\*NO位点 |
| ... | ... | ... | ... | ... |

### 7.4 映射到SPARK

```python
# 在模型定义中直接使用DFT计算的ε值
model.add_lateral_interaction("NO", "NO", epsilon_1NN=0.150, epsilon_2NN=0.050)
model.add_lateral_interaction("NO", "H",  epsilon_1NN=0.060, epsilon_2NN=0.015)
model.add_lateral_interaction("H",  "H",  epsilon_1NN=0.045, epsilon_2NN=0.010)
model.add_lateral_interaction("OH", "NO", epsilon_1NN=0.120, epsilon_2NN=0.035)
# ...

# KMC引擎自动计算每个位点的速率修正：
# E_lat = Σ_j ε(sp_i, sp_j, d_ij)
# k(site) = k_base × exp(+E_lat / kBT)
```

### 7.5 自洽性检查

**必须验证以下几点：**

1. **晶格常数收敛**：bulk Cu优化的a与实验值(3.615 Å)偏差 < 2%
2. **slab层数收敛**：比较4层和5层slab的E(clean)，差值 < 5 meV/atom
3. **k-points收敛**：p(3×3)比较3×3×1和5×5×1的E，差值 < 10 meV
4. **吸附位点确认**：检查弛豫后吸附物是否仍在原位点（没有迁移）
5. **对称性检查**：等价位点的E值应 < 5 meV差异
6. **镜像检查**：比较p(3×3)和p(4×4)的单吸附物能量，差值即为镜像误差

---

## 8. 预期结果与验证

### 8.1 文献参考值（类似体系）

| 吸附物对 | 表面 | ε_1NN (eV) | 来源 |
|---------|------|-----------|------|
| O\*–O\* | Pt(111) | +0.140 | Schneider, J. Catal. 286, 88 (2012) |
| CO\*–CO\* | Pt(111) | +0.560 | Zacros Tutorial 4 |
| CO\*–CO\* | Cu上 | +0.15~0.32 | 文献综合 |
| H\*–H\* | Pt(111) | +0.05~0.10 | 文献综合 |

你的Cu(111)上\*NO–\*NO结果应该在 **+0.05 ~ +0.30 eV** 范围内。如果超出此范围，检查计算设置。

### 8.2 物理合理性检查

| 检查项 | 预期 | 异常信号 |
|-------|------|---------|
| 同种排斥 | ε > 0 (排斥为主) | ε < −0.1 eV (过强吸引，检查结构) |
| 距离衰减 | \|ε_2NN\| < \|ε_1NN\| | ε_2NN > ε_1NN (不物理) |
| 异种交叉 | 通常排斥 | 强吸引可能暗示化学反应发生 |
| NH₃弱吸附 | ε ≈ 0 | \|ε\| > 0.1 eV (NH₃太弱不该有强相互作用) |

### 8.3 三体项检测

如果Tier 1的ε_1NN中有 > 0.15 eV的，建议补做三体项检查：

```
ω(A,B,C) = E(A+B+C) − E(A+B) − E(A+C) − E(B+C) + E(A) + E(B) + E(C) − E(clean)
```

如果 |ω| < 0.03 eV，三体项可忽略（pairwise近似成立）。

---

## 9. 参考文献

### 方法论
1. Curulla Ferré et al. "Testing the Pairwise Additive Potential Approximation" — *ChemPhysChem* **6**, 1009 (2005). DOI: 10.1002/cphc.200400399
2. Wu, Schmidt, Wolverton, Schneider. "Coverage-dependence in first-principles kinetic models: NO oxidation on Pt(111)" — *J. Catal.* **286**, 88 (2012). DOI: 10.1016/j.jcat.2011.10.020
3. Getman & Schneider. "DFT-Based Coverage-Dependent Model of Pt-Catalyzed NO Oxidation" — *ChemCatChem* **2**, 1450 (2010). DOI: 10.1002/cctc.201000146
4. Hess, Sander, Stamatakis. "Efficient Implementation of Cluster Expansion Models in Surface KMC" — *J. Comput. Chem.* **40**, 2664 (2019). DOI: 10.1002/jcc.26041
5. Schmidt et al. "Comparison of cluster expansion fitting algorithms" — *Surf. Sci.* (2015). DOI: 10.1016/j.susc.2015.01.019
6. Mason & Rappé. "Adsorbate-Adsorbate Interactions and Chemisorption at Different Coverages" — *J. Phys. Chem. B* **110**, 3816 (2006)

### Cu表面吸附
7. Xu, Lin, Mavrikakis. "Atomic and Molecular Adsorption on Cu(111)" — *Topics Catal.* **61**, 736 (2018). DOI: 10.1007/s11244-018-0943-0
8. Al Fauzan et al. "Coverage-dependent adsorption of NO on Cu(100)" — (2024). DOI: 10.1016/j.mseb.2024.117299

### KMC方法论
9. Zacros Tutorial 3 — Cluster Expansion for O/Pt(111) (zacros.org)
10. Zacros Tutorial 4 — Mapping DFT Energies to Zacros Input (zacros.org)
11. Andersen, Panosetti, Reuter. "A Practical Guide to Surface KMC Simulations" — *Front. Chem.* (2019). DOI: 10.3389/fchem.2019.00202

### NO₃RR
12. J. Phys. Chem. C (2025) — Cu(111) GC-DFT NO₃RR. DOI: 10.1021/acs.jpcc.5c02461
13. ACS Catal. 9, 7052 (2019) — NO₃RR MKM volcano. DOI: 10.1021/acscatal.9b02179
14. EES Catalysis (2024) — Cu-based NO₃RR review

---

*本文档作为SPARK Phase 2 DFT计算的操作指南*
*Last updated: 2026-03-22*
