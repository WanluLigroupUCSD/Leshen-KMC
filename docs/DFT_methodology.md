# NO₃RR on Cu(100): DFT 计算方法论

> SPARK Phase 2 | 2026-03-23 (D4 更新)
> 本文档描述**如何从DFT计算得到KMC所需的每一个能垒**
> 操作手册(VASP输入文件模板)见 `DFT_workflow_operational.md`
> GC-DFT详细教程见用户已有文档: `project4-gs-colab/Mo/transition/H/transition-energy/GC-DFT_tutorial.md`

---

## 目录

1. [计算框架总览](#1-计算框架总览)
2. [电子结构方法](#2-电子结构方法)
3. [表面模型](#3-表面模型)
4. [吸附能计算](#4-吸附能计算)
5. [过渡态搜索: CI-NEB](#5-过渡态搜索-ci-neb)
6. [恒电势能垒: GC-DFT](#6-恒电势能垒-gc-dft)
7. [PCET步骤的处理策略](#7-pcet步骤的处理策略)
8. [自由能修正](#8-自由能修正)
9. [横向相互作用计算](#9-横向相互作用计算)
10. [完整反应步骤计算矩阵](#10-完整反应步骤计算矩阵)
11. [从DFT到KMC的参数映射](#11-从dft到kmc的参数映射)

---

## 1. 计算框架总览

### 1.1 KMC需要什么参数？

KMC引擎中每个基元步骤的速率常数为：

```
k = (kBT/h) × exp(−Ea / kBT)       热力学步骤(TST)
k = (kBT/h) × exp(−Ea(U) / kBT)    电化学步骤(恒电势活化能)
```

因此每个步骤需要：
- **正向活化能 Ea,fwd(U)** — 恒电势下的反应速率
- **逆向活化能 Ea,rev(U)** — 可逆性 (Ea,rev = Ea,fwd − ΔΩ)
- **反应能 ΔE / ΔΩ(U)** — 热力学驱动力

### 1.2 核心方法选择

| 步骤类型 | 过渡态搜索 | 能垒计算 | 电位依赖 |
|---------|-----------|---------|---------|
| **热力学反应** (扩散、偶联、Tafel) | CI-NEB | GC-DFT | Ω(U)拟合 |
| **电化学反应 (PCET)** | CI-NEB (热力学类比) | GC-DFT | Ω(U)拟合 |
| **PCET (NEB不收敛时)** | — | BEP 备选 | Ea(U)=α×ΔG+Ea0+αeU |

### 1.3 计算流程总图

```
                    DFT计算流程 (D4)
                        │
         ┌──────────────┼──────────────┐
         │              │              │
    Phase 0-2:      Phase TS:       Phase GC-DFT:
    吸附能/ΔE        过渡态搜索       恒电势能垒
         │              │              │
    结构优化          CI-NEB         VASPsol +
    p(4×4) Cu(100)   (所有步骤)      变电荷单点
    (IBRION=2)        (VTST)        (NELECT±)
         │              │              │
         ▼              ▼              ▼
    E(slab+A)        E(TS)         Ω(U)二次拟合
    E(clean)         IS, TS结构       │
    E(gas)              │              │
         │              │              ▼
         │         频率计算         Ea(U) = Ω_TS(U)
         │         (ZPE)            − Ω_IS(U) + ΔZPE
         │              │              │
         └──────────────┼──────────────┘
                        │
                   KMC参数: Ea(U)
```

---

## 2. 电子结构方法

### 2.1 交换关联泛函

**PBE (Perdew-Burke-Ernzerhof) + DFT-D3(BJ) 色散修正**

```fortran
GGA  = PE          ! PBE 泛函
IVDW = 12          ! DFT-D3(BJ) 色散修正
```

**理由**: 横向相互作用ε和活化能Ea都是能量差，PBE的系统误差在取差时大部分抵消。D3(BJ)改善弱吸附物种(*NH₃, *N₂O)的描述。

### 2.2 赝势与基组

```
PAW方法 (Projector Augmented Wave)
POTCAR: PAW_PBE
  Cu: PAW_PBE Cu 22Jun2005      (ENMAX = 295 eV, 11价电子: 3d¹⁰4s¹)
  N:  PAW_PBE N  08Apr2002      (ENMAX = 400 eV)
  O:  PAW_PBE O  08Apr2002      (ENMAX = 400 eV)
  H:  PAW_PBE H  15Jun2001      (ENMAX = 250 eV)

ENCUT = 450 eV  (≥ max(ENMAX) × 1.125, 即400×1.125=450)
```

### 2.3 k-points采样

| 超胞 | k-points | 方法 |
|------|----------|------|
| Bulk Cu (常规胞) | 11×11×11 | Γ-centered |
| **p(4×4) slab** | **2×2×1** | **Γ-centered** |
| 气相分子 | 1×1×1 (Γ only) | — |

### 2.4 自旋极化

```fortran
ISPIN = 2  (所有计算统一开启)
```

**原因**: *NO是自由基(未配对电子)，吸附后磁矩可能淬灭也可能残留。统一开ISPIN=2保证一致性。

### 2.5 展宽方法

```fortran
ISMEAR = 0          ! Gaussian smearing (所有计算统一)
SIGMA  = 0.05       ! eV
```

**所有计算（slab、气相、GC-DFT单点）统一使用 ISMEAR=0, SIGMA=0.05。**

**检查**: entropy T×S应 < 1 meV/atom，否则SIGMA过大。使用 `energy(sigma->0)` 即E₀值。

---

## 3. 表面模型

### 3.1 Slab构建

```
Cu(100), fcc结构, 正方晶格
├── 晶格常数 a₀: PBE优化值 (先做bulk优化, 预期~3.63 Å)
├── 表面晶格常数: a_surf = a₀/√2 ≈ 2.566 Å
├── 层数: 4层 (底2层固定, 顶2层弛豫)
├── 真空层: 15 Å
└── 偶极修正: LDIPOL=.TRUE., IDIPOL=3 (非对称slab必须)
```

### 3.2 超胞

**统一使用 p(4×4) Cu(100)**，不做 p(3×3) 位点筛选。

- 每层 16 个 Cu 原子，4 层 = 64 Cu
- 超胞尺寸 ≈ 10.26 × 10.26 Å
- 对已知稳定位点，直接放置吸附物在 p(4×4) 上优化

### 3.3 吸附位点

Cu(100) 有 **3种** 高对称位点（正方晶格）：

```
  atop     bridge     hollow
   ●         ●         ●
   |        / \       /|\
   Cu     Cu   Cu   Cu Cu
                     Cu Cu
```

| 位点 | 配位数 | 描述 |
|------|--------|------|
| **Hollow** (4-fold) | 4 | 四个 Cu 原子围成的正方形中心 |
| **Bridge** | 2 | 两个相邻 Cu 之间 |
| **Atop** | 1 | 单个 Cu 正上方 |

**注意**: Cu(100) 只有一种 hollow 位点（不像 Cu(111) 有 fcc 和 hcp 两种）。

### 3.4 邻居距离

| 邻居 | Cu(100) 表面距离 | 说明 |
|------|----------------|------|
| 1NN | a₀/√2 ≈ 2.566 Å | 最近邻 (沿 [110] 方向) |
| 2NN | a₀ ≈ 3.63 Å | 次近邻 (沿 [100] 方向) |
| 3NN | a₀√(3/2) ≈ 4.44 Å | 通常忽略 |

---

## 4. 吸附能计算

### 4.1 定义

```
BE(A) = E(slab+A) − E(clean slab) − E(A, gas)
```

- BE < 0: 放热吸附（稳定）
- BE > 0: 吸热（不稳定）

### 4.2 反应能ΔE

对于表面反应 *A → *B + gas_product:

```
ΔE = [E(slab+B) + E(gas_product)] − E(slab+A)
```

对于 *A + *B(adj) → *C + *（双分子反应）:

```
ΔE = E(slab+C, after) − E(slab+A+B, before)
(在共吸附构型中直接计算)
```

### 4.3 电化学步骤的反应自由能 (CHE框架)

对于PCET步骤 *A + H⁺ + e⁻ → *B:

```
ΔG(U=0) = G(*B) − G(*A) − ½G(H₂, gas)    (CHE: μ(H⁺+e⁻) = ½G(H₂))
ΔG(U)   = ΔG(U=0) + eU                    (线性修正)
```

**注意**: CHE 仅用于快速估算 ΔG。活化能 Ea(U) 通过 GC-DFT 获得（见第6节）。

### 4.4 涉及H₂O的步骤

对于产生水的步骤 *A + H⁺ + e⁻ → *B + H₂O:

```
ΔG(U=0) = [G(*B) + G(H₂O, liq)] − G(*A) − ½G(H₂, gas)
```

G(H₂O, liq) 可用: E(H₂O, gas) + 0.035 eV (经验修正)

---

## 5. 过渡态搜索: CI-NEB

### 5.1 适用范围

**所有基元步骤均用 CI-NEB 搜索过渡态**，包括：
- 热力学步骤: 扩散、N-N偶联、Tafel
- PCET步骤: 通过热力学类比建模（*A + *H → *B，见第7节）

### 5.2 CI-NEB 方法

**原理**:
在初态(IS)和末态(FS)之间构建一系列中间构型(images)，用弹性弹簧连接，同时优化所有images使其收敛到最小能量路径(MEP)。Climbing image让能量最高的image向鞍点攀升，精确定位过渡态。

**步骤**:

```
Step 1: 优化初态(IS) → CONTCAR → 00/POSCAR
Step 2: 优化末态(FS) → CONTCAR → 0N/POSCAR (N=images+1)
Step 3: 线性插值 → nebmake.pl IS/POSCAR FS/POSCAR N_images
Step 4: 运行CI-NEB (VTST版本VASP)
Step 5: 提取 Ea,fwd = E(TS) − E(IS), Ea,rev = E(TS) − E(FS)
```

### 5.3 CI-NEB 的 VASP 参数

```fortran
! 基本设置
PREC    = Accurate
ENCUT   = 450
EDIFF   = 1E-5
ISIF    = 2
ISMEAR  = 0
SIGMA   = 0.05
ISPIN   = 2
LDIPOL  = .TRUE.
IDIPOL  = 3
LREAL   = Auto
GGA     = PE
IVDW    = 12          ! D3(BJ) — 与结构优化一致

! NEB专用
IBRION  = 3          ! VTST优化器
POTIM   = 0          ! 必须=0, VTST自己控制步长
ICHAIN  = 0          ! NEB方法
IMAGES  = 5          ! 中间image数 (不含IS和FS)
SPRING  = -5         ! 弹性常数 (eV/Å²), 负值=变弹簧
LCLIMB  = .TRUE.     ! 开启climbing image
IOPT    = 1          ! 优化器: 1=L-BFGS

! 收敛
EDIFFG  = -0.03      ! NEB力收敛 (eV/Å)
NSW     = 500        ! 最大离子步数
```

**⚠️ 必须使用VTST编译的VASP**: `grep "VTST" OUTCAR` 验证

### 5.4 各步骤的 CI-NEB 设置

#### 热力学步骤

| 步骤 | IS | FS | IMAGES | 预计Ea | 注意 |
|------|----|----|--------|--------|------|
| R12: *NO+*NO→*N₂O+* | 两个*NO共吸附1NN | *N₂O + 空位 | 5-7 | 0.05-0.80 eV | **最关键! N₂选择性核心** |
| R14: *N+*N→N₂(g)+2* | N_N_1NN | N₂(g) + 空 | 5 | ~0.45 eV | |
| R16: *H+*H→H₂+2* (Tafel) | H_H_1NN | H₂(g) + 空 | 3-5 | 0.8-1.0 eV | |
| R18: *NO扩散 | NO@hollow_1 | NO@hollow_2 | 3 | 0.1-0.3 eV | 经bridge |
| R19: *H扩散 | H@hollow_1 | H@hollow_2 | 3 | 0.1-0.2 eV | |

#### PCET步骤 (热力学类比 → CI-NEB)

| 步骤 | IS (类比) | FS (类比) | IMAGES | 理由 |
|------|----------|----------|--------|------|
| R2: *NO₃→*NO₂+*OH | *NO₃+*H | *NO₂+*OH | 5 | 脱氧 |
| R4: *NO₂→*NO+*OH | *NO₂+*H | *NO+*OH | 5 | 脱氧 |
| R5: *OH→H₂O | *OH+*H | H₂O+2* | 5 | OH清除 |
| **R6: *NO→*NOH** | ***NO+*H** | ***NOH+*** | **5** | **PLS! 必须精确** |
| R7: *NOH→*NHOH | *NOH+*H | *NHOH+* | 5 | 加氢 |
| R8: *NHOH→*NH+H₂O | *NHOH+*H | *NH+H₂O | 5 | 脱水，路径复杂 |
| R9: *NH→*NH₂ | *NH+*H | *NH₂+* | 5 | 加氢 |
| R10: *NH₂→*NH₃ | *NH₂+*H | *NH₃+* | 5 | 加氢 |
| R15: Volmer | H₂O+* | *H+OH⁻ | 5 | HER |

### 5.5 TS验证

**必须验证过渡态的正确性**:

```
1. 频率分析: IBRION=5 → 有且仅有一个虚频
2. 能量检查: E(TS) > E(IS) 且 E(TS) > E(FS)
3. 结构合理性: 可视化TS, 检查键长/键角
```

频率计算INCAR:
```fortran
IBRION  = 5          ! 有限差分法
NSW     = 1
NFREE   = 2          ! 中心差分
POTIM   = 0.015      ! 位移步长 (Å)
EDIFF   = 1E-7       ! 频率计算需更严格的电子收敛
ISMEAR  = 0
SIGMA   = 0.05
! 仅振动吸附物原子 (slab F F F)
```

### 5.6 备选: Dimer方法

若NEB不收敛, 可从NEB最高image出发用Dimer方法:

```fortran
ICHAIN  = 2          ! Dimer方法
IOPT    = 2          ! CG优化器
DdR     = 0.005      ! 有限差分步长 (Å)
DRotMax = 4          ! 每步最大旋转次数
```

---

## 6. 恒电势能垒: GC-DFT

### 6.1 原理

**所有步骤**（热力学和电化学）均通过 GC-DFT 获取恒电势活化能 Ea(U)。

常规DFT在固定电子数(CE)下计算，但电催化反应发生在恒电势条件下。GC-DFT将体系从正则系综转换到巨正则系综：

```
巨正则自由能:
Ω(U) = E(N_e) + ΔN_e × Φ

其中:
  E(N_e) — 给定电子数下的DFT总能 (energy without entropy)
  ΔN_e = N_e − N_e^0 — 相对于中性体系的多余电子数
  Φ — 功函数 = E_vacuum − E_Fermi
  U_SHE = Φ − φ_SHE (φ_SHE = 4.6 eV)
```

### 6.2 工作流程

对每个反应步骤的 IS 和 TS 结构分别执行:

```
Step 1: 确定 NELECT_neutral (中性体系电子数)
        ← grep NELECT OUTCAR 或 从POTCAR的ZVAL计算
           Cu: ZVAL=11, N: 5, O: 6, H: 1

Step 2: 变电荷单点能计算 (7-9个点)
        NELECT = N_e^0 + offset
        offset = -1.5, -1.0, -0.5, 0, +0.5, +1.0, +1.5
        (资源充足可扩展到 ±2.0, 共9个点)

Step 3: 提取功函数 Φ
        vaspkit 426 或 从LOCPOT/PLANAR_AVERAGE.dat

Step 4: 计算 Ω(U) 并二次拟合
        Ω = E + ΔN_e × Φ
        U_SHE = Φ − 4.6
        拟合: Ω(U) = aU² + bU + c (要求 R² > 0.99)

Step 5: 恒电势活化能
        Ea(U) = Ω_TS(U) − Ω_IS(U) + ΔZPE
```

### 6.3 GC-DFT 单点能 INCAR

```fortran
! ===== 电子步 =====
ALGO   = F          ! Fast算法
EDIFF  = 1e-7       ! 比结构优化更严格
NELM   = 500        ! 带电体系可能更难收敛
NELMIN = 6
PREC   = Accurate

! ===== 离子步 =====
IBRION = -1         ! 不弛豫
NSW    = 0          ! 单点能

! ===== 电子占据 =====
ISMEAR = 0          ! Gaussian smearing
SIGMA  = 0.05

! ===== 自旋与对称性 =====
ISPIN  = 2
ISYM   = -1         ! 关闭对称性

! ===== 泛函与色散 =====
GGA    = PE         ! PBE
ENCUT  = 450        ! 与结构优化一致!
IVDW   = 12         ! DFT-D3(BJ)

! ===== VASPsol 隐式溶剂化 =====
LSOL       = .TRUE.     ! 开启溶剂化
EB_K       = 78.4        ! 水的介电常数
TAU        = 0           ! 腔体表面张力
LAMBDA_D_K = 3.0         ! Debye 屏蔽长度 (Å)

! ===== 关键输出 =====
LVHAR  = .TRUE.     ! 输出 LOCPOT (计算功函数!)
LCHARG = .TRUE.
LORBIT = 11
LREAL  = Auto

! ===== 偶极校正 =====
DIPOL = 0.5 0.5 0.5

! ===== 变量 =====
NELECT = {nelect}   ! ← 每个电荷态不同
```

### 6.4 GC-DFT 目录结构

```
charge/
├── IS/
│   ├── R6_NO-NOH/          ← 反应步骤名
│   │   ├── charge_-1.5/    ← NELECT = N_e^0 − 1.5
│   │   │   ├── INCAR
│   │   │   ├── POSCAR      ← IS优化后的CONTCAR
│   │   │   ├── POTCAR
│   │   │   └── KPOINTS
│   │   ├── charge_-1.0/
│   │   ├── charge_-0.5/
│   │   ├── charge_0/       ← 中性参考
│   │   ├── charge_0.5/
│   │   ├── charge_1.0/
│   │   └── charge_1.5/
│   └── ...
├── TS/
│   ├── R6_NO-NOH/
│   │   ├── charge_-1.5/
│   │   │   ├── POSCAR      ← CI-NEB的TS结构(CONTCAR)!
│   │   │   └── ...
│   │   └── ...
│   └── ...
```

### 6.5 数据提取与拟合

```python
import numpy as np
import re

phi_SHE = 4.6  # eV, SHE参考电势

def get_energy_without_entropy(outcar_path):
    """从OUTCAR提取 energy without entropy"""
    energy = None
    with open(outcar_path) as f:
        for line in f:
            if "energy  without entropy" in line:
                m = re.search(r"energy\s+without entropy=\s*([-\d.]+)", line)
                if m:
                    energy = float(m.group(1))
    return energy

# 对每个电荷态:
# Ω = E + ΔN_e × Φ
# U_SHE = Φ − φ_SHE

# 二次拟合:
coeffs = np.polyfit(U_array, Omega_array, 2)  # Ω(U) = aU² + bU + c
# 要求 R² > 0.99

# 恒电势活化能:
# Ea(U) = Ω_TS(U) − Ω_IS(U) + ΔZPE
```

### 6.6 计算量估算

每个反应步骤的 GC-DFT:
- IS: 7 个变电荷单点能
- TS: 7 个变电荷单点能
- 小计: **14 个单点能/步骤**
- 每个单点 ~0.5-2 小时

全部 ~15 个反应步骤: **~210 个单点能** (可全部并行提交)

---

## 7. PCET步骤的处理策略

### 7.1 首选: 热力学类比 + CI-NEB + GC-DFT

将 PCET 步骤 `*A + H⁺ + e⁻ → *B` 建模为热力学类比:

```
*A + *H(邻位) → *B + *(空位)
```

即把质子当作已在表面上的*H，做热力学表面反应的CI-NEB。

**然后对 IS 和 TS 做 GC-DFT**，直接得到 Ea(U)。

**优点**: 真实TS几何 + 恒电势活化能，最高精度
**缺点**: PCET的过渡态较难找到（势能面平坦，NEB可能不收敛）

### 7.2 备选: BEP 关系

若某个 PCET 步骤的 CI-NEB 始终不收敛，回退到 BEP:

```
Ea(U) = α × ΔG(U) + Ea0
      = α × [ΔG(U=0) + eU] + Ea0
      = [α×ΔG(U=0) + Ea0] + α×eU

其中:
  α ≈ 0.5 (BEP斜率 ≈ 转移系数 β)
  Ea0 ≈ 0.5-0.8 eV (需校准)
  ΔG(U=0) 从 CHE 计算:
    ΔG(U=0) = E(*B) − E(*A) − ½E(H₂) + ΔZPE − TΔS
```

**BEP 精度有限** (~0.1-0.3 eV 误差)，仅作为不得已的备选方案。

### 7.3 各 PCET 步骤策略

| 步骤 | 首选 | 备选(若NEB不收敛) | 优先级 |
|------|------|-------------------|--------|
| **R6: *NO→*NOH** | **CI-NEB + GC-DFT** | — (必须成功) | **★★★★★ PLS** |
| R2: *NO₃→*NO₂+*OH | CI-NEB + GC-DFT | BEP | ★★★ |
| R4: *NO₂→*NO+*OH | CI-NEB + GC-DFT | BEP | ★★★ |
| R5: *OH→H₂O | CI-NEB + GC-DFT | BEP | ★★★ |
| R7: *NOH→*NHOH | CI-NEB + GC-DFT | BEP | ★★★ |
| R8: *NHOH→*NH+H₂O | CI-NEB + GC-DFT | BEP | ★★★ |
| R9: *NH→*NH₂ | CI-NEB + GC-DFT | BEP | ★★ |
| R10: *NH₂→*NH₃ | CI-NEB + GC-DFT | BEP | ★★ |
| R15: Volmer | CI-NEB + GC-DFT | BEP | ★★★ |

---

## 8. 自由能修正

### 8.1 零点能修正 (ZPE)

```
ZPE = ½ Σᵢ hνᵢ    (对所有实数频率求和)
```

**计算方法**: 在优化构型上做频率计算(IBRION=5), 只振动吸附物原子(固定slab)

```fortran
IBRION = 5
NSW    = 1
NFREE  = 2
POTIM  = 0.015
EDIFF  = 1E-7       ! 频率计算需更严格电子收敛
ISMEAR = 0
SIGMA  = 0.05
```

**频率计算在真空条件下进行即可**（不加溶剂化），因为溶剂效应主要影响电子能量而非振动频率。

**ZPE 用于 GC-DFT**: ΔZPE = ZPE(TS) − ZPE(IS)

典型ZPE值 (文献参考):

| 物种 | ZPE (eV) |
|------|----------|
| *H | ~0.17 |
| *OH | ~0.35 |
| *NO | ~0.09 |
| ½H₂(gas) | 0.135 |
| H₂O(gas) | 0.56 |

### 8.2 熵修正

```
G = E + ZPE − TS
```

**表面吸附物**: 振动熵通常很小 (~0.01-0.05 eV at 298K), 可近似忽略。

**气相分子**: 平动+转动+振动熵 (从NIST热力学数据库获取)

| 物种 | −TS at 298K, 1 bar (eV) |
|------|----------------------|
| H₂(g) | −0.40 |
| H₂O(g) | −0.67 |
| NH₃(g) | −0.60 |
| N₂(g) | −0.59 |
| N₂O(g) | −0.66 |
| NO(g) | −0.63 |

### 8.3 溶剂效应

**GC-DFT 计算中已包含 VASPsol 隐式溶剂化**（LSOL=.TRUE.）。
溶剂效应在 GC-DFT 变电荷单点能中自动计入，无需额外修正。

---

## 9. 横向相互作用计算

### 9.1 公式

```
ε(A,B,d) = E(slab+A+B@d) − E(slab+A) − E(slab+B) + E(clean)
```

4个结构优化得到1个ε值。单吸附能和clean slab能量可在不同对之间复用。

### 9.2 需要计算的对

**Tier 1** (必须, 4对×2距离=8个共吸附计算):
- *NO–*NO, *NO–*H, *H–*H, *OH–*NO

**Tier 2** (建议, 6对×2距离=12个):
- *OH–*OH, *OH–*H, *N–*NO, *N–*N, *O–*NO, *N–*H

### 9.3 Cu(100)上的距离

| 邻居 | 距离 | 位点排列 |
|------|------|---------|
| **1NN** | **a₀/√2 ≈ 2.566 Å** | 相邻hollow位点 (沿[110]) |
| **2NN** | **a₀ ≈ 3.63 Å** | 对角hollow位点 (沿[100]) |

共吸附时检查:
1. 吸附物未迁移/解离
2. 磁矩变化
3. 构型稳定性

### 9.4 三体项检测

如果Tier 1中任何 |ε_1NN| > 0.15 eV, 建议检测三体项:

```
ω(A,B,C) = E(A+B+C) − E(A+B) − E(A+C) − E(B+C) + E(A) + E(B) + E(C) − E(clean)

如果 |ω| < 0.03 eV → 成对可加性近似成立
```

---

## 10. 完整反应步骤计算矩阵

### 10.1 总计算量

| 类别 | 内容 | DFT计算数 | 方法 |
|------|------|----------|------|
| Phase 0 | bulk + 8气相 | 9 | 结构优化 |
| Phase 1 | 稳定位点吸附 (p(4×4)直接优化) | ~10-13 | 结构优化 |
| Phase 2 | Tier 1共吸附 | 8 | 结构优化 |
| Phase 3 | Tier 2共吸附 | 12 | 结构优化 |
| Phase TS | 所有步骤CI-NEB | ~15×(5 images+2) ≈ 105 | CI-NEB |
| Phase FREQ | IS+TS频率 | ~30 | IBRION=5 |
| **Phase GC-DFT** | **IS+TS 变电荷单点** | **~15×14 = 210** | **GC-DFT** |
| **总计** | | **~390-400** | |

### 10.2 每步最终参数表模板

```
步骤        | Ea(U) 方法        | 过渡态  | 关键性
──────────────────────────────────────────────────
R1  NO₃吸附  | 无垒              | —      |
R2  脱氧1    | CI-NEB+GC-DFT    | NEB    | ★★★
R4  脱氧2    | CI-NEB+GC-DFT    | NEB    | ★★★
R5  OH清除   | CI-NEB+GC-DFT    | NEB    | ★★★
R6  *NO→*NOH | CI-NEB+GC-DFT    | NEB    | ★★★★★ PLS
R7  加氢     | CI-NEB+GC-DFT    | NEB    | ★★★
R8  脱水     | CI-NEB+GC-DFT    | NEB    | ★★★
R9  加氢     | CI-NEB+GC-DFT    | NEB    | ★★
R10 加氢     | CI-NEB+GC-DFT    | NEB    | ★★
R11 NH₃脱附  | |BE|             | —      |
R12 N-N偶联  | CI-NEB+GC-DFT    | NEB    | ★★★★★ N₂选择性
R13 N₂O脱附  | |BE|             | —      |
R14 N+N→N₂  | CI-NEB+GC-DFT    | NEB    | ★★★
R15 Volmer   | CI-NEB+GC-DFT    | NEB    | ★★★
R16 Tafel    | CI-NEB+GC-DFT    | NEB    | ★★★
R18 NO扩散   | CI-NEB+GC-DFT    | NEB    | ★★★
R19 H扩散    | CI-NEB+GC-DFT    | NEB    | ★★

若PCET步骤NEB不收敛 → 回退BEP: Ea=α×ΔG+Ea0, α≈0.5
```

### 10.3 推荐计算顺序

```
阶段1: Phase 0 (bulk Cu + 气相分子)
  ↓ 得到 a₀ 和气相参考能量
阶段2: Phase 1 (p(4×4) 吸附物结构优化, 直接放稳定位点)
  ↓ 得到所有 E(slab+A) 和 E(clean)
阶段3 (并行): 共吸附(ε) + NEB末态(FS)优化
  ↓
阶段4: CI-NEB (所有步骤)
  ↓ 得到 TS 结构
阶段5: 频率计算 (IS + TS)
  ↓ 得到 ZPE
阶段6: GC-DFT 变电荷单点 (IS + TS, 每个7点)
  ↓ 得到 Ω(U) → Ea(U)
阶段7: 数据整合 → KMC参数
```

---

## 11. 从DFT到KMC的参数映射

### 11.1 恒电势速率常数

```python
# 所有步骤统一用 GC-DFT 的 Ea(U):
# k(U) = (kBT/h) × exp(−Ea(U) / kBT)

# Ea(U) = Ω_TS(U) − Ω_IS(U) + ΔZPE
# 其中 Ω(U) = aU² + bU + c (二次拟合)

# KMC中，对给定电位U:
# 直接代入拟合公式得到 Ea(U) → k(U)
```

### 11.2 脱附步骤

```python
# 脱附无需TS搜索:
# Ea_des ≈ |BE(A)|
# 例: *NH₃ → NH₃(aq), Ea_des ≈ 0.37 eV
```

### 11.3 横向相互作用

```python
# KMC速率修正:
# k(site) = k_base × exp(+E_lat / kBT)
# E_lat = Σ_j ε(species_i, species_j, d_ij)

# ε值从共吸附计算得到
# 正值(排斥) → 加速反应/脱附
# 负值(吸引) → 减慢脱附
```

---

*本文档为DFT计算的方法论框架 (D4 更新: Cu(100) + GC-DFT)*
*GC-DFT教程见: `project4-gs-colab/Mo/transition/H/transition-energy/GC-DFT_tutorial.md`*
*VASP输入文件模板和操作脚本见: `DFT_workflow_operational.md`*
*整体研究计划见: `research_plan_NO3RR.md`*
*Last updated: 2026-03-23 (D4)*
