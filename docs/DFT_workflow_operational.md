# NO₃RR on Cu(111): DFT 计算操作手册

> SPARK Phase 2 | 2026-03-23
> 配合 `DFT_lateral_interactions_guide.md` 使用（理论公式），本文档为**逐步操作指南**

---

## 目录

1. [项目目录结构](#1-项目目录结构)
2. [Phase 0: Bulk Cu + 气相分子](#2-phase-0-bulk-cu--气相分子)
3. [Phase 1: 位点偏好测试](#3-phase-1-位点偏好测试)
4. [Phase 2: 精确单吸附物](#4-phase-2-精确单吸附物)
5. [Phase 3-4: 成对横向相互作用](#5-phase-3-4-成对横向相互作用)
6. [Phase TS: 过渡态搜索](#6-phase-ts-过渡态搜索)
7. [数据提取与后处理脚本](#7-数据提取与后处理脚本)
8. [收敛性检查清单](#8-收敛性检查清单)
9. [常见问题排查](#9-常见问题排查)
10. [附录: 坐标生成参考](#10-附录-坐标生成参考)

---

## 1. 项目目录结构

```
DFT_NO3RR_Cu111/
├── phase0_baseline/
│   ├── bulk_Cu/               # bulk Cu fcc优化
│   │   ├── INCAR
│   │   ├── POSCAR
│   │   ├── KPOINTS
│   │   ├── POTCAR
│   │   └── run.sh
│   └── gas_phase/             # 气相分子
│       ├── NO/
│       ├── H2/
│       ├── H2O/
│       ├── NH3/
│       ├── N2/
│       ├── N2O/
│       ├── NO2/
│       └── OH/
│
├── phase1_site_test/
│   ├── clean_slab_3x3/       # 清洁p(3×3) slab
│   ├── H_fcc/                 # *H在fcc hollow
│   ├── H_hcp/                 # *H在hcp hollow
│   ├── H_bridge/              # *H在bridge
│   ├── H_atop/                # *H在atop
│   ├── NO_fcc/
│   ├── NO_hcp/
│   ├── NO_bridge/
│   ├── NO_atop/
│   ├── ... (每个物种×4位点)
│   └── summary_site_test.csv  # 汇总结果
│
├── phase2_single_ads/
│   ├── clean_slab_4x4/       # 清洁p(4×4) slab
│   ├── NO_fcc_4x4/           # *NO在最优位点(fcc), p(4×4)
│   ├── H_fcc_4x4/
│   ├── OH_fcc_4x4/
│   ├── O_fcc_4x4/
│   ├── N_fcc_4x4/
│   ├── NH_fcc_4x4/
│   ├── NH2_bridge_4x4/
│   └── summary_single_ads.csv
│
├── phase3_tier1_pairs/
│   ├── NO_NO_1NN/             # *NO+*NO共吸附, 1NN距离
│   ├── NO_NO_2NN/
│   ├── NO_H_1NN/
│   ├── NO_H_2NN/
│   ├── H_H_1NN/
│   ├── H_H_2NN/
│   ├── OH_NO_1NN/
│   ├── OH_NO_2NN/
│   └── summary_lateral_tier1.csv
│
├── phase4_tier2_pairs/
│   ├── OH_OH_1NN/
│   ├── ... (6对 × 2距离 = 12个)
│   └── summary_lateral_tier2.csv
│
├── phase_TS/
│   ├── NEB_NO_NO_coupling/    # *NO+*NO→*N₂O CI-NEB
│   ├── NEB_H_H_tafel/        # Tafel步骤
│   ├── NEB_NO_diffusion/      # *NO扩散
│   ├── NEB_H_diffusion/       # *H扩散
│   └── summary_barriers.csv
│
├── scripts/
│   ├── extract_energy.sh      # 批量提取能量
│   ├── calc_epsilon.py        # 计算横向相互作用ε
│   ├── check_convergence.sh   # 检查收敛
│   ├── gen_poscar_slab.py     # 生成slab POSCAR
│   └── gen_poscar_coads.py    # 生成共吸附POSCAR
│
└── results/
    ├── all_energies.csv       # 所有DFT总能量汇总
    ├── lateral_interactions.csv # 所有ε值
    ├── barriers.csv           # 所有活化能
    └── kmc_params.json        # 最终KMC输入参数
```

---

## 2. Phase 0: Bulk Cu + 气相分子

### 2.1 Bulk Cu 优化

**目标**: 获取PBE优化的Cu晶格常数 a₀

**POSCAR** (bulk_Cu/POSCAR):
```
Cu bulk fcc
1.0
  3.615  0.000  0.000
  0.000  3.615  0.000
  0.000  0.000  3.615
Cu
4
Direct
  0.000  0.000  0.000
  0.500  0.500  0.000
  0.500  0.000  0.500
  0.000  0.500  0.500
```

**INCAR** (bulk_Cu/INCAR):
```fortran
SYSTEM  = bulk_Cu_fcc
PREC    = Accurate
ENCUT   = 450
EDIFF   = 1E-6
EDIFFG  = -0.005
IBRION  = 2
NSW     = 50
ISIF    = 3          ! ← 关键: ISIF=3 允许晶胞形状和体积变化
ISMEAR  = 1
SIGMA   = 0.15
ALGO    = Fast
ISPIN   = 1          ! Cu无磁性
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

**KPOINTS** (bulk_Cu/KPOINTS):
```
K-Points
0
Gamma
11 11 11
0  0  0
```

**POTCAR**: Cu PAW_PBE

**操作步骤**:
```bash
# 1. 创建目录
mkdir -p DFT_NO3RR_Cu111/phase0_baseline/bulk_Cu
cd DFT_NO3RR_Cu111/phase0_baseline/bulk_Cu

# 2. 准备输入文件 (INCAR, POSCAR, KPOINTS 如上)
# 3. 生成POTCAR
cat /path/to/potpaw_PBE/Cu/POTCAR > POTCAR

# 4. 提交任务
sbatch run.sh   # 或你的调度系统

# 5. 完成后提取晶格常数
grep "ALAT" OUTCAR            # 如果有
# 或从CONTCAR读取:
head -5 CONTCAR               # 第2-4行是晶胞向量
# a₀ = 第2行第1列的值
# PBE预期: a₀ ≈ 3.63 Å (实验值 3.615 Å, 偏差 <0.5%)
```

**验证**: |a₀(DFT) − 3.615| / 3.615 < 2%

---

### 2.2 气相分子优化

**目标**: 获取各气相分子的总能量 E(gas), 用于计算吸附能

**通用INCAR** (gas_phase/INCAR_template):
```fortran
SYSTEM  = gas_phase_molecule
PREC    = Accurate
ENCUT   = 450
EDIFF   = 1E-6
EDIFFG  = -0.01
IBRION  = 2
NSW     = 100
ISIF    = 2          ! 不优化晶胞
ISMEAR  = 0          ! ← 分子用Gaussian展宽
SIGMA   = 0.01       ! ← 小展宽
ALGO    = Fast
NELM    = 200
ISPIN   = 2          ! ← 默认开自旋, 安全
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

**KPOINTS** (所有气相分子共用):
```
Gamma only
0
Gamma
1 1 1
0 0 0
```

**各分子的POSCAR**:

**NO** (gas_phase/NO/POSCAR) — **自旋极化！**
```
NO molecule
1.0
  15.000   0.000   0.000
   0.000  16.000   0.000
   0.000   0.000  17.000
N O
1 1
Cartesian
  7.500  8.000  8.500    ! N
  7.500  8.000  9.660    ! O (N-O键长 ~1.16 Å)
```
> 注意: box三个方向不同(15×16×17)，避免对称性问题

**H₂** (gas_phase/H2/POSCAR):
```
H2 molecule
1.0
  15.000   0.000   0.000
   0.000  16.000   0.000
   0.000   0.000  17.000
H
2
Cartesian
  7.500  8.000  8.125    ! H
  7.500  8.000  8.875    ! H (H-H键长 ~0.75 Å)
```

**H₂O** (gas_phase/H2O/POSCAR):
```
H2O molecule
1.0
  15.000   0.000   0.000
   0.000  16.000   0.000
   0.000   0.000  17.000
H O
2 1
Cartesian
  7.500  8.000  8.500    ! O
  7.500  8.763  9.058    ! H (O-H=0.97 Å, 角度104.5°)
  7.500  7.237  9.058    ! H
```

**NH₃** (gas_phase/NH3/POSCAR):
```
NH3 molecule
1.0
  15.000   0.000   0.000
   0.000  16.000   0.000
   0.000   0.000  17.000
N H
1 3
Cartesian
  7.500  8.000  8.500    ! N
  7.500  8.942  8.820    ! H (N-H=1.02 Å, 三角锥)
  7.500  7.529  9.315    ! H
  8.315  7.529  8.185    ! H
```

**N₂** (gas_phase/N2/POSCAR):
```
N2 molecule
1.0
  15.000   0.000   0.000
   0.000  16.000   0.000
   0.000   0.000  17.000
N
2
Cartesian
  7.500  8.000  7.950    ! N
  7.500  8.000  9.050    ! N (N≡N键长 ~1.10 Å)
```

**N₂O** (gas_phase/N2O/POSCAR) — 线性分子 N=N=O:
```
N2O molecule
1.0
  15.000   0.000   0.000
   0.000  16.000   0.000
   0.000   0.000  17.000
N O
2 1
Cartesian
  7.500  8.000  7.370    ! N (端N)
  7.500  8.000  8.500    ! N (中心N, N-N=1.13 Å)
  7.500  8.000  9.686    ! O (N-O=1.186 Å)
```

**NO₂** (gas_phase/NO2/POSCAR) — 弯曲分子, **自旋极化**:
```
NO2 molecule
1.0
  15.000   0.000   0.000
   0.000  16.000   0.000
   0.000   0.000  17.000
N O
1 2
Cartesian
  7.500  8.000  8.500    ! N
  7.500  8.998  9.100    ! O (N-O=1.20 Å, O-N-O=134°)
  7.500  7.002  9.100    ! O
```

**OH** (gas_phase/OH/POSCAR) — **自由基, 自旋极化!**
```
OH radical
1.0
  15.000   0.000   0.000
   0.000  16.000   0.000
   0.000   0.000  17.000
O H
1 1
Cartesian
  7.500  8.000  8.500    ! O
  7.500  8.000  9.470    ! H (O-H=0.97 Å)
```

**操作步骤**:
```bash
# 批量创建目录
for mol in NO H2 H2O NH3 N2 N2O NO2 OH; do
  mkdir -p DFT_NO3RR_Cu111/phase0_baseline/gas_phase/$mol
  cp INCAR_template DFT_NO3RR_Cu111/phase0_baseline/gas_phase/$mol/INCAR
  cp KPOINTS_gamma  DFT_NO3RR_Cu111/phase0_baseline/gas_phase/$mol/KPOINTS
  # 手动放置各分子的POSCAR
  # 生成对应的POTCAR (注意元素顺序与POSCAR一致!)
done

# 批量提交
for mol in NO H2 H2O NH3 N2 N2O NO2 OH; do
  cd DFT_NO3RR_Cu111/phase0_baseline/gas_phase/$mol
  sbatch run.sh
  cd -
done

# 批量提取能量
for mol in NO H2 H2O NH3 N2 N2O NO2 OH; do
  E=$(grep "energy(sigma->0)" \
    DFT_NO3RR_Cu111/phase0_baseline/gas_phase/$mol/OUTCAR \
    | tail -1 | awk '{print $NF}')
  echo "$mol  $E"
done > phase0_gas_energies.csv
```

**验证**:
- H₂: E ≈ −6.77 eV (PBE)
- H₂O: E ≈ −14.22 eV
- N₂: E ≈ −16.64 eV
- NO: 检查ISPIN=2有效, OUTCAR中 `magnetization` 应 ≈ 1 μ_B

---

## 3. Phase 1: 位点偏好测试

### 3.1 构建清洁 p(3×3) Cu(111) slab

**目标**: 用Phase 0得到的a₀构建slab, 获取E(clean slab)

**POSCAR生成方法** (推荐用ASE):

```python
# scripts/gen_poscar_slab.py
from ase.build import fcc111
from ase.io import write
from ase.constraints import FixAtoms
import numpy as np

# --- 用户修改区 ---
a0 = 3.63       # ← 填入Phase 0优化的晶格常数!
size = (3, 3)   # p(3×3)
nlayers = 4
vacuum = 15.0   # Å
# -------------------

slab = fcc111('Cu', size=(*size, nlayers), a=a0, vacuum=vacuum,
              periodic=True)

# 固定底2层
z_coords = slab.positions[:, 2]
z_sorted = np.sort(np.unique(np.round(z_coords, 2)))
z_cutoff = (z_sorted[1] + z_sorted[2]) / 2  # 第2层和第3层之间
fix = FixAtoms(indices=[i for i, z in enumerate(z_coords) if z < z_cutoff])
slab.set_constraint(fix)

write('POSCAR', slab, format='vasp', vasp5=True, sort=True)
print(f"Slab: {len(slab)} atoms, cell: {slab.cell.lengths()}")
print(f"Fixed atoms: {sum(1 for z in z_coords if z < z_cutoff)}")

# 打印fcc hollow位点坐标 (用于后续放置吸附物)
top_z = z_sorted[-1]
ads_z = top_z + 1.8  # 吸附物大约在表面上方1.8 Å
print(f"\n吸附物放置高度: z ≈ {ads_z:.3f} Å")

# 生成fcc hollow位点坐标
a_nn = a0 / np.sqrt(2)  # 最近邻距离
print(f"最近邻距离: {a_nn:.3f} Å")
```

**INCAR** (所有slab计算通用):
```fortran
SYSTEM  = Cu111_3x3_clean
PREC    = Accurate
ENCUT   = 450
EDIFF   = 1E-5
EDIFFG  = -0.02
IBRION  = 2
NSW     = 300
ISIF    = 2
ISMEAR  = 1
SIGMA   = 0.15
ALGO    = Fast
NELM    = 200
ISPIN   = 2          ! 全部开自旋（即使清洁slab）
LDIPOL  = .TRUE.     ! 偶极修正
IDIPOL  = 3
LREAL   = Auto       ! 大超胞加速
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

**KPOINTS**:
```
K-Points
0
Gamma
3 3 1
0 0 0
```

---

### 3.2 在各位点放置吸附物

**Cu(111) 四种高对称位点坐标** (相对于表面Cu原子):

```
俯视图 (●=顶层Cu, ○=第2层Cu):

       ●───────●───────●
      / \  hcp/ \  fcc/ \
     /   \ ▲ /   \ ▼ /   \
    /     \/     \/     \
   ●───br──●───br──●───br──●
    \     /\     /\     /
     \   / ▼\   / ▲\   /
      \ / fcc\ / hcp\ /
       ●───────●───────●

▼ fcc hollow: 正下方无第2层原子
▲ hcp hollow: 正下方有第2层原子
br: bridge — 两个顶层Cu之间
atop: 直接在顶层Cu原子正上方
```

**用ASE放置吸附物** (推荐):

```python
# scripts/add_adsorbate.py
from ase.io import read, write
from ase.build import add_adsorbate
from ase import Atoms

slab = read('clean_slab_3x3/CONTCAR')  # 用弛豫后的清洁slab

# === 放置 *H 在 fcc hollow ===
# fcc hollow位点: 用fractional坐标 (1/3, 1/3) 相对于单胞
add_adsorbate(slab, 'H', height=1.0, position='fcc')
write('H_fcc/POSCAR', slab, format='vasp', vasp5=True)
# height: H距表面顶层Cu的距离 (Å), H一般~1.0 Å

# === 放置 *NO 在 fcc hollow ===
slab2 = read('clean_slab_3x3/CONTCAR')
# NO分子: N朝下(朝表面), O在上
no_mol = Atoms('NO', positions=[[0, 0, 0], [0, 0, 1.16]])  # N在下
add_adsorbate(slab2, no_mol, height=1.3, position='fcc')
write('NO_fcc/POSCAR', slab2, format='vasp', vasp5=True)
# height: N距表面的距离, NO约1.3 Å
```

**手动放置方法** (不用ASE时):

对于p(3×3) Cu(111), 晶胞向量为:
```
a1 = (3×a_nn,    0,       0) = (7.669, 0, 0)
a2 = (1.5×a_nn, 1.5×a_nn×√3, 0) = (3.835, 6.641, 0)
```

各位点的**分数坐标** (x, y相对于a1, a2):

| 位点 | 分数坐标 (x, y) | 说明 |
|------|----------------|------|
| atop | (0, 0) | 顶层Cu原子正上方 |
| fcc hollow | (1/3, 1/3) | 三个Cu围成的三角, 下方无第2层Cu |
| hcp hollow | (2/3, 2/3) | 三个Cu围成的三角, 下方有第2层Cu |
| bridge | (1/6, 1/6) | 两个Cu之间 (近似) |

**吸附物距表面的初始高度** (z方向, 粗略参考):

| 物种 | 朝向 | 距表面高度 (Å) | 注意 |
|------|------|---------------|------|
| *H | H朝下 | ~1.0 | 最简单 |
| *O | O朝下 | ~1.2 | — |
| *N | N朝下 | ~1.2 | — |
| *OH | O朝下, H斜上 | ~1.3 (O到表面) | O-H键 ~0.97 Å, 倾斜 ~30° |
| *NO | **N朝下**, O在上 | ~1.3 (N到表面) | **分子近垂直**, N-O=1.16 Å |
| *NOH | N朝下, O-H斜上 | ~1.3 | 在NO基础上加H |
| *NH | N朝下, H斜上 | ~1.2 | — |
| *NH₂ | N朝下, 2H斜上 | ~1.4 (bridge时) | 平面构型 |
| *NH₃ | N朝下, 3H上方 | ~2.1 (atop) | 弱吸附, 距离远 |
| *NO₂ | 两个O朝下 | ~1.8 | bidentate, 两个O桥接两个Cu |
| *NO₃ | 三个O朝下 | ~2.0 | tridentate |
| *NHOH | N朝下 | ~1.3 | — |
| *N₂O | N-N-O线性 | ~2.0 (atop) | 可能弱吸附 |

---

### 3.3 批量计算工作流

```bash
# === Step 1: 准备所有位点测试目录 ===
SPECIES="H O N OH NO NOH NH NH2 NH3 NO2 NO3 NHOH N2O"
SITES="fcc hcp bridge atop"

for sp in $SPECIES; do
  for site in $SITES; do
    dir="phase1_site_test/${sp}_${site}"
    mkdir -p $dir
    # 复制通用INCAR和KPOINTS
    cp templates/INCAR_slab_3x3 $dir/INCAR
    cp templates/KPOINTS_3x3    $dir/KPOINTS
    # POSCAR: 需要用ASE脚本或手动生成
    # POTCAR: 按物种组合生成
  done
done

# === Step 2: 生成POTCAR ===
# POTCAR元素顺序必须与POSCAR一致!
# 例如: *NO在Cu slab → POTCAR = Cu + N + O
for sp in $SPECIES; do
  for site in $SITES; do
    dir="phase1_site_test/${sp}_${site}"
    # 根据物种确定需要的元素
    case $sp in
      H)    cat Cu/POTCAR H/POTCAR > $dir/POTCAR ;;
      O)    cat Cu/POTCAR O/POTCAR > $dir/POTCAR ;;
      N)    cat Cu/POTCAR N/POTCAR > $dir/POTCAR ;;
      OH)   cat Cu/POTCAR O/POTCAR H/POTCAR > $dir/POTCAR ;;
      NO)   cat Cu/POTCAR N/POTCAR O/POTCAR > $dir/POTCAR ;;
      NOH)  cat Cu/POTCAR N/POTCAR O/POTCAR H/POTCAR > $dir/POTCAR ;;
      NH)   cat Cu/POTCAR N/POTCAR H/POTCAR > $dir/POTCAR ;;
      NH2)  cat Cu/POTCAR N/POTCAR H/POTCAR > $dir/POTCAR ;;
      NH3)  cat Cu/POTCAR N/POTCAR H/POTCAR > $dir/POTCAR ;;
      NO2)  cat Cu/POTCAR N/POTCAR O/POTCAR > $dir/POTCAR ;;
      NO3)  cat Cu/POTCAR N/POTCAR O/POTCAR > $dir/POTCAR ;;
      NHOH) cat Cu/POTCAR N/POTCAR O/POTCAR H/POTCAR > $dir/POTCAR ;;
      N2O)  cat Cu/POTCAR N/POTCAR O/POTCAR > $dir/POTCAR ;;
    esac
  done
done

# === Step 3: 批量提交 ===
for sp in $SPECIES; do
  for site in $SITES; do
    dir="phase1_site_test/${sp}_${site}"
    cd $dir
    sbatch run.sh
    cd -
  done
done

# === Step 4: 监控任务 ===
# 检查哪些任务完成
for sp in $SPECIES; do
  for site in $SITES; do
    dir="phase1_site_test/${sp}_${site}"
    if grep -q "reached required accuracy" $dir/OUTCAR 2>/dev/null; then
      echo "DONE: $sp @ $site"
    else
      echo "RUNNING/FAILED: $sp @ $site"
    fi
  done
done
```

### 3.4 分析位点测试结果

```bash
# 提取所有能量
echo "species,site,energy_eV,converged" > summary_site_test.csv
for sp in $SPECIES; do
  for site in $SITES; do
    dir="phase1_site_test/${sp}_${site}"
    E=$(grep "energy(sigma->0)" $dir/OUTCAR 2>/dev/null | tail -1 | awk '{print $NF}')
    CONV=$(grep -c "reached required accuracy" $dir/OUTCAR 2>/dev/null)
    echo "$sp,$site,${E:-NA},${CONV:-0}" >> summary_site_test.csv
  done
done

# 找每种物种的最低能量位点
python3 << 'PYEOF'
import csv
from collections import defaultdict

data = defaultdict(list)
with open('summary_site_test.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['energy_eV'] != 'NA' and row['converged'] == '1':
            data[row['species']].append((row['site'], float(row['energy_eV'])))

print(f"{'物种':<8} {'最优位点':<10} {'E(eV)':<14} {'次优位点':<10} {'ΔE(meV)':<10}")
print("-" * 55)
for sp, entries in sorted(data.items()):
    entries.sort(key=lambda x: x[1])  # 按能量排序（最低最稳定）
    best_site, best_E = entries[0]
    if len(entries) > 1:
        second_site, second_E = entries[1]
        dE = (second_E - best_E) * 1000  # meV
        print(f"{sp:<8} {best_site:<10} {best_E:<14.6f} {second_site:<10} {dE:<10.1f}")
    else:
        print(f"{sp:<8} {best_site:<10} {best_E:<14.6f}")
PYEOF
```

**⚠️ 必查项**:
1. 弛豫后吸附物是否仍在初始位点？（可能迁移到更稳定位点）
   ```bash
   # 比较初始POSCAR和最终CONTCAR中吸附物的xy坐标
   # 如果偏移 > 0.5 Å, 说明发生了迁移
   ```
2. *NO的磁矩: `grep "magnetization (x)" OUTCAR | tail -1`，应接近0（吸附后淬灭）
3. 结构是否合理: 可视化CONTCAR (用VESTA或ASE的`view`)

---

## 4. Phase 2: 精确单吸附物

### 4.1 构建 p(4×4) slab

```python
# 修改gen_poscar_slab.py中的size参数
size = (4, 4)   # p(4×4), 64 Cu atoms (4 layers)
# KPOINTS改为 2×2×1
```

**KPOINTS** (phase2用):
```
K-Points
0
Gamma
2 2 1
0 0 0
```

### 4.2 工作流

```bash
# 1. 先优化清洁slab
cd phase2_single_ads/clean_slab_4x4
# (POSCAR: p(4×4) clean slab, INCAR同上但KPOINTS改2×2×1)
sbatch run.sh

# 2. 等清洁slab完成后，用CONTCAR为基础放置吸附物
# 只在 Phase 1 确定的最优位点放置, 共7个关键物种
BEST_SITES="NO:fcc H:fcc OH:fcc O:fcc N:fcc NH:fcc NH2:bridge"
# ↑ 根据Phase 1结果填写!

# 3. 对每个物种生成p(4×4)的POSCAR并提交

# 4. 提取精确能量
echo "species,E_slab_plus_A" > summary_single_ads.csv
E_clean=$(grep "energy(sigma->0)" clean_slab_4x4/OUTCAR | tail -1 | awk '{print $NF}')
echo "clean,$E_clean" >> summary_single_ads.csv
for entry in NO H OH O N NH NH2; do
  dir="${entry}_*_4x4"
  E=$(grep "energy(sigma->0)" $dir/OUTCAR | tail -1 | awk '{print $NF}')
  echo "$entry,$E" >> summary_single_ads.csv
done

# 5. 计算吸附能 (参考用, 非ε所需)
# BE(A) = E(slab+A) - E(clean) - E(A_gas)
```

### 4.3 镜像检查

**关键验证**: 比较p(3×3)和p(4×4)的单吸附物能量差

```python
# 镜像误差 = [E(3x3,slab+A) - E(3x3,clean)] - [E(4x4,slab+A) - E(4x4,clean)]
# 如果 |误差| < 20 meV, p(3×3)的位点测试结果可信
# 如果 |误差| > 50 meV, 说明p(3×3)位点测试有较大镜像干扰
```

---

## 5. Phase 3-4: 成对横向相互作用

### 5.1 共吸附结构的构建

**核心**: 在p(4×4) slab上放置两个吸附物，控制它们的间距

**Cu(111) p(4×4) fcc hollow位点编号** (俯视图):

```
行4:  ×13  ×14  ×15  ×16
行3:  ×9   ×10  ×11  ×12
行2:  ×5   ×6   ×7   ×8
行1:  ×1   ×2   ×3   ×4

位点间距离关系 (a_nn = a₀/√2):
  ×1 → ×2 : 1NN = a_nn = 2.566 Å
  ×1 → ×6 : 2NN = a_nn×√3 = 4.443 Å
  ×1 → ×3 : 3NN = 2×a_nn = 5.131 Å
```

**用ASE构建共吸附结构**:

```python
# scripts/gen_poscar_coads.py
from ase.io import read, write
from ase import Atoms
import numpy as np

slab = read('phase2_single_ads/clean_slab_4x4/CONTCAR')
a_nn = 3.63 / np.sqrt(2)  # 用你的a₀

# --- *NO + *NO at 1NN (fcc1 + fcc2) ---
# fcc1的Cartesian坐标 (根据slab几何计算)
# fcc2 = fcc1 + (a_nn, 0, 0) 沿a1方向
z_top = max(slab.positions[:, 2])  # 顶层Cu的z

# 位点1: fcc hollow #1
pos1_xy = np.array([a_nn/3, a_nn*np.sqrt(3)/3])  # 分数→Cartesian
# 位点2: fcc hollow #2 (1NN)
pos2_xy = pos1_xy + np.array([a_nn, 0.0])

# 放置 NO #1 (N朝下)
h_N = z_top + 1.3  # N距表面1.3 Å
h_O = h_N + 1.16   # N-O键长
slab.append(Atoms('N', positions=[[pos1_xy[0], pos1_xy[1], h_N]]))
slab.append(Atoms('O', positions=[[pos1_xy[0], pos1_xy[1], h_O]]))

# 放置 NO #2
slab.append(Atoms('N', positions=[[pos2_xy[0], pos2_xy[1], h_N]]))
slab.append(Atoms('O', positions=[[pos2_xy[0], pos2_xy[1], h_O]]))

write('phase3_tier1_pairs/NO_NO_1NN/POSCAR', slab, format='vasp', vasp5=True)
print(f"NO-NO distance: {np.linalg.norm(pos2_xy - pos1_xy):.3f} Å (should be {a_nn:.3f})")
```

### 5.2 ε的计算流程

完成Phase 2和Phase 3的所有计算后：

```python
# scripts/calc_epsilon.py
"""从DFT总能量计算成对横向相互作用ε"""

import json

# === 从VASP输出收集的能量 (eV) ===
# Phase 2 精确值
E_clean = -XXX.XXX  # ← 填入 phase2_single_ads/clean_slab_4x4
E_NO    = -XXX.XXX  # ← 填入 phase2_single_ads/NO_fcc_4x4
E_H     = -XXX.XXX  # ← 填入 phase2_single_ads/H_fcc_4x4
E_OH    = -XXX.XXX  # ← 填入 phase2_single_ads/OH_fcc_4x4
E_O     = -XXX.XXX  # ← 填入 phase2_single_ads/O_fcc_4x4
E_N     = -XXX.XXX  # ← 填入 phase2_single_ads/N_fcc_4x4

# Phase 3 共吸附能量
E_NO_NO_1NN = -XXX.XXX
E_NO_NO_2NN = -XXX.XXX
E_NO_H_1NN  = -XXX.XXX
E_NO_H_2NN  = -XXX.XXX
E_H_H_1NN   = -XXX.XXX
E_H_H_2NN   = -XXX.XXX
E_OH_NO_1NN = -XXX.XXX
E_OH_NO_2NN = -XXX.XXX

# === 计算ε ===
# ε(A,B) = E(slab+A+B) - E(slab+A) - E(slab+B) + E(clean)

results = {}

# *NO + *NO (同种, 简化公式)
eps_NO_NO_1NN = E_NO_NO_1NN - 2*E_NO + E_clean
eps_NO_NO_2NN = E_NO_NO_2NN - 2*E_NO + E_clean
results['NO-NO'] = {'1NN': eps_NO_NO_1NN, '2NN': eps_NO_NO_2NN}

# *NO + *H (异种)
eps_NO_H_1NN = E_NO_H_1NN - E_NO - E_H + E_clean
eps_NO_H_2NN = E_NO_H_2NN - E_NO - E_H + E_clean
results['NO-H'] = {'1NN': eps_NO_H_1NN, '2NN': eps_NO_H_2NN}

# *H + *H (同种)
eps_H_H_1NN = E_H_H_1NN - 2*E_H + E_clean
eps_H_H_2NN = E_H_H_2NN - 2*E_H + E_clean
results['H-H'] = {'1NN': eps_H_H_1NN, '2NN': eps_H_H_2NN}

# *OH + *NO (异种)
eps_OH_NO_1NN = E_OH_NO_1NN - E_OH - E_NO + E_clean
eps_OH_NO_2NN = E_OH_NO_2NN - E_OH - E_NO + E_clean
results['OH-NO'] = {'1NN': eps_OH_NO_1NN, '2NN': eps_OH_NO_2NN}

# === 输出 ===
print(f"{'Pair':<10} {'ε_1NN (eV)':<14} {'ε_2NN (eV)':<14} {'性质':<8} {'2NN/1NN'}")
print("-" * 60)
for pair, eps in results.items():
    sign = '排斥' if eps['1NN'] > 0 else '吸引'
    ratio = abs(eps['2NN']/eps['1NN'])*100 if abs(eps['1NN']) > 0.001 else 'N/A'
    print(f"{pair:<10} {eps['1NN']:+.4f}       {eps['2NN']:+.4f}       {sign:<8} "
          f"{ratio:.0f}%" if isinstance(ratio, float) else f"{ratio}")

# === 物理合理性检查 ===
for pair, eps in results.items():
    if eps['1NN'] < -0.10:
        print(f"⚠️ {pair} 1NN过强吸引({eps['1NN']:.3f} eV), 检查结构!")
    if abs(eps['2NN']) > abs(eps['1NN']):
        print(f"⚠️ {pair} |ε_2NN| > |ε_1NN|, 不物理!")
    if abs(eps['1NN']) > 0.50:
        print(f"⚠️ {pair} |ε_1NN|过大({eps['1NN']:.3f} eV), 可能有化学反应发生!")

# === 保存为KMC输入格式 ===
with open('results/lateral_interactions.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\n已保存到 results/lateral_interactions.json")
```

### 5.3 共吸附计算的关键检查

每个共吸附计算完成后**必须检查**:

```bash
# 1. 结构检查: 吸附物没有迁移或解离
# 比较POSCAR和CONTCAR中吸附物坐标
python3 -c "
from ase.io import read
ini = read('POSCAR')
fin = read('CONTCAR')
# 比较最后几个原子(吸附物)的坐标
n_ads = len(ini) - 36  # 36=Cu原子数(p3x3, 4层); p4x4为64
for i in range(-n_ads, 0):
    d = fin.positions[i] - ini.positions[i]
    print(f'Atom {i}: Δx={d[0]:.3f} Δy={d[1]:.3f} Δz={d[2]:.3f} |Δ|={sum(d**2)**0.5:.3f} Å')
"

# 2. 磁矩检查 (含NO的计算)
grep "magnetization (x)" OUTCAR | tail -1

# 3. 两吸附物间距检查
python3 -c "
from ase.io import read
atoms = read('CONTCAR')
# 找吸附物原子(如两个N原子)
n_atoms = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == 'N']
if len(n_atoms) >= 2:
    d = atoms.get_distance(n_atoms[0], n_atoms[1], mic=True)
    print(f'N-N distance: {d:.3f} Å')
"

# 4. 检查力收敛
grep "EDIFFG" OUTCAR
grep "FORCES:" OUTCAR | tail -5
# 最大力应 < EDIFFG的绝对值 (0.02 eV/Å)
```

---

## 6. Phase TS: 过渡态搜索

### 6.1 CI-NEB 方法 (推荐用于热力学步骤)

**适用**: *NO+*NO→*N₂O, *H+*H→H₂, *NO扩散, *H扩散

**工作流**:

```
1. 优化初态(IS)和末态(FS) → 获得CONTCAR
2. 用VTST工具生成中间images
3. 运行CI-NEB
4. 从saddle point提取Ea
```

**Step 1: 准备初态和末态**

以 *NO+*NO→*N₂O+* 为例:
- **初态(IS)**: 两个*NO共吸附在相邻fcc位点 (已有, 就是Phase 3的结构)
- **末态(FS)**: *N₂O吸附 + 空位 (需要单独优化)

```bash
# IS: 复制Phase 3的共吸附CONTCAR
cp phase3_tier1_pairs/NO_NO_1NN/CONTCAR phase_TS/NEB_NO_NO_coupling/IS/POSCAR

# FS: 需要构建*N₂O结构
# N₂O在atop位点, 端N朝下, N-N-O线性/弯曲
# 手动构建或从文献参考
```

**Step 2: 生成中间images**

```bash
# 需要VTST scripts (vtst.eng.ucsd.edu)
# 安装: 下载nebmake.pl等脚本

cd phase_TS/NEB_NO_NO_coupling

# 线性插值生成5个中间image
nebmake.pl IS/POSCAR FS/POSCAR 5
# 生成: 00/ (IS), 01/ 02/ 03/ 04/ 05/ (images), 06/ (FS)

# 检查插值合理性
# 可视化每个image的POSCAR, 确保原子没有重叠或不物理的位置
```

**Step 3: NEB的INCAR**

```fortran
SYSTEM  = NEB_NO_NO_coupling
PREC    = Accurate
ENCUT   = 450
EDIFF   = 1E-5
EDIFFG  = -0.03       ! NEB可以稍宽松
NSW     = 300
IBRION  = 3            ! ← VTST的quasi-Newton
POTIM   = 0            ! ← 必须! VTST自己控制步长
ISIF    = 2

! === NEB 特定参数 (需VTST编译的VASP) ===
ICHAIN  = 0            ! NEB方法
IMAGES  = 5            ! 中间image数
SPRING  = -5           ! 弹性常数 (eV/Å²)
LCLIMB  = .TRUE.       ! CI-NEB (收敛后自动开启climbing image)
IOPT    = 1            ! 优化器: 1=LBFGS (推荐), 2=CG

ISMEAR  = 1
SIGMA   = 0.15
ALGO    = Fast
ISPIN   = 2
LDIPOL  = .TRUE.
IDIPOL  = 3
LREAL   = Auto
LWAVE   = .FALSE.
LCHARG  = .FALSE.
```

> **重要**: CI-NEB需要用VTST编译的VASP! 检查: `grep VTST OUTCAR` 应有输出

**Step 4: 提取活化能**

```bash
# VTST工具提取MEP (minimum energy path)
nebef.pl    # 输出每个image的能量和力

# 或手动:
for i in 00 01 02 03 04 05 06; do
  E=$(grep "energy(sigma->0)" $i/OUTCAR | tail -1 | awk '{print $NF}')
  echo "$i  $E"
done

# Ea_fwd = E(saddle) - E(IS)
# Ea_rev = E(saddle) - E(FS)
# ΔH = E(FS) - E(IS)
```

---

### 6.2 PCET步骤的活化能处理

**问题**: 电化学步骤涉及质子-电子转移, 不能用标准NEB

**推荐方法: CHE (Computational Hydrogen Electrode)**

```
Ea(U) = Ea(U=0) + e × β × U

其中:
- Ea(U=0): 在U=0V vs SHE的活化能 (从NEB计算或文献获取)
- β: 对称因子, 通常 ≈ 0.5
- U: 施加电位 (V vs RHE)
- e: 电子电荷 = 1 (以eV为单位时)
```

**实际操作**:

```python
# 对于PCET步骤, 可以用两种方法获取Ea(U=0):

# 方法A: 从反应自由能用BEP关系估算 (最简单)
# Ea = alpha * dG + Ea0
# 典型参数: alpha=0.5, Ea0=0.5 eV (需校准)

# 方法B: 做表面加氢的NEB (不含电化学, 纯热力学)
# 将PCET步骤拆分为:
#   *NO + *H → *NOH + *     (热力学, 可以做NEB)
# 然后用CHE修正:
#   Ea(U) = Ea(thermal) + β×e×U

# 方法C: 用文献值 (最快)
# 从 ACS Catal. 9, 7052 (2019) 或 J. Phys. Chem. C (2025) 提取
```

**文献Ea估计值** (可作为初始参数, 后续替换为自己的DFT):

```python
# NO3RR PCET步骤活化能参考 (单位 eV, U=0 vs RHE)
PCET_barriers = {
    'NO3_to_NO2':  0.40,   # *NO₃ + H⁺+e⁻ → *NO₂ + *OH
    'NO2_to_NO':   0.35,   # *NO₂ + H⁺+e⁻ → *NO + *OH
    'NO_to_NOH':   0.70,   # *NO + H⁺+e⁻ → *NOH  ← PLS!
    'NOH_to_NHOH': 0.30,   # *NOH + H⁺+e⁻ → *NHOH
    'NHOH_to_NH':  0.30,   # *NHOH + H⁺+e⁻ → *NH + H₂O
    'NH_to_NH2':   0.30,   # *NH + H⁺+e⁻ → *NH₂
    'NH2_to_NH3':  0.20,   # *NH₂ + H⁺+e⁻ → *NH₃
    'OH_removal':  0.10,   # *OH + H⁺+e⁻ → H₂O + *
    'Volmer':      0.70,   # H₂O + * + e⁻ → *H + OH⁻  (碱性)
    # 热力学步骤 (需自己NEB)
    'NO_NO_coupling': None,  # *NO+*NO→*N₂O (估计 0.5-0.8 eV)
    'Tafel':          None,  # *H+*H→H₂ (估计 0.8-1.0 eV)
    'NO_diffusion':   None,  # *NO跳跃 (估计 0.1-0.3 eV)
    'H_diffusion':    None,  # *H跳跃 (估计 0.1-0.2 eV)
}

# 电位依赖: Ea(U) = Ea(U=0) + 0.5 × U  (β=0.5)
# 例: *NO→*NOH, U=-0.6V: Ea = 0.70 + 0.5×(-0.6) = 0.40 eV
```

---

### 6.3 Dimer方法 (CI-NEB的替代)

当初态/末态不明确或NEB不收敛时使用:

```fortran
! Dimer INCAR (需VTST编译的VASP)
IBRION = 3
POTIM  = 0
ICHAIN = 2          ! Dimer方法
IOPT   = 2          ! CG优化器
DdR    = 0.005      ! Dimer有限差分步长
DRotMax = 4         ! 最大旋转步数
```

**操作**:
```bash
# 1. 从NEB的最高image开始 (近似TS)
# 2. 或手动构建近似TS结构
# 3. 需要提供MODECAR文件 (位移方向, 可从NEB获得)
# 4. Dimer会沿虚频方向搜索saddle point
```

---

### 6.4 扩散势垒 (简单NEB)

**\*NO在fcc→fcc跳跃** (经过bridge位点):

```
IS: *NO @ fcc site #1
TS: *NO @ bridge (过渡态)
FS: *NO @ fcc site #2 (相邻)
```

```bash
# 只需3个images (路径短)
nebmake.pl IS/POSCAR FS/POSCAR 3

# INCAR中 IMAGES = 3
# 典型结果: Ea_diff(NO) ≈ 0.1-0.3 eV
#           Ea_diff(H)  ≈ 0.1-0.2 eV
```

---

## 7. 数据提取与后处理脚本

### 7.1 批量能量提取

```bash
#!/bin/bash
# scripts/extract_energy.sh
# 用法: ./extract_energy.sh <directory>

echo "Directory,E_sigma0(eV),E_without_entropy(eV),Converged,NSW_used,Magnetization"

find "$1" -name OUTCAR -type f | sort | while read outcar; do
  dir=$(dirname "$outcar")
  E0=$(grep "energy(sigma->0)" "$outcar" | tail -1 | awk '{print $NF}')
  Ewe=$(grep "energy  without entropy" "$outcar" | tail -1 | awk '{print $4}')
  conv=$(grep -c "reached required accuracy" "$outcar")
  nsw=$(grep -c "LOOP+" "$outcar")
  mag=$(grep "magnetization (x)" "$outcar" | tail -1 | awk '{print $NF}')
  echo "${dir},$E0,$Ewe,$conv,$nsw,${mag:-0}"
done
```

### 7.2 完整KMC参数整合

```python
# scripts/compile_kmc_params.py
"""将所有DFT结果整合为SPARK可用的参数文件"""

import json

# === 从DFT提取的参数 ===

params = {
    "system": "NO3RR on Cu(111)",
    "lattice_constant_A": 3.63,       # ← Phase 0
    "nearest_neighbor_A": 2.566,      # a₀/√2
    "temperature_K": 298,

    "species": [
        "empty", "NO3", "NO2", "NO", "NOH", "NHOH",
        "NH", "NH2", "NH3", "N2O", "OH", "H"
    ],

    "adsorption_energies_eV": {
        # BE = E(slab+A) - E(clean) - E(A_gas)
        "NO":  None,   # ← 填入
        "H":   None,
        "OH":  None,
        "O":   None,
        "N":   None,
        "NH":  None,
        "NH2": None,
        "NH3": None,   # 预期 ~ -0.37 eV (弱)
        "NO3": None,
        "NO2": None,
        "N2O": None,   # 预期弱吸附
    },

    "preferred_sites": {
        # ← 从Phase 1填入
        "NO": "fcc",
        "H":  "fcc",
        "OH": "fcc",   # or bridge
        "O":  "fcc",
        "N":  "fcc",
        "NH": "fcc",
        "NH2": "bridge",
        "NH3": "atop",
        "NO3": "bidentate",
        "NO2": "bidentate",
    },

    "lateral_interactions_eV": {
        # ε(A,B,distance) — 从Phase 3-4填入
        "NO-NO": {"1NN": None, "2NN": None},
        "NO-H":  {"1NN": None, "2NN": None},
        "H-H":   {"1NN": None, "2NN": None},
        "OH-NO": {"1NN": None, "2NN": None},
        "OH-OH": {"1NN": None, "2NN": None},
        "OH-H":  {"1NN": None, "2NN": None},
        "N-NO":  {"1NN": None, "2NN": None},
        "N-N":   {"1NN": None, "2NN": None},
        "O-NO":  {"1NN": None, "2NN": None},
        "N-H":   {"1NN": None, "2NN": None},
    },

    "activation_energies_eV": {
        # 热力学步骤 — 从NEB填入
        "NO_NO_coupling": {"Ea_fwd": None, "Ea_rev": None},
        "Tafel":          {"Ea_fwd": None, "Ea_rev": None},
        "NO_diffusion":   {"Ea_fwd": None},
        "H_diffusion":    {"Ea_fwd": None},
        # PCET步骤 — Ea(U) = Ea0 + β×U
        "NO3_to_NO2":  {"Ea0": 0.40, "beta": 0.5, "n_electron": 1},
        "NO2_to_NO":   {"Ea0": 0.35, "beta": 0.5, "n_electron": 1},
        "NO_to_NOH":   {"Ea0": 0.70, "beta": 0.5, "n_electron": 1},
        "NOH_to_NHOH": {"Ea0": 0.30, "beta": 0.5, "n_electron": 1},
        "NHOH_to_NH":  {"Ea0": 0.30, "beta": 0.5, "n_electron": 1},
        "NH_to_NH2":   {"Ea0": 0.30, "beta": 0.5, "n_electron": 1},
        "NH2_to_NH3":  {"Ea0": 0.20, "beta": 0.5, "n_electron": 1},
        "OH_removal":  {"Ea0": 0.10, "beta": 0.5, "n_electron": 1},
        "Volmer":      {"Ea0": 0.70, "beta": 0.5, "n_electron": 1},
    },

    "desorption_energies_eV": {
        "NH3": 0.37,     # |BE(NH₃)|
        "N2O": 0.30,     # |BE(N₂O)|
    },

    "product_electrons": {
        "NH3": 8,   # NO₃⁻ + 9H⁺ + 8e⁻ → NH₃ + 3H₂O
        "N2":  10,  # 2NO₃⁻ + 12H⁺ + 10e⁻ → N₂ + 6H₂O (per N₂)
        "N2O": 8,   # 2NO₃⁻ + 10H⁺ + 8e⁻ → N₂O + 5H₂O (per N₂O)
        "H2":  2,   # 2H⁺ + 2e⁻ → H₂
    },

    "references": [
        "Phase 0-4: this work (PBE/PAW, 450 eV, p(4x4) Cu(111) 4-layer slab)",
        "PCET Ea0: ACS Catal. 9, 7052 (2019); J. Phys. Chem. C (2025)",
    ]
}

with open('results/kmc_params.json', 'w') as f:
    json.dump(params, f, indent=2, ensure_ascii=False)

print("KMC参数文件已生成: results/kmc_params.json")
print("需要填入的None值:", sum(1 for v in str(params).split('null')))
```

---

## 8. 收敛性检查清单

### 8.1 每个VASP计算完成后的必查项

```bash
# scripts/check_convergence.sh
#!/bin/bash
# 用法: ./check_convergence.sh <OUTCAR_directory>

DIR=${1:-.}
OUTCAR="$DIR/OUTCAR"

echo "=== 收敛性检查: $DIR ==="

# 1. 电子收敛
echo -n "电子收敛: "
if grep -q "reached required accuracy" "$OUTCAR"; then
  echo "✅ PASS"
else
  echo "❌ FAIL - 未达到收敛精度"
  echo "  最后EDIFF: $(grep 'EDIFF ' $OUTCAR | head -1)"
  echo "  最后变化: $(grep 'total energy' $OUTCAR | tail -3)"
fi

# 2. 离子收敛 (力)
echo -n "力收敛: "
NSW_ACTUAL=$(grep -c "LOOP+" "$OUTCAR")
NSW_MAX=$(grep "NSW" "$DIR/INCAR" | awk -F= '{print $2}' | tr -d ' ')
if [ "$NSW_ACTUAL" -lt "${NSW_MAX:-300}" ]; then
  echo "✅ PASS (步数: $NSW_ACTUAL)"
else
  echo "⚠️ 可能未收敛 (用完$NSW_ACTUAL步)"
fi

# 3. 最大残余力
echo -n "最大残余力: "
FMAX=$(grep "FORCES:" "$OUTCAR" | tail -1 | awk '{print $NF}')
echo "${FMAX:-N/A} eV/Å (目标 < 0.02)"

# 4. 能量变化
echo "最后5步能量变化:"
grep "energy(sigma->0)" "$OUTCAR" | tail -5 | awk '{print NR": "$NF" eV"}'

# 5. 磁矩 (如果ISPIN=2)
MAG=$(grep "magnetization (x)" "$OUTCAR" | tail -1 | awk '{print $NF}')
if [ -n "$MAG" ]; then
  echo "总磁矩: $MAG μ_B"
fi

# 6. SIGMA校正
echo -n "SIGMA校正: "
ENTROPY=$(grep "entropy T\*S" "$OUTCAR" | tail -1 | awk '{print $NF}')
echo "${ENTROPY:-N/A} eV (应 < 1 meV/atom)"
```

### 8.2 全局收敛性验证 (Phase完成后)

| 检查项 | 方法 | 标准 |
|-------|------|------|
| 晶格常数 | Phase 0 a₀ vs 实验 3.615 Å | 偏差 < 2% |
| k-points收敛 | p(3×3): 比较3×3×1和5×5×1 | ΔE < 10 meV |
| slab层数 | 比较4层和5层slab E(clean) | ΔE < 5 meV/atom |
| 超胞大小 | 比较p(3×3)和p(4×4)单吸附 | ΔE_ads < 30 meV |
| 位点稳定性 | CONTCAR vs POSCAR吸附物坐标 | 偏移 < 0.5 Å |
| 对称性 | 等价位点的E差异 | < 5 meV |
| SIGMA | entropy T*S < 1 meV/atom | 所有计算 |

---

## 9. 常见问题排查

### 9.1 计算不收敛

| 症状 | 可能原因 | 解决方案 |
|------|---------|---------|
| 电子不收敛 (NELM耗尽) | 初始结构太差 / SIGMA太小 | 增大SIGMA到0.2; 用ALGO=All; 增大NELM |
| 离子步达到NSW上限 | 初始力太大 | 检查初始结构; 减小POTIM; 改IBRION=1(RMM-DIIS) |
| 能量振荡 | IBRION=2(CG)困在鞍点 | 改IBRION=1; 或加小随机扰动 |
| 磁矩异常 | 自旋态错误 | 设MAGMOM初始值; 检查ISPIN=2 |

### 9.2 吸附物迁移/解离

| 症状 | 原因 | 解决方案 |
|------|------|---------|
| 吸附物从bridge跑到hollow | bridge不是最优位点 | 正常! 记录结果 |
| NO解离为N+O | N-O键断裂, 能量更低 | 检查是否正确; 可能需要约束N-O键 |
| NH₃脱附离开表面 | 结合太弱 | 减小EDIFFG; 或用更大SIGMA |
| 两个共吸附物合并 | 初始距离太近 | 确认1NN距离正确; 用p(4×4) |

### 9.3 POTCAR元素顺序

**⚠️ VASP最常见错误!**

POTCAR中元素的顺序**必须**与POSCAR第6行(元素行)一致:
```
POSCAR:  Cu  N  O  H    ← 第6行
POTCAR:  Cu + N + O + H  ← 必须同顺序cat
```

验证:
```bash
# 检查POTCAR中的元素
grep "TITEL" POTCAR
# 应该按POSCAR中的顺序出现
```

---

## 10. 附录: 坐标生成参考

### 10.1 Cu(111) fcc hollow位点的Cartesian坐标

对于p(3×3), 晶格常数a₀, a_nn = a₀/√2:

```python
import numpy as np

a0 = 3.63  # 你的晶格常数
a_nn = a0 / np.sqrt(2)

# p(3×3)晶胞向量
a1 = np.array([3*a_nn, 0, 0])
a2 = np.array([1.5*a_nn, 1.5*a_nn*np.sqrt(3), 0])

# fcc hollow位点 (分数坐标)
# 在单胞中: (1/3, 1/3)
# 在p(3×3)中有9个等价fcc位点:
fcc_sites = []
for i in range(3):
    for j in range(3):
        frac = np.array([(i + 1/3)/3, (j + 1/3)/3])
        cart = frac[0]*a1 + frac[1]*a2
        fcc_sites.append(cart)
        print(f"fcc({i},{j}): x={cart[0]:.3f}, y={cart[1]:.3f}")

# 相邻fcc位点间距
for i, s1 in enumerate(fcc_sites):
    for j, s2 in enumerate(fcc_sites):
        if i < j:
            d = np.linalg.norm(s1 - s2)
            if d < a_nn * 1.1:
                print(f"1NN: fcc{i}-fcc{j}, d={d:.3f} Å")
            elif d < a_nn * 1.8:
                print(f"2NN: fcc{i}-fcc{j}, d={d:.3f} Å")
```

### 10.2 Selective Dynamics 标记生成

```python
from ase.io import read, write
import numpy as np

slab = read('POSCAR')
z = slab.positions[:, 2]
z_layers = np.sort(np.unique(np.round(z, 1)))

# 底2层固定, 顶2层弛豫
z_cut = (z_layers[1] + z_layers[2]) / 2

selective = []
for atom_z in z:
    if atom_z < z_cut:
        selective.append([False, False, False])  # F F F
    else:
        selective.append([True, True, True])      # T T T

# 写入带Selective Dynamics的POSCAR
# ASE方式:
from ase.constraints import FixAtoms
fix = FixAtoms(indices=[i for i, s in enumerate(selective) if not s[0]])
slab.set_constraint(fix)
write('POSCAR', slab, format='vasp', vasp5=True, direct=True)
```

### 10.3 典型 run.sh (SLURM)

```bash
#!/bin/bash
#SBATCH --job-name=Cu111_DFT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --output=slurm-%j.out

module load vasp/6.4.1   # 或你的版本
# 如果需要VTST (NEB):
# module load vasp/6.4.1-vtst

srun vasp_std
# 或 mpirun -np $SLURM_NTASKS vasp_std
```

---

*本文档为 `DFT_lateral_interactions_guide.md` 的操作手册补充*
*Last updated: 2026-03-23*
