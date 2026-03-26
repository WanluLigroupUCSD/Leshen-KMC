# SPARK Tutorial

SPARK (SPatial Atomistic Reaction Kinetics) 完整教程。

---

## 目录

1. [安装](#1-安装)
2. [SPARKIN 输入文件](#2-sparkin-输入文件)
3. [快速入门：CO 吸脱附](#3-快速入门co-吸脱附)
4. [完整示例：HER on Pt(111)](#4-完整示例her-on-pt111)
5. [运行 KMC 仿真](#5-运行-kmc-仿真)
6. [均场微观动力学模型](#6-均场微观动力学模型)
7. [极化曲线](#7-极化曲线)
8. [高级功能](#8-高级功能)
9. [Python API](#9-python-api)
10. [常见问题](#10-常见问题)

---

## 1. 安装

```bash
git clone https://github.com/WanluLigroupUCSD/SPARK.git
cd SPARK
pip install numpy scipy matplotlib pyyaml
```

验证安装：

```python
import spark
print(spark.__version__)  # 0.4.0
```

---

## 2. SPARKIN 输入文件

SPARK 使用类似 VASP INCAR 的 `TAG = VALUE` 格式。所有模型定义和仿真参数写在一个 `.sparkin` 文件中。

### 2.1 基本语法

```
# 这是注释（整行）
SYSTEM    = HER_Pt111           (这是行内注释)
TEMP      = 298                 (K)
```

规则：
- `TAG = VALUE`，等号两边空格可选
- `#` 开头为注释行
- `(...)` 为行内注释（类似 VASP），解析时自动忽略
- 标签名**不区分大小写**（推荐全大写）
- 空行自动忽略

### 2.2 所有标签

| 标签 | 用途 | 可否重复 | 示例 |
|------|------|----------|------|
| `SYSTEM` | 模型名称 | 否 | `SYSTEM = HER_Pt111` |
| `SURFACE` | 表面预设 | 否 | `SURFACE = Pt(111)` |
| `LATTICE` | 自定义晶格 (Å) | 否 | `LATTICE = 3.5 3.5 10.0` |
| `DIMENSION` | 维度 (1/2/3) | 否 | `DIMENSION = 2` |
| `TEMP` | 温度 (K) | 否 | `TEMP = 298` |
| `POTENTIAL` | 电位 (V vs RHE) | 否 | `POTENTIAL = -0.2` |
| `SPECIES` | 表面物种 | 否 | `SPECIES = H CO` |
| `PARAMETER` | 自定义参数 | 是 | `PARAMETER = p_CO 1.0` |
| `REACTION` | 基元反应 | 是 | 见下文 |
| `DIFFUSION` | 扩散 | 是 | `DIFFUSION = H : 0.10` |
| `LATERAL` | 侧向相互作用 | 是 | `LATERAL = H-H : 0.10` |
| `BEP` | BEP 关系 | 是 | `BEP = Volmer_fwd : 0.5` |
| `LSIZE` | 晶格大小 | 否 | `LSIZE = 20 20` |
| `NEQUIL` | 平衡步数 | 否 | `NEQUIL = 100000` |
| `NSTEPS` | 采样步数 | 否 | `NSTEPS = 100000` |

### 2.3 SURFACE 预设

`SURFACE` 自动配置晶格常数和位点类型，无需手动定义晶格：

| 预设 | 晶格常数 (Å) | 位点 | 预设 | 晶格常数 (Å) | 位点 |
|------|-------------|------|------|-------------|------|
| Pt(111) | 2.775 | top | Cu(111) | 2.556 | top |
| Pt(100) | 2.775 | hollow | Cu(100) | 2.556 | hollow |
| Pd(111) | 2.751 | top | Au(111) | 2.884 | top |
| Pd(100) | 2.751 | hollow | Au(100) | 2.884 | hollow |
| Ni(111) | 2.492 | top | Ag(111) | 2.889 | top |
| Ru(0001) | 2.706 | top | Ir(111) | 2.715 | top |
| Mo(110) | 2.725 | top | Fe(110) | 2.482 | top |

如果你的表面不在预设中，使用 `LATTICE`：

```
LATTICE = 2.55 2.55 15.0
```

### 2.4 SPECIES

声明表面吸附物种。**空位 `*` 始终隐式存在**，不需要声明。

```
SPECIES = H          (单物种)
SPECIES = H CO OH    (多物种，空格分隔)
```

### 2.5 REACTION

反应定义是 SPARKIN 的核心。格式为冒号 `:` 分隔的字段：

```
REACTION = 名称 : 方程式 : Ea(eV) [: 类型 : beta]
```

#### 热力学反应（默认）

速率自动生成为 TST：k = k_B T/h × exp(-Ea / k_B T)

```
REACTION = Tafel : H + H -> * + * : 0.85
```

含义：H\* + H\*(邻) → H₂(g) + 2\*，Ea = 0.85 eV

#### 电化学反应（PCET）

速率：k = k_B T/h × exp(-(Ea + β×U) / k_B T)

```
REACTION = Volmer_fwd : * -> H : 0.67 : electrochemical : 0.5
```

含义：H⁺(aq) + e⁻ + \* → H\*，Ea = 0.67 eV，β = 0.5

β > 0 用于正向反应，β < 0 用于逆向反应：

```
REACTION = Volmer_rev : H -> * : 0.62 : electrochemical : -0.5
```

#### 自定义速率

用 `rate=` 覆盖自动生成的速率表达式：

```
REACTION = CO_ads : * -> CO : rate=p_CO*bar*A/sqrt(2*pi*m_CO*umass/(kB*T))
```

可用变量：`kB`, `h`, `T`, `U`, `eV`, `bar`, `angstrom`, `umass`, `pi`, `exp`, `sqrt`, `log`, 以及所有 `PARAMETER` 和分子质量 `m_CO`, `m_H2` 等。

#### 方程式语法

| 写法 | 含义 |
|------|------|
| `* -> H` | 单位点：空位变为 H |
| `H -> *` | 单位点：H 脱附 |
| `H + H -> * + *` | 双位点：`+` 后的第二物种默认在近邻位 |
| `H + H(nn) -> * + *(nn)` | 同上，`(nn)` 显式标注近邻 |
| `CO@hollow -> *@hollow` | 指定位点名称 |

### 2.6 DIFFUSION

自动展开为所有近邻方向的扩散过程（2D 生成 4 个，3D 生成 6 个）：

```
DIFFUSION = H : 0.10
```

等价于定义 4 个反应：
```
H_diff_right: H@(0,0,0) + *@(1,0,0) → *@(0,0,0) + H@(1,0,0)
H_diff_left:  H@(0,0,0) + *@(-1,0,0) → *@(0,0,0) + H@(-1,0,0)
H_diff_up:    H@(0,0,0) + *@(0,1,0) → *@(0,0,0) + H@(0,1,0)
H_diff_down:  H@(0,0,0) + *@(0,-1,0) → *@(0,0,0) + H@(0,-1,0)
```

### 2.7 LATERAL 和 BEP

```
LATERAL = H-H : 0.10     (eV, 正值=排斥, 负值=吸引)
LATERAL = CO-CO : -0.05

BEP = Volmer_fwd : 0.5   (反应名 : alpha)
```

**侧向相互作用**：近邻吸附物之间的相互作用能。影响含该吸附物的所有反应速率。排斥使脱附速率增大，吸引使脱附速率减小。

**BEP 关系**：Ea(env) = Ea(0) + α × ΔΔH，其中 ΔΔH 来自侧向相互作用引起的反应焓变化。

### 2.8 仿真控制

```
LSIZE  = 20 20       (晶格大小，20×20 = 400 个位点)
NEQUIL = 100000      (平衡步数)
NSTEPS = 100000      (采样步数)
```

这些参数存储在 `pt.meta` 中，Python 脚本可以读取：

```python
pt = load_model('model.sparkin')
size = pt.meta.get('lsize', [20, 20])
```

---

## 3. 快速入门：CO 吸脱附

最简单的 KMC 模型：CO 在 Pd(100) 上的吸附/脱附。

### 3.1 输入文件

创建 `co_model.sparkin`：

```
SYSTEM    = CO_on_Pd100
SURFACE   = Pd(100)
TEMP      = 600
SPECIES   = CO

PARAMETER = p_CO    1.0
PARAMETER = A       (3.5*angstrom)**2
PARAMETER = deltaG  -0.5

REACTION  = CO_ads : * -> CO : rate=p_CO*bar*A/sqrt(2*pi*m_CO*umass/(kB*T))
REACTION  = CO_des : CO -> * : rate=p_CO*bar*A/sqrt(2*pi*m_CO*umass/(kB*T))*exp(deltaG*eV/(kB*T))

LSIZE     = 20 20
NSTEPS    = 200000
```

### 3.2 运行

```python
from spark import load_model, KMCEngine

# 加载模型
pt = load_model('co_model.sparkin')
pt.summary()

# 创建引擎
size = pt.meta.get('lsize', [20, 20])
engine = KMCEngine(pt, size=size, print_rates=True)

# 平衡
engine.do_steps(100000)

# 重置 TOF 计数器
engine.get_tof()

# 采样
engine.do_steps(100000)

# 结果
print('Coverage:', engine.get_coverage())
print('TOF:', engine.get_tof())
```

### 3.3 Langmuir 等温线验证

对简单吸脱附体系，KMC 结果应与 Langmuir 解析解一致：

```python
import numpy as np
from spark import load_model, KMCEngine
from spark.units import kB, eV

pt = load_model('co_model.sparkin')

# 扫描 deltaG
deltaG_range = np.linspace(-1.0, 0.5, 16)
theta_kmc = []

for dG in deltaG_range:
    # 修改参数
    for p in pt.parameter_list:
        if p.name == 'deltaG':
            p.value = dG

    engine = KMCEngine(pt, size=[20, 20], print_rates=False, banner=False)
    engine.do_steps(200000)
    cov = engine.get_coverage()
    theta_kmc.append(cov.get('CO', 0.0))

# 解析解
T = 600.0
K_eq = np.exp(-deltaG_range * eV / (kB * T))
theta_langmuir = K_eq / (1 + K_eq)

# 比较
for dG, kmc, ana in zip(deltaG_range, theta_kmc, theta_langmuir):
    print(f'deltaG={dG:+.2f}  KMC={kmc:.3f}  Langmuir={ana:.3f}')
```

---

## 4. 完整示例：HER on Pt(111)

### 4.1 输入文件

创建 `her.sparkin`：

```
# ============================================================
#  HER on Pt(111) in acidic media
# ============================================================
#
# Volmer:    H+(aq) + e- + * <-> H*        (PCET)
# Tafel:     H* + H*(nn) -> H2(g) + 2*     (thermal)
# Heyrovsky: H* + H+(aq) + e- -> H2 + *    (PCET)
# Diffusion: H*@s1 + *@s2 <-> *@s1 + H*@s2
#
# DFT references:
#   Li et al., ACS Catal. 14, 2696 (2024)
#   Skulason et al., J. Phys. Chem. C 114, 18182 (2010)

SYSTEM    = HER_Pt111
SURFACE   = Pt(111)
TEMP      = 298                 (K)
POTENTIAL = -0.2                (V vs RHE)
SPECIES   = H

# Elementary reactions
REACTION  = Volmer_fwd : * -> H : 0.67 : electrochemical : 0.5
REACTION  = Volmer_rev : H -> * : 0.62 : electrochemical : -0.5
REACTION  = Tafel : H + H -> * + * : 0.85
REACTION  = Heyrovsky : H -> * : 0.70 : electrochemical : 0.5
DIFFUSION = H : 0.10

# Interactions
LATERAL   = H-H : 0.10         (eV, repulsive NN)
BEP       = Volmer_fwd : 0.5

# Simulation
LSIZE     = 20 20
NEQUIL    = 100000
NSTEPS    = 500000
```

### 4.2 运行仿真

```python
from spark import load_model, KMCEngine
from spark.analysis import TrajectoryRecorder

pt = load_model('her.sparkin')
size = pt.meta.get('lsize', [20, 20])
nequil = pt.meta.get('nequil', 100000)
nsteps = pt.meta.get('nsteps', 500000)

engine = KMCEngine(pt, size=size, print_rates=True)

# 平衡阶段
engine.do_steps(nequil)
engine.get_tof()  # 重置 TOF

# 采样阶段：每 5000 步记录一次
recorder = TrajectoryRecorder(engine)
n_samples = 100
steps_per_sample = nsteps // n_samples

for i in range(n_samples):
    engine.do_steps(steps_per_sample)
    recorder.record()

# 查看结果
print('Final coverage:', engine.get_coverage())
print('TOF:', engine.get_tof())
print('Process stats:', engine.get_process_stats())
```

### 4.3 绘制覆盖度随时间演化

```python
import matplotlib.pyplot as plt

times = recorder.get_times()
theta_H = recorder.get_coverage_array('H')

plt.figure(figsize=(8, 5))
plt.plot(times * 1e6, theta_H, 'b-')
plt.xlabel('Time (μs)')
plt.ylabel('θ_H')
plt.title('H coverage evolution on Pt(111)')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('coverage_evolution.png', dpi=150)
```

### 4.4 电位扫描（极化曲线）

```python
import numpy as np
from spark import load_model, KMCEngine

pt = load_model('her.sparkin')
U_range = np.linspace(-0.5, 0.0, 11)

results = {'U': [], 'theta_H': [], 'j_Heyrovsky': []}

for U in U_range:
    # 设置电位
    engine = KMCEngine(pt, size=[20, 20], print_rates=False, banner=False)
    engine.parameters.U = U

    # 平衡 + 采样
    engine.do_steps(200000)
    engine.get_tof()
    engine.do_steps(500000)

    cov = engine.get_coverage()
    tof = engine.get_tof()
    results['U'].append(U)
    results['theta_H'].append(cov.get('H', 0.0))
    results['j_Heyrovsky'].append(tof.get('Heyrovsky', 0.0))

    print(f'U={U:.3f} V  theta_H={cov.get("H",0):.4f}  '
          f'TOF_Hey={tof.get("Heyrovsky",0):.4e}')
```

---

## 5. 运行 KMC 仿真

### 5.1 KMCEngine 基本用法

```python
from spark import load_model, KMCEngine

pt = load_model('model.sparkin')
engine = KMCEngine(pt, size=[30, 30], print_rates=True)
```

**参数说明：**
- `size`：晶格大小，如 `[20, 20]` 为 400 位点
- `print_rates=True`：打印初始速率常数
- `banner=True`：打印模型信息横幅

### 5.2 运行步骤

```python
# 运行固定步数
engine.do_steps(100000)

# 单步运行
success = engine.do_kmc_step()  # 返回 True/False

# 查看状态
print(engine.kmc_time)   # 仿真时间 [s]
print(engine.kmc_step)   # 已执行步数
```

### 5.3 读取结果

```python
# 覆盖度
cov = engine.get_coverage()
# {'empty': 0.66, 'H': 0.34}

# 转化频率 (TOF)
tof = engine.get_tof()
# {'Volmer_fwd': 1234.5, 'Heyrovsky': 567.8}
# 注意：get_tof() 返回上次调用以来的平均值，并重置计数器

# 过程统计
stats = engine.get_process_stats()
# {'Volmer_fwd': 123456, 'Tafel': 789}

# 2D 晶格快照
lattice_2d = engine.get_lattice_2d()  # shape: (ny, nx)
```

### 5.4 动态修改参数

```python
# 改变温度（自动重新计算所有速率）
engine.parameters.T = 350

# 改变电位
engine.parameters.U = -0.3

# 查看所有参数
print(engine.parameters)
```

### 5.5 轨迹记录

```python
from spark.analysis import TrajectoryRecorder

recorder = TrajectoryRecorder(engine)

for _ in range(200):
    engine.do_steps(5000)
    recorder.record()

# 提取数据
times = recorder.get_times()           # shape: (200,)
theta_H = recorder.get_coverage_array('H')  # shape: (200,)
all_cov = recorder.get_all_coverages()      # dict of arrays
```

### 5.6 稳态检测

```python
from spark.analysis import run_to_steady_state

recorder = run_to_steady_state(
    engine,
    max_steps=1e7,
    sample_interval=10000,
    rtol=0.05,       # 相对容差 5%
    min_samples=100,
    verbose=True,
)
print('Steady-state coverage:', engine.get_coverage())
```

### 5.7 直接操作晶格

```python
# 在指定位点放置物种
engine.put((10, 10), 'H')

# 读取指定位点
sp = engine.get((10, 10))  # 'H' or 'empty'

# 重置为初始状态
engine.reset()

# 保存/恢复配置
config = engine.get_configuration()
engine.set_configuration(config)
```

---

## 6. 均场微观动力学模型

对于不需要空间分辨的快速计算，使用均场 ODE 求解器。

### 6.1 基本用法

```python
from spark import MicroKineticModel
from spark.rates import tst_rate, electrochemical_rate

mkm = MicroKineticModel()

# 添加物种（不需要添加 empty）
mkm.add_species('H')

# 添加反应
mkm.add_reaction(
    name='Volmer_fwd',
    reactants={},           # 空位吸附，reactants 为空
    products={'H': 1},
    rate_fwd=lambda p: electrochemical_rate(0.67, p['T'], p['U'], beta_bv=0.5),
    rate_rev=lambda p: electrochemical_rate(0.62, p['T'], p['U'], beta_bv=-0.5),
    tof_count={'Volmer': 1},
)

mkm.add_reaction(
    name='Heyrovsky',
    reactants={'H': 1},
    products={},
    rate_fwd=lambda p: electrochemical_rate(0.70, p['T'], p['U'], beta_bv=0.5),
    tof_count={'H2_production': 1},
)

# 设置参数
mkm.parameters = {'T': 298, 'U': -0.2}

# 求稳态
ss = mkm.solve_steady_state()
mkm.print_summary(ss)
```

### 6.2 瞬态 ODE 求解

```python
import numpy as np

# 求解时间演化
sol = mkm.solve_ode(
    t_span=(0, 1.0),              # 0 到 1 秒
    t_eval=np.logspace(-9, 0, 500),
    method='BDF',                 # 刚性 ODE 推荐
)

if sol.success:
    for i, sp in enumerate(mkm.species):
        print(f'{sp}: {sol.y[i, -1]:.6e}')
```

### 6.3 参数扫描

```python
import numpy as np

# 电位扫描
results = mkm.scan_parameter('U', np.linspace(-0.5, 0.0, 21),
                              observable='H2_production')

for U, tof, cov in zip(results['values'], results['tofs'], results['coverages']):
    print(f'U={U:.2f} V  TOF={tof:.4e}  theta_H={cov.get("H", 0):.4f}')
```

### 6.4 速率控制度分析

```python
# Campbell 速率控制度 (DRC)
ss = mkm.solve_steady_state()
drc = mkm.degree_of_rate_control(ss, observable='H2_production')

print('Degree of Rate Control:')
for rxn, xrc in sorted(drc.items(), key=lambda x: -abs(x[1])):
    print(f'  {rxn}: {xrc:+.4f}')
```

---

## 7. 极化曲线

从 KMC 或 MKM 结果计算电化学极化曲线 (j-U)。

### 7.1 从 MKM 计算

```python
import numpy as np
from spark.polarization import PolarizationCurve

# mkm 是已定义好的 MicroKineticModel

pc = PolarizationCurve(
    mkm=mkm,
    n_electrons={'H2_production': 2},  # 每个 H2 涉及 2 个电子
    A_site=(2.775e-10)**2,             # 位点面积 [m²]
)

results = pc.compute(
    U_range=np.linspace(-0.5, 0.0, 21),
    T=298,
    verbose=True,
)

pc.print_results()
pc.save('polarization.dat')
```

### 7.2 从 DFT 能量数据

如果有恒电位 DFT 数据（不同电位下的能垒）：

```python
from spark.polarization import EnergyLandscape, load_energy_data

# 加载 JSON 格式的 DFT 数据
data = load_energy_data('dft_energies.json')

landscape = EnergyLandscape(
    potentials=data['potentials'],
    state_energies=data['state_energies'],
    ts_energies=data['ts_energies'],
)

landscape.summary(U=-0.3)  # 打印 U=-0.3V 下的能垒
```

### 7.3 TOF 转电流密度

```python
from spark.polarization import tof_to_current_density

j = tof_to_current_density(
    tof=1e4,          # s⁻¹ per site
    n_electrons=2,    # H2: 2e⁻
    A_site=7.7e-20,   # m²
)
print(f'j = {j:.2f} mA/cm²')
```

---

## 8. 高级功能

### 8.1 侧向相互作用机制

添加 `LATERAL` 后，每个位点的反应速率不再是全局常数，而是取决于局部环境：

```
LATERAL = H-H : 0.10   (repulsive)
```

物理含义：每个 H 邻居贡献 +0.10 eV 的不稳定化能。对于 4 个近邻的 2D 晶格：
- 0 个 H 邻居：Ea_eff = Ea
- 1 个 H 邻居：Ea_eff = Ea + 0.10（脱附加速）
- 4 个 H 邻居：Ea_eff = Ea + 0.40

这导致覆盖度存在上限（~0.34 ML for HER），远低于无相互作用时的值。

### 8.2 BEP 关系

BEP (Brønsted-Evans-Polanyi) 关系将活化能修正与反应焓变化关联：

```
BEP = Volmer_fwd : 0.5
```

含义：Ea(env) = Ea(0) + α × (ΔH(env) - ΔH(0))

其中 α = 0.5 为接近因子。当侧向相互作用改变反应焓时，活化能也相应调整。

### 8.3 多位点类型

支持异质表面（如台阶、缺陷）：

```python
# 在 Python API 中设置位点类型
engine.set_site_type((0, 0), 1)     # 位点 (0,0) 设为类型 1
engine.set_site_types_region(
    lambda x, y: x < 5,             # x < 5 的区域
    site_type=2,                     # 设为类型 2
)
```

反应可以指定只在特定位点类型上发生：

```python
pt.add_process(name='step_adsorption', ..., site_type=1)
```

### 8.4 选择性分析

```python
from spark.analysis import compute_selectivity

tof = engine.get_tof()
sel = compute_selectivity(tof, target='NH3_production',
                          total_products=['NH3_production', 'N2H4_production'])
print(f'NH3 selectivity: {sel:.1%}')
```

### 8.5 表观活化能

```python
from spark.analysis import apparent_activation_energy

temperatures = [280, 290, 300, 310, 320]
tofs = []

for T in temperatures:
    engine = KMCEngine(pt, size=[20, 20], print_rates=False, banner=False)
    engine.parameters.T = T
    engine.do_steps(500000)
    engine.get_tof()
    engine.do_steps(500000)
    tofs.append(engine.get_tof()['Heyrovsky'])

Ea_app = apparent_activation_energy(temperatures, tofs)
print(f'Apparent Ea = {Ea_app:.3f} eV')
```

---

## 9. Python API

除了 SPARKIN 输入文件，也可以完全用 Python 构建模型。

### 9.1 从零构建

```python
from spark import Project, Site, Condition, Action
import numpy as np

pt = Project()
pt.set_meta(model_name='my_model', model_dimension=2)

# 物种
pt.add_species(name='empty')
pt.add_species(name='A')
pt.add_species(name='B')

# 晶格
layer = pt.add_layer(name='surface')
layer.sites.append(Site(name='top', pos=(0.5, 0.5, 0.5)))
pt.lattice.cell = np.diag([3.0, 3.0, 15.0])

# 参数
pt.add_parameter(name='T', value=300)
pt.add_parameter(name='Ea_fwd', value=0.5)

# 反应
coord = pt.lattice.generate_coord('top')
pt.add_process(
    name='A_to_B',
    conditions=[Condition(coord, 'A')],
    actions=[Action(coord, 'B')],
    rate_constant='kB*T/h*exp(-Ea_fwd*eV/(kB*T))',
    tof_count={'A_to_B': 1},
)

# 扩散
pt.add_diffusion('A', rate_constant='kB*T/h*exp(-0.1*eV/(kB*T))')

# 侧向相互作用
pt.add_lateral_interaction('A', 'A', energy=0.15)

# BEP
pt.add_bep_relation('A_to_B', alpha=0.5)
```

### 9.2 模型导出

```python
from spark.io import project_to_sparkin, project_to_yaml

# 导出为 SPARKIN 格式
project_to_sparkin(pt, 'my_model.sparkin')

# 导出为 YAML
project_to_yaml(pt, 'my_model.yaml')

# 导出为 JSON
pt.save('my_model.json')
```

### 9.3 模型加载

```python
from spark import load_model

# 自动检测格式
pt = load_model('model.sparkin')   # VASP-style
pt = load_model('model.yaml')     # YAML
pt = load_model('model.json')     # JSON
```

---

## 10. 常见问题

### Q: 覆盖度不变 / TOF 为零

**原因**：步数不够或速率太低。

**解决**：
1. 检查 `print_rates=True` 输出的速率常数是否合理
2. 增加步数：`engine.do_steps(1000000)`
3. 检查温度和活化能：Ea 过高会导致速率接近零

### Q: 仿真很慢

**原因**：扩散速率远高于反应速率，大部分步骤在做无效扩散。

**解决**：
1. 降低扩散 Ea 或完全移除扩散
2. 减小晶格尺寸
3. 使用 Rust 引擎 (`spark-rs`) 获得 20-50x 加速

### Q: SPARKIN 文件报错

**常见错误**：
```
ValueError: model.sparkin:12: Unknown surface 'Cu111'
```
→ 表面名称需要括号：`Cu(111)` 不是 `Cu111`

```
ValueError: Reaction must contain '->'
```
→ 反应方程式必须包含 `->` 箭头

```
ValueError: Mismatch in '...'
```
→ 反应物和产物数量不匹配（`+` 分隔的项数需相等）

### Q: KMC 和 MKM 结果不一致

这是**正常的**。KMC 包含空间关联效应，而 MKM 是均场近似。差异在高覆盖度时尤为明显。一般来说：
- 低覆盖度 (θ < 0.1)：KMC ≈ MKM
- 高覆盖度 (θ > 0.3)：KMC 和 MKM 可能有显著差异
- 有侧向相互作用时：差异更大

### Q: 如何添加新的表面预设？

在 Python 中：

```python
from spark.io import SURFACE_PRESETS
SURFACE_PRESETS['Rh(111)'] = (2.689, 15.0, 'top')
```

或修改 `spark/io.py` 中的 `SURFACE_PRESETS` 字典。

### Q: 速率表达式中可以用什么变量？

| 变量 | 含义 | 值 |
|------|------|----|
| `kB` | Boltzmann 常数 | 1.381e-23 J/K |
| `h` | Planck 常数 | 6.626e-34 J·s |
| `eV` | 电子伏特 | 1.602e-19 J |
| `T` | 温度 | 来自 TEMP 标签 |
| `U` | 电位 | 来自 POTENTIAL 标签 |
| `beta` | 1/(kB×T) | 自动计算 |
| `bar` | 巴 | 1e5 Pa |
| `angstrom` | 埃 | 1e-10 m |
| `umass` | 原子质量单位 | 1.661e-27 kg |
| `pi` | 圆周率 | 3.14159... |
| `m_CO`, `m_H2`, ... | 分子质量 (amu) | 自动查表 |
| `exp`, `sqrt`, `log` | 数学函数 | numpy 实现 |
| 所有 `PARAMETER` | 用户参数 | 来自输入文件 |
