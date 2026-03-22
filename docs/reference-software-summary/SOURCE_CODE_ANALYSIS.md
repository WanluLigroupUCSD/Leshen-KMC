# 参考软件源码深度分析 —— Leshen-KMC 功能移植指南

> 分析日期: 2026-03-21
> 已分析软件: MonteCoffee, kmos3/kmcos, KMCLib, SPPARKS, Zacros 4.0

---

## 1. 各软件核心架构对比

| 维度 | MonteCoffee | kmos3 | KMCLib | SPPARKS | Zacros 4.0 |
|------|------------|-------|--------|---------|------------|
| 语言 | 纯 Python | Python + 生成 Fortran | Python + C++ (SWIG) | C++ | Fortran 2003 |
| KMC算法 | FRM | VSSM (Direct) | VSSM (Direct) | VSSM + 多种solver | FRM + Direct |
| 邻居定义 | 距离截断 neighbor list | 编译期 offset 生成 | cell-based + shells | cmap连接表 | 图邻接表 |
| 横向作用 | get_rate()中用户计算 | OTF后端 bystander扫描 | RateCalculator插件 | site_propensity()中用户计算 | Cluster Expansion |
| 并行 | 无 | 无 | MPI(数据并行匹配) | MPI(sector空间分解) | MPI(Time-Warp) + OpenMP |
| 刚度处理 | α缩放(temporal accel.) | 无 | 无 | Group solver(composition-rejection) | 准平衡检测+缩放 |

---

## 2. 功能移植方案 —— 按优先级排序

### P1: 邻居列表 + 空间依赖事件 🔴

**最佳参考: MonteCoffee + SPPARKS**

#### 数据结构 (参考 MonteCoffee `base/system.py`)

```python
# Python 端
neighbors: list[list[int]]  # neighbors[site] = [site1, site2, ...]

# Rust 端 (CSR格式，缓存友好)
struct NeighborList {
    offsets: Vec<usize>,    # offsets[i]..offsets[i+1] 是 site i 的邻居范围
    neighbors: Vec<usize>,  # 扁平存储所有邻居
}
```

**构建算法** (参考 MonteCoffee `user_system.py`):
```
对于周期性格点 (size_x, size_y):
  FOR each site i at (x, y):
    FOR each offset (dx, dy) in [(-1,0),(1,0),(0,-1),(0,1),...]:
      j = site_at((x+dx) % size_x, (y+dy) % size_y)
      neighbors[i].append(j)
```

**邻居依赖事件** (参考 MonteCoffee `base/events.py`):
```python
# 单位点事件: 只检查自身状态
def possible(system, site, other_site=None):
    return system.sites[site].covered == EMPTY

# 双位点事件: 检查自身+邻居
def possible(system, site, other_site):
    return (system.sites[site].covered == CO and
            system.sites[other_site].covered == O)
```

**事件更新范围** (参考 SPPARKS `app_lattice.h`):
```
delpropensity = 1  # 计算位点速率需要看1层邻居
delevent = 1       # 事件执行影响1层邻居
→ 更新区域 = 事件涉及位点的 (delpropensity + delevent) 层邻居
```

**关键源文件:**
- `MonteCoffee/NeighborKMC/base/system.py` → 邻居构建 + 查询
- `MonteCoffee/NeighborKMC/base/events.py` → 事件接口(possible/get_rate/do_event)
- `MonteCoffee/NeighborKMC/base/kmc.py:284-419` → FRM初始化 + 局部更新
- `SPPARKS/src/app_erbium.cpp:308-362` → site_propensity() 邻居匹配
- `SPPARKS/src/app_erbium.cpp:369-461` → site_event() 执行+更新

---

### P2: 横向相互作用 (Pairwise) 🔴

**最佳参考: kmos3 OTF后端 + MonteCoffee**

#### 方案A: MonteCoffee 式 —— 在 get_rate() 中直接计算

```python
# MonteCoffee/NeighborKMC/user_events.py 示例
def get_rate(self, system, site, other_site):
    # 获取所有邻居的覆盖状态
    ncovs = [system.sites[n].covered for n in system.neighbors[site]]
    # Pairwise 相互作用矩阵
    H_int = [[0, 0, 0],      # empty与任何物种无相互作用
             [0, E_CO_CO, E_CO_O],  # CO-CO, CO-O
             [0, E_O_CO, E_O_O]]    # O-CO, O-O
    repulsion = sum(H_int[system.sites[site].covered][j] for j in ncovs)
    E_ads = E_ads_0 - repulsion
    return (kB*T/h) * exp(-max(E_ads, 0) / (kB*T))
```

**优点:** 实现最简单，完全灵活
**缺点:** 每次评估速率都要扫描邻居，性能开销

#### 方案B: kmos3 OTF式 —— Bystander 计数 + 参数化速率

**核心概念** (参考 `kmos3/kmcos/fortran_src/base_otf.f90`):

```
每个Process定义:
  - conditions: [(site, species_before)]  # 反应物
  - actions: [(site, species_after)]      # 产物
  - bystanders: [(neighbor_coord, allowed_species, flag)]  # 旁观者
  - base_rate: 零覆盖度极限速率
  - otf_rate: "base_rate * exp(beta * nr_CO_1nn * E_CO_nn * eV)"

运行时 (gr_<proc>函数):
  1. 扫描所有 bystander 位点
  2. 按 flag 分组计数各物种数目 → nr_CO_1nn, nr_O_1nn, ...
  3. 代入 otf_rate 表达式计算实际速率
```

**关键数据结构 — rates_matrix** (kmos3独有):
```
rates_matrix[proc][site] = 该process在该site的实际速率 (不再是全局统一)
```
这是OTF后端与标准后端的核心差异。

**事件选择变为两级:**
```
1. 二分搜索选择 process (按 sum(rates_matrix[proc]) 加权)
2. 在选中的 process 内二分搜索选择 site (按 rates_matrix[proc][site] 加权)
```

**关键源文件:**
- `kmos3/kmcos/fortran_src/base_otf.f90:1186-1260` → determine_procsite 两级选择
- `kmos3/kmcos/io.py:2117-2942` → OTF代码生成 (bystander扫描、gr_函数)
- `kmos3/examples/pairwise_interaction_otf__build.py` → OTF用户接口
- `MonteCoffee/NeighborKMC/user_energy.py` → Pairwise相互作用矩阵

**推荐:** 先实现方案A(简单直接)，验证正确性后优化为方案B(更高效)。

---

### P3: 多位点类型 🟡

**最佳参考: kmos3 + SPPARKS**

```python
# kmos3 的 Layer + Site 系统 (types.py)
class Site:
    name: str      # "bridge", "cus", "hollow", "top"
    pos: ndarray   # 单元胞内的分数坐标
    tags: str      # 附加标签

class Layer:
    name: str      # "ruo2", "default"
    sites: list[Site]

# SPPARKS 的做法更简单 (app_erbium.cpp)
type[i] = FCC / OCTA / TETRA  # 每个位点一个整数标签
# 反应匹配时检查: if (type[i] != required_type) continue;
```

**Leshen-KMC 移植方案:**
```rust
struct Site {
    site_type: u16,       // 0=top, 1=bridge, 2=hollow, ...
    species: u16,         // 当前占据物种
    pos: [f64; 3],        // 坐标
}

struct Process {
    conditions: Vec<(Offset, u16, Option<u16>)>,  // (offset, species, site_type)
    // site_type = None 表示不限制
}
```

**关键源文件:**
- `kmos3/kmcos/types.py:1645-1710` → Layer/Site 定义
- `kmos3/kmcos/fortran_src/lattice.mpy` → spuck (Sites Per Unit Cell)
- `SPPARKS/src/app_erbium.cpp:315-320` → type匹配逻辑

---

### P4: 表面扩散 🟡

**最佳参考: MonteCoffee**

扩散本质上就是一个双位点Process:
```
条件: A@site + empty@neighbor
动作: empty@site + A@neighbor
```

**MonteCoffee 的扩散速率计算** (`user_events.py:183-207`):
```python
def get_rate(self, system, site, other_site):
    # 计算初态能量 (A在site)
    E_initial = get_energy(system, site)
    # 临时交换状态计算末态能量 (A在other_site)
    system.sites[site].covered = EMPTY
    system.sites[other_site].covered = A
    E_final = get_energy(system, other_site)
    # 恢复
    system.sites[site].covered = A
    system.sites[other_site].covered = EMPTY
    # 扩散壁垒
    E_diff = max(E_final - E_initial, 0) + E_diff_barrier
    return prefactor * exp(-E_diff / (kB*T))
```

**注意:** 扩散频率通常远高于反应频率，需要配合P5(刚度缩放)。

---

### P5: 刚度缩放 🟡

**三种方案可选，按复杂度递增:**

#### 方案A: MonteCoffee α缩放 (最简单)

**核心思想** (参考 `MonteCoffee/NeighborKMC/base/basin.py:145-216`):
```
超盆(Superbasin)检测:
1. 追踪每种事件的执行次数 nem[event_type]
2. 对正/逆反应对, 计算可逆性指标: rev = |nem[fwd] - nem[rev]| / (nem[fwd] + nem[rev])
3. 若 rev < δ(默认0.2), 标记为准平衡

缩放:
4. 对准平衡事件: α = min(2 * R_nonequil / (R_equil_fwd + R_equil_rev) * Nf, 1)
5. k_scaled = α * k_original
6. 非可逆事件触发时, 退出superbasin, 重置所有α=1
```

**参数:**
- `delta = 0.2` — 可逆性阈值
- `Nf = 1` — 缩放因子
- `Ns = 100` — 每Ns步检查一次
- `ne = 100` — 最小执行次数

**关键源文件:**
- `MonteCoffee/NeighborKMC/base/basin.py:71-142` → scale_rate / scale_rate_constant / scale_constant
- `MonteCoffee/NeighborKMC/base/basin.py:145-216` → superbasin检测逻辑
- `MonteCoffee/NeighborKMC/base/basin.py:11-44` → FRM时间重缩放

#### 方案B: SPPARKS Group Solver (O(1)选择)

**核心思想** (参考 `SPPARKS/src/groups.cpp`):
```
将所有事件按速率分组到对数尺度的桶中:
  Group 0: [hi/2, hi]
  Group 1: [hi/4, hi/2]
  Group 2: [hi/8, hi/4]
  ...

选择事件:
  1. 按组总速率选组 (线性扫描，组数少)
  2. 在组内随机选事件，用 rejection 接受/拒绝
  → 接受概率 = rate[i] / group_upper_bound ≥ 0.5 (对数分组保证)

更新: O(1) — 从旧组移除(swap-with-last)，加入新组
```

**关键源文件:**
- `SPPARKS/src/groups.cpp:96-130` → 对数分组逻辑
- `SPPARKS/src/groups.cpp:203-213` → composition-rejection采样
- `SPPARKS/src/solve_group.cpp:120-136` → propensity clamping(近似)

#### 方案C: Zacros 准平衡检测 (最精确)

```
周期性检查 (每 icheckeverynoccur 事件):
  FOR each 正/逆反应对 (i, j):
    检查 ratio = rate_fwd * N_fwd / (rate_rev * N_rev)
    若 |1 - ratio| < threshold:
      同步缩小: rate_fwd *= factor, rate_rev *= factor
      确保 time_scale_separation > min_threshold
```

**关键源文件:**
- `zacros_4.0/stiffness_scaling_module.f90` → 完整实现

**推荐:** 先实现方案A (MonteCoffee α缩放)，简单有效。后续可升级为方案B。

---

### P6: BEP 关系 🟢

**最佳参考: kmos3**

```python
# kmos3/kmcos/interactions/BEPmodule.py
class BEPRelation:
    alpha: float = 0.5  # BEP斜率
    beta: float = 0.0   # BEP截距

    def Ea_from_DeltaH(self, DeltaH):
        return self.alpha * DeltaH + self.beta

    def change_in_Ea(self, change_in_DeltaH):
        return self.alpha * change_in_DeltaH

    def reverse(self):
        return BEPRelation(1 - self.alpha, self.beta)
```

**与横向作用整合** (`configurationsAndInteractions.py:592`):
```
1. 计算反应物态的总相互作用能: E_reactant = Σ pairwise(conditions, bystanders)
2. 计算产物态的总相互作用能: E_product = Σ pairwise(actions, bystanders)
3. ΔΔH = E_product - E_reactant
4. ΔEa = α * ΔΔH
5. k_new = k_base * exp(-ΔEa / kBT)
```

**注意:** kmos3 的 OTF 后端与 BEP 不兼容 (因为OTF只知道反应物态的bystander，不知产物态)。Leshen-KMC 可以解决这个问题，因为 action_list 已知，产物态的bystander配置可以推断。

---

### P7: 高效事件选择数据结构 🟢

**三种可选:**

| 数据结构 | 选择 | 更新 | 参考实现 |
|---------|------|------|---------|
| 线性扫描 | O(N) | O(1) | SPPARKS solve_linear |
| 二叉求和树 | O(log N) | O(log N) | SPPARKS solve_tree |
| 索引二叉堆 | O(1) | O(log N) | Zacros execution_queue_binary_heap |
| Composition-rejection | O(1)摊销 | O(1) | SPPARKS solve_group |

**SPPARKS 二叉求和树** (`solve_tree.cpp`):
```
数组布局: tree[0..2N-1]
  - tree[N-1..2N-2] = 叶子 (各事件的propensity)
  - tree[k] = tree[2k+1] + tree[2k+2] (内部节点 = 子节点之和)
  - tree[0] = 总速率

更新 set(i, value):
  tree[N-1+i] = value
  k = N-1+i
  while k > 0:
    parent = (k-1)/2
    tree[parent] = tree[2*parent+1] + tree[2*parent+2]
    k = parent

选择 find(target):
  k = 0
  while k < N-1:
    if target <= tree[2k+1]: k = 2k+1
    else: target -= tree[2k+1]; k = 2k+2
  return k - (N-1)
```

**关键源文件:**
- `SPPARKS/src/solve_tree.cpp:174-224` → set() + find()
- `SPPARKS/src/solve_group.cpp` + `SPPARKS/src/groups.cpp` → composition-rejection
- `zacros_4.0/execution_queue_binary_heap_class.f90` → 索引堆

---

### P8: 检查点/重启 🟢

**最佳参考: Zacros 4.0**

Zacros 每3600秒写 `restart.inf`，包含完整状态:
- 当前时间、步数、随机数种子
- 格点状态 (每个位点的物种)
- 事件队列
- 统计数据

**Leshen-KMC 方案:** JSON序列化 (已有基础设施)
```json
{
  "time": 1.23e-5,
  "step": 100000,
  "rng_state": "...",
  "lattice": [0, 1, 3, 0, 2, ...],
  "process_stats": {...},
  "tof_counts": {...}
}
```

---

### P9: 并行化 (长期) 🟢

**最佳参考: SPPARKS sector方法**

```
将格点划分为 2×2=4 个sector:
  FOR each sector:
    1. 交换边界ghost数据
    2. 独立运行KMC dt_kmc时间
    3. 回传更新状态

dt_kmc = nstop / pmax_avg  (自适应)
```

**Rust实现:** 使用 Rayon 的 `par_iter` 处理sector内并行，无需MPI。

**关键源文件:**
- `SPPARKS/src/app_lattice.cpp:517-626` → iterate_kmc_sector
- `SPPARKS/src/comm_lattice.h` → ghost数据交换

---

## 3. Leshen-KMC 现有代码的改动点

### engine.py / engine.rs 改动

```
现有结构:
  lattice: ndarray[nsites]         → species per site
  _avail_sites: dict[proc_id → list[site_id]]
  _site_in_avail: dict[proc_id → dict[site_id → index]]

需要新增:
  neighbors: list[list[int]]       → 邻居列表 (P1)
  site_types: ndarray[nsites]      → 位点类型 (P3)
  interaction_matrix: dict[(sp1,sp2) → float]  → pairwise相互作用 (P2)

需要修改:
  _check_available(): 当前只检查本位点 → 需要检查邻居状态 (P1)
  _update_after_event(): 当前更新 max_offset 范围 → 需要更新邻居链 (P1)
  _compute_rate(): 当前返回固定值 → 需要考虑邻居环境 (P2)
```

### types.py 改动

```
Site 类新增:
  site_type: int = 0

Process 类新增:
  neighbor_conditions: list[(neighbor_index, species)]  # 邻居条件
  bystanders: list[(offset, species_flag)]  # OTF bystander
  otf_rate_expr: str  # OTF速率表达式
  reverse_process: Optional[int]  # 逆反应ID (用于刚度缩放)
```

---

## 4. 推荐实现顺序

```
Week 1-2: P1 邻居列表 + 空间依赖事件
  → 先在 Python 端实现, 用 CO氧化/Pt 验证
  → 对比 well-mixed vs spatial 的覆盖度差异

Week 3-4: P2 Pairwise 横向相互作用 (方案A)
  → 在 get_rate() 中加入邻居计数
  → 验证: 覆盖度 vs 相互作用强度的等温线

Week 5: P3 多位点类型
  → Site 加 type 字段, Process 加 type 匹配
  → 测试: bridge + top 双位点 CO吸附

Week 6: P4 表面扩散
  → 实现 hop 过程
  → 验证: 扩散系数 vs 解析解

Week 7-8: P5 刚度缩放 (方案A: α缩放)
  → 追踪正逆反应执行次数
  → 验证: 有/无缩放的稳态一致性

Week 9: 移植到 Rust
  → 核心数据结构: CSR邻居表, site_type, interaction_matrix
  → 关键算法: 邻居依赖 check_available, OTF get_rate
```

---

## 5. 各软件源码路径索引

| 功能 | 软件 | 关键文件 |
|------|------|---------|
| 邻居列表构建 | MonteCoffee | `NeighborKMC/base/system.py:37-58` |
| 邻居依赖事件 | MonteCoffee | `NeighborKMC/base/events.py` (全文) |
| FRM + 局部更新 | MonteCoffee | `NeighborKMC/base/kmc.py:284-419` |
| OTF横向作用 | kmos3 | `kmcos/fortran_src/base_otf.f90:1186-1260` |
| Bystander定义 | kmos3 | `kmcos/types.py:2128`, `kmcos/io.py:2117-2942` |
| BEP关系 | kmos3 | `kmcos/interactions/BEPmodule.py` |
| 模式匹配 | KMCLib | `c++/src/matcher.cpp:41-198` |
| 自定义速率插件 | KMCLib | `c++/src/ratecalculator.h`, `python/.../KMCRateCalculatorPlugin.py` |
| 多粒子TypeBucket | KMCLib | `c++/src/typebucket.h` |
| 二叉求和树 | SPPARKS | `src/solve_tree.cpp:174-224` |
| Group solver | SPPARKS | `src/groups.cpp:96-213` |
| Sector并行 | SPPARKS | `src/app_lattice.cpp:517-626` |
| 催化反应匹配 | SPPARKS | `src/app_erbium.cpp:308-461` |
| α缩放 (temporal accel.) | MonteCoffee | `NeighborKMC/base/basin.py:71-216` |
| 准平衡刚度缩放 | Zacros | `stiffness_scaling_module.f90` |
| Cluster Expansion | Zacros | `energetics_parser_module.f90`, `energetics_handle_module.F90` |
| 子图同构 | Zacros | `SubIsoSolver.f90`, `SubIsoSolver_VF2.f90` |
