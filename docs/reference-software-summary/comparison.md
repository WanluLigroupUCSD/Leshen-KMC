# SPARK 与现有动力学蒙特卡洛软件的系统对比分析

**SPARK: Kinetic Monte Carlo & Microkinetic Modeling for Heterogeneous Electrocatalysis**

*Wanlu Li Group, UC San Diego*

---

## 1. 引言：现有 KMC 软件生态概览

动力学蒙特卡洛（KMC）模拟结合第一性原理计算，已成为桥接多相催化中多尺度时间和空间的黄金标准计算框架。过去二十年中，涌现了多个通用 KMC 软件包，主要包括：

| 软件 | 开发方 | 语言 | 开源 | 首次发布 | 主要应用领域 |
|------|--------|------|------|----------|-------------|
| **CARLOS** | J. Lukkien | C | 否 | ~2000 | 催化表面反应（最早期） |
| **SPPARKS** | Sandia National Labs | C++ | GPL | 2009 | 材料科学（晶粒生长、薄膜、增材制造） |
| **Zacros** | UCL / Oxford (Stamatakis) | Fortran 2003 | 学术免费 | 2013 | 催化表面反应（最活跃） |
| **kmos/kmos3** | Hoffmann, Reuter 等 | Python + Fortran | GPL | 2014 | 催化表面反应 |
| **KMCLib** | Leetmaa, Skorodumova | Python + C++ | GPL | 2014 | 通用格点扩散/反应 |
| **MoCKA** | KIT (Deutschmann) | — | — | 2015 | 纳米颗粒催化 |
| **MonteCoffee** | Chalmers (Grönbeck) | Python | MIT | 2018 | 纳米颗粒催化 |
| **SPARK** | UCSD (Wanlu Li) | Python + Rust | MIT | 2025 | **电催化（唯一）** |

**关键事实：上述所有软件（除 SPARK 外）均为热催化设计，没有任何一个原生支持电催化功能。**

---

## 2. 核心算法对比

### 2.1 KMC 算法选择

| 软件 | 算法 | 事件选择复杂度 | 拒绝率 |
|------|------|---------------|--------|
| **SPARK** | BKL / VSSM (Variable Step-Size Method) | O(log N_proc) binary search | 零拒绝（rejection-free） |
| **Zacros** | FRM (First Reaction Method) | O(N_events) 或优化后更低 | 零拒绝 |
| **kmos3** | Direct Method (BKL/Gillespie) | O(1)（编译期生成优化代码） | 零拒绝 |
| **SPPARKS** | 多种可选（VSSM, null-event, rejection-free） | O(1) ~ O(N log N) | 取决于算法 |
| **KMCLib** | Direct Method | 标准实现 | 零拒绝 |
| **MonteCoffee** | FRM | O(N_events) | 零拒绝 |

**分析：**

- SPARK 采用 BKL/VSSM 算法，与 kmos 同族。核心区别在于实现方式：kmos 在"导出"（export）阶段根据具体模型编译生成优化的 Fortran 代码，使得运行时局部更新达到 O(1)——即运行时间与格点大小完全无关。SPARK 使用通用引擎（Python/Rust），通过 swap-with-last 技巧实现 O(1) 的位点增删操作，过程选择通过二分搜索完成，复杂度为 O(log N_proc)。对于催化反应（通常 N_proc ≤ 20），这个差异可以忽略。

- Zacros 的 FRM 为每个可行事件独立生成发生时间，选取最早者执行。对于事件频繁增删的复杂系统（如含多齿物种和长程横向相互作用），FRM 的优势更明显。

- **SPARK 的独特之处：Rust 高性能后端。** 所有现有 KMC 软件的性能瓶颈后端都是 Fortran 或 C++，SPARK 是第一个提供 Rust 实现的 KMC 框架，实测比 Python 版本加速约 24 倍，且编译产物仅 1.3 MB，无外部依赖。

### 2.2 格点表示与事件匹配

| 软件 | 格点类型 | 事件匹配方法 | 空间关联性 |
|------|---------|-------------|-----------|
| **SPARK** | 均匀格点（well-mixed） | 直接查表 | ✗ 无 |
| **Zacros** | 任意拓扑格点（图表示） | 子图同构匹配（VF2/RI 算法） | ✅ 完整空间关联 |
| **kmos3** | 周期性格点（1D-3D） | 编译期代码生成 | ✅ 邻居依赖事件 |
| **SPPARKS** | 格点 + 非格点 | 用户自定义 | ✅ 空间分解 |
| **KMCLib** | 周期性格点（1D-3D） | 标准匹配 | ✅ 支持 |
| **MonteCoffee** | 非格点（邻居列表） | 邻居列表查询 | ✅ 位点连通性 |

**分析：**

Zacros 的 Graph-Theoretical（GT）方法是该领域最先进的事件匹配框架：将格点和基元步骤都表示为图（graph），通过子图同构问题的求解来自动识别格点上所有可行事件。这使其能天然处理多齿物种（如 CO* 占据 bridge 位需两个相邻位点）、复杂邻居模式（如要求反应物 A* 和 B* 相邻且第三邻位为空）等情况。

**SPARK 目前使用 well-mixed 格点模型，即每个位点独立选择过程，不考虑邻居状态。** 这是一个重要的简化假设——等价于假设表面扩散无限快（mean-field 极限）。对于覆盖度主导的体系，这个近似合理；对于岛形成、pattern formation、扩散限制等现象，则无法捕捉。

### 2.3 平均场微动力学（MKM）求解器

| 软件 | 内置 MKM | ODE 求解器 | 稳态求解 |
|------|---------|-----------|---------|
| **SPARK** | ✅ 内置（核心功能） | SciPy solve_ivp（Python）, RK4+自适应 Euler（Rust） | fsolve（Python）, 自适应步长收敛（Rust） |
| **Zacros** | ✗（需外部 MKMCXX） | — | — |
| **kmos3** | ✗（需外部 CatMAP） | — | — |
| **SPPARKS** | ✗ | — | — |
| **KMCLib** | ✗ | — | — |
| **MonteCoffee** | ✗ | — | — |

**分析：**

**这是 SPARK 的核心差异化优势之一。** 所有现有 KMC 软件都是纯 KMC 引擎，要进行平均场微动力学分析需要：(1) 将反应模型手动转换到 CatMAP 或 MKMCXX 的格式，(2) 分别运行两个独立软件，(3) 手动对比结果。

SPARK 将 MKM 和 KMC 集成在同一框架内，使用完全相同的反应模型定义（Species, Process, rate expressions）。用户可以对同一个模型：
- 运行 MKM 快速获得稳态覆盖度和 TOF
- 运行 KMC 获得随机动力学和时间演化
- 直接对比两者结果验证 mean-field 近似的有效性

这对方法论验证和快速筛选-精细模拟的工作流非常有价值。

---

## 3. 功能特性全面对比

### 3.1 电催化功能（核心差异化）

| 电催化功能 | **SPARK** | Zacros | kmos3 | SPPARKS | KMCLib | MonteCoffee |
|-----------|:-----------:|:------:|:-----:|:-------:|:-----:|:-----------:|
| 电位依赖速率常数（PCET） | ✅ Butler-Volmer 原生 | ✗ | ✗ | ✗ | ✗ | ✗ |
| 极化曲线输出（j-V） | ✅ 直接计算 | ✗ | ✗ | ✗ | ✗ | ✗ |
| DFT 恒电位数据插值 | ✅ 三次样条 | ✗ | ✗ | ✗ | ✗ | ✗ |
| 电化学步/热反应步区分 | ✅ PCET vs Thermal | ✗ | ✗ | ✗ | ✗ | ✗ |
| 对称因子 β / 电荷转移系数 | ✅ 输入参数 | ✗ | ✗ | ✗ | ✗ | ✗ |
| 电位扫描（内置） | ✅ CLI 一行命令 | ✗ | ✗ | ✗ | ✗ | ✗ |
| 法拉第效率追踪 | ✅（通过 TOF 分产物） | ✗ | ✗ | ✗ | ✗ | ✗ |

**详细说明：**

**(1) Butler-Volmer PCET 速率**

SPARK 原生区分质子耦合电子转移（PCET）步骤和纯热化学步骤。对于 PCET 步骤：

```
k(U) = (kB·T/h) · exp(−(Ea + β_BV·U) · eV / (kB·T))
```

对于热化学步骤（如 N-N 键断裂）：

```
k = (kB·T/h) · exp(−Ea · eV / (kB·T))
```

其他 KMC 软件中的速率常数仅依赖温度，不包含电位项。文献中做电催化 KMC 的少数工作（如 Wei et al., JPCL 2025 的 CO₂RR/Cu KMC 预测线性扫描伏安法；Li et al., ACS Catal. 2024 的 HER/Pt KMC 预测极化曲线）均使用各组自行编写的一次性代码，而非通用软件。

**(2) 极化曲线直接输出**

SPARK 提供两种模式计算极化曲线：

- **模式 1（Butler-Volmer）：** 使用 MKM 或 KMC 在不同电位下运行，获得 TOF，通过 j = TOF × n_e × e × N_site / A_geo 转换为电流密度。
- **模式 2（DFT 恒电位数据）：** 输入多个电位下的 DFT 态能量和过渡态能量（JSON 格式），用三次样条自动插值出任意电位下的能垒，再运行 MKM 获得 j-V 曲线。

这对接了 CP-VASP / TPOT / GC-DFT 等恒电位 DFT 方法的输出，实现了从 DFT 到极化曲线的完整计算管道。

**(3) 电化学/热反应步骤分类**

在 NRR 反应网络中，N₂* → NNH* 是 PCET 步（速率与电位相关），而 NNH₂* → N* + NH₃ 是纯热步（N-N 键断裂，速率与电位无关）。SPARK 的模型定义中明确标记每个步骤的类型，引擎自动应用正确的速率表达式。

### 3.2 反应模型定义与灵活性

| 功能 | **SPARK** | Zacros | kmos3 | MonteCoffee |
|------|:-----------:|:------:|:-----:|:-----------:|
| 速率表达式解析器 | ✅ 字符串表达式，运行时解析 | ✗（keyword 文件） | ✅ 类似语法 | ✗（用户代码） |
| kmos 兼容 API | ✅ Project/Species/Process/Condition/Action | — | 原生 | — |
| JSON 模型 I/O | ✅ | ✗ | ✗ (.ini/.xml) | ✗ |
| 物理常数库 | ✅ kB, h, eV, bar, Å, umass, 分子质量 | 内置 | ✅ | ✗ |
| 双模式（Full/Cycle） | ✅ | ✗ | ✗ | ✗ |
| 预置模型 | ✅ NRR/Mo + CO/Pd(100) | ✗ | ✗ | CO/Pt 示例 |
| 教程 | ✅ 7 部分完整教程 | ✅ 4 个教程 | ✅ 文档 | 基本文档 |

**分析：**

SPARK 采用了与 kmos 兼容的 API 概念（Project, Species, Site, Layer, Process, Condition, Action），降低了已有 kmos 用户的学习成本。

JSON 模型 I/O 是一个实用创新——模型可以在 Python 和 Rust 之间无缝共享，也便于版本控制和与其他工具集成。Zacros 使用自定义 keyword 文件格式，kmos 使用 .ini 或 .xml，均不如 JSON 通用。

双模式（Full + Cycle）设计也是独有的：Full mode 包含完整的吸附/脱附动力学，适合研究覆盖度效应；Cycle mode 假设所有位点始终被占据，仅考虑催化循环中的化学转化，适合快速机理筛选。两种模式共享同一套反应参数，可直接对比。

### 3.3 表面物理与横向相互作用

| 功能 | **SPARK** | Zacros | kmos3 | SPPARKS | MonteCoffee |
|------|:-----------:|:------:|:-----:|:-------:|:-----------:|
| Cluster Expansion Hamiltonian | ✗ | ✅ 多体项 | 有限 | ✗ | ✗ |
| BEP 关系（on-the-fly） | ✗ | ✅ | ✅ OTF 后端 | ✗ | 用户代码 |
| 吸附物间相互作用 | ✗ | ✅（pair, trio, ...） | ✅（pair） | 用户定义 | 用户定义 |
| 表面扩散 | ✗ | ✅ | ✅ | ✅ | ✅ |
| 多齿物种 | ✗ | ✅ GT 框架原生 | ✅ | 需自定义 | 用户代码 |
| 空间关联/岛形成 | ✗ | ✅ | ✅ | ✅ | ✅ |
| 纳米颗粒多面 | ✗ | 需自定义格点 | ✗ | 部分 | ✅（CN-based） |

**这是 SPARK 当前最大的局限。** 横向相互作用和空间关联效应对于高覆盖度体系（如 CO₂RR 中的高 CO* 覆盖度）至关重要。Zacros 的 Cluster Expansion Hamiltonian 可以系统地处理多体相互作用对吸附能和活化能的影响，而 SPARK 的 well-mixed 格点模型隐含假设吸附物间无相互作用、扩散无限快。

### 3.4 模拟分析工具

| 功能 | **SPARK** | Zacros | kmos3 | MonteCoffee |
|------|:-----------:|:------:|:-----:|:-----------:|
| 覆盖度统计 | ✅ | ✅ | ✅ | ✅ |
| TOF 计算 | ✅ | ✅ | ✅ | ✅ |
| 轨迹记录 | ✅ TrajectoryRecorder | ✅ | ✅ | ✅ |
| 参数扫描（T, U） | ✅ 内置 CLI | 需脚本 | 需脚本 | 需脚本 |
| 稳态检测 | ✅ 自动（dθ/dt 收敛） | 用户判断 | 用户判断 | 用户判断 |
| 灵敏度分析 (DRC) | ✗ | ✅ | ✗ | ✗ |
| TPD/TPR 模拟 | ✗ | ✅ 内置 | 可实现 | 可实现 |
| 可视化 | 基本（数据输出） | ✅ Zacros-post GUI | ✅ 内置 GUI | ✗ |
| 吸附层构型快照 | ✗ | ✅ | ✅ | ✅ |

### 3.5 性能与可扩展性

| 特征 | **SPARK** | Zacros | kmos3 | SPPARKS | KMCLib |
|------|:-----------:|:------:|:-----:|:-------:|:-----:|
| 分布式并行 (MPI) | ✗ | ✅ Time-Warp | ✗ | ✅ 空间分解 | ✅ |
| 多线程 | ✗ | ✗ | ✗ | ✗ | ✗ |
| 高性能后端 | ✅ Rust（24× vs Python） | Fortran 2003 | Fortran（编译期优化） | C++ | C++ |
| 最大格点规模 | ~10⁴ 位点（单核） | ~10⁷ 位点（分布式） | ~10⁴ 位点 | ~10⁸ 位点 | ~10⁶ 位点 |
| 二进制大小 | 1.3 MB | ~数 MB | 取决于模型 | ~数 MB | 需编译 |
| 验证测试 | ✅ Langmuir 等温线 | ✅ 多项 | ✅ | ✅ | ✅ |

---

## 4. SPARK 的核心优势总结

### 优势 1：首个面向电催化的通用 KMC 框架

这是 SPARK 最重要的差异化特征。所有现有通用 KMC 软件都是为热催化设计的，要做电催化 KMC 需要：(a) 自行编写电位依赖的速率函数，(b) 手动在不同电位下扫描运行，(c) 自行编写 TOF 到电流密度的转换。SPARK 将这些作为一等公民内置，提供从 DFT 数据到极化曲线的一站式工作流。

### 优势 2：KMC + MKM 同框对比

唯一一个在同一框架内同时支持格点 KMC 模拟和平均场微动力学求解的软件。这使得：
- 快速 MKM 筛选 → 精细 KMC 验证的工作流无缝衔接
- 直接定量评估 mean-field 近似对特定体系的误差
- 同一反应模型无需重复定义

### 优势 3：双语言架构（Python + Rust）

Python 层提供灵活的研究接口和丰富的科学计算生态（NumPy, SciPy, Matplotlib），适合探索性研究和分析。Rust 层提供约 24 倍加速，适合参数扫描和生产级模拟。两层通过 JSON 模型格式互通。这是 KMC 领域首个提供 Rust 后端的软件。

### 优势 4：DFT 恒电位数据直接对接

支持输入多个电位下的 DFT 态能量和过渡态能量（JSON 格式），用三次样条自动插值出连续的电位-能垒关系。这直接对接了 CP-VASP、TPOT、GC-DFT 等恒电位 DFT 方法的输出格式。

### 优势 5：低门槛与即时可用性

- JSON 模型 I/O 便于版本控制和跨工具集成
- kmos 兼容 API 降低学习成本
- 内置 NRR 和 CO 吸附预置模型
- 7 部分完整教程
- CLI 命令一行即可运行模拟、参数扫描、极化曲线计算
- Rust 二进制无外部依赖，仅 1.3 MB

---

## 5. SPARK 的不足与未来发展方向

### 不足 1：缺乏空间关联效应（最关键）

**现状：** SPARK 使用 well-mixed 格点模型，每个位点独立选择过程，不考虑邻居状态。这等价于 mean-field 近似在空间维度上的退化。

**影响：** 无法捕捉以下重要物理现象：
- 吸附物岛状聚集（island formation）
- 表面 pattern formation 和空间振荡
- 扩散限制导致的反应速率偏差
- 邻居依赖的反应事件（如 A* + B* → AB*，要求 A 和 B 相邻）

**对比：** Zacros、kmos、SPPARKS 等均支持真实空间格点上的邻居依赖事件和表面扩散模拟。

**改进方向：** 引入邻居列表和空间依赖事件选择机制。可参考 MonteCoffee 的邻居列表方法（相对简单）或 Zacros 的 GT 框架（更通用但实现复杂）。

### 不足 2：缺乏横向相互作用处理

**现状：** 速率常数独立于表面构型。

**影响：** 在高覆盖度体系中（如 CO₂RR 中 CO* 覆盖度可达 >0.5 ML），吸附物间的排斥/吸引相互作用会显著影响吸附能和活化能。忽略这一效应可能导致定量预测的系统偏差。

**对比：** Zacros 的 Cluster Expansion Hamiltonian 可系统处理 pair、trio 及更高阶相互作用；kmos 支持 on-the-fly BEP 关系修正。

**改进方向：** 最简方案是引入基于覆盖度的线性修正（类似 CatMAP 的 cross-interaction 参数）；更完善的方案是实现 pair-wise lateral interaction 框架。

### 不足 3：无分布式并行能力

**现状：** 单核运行（Rust 后端虽快但仍为串行）。

**影响：** 格点规模受限于约 10⁴ 位点。对于需要大格点的问题（如模拟纳米颗粒上的长程空间效应、μm 尺度的 pattern formation），当前架构不足。

**对比：** Zacros 的 Time-Warp 分布式并行可处理 10⁷ 位点；SPPARKS 的空间分解并行可达 10⁸ 位点级别。

**改进方向：** Rust 生态有良好的并发支持（Rayon, tokio），可相对容易地引入域分解并行。但对于 KMC 的并行化，需要处理因果性违反问题（Time-Warp 或保守同步算法）。

### 不足 4：无 TPD/TPR 模拟能力

**现状：** 不支持模拟升温脱附/升温反应谱。

**影响：** TPD/TPR 是表面科学中验证动力学模型的核心实验手段。无法直接模拟 TPD 谱意味着模型参数的实验验证受限。

**对比：** Zacros 原生支持 TPD/TPR 模拟。

**改进方向：** 实现时间依赖的温度斜坡功能（T(t) = T₀ + β·t），在每步 KMC 中更新温度和对应的速率常数。

### 不足 5：灵敏度分析工具不足

**现状：** 未提供系统的灵敏度分析（如 Degree of Rate Control, DRC）工具。

**影响：** 灵敏度分析对于识别速率决定步骤和关键参数至关重要，是 KMC 模拟结果物理解释的核心工具。

**对比：** Zacros 提供内置的 DRC 分析框架。

**改进方向：** 实现 Campbell 的 DRC 分析方法，或基于 likelihood ratio 方法的灵敏度分析（Vlachos 组 Zacros-Wrapper 中有实现）。

### 不足 6：格点可视化缺失

**现状：** 输出为数值数据，无图形化格点构型展示。

**影响：** 空间分辨的覆盖度快照对于理解表面动力学（如吸附层有序结构、相分离）非常有价值。

**对比：** Zacros 有 Zacros-post GUI 展示格点构型；kmos 内置格点可视化。

**改进方向：** 可利用 Matplotlib 或 CatGo 项目中的 MatterViz 实现格点构型的可视化。

### 不足 7：仅验证了 Langmuir 等温线

**现状：** 验证测试仅包含 CO/Pd(100) 的 Langmuir 吸附-脱附平衡。

**影响：** 单一的验证测试不足以充分证明引擎的正确性，尤其是对于复杂反应网络和电化学速率函数的实现。

**对比：** 成熟 KMC 软件通常包含多种验证：解析可解模型（如 ZGB 模型的相变点）、与 mean-field 解的系统对比、不同算法间的交叉验证。

**改进方向：** 添加更多验证测试，如 ZGB 模型的相变点、与 CatMAP mean-field 解的对比、已发表 KMC 文献结果的复现。

---

## 6. 竞争定位图

```
         电催化支持
           ↑
      高   │  ★ SPARK     （唯一占据此象限）
           │
           │
           │
      低   │  Zacros  kmos3  MonteCoffee  SPPARKS  KMCLib
           └──────────────────────────────────────────→
              简单                                 复杂
                     表面物理模型复杂度
```

SPARK 填补的是「电催化 + KMC」这个交叉领域的工具空白。现有软件在热催化格点 KMC 方面更加成熟（尤其是 Zacros 和 kmos），但在电催化方向没有通用工具。SPARK 的战略定位是成为**电催化 KMC 的默认选择**，同时逐步补齐空间关联和横向相互作用等表面物理功能。

---

## 7. 主要参考文献

### SPARK
- GitHub: https://github.com/WanluLigroupUCSD/SPARK

### Zacros
- Stamatakis & Vlachos, J. Chem. Phys. 134, 214115 (2011). DOI: 10.1063/1.3596751
- Nielsen et al., J. Chem. Phys. 139, 224706 (2013). DOI: 10.1063/1.4840395
- Ravipati et al., Comput. Phys. Commun. 270, 108148 (2022). DOI: 10.1016/j.cpc.2021.108148
- Prats, J. Phys. Chem. A 129, 6608–6614 (2025). DOI: 10.1021/acs.jpca.5c02802
- https://zacros.org/

### kmos / kmos3
- Hoffmann, Matera & Reuter, Comput. Phys. Commun. 185, 2138–2150 (2014). DOI: 10.1016/j.cpc.2014.04.003
- Andersen, Panosetti & Reuter, Front. Chem. 7, 202 (2019). DOI: 10.3389/fchem.2019.00202
- https://www.kmos3.org/doc/

### SPPARKS
- Plimpton et al., Sandia Report SAND2009-6226 (2009)
- https://spparks.github.io/

### KMCLib
- Leetmaa & Skorodumova, Comput. Phys. Commun. 185, 2340–2349 (2014). DOI: 10.1016/j.cpc.2014.04.017
- https://github.com/leetmaa/KMCLib

### MonteCoffee
- Jørgensen & Grönbeck, J. Chem. Phys. 149, 114101 (2018). DOI: 10.1063/1.5046635
- https://gitlab.com/ChemPhysChalmers/MonteCoffee

### MoCKA
- Kunz, Kuhn & Deutschmann, J. Chem. Phys. 143, 044108 (2015). DOI: 10.1063/1.4926924

### 综述
- Pineda & Stamatakis, J. Chem. Phys. 156, 120902 (2022). DOI: 10.1063/5.0083251
- Stamatakis & Vlachos, ACS Catal. 2, 2648–2663 (2012). DOI: 10.1021/cs3005709

### 电催化 KMC 文献（非通用软件，自编代码）
- Wei et al., J. Phys. Chem. Lett. 16, 2896–2904 (2025). DOI: 10.1021/acs.jpclett.4c03426 — CO₂RR on Cu KMC
- Li et al., ACS Catal. (2024) — HER on Pt KMC with polarization curve
