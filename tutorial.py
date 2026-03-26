#!/usr/bin/env python3
"""
=============================================================================
  spark 完整使用教程
  ——从零开始构建动力学蒙特卡洛与微观动力学模型
=============================================================================

本教程分为 7 个部分：
  第1部分：构建简单模型（CO 吸附/脱附）
  第2部分：运行 KMC 模拟
  第3部分：参数动态调节与覆盖度分析
  第4部分：构建复杂反应网络（N₂ 还原）
  第5部分：平均场微观动力学 ODE 求解
  第6部分：参数扫描（电位、温度）与速率控制度分析
  第7部分：结果可视化

运行方式：
  python3.11 tutorial.py           # 运行全部教程
  python3.11 tutorial.py --part 1  # 只运行第1部分
  python3.11 tutorial.py --part 5  # 只运行第5部分
"""

import sys
import os
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ═══════════════════════════════════════════════════════════════════════════
#  第1部分：构建简单模型 —— CO 在 Pd(100) 上的吸附/脱附
# ═══════════════════════════════════════════════════════════════════════════

def part1_build_simple_model():
    """
    演示如何从零构建一个 KMC 模型。
    这是 kmos 官方教程的 Python 复现版。
    """
    print("\n" + "=" * 70)
    print("  第1部分：构建简单模型 —— CO/Pd(100) 吸附脱附")
    print("=" * 70)

    # --- 步骤 1：导入必要模块 ---
    from spark import (
        Project, Species, Site, Layer,
        Condition, Action, Parameter, Coord,
    )

    # --- 步骤 2：创建项目并设置元数据 ---
    pt = Project()
    pt.set_meta(
        author='Tutorial',
        email='user@example.com',
        model_name='CO_on_Pd100',    # 模型名称
        model_dimension=2,            # 2D 表面模型
    )
    print("\n[步骤2] 项目创建完成")

    # --- 步骤 3：定义物种 ---
    # 重要：第一个添加的物种是默认物种（所有位点的初始状态）
    pt.add_species(name='empty', color='#ffffff')  # 空位
    pt.add_species(name='CO', color='#ff0000')     # CO 吸附态

    print(f"[步骤3] 物种定义完成：{[s.name for s in pt.species_list]}")

    # --- 步骤 4：定义晶格 ---
    # 创建一个简立方层，包含一个 hollow 位点
    layer = pt.add_layer(name='simple_cubic')
    layer.sites.append(
        Site(name='hollow',
             pos=(0.5, 0.5, 0.5),          # 分数坐标
             default_species='empty')        # 初始为空
    )

    # 设置单胞大小（Å）
    pt.lattice.cell = np.diag([3.5, 3.5, 10.0])

    print(f"[步骤4] 晶格定义完成：{layer.name}, 位点 = {[s.name for s in layer.sites]}")

    # --- 步骤 5：定义参数 ---
    # adjustable=True 的参数可以在运行时动态修改
    pt.add_parameter(name='T', value=600.0,
                     adjustable=True, min=400, max=800)
    pt.add_parameter(name='p_CO', value=1.0,
                     adjustable=True, min=1e-10, max=100)
    pt.add_parameter(name='A', value='(3.5*angstrom)**2')  # 位点面积
    pt.add_parameter(name='deltaG', value=-0.5,            # 吸附自由能 [eV]
                     adjustable=True, min=-1.3, max=0.3)

    print(f"[步骤5] 参数定义完成：T={600}K, p_CO={1}bar, ΔG={-0.5}eV")

    # --- 步骤 6：定义基元过程 ---
    # 首先生成坐标引用
    coord = pt.lattice.generate_coord('hollow')

    # 过程 1：CO 吸附（* → CO*）
    # 速率 = Hertz-Knudsen 碰撞频率
    # k_ads = p * A / sqrt(2π * m * kB * T)
    pt.add_process(
        name='CO_adsorption',
        conditions=[Condition(coord, 'empty')],   # 前提：位点为空
        actions=[Action(coord, 'CO')],             # 结果：位点变为 CO
        rate_constant='p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)',
        tof_count={'CO_adsorption': 1},            # TOF 计数器
    )

    # 过程 2：CO 脱附（CO* → *）
    # 速率 = k_ads * exp(ΔG / kBT)   （细致平衡）
    # 注意：ΔG < 0（放热吸附）→ 脱附速率 < 吸附速率
    pt.add_process(
        name='CO_desorption',
        conditions=[Condition(coord, 'CO')],
        actions=[Action(coord, 'empty')],
        rate_constant='p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)*exp(beta*deltaG*eV)',
        tof_count={'CO_desorption': 1},
    )

    print(f"[步骤6] 过程定义完成：{[p.name for p in pt.process_list]}")

    # --- 步骤 7：打印模型摘要 ---
    print("\n--- 模型摘要 ---")
    pt.summary()

    # --- 步骤 8：保存模型 ---
    pt.filename = 'co_model.json'
    pt.save()

    return pt


# ═══════════════════════════════════════════════════════════════════════════
#  第2部分：运行 KMC 模拟
# ═══════════════════════════════════════════════════════════════════════════

def part2_run_kmc():
    """
    演示如何用 KMCEngine 运行蒙特卡洛模拟。
    """
    print("\n" + "=" * 70)
    print("  第2部分：运行 KMC 模拟")
    print("=" * 70)

    from spark import KMCEngine
    from spark.analysis import TrajectoryRecorder

    # --- 步骤 1：构建模型 ---
    pt = part1_build_simple_model()

    # --- 步骤 2：创建 KMC 引擎 ---
    # size = [Lx, Ly]，定义格点大小
    engine = KMCEngine(
        pt,
        size=[30, 30],           # 30×30 = 900 个位点
        print_rates=True,        # 打印速率常数
        banner=True,             # 打印欢迎信息
    )

    print(f"\n格点大小: {engine.lattice_size}")
    print(f"总位点数: {engine.nsites}")

    # --- 步骤 3：运行模拟 ---
    print("\n[运行] 执行 500,000 步 KMC ...")
    engine.do_steps(500000, progress=True)

    print(f"\n模拟结果：")
    print(f"  KMC 步数: {engine.kmc_step}")
    print(f"  模拟时间: {engine.kmc_time:.6e} s")

    # --- 步骤 4：获取覆盖度 ---
    coverage = engine.get_coverage()
    print(f"\n表面覆盖度：")
    for sp, theta in coverage.items():
        print(f"  θ({sp}) = {theta:.4f}")

    # --- 步骤 5：获取 TOF ---
    tof = engine.get_tof()
    print(f"\n转化频率 (TOF)：")
    for obs, val in tof.items():
        print(f"  {obs}: {val:.4e} site⁻¹·s⁻¹")

    # --- 步骤 6：查看过程统计 ---
    stats = engine.get_process_stats()
    print(f"\n各过程执行次数：")
    for name, count in stats.items():
        print(f"  {name}: {count}")

    # --- 步骤 7：查看/修改单个位点 ---
    sp_at_origin = engine.get(site_coord=(0, 0))
    print(f"\n位点 (0,0) 的物种: {sp_at_origin}")

    # 手动设置某个位点
    engine.put(site_coord=(5, 5), species_name='CO')
    print(f"手动设置 (5,5) → CO，当前物种: {engine.get((5, 5))}")

    # --- 步骤 8：使用轨迹记录器 ---
    print("\n[轨迹记录] 继续运行并记录覆盖度变化 ...")
    recorder = TrajectoryRecorder(engine)

    for i in range(50):
        engine.do_steps(10000)
        recorder.record()

    times = recorder.get_times()
    cov_CO = recorder.get_coverage_array('CO')
    print(f"  记录了 {len(times)} 个数据点")
    print(f"  CO 覆盖度范围: {cov_CO.min():.4f} ~ {cov_CO.max():.4f}")
    print(f"  CO 覆盖度均值: {cov_CO.mean():.4f}")

    # --- 步骤 9：Context Manager 用法（推荐）---
    print("\n[Context Manager] 使用 with 语句：")
    with KMCEngine(pt, size=[20, 20], print_rates=False, banner=False) as model:
        model.do_steps(100000)
        cov = model.get_coverage()
        print(f"  θ(CO) = {cov['CO']:.4f}")
    # 退出 with 后自动清理

    return engine, recorder


# ═══════════════════════════════════════════════════════════════════════════
#  第3部分：参数动态调节与覆盖度分析
# ═══════════════════════════════════════════════════════════════════════════

def part3_parameter_adjustment():
    """
    演示如何在运行时修改参数，以及与解析解（Langmuir 等温线）的对比。
    """
    print("\n" + "=" * 70)
    print("  第3部分：参数动态调节与 Langmuir 等温线验证")
    print("=" * 70)

    from spark import KMCEngine
    from spark.units import kB, eV

    pt = _build_co_project()

    # Langmuir 等温线解析解
    def langmuir(T, deltaG, p=1.0):
        K = np.exp(-deltaG * eV / (kB * T))
        return K * p / (1.0 + K * p)

    # --- 实验 1：扫描温度 ---
    print("\n--- 实验 1：固定 ΔG=-0.5 eV，扫描温度 ---")
    print(f"{'T [K]':>8s} {'θ_KMC':>8s} {'θ_Langmuir':>12s} {'误差':>8s}")
    print("-" * 40)

    for T in [400, 500, 600, 700, 800]:
        theta_exact = langmuir(T, -0.5)

        engine = KMCEngine(pt, size=[30, 30], print_rates=False, banner=False)
        engine.parameters.T = T  # 修改温度 → 速率自动更新！

        engine.do_steps(300000)   # 平衡
        samples = []
        for _ in range(10):
            engine.do_steps(10000)
            samples.append(engine.get_coverage()['CO'])
        theta_kmc = np.mean(samples)

        err = abs(theta_kmc - theta_exact) / max(theta_exact, 1e-10)
        print(f"{T:>8.0f} {theta_kmc:>8.4f} {theta_exact:>12.4f} {err:>7.2%}")

    # --- 实验 2：扫描 ΔG ---
    print("\n--- 实验 2：固定 T=600 K，扫描吸附自由能 ΔG ---")
    print(f"{'ΔG [eV]':>8s} {'θ_KMC':>8s} {'θ_Langmuir':>12s}")
    print("-" * 32)

    for dG in np.arange(-1.0, 0.4, 0.2):
        theta_exact = langmuir(600, dG)

        engine = KMCEngine(pt, size=[30, 30], print_rates=False, banner=False)
        engine.parameters.deltaG = dG

        engine.do_steps(300000)
        samples = []
        for _ in range(10):
            engine.do_steps(10000)
            samples.append(engine.get_coverage()['CO'])
        theta_kmc = np.mean(samples)

        print(f"{dG:>8.1f} {theta_kmc:>8.4f} {theta_exact:>12.4f}")

    # --- 实验 3：运行中动态修改参数 ---
    print("\n--- 实验 3：运行中动态修改参数 ---")
    engine = KMCEngine(pt, size=[30, 30], print_rates=False, banner=False)

    print("初始: T=600K, ΔG=-0.5eV")
    engine.do_steps(200000)
    print(f"  θ(CO) = {engine.get_coverage()['CO']:.4f}")

    print("修改: T=600K, ΔG=0.0eV （降低吸附强度）")
    engine.parameters.deltaG = 0.0   # 实时修改！
    engine.do_steps(200000)
    print(f"  θ(CO) = {engine.get_coverage()['CO']:.4f}")

    print("修改: T=800K, ΔG=0.0eV （升温）")
    engine.parameters.T = 800.0       # 实时修改！
    engine.do_steps(200000)
    print(f"  θ(CO) = {engine.get_coverage()['CO']:.4f}")

    # --- 查看当前参数 ---
    print(f"\n当前参数状态：")
    print(engine.parameters)

    # --- 查看当前速率 ---
    engine.print_rates()


# ═══════════════════════════════════════════════════════════════════════════
#  第4部分：构建复杂反应网络 —— Mo 表面 N₂ 还原
# ═══════════════════════════════════════════════════════════════════════════

def part4_complex_model():
    """
    演示如何构建包含多物种、多步骤的复杂反应网络。
    """
    print("\n" + "=" * 70)
    print("  第4部分：构建复杂反应网络 —— Mo 表面 N₂ 还原")
    print("=" * 70)

    from spark import Project, Site, Condition, Action

    pt = Project()
    pt.set_meta(model_name='N2_reduction_Mo', model_dimension=2)

    # --- 定义 11 个物种 ---
    species_list = [
        ('empty',   '#ffffff'),    # 空位（默认）
        ('N2',      '#0000ff'),    # 吸附态 N₂
        ('NNH',     '#0044ff'),    # 第一步加氢
        ('HNNH',    '#0088ff'),    # 交替路径中间体
        ('NNH2',    '#00ccff'),    # 远端路径中间体
        ('HNNH2',   '#00ffcc'),    # 交替路径
        ('H2NNH2',  '#00ff88'),    # 肼中间体
        ('N',       '#ff0000'),    # 原子 N
        ('NH',      '#ff4400'),    # NH
        ('NH2',     '#ff8800'),    # NH₂
        ('NH3',     '#ffcc00'),    # 产物 NH₃
    ]

    for name, color in species_list:
        pt.add_species(name=name, color=color)

    print(f"物种数量: {len(pt.species_list)}")
    print(f"物种列表: {[s.name for s in pt.species_list]}")

    # --- 晶格 ---
    layer = pt.add_layer(name='Mo_surface')
    layer.sites.append(Site(name='top', default_species='empty'))
    pt.lattice.cell = np.diag([3.2, 3.2, 15.0])

    # --- 参数 ---
    pt.add_parameter(name='T', value=300.0, adjustable=True)
    pt.add_parameter(name='p_N2', value=1.0, adjustable=True)
    pt.add_parameter(name='U', value=-0.5, adjustable=True)   # 电位 V vs RHE
    pt.add_parameter(name='A_site', value='(3.2*angstrom)**2')
    pt.add_parameter(name='E_bind_N2', value=0.8)
    pt.add_parameter(name='E_bind_NH3', value=0.5)
    pt.add_parameter(name='beta_BV', value=0.5)

    # DFT 活化能参数 [eV]
    barriers = {
        'Ea_N2_to_NNH': 2.82,       # N₂* → NNH*
        'Ea_NNH_to_HNNH': 1.14,     # NNH* → HNNH* （交替）
        'Ea_NNH_to_NNH2': 3.22,     # NNH* → NNH₂* （远端）
        'Ea_HNNH_to_HNNH2': 3.32,   # HNNH* → HNNH₂*
        'Ea_NNH2_to_N': 2.94,       # NNH₂* → N* + NH₃(g)
        'Ea_HNNH2_to_H2NNH2': 3.36, # HNNH₂* → H₂NNH₂*
        'Ea_H2NNH2_to_NH2': 5.14,   # H₂NNH₂* → NH₂*（N-N 断裂）
        'Ea_N_to_NH': 2.01,         # N* → NH*
        'Ea_NH_to_NH2': 2.68,       # NH* → NH₂*
        'Ea_NH2_to_NH3': 4.35,      # NH₂* → NH₃*
    }
    for name, val in barriers.items():
        pt.add_parameter(name=name, value=val)

    # --- 定义过程 ---
    coord = pt.lattice.generate_coord('top')

    # 辅助函数：添加电化学加氢步骤（PCET）
    def add_pcet_process(name, reactant, product, Ea_param, tof=None):
        """
        质子耦合电子转移步骤：
          k = (kB*T/h) * exp(-(Ea + β_BV * U) * eV / (kB*T))
        施加负电位 → 降低有效活化能 → 加速反应
        """
        pt.add_process(
            name=name,
            conditions=[Condition(coord, reactant)],
            actions=[Action(coord, product)],
            rate_constant=f'kB*T/h*exp(-({Ea_param} + beta_BV*U)*eV/(kB*T))',
            tof_count=tof or {},
        )

    # 辅助函数：添加热反应步骤（不受电位影响）
    def add_thermal_process(name, reactant, product, Ea_param, tof=None):
        """
        热反应步骤（如 N-N 键断裂）：
          k = (kB*T/h) * exp(-Ea * eV / (kB*T))
        """
        pt.add_process(
            name=name,
            conditions=[Condition(coord, reactant)],
            actions=[Action(coord, product)],
            rate_constant=f'kB*T/h*exp(-{Ea_param}*eV/(kB*T))',
            tof_count=tof or {},
        )

    # 1. N₂ 吸附
    pt.add_process(
        name='N2_adsorption',
        conditions=[Condition(coord, 'empty')],
        actions=[Action(coord, 'N2')],
        rate_constant='p_N2*bar*A_site/sqrt(2*pi*m_N2*umass/beta)',
    )

    # 2. N₂ 脱附
    pt.add_process(
        name='N2_desorption',
        conditions=[Condition(coord, 'N2')],
        actions=[Action(coord, 'empty')],
        rate_constant='p_N2*bar*A_site/sqrt(2*pi*m_N2*umass/beta)*exp(-beta*E_bind_N2*eV)',
    )

    # 3-12. 加氢与断键步骤
    add_pcet_process('N2_to_NNH', 'N2', 'NNH', 'Ea_N2_to_NNH',
                     tof={'N2_consumption': 1})
    add_pcet_process('NNH_to_HNNH', 'NNH', 'HNNH', 'Ea_NNH_to_HNNH')
    add_pcet_process('NNH_to_NNH2', 'NNH', 'NNH2', 'Ea_NNH_to_NNH2')
    add_pcet_process('HNNH_to_HNNH2', 'HNNH', 'HNNH2', 'Ea_HNNH_to_HNNH2')
    add_thermal_process('NNH2_to_N_NH3', 'NNH2', 'N', 'Ea_NNH2_to_N',
                        tof={'NH3_production': 1})
    add_pcet_process('HNNH2_to_H2NNH2', 'HNNH2', 'H2NNH2',
                     'Ea_HNNH2_to_H2NNH2')
    add_thermal_process('H2NNH2_to_NH2', 'H2NNH2', 'NH2',
                        'Ea_H2NNH2_to_NH2')
    add_pcet_process('N_to_NH', 'N', 'NH', 'Ea_N_to_NH')
    add_pcet_process('NH_to_NH2', 'NH', 'NH2', 'Ea_NH_to_NH2')
    add_pcet_process('NH2_to_NH3', 'NH2', 'NH3', 'Ea_NH2_to_NH3')

    # 13. NH₃ 脱附
    pt.add_process(
        name='NH3_desorption',
        conditions=[Condition(coord, 'NH3')],
        actions=[Action(coord, 'empty')],
        rate_constant='kB*T/h*exp(-E_bind_NH3*eV/(kB*T))',
        tof_count={'NH3_production': 1},
    )

    # --- 打印模型 ---
    print(f"\n过程数量: {len(pt.process_list)}")
    pt.summary()

    return pt


# ═══════════════════════════════════════════════════════════════════════════
#  第5部分：平均场微观动力学 ODE 求解
# ═══════════════════════════════════════════════════════════════════════════

def part5_microkinetic_ode():
    """
    演示平均场微观动力学模型的构建和求解。
    比 KMC 快得多，适合参数扫描和机理分析。
    """
    print("\n" + "=" * 70)
    print("  第5部分：平均场微观动力学 ODE 求解")
    print("=" * 70)

    from spark import MicroKineticModel
    from spark.rates import tst_rate, electrochemical_rate, hertz_knudsen

    # --- 步骤 1：构建模型 ---
    mkm = MicroKineticModel()

    # 添加表面物种（不需要添加 'empty'，它是隐式的）
    for sp in ['N2', 'NNH', 'HNNH', 'NNH2', 'HNNH2',
               'H2NNH2', 'N', 'NH', 'NH2', 'NH3']:
        mkm.add_species(sp)

    # 设置参数
    mkm.parameters = {
        'T': 300.0,         # 温度 [K]
        'p_N2': 1.0,        # N₂ 分压 [bar]
        'U': -1.0,          # 电位 [V vs RHE]
        'beta_BV': 0.5,     # Butler-Volmer 对称因子
        'E_bind_N2': 0.8,   # N₂ 吸附能 [eV]
        'E_bind_NH3': 0.5,  # NH₃ 吸附能 [eV]
    }

    # DFT 活化能
    Ea = {
        'N2_to_NNH': 2.82,
        'NNH_to_HNNH': 1.14,
        'NNH_to_NNH2': 3.22,
        'HNNH_to_HNNH2': 3.32,
        'NNH2_to_N': 2.94,
        'HNNH2_to_H2NNH2': 3.36,
        'H2NNH2_to_NH2': 5.14,
        'N_to_NH': 2.01,
        'NH_to_NH2': 2.68,
        'NH2_to_NH3': 4.35,
    }

    # 速率函数生成器
    def echem(key):
        """电化学步骤速率（受电位影响）"""
        def f(p):
            return electrochemical_rate(
                Ea[key], p['T'], p['U'], U0=0.0, beta_bv=p['beta_BV'])
        return f

    def thermal(key):
        """热反应速率（不受电位影响）"""
        def f(p):
            return tst_rate(Ea[key], p['T'])
        return f

    def k_ads_N2(p):
        """N₂ 吸附速率（Hertz-Knudsen）"""
        return hertz_knudsen(
            p['p_N2'] * 1e5, p['T'],
            28.014 * 1.66054e-27, (3.2e-10)**2)

    def k_des_N2(p):
        """N₂ 脱附速率"""
        return tst_rate(p['E_bind_N2'], p['T'])

    def k_des_NH3(p):
        """NH₃ 脱附速率"""
        return tst_rate(p['E_bind_NH3'], p['T'])

    # --- 步骤 2：添加反应 ---

    mkm.add_reaction('N2_adsorption',
                     reactants={'empty': 1}, products={'N2': 1},
                     rate_fwd=k_ads_N2, rate_rev=k_des_N2)

    mkm.add_reaction('N2_to_NNH',
                     reactants={'N2': 1}, products={'NNH': 1},
                     rate_fwd=echem('N2_to_NNH'),
                     tof_count={'N2_consumption': 1})

    mkm.add_reaction('NNH_to_HNNH',
                     reactants={'NNH': 1}, products={'HNNH': 1},
                     rate_fwd=echem('NNH_to_HNNH'))

    mkm.add_reaction('NNH_to_NNH2',
                     reactants={'NNH': 1}, products={'NNH2': 1},
                     rate_fwd=echem('NNH_to_NNH2'))

    mkm.add_reaction('HNNH_to_HNNH2',
                     reactants={'HNNH': 1}, products={'HNNH2': 1},
                     rate_fwd=echem('HNNH_to_HNNH2'))

    mkm.add_reaction('NNH2_to_N',
                     reactants={'NNH2': 1}, products={'N': 1},
                     rate_fwd=thermal('NNH2_to_N'),
                     tof_count={'NH3_production': 1})  # 释放 1 个 NH₃

    mkm.add_reaction('HNNH2_to_H2NNH2',
                     reactants={'HNNH2': 1}, products={'H2NNH2': 1},
                     rate_fwd=echem('HNNH2_to_H2NNH2'))

    mkm.add_reaction('H2NNH2_to_NH2',
                     reactants={'H2NNH2': 1}, products={'NH2': 1},
                     rate_fwd=thermal('H2NNH2_to_NH2'))

    mkm.add_reaction('N_to_NH',
                     reactants={'N': 1}, products={'NH': 1},
                     rate_fwd=echem('N_to_NH'))

    mkm.add_reaction('NH_to_NH2',
                     reactants={'NH': 1}, products={'NH2': 1},
                     rate_fwd=echem('NH_to_NH2'))

    mkm.add_reaction('NH2_to_NH3',
                     reactants={'NH2': 1}, products={'NH3': 1},
                     rate_fwd=echem('NH2_to_NH3'))

    mkm.add_reaction('NH3_desorption',
                     reactants={'NH3': 1}, products={'empty': 1},
                     rate_fwd=k_des_NH3,
                     tof_count={'NH3_production': 1})

    print(f"\n模型构建完成：")
    print(f"  物种: {mkm.species}")
    print(f"  反应: {[r['name'] for r in mkm.reactions]}")
    print(f"  参数: {mkm.parameters}")

    # --- 步骤 3：求解稳态 ---
    print("\n--- 3.1 稳态求解 ---")
    ss = mkm.solve_steady_state()
    mkm.print_summary(ss)

    # --- 步骤 4：ODE 瞬态求解 ---
    print("\n--- 3.2 ODE 瞬态求解 ---")
    t_eval = np.logspace(-10, 6, 300)
    sol = mkm.solve_ode((0, 1e6), t_eval=t_eval, method='BDF')

    if sol.success:
        print(f"积分成功: {len(sol.t)} 个时间点")

        # 打印几个关键时刻的覆盖度
        print(f"\n{'时间 [s]':>12s}", end='')
        for sp in ['N2', 'NNH', 'HNNH', 'N', 'NH', 'empty']:
            print(f' {sp:>10s}', end='')
        print()

        for idx in [0, 20, 50, 100, 150, 200, 250, 299]:
            if idx < len(sol.t):
                t = sol.t[idx]
                y = sol.y[:, idx]
                empty = max(1.0 - np.sum(y), 0)
                print(f'{t:>12.3e}', end='')
                for sp in ['N2', 'NNH', 'HNNH', 'N', 'NH']:
                    i = mkm.species.index(sp)
                    print(f' {max(y[i],0):>10.3e}', end='')
                print(f' {empty:>10.3e}')
    else:
        print(f"积分失败: {sol.message}")

    # --- 步骤 5：查看各反应速率 ---
    print("\n--- 3.3 各反应速率详情 ---")
    rates = mkm.get_reaction_rates(ss)
    for r in rates:
        if r['net'] != 0:
            print(f"  {r['name']:<30s} "
                  f"正向={r['rate_fwd']:.3e}  "
                  f"净速率={r['net']:.3e}")

    return mkm


# ═══════════════════════════════════════════════════════════════════════════
#  第6部分：参数扫描与速率控制度分析
# ═══════════════════════════════════════════════════════════════════════════

def part6_parameter_scan():
    """
    演示电位扫描、温度扫描、以及 Campbell 速率控制度分析。
    """
    print("\n" + "=" * 70)
    print("  第6部分：参数扫描与速率控制度分析")
    print("=" * 70)

    from models.n2_reduction_Mo import build_microkinetic_model

    mkm = build_microkinetic_model()

    # --- 6.1 电位扫描 ---
    print("\n--- 6.1 电位扫描 (T=300K, p_N2=1bar) ---")
    mkm.parameters['T'] = 300.0
    mkm.parameters['p_N2'] = 1.0

    U_values = np.linspace(-0.5, -2.0, 16)
    results = mkm.scan_parameter('U', U_values, observable='NH3_production')

    print(f"{'U [V]':>8s} {'TOF_NH3 [s⁻¹]':>15s} {'主要物种':>12s} {'覆盖度':>8s}")
    print("-" * 48)
    for i, U in enumerate(results['values']):
        cov = results['coverages'][i]
        dominant = max(cov, key=cov.get)
        tof = results['tofs'][i]
        print(f"{U:>8.2f} {tof:>15.4e} {dominant:>12s} {cov[dominant]:>8.4f}")

    # --- 6.2 温度扫描 ---
    print("\n--- 6.2 温度扫描 (U=-1.0V, p_N2=1bar) ---")
    mkm.parameters['U'] = -1.0
    T_values = np.linspace(250, 500, 11)
    results_T = mkm.scan_parameter('T', T_values, observable='NH3_production')

    print(f"{'T [K]':>8s} {'TOF_NH3 [s⁻¹]':>15s} {'主要物种':>12s}")
    print("-" * 40)
    for i, T in enumerate(results_T['values']):
        cov = results_T['coverages'][i]
        dominant = max(cov, key=cov.get)
        tof = results_T['tofs'][i]
        print(f"{T:>8.1f} {tof:>15.4e} {dominant:>12s}")

    # --- 6.3 Campbell 速率控制度 (DRC) ---
    print("\n--- 6.3 速率控制度分析 (T=300K, U=-1.0V) ---")
    print("  X_RC > 0: 加速该步骤可提高整体 TOF")
    print("  X_RC < 0: 加速该步骤会降低整体 TOF")
    print()

    mkm.parameters['T'] = 300.0
    mkm.parameters['U'] = -1.0
    ss = mkm.solve_steady_state()
    tof_check = mkm.get_tof(ss)

    if tof_check.get('NH3_production', 0) > 1e-100:
        drc = mkm.degree_of_rate_control(ss, 'NH3_production', dk=0.01)
        print(f"{'反应步骤':<30s} {'X_RC':>10s}")
        print("-" * 42)
        for name, val in sorted(drc.items(), key=lambda x: -abs(x[1])):
            if abs(val) > 0.001:
                print(f"  {name:<28s} {val:>+10.4f}")
    else:
        print("  TOF 太小，无法计算有意义的 DRC")
        print("  提示：尝试更负的电位（如 U=-5V）或更高温度")

    return results


# ═══════════════════════════════════════════════════════════════════════════
#  第7部分：结果可视化
# ═══════════════════════════════════════════════════════════════════════════

def part7_visualization():
    """
    使用 matplotlib 绘制分析结果图。
    """
    print("\n" + "=" * 70)
    print("  第7部分：结果可视化")
    print("=" * 70)

    import matplotlib
    matplotlib.use('Agg')  # 非交互式后端
    import matplotlib.pyplot as plt

    from spark import KMCEngine
    from spark.analysis import TrajectoryRecorder
    from models.n2_reduction_Mo import build_microkinetic_model

    # --- 图 1：KMC 覆盖度随时间演化 (CO 模型) ---
    print("\n绘制图1: KMC 覆盖度随时间演化 ...")
    pt = _build_co_project()
    engine = KMCEngine(pt, size=[30, 30], print_rates=False, banner=False)
    engine.parameters.deltaG = -0.3

    recorder = TrajectoryRecorder(engine)
    for _ in range(200):
        engine.do_steps(5000)
        recorder.record()

    fig, ax = plt.subplots(figsize=(8, 5))
    times = recorder.get_times()
    cov_CO = recorder.get_coverage_array('CO')
    cov_empty = recorder.get_coverage_array('empty')

    ax.plot(times, cov_CO, 'r-', label=r'$\theta_{CO}$', linewidth=2)
    ax.plot(times, cov_empty, 'b--', label=r'$\theta_{empty}$', linewidth=2)
    ax.set_xlabel('Time [s]', fontsize=12)
    ax.set_ylabel('Coverage', fontsize=12)
    ax.set_title('KMC: CO Coverage Evolution on Pd(100)\n'
                 r'T=600K, $\Delta G$=-0.3eV, p$_{CO}$=1bar', fontsize=13)
    ax.legend(fontsize=12)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig('fig1_kmc_coverage_evolution.png', dpi=150)
    print("  保存为 fig1_kmc_coverage_evolution.png")

    # --- 图 2：电位扫描 (N₂ 还原) ---
    print("\n绘制图2: 电位扫描 TOF ...")
    mkm = build_microkinetic_model()
    mkm.parameters['T'] = 300.0

    U_arr = np.linspace(-0.5, -2.0, 31)
    tof_arr = []
    cov_N2 = []
    cov_empty_arr = []

    for U in U_arr:
        mkm.parameters['U'] = float(U)
        ss = mkm.solve_steady_state()
        tof = mkm.get_tof(ss)
        tof_arr.append(tof.get('NH3_production', 0))
        cov_N2.append(ss.get('N2', 0))
        cov_empty_arr.append(1.0 - sum(ss.values()))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    ax1.semilogy(U_arr, [max(t, 1e-100) for t in tof_arr],
                 'ko-', linewidth=2, markersize=4)
    ax1.set_ylabel(r'TOF$_{NH_3}$ [s$^{-1}$]', fontsize=12)
    ax1.set_title(r'N$_2$ Reduction on Mo: Potential Scan (T=300K)',
                  fontsize=13)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(1e-100, 1e2)

    ax2.plot(U_arr, cov_N2, 'b-', label=r'$\theta_{N_2}$', linewidth=2)
    ax2.plot(U_arr, cov_empty_arr, 'k--', label=r'$\theta_{empty}$',
             linewidth=2)
    ax2.set_xlabel('U [V vs RHE]', fontsize=12)
    ax2.set_ylabel('Coverage', fontsize=12)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig('fig2_potential_scan.png', dpi=150)
    print("  保存为 fig2_potential_scan.png")

    # --- 图 3：ODE 瞬态覆盖度演化 ---
    print("\n绘制图3: ODE 瞬态覆盖度演化 ...")
    mkm.parameters['U'] = -1.0
    t_eval = np.logspace(-10, 4, 500)
    sol = mkm.solve_ode((0, 1e4), t_eval=t_eval, method='BDF')

    if sol.success:
        fig, ax = plt.subplots(figsize=(8, 5))
        colors = plt.cm.tab10(np.linspace(0, 1, len(mkm.species)))

        for i, sp in enumerate(mkm.species):
            y = np.clip(sol.y[i], 1e-30, 1)
            if np.max(y) > 1e-20:
                ax.semilogy(sol.t, y, color=colors[i],
                            label=sp, linewidth=1.5)

        # empty coverage
        empty = np.clip(1.0 - np.sum(sol.y, axis=0), 1e-30, 1)
        ax.semilogy(sol.t, empty, 'k--', label='empty', linewidth=1.5)

        ax.set_xscale('log')
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Coverage', fontsize=12)
        ax.set_title(r'ODE: Coverage Evolution (T=300K, U=-1.0V)',
                     fontsize=13)
        ax.legend(fontsize=9, ncol=3, loc='lower left')
        ax.set_ylim(1e-30, 2)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig('fig3_ode_coverage_evolution.png', dpi=150)
        print("  保存为 fig3_ode_coverage_evolution.png")

    # --- 图 4：能量图 (Reaction Energy Diagram) ---
    print("\n绘制图4: 反应能量图 ...")
    fig, ax = plt.subplots(figsize=(12, 6))

    # 交替路径
    alt_path = ['N₂*', 'NNH*', 'HNNH*', 'HNNH₂*', 'H₂NNH₂*',
                'NH₂*', 'NH₃*', '*']
    alt_barriers = [0, 2.82, 1.14, 3.32, 3.36, 5.14, 4.35, 0.5]

    # 远端路径
    dis_path = ['N₂*', 'NNH*', 'NNH₂*', 'N*+NH₃', 'NH*',
                'NH₂*', 'NH₃*', '*']
    dis_barriers = [0, 2.82, 3.22, 2.94, 2.01, 2.68, 4.35, 0.5]

    # 画阶梯图
    x = np.arange(len(alt_path))
    cum_alt = np.cumsum([0] + alt_barriers[1:])
    cum_dis = np.cumsum([0] + dis_barriers[1:])

    for i in range(len(x)):
        # 交替路径
        ax.plot([x[i]-0.3, x[i]+0.3], [cum_alt[i], cum_alt[i]],
                'b-', linewidth=2.5)
        if i < len(x) - 1:
            ax.plot([x[i]+0.3, x[i+1]-0.3],
                    [cum_alt[i], cum_alt[i+1]], 'b:', alpha=0.5)
        # 远端路径
        ax.plot([x[i]-0.3, x[i]+0.3], [cum_dis[i], cum_dis[i]],
                'r-', linewidth=2.5)
        if i < len(x) - 1:
            ax.plot([x[i]+0.3, x[i+1]-0.3],
                    [cum_dis[i], cum_dis[i+1]], 'r:', alpha=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(alt_path, fontsize=9, rotation=30, ha='right')
    ax.set_ylabel('Cumulative Barrier [eV]', fontsize=12)
    ax.set_title('N₂ Reduction on Mo: Reaction Energy Diagram', fontsize=13)
    ax.plot([], [], 'b-', linewidth=2.5, label='Alternating pathway')
    ax.plot([], [], 'r-', linewidth=2.5, label='Distal pathway')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig('fig4_energy_diagram.png', dpi=150)
    print("  保存为 fig4_energy_diagram.png")

    print("\n所有图片已保存到当前目录。")


# ═══════════════════════════════════════════════════════════════════════════
#  辅助函数
# ═══════════════════════════════════════════════════════════════════════════

def _build_co_project():
    """快速构建 CO 模型（供多个部分复用）。"""
    from spark import Project, Site, Condition, Action

    pt = Project()
    pt.set_meta(model_name='CO_on_Pd100', model_dimension=2)
    pt.add_species(name='empty')
    pt.add_species(name='CO')
    layer = pt.add_layer(name='simple_cubic')
    layer.sites.append(Site(name='hollow', default_species='empty'))
    pt.lattice.cell = np.diag([3.5, 3.5, 10.0])
    pt.add_parameter(name='T', value=600.0, adjustable=True)
    pt.add_parameter(name='p_CO', value=1.0, adjustable=True)
    pt.add_parameter(name='A', value='(3.5*angstrom)**2')
    pt.add_parameter(name='deltaG', value=-0.5, adjustable=True)
    coord = pt.lattice.generate_coord('hollow')
    pt.add_process(name='CO_adsorption',
                   conditions=[Condition(coord, 'empty')],
                   actions=[Action(coord, 'CO')],
                   rate_constant='p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)')
    pt.add_process(name='CO_desorption',
                   conditions=[Condition(coord, 'CO')],
                   actions=[Action(coord, 'empty')],
                   rate_constant='p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)*exp(beta*deltaG*eV)')
    return pt


# ═══════════════════════════════════════════════════════════════════════════
#  主入口
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='spark 完整使用教程')
    parser.add_argument('--part', type=int, default=0,
                        help='运行指定部分 (1-7)，默认运行全部')
    args = parser.parse_args()

    parts = {
        1: ('构建简单模型', part1_build_simple_model),
        2: ('运行 KMC 模拟', part2_run_kmc),
        3: ('参数调节与 Langmuir 验证', part3_parameter_adjustment),
        4: ('构建复杂反应网络', part4_complex_model),
        5: ('平均场微观动力学', part5_microkinetic_ode),
        6: ('参数扫描与速率控制度', part6_parameter_scan),
        7: ('结果可视化', part7_visualization),
    }

    if args.part > 0:
        if args.part in parts:
            name, func = parts[args.part]
            print(f"\n运行第 {args.part} 部分：{name}")
            func()
        else:
            print(f"错误：没有第 {args.part} 部分，可选 1-7")
    else:
        print("\n" + "═" * 70)
        print("  spark 完整使用教程 —— 运行全部 7 个部分")
        print("═" * 70)
        for num, (name, func) in parts.items():
            try:
                func()
            except Exception as e:
                print(f"\n第 {num} 部分出错: {e}")
                import traceback
                traceback.print_exc()

    print("\n" + "═" * 70)
    print("  教程完成！")
    print("═" * 70)


if __name__ == '__main__':
    main()
