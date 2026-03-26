# Molecular Adsorption and Desorption in Electrocatalytic KMC Simulations
# 电催化KMC模拟中的分子吸附与脱附

## Literature Review for SPARK Development

**Date:** 2026-03-22
**Scope:** Treatment of adsorption/desorption elementary steps in lattice KMC for electrocatalysis

---

## Table of Contents

1. [Adsorption in KMC](#1-adsorption-in-kmc)
   - 1.1 Non-Activated Adsorption (Barrierless)
   - 1.2 Activated Adsorption (Dissociative)
   - 1.3 Electrochemical Adsorption (PCET)
2. [Desorption in KMC](#2-desorption-in-kmc)
   - 2.1 Thermal Desorption
   - 2.2 Electrochemical Desorption
   - 2.3 Molecular vs. Dissociative Desorption
3. [Specific Challenges for Multi-Product Systems](#3-specific-challenges-for-multi-product-systems)
   - 3.1 Multiple Gas-Phase Products
   - 3.2 Solution-Phase Products
   - 3.3 Competitive Adsorption
4. [Key Implementation Details](#4-key-implementation-details)
   - 4.1 Software Comparison
   - 4.2 Rate Constant Units and Conversions
   - 4.3 Detailed Balance and Microscopic Reversibility
5. [Recommendations for SPARK](#5-recommendations-for-spark)
6. [Sources](#6-sources)

---

## 1. Adsorption in KMC
## 1. KMC中的吸附

### 1.1 Non-Activated Adsorption (Barrierless / 无势垒吸附)

#### 1.1.1 Hertz-Knudsen Equation (赫兹-克努森方程)

For non-activated (barrierless) molecular adsorption from the gas phase, the adsorption rate
is derived from kinetic gas theory via the Hertz-Knudsen equation. The molecular flux
impinging on a surface is:

```
Φ = P / √(2π m k_B T)
```

where:
- `P` = partial pressure of the adsorbate species (Pa)
- `m` = mass of one molecule (kg)
- `k_B` = Boltzmann constant (1.381 × 10⁻²³ J/K)
- `T` = gas temperature (K)

The **per-site adsorption rate constant** (units: s⁻¹) for KMC is obtained by multiplying
by the site area and sticking coefficient:

```
k_ads = (P · A_site · S) / √(2π m k_B T)
```

where:
- `A_site` = area of one adsorption site (m²), typically the unit cell area
- `S` = sticking coefficient (dimensionless, 0 < S ≤ 1)

**Key point (关键要点):** In KMC, this rate constant has units of s⁻¹ (events per site per
second). It is already a "per site" quantity because the impinging flux is multiplied by the
area of one site. The KMC algorithm selects which specific empty site the molecule lands on.

**Example from kmos/kmcos implementation:**
In the kmos software, the CO adsorption rate constant is coded as:
```python
rate_constant = 'p_CO * bar * A / sqrt(2 * pi * umass * m_CO / beta)'
```
where `A = (3.5*angstrom)**2` is the site area, `bar` converts pressure units, `umass` is
atomic mass unit, and `beta = 1/(k_B * T)`.

**Numerical example (数值举例):**
For CO adsorption at P = 1 bar, T = 600 K, on a Pt(111) site with A_site ≈ 6.67 × 10⁻²⁰ m²:
```
m_CO = 28 × 1.66 × 10⁻²⁷ = 4.65 × 10⁻²⁶ kg
√(2π m k_B T) = √(2π × 4.65e-26 × 1.381e-23 × 600) = 4.42 × 10⁻²³ kg·m/s
k_ads = (1e5 × 6.67e-20 × S) / 4.42e-23 ≈ 1.5 × 10⁸ × S   s⁻¹
```
This gives k_ads on the order of 10⁶–10⁸ s⁻¹ for typical conditions (depending on S),
consistent with literature values (Stamatakis, 2019; Reuter & Scheffler, 2006).

#### 1.1.2 Sticking Coefficient Models (粘附系数模型)

The sticking coefficient `S` accounts for the probability that an impinging molecule actually
adsorbs rather than scattering back to the gas phase:

**(a) Langmuir (direct) model (朗缪尔直接吸附模型):**
```
S(θ) = S₀ · (1 - θ)
```
- `S₀` = initial sticking coefficient on a clean surface (e.g., S₀ = 0.84 for O₂/Pt(111))
- `θ` = fractional surface coverage
- Assumes a molecule hitting an occupied site always bounces off
- In KMC, this coverage dependence is handled **automatically**: the algorithm only
  attempts adsorption on empty sites, so the factor `(1-θ)` is implicitly built in.
  You use S₀ as the rate constant prefactor.

**(b) Precursor-mediated (Kisliuk) model (前驱体介导的Kisliuk模型):**
```
S(θ) = S₀ / (1 + K · θ/(1-θ))
```
- `K` = Kisliuk parameter (ratio of desorption rate from extrinsic precursor to chemisorption
  rate from intrinsic precursor)
- When K → ∞: reduces to Langmuir behavior S = S₀(1-θ)
- When K → 0: S = S₀ regardless of coverage (strong precursor trapping)
- CO on Pt: S₀ = 0.68, well-described by Kisliuk kinetics
- The precursor model gives higher sticking at intermediate coverages because molecules
  can "hop" over occupied sites in the precursor state before finding empty sites

**In KMC (KMC中的处理):** Since KMC already handles site availability explicitly, the
Langmuir (1-θ) factor is automatic. The precursor effect is harder to capture in standard
lattice KMC. Two approaches:
1. Use a coverage-dependent S(θ) evaluated at local coverage
2. Add an explicit precursor state as an intermediate species that can diffuse and either
   chemisorb or desorb

#### 1.1.3 Adsorption from Solution Phase (溶液相吸附)

For electrochemical systems, reactants come from solution, not gas phase. The Hertz-Knudsen
equation must be replaced with a concentration-dependent expression:

```
k_ads = k°_ads · (c / c°)
```

where:
- `c` = bulk concentration of the reactant species in solution (mol/L)
- `c°` = reference concentration (typically 1 M)
- `k°_ads` = reference adsorption rate constant at c = c° (s⁻¹)

**Physical basis:** Instead of collision theory, the rate of arrival at the surface depends on:
1. Diffusion from bulk to the electrode surface (mass transport)
2. Encounter rate at the surface

In the diffusion-limited regime, the flux to the surface is:
```
J = D · c / δ
```
where `D` is the diffusion coefficient and `δ` is the diffusion layer thickness.

For a site-based KMC rate constant:
```
k_ads = (D · c · A_site) / (δ · 1 molecule)
```

**Practical approach for KMC (实用方法):** Most electrocatalytic KMC studies assume the
solution concentration at the surface is constant (well-stirred assumption, no mass transport
limitation). The adsorption rate is then:
```
k_ads(c, U) = k⁰ · (c/c°) · exp(-βeΔU/k_BT)    [if charge transfer is involved]
k_ads(c)    = k⁰ · (c/c°)                         [if purely non-Faradaic]
```

**Example for NO₃⁻ adsorption (硝酸根吸附示例):**
```
NO₃⁻(aq) + * → *NO₃ + e⁻

k_ads = k⁰_ads · (c_NO₃⁻ / c°) · exp(α·e·(U - U°) / k_BT)
```
This combines concentration dependence (1st order in NO₃⁻) with Butler-Volmer
electrochemistry (see Section 1.3).


### 1.2 Activated Adsorption (Dissociative / 活化吸附/解离吸附)

#### 1.2.1 Dissociative Adsorption Mechanism

Dissociative adsorption involves breaking a molecular bond upon adsorption:
```
A₂(g) + 2* → 2A*
```
Examples:
- O₂ + 2* → 2O*   (oxygen dissociation)
- N₂ + 2* → 2N*   (nitrogen dissociation)
- H₂ + 2* → 2H*   (hydrogen dissociation)

The rate constant includes an activation barrier E_a:
```
k_diss = (P · A_site · S₀) / √(2π m k_B T) · exp(-E_a / k_BT)
```

Or equivalently using transition state theory (TST / 过渡态理论):
```
k_diss = ν · exp(-E_a / k_BT)
```
where `ν` is an effective prefactor that includes the entropic contribution.

#### 1.2.2 Multi-Site Requirement in Lattice KMC (多位点需求)

**This is a critical KMC implementation detail (这是KMC实现的关键细节).**

In lattice KMC, dissociative adsorption requires **two adjacent empty sites**:

```
Process definition:
  Condition: site_i = empty AND site_j = empty AND (i,j) are neighbors
  Action:    site_i = A*   AND site_j = A*
  Rate:      k_diss
```

**How KMC codes handle this:**

**(a) Zacros approach:**
- All pairs of nearest-neighbor sites are enumerated based on Euclidean distance
- Each pair of adjacent empty sites constitutes a possible event
- The total propensity sums over all such valid pairs
- Example: O₂ dissociative adsorption on two empty next-nearest neighbor hollow sites

**(b) kmos approach:**
- Process definition specifies conditions on multiple sites:
```python
pt.add_process(
    name='O2_dissociative_adsorption',
    conditions=[
        Condition(coord=coord_A, species='empty'),
        Condition(coord=coord_B, species='empty')
    ],
    actions=[
        Action(coord=coord_A, species='O'),
        Action(coord=coord_B, species='O')
    ],
    rate_constant='p_O2*bar*A/(2*sqrt(2*pi*umass*m_O2/beta))*exp(-E_diss*beta)'
)
```

**(c) Geometry factor (几何因子):**
For a hexagonal (111) surface, each site has 6 nearest neighbors, so there are 3 distinct
pairs per site. The factor of 1/2 avoids double counting:
```
N_pairs = (N_sites × z) / 2
```
where z = coordination number (6 for (111), 4 for (100)).

Some codes (e.g., kmos) handle the pair counting automatically when processes are defined.
Others require the user to include a geometry factor in the rate constant.

**Important subtlety (重要细节):** For dissociative adsorption from the gas phase, the
molecule arrives with a certain impingement rate regardless of surface configuration. Only
molecules that land near a valid pair of empty sites can dissociate. The effective rate
per valid pair is:
```
k_diss_per_pair = (P · A_site · S₀ / √(2πmk_BT)) · (1/z) · exp(-E_a/k_BT)
```
The factor 1/z accounts for the fact that the molecule can dissociate onto any of z
neighbor pairs, distributing the flux.

#### 1.2.3 Lattice KMC vs. Off-Lattice (格子KMC vs 非格子)

| Aspect | Lattice KMC | Off-lattice KMC |
|--------|------------|-----------------|
| Adjacent sites | Pre-defined by lattice | Determined by neighbor lists |
| Pair enumeration | Straightforward | Requires distance cutoff |
| Coverage effects | Through pair availability | More flexible geometry |
| Code examples | Zacros, kmos, SPPARKS | MonteCoffee |
| Typical use | Flat surfaces | Nanoparticles, defects |

MonteCoffee (off-lattice, Python) uses neighbor lists to represent site connectivity rather
than mapping onto a lattice, making it well-suited for nanoparticle simulations where
surface sites have varying coordination environments.


### 1.3 Electrochemical Adsorption (PCET / 质子耦合电子转移)

#### 1.3.1 Butler-Volmer Formulation for PCET Steps

For electrochemical steps involving proton-coupled electron transfer (PCET), the rate
constant depends on the applied electrode potential U:

**Forward (cathodic) rate constant (正向/阴极速率常数):**
```
k_f(U) = k⁰ · exp[-α_c · e · (U - U°) / (k_B T)]
```

**Reverse (anodic) rate constant (逆向/阳极速率常数):**
```
k_r(U) = k⁰ · exp[α_a · e · (U - U°) / (k_B T)]
```

where:
- `k⁰` = exchange rate constant at equilibrium potential (s⁻¹)
- `α_c` = cathodic transfer coefficient (typically 0.5 for symmetric barrier)
- `α_a` = anodic transfer coefficient (α_a = 1 - α_c for single-electron transfer)
- `e` = elementary charge (1.602 × 10⁻¹⁹ C)
- `U` = applied electrode potential (V vs. reference)
- `U°` = equilibrium potential for the elementary step (V)

**Relationship to CHE (与计算氢电极的关系):**
Using the Computational Hydrogen Electrode (CHE) developed by Norskov and coworkers,
the free energy of each PCET step shifts by:
```
ΔG(U) = ΔG(U=0) + n·e·U
```
where n = number of electrons transferred.

For the activation barrier, using a Bronsted-Evans-Polanyi (BEP) relation:
```
E_a(U) = E_a(U=0) + α · e · (U - U°)
```
or equivalently:
```
E_a(U) = E_a⁰ + α · ΔG_rxn(U)
```

This gives the full KMC rate constant for an electrochemical step:
```
k(U) = ν · exp[-E_a(U) / (k_B T)]
     = ν · exp[-(E_a⁰ + α·e·(U-U°)) / (k_B T)]
     = [ν · exp(-E_a⁰/(k_BT))] · exp[-α·e·(U-U°)/(k_BT)]
     = k⁰ · exp[-α·e·η / (k_BT)]
```
where η = U - U° is the overpotential.

#### 1.3.2 Implementation in KMC (KMC中的实现)

For each electrochemical elementary step, the KMC rate constant is computed as a function
of the applied potential:

```python
# Example: Volmer step in HER
# H⁺(aq) + e⁻ + * → H*
def rate_volmer(U, T, c_H, E_a0, U0, alpha=0.5):
    """
    U: applied potential (V vs. RHE)
    T: temperature (K)
    c_H: proton concentration (mol/L), normalized to 1 M
    E_a0: activation energy at equilibrium potential (eV)
    U0: equilibrium potential for the step (V vs. RHE)
    alpha: transfer coefficient
    """
    kB = 8.617e-5  # eV/K
    eta = U - U0   # overpotential
    E_a = E_a0 + alpha * eta  # potential-dependent barrier
    # Ensure barrier doesn't go negative
    E_a = max(E_a, 0.0)
    nu = 1e13  # prefactor, s^-1
    k = nu * (c_H / 1.0) * np.exp(-E_a / (kB * T))
    return k
```

#### 1.3.3 Specific Example: NO₃⁻(aq) Adsorption (硝酸根吸附举例)

For the initial adsorption of nitrate from solution:
```
NO₃⁻(aq) + * → *NO₃ + e⁻
```

This is simultaneously:
- An adsorption event (species goes from solution to surface)
- An electrochemical event (electron transfer to the electrode)

The KMC rate constant combines both effects:
```
k_ads_NO3 = k⁰ · (c_NO₃⁻ / c°) · exp[-α·e·(U - U°_NO₃) / (k_BT)]
```

Typical parameters:
- c_NO₃⁻ ≈ 0.1–1.0 M for NO₃RR experiments
- U°_NO₃ from DFT: related to the free energy of *NO₃ formation
- α ≈ 0.5 (symmetric barrier assumption)
- k⁰ estimated from DFT barrier at equilibrium potential

#### 1.3.4 Coupled Ion-Electron Transfer (CIET) Theory

Recent work by Bazant (2023) presents a unified quantum theory of coupled ion-electron
transfer (CIET) that connects Marcus kinetics of electron transfer with Butler-Volmer
kinetics of ion transfer. In the limit of large ion transfer energies, the theory predicts
Butler-Volmer kinetics of "ion-coupled electron transfer" (ICET).

For KMC purposes, the standard Butler-Volmer form is generally sufficient, but CIET theory
provides a more rigorous treatment when both ion and electron reorganization energies are
important.


---

## 2. Desorption in KMC
## 2. KMC中的脱附

### 2.1 Thermal Desorption (热脱附)

#### 2.1.1 Arrhenius Rate Expression

For thermal (non-electrochemical) desorption:
```
k_des = ν · exp(-E_des / k_B T)
```

where:
- `ν` = pre-exponential factor (attempt frequency, s⁻¹)
- `E_des` = desorption energy (= binding energy for non-activated adsorption) (eV or kJ/mol)

**Pre-exponential factor estimation (指前因子估计):**

From transition state theory:
```
ν = (k_B T / h) · (q_TS / q_IS)
```
where q_TS and q_IS are partition functions of the transition and initial states.

| Estimation Method | ν (s⁻¹) | Notes (备注) |
|-------------------|---------|------|
| TST: k_BT/h | ~6 × 10¹² at 300K | Simple upper bound estimate |
| Typical range | 10¹² – 10¹³ | Valid for many gas/metal systems |
| Semiconductors | 10⁸ – 10¹⁸ | Much wider range observed |
| Common assumption | 10¹³ | ~60% of cases within 1 order of magnitude |

**Warning (注意):** Using 10¹³ s⁻¹ universally can introduce significant errors. The actual
prefactor depends on the entropy difference between the transition state and the adsorbed
state. For tightly bound species (low entropy adsorbed state), ν can be much larger than
10¹³. For loosely bound species (high entropy adsorbed state, approaching 2D gas), ν can
be smaller. Proper TST calculation from DFT vibrational frequencies is strongly recommended.

**Better estimate using vibrational frequencies from DFT:**
```
ν = Π(ν_i,IS) / Π(ν_j,TS)
```
where the products are over real vibrational frequencies of the initial state (IS, all modes)
and transition state (TS, all modes except the reaction coordinate).

#### 2.1.2 Coverage-Dependent Desorption Energy (覆盖度依赖的脱附能)

Lateral interactions between adsorbates modify the desorption energy:

**(a) Pairwise interaction model (对相互作用模型):**
```
E_des(config) = E_des(0) - Σ_j ε_ij · n_j
```
where:
- `E_des(0)` = desorption energy at zero coverage (isolated adsorbate)
- `ε_ij` = pairwise interaction energy between species i and neighbor j
  (ε > 0 for repulsive interactions that weaken binding)
- `n_j` = number of nearest-neighbor j-type adsorbates

**(b) Cluster Expansion Hamiltonian (团簇展开哈密顿量):**

Zacros uses a more sophisticated approach: the total energy of the adsorbate layer is
decomposed into single-body, two-body, and many-body "cluster" contributions:
```
E_total = Σ_i E_i + Σ_{i<j} J_ij · σ_i · σ_j + Σ_{i<j<k} K_ijk · σ_i · σ_j · σ_k + ...
```

The **desorption barrier for a specific local environment** is then:
```
E_a,des(config) = E_a,des(isolated) + ΔE_lateral(config)
```

where ΔE_lateral accounts for how the specific arrangement of neighbors shifts the
energy of the initial state (and possibly the transition state).

**Effect on KMC rates (对KMC速率的影响):**
- Repulsive lateral interactions (ε > 0): lower E_des → faster desorption
- Attractive lateral interactions (ε < 0): higher E_des → slower desorption

In the cluster expansion approach, the rate constant for each individual desorption event
depends on the **local configuration** around the desorbing molecule. This means rates
must be recalculated whenever the local environment changes (when a neighbor adsorbs,
desorbs, or diffuses).

**Typical lateral interaction strengths (典型横向相互作用强度):**
| System | ε (eV) | Type |
|--------|--------|------|
| CO-CO on Pd(111) | +0.19 | Repulsive (nearest neighbor) |
| O-O on Pt(111) | +0.30 | Repulsive |
| CO-O on Ru(0001) | −0.10 | Attractive |
| H-H on Pt(111) | +0.04 | Weakly repulsive |


### 2.2 Electrochemical Desorption (电化学脱附)

#### 2.2.1 Potential-Dependent Desorption

For desorption steps involving electron transfer:
```
*A + H⁺ + e⁻ → AH(aq)   or   *A → A⁺(aq) + e⁻
```

The rate constant has Butler-Volmer potential dependence:
```
k_des(U) = ν · exp[-(E_des⁰ + (1-α)·e·(U - U°)) / (k_BT)]
```

Note the sign: for a cathodic (reductive) desorption, more negative potential increases
the rate. For an anodic (oxidative) desorption, more positive potential increases the rate.

#### 2.2.2 Desorption vs. Reaction: The OH* Example (脱附vs反应: OH*的例子)

Consider:
```
*OH + H⁺ + e⁻ → H₂O(l) + *
```

**Is this desorption or reaction?**

From the KMC perspective, this is a **surface reaction that produces a desorbed product**.
It should be treated as:
- A single elementary step
- Involving one surface species (*OH)
- With a potential-dependent rate constant (Butler-Volmer)
- The product (H₂O) leaves the surface, freeing the site

In KMC implementation:
```python
process = Process(
    name='OH_reduction',
    conditions=[Condition(coord=site, species='OH')],
    actions=[Action(coord=site, species='empty')],
    rate_constant='nu * exp(-(Ea_OH + (1-alpha)*e*(U-U0_OH)) / (kB*T))'
)
```

**The distinction matters for bookkeeping (区分对于簿记很重要):**
- Count this as a Faradaic event (contributes to current)
- Count H₂O as a product
- The site becomes free for subsequent adsorption

#### 2.2.3 Product Release Steps That Are Purely Thermal

Some product desorption steps do not involve electron transfer:
```
*NH₃ → NH₃(aq) + *     (thermal, non-Faradaic)
*N₂O → N₂O(g) + *      (thermal, non-Faradaic)
```

These are treated as standard thermal desorption (Section 2.1) with no potential dependence.
The rate constant is purely Arrhenius:
```
k_des = ν · exp(-E_des / k_BT)
```


### 2.3 Molecular vs. Dissociative Desorption (分子vs解离脱附)

#### 2.3.1 Associative (Recombinative) Desorption (结合脱附)

The reverse of dissociative adsorption:
```
2H* → H₂(g) + 2*     (Tafel step in HER)
2N* → N₂(g) + 2*
2CO* → (CO)₂ is not typical; but CO* + O* → CO₂(g) + 2* is
```

**KMC implementation (KMC实现):**

Associative desorption requires **two adjacent occupied sites** with the correct species:

```python
process = Process(
    name='Tafel_H2_desorption',
    conditions=[
        Condition(coord=coord_A, species='H'),
        Condition(coord=coord_B, species='H'),
        # coord_A and coord_B must be neighbors
    ],
    actions=[
        Action(coord=coord_A, species='empty'),
        Action(coord=coord_B, species='empty'),
    ],
    rate_constant='nu * exp(-E_a_Tafel / (kB * T))'
)
```

**Critical requirement:** Both sites must be occupied by the correct species. KMC naturally
handles this — the event is only possible when the local configuration is satisfied.

**Rate expression considerations:**
- The rate constant is per valid pair (each pair of adjacent H* atoms)
- Barrier E_a_Tafel from DFT NEB calculation
- For the Tafel step: E_a typically 0.5–0.9 eV on various metals
- The associative desorption barrier includes both the recombination barrier and the
  desorption energy

#### 2.3.2 Electrochemical Associative Desorption (Heyrovsky Step)

```
H* + H⁺(aq) + e⁻ → H₂(g) + *     (Heyrovsky step in HER)
```

This is a special case: one species from the surface + one from solution + electron transfer.

KMC rate constant:
```
k_Hey(U) = ν · (c_H⁺/c°) · exp[-(E_a⁰ + α·e·(U-U°)) / (k_BT)]
```

Implementation:
```python
process = Process(
    name='Heyrovsky',
    conditions=[Condition(coord=site, species='H')],
    actions=[Action(coord=site, species='empty')],
    rate_constant='nu * (c_H/c_ref) * exp(-(Ea_Hey + alpha*e*(U-U0)) / (kB*T))'
)
```
Note: only one surface site is needed (the H*). The proton comes from solution (captured
in the concentration-dependent prefactor).

#### 2.3.3 Detailed Balance: Linking Adsorption and Desorption

For the adsorption/desorption pair:
```
A(g) + * ⇌ A*
```

Microscopic reversibility requires:
```
k_ads / k_des = K_eq = exp(-ΔG_ads / k_BT)
```

This means if you set k_ads from the Hertz-Knudsen equation, then k_des is **not** a free
parameter — it must be:
```
k_des = k_ads / K_eq = k_ads · exp(ΔG_ads / k_BT)
```

In practice, you can set either:
1. k_ads from Hertz-Knudsen and compute k_des from detailed balance, or
2. k_des from Arrhenius with DFT-computed binding energy, and compute k_ads from
   detailed balance

**Zacros enforces this automatically** through its cluster expansion Hamiltonian: the
energetics model defines the thermodynamic landscape, and rate constants for forward
and reverse processes are constructed to be thermodynamically consistent.


---

## 3. Specific Challenges for Multi-Product Systems
## 3. 多产物体系的特殊挑战

### 3.1 Multiple Gas-Phase Products (多种气相产物)

In NO₃RR, multiple products can form:
```
NO₃⁻ → NH₃ (desired, 8e⁻)
NO₃⁻ → N₂ (undesired, 10e⁻ total for 2 NO₃⁻)
NO₃⁻ → N₂O (undesired, 8e⁻ total for 2 NO₃⁻)
2H⁺ + 2e⁻ → H₂ (competing HER)
```

#### 3.1.1 Tracking Production Rates

In KMC, each desorption event is logged with:
- Product identity (NH₃, N₂, N₂O, H₂, etc.)
- Time of event
- Number of electrons transferred in the pathway

**Turnover Frequency (TOF / 转换频率) for each product:**
```
TOF_i = N_des_i / (N_sites · t_total)
```
where N_des_i = number of desorption events producing species i.

**Instantaneous production rate:**
Use a time window Δt and count desorption events within that window:
```
r_i(t) = ΔN_des_i / (N_sites · Δt)
```

#### 3.1.2 Partial Current Density from KMC (从KMC计算分电流密度)

The partial current density for product i:
```
j_i = (n_i · e · N_des_i) / (A_total · t_total)
```
where:
- n_i = number of electrons transferred per molecule of product i
- A_total = total surface area simulated = N_sites × A_site

**Faradaic Efficiency (法拉第效率):**
```
FE_i = j_i / j_total = (n_i · N_des_i) / Σ_k (n_k · N_des_k)
```

**Example for NO₃RR products (NO₃RR产物示例):**

| Product | n_i (electrons) | If 100 events each | j_i contribution |
|---------|----------------|---------------------|------------------|
| NH₃     | 8              | 800e                | 800e/(A·t)       |
| N₂      | 10 (for 2 NO₃⁻)| 1000e               | 1000e/(A·t)      |
| N₂O     | 8 (for 2 NO₃⁻) | 800e                | 800e/(A·t)       |
| H₂      | 2              | 200e                | 200e/(A·t)       |

### 3.2 Solution-Phase Products (溶液相产物)

#### 3.2.1 NH₃(aq) vs NH₃(g): Does It Matter?

In aqueous electrolysis, NH₃ dissolves in solution. For KMC:

**In most cases, it does NOT matter** whether the product is gas-phase or stays dissolved:
- The KMC simulation models surface kinetics
- Once a product desorbs from the surface, it is "gone" from the KMC perspective
- The site becomes available for new reactions
- Whether NH₃ then dissolves or escapes to gas phase is a mass-transport question
  outside the scope of surface KMC

**When it DOES matter:**
- If a dissolved product can **re-adsorb** (e.g., NH₃ readsorption and further oxidation)
- If solution-phase concentration buildup affects thermodynamics
- If pH changes locally due to product formation (NH₃ is basic: NH₃ + H₂O → NH₄⁺ + OH⁻)

**Practical treatment (实际处理方法):**
1. Treat product desorption as irreversible (product leaves and doesn't return)
2. If re-adsorption is important, include reverse adsorption with a rate proportional
   to the local product concentration
3. For NH₃ in acidic media, the protonation NH₃* + H⁺ → NH₄⁺(aq) + * is often barrierless
   and fast, making re-adsorption negligible

#### 3.2.2 Treatment of NO₂⁻ as a Dissolved Intermediate

In NO₃RR, NO₂⁻ can desorb into solution and re-adsorb:
```
*NO₂ → NO₂⁻(aq) + *    (desorption)
NO₂⁻(aq) + * → *NO₂    (re-adsorption)
```

This creates a "shuttle" mechanism. In KMC, you need:
- A variable tracking solution-phase NO₂⁻ concentration
- Re-adsorption rate proportional to [NO₂⁻]
- Update [NO₂⁻] after each desorption/adsorption event

This is analogous to the "reservoir" concept in grand-canonical KMC.


### 3.3 Competitive Adsorption (竞争吸附)

#### 3.3.1 Multiple Species Competing for Same Sites

In NO₃RR, several species compete for surface sites:
- NO₃⁻, NO₂⁻ (reactant/intermediate adsorption)
- H⁺ / H₂O (for hydrogen supply)
- Reaction intermediates (*NO, *N, *NH, *NH₂, *NH₃, *OH, etc.)

**KMC's natural advantage (KMC的天然优势):**
Unlike mean-field microkinetic models that use average coverages, KMC naturally handles
competitive adsorption because:
1. Each site is explicitly tracked as occupied or empty
2. Adsorption only occurs on available (empty) sites
3. The stochastic selection automatically weights by the number of available sites
4. Spatial correlations (clustering, island formation) are captured

#### 3.3.2 Mean-Field vs KMC for Competitive Adsorption

In mean-field models, competitive Langmuir adsorption gives:
```
θ_A = K_A · P_A / (1 + K_A · P_A + K_B · P_B)
θ_B = K_B · P_B / (1 + K_A · P_A + K_B · P_B)
```

This assumes **random mixing** of adsorbates. KMC goes beyond this:

| Phenomenon | Mean-Field | KMC |
|-----------|-----------|-----|
| Average coverage | Correct | Correct |
| Spatial correlations | Ignored | Captured |
| Island formation | Cannot capture | Naturally captured |
| Site blocking by multi-dentate species | Approximated | Exact |
| Local coverage effects on rates | Uses average θ | Uses actual local config |
| Phase transitions in adsorbate layer | Discontinuous | Smooth/realistic |

#### 3.3.3 Site Blocking Effects (位点阻塞效应)

Multi-dentate adsorbates (species that occupy multiple sites) create interesting
competitive effects:

Example: On a Pt surface, if C and CH species adsorb on hollow sites, they block all four
surrounding Pt sites from adsorbing other species. This creates:
- Larger "exclusion zones" around big adsorbates
- Non-trivial packing effects at high coverage
- Coverage-dependent selectivity changes

In KMC, this is handled by defining multi-site species:
```python
# Example: *NO₃ occupying a top site plus blocking 3 neighbor bridge sites
process = Process(
    name='NO3_adsorption',
    conditions=[
        Condition(coord=top_site, species='empty'),
        Condition(coord=bridge_1, species='empty'),
        Condition(coord=bridge_2, species='empty'),
        Condition(coord=bridge_3, species='empty'),
    ],
    actions=[
        Action(coord=top_site, species='NO3'),
        Action(coord=bridge_1, species='blocked'),
        Action(coord=bridge_2, species='blocked'),
        Action(coord=bridge_3, species='blocked'),
    ],
    rate_constant='k_ads_NO3'
)
```

#### 3.3.4 Kinetic Selectivity vs Thermodynamic Selectivity (动力学选择性vs热力学选择性)

KMC reveals that the species with weaker binding but faster adsorption rate can initially
dominate the surface coverage (kinetic selectivity), before eventually being displaced by
the more strongly binding species at equilibrium (thermodynamic selectivity).

This is particularly important at early reaction times or under transient conditions, and
KMC is the right tool to capture these effects.


---

## 4. Key Implementation Details
## 4. 关键实现细节

### 4.1 Software Comparison (软件比较)

| Feature | Zacros | kmos/kmcos | SPPARKS | MonteCoffee |
|---------|--------|-----------|---------|-------------|
| **Language** | Fortran | Python/Fortran | C++ | Python |
| **Lattice type** | On-lattice (graph-theoretical) | On-lattice | On-lattice | Off-lattice (neighbor lists) |
| **Lateral interactions** | Cluster expansion Hamiltonian | BEP + on-the-fly | Pairwise | User-defined |
| **Adsorption** | Non-activated + activated | Non-activated + activated | User-defined events | User-defined |
| **Desorption** | Arrhenius with cluster expansion | Arrhenius with BEP corrections | User-defined events | User-defined |
| **Multi-site processes** | Yes (graph patterns) | Yes (multi-coord conditions) | Yes (2-site reactions) | Yes (neighbor lists) |
| **Electrochemistry** | Not built-in | Not built-in | Not built-in | Not built-in |
| **Algorithm** | Variable step size method (VSSM) | First-reaction method (FRM) or VSSM | Rejection KMC + KMC | First-reaction method |
| **Parallelization** | Time-warp (optimistic) | Limited | Sublattice algorithm | Limited |
| **Best for** | Detailed surface chemistry | Rapid prototyping, DFT integration | Large-scale parallel | Nanoparticles |

**Key observation for SPARK (关键观察):** None of the existing KMC codes have
built-in support for electrochemistry (potential-dependent rates, Butler-Volmer kinetics).
This is a major opportunity for SPARK to fill this gap.

#### 4.1.1 How Zacros Handles Adsorption/Desorption

Zacros defines all elementary steps in a `mechanism_input.dat` file:
```
reversible_step CO_ads
  sites 1
  initial
    1 *          # empty site
  final
    1 CO*        # CO adsorbed
  site_types 1 top
  pre_expon  7.6e7    # ≈ P·A/√(2πmkT), s⁻¹
  activ_eng  0.0      # non-activated (barrierless)
end_reversible_step
```

For dissociative adsorption:
```
reversible_step O2_diss_ads
  sites 2
  neighboring 1-2
  initial
    1 *    2 *    # two adjacent empty sites
  final
    1 O*   2 O*   # two O atoms
  site_types 1 fcc  2 fcc
  pre_expon  2.5e8
  activ_eng  0.26     # activation barrier in eV
end_reversible_step
```

Zacros automatically computes the reverse rate from detailed balance using the
energetics model (cluster expansion).

#### 4.1.2 How kmos Handles Adsorption/Desorption

In kmos, the user defines processes with explicit rate constant expressions:
```python
# CO adsorption
pt.add_process(
    name='CO_adsorption',
    conditions=[Condition(coord=coord, species='empty')],
    actions=[Action(coord=coord, species='CO')],
    rate_constant='p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)'
)

# CO desorption (reverse, using detailed balance)
pt.add_process(
    name='CO_desorption',
    conditions=[Condition(coord=coord, species='CO')],
    actions=[Action(coord=coord, species='empty')],
    rate_constant='p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)*exp(beta*deltaG)'
)
```

The `beta*deltaG` term in the desorption rate ensures detailed balance:
k_des = k_ads × exp(ΔG/k_BT), where ΔG > 0 for desorption.


### 4.2 Rate Constant Units and Conversions (速率常数单位与转换)

#### 4.2.1 KMC Rate Constant Units

**All KMC rate constants must have units of s⁻¹ (per site per second).**

This is because in the Variable Step Size Method (VSSM / Gillespie algorithm):
```
Total propensity: R_total = Σ_i k_i · N_i
```
where k_i is in s⁻¹ and N_i is the number of sites where event i can occur.
The time step is: Δt = -ln(u) / R_total, where u ∈ (0,1) is a random number.

#### 4.2.2 Conversion from Different Rate Expressions

**(a) From pressure-dependent gas-phase adsorption rate (从气相吸附速率转换):**
```
r_ads = P·S / √(2πmk_BT)                    [m⁻²·s⁻¹]  (flux)
k_ads = P·S·A_site / √(2πmk_BT)             [s⁻¹]       (per site)
```

**(b) From TST-derived rate constants (从TST导出的速率常数转换):**
```
k_TST = (k_BT/h) · exp(-ΔG‡/k_BT)           [s⁻¹]       (already per site)
```

**(c) From experimental TOF (从实验TOF转换):**
```
k = TOF_experimental                          [s⁻¹]       (if per active site)
```

**(d) From bimolecular rate constants (从双分子速率常数转换):**
For a Langmuir-Hinshelwood reaction A* + B* → product:
- Mean-field: rate = k_LH · θ_A · θ_B · N_sites
- KMC: k_LH (s⁻¹) is the rate per valid A*-B* pair; the number of such pairs is
  counted automatically by KMC

**(e) From solution-phase concentration (从溶液相浓度转换):**
```
k_ads = k° · (c/c°)                          [s⁻¹]       (1st order in concentration)
```
where k° incorporates the mass-transport-limited flux.

#### 4.2.3 Common Pitfalls (常见陷阱)

1. **Forgetting the site area factor:** The Hertz-Knudsen flux is per m²; you must multiply
   by A_site to get a per-site rate.

2. **Double-counting in pair processes:** For A* + B* → products where A and B are on
   neighbor pairs, make sure to count each pair once, not twice.

3. **Mixing units:** Ensure all energies are in consistent units (eV with k_B in eV/K,
   or J with k_B in J/K). A common source of bugs.

4. **Pressure-to-concentration conversion:** For electrochemistry,
   P/(k_BT) gives number density (m⁻³). For solution: c(mol/L) × N_A gives m⁻³.

5. **Forgetting concentration dependence:** For solution-phase reactants, the adsorption
   rate must be proportional to concentration. This is NOT automatic in KMC codes
   designed for gas-phase catalysis.


### 4.3 Detailed Balance and Microscopic Reversibility (细致平衡与微观可逆性)

#### 4.3.1 Requirements

For every pair of forward and reverse elementary steps:
```
k_forward / k_reverse = K_eq = exp(-ΔG_rxn / k_BT)
```

This is **mandatory** for thermodynamic consistency. If violated, the KMC simulation will
not reach the correct equilibrium state (if isolated) and may produce unphysical results
under steady-state conditions.

#### 4.3.2 How to Ensure It in Practice

**Method 1: Set both rates independently from DFT, then check consistency**
- Calculate E_a,forward and E_a,reverse from NEB
- Verify: E_a,forward - E_a,reverse = ΔE_rxn (to within DFT accuracy)
- Adjust one barrier if inconsistent

**Method 2: Set one rate from physical model, derive the other (推荐方法)**
- Set k_ads from Hertz-Knudsen equation
- Set ΔG_ads from DFT
- Compute k_des = k_ads × exp(ΔG_ads / k_BT)

**Method 3: Use energetics-based approach (Zacros style)**
- Define the energy landscape (initial state, transition state, final state energies)
- Both forward and reverse rate constants are computed from the same energy landscape
- Detailed balance is automatically satisfied

#### 4.3.3 Special Considerations for Electrochemistry

For electrochemical steps, the equilibrium constant depends on potential:
```
K_eq(U) = exp[-(ΔG°_rxn + n·e·U) / (k_BT)]
```

The forward and reverse rate constants must satisfy:
```
k_f(U) / k_r(U) = K_eq(U)
```

If k_f = k⁰_f · exp[-α·e·η/(k_BT)] and k_r = k⁰_r · exp[(1-α)·e·η/(k_BT)], then:
```
k_f / k_r = (k⁰_f/k⁰_r) · exp[-e·η/(k_BT)]
```
which is consistent with K_eq(U) only if k⁰_f/k⁰_r = K_eq(U°).

**This is the fundamental constraint (这是根本性约束):** once you set the forward rate
constant and the thermodynamics (from DFT/CHE), the reverse rate constant is fixed.

#### 4.3.4 Cycle Consistency

In a multi-step mechanism like NO₃RR, you may have thermodynamic cycles. For example:
```
*NO₃ → *NO₂ + *O    (dissociation on surface)
*NO₂ → NO₂⁻(aq)     (desorption)
NO₂⁻(aq) → *NO₂     (re-adsorption)
```

The free energies around any such cycle must sum to zero. This is an additional constraint
beyond pairwise detailed balance.


---

## 5. Recommendations for SPARK Implementation
## 5. SPARK实现建议

Based on this literature review, here are specific design recommendations:

### 5.1 Adsorption Module Design

```python
class AdsorptionEvent:
    """Base class for all adsorption events in SPARK"""

    # Type 1: Non-activated gas-phase adsorption
    # k = P * A_site * S / sqrt(2*pi*m*kB*T)
    GAS_PHASE_MOLECULAR = 'gas_molecular'

    # Type 2: Activated (dissociative) gas-phase adsorption
    # k = P * A_site * S / sqrt(2*pi*m*kB*T) * exp(-Ea/kBT)
    # Requires 2+ adjacent empty sites
    GAS_PHASE_DISSOCIATIVE = 'gas_dissociative'

    # Type 3: Solution-phase non-Faradaic adsorption
    # k = k0 * (c/c_ref)
    SOLUTION_MOLECULAR = 'solution_molecular'

    # Type 4: Electrochemical adsorption (PCET)
    # k = k0 * (c/c_ref) * exp(-alpha*e*(U-U0)/(kB*T))
    ELECTROCHEMICAL = 'electrochemical'
```

### 5.2 Desorption Module Design

```python
class DesorptionEvent:
    """Base class for all desorption events"""

    # Type 1: Thermal desorption (non-Faradaic)
    # k = nu * exp(-E_des(config) / kBT)
    THERMAL = 'thermal'

    # Type 2: Electrochemical desorption (Faradaic)
    # k = nu * exp(-(Ea + (1-alpha)*e*(U-U0)) / kBT)
    ELECTROCHEMICAL = 'electrochemical'

    # Type 3: Associative desorption (2 surface species → 1 gas product)
    # k = nu * exp(-Ea / kBT), requires adjacent pair
    ASSOCIATIVE = 'associative'

    # Type 4: Electrochemical associative (Heyrovsky-type)
    # k = nu * (c/c_ref) * exp(-(Ea + alpha*e*(U-U0)) / kBT)
    ELECTROCHEMICAL_ASSOCIATIVE = 'electrochemical_associative'
```

### 5.3 Product Tracking System

```python
class ProductTracker:
    """Track desorption events to calculate FE, TOF, partial current density"""

    def record_desorption(self, product_name, n_electrons, time, site_id):
        """Record a desorption event"""
        self.events[product_name].append({
            'n_electrons': n_electrons,
            'time': time,
            'site_id': site_id
        })

    def get_faradaic_efficiency(self, product_name):
        """FE_i = n_i * N_i / Σ(n_k * N_k)"""
        total_charge = sum(
            p['n_electrons'] * len(self.events[p_name])
            for p_name in self.events
        )
        product_charge = (
            self.events[product_name][0]['n_electrons']
            * len(self.events[product_name])
        )
        return product_charge / total_charge

    def get_partial_current_density(self, product_name, total_area, total_time):
        """j_i = n_i * e * N_i / (A * t)"""
        e = 1.602e-19  # C
        n_events = len(self.events[product_name])
        n_e = self.events[product_name][0]['n_electrons']
        return n_e * e * n_events / (total_area * total_time)
```

### 5.4 Detailed Balance Enforcement

```python
class ThermodynamicConsistency:
    """Enforce detailed balance across all elementary steps"""

    def compute_reverse_rate(self, k_forward, delta_G, T):
        """Given forward rate and free energy change, compute reverse rate"""
        kB = 8.617e-5  # eV/K
        return k_forward * np.exp(delta_G / (kB * T))

    def check_cycle_consistency(self, cycle_delta_Gs):
        """Verify that free energies around any cycle sum to zero"""
        total = sum(cycle_delta_Gs)
        if abs(total) > 0.01:  # tolerance in eV
            warnings.warn(
                f"Thermodynamic cycle inconsistency: "
                f"ΣΔG = {total:.4f} eV (should be 0)"
            )
```

### 5.5 Key Advantages of SPARK Over Existing Codes

1. **Built-in electrochemistry support** — None of the existing codes (Zacros, kmos,
   SPPARKS, MonteCoffee) natively support potential-dependent rate constants. SPARK
   should make Butler-Volmer / CHE-based rates a first-class feature.

2. **Concentration-dependent adsorption** — Solution-phase reactant adsorption with
   concentration dependence should be as easy to specify as gas-phase adsorption.

3. **Product tracking and FE calculation** — Built-in tools for tracking desorption events,
   computing Faradaic efficiency, partial current density, and selectivity.

4. **Thermodynamic consistency enforcement** — Automatic checking and optional enforcement
   of detailed balance and cycle consistency.


---

## 6. Sources

### Key Review Articles
- [A Practical Guide to Surface Kinetic Monte Carlo Simulations (Stamatakis, 2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6465329/)
- [Kinetic Monte Carlo simulations for heterogeneous catalysis: Fundamentals, current status, and challenges (2022)](https://pubs.aip.org/aip/jcp/article/156/12/120902/2840948/Kinetic-Monte-Carlo-simulations-for-heterogeneous)
- [DFT-Based Multiscale Modeling of Heterogeneous (Electro)Catalytic Reactions (2025)](https://pubs.acs.org/doi/10.1021/acscatal.5c07967)
- [Microkinetic modeling in electrocatalysis: Applications, limitations, and recommendations](https://www.sciencedirect.com/science/article/abs/pii/S0021951721003523)

### Adsorption Kinetics
- [Kinetics of Adsorption - Chemistry LibreTexts](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Surface_Science_(Nix)/02:_Adsorption_of_Molecules_on_Surfaces/2.03:_Kinetics_of_Adsorption)
- [Hertz-Knudsen equation - Wikipedia](https://en.wikipedia.org/wiki/Hertz%E2%80%93Knudsen_equation)
- [Sticking coefficient - Wikipedia](https://en.wikipedia.org/wiki/Sticking_coefficient)
- [The Role of Precursor States in Adsorption, Surface Reactions and Catalysis](https://link.springer.com/article/10.1007/s11244-016-0538-6)

### First-Principles KMC
- [First-principles kinetic Monte Carlo simulations: CO oxidation at RuO₂(110) (Reuter & Scheffler, 2006)](https://link.aps.org/doi/10.1103/PhysRevB.73.045433)
- [CO Oxidation on Pd(111): A First-Principles-Based Kinetic Monte Carlo Study](https://pubs.acs.org/doi/10.1021/cs500377j)
- [H2 Thermal Desorption Spectra on Pt(111): A DFT and KMC Study](https://www.mdpi.com/2073-4344/8/10/450)

### Lateral Interactions and Cluster Expansion
- [Efficient Implementation of Cluster Expansion Models in Surface KMC Simulations](https://researchgate.net/publication/335079263_Efficient_Implementation_of_Cluster_Expansion_Models_in_Surface_Kinetic_Monte_Carlo_Simulations_with_Lateral_Interactions_Subtraction_Schemes_Supersites_and_the_Supercluster_Contraction)
- [Speeding up the Detection of Adsorbate Lateral Interactions in Graph-Theoretical KMC](https://pubs.acs.org/doi/10.1021/acs.jpca.3c05581)
- [Parallel kinetic Monte Carlo simulation framework with accurate lateral interactions](https://pubs.aip.org/aip/jcp/article/139/22/224706/193591/Parallel-kinetic-Monte-Carlo-simulation-framework)

### Electrochemical Kinetics
- [Unified quantum theory of electrochemical kinetics by coupled ion-electron transfer (Bazant, 2023)](https://pubs.rsc.org/en/content/articlehtml/2023/fd/d3fd00108c)
- [Butler-Volmer equation - Wikipedia](https://en.wikipedia.org/wiki/Butler%E2%80%93Volmer_equation)
- [Butler-Volmer equation lecture (MIT OCW, Bazant)](https://ocw.mit.edu/courses/10-626-electrochemical-energy-systems-spring-2014/56cfa6e0f28bc8fc1a647cbe679384d1_MIT10_626S14_S11lec13.pdf)
- [A simple method to approximate electrode potential-dependent activation energies using DFT](https://www.sciencedirect.com/science/article/abs/pii/S0920586117300597)
- [Using BEP relations to predict electrode potential-dependent activation energies](https://www.sciencedirect.com/science/article/abs/pii/S0920586118303122)

### BEP Relations
- [BEP relation for hydrogen evolution reaction from first-principles](https://www.nature.com/articles/s41524-024-01244-3)
- [Direct Demonstration of Unified BEP Relationships for PCET Reactions on Transition Metal Surfaces](https://www.osti.gov/biblio/1851368)
- [Brønsted-Evans-Polanyi Relation of Multistep Reactions and Volcano Curves](https://pubs.acs.org/doi/10.1021/jp711191j)

### Pre-exponential Factors
- [Estimating Pre-Exponential Factors for Desorption from Semiconductors](https://www.sciencedirect.com/science/article/abs/pii/S0169433201003828)
- [Adsorption and desorption equilibria from statistical thermodynamics and rates from TST](https://acp.copernicus.org/articles/21/15725/2021/)

### KMC Software
- [Zacros - Home](https://www.zacros.org/home)
- [kmos documentation](https://kmos.readthedocs.io/en/latest/tutorials/first_model_api.html)
- [SPPARKS Kinetic Monte Carlo Simulator](https://spparks.github.io/)
- [MonteCoffee: A programmable kinetic Monte Carlo framework](https://research.chalmers.se/publication/506416/file/506416_Fulltext.pdf)
- [Zacros: Mapping DFT Energies to Input](https://zacros.org/resources/tutorials/10-tutorial-4-mapping-dft-energies-to-zacros-input?limitstart=&showall=1)

### NO₃RR Specific
- [Investigating High-Performance Non-Precious TMO Catalysts for NRR: DFT-kMC-LSTM Approach](https://pubs.acs.org/doi/10.1021/acscatal.3c01360)
- [Recent advances in mechanistic studies for electrochemical NO₃RR to ammonia](https://www.nature.com/articles/s42004-025-01864-w)
- [Electrocatalytic nitrate to ammonia conversion: mechanistic insights to applications](https://link.springer.com/article/10.1007/s44422-025-00010-w)
- [Kinetics of NH₃ Desorption and Diffusion on Pt: Implications for the Ostwald Process](https://pubs.acs.org/doi/10.1021/jacs.1c09269)

### Selectivity and Faradaic Efficiency
- [Reliable reporting of Faradaic efficiencies for electrocatalysis research](https://www.nature.com/articles/s41467-023-36880-8)
- [The Significance of Properly Reporting Turnover Frequency in Electrocatalysis Research](https://pmc.ncbi.nlm.nih.gov/articles/PMC8596788/)
- [Achieving Theory-Experiment Parity for Activity and Selectivity in Heterogeneous Catalysis Using Microkinetic Modeling](https://pmc.ncbi.nlm.nih.gov/articles/PMC9069691/)

### Competitive Adsorption and Multi-product Systems
- [Multiscale Investigation of the Mechanism and Selectivity of CO₂ Hydrogenation over Rh(111)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11036393/)
- [Kinetic Understanding of Catalytic Selectivity and Product Distribution of Electrochemical CO₂RR](https://pmc.ncbi.nlm.nih.gov/articles/PMC10052237/)
- [pH Effects in a Model Electrocatalytic Reaction Disentangled](https://pubs.acs.org/doi/10.1021/jacsau.2c00662)
