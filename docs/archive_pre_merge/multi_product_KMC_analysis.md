# Multi-Product KMC Analysis Methodology for Electrocatalysis

## Literature Review for SPARK Development

---

## 1. Multi-Product Selectivity Analysis in KMC

### 1.1 Tracking Product Formation in KMC

In lattice KMC simulations for heterogeneous catalysis, product formation is tracked by counting **desorption events** for each product species. Every time a product molecule (e.g., CH4, C2H4, CO, H2, NH3, N2, H2O2) desorbs from the surface, the event is logged with a timestamp and the identity of the desorbing species.

**Turnover Frequency (TOF) Calculation:**

The turnover frequency is defined as the average rate of production of a certain molecule per second per surface site (Andersen, Panosetti, and Reuter, *Frontiers in Chemistry*, 2019):

```
TOF_i = N_desorption_i / (N_sites * t_simulation)
```

where:
- `N_desorption_i` = total count of desorption events for product i during steady-state sampling
- `N_sites` = total number of surface sites in the lattice
- `t_simulation` = total simulated time during the steady-state sampling window

The TOF has units of s^{-1} site^{-1}, or equivalently molecules/(site * s). For surface-area-normalized rates, multiply by site density (e.g., ~1.5 x 10^{15} sites/cm^2 for close-packed metal surfaces).

**Steady-State Determination:**

Reaching steady state is critical for reliable TOF and selectivity measurements. The recommended approach uses an **iterative batch-means strategy** (Hoffmann and Matera, kmos documentation; Andersen et al., 2019):

1. Divide the simulation trajectory into successive batches of equal length
2. Compute the TOF within each batch
3. Check convergence criteria:
   - The variance of TOF estimates across batches should be sufficiently low
   - Adjacent batches should be statistically uncorrelated (verify via autocorrelation)
4. If not converged, double the trajectory length and repeat

A practical equilibration protocol:
- Run the simulation for an initial "warm-up" period (discard this data)
- Monitor surface coverages theta_i vs. time; steady state is indicated when coverages fluctuate around constant mean values
- Begin sampling only after coverages have stabilized
- Use EWMA (Exponentially Weighted Moving Average) for online detection of steady state

**Statistical Sampling Requirements:**

For reliable selectivity, the number of product desorption events for **each** product must be statistically significant:

- **Minimum events:** At least ~100-1000 desorption events per product species for ~10% relative error in TOF
- **For selectivity ratios:** The minority product must have at least ~100 events to establish its TOF to within ~10% relative precision
- **Error estimation:** Standard error of the mean TOF is estimated from batch means:
  ```
  sigma_TOF = std(TOF_batches) / sqrt(N_batches)
  ```
- **95% confidence interval:** TOF +/- 1.96 * sigma_TOF

For a system producing 95% product A and 5% product B, if 10,000 total product events are observed, only ~500 will be product B. The relative error in S_B is then ~1/sqrt(500) ~ 4.5%. For 1% selectivity products, one would need ~100,000 total events to get ~100 minority events.

**Multiple Replica Simulations:**

Running several independent KMC replicas (different random seeds) and averaging results is the gold standard for error estimation. The acceleration and sensitivity analysis work by Nielsen, Plaisance, and Vlachos (*J. Chem. Phys.*, 2013) and Dybeck et al. (*J. Chem. Theory Comput.*, 2017) recommend:
- At minimum 5-10 independent replicas
- Compute mean and standard deviation of TOF across replicas
- Parallel processing of replicas can decrease the sampling needed from a single trajectory

### 1.2 Selectivity Metrics

**Selectivity (mole-fraction based):**

```
S_i = TOF_i / sum_j(TOF_j)
```

where the sum runs over all product species j. This gives the fraction of all product molecules that are species i.

**Faradaic Efficiency from KMC:**

For electrocatalytic systems, Faradaic efficiency (FE) measures the fraction of total current going to each product:

```
FE_i = (n_i * TOF_i) / sum_j(n_j * TOF_j)
```

where:
- `n_i` = number of electrons transferred to make one molecule of product i
- `TOF_i` = turnover frequency for product i (from desorption event counting)

Equivalently, in terms of partial current densities:

```
FE_i = j_i / j_total
```

Common electron counts for CO2RR products:
| Product | Formula | n_i (electrons) |
|---------|---------|-----------------|
| CO | CO | 2 |
| Formate | HCOO^- | 2 |
| Methane | CH4 | 8 |
| Ethylene | C2H4 | 12 |
| Ethanol | C2H5OH | 12 |
| Methanol | CH3OH | 6 |
| H2 (HER) | H2 | 2 |

Common electron counts for NO3RR products:
| Product | Formula | n_i (electrons) |
|---------|---------|-----------------|
| NO2^- | NO2^- | 2 |
| NO | NO | 3 |
| N2O | N2O | 8 (per N2O) |
| N2 | N2 | 10 (per N2) |
| NH3 | NH3 | 8 |
| NH2OH | NH2OH | 6 |
| H2 (HER) | H2 | 2 |

**Computing FE from KMC Event Logs:**

Implementation in a KMC code:
1. Maintain a counter array `count[i]` for each product species i
2. At each product desorption event, increment the appropriate counter
3. Record the simulation time at the start (`t_start`) and end (`t_end`) of the steady-state window
4. Compute:
   ```python
   dt = t_end - t_start
   TOF_i = count[i] / (N_sites * dt)
   j_i = n_i * e * TOF_i * N_sites / A_geometric  # partial current density
   j_total = sum(j_i for all products)
   FE_i = j_i / j_total
   ```

### 1.3 KMC vs Mean-Field Microkinetic Models (MKM) for Selectivity

This is a critical topic with substantial literature demonstrating that KMC and MKM can give **dramatically different** selectivity predictions.

**When KMC and MKM Agree:**

In the absence of both (a) surface diffusion limitations and (b) adsorbate-adsorbate interactions, MKM and KMC approaches yield identical trends (Nunez, Filot, and Hensen, *Catal. Sci. Technol.*, 2021). The mean-field approximation is valid when:
- Adsorbates are well-mixed on the surface (fast diffusion)
- No lateral interactions modify rate constants
- No site-blocking creates spatial correlations

**When KMC Gives Different Selectivity from MKM:**

**Case 1: CO Oxidation on RuO2(110) - Correlation Effects Without Lateral Interactions**

Matera, Meskine, and Reuter (*J. Chem. Phys.*, 2011) demonstrated that even in a system **without** appreciable lateral interactions, mean-field rate equations fail to predict catalytic activity by **orders of magnitude**. The discrepancy arises from the inability of mean-field models to account for **vacancy pair formation** that is kinetically driven by ongoing reactions. The adlayer becomes spatially inhomogeneous purely due to the interplay between adsorption, desorption, and reaction rates -- an effect invisible to mean-field models.

**Case 2: Selective Oxidation of NH3 on RuO2(110) - Multi-Product Selectivity**

Hong, Rahman, Jacobi, and Bhattacharya (*J. Phys. Chem. C*, 2007; arxiv:1006.4297, 2010) used combined DFT+KMC to study NH3 oxidation on RuO2(110), which produces NO, N2, and N2O. Key findings:
- KMC correctly predicted **93% selectivity** for NO over N2, in very close agreement with experimental value of **95%**
- The selectivity is controlled by **spatial effects**: N+N recombination (to form N2) requires two adjacent N atoms, which is strongly inhibited by the presence of other surface intermediates that reduce N diffusion
- N+O recombination (to form NO) is far less affected because O from O2 dissociation is readily available nearby
- At ambient pressures, NO selectivity disappears due to abundant N species from active NHx decomposition -- a coverage-dependent spatial effect
- Mean-field models would dramatically overestimate N2 formation by ignoring these spatial constraints on bimolecular steps

**Case 3: CO2RR on Cu - Surface Diffusion and C-C Coupling**

Li and Maresi (*J. Chem. Phys.*, 2021) developed a lattice KMC approach for electrocatalytic CO2 reduction on Cu considering surface diffusion effects. Key findings:
- In the mean-field limit (infinitely fast diffusion), intermediates ^{12}CO* and ^{13}CO* would be well-mixed and no isotopic site selectivity could be observed
- KMC with finite diffusion rates revealed that the size of active sites and total adsorbate coverage strongly influence C2 product selectivity
- The spatial distribution of CO* intermediates on the surface determines the probability of C-C coupling events

**Case 4: Evaluating KMC vs MKM for Catalyst Design with Lateral Interactions**

Nunez, Filot, and Hensen (*Catal. Sci. Technol.*, 2021) systematically compared KMC and MKM for catalyst design studies:
- MKMs are ~**1000x faster** computationally than KMC
- Mean-field MKMs **cannot properly account for local coverage effects**
- When lateral interactions are included via cluster expansion in KMC, it can differentiate among highly active metals, but results are very sensitive to the set of included interaction parameters
- Mean-field implementation of lateral interactions causes **artificial overprediction** of activity of strongly binding metals
- Conclusion: KMC is essential when lateral interactions and spatial correlations are important for selectivity

**The Role of Spatial Correlations in Bimolecular Selectivity-Determining Steps:**

Bimolecular surface reactions (A* + B* -> C* or A* + B* -> AB(g)) are the steps most sensitive to spatial effects:
- The reaction rate depends on the probability of finding A* and B* on adjacent sites
- Mean-field approximation: rate = k * theta_A * theta_B * N_sites (assumes random mixing)
- KMC reality: the actual pair correlation function <n_A * n_B>_nn can differ significantly from theta_A * theta_B
- For selectivity-determining bimolecular steps (e.g., CO* + CO* -> OCCO* for C2 products vs. CO* + H* -> CHO* for C1 products), the spatial distribution of intermediates on the lattice directly controls the product ratio
- Island formation, surface segregation, and clustering of adsorbates all affect these pair correlations

---

## 2. Polarization Curves from KMC

### 2.1 Current Density Calculation

**Total Current Density:**

The total current density from an electrocatalytic KMC simulation is computed by summing the electron flux from all electrochemical (PCET) steps:

```
j_total = e * sum_i(n_i * TOF_i_desorption) * rho_sites
```

where:
- `e` = elementary charge (1.602 x 10^{-19} C)
- `n_i` = electrons transferred per molecule of product i
- `TOF_i_desorption` = desorption rate of product i (molecules / site / s)
- `rho_sites` = surface site density (sites/cm^2)

Alternatively, in SI units:
```
j [A/cm^2] = (e * rho_sites) * sum_i(n_i * TOF_i)
```

For a Cu(111) surface with rho_sites ~ 1.5 x 10^{15} cm^{-2}:
- If TOF_total_electrons = 10 s^{-1} site^{-1}, then j ~ 2.4 mA/cm^2

**Partial Current Density for Each Product:**

```
j_i = n_i * e * TOF_i * rho_sites
```

This allows direct comparison with experimental partial current densities measured by online product detection (GC, DEMS, etc.).

**Counting Electrons: At Desorption vs. At Each PCET Step**

There are two valid approaches, and they must give the same total current under steady state:

**Method 1: Count at product desorption (recommended for simplicity)**
- When a product molecule desorbs, credit all n_i electrons to that desorption event
- j_i = n_i * e * TOF_i_desorption * rho_sites
- Advantage: Simple, unambiguous product assignment
- Requirement: Must correctly assign n_i for each product

**Method 2: Count at each PCET step (more detailed)**
- Every time an electrochemical elementary step fires (e.g., CO2* + H+ + e- -> COOH*), count 1 electron
- j_total = e * sum_k(TOF_k_PCET) * rho_sites where k runs over all PCET steps
- Advantage: Can separately analyze current from each intermediate step
- Complication: Assigning partial current to a specific product requires knowing the eventual fate of each intermediate

**Under steady state, both methods give identical j_total**, because the flux through each PCET step must equal the sum of product desorption fluxes that pass through that step. Method 1 is generally preferred for computing FE and product-specific partial currents.

**Handling Intermediate PCET Steps:**

Some care is needed for branching mechanisms:
```
CO2 -> COOH* -> CO* (2e- total to CO)
                  \-> CHO* -> ... -> CH4 (8e- total to CH4)
```
If counting at desorption:
- Each CO desorption: 2 electrons
- Each CH4 desorption: 8 electrons
- The total is correct regardless of branching

If counting at each PCET step:
- The PCET step CO2 -> COOH* fires for **both** CO and CH4 pathways
- You cannot assign that electron to a specific product without tracking the trajectory of each intermediate

### 2.2 Potential Sweep in KMC

**Independent KMC Simulations at Each Potential:**

The standard approach for generating polarization curves from KMC is to run **independent** KMC simulations at each applied potential, as demonstrated in:
- Liu et al. (*ChemCatChem*, 2018): ORR on Ag(111), reproduced experimental polarization curves
- Voltage-dependent CO2RR KMC (*J. Phys. Chem. Lett.*, 2024): 178 elementary reactions on Cu(111)/(100)

Protocol:
1. Define a set of applied potentials: V = {V_1, V_2, ..., V_M}
2. At each potential V_m:
   a. Compute all potential-dependent rate constants k_j(V_m) using Butler-Volmer or Marcus theory
   b. Run KMC simulation to steady state
   c. Extract TOF_i for each product
   d. Compute j_i(V_m) and j_total(V_m)
3. Plot j vs. V to obtain the polarization curve

**Potential-Dependent Rate Constants:**

For electrochemical (PCET) elementary steps, the rate constant depends on the applied potential:

**Butler-Volmer formulation:**
```
k_forward(V) = k_0 * exp(-beta * e * (V - V_eq) / (k_B * T))
k_backward(V) = k_0 * exp((1 - beta) * e * (V - V_eq) / (k_B * T))
```
where beta is the symmetry factor (often assumed 0.5), V_eq is the equilibrium potential, and k_0 is the rate constant at zero overpotential.

**Marcus theory formulation (more accurate at large overpotentials):**
```
Delta_G_act(V) = (lambda + Delta_G_rxn(V))^2 / (4 * lambda)
```
where lambda is the reorganization energy and Delta_G_rxn(V) = Delta_G_0 + e * V. The Marcus-Hush-Chidsey (MHC) model further integrates over the electronic density of states of the electrode, preventing the unphysical inverted region prediction.

For non-electrochemical steps (thermal surface reactions, diffusion), rate constants are independent of potential.

**Equilibration Time Determination:**

How to know steady state is reached:
1. **Coverage monitoring:** Plot theta_i vs. time for all surface species. Steady state when all theta_i fluctuate around constant values.
2. **Running average of TOF:** Compute TOF in successive time windows. Steady state when the running average converges.
3. **Practical criterion:** Discard the first ~10-30% of the total simulated time as equilibration, then verify that TOF statistics in the remaining period are stationary.
4. **Automated approach:** Use the batch-means method with statistical tests for stationarity.

Typical equilibration times range from 10^3 to 10^8 KMC steps depending on the system stiffness. Systems with very fast reversible steps (e.g., fast diffusion) relative to slow reaction steps require longer simulated times.

**Statistical Error Estimation:**

For each potential point:
- Run at least 3-5 independent replicas
- Report mean j +/- standard error of the mean
- Alternatively, use batch-means within a single long trajectory
- Error bars on FE propagate from errors in individual j_i values

**Number of Potential Points:**

For a well-resolved polarization curve:
- **Minimum:** 10-15 points spanning the potential range of interest
- **For Tafel analysis:** At least 5-8 points in the linear log(j) vs. V region, with spacing of ~25-50 mV
- **Near onset potential:** Higher density of points (every 25 mV) to capture the transition
- **At high overpotential:** Wider spacing acceptable (50-100 mV)
- **Full range:** Typically 0 to -1.5 V vs. RHE for CO2RR, or 0 to +1.0 V vs. RHE for OER

### 2.3 Tafel Analysis from KMC

**Extracting Tafel Slope:**

The Tafel slope b is extracted from the linear region of the log(j) vs. V plot:

```
log10(j) = log10(j_0) + V / b
```

or equivalently:
```
b = dV / d(log10 j)    [units: mV/dec]
```

From KMC:
1. Compute j_total (or j_i for product-specific Tafel slopes) at each potential
2. Plot log10(j) vs. V
3. Identify the linear region (avoid mass-transport-limited or coverage-transition regions)
4. Fit a straight line: slope = 1/b (if plotting V on x-axis) or b (if plotting V on y-axis)

**Theoretical Tafel Slopes:**

For a mechanism with rate-determining step involving n_alpha electrons before the RDS and the RDS having symmetry factor beta:
```
b = 2.303 * R * T / ((n_alpha + beta) * F)
```

At 298 K:
| Mechanism | Tafel slope (mV/dec) |
|-----------|---------------------|
| 1st PCET is RDS (beta=0.5) | 120 |
| 2nd PCET is RDS | 40 |
| Chemical step after 1st PCET | 60 (if 1st PCET is quasi-equilibrated) |
| 1st PCET, coverage dependent | 60-120 (varies with theta) |

**Important caveat** from recent Bayesian analysis (Limaye et al., *Nature Communications*, 2021): Experimental and simulated Tafel slopes may not necessarily correspond to "cardinal" values (40, 60, 120 mV/dec) and can be coverage-dependent. KMC simulations naturally capture this coverage dependence.

**Exchange Current Density from KMC:**

The exchange current density j_0 is the current density at zero overpotential (V = V_eq):
```
j_0 = j_total(V = V_eq)
```

In practice, j_0 is obtained by extrapolating the linear Tafel region to V = V_eq. From KMC, one can directly simulate at V = V_eq and measure the (equal) forward and reverse current.

**Comparison with Experimental Tafel Slopes:**

The voltage-dependent KMC study of CO2RR on Cu (*J. Phys. Chem. Lett.*, 2024) demonstrated that KMC simulations with 178 elementary reactions could reproduce experimentally observed:
- Linear sweep voltammetry (LSV) profiles
- Potential-dependent product distributions
- Transition in rate-determining step with potential on Cu(111): from CO hydrogenation to CHO* at low overpotential to COH* formation at high overpotential

The ORR study on Ag(111) (Liu et al., *ChemCatChem*, 2018) showed that combined DFT+KMC simulations reproduced experimental polarization curves in both low and high potential regions.

---

## 3. Specific Multi-Product Electrocatalytic Systems Studied by KMC

### 3.1 CO2 Reduction (CO2RR) on Cu

**Study 1: Surface Diffusion Effects in CO2RR**

- **Paper:** "Effects of surface diffusion in electrocatalytic CO2 reduction on Cu revealed by kinetic Monte Carlo simulations"
- **Authors:** Li, Maresi
- **Journal:** *J. Chem. Phys.* 155, 164701 (2021)
- **DOI:** 10.1063/5.0065348
- **System:** CO2RR on Cu, lattice KMC with surface diffusion
- **Key features:**
  - Lattice KMC model considering CO* diffusion effects on C-C coupling
  - Motivated by isotopic labeling experiments (12CO2 + 13CO)
  - Demonstrated that mean-field limit (infinitely fast diffusion) fails to explain isotope scrambling data
  - Site size and adsorbate coverage strongly influence C2 selectivity
  - Delta-13C is sensitive to active site size but less to CO* diffusion rate within estimated ranges
- **KMC vs MKM:** Clear differences -- mean-field cannot capture the spatial effects that determine C2 selectivity

**Study 2: Voltage-Dependent CO2RR Mechanism**

- **Paper:** "Voltage-Dependent Electrochemical Carbon Dioxide Reduction Mechanism Unveiled by Kinetic Monte Carlo Simulation"
- **Journal:** *J. Phys. Chem. Lett.* (2024)
- **DOI:** 10.1021/acs.jpclett.4c03426
- **System:** CO2RR on Cu(111) and Cu(100)
- **Key features:**
  - 178 elementary reactions from DFT
  - Predicted LSV and potential-dependent product distribution
  - Cu(111): primarily CH4 (C1); RDS shifts from CO->CHO* to COH* formation with increasing overpotential
  - Cu(100): more active for C2H4 and C2H5OH; RDS is CO*+CO*->OCCO* (symmetric coupling)
  - Multiple products tracked simultaneously: CO, CH4, C2H4, C2H5OH, H2
  - Selectivity varies with potential, correctly reproducing experimental trends

### 3.2 Selective Oxidation of NH3 on RuO2(110)

**Study 3: NH3 Oxidation - NO vs N2 Selectivity**

- **Paper:** "Selective Oxidation of Ammonia on RuO2(110): a combined DFT and KMC study"
- **Authors:** Hong, Rahman, Jacobi, Bhattacharya (and related: Reuter group)
- **Journal:** *J. Phys. Chem. C* (2007); arXiv:1006.4297 (2010)
- **System:** NH3 + O2 on RuO2(110), products: NO, N2, N2O
- **Elementary steps:** ~20+ reactions including NH3 adsorption, sequential dehydrogenation (NH3->NH2->NH->N), N+O->NO, N+N->N2, NO+N->N2O, desorption
- **Key barriers:** NH3+O->NH+H2O (0.56 eV), N+N->N2 (0.27 eV), N+O->NO (0.14 eV)
- **Selectivity results:**
  - KMC predicted 93% NO selectivity, experiment: 95%
  - N2 formation inhibited by reduced N diffusion due to surface crowding
  - At ambient pressures: NO selectivity disappears due to enhanced NHx decomposition
- **KMC vs MKM:** Essential for this system -- spatial distribution of N* atoms determines N+N recombination probability. Mean-field would overestimate N2 production.

### 3.3 Oxygen Reduction Reaction (ORR) on Ag(111)

**Study 4: ORR Polarization Curve from DFT+KMC**

- **Paper:** "Oxygen Reduction Reaction on Ag(111) in Alkaline Solution: A Combined Density Functional Theory and Kinetic Monte Carlo Study"
- **Authors:** Liu et al.
- **Journal:** *ChemCatChem* (2018)
- **DOI:** 10.1002/cctc.201701539
- **System:** ORR on Ag(111) in alkaline solution, products: OH^- (4e-) and HO2^- (2e-)
- **Key features:**
  - Multiple pathways: 2e- (H2O2 pathway) and 4e- (H2O pathway)
  - Reproduced experimentally measured polarization curves in both low and high potential regions
  - Identified dominant pathway: *H2O-mediated 4e- associative pathway
  - *OH coverage inhibits O2 activation at high potentials
  - Selectivity between 2e- and 4e- products depends on potential
- **KMC advantage:** Captured coverage-dependent pathway switching that determines H2O2 vs H2O selectivity

### 3.4 Nitrogen Reduction Reaction (NRR) on Transition Metal Oxides

**Study 5: DFT-kMC-LSTM for NRR on V2O3**

- **Paper:** "Investigating High-Performance Non-Precious Transition Metal Oxide Catalysts for Nitrogen Reduction Reaction: A Multifaceted DFT-kMC-LSTM Approach"
- **Authors:** Lee, Pahari, et al.
- **Journal:** *ACS Catalysis* 13, 8327-8343 (2023)
- **DOI:** 10.1021/acscatal.3c01360
- **System:** N2 reduction on transition metal oxides (V2O3, Cr2O3, Fe2O3, etc.)
- **Key features:**
  - Integrated DFT + kMC + LSTM (machine learning for long-term degradation)
  - V2O3 predicted TOF 1000x higher than Ru catalyst
  - Complementary H transfer between TM and O sites accelerates NRR
  - kMC used to simulate the entire NRR pathway on candidate materials
  - Selectivity between NRR and competing HER explicitly tracked
  - Stable performance predicted for >10,000 hours

**Study 6: Bimetallic NRR Catalysts with HER Suppression**

- **Paper:** "DFT-kMC Analysis for Identifying Novel Bimetallic Electrocatalysts for Enhanced NRR Performance by Suppressing HER at Ambient Conditions Via Active-Site Separation"
- **Authors:** Lee, Chi Ho; Pahari, Silabrata; Sitapure, Niranjan; Barteau, Mark A.; Kwon, Joseph Sang-Il
- **Journal:** *ACS Catalysis* (2022)
- **DOI:** 10.1021/acscatal.2c04797
- **System:** Bimetallic alloys (RuTi, RuV2, Ru3W, RuZn3, RuZr) for NRR
- **Key features:**
  - Separate active sites for N2 and H adsorption designed to suppress HER
  - RuV2 showed superior NRR: TOF = 1.1 x 10^{-4} s^{-1} (1000x greater than Ru)
  - kMC explicitly modeled competition between NRR and HER pathways
  - Selectivity toward NH3 vs H2 was the primary metric
  - Active-site separation as a design principle requires spatial KMC

### 3.5 CO Oxidation on RuO2(110) and Pd(111)

**Study 7: CO Oxidation - Mean-Field Failure**

- **Paper:** "Adlayer inhomogeneity without lateral interactions: Rationalizing correlation effects in CO oxidation at RuO2(110) with first-principles kinetic Monte Carlo"
- **Authors:** Matera, Meskine, Reuter
- **Journal:** *J. Chem. Phys.* 134, 064713 (2011)
- **System:** CO + O2 on RuO2(110)
- **Key finding:** Even without lateral interactions, MKM fails by **orders of magnitude** due to kinetically-driven adlayer inhomogeneity (vacancy pair correlations)

### 3.6 Bifunctional Catalysis: Water-Gas Shift on Au/MoC

**Study 8: Synergic Effects on Bifunctional Catalysts**

- **Paper:** "Kinetic Monte Carlo Simulations Unveil Synergic Effects at Work on Bifunctional Catalysts"
- **Journal:** *ACS Catalysis* 9(10), 9117-9126 (2019)
- **DOI:** 10.1021/acscatal.9b02813
- **System:** Water-gas shift reaction on Au/MoC
- **Key features:**
  - Different catalyst regions have different functions (MoC: H2O dissociation; Au: COOH formation)
  - KMC essential to capture the spatial cooperativity between regions
  - Selectivity emerged from the spatial proximity of different functional sites
  - Mean-field models would miss the synergistic effect entirely

### 3.7 Single-Cluster Catalysis: CO2 + CH4 on Pt/HfC

**Study 9: First-Principles KMC for Single-Cluster Catalysis**

- **Paper:** "First-Principles Kinetic Monte Carlo Simulations for Single-Cluster Catalysis: Study of CO2 and CH4 Conversion on Pt/HfC"
- **Journal:** *ACS Catalysis* (2025)
- **DOI:** 10.1021/acscatal.4c07877
- **System:** Dry reforming of methane (CO2 + CH4) on Pt/HfC
- **Key features:**
  - Multi-product system with simultaneous CO2 and CH4 conversion
  - Investigated interplay between different catalytic sites
  - Evaluated activity, selectivity, and adlayer composition across operating conditions

### 3.8 NH3-SCR on Cu-Zeolites (Cu-CHA)

**Study 10: Cu-CHA Zeolite for NOx Reduction**

- **Paper:** "Kinetic Monte Carlo Simulations of Low-Temperature NH3-SCR over Cu-Exchanged Chabazite"
- **Journal:** *ChemPhysChem* (2024)
- **DOI:** 10.1002/cphc.202400558
- **System:** NH3-SCR of NOx on Cu-CHA zeolite
- **Key features:**
  - Mobile [Cu(NH3)2]+ complexes require pairing to form [Cu2(NH3)4O2]2+ peroxo-species for O2 activation
  - KMC essential: the Al-distribution in the zeolite framework controls Cu pairing probability
  - Mean-field models cannot capture the dependence on spatial distribution of Cu sites
  - Selectivity between different NOx reduction products depends on local Cu-Cu distances

---

## 4. Best Practices for Multi-Product KMC

### 4.1 Lattice Size Convergence for Selectivity

Lattice size is a critical parameter. Too small a lattice introduces artificial correlations from periodic boundary conditions.

**Convergence Testing Protocol:**
1. Start with a small lattice (e.g., 10 x 10 = 100 sites)
2. Double in each dimension: 20x20, 40x40, etc.
3. For each size, run to steady state and measure TOF_i for all products
4. Convergence is reached when doubling the lattice changes TOF by < some threshold (typically 5-10%)

**Literature values for converged lattice sizes:**
- CO oxidation on RuO2(110): fully converged at (40 x 20) = 1600 sites (Reuter and Scheffler, *Phys. Rev. B*, 2006)
- Typical catalysis models: 20x20 to 50x50 sites (400-2500 sites) are often sufficient for single-product reactions
- For multi-product systems with bimolecular selectivity-determining steps: larger lattices (50x50 to 100x100) may be needed to avoid finite-size effects on pair correlations
- For CO2RR C-C coupling: lattice must be large enough to contain multiple CO* "islands" to properly sample C-C coupling probabilities

**Practical recommendation:**
- Always perform a lattice size convergence study for selectivity, not just activity
- Selectivity can be more sensitive to lattice size than total TOF
- Report the lattice size used and demonstrate convergence in publications

### 4.2 Statistical Convergence Criteria

**For Total TOF:**
- Standard error of the mean < 5% of the mean TOF
- At least 1000 product desorption events (total across all products)

**For Selectivity:**
- Each product must have sufficient events
- For product with selectivity S_i, the relative error scales as:
  ```
  delta_S_i / S_i ~ 1 / sqrt(N_i)
  ```
  where N_i is the number of desorption events for product i
- For S_i = 0.01 (1% FE), need N_i ~ 10,000 for 1% relative error, which requires ~10^6 total events

**Batch-means convergence check:**
1. Divide the steady-state trajectory into K equal batches (K >= 20)
2. Compute TOF_i in each batch
3. Check that the batch-to-batch variance is consistent with the expected statistical noise
4. Verify zero autocorrelation between adjacent batches (if correlated, batches are too short)

### 4.3 Handling Very Rare Products (Low FE)

This is one of the most challenging aspects of multi-product KMC. For a product with FE < 1%, the number of desorption events in a typical simulation may be far too few for reliable statistics.

**Strategies:**

1. **Longer simulations:** Simply run for more KMC time. If computational cost is dominated by fast processes (diffusion, adsorption/desorption equilibria), this may be impractical.

2. **Rate constant rescaling (quasi-equilibrium approximation):**
   - Fast reversible processes (e.g., diffusion, adsorption/desorption) are rescaled to be slower while maintaining their equilibrium ratios
   - This accelerates simulated time per CPU time
   - Dybeck, Plaisance, and Vlachos (*J. Chem. Theory Comput.*, 2017) developed statistical criteria to ensure sufficient sampling under rescaling
   - Caution: rescaling must not alter the selectivity by changing relative rates of competing slow steps

3. **Parallel replica simulations:**
   - Run many independent replicas and aggregate statistics
   - Particularly effective for rare products: if 100 replicas each produce 10 minority events, the aggregate 1000 events gives ~3% relative error

4. **Transition path sampling or enhanced sampling:**
   - For extremely rare events, specialized techniques like forward flux sampling can be used
   - However, these are rarely needed for product selectivity (more for nucleation events)

5. **Sensitivity analysis to bound the selectivity:**
   - If a product is too rare to sample, use sensitivity analysis to determine which rate constants most influence its formation
   - This can bracket the expected selectivity even without direct sampling

### 4.4 Sensitivity Analysis: Degree of Rate Control in KMC

**Campbell's Degree of Rate Control (DRC):**

The DRC quantifies the sensitivity of the overall reaction rate to perturbations in individual elementary step rate constants (Campbell, *J. Catal.*, 1994, 2001; Campbell, *ACS Catal.*, 2017):

```
X_RC,i = (k_i / r) * (partial r / partial k_i)_{k_j!=i, K_eq,i}
```

where:
- X_RC,i = degree of rate control for step i
- k_i = rate constant of step i
- r = overall reaction rate (TOF)
- The derivative is taken holding all other rate constants and the equilibrium constant of step i fixed

Properties:
- Sum rule: sum_i(X_RC,i) = 1
- X_RC,i > 0: rate-limiting step (increasing k_i increases rate)
- X_RC,i < 0: rate-inhibiting step (increasing k_i decreases rate)
- X_RC,i ~ 0: kinetically unimportant step

**Degree of Selectivity Control (DSC):**

Campbell extended DRC to selectivity (ACS Catal., 2017):

```
X_SC,i^{A/B} = (k_i / S_{A/B}) * (partial S_{A/B} / partial k_i)
```

where S_{A/B} = TOF_A / TOF_B is the selectivity ratio.

This identifies which elementary steps most influence the selectivity between products A and B -- critical for rational catalyst design.

**Implementing DRC/DSC in KMC:**

The challenge: KMC is stochastic, so numerical derivatives are noisy.

**Three-stage approach** (Hoffmann, Bligaard, et al., 2017):

1. **Screening stage:** Perturb each rate constant by a large factor (e.g., 10x). Run short KMC simulations. Identify the ~5-10 steps with largest TOF change.

2. **Refinement stage:** For the important steps, use smaller perturbations (e.g., 1.1x and 0.9x) with longer simulations to get better estimates.

3. **Final stage:** Compute DRC as finite difference:
   ```
   X_RC,i ~ (k_i / TOF) * (Delta TOF / Delta k_i)
   ```
   Use centered finite differences with perturbation factor epsilon:
   ```
   X_RC,i ~ [ln(TOF(k_i * (1+epsilon))) - ln(TOF(k_i * (1-epsilon)))] / (2 * epsilon)
   ```

**Computational cost:** For N elementary steps, computing all DRC values requires ~2N additional KMC simulations (forward and backward perturbation for each step). With ~50 steps, this means ~100 KMC runs -- feasible but expensive.

**Alternative: Likelihood ratio method**
- Can extract sensitivity information from a **single** KMC trajectory
- But suffers from variance explosion for stiff systems (large time-scale separation)
- Works well when rate constants span < 3-4 orders of magnitude

### 4.5 Summary of Recommended Workflow for Multi-Product KMC

1. **Define the reaction network:** All elementary steps, rate constants, potential dependence
2. **Set up the lattice:** Choose lattice type (square, hexagonal, etc.), define site types
3. **Convergence tests:**
   a. Lattice size: Run at 20x20, 40x40, 80x80 and compare selectivities
   b. Time: Ensure steady state (batch-means test)
   c. Statistical: Ensure sufficient events for all products of interest
4. **Potential sweep:** Run independent simulations at 10-20 potential points
5. **At each potential:**
   a. Equilibrate (discard initial transient)
   b. Sample TOF_i for all products
   c. Compute j_i, j_total, FE_i
6. **Post-processing:**
   a. Plot polarization curves: j_total vs. V, j_i vs. V
   b. Plot FE_i vs. V
   c. Extract Tafel slopes
   d. Perform DRC/DSC analysis at key potentials
7. **Compare with experiment and MKM:**
   a. Overlay with experimental data
   b. Run equivalent MKM to identify where spatial effects matter

---

## 5. Key Software Frameworks for Multi-Product KMC

### 5.1 Zacros (Graph-Theoretical KMC)

- **Developer:** Stamatakis group (Oxford/UCL)
- **Language:** Fortran
- **Features:** Graph-theoretical cluster expansion for lateral interactions, arbitrary lattice geometries, general reaction mechanisms
- **Selectivity:** Native support for tracking multiple product TOFs
- **Companion tool:** ZacrosTools (Python) for automated analysis and visualization of selectivity, TOF, coverage
- **URL:** https://zacros.org

### 5.2 kmos

- **Developer:** Hoffmann, Matera, Reuter
- **Language:** Python frontend, Fortran backend
- **Features:** API-based model definition, automatic Fortran code generation, real-time visualization
- **Selectivity:** TOF for arbitrary product species, batch-means convergence
- **URL:** https://github.com/mhoffman/kmos

### 5.3 MKMCXX + Zacros (via AMS)

- **Developer:** SCM (Software for Chemistry & Materials)
- **Features:** Integrated MKM and KMC within the Amsterdam Modeling Suite, allowing direct comparison of mean-field and KMC results
- **URL:** https://www.scm.com/amsterdam-modeling-suite/kinetic-monte-carlo-and-microkinetics/

---

## 6. Critical Equations Summary

### Current Density and Faradaic Efficiency

```
j_i = n_i * e * TOF_i * rho_sites              [A/cm^2]
j_total = sum_i(j_i)                             [A/cm^2]
FE_i = j_i / j_total = (n_i * TOF_i) / sum_j(n_j * TOF_j)
```

### Selectivity

```
S_i = TOF_i / sum_j(TOF_j)                      [mole fraction]
```

### TOF from KMC

```
TOF_i = N_desorption_i / (N_sites * t_steady_state)   [s^{-1} site^{-1}]
```

### Tafel Slope

```
b = dV / d(log10 j) = 2.303 * k_B * T / (alpha * e)   [V/dec]
```
where alpha = n_alpha + beta (number of electrons before RDS + symmetry factor of RDS)

### Statistical Error

```
sigma_{TOF} = std(TOF_batches) / sqrt(N_batches)
sigma_{FE_i} ~ FE_i * sqrt((sigma_{TOF_i}/TOF_i)^2 + (sigma_{j_total}/j_total)^2)
```

### Degree of Rate Control

```
X_{RC,i} = (k_i / TOF) * (partial TOF / partial k_i)
```

### Butler-Volmer Rate Constant

```
k(V) = (k_B T / h) * exp(-Delta G_act^0 / (k_B T)) * exp(-beta * e * eta / (k_B T))
```
where eta = V - V_eq is the overpotential.

---

## 7. Key Gaps and Opportunities for SPARK

Based on this literature review, several opportunities exist for the SPARK framework:

1. **No existing general-purpose electrocatalytic KMC code** exists that natively handles potential-dependent rate constants, Faradaic efficiency tracking, and polarization curve generation. Zacros and kmos are thermal catalysis tools that can be adapted but lack built-in electrochemistry support.

2. **Multi-product tracking with automatic FE computation** would be a unique feature. The user defines n_i for each product, and the code automatically computes j_i(V) and FE_i(V).

3. **KMC for urea synthesis (CO2 + NO3RR co-reduction)** has **not yet been reported** in the literature. This represents a completely open field for SPARK.

4. **Integrated DRC/DSC analysis** within the KMC framework would be valuable -- most current implementations require manual scripting.

5. **Automatic steady-state detection** with batch-means convergence criteria is available in kmos but should be a standard feature.

6. **Rare product handling** via rate constant rescaling with selectivity-preserving constraints would be an advanced feature addressing a known pain point.

---

Sources:
- [Andersen, Panosetti, Reuter - A Practical Guide to Surface KMC Simulations (2019)](https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2019.00202/full)
- [Reuter, Scheffler - First-principles KMC for CO oxidation on RuO2(110) (2006)](https://link.aps.org/doi/10.1103/PhysRevB.73.045433)
- [Matera, Meskine, Reuter - Adlayer inhomogeneity without lateral interactions (2011)](https://pubs.aip.org/aip/jcp/article-abstract/134/6/064713/645472)
- [Li, Maresi - Surface diffusion in CO2RR on Cu by KMC (2021)](https://pubs.aip.org/aip/jcp/article/155/16/164701/199778)
- [Voltage-dependent CO2RR KMC on Cu (2024)](https://pubs.acs.org/doi/10.1021/acs.jpclett.4c03426)
- [Liu et al. - ORR on Ag(111) DFT+KMC (2018)](https://chemistry-europe.onlinelibrary.wiley.com/doi/abs/10.1002/cctc.201701539)
- [Hong et al. - NH3 oxidation on RuO2(110) DFT+KMC (2007-2010)](https://arxiv.org/abs/1006.4297)
- [Lee et al. - DFT-kMC-LSTM for NRR on TMOs (2023)](https://pubs.acs.org/doi/10.1021/acscatal.3c01360)
- [Lee et al. - DFT-kMC NRR bimetallic catalysts (2022)](https://pubs.acs.org/doi/10.1021/acscatal.2c04797)
- [Nunez, Filot, Hensen - KMC vs MKM with lateral interactions (2021)](https://www.sciencedirect.com/science/article/abs/pii/S092058612100119X)
- [KMC Simulations Unveil Synergic Effects on Bifunctional Catalysts (2019)](https://pubs.acs.org/doi/abs/10.1021/acscatal.9b02813)
- [Campbell - Degree of Rate Control review (2017)](https://pubs.acs.org/doi/full/10.1021/acscatal.7b00115)
- [Hoffmann et al. - Sensitivity analysis for KMC (2017)](https://arxiv.org/abs/1611.07554)
- [Dybeck et al. - Acceleration and sensitivity of lattice KMC (2017)](https://www.osti.gov/biblio/1512927)
- [Microkinetic modeling in electrocatalysis review (2021)](https://www.sciencedirect.com/science/article/abs/pii/S0021951721003523)
- [Gholizadeh et al. - Multiscale CO2RR modeling review (2025)](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202400898)
- [JACS Au - Kinetic Understanding of CO2RR Selectivity (2023)](https://pubs.acs.org/doi/10.1021/jacsau.3c00002)
- [KMC Review: Fundamentals and Challenges (2022)](https://pubs.aip.org/aip/jcp/article/156/12/120902/2840948)
- [Nellis et al. - Generalized Degree of Rate Control in Electrocatalysis (2021)](https://www.sciencedirect.com/science/article/abs/pii/S0021951721001184)
- [DFT-Based Multiscale Modeling of (Electro)Catalytic Reactions (2025)](https://pubs.acs.org/doi/10.1021/acscatal.5c07967)
- [Zacros KMC Software](https://zacros.org/)
- [kmos Framework](https://github.com/mhoffman/kmos)
- [SCM KMC and Microkinetics](https://www.scm.com/amsterdam-modeling-suite/kinetic-monte-carlo-and-microkinetics/)
- [KMC NH3-SCR on Cu-CHA (2024)](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cphc.202400558)
- [Limaye et al. - Bayesian Tafel analysis (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7846806/)
- [First-Principles KMC for Single-Cluster Catalysis on Pt/HfC (2025)](https://pubs.acs.org/doi/10.1021/acscatal.4c07877)
- [ZacrosTools Python Library](https://pubmed.ncbi.nlm.nih.gov/40658375/)
