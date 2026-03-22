---
name: HER Pt(111) validation — Phase 1 COMPLETE with key findings
description: Phase 1 complete. HER validated with 6 tests + polarization curves. Key finding: MKM overestimates θ_H by 2x vs spatial KMC due to lateral interactions.
type: project
---

**Phase 1 — HER on Pt(111) validation: COMPLETE (2026-03-21)**

**Model:** 2 species (H*, empty), 5 processes (Volmer fwd/rev, Tafel, Heyrovsky, diffusion)
**DFT:** Li ACS Catal 2024, Skulason 2010, Karlberg 2007, Greeley 2004
**Files:**
- `models/her_Pt111.py` — KMC + MKM model definitions
- `validate_her.py` — 6 validation tests (all PASSED)
- `run_polarization_her.py` — T4 polarization curve script
- `results/` — data (.dat) + figures (.png)

**Polarization curve results (T4):**
- j₀ = 0.45 mA/cm² (exp ~0.5-10 mA/cm²)
- Tafel slope = 81 mV/dec
- Heyrovsky pathway ~100% dominant (Tafel negligible)
- **Key finding:** θ_H MKM=0.76 vs KMC=0.33 at U=-0.2V
  - Mean-field overestimates coverage by ~2x when H*-H* repulsion (+0.10 eV) is present
  - Spatial correlations from lateral interactions break mean-field assumption

**Performance:** Pure Python engine ~1500 steps/s. KMC 11-point polarization took ~13 min (20×20, 50k eq + 50k prod).

**Why:** Validates Leshen-KMC for electrochemistry. The MKM/KMC divergence demonstrates the unique value of spatial KMC.
**How to apply:** Use as reference for Phase 2. The validated model confirms Butler-Volmer, lateral, BEP implementations work correctly.
