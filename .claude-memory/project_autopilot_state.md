---
name: AutoPilot progress — Phase 1 complete, ready for Phase 2
description: T1-T4 done (HER Pt111 validated). Phase 1 complete. Next T5 CO2RR parameters. 5/11 tasks done.
type: project
---

**Status (2026-03-21): Phase 1 (HER Pt111) COMPLETE — 5/11 tasks done.**

Completed tasks:
- T1: DFT parameters collected (Li 2024, Skulason 2010, Karlberg 2007)
- T2: Python KMC + MKM model implemented (models/her_Pt111.py)
- T3: 6 validation tests PASSED (validate_her.py, ~1500 steps/s pure Python)
- T4: Polarization curves generated (run_polarization_her.py + results/)
- T11: Research plan documents written

**T4 key results:**
- j₀ = 0.45 mA/cm² (exp: 0.5-10 mA/cm²) — close to experimental
- Tafel slope = 81 mV/dec (mixed Volmer/Heyrovsky)
- Heyrovsky pathway dominates ~100%
- **MKM vs KMC divergence**: θ_H MKM=0.76 vs KMC=0.33 (~2x difference)
  - Lateral interactions (H*-H* +0.10 eV) strongly suppress coverage in spatial KMC
  - This validates the unique value of spatial KMC over mean-field
- Output: results/fig_summary.png (4-panel), results/her_mkm_polarization.dat, her_kmc_polarization.dat

**Next: T5 — CO2RR on Cu(100) DFT parameter collection (Phase 2 starts)**

**Why:** Phase 1 validates the software works correctly. Phase 2 is the actual research contribution.
**How to apply:** Phase 1 is reference/validation. Move to T5-T10 for CO2RR research.
