---
name: τ-leaping algorithm plan
description: User wants to attempt τ-leaping implementation later for stiff electrochemical KMC systems (10²-10⁵× speedup potential)
type: project
---

**τ-leaping for stiff KMC systems** — deferred, user will attempt later.

Algorithm outline:
1. Detect fast/slow process separation (rate ratio > 10³)
2. Fast processes: execute Poisson(k_fast × τ) events per leap step
3. Slow processes: keep exact BKL
4. Adaptive τ selection to maintain accuracy (Cao-Gillespie-Petzold criterion)

**Why:** Electrochemical KMC is inherently stiff — diffusion ~10⁸ s⁻¹ vs PCET ~10⁻² s⁻¹. 99.999999% of BKL steps are wasted on fast diffusion. τ-leaping skips these.

**How to apply:** Implement as a new method `do_kmc_step_leaping()` in engine.rs alongside existing `do_kmc_step()`. No existing KMC software (Zacros, kmos, SPPARKS, MonteCoffee) has this — would be a unique competitive advantage.

**Key references:**
- Gillespie, J. Chem. Phys. 115, 1716 (2001) — original τ-leaping
- Cao, Gillespie & Petzold, J. Chem. Phys. 124, 044109 (2006) — adaptive τ selection
- Chatterjee & Vlachos, J. Chem. Phys. 122, 024112 (2005) — spatial τ-leaping for lattice KMC
