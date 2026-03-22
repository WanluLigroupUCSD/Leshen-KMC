---
name: KMC project development context
description: User is developing Leshen-KMC and using Zacros 4.0 as reference software for feature comparison
type: project
---

User is working on KMC (Kinetic Monte Carlo) simulation software development.

- **Project repo**: cloned from https://github.com/WanluLigroupUCSD/Leshen-KMC to `/ibex/user/reny0b/zls/KMC/Leshen-KMC/`
- **User's software**: Leshen-KMC, has both Python (`mykmc/`) and Rust (`mykmc-rs/`) implementations
- **Reference software**: Zacros 4.0 (Fortran), extracted at `reference-software/zacros_4.0/`
- **Zacros build**: CMake configured and compilation started at `reference-software/zacros_4.0/build/` using gfortran 11.4.1
- **Goal**: Compare Zacros features with Leshen-KMC to guide development

**Why:** User wants to understand the feature gap between their KMC software and the established Zacros package to plan development priorities.
**How to apply:** When discussing features or planning tasks, frame in context of Zacros as the reference standard.
