# SPARK Multi-Lattice kMC — Design Note

Implementation plan for adding Hoffmann-Reuter-Scheffler 2015 multi-lattice kMC capability to `spark.engine`. Tracks the algorithm of:

> Hoffmann, M. J.; Scheffler, M.; Reuter, K. Multi-Lattice Kinetic Monte Carlo Simulations from First Principles: Reduction of the Pd(100) Surface Oxide by CO. *ACS Catal.* **2015**, *5*, 1199-1209. DOI 10.1021/cs501352t

Reference implementation surveyed: `kmcos` (https://github.com/kmcos/kmcos, GPL v3). Our SPARK reimplementation is independent and reads the kmcos source for **algorithmic ideas only** — no code lifted, MIT/Apache compatibility preserved.

---

## 1. What multi-lattice kMC is

Two crystal phases (e.g. Pd(100) metal lattice and √5×√5-PdO surface oxide) coexist in space. Atoms migrate between them as the system evolves. Conventional single-lattice kMC cannot represent this because the lattice geometry itself changes.

**Hoffmann's solution (paper §2.1):**

1. **Super-lattice** — work on a virtual lattice whose unit cell embeds the primitive cells of *all* participating sub-lattices simultaneously. Every site of every sub-lattice exists at every cell.
2. **`null` species** (paper concept) — sites belonging to a currently-disabled sub-lattice are flagged. No process can fire at a `null`-occupied site.
3. **Lattice-swap elementary processes** — local lattice transitions are encoded as ordinary kMC processes that *place or remove* `null` species at sites of the relevant sub-lattices. From the engine's point of view they are just multi-site processes (10-20 sites each, vs. 1-8 for surface chemistry).

**kmcos's actual realization** (`kmcos/types.py`, `kmcos/fortran_src/lattice.mpy`):

* `null` is **not** a hardcoded sentinel species — it can be any ordinary species the user chooses (in `sqrt5PdO/kmc_settings.py` the role is played by `'Pd'` for oxide-occupied metal sites). The "no process available" semantic is achieved simply by *not writing* any process whose conditions accept that species at that site.
* Each sub-lattice is a `Layer`. The `Project` carries a `layer_list` with `default_layer` and `substrate_layer` attributes.
* All layers share a single physical `cell`; each `Site` lives in exactly one `Layer` but every `(cell, layer, site_in_layer)` tuple maps to a unique global 1D site index via the `spuck` constant (Sites Per Unit Cell, summed across all layers):

  ```
  spuck = sum(len(layer.sites) for layer in project.layers)
  nr = spuck * (cx + Lx*cy + Lx*Ly*cz) + site_within_cell        # forward
  cz = (nr-1) // (Lx*Ly*spuck)                                    # inverse
  cy = (nr-1 - Lx*Ly*spuck*cz) // (Lx*spuck)
  cx = (nr-1 - spuck*(Lx*Ly*cz + Lx*cy)) // spuck
  s  = nr - spuck*(Lx*Ly*cz + Lx*cy + cx)
  ```

  Within `site_within_cell`, the first N₁ entries belong to layer 0, the next N₂ to layer 1, etc. — i.e. **layer is encoded inside `site_within_cell`, not stored in a separate dimension.** This is the algorithmic heart of "multi-lattice = single-lattice with extra sites" (paper §2.2).

---

## 2. SPARK current state vs target

`spark.engine` already has the foundations:

| Need | Current SPARK | Action |
|---|---|---|
| `Layer` data class | ✅ `spark/types.py:30 Layer` | Reuse |
| `Coord(layer=...)` layer-aware coordinate | ✅ `spark/types.py:43` | Wire into engine (currently ignored) |
| `Project.add_layer()` | ✅ | Extend to track `default_layer` / `substrate_layer` |
| Pairwise lateral interactions | ✅ `spark/engine.py:_setup_lateral_interactions` | No change — paper §2.3 lateral-interaction handling already matches |
| Site-type filtering on processes | ✅ `Process.site_type` | No change |
| Multi-layer 1D site indexing | ❌ Engine assumes single `nsites` array per cell, no `spuck` concept | **Build** |
| Cross-layer `Action` / `Condition` resolution | ❌ Engine ignores `Coord.layer` | **Build** |
| `default_layer` / `substrate_layer` registry | ❌ | **Add** |

---

## 3. Design decisions (locked in)

### 3.1 Indexing scheme — "shared nsites array, layer encoded in site_within_cell" (option 2a)

* Single 1D `species` array of length `Lx · Ly · Lz · spuck`.
* `spuck` is computed at engine setup as the total sites across all layers.
* Layer membership of a site index is derived by table lookup: `_site_to_layer[s_in_cell] -> layer_id`. No per-layer separate array — kmcos has run this scheme for 10+ years; performance and code simplicity both win.
* All existing single-layer models reduce to `nlayers = 1` with no behavior change.

### 3.2 No hardcoded `null` species

Follow kmcos: any user-chosen ordinary species can play the disabled-site role. We don't reserve a name. The `Process` system already supports this — just write processes that don't have `Condition` matching the disabled-site occupancy.

### 3.3 Layers share one `cell`

`Lattice.cell` is one 3×3 matrix. Each `Layer` carries only its sites' fractional positions. Two sub-lattices with different primitive cells must be embedded into a common commensurate super-cell *before* declaring them as layers (this is the user's job — paper §2.1 calls it "a suitable, large enough commensurate cell"). SPARK does not do super-cell construction automatically.

### 3.4 Cross-layer processes are first-class

A `Process` whose `conditions` / `actions` reference coords with different `layer` fields is treated identically to single-layer processes by the engine — they just touch a wider site-stencil. Paper §2.2: "From a purely algorithmic point of view ... multi-lattice kMC is identical to single-lattice kMC."

### 3.5 GPL boundary

Read kmcos source for ideas / verification only. **No code copied.** All new SPARK code is original. License remains permissive.

---

## 4. Phase plan

### Phase A — Data layer (`spark/types.py`), ~2 days
- A.1 `Lattice` gains `default_layer`, `substrate_layer` attrs (init when first `add_layer` call)
- A.1 `Project` maintains `_layer_site_offsets[layer_id] -> first_site_in_cell_id` so that `(layer, site_in_layer) -> site_in_cell_id` is O(1)
- A.2 Helper functions `lattice_to_nr(cell, site_in_cell, system_size, spuck)` and `nr_to_lattice(nr, system_size, spuck)` in `spark/engine.py` (or new `spark/indexing.py`)
- A.2 Unit test `tests/test_multi_lattice_indexing.py` — round-trip nr ↔ (cell, site) for both 1-layer and 3-layer toy projects; assert single-layer test from existing engine tests still passes byte-for-byte

**DoD:** All existing tests green. New indexing test green. No engine changes yet.

### Phase B — Engine (`spark/engine.py`), ~3 days
- B.1 `_setup_lattice` rewritten to use `spuck` and the global 1D `species` array. Introspects `len(project.layers)`; falls through to existing single-layer path when `nlayers == 1` (this is automatic since `spuck = nsites_layer0` then).
- B.1 `set_site_types_region` updated to accept layer-aware predicates.
- B.2 `_setup_processes` resolves each `Condition.coord` and `Action.coord` to a site-in-cell offset using the layer registry. Cross-layer offsets handled by simply pointing into the appropriate slot of `site_within_cell`.
- B.2 New unit test `tests/test_cross_layer_process.py` — declare 2 layers, write a single cross-layer process (e.g. "remove species X from layer 0, place species Y on layer 1"), step the engine, assert both species changes occurred at the correct global site indices.

**DoD:** All existing tests green. Cross-layer process test green. `git diff --stat spark/engine.py` shows additions only at `_setup_lattice` and `_setup_processes`; no algorithmic changes elsewhere.

### Phase C — Validation example (`examples/`), ~2 days
- C.1 `examples/multi_lattice_PdO_reduction.py` — implements a 2-layer Pd(100)+√5-PdO toy model in pure SPARK API
  - Layer 0 = `Pd100`: 5 hollow + 10 bridge sites per cell
  - Layer 1 = `PdO`: 4 hollow + 2 bridge + 4 metal-Pd sites per cell
  - Species: `empty`, `CO`, `O`, `Pd_atom` (this last plays the "disabled" role on Pd100 sites when oxide is intact)
  - Processes: `CO_ads` / `CO_des` / `O2_ads` / `O2_des` on each layer + diffusion + ~5 `destruct_*` processes that flip a region from PdO → Pd100 occupancy + 1 `oxidize` reverse
  - DFT-style placeholder rates (no need to be physically tight; structural validation only)
- C.1 Run 1000 steps; collect trajectory; visually verify a destruct cluster propagates through the lattice (basic sanity, not paper-level reproduction)

**DoD:** Example runs to completion with `python examples/multi_lattice_PdO_reduction.py`. Lattice-swap process fires at least once. Trajectory dump shows species correctly flipped on both layers.

### Phase D — Documentation + commit, ~1 day
- D `spark/__init__.py` top-level docstring grows a `Multi-Lattice kMC` section (alongside Lattice / Off-Lattice / Dynamic)
- D `README.md` — short paragraph + reference to this design note
- D `EXECUTION_STATE.json` `software_status.layers` — new entry `multi_lattice_kmc`
- D Commit + push `feat(engine): multi-lattice kMC support per Hoffmann-Reuter 2015`

**DoD:** GitHub `multi-function` HEAD updated, three deployments (local + GitHub + Shaheen) re-synced.

**Total:** 8 working days, with checkpoints at end of each phase.

---

## 5. What we are NOT doing in this round

- Not implementing *automatic* super-cell construction — user supplies the commensurate cell.
- Not implementing the full DFT data pipeline of the 2015 paper (their Pd(100)/√5-PdO model has 40+ processes with specific DFT barriers; we ship a structural toy).
- Not touching `spark.offlattice` or `spark.dynamic`. Multi-lattice is a `spark.engine`-only addition.
- Not bringing in any kmcos-licensed Fortran backend code. Pure Python implementation that reuses existing engine vectorization.

---

## 6. Risk register

| Risk | Mitigation |
|---|---|
| Single-layer regression after engine refactor | Phase A.2 + B.1 keep single-layer path identical via `nlayers == 1` shortcut; full existing test suite must pass at every commit |
| Performance regression (Python overhead per cross-layer event) | Use `np.ndarray` for `_layer_site_offsets`; per-step lookup is `O(1)` array index, no dict |
| Subtle off-by-one in `nr ↔ (cell, layer, site)` map | Phase A.2 indexing test does forward+inverse round-trip on every site of a 3×3×1 × 3-layer toy; CI will catch any drift |
| Lateral interaction across layers semantically ambiguous | Defer to Phase B.2 review; if ambiguous, scope cross-layer LIs to 0 in v1 (pairwise within-layer only — matches paper §2.3) |

---

## 7. Glossary

* **Layer / sub-lattice**: one of N coexisting commensurate lattices in the super-cell.
* **`spuck`**: Sites Per Unit Cell — total sites across all layers in a single physical cell.
* **`site_within_cell`** ∈ `[0, spuck)`: flat index that encodes (layer_id, site_in_layer).
* **`nr`** ∈ `[0, Lx·Ly·Lz·spuck)`: global 1D site index used by the engine.
* **Lattice-swap process**: an elementary kMC process whose conditions/actions cross layer boundaries.
* **Disabled site**: a site currently in a non-active sub-lattice; occupied by a chosen species (no `null` reserved name).
