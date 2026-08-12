# Test specification — spark-atomistic (Python reference)

Normative inputs, and the only design inputs consulted:

| Document | Self-excluding SHA-256 | Verified |
|---|---|---|
| `OFFLATTICE_OTF_KMC_SPEC_V1.md` | `8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84` | recomputed, matches |
| `..._ERRATA_1.md` | `52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40` | recomputed, matches |
| `..._ERRATA_2_PARITY.md` | `eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995` | recomputed, matches |
| `..._ERRATA_3.md` | `eba384af3694c5f3997caf28829e56d188ef9929f29d99a2116520f0067d8a96` | recomputed, matches |

## How to run

```text
python3 -m unittest discover -s tests -p 'test_*.py'
```

`tests/fixtures.py` puts `src/` on `sys.path`, so no install step and no `PYTHONPATH`
is required; every test module imports it first. Each module also runs standalone.

## Run record

| Date | Interpreter | Result |
|---|---|---|
| 2026-08-11 | CPython 3.12.3 (GCC 13.3.0), Linux x86-64 | **68 tests, 68 passed, 0 failed, 0 skipped** |
| 2026-08-11 (cross-language round) | CPython 3.12.3, Linux x86-64 | **77 tests, 77 passed, 1 skipped without `SPARK_XLANG_OUT`; 77 passed, 0 skipped with it** |
| 2026-08-11 (defect-repair round, F7/F8/F9) | CPython 3.12.3, Linux x86-64 | **100 tests: 99 passed, 0 failed, 1 skipped without `SPARK_XLANG_OUT`** |
| 2026-08-12 (normative Errata 3) | CPython 3.12.3, Linux x86-64 | **85/85 cross-language cases, 46/46 fixtures; core 83, adapter 2** |

Per module: `test_spec` 24, `test_checkpoint` 22, `test_engine_output` 12,
`test_saddle_domain` 10, `test_cross_language_parity` 9,
`test_rate_detailed_balance` 8, `test_event_application` 8, `test_saddle_order_gate` 7.

**A green suite is not conformance.** These are behavioural tests.
The normative Errata 3 tiered comparison passes all 46 fixtures. The implementation
status stays `implemented_unvalidated` because parity does not establish physical accuracy.

## Rules every test in this suite obeys

1. **Every discrepancy is decided by quoting a requirement ID.** No test is adjusted to
   make it pass.
2. **A fixture asserted to be REJECTED is paired with a baseline that is ACCEPTED.**
   A rejected fixture also satisfies a rejection assertion when it is rejected for an
   unrelated reason, which is how `corpus/minimal_model.json` turned the retry/alpha
   test into a false positive before it was deleted.
3. **Golden files are referenced, never copied.** `E2-PAR-001` requires both backends to
   consume the same canonical corpus; the single copy lives at
   `../spark-atomistic-rs/tests/corpus/` and is reached through `fixtures.CORPUS` and
   `xlang_harness.XLANG`.
4. **Defects are quantified**, not described.
5. **A divergence is reported, never normalised.** The comparator prints the exact first
   differing byte offset and gates every core artifact; adapter artifacts are reported
   separately under `E3-PAR-001`.

## Coverage against `E2-PAR-002`

`E2-PAR-002` enumerates eleven mandatory categories. `NOT COVERED` below means no
executed Python test exists yet; it does not mean the behaviour is absent.

| # | `E2-PAR-002` category | Status | Test |
|---|---|---|---|
| 1 | capability request and exact capability value | covered | `test_capability_response_matches_shared_parity_golden_byte_for_byte`, `test_capabilities_are_conservative_and_immutable` |
| 2 | valid minimal model; metadata-only variants share a config digest | covered | `test_schema_fixture_and_explicit_run_gate`, `test_metadata_is_free_form_and_excluded_from_the_config_digest` |
| 3 | duplicate key, malformed UTF-8, lone surrogate, nonfinite, integer boundary, out-of-domain integer, unknown key, adapter extension | covered | `test_duplicate_key_golden_file_is_refused_at_any_depth`, `test_nonfinite_tokens_and_lone_surrogates_are_refused`, `test_portable_integer_boundary_and_out_of_domain`, `test_unknown_keys_are_refused_in_every_normative_object`, `test_extension_is_excluded_from_the_config_digest` |
| 4 | canonical numbers: negative zero, `1e-6` and neighbours, `1e21` and neighbours | covered (7/7) | `test_canonical_number_cases_match_the_shared_golden` |
| 5 | translated, periodic-image, same-species-permuted identities | partial | `test_state_identity_translation_image_and_same_species_permutation`, `test_candidate_hash_never_blocks_tolerance_geometry_match` — a dedicated periodic-image wrap case is **NOT COVERED** |
| 6 | Errata-1 boundaries, Philox zero key/counter, trajectory and saddle derivation, state restore, counter carry | covered | `PhiloxCorpusTests` (5 tests) |
| 7 | reversible pair, reciprocal mapping, tolerated and rejected negative barrier, detailed-balance failure, parallel saddles | partial | reversible pair and reciprocal mapping are exercised through the checkpoint fixture and `test_reciprocal_event_pair_is_recomputed_not_trusted`; the four negative-barrier / detailed-balance / parallel-saddle cases are **NOT COVERED** |
| 8 | deterministic discovery class sequence, out-of-order completion commit | **NOT COVERED** | — |
| 9 | two-event KMC selection, application rollback, lost-rate bound, historical snapshot replay after catalog growth | partial | `HistoricalRateSnapshotReplayTests` (7 tests) covers two-event selection and the snapshot-after-growth case; rollback and the lost-rate bound are **NOT COVERED** |
| 10 | clean checkpoint, corrupted hash, incompatible digest, recursive nested corruption, cancellation, resource limit, exact next-event resume | partial | `DeepCheckpointInvariantTests` (14 tests) covers the first four; cancellation, resource-limit checkpointing, and exact next-event resume are **NOT COVERED** |
| 11 | basin-disabled capability, input rejection, checkpoint record, serial fallback | partial | `test_basin_enabled_is_schema_valid_but_never_advertised`, `test_basin_record_is_frozen_disabled`; serial fallback after a `BASIN_DISABLED` outcome is **NOT COVERED** |

## Module map

### `test_spec.py` — wire contract (24 tests)

`ContractSmokeTests` (12) — capability value shape and byte-identity against the shared
golden, machine-readable spec SHAs, request envelope and the `E2-API-004` trust gate,
duplicate-key rejection layer, exact severity strings, retry/alpha rejection against an
accepted baseline, Errata-1 uniform boundaries, triclinic closest image, state identity
under translation and same-species permutation, rate-table snapshot envelope shape.

`ParityCorpusTests` (7) — `E2-CAN-004` canonical numbers 7/7, `E2-JSON-002` integer
domain, `E2-JSON-001`/`E2-JSON-003` malformed input, nested duplicate keys,
`E2-CAN-007` metadata digest invariance paired with a behavioural field that must change
the digest, `E2-SCOPE-003` unknown keys in five normative objects, `E2-SCHEMA-011`
basin `true` validating while the capability stays `disabled`.

`PhiloxCorpusTests` (5) — `E2-RNG-001..007`: zero-key/zero-counter block against the
shared golden, both boundary records and their normative SHA-256 line hashes,
128-bit counter carry, the buffered-block/`next_pair` state machine and exact restore
round trip, and scheduling-independent substream derivation including the length-prefix
separation of state ID from search class.

### `test_checkpoint.py` — sections 10 and 11 (22 tests)

`HistoricalRateSnapshotReplayTests` (7). The fixture commits two events from one origin,
takes one KMC step, then commits a third event from the same origin, then checkpoints.
Measured discriminators:

| Quantity | Historical snapshot (2 events) | Grown catalog (3 events) |
|---|---|---|
| total rate | `17.762303659592` s⁻¹ | `138.121865606364` s⁻¹ |
| Δt | `0.018752335929856312` s | `0.0024115275568463812` s |
| selected event | differs | differs |

The Δt disagreement is 87.1 % relative, about 1.7 × 10¹⁴ times the `PAR-005` replay
tolerance of relative 5 × 10⁻¹⁵.

- ACCEPTED: `test_replay_accepts_the_historical_snapshot_after_catalog_growth`.
  **Mutation-verified**: replacing `_restore_rate_snapshot` with a version that rebuilds
  the table from the restored catalog makes this test fail with `CHECKPOINT_CORRUPT` /
  `E2-KMC-005`. The test therefore genuinely kills the defect it is written against.
- REJECTED counterparts: the current-catalog time increment, the current-catalog
  selection, and the current-catalog snapshot substituted into the record — each must
  raise `CHECKPOINT_CORRUPT` / `E2-KMC-005`.
- `E2-KMC-002` snapshot envelope shape, hash, and rate/log-rate agreement.

`DeepCheckpointInvariantTests` (14) — one accepted baseline re-verified before every
rejection: `E2-CKPT-001` envelope hash and canonical bytes, `E2-CKPT-002` exact
fifteen-key payload, `E2-CKPT-003` digest mismatch as `CHECKPOINT_INCOMPATIBLE` rather
than `CHECKPOINT_CORRUPT`, `E2-BASIN-001` frozen disabled record, `E2-CKPT-004` flag key
set and derived consistency, `E2-CKPT-005` resource counters with empty retry history
and `saddle_attempts_by_state` equal to the discovery attempt counts, `E2-CKPT-006`
catalog digest over the object without its own digest and multiplicity symmetry,
`E2-EVENT-006` recomputed pair and event identifiers, `E2-RATE-001` recomputed raw
barriers, `E2-CKPT-008` step/log/trajectory length identity and `2 × step_index`
consumed uniforms, `E2-RNG-005` complete and exact substream map, `E2-DISC-004` counter
identities, `E2-CKPT-009` canonical bytes with no trailing newline and no temporary
left behind.

`CheckpointRestoreOrderTests` (1) — records a **known deviation** from the
`E2-CKPT-007` verification order; see "Findings" below.

### `test_engine_output.py` — output transaction (12 tests)

`OutputReservationTests` (6) — reservation strictly precedes `_atomic_write`, reserved
bytes equal written bytes exactly (the self-referential `output_bytes` field is driven
to a fixed point first), the stored counter equals the on-disk file size, an output-byte
limit refuses the write with `RESOURCE_LIMIT` / `RES-002` leaving the output directory
empty and `checkpoint_sequence` rolled back, and failures injected at both the content
flush and the atomic replacement leave the previous artifact byte-identical with no
`*.tmp-*` sibling.

`OutputCollisionPolicyTests` (4) — `IO-004` collision policy, `IO-003` path validation.

`AdapterBoundaryTests` (2) — `E2-SCOPE-004`: the run-request `extension` never enters
the config digest, and a relative `calculator_command` is refused.

### `test_saddle_domain.py` — solver domain (10 tests)

Characterisation of `DirectionalDimerSearcher` on the analytic double well
`V = Σᵢ (xᵢ²−1)² + yᵢ² + zᵢ²`. See "Findings".

### `test_rate_detailed_balance.py` — `E2-EVENT-003` residual (8 tests)

`CommitTimeResidualTests` (4) — the residual is an algebraic identity of
`E2-RATE-001`/`E2-RATE-002`, so no commit-time gate may be applied to it; includes the
anti-vacuity control that the sweep reaches and passes through the rate block.
`RestoreTimeDetailedBalanceTests` (4) — the `E2-CKPT-007`(5) recomputation, where the
stored residual IS independent, with an accepted baseline, resealed tampers, and a
control proving the tamper is not answered by the catalog digest. See F7.

### `test_saddle_order_gate.py` — `SADDLE-005` orthogonal evidence (7 tests)

Exact quadratic fixture with a prescribed Hessian spectrum in a Householder-rotated
eigenbasis: index-1 accepted, index-2 and index-7 rejected, local maximum rejected, the
reported value bounded below by the true lowest orthogonal eigenvalue and above by the
random start it came from, and the rotation budget bounded. See F8.

### `test_event_application.py` — `EVENT-004` application (8 tests)

Genuine A→B double-well event accepted; single-field fabrications of
`destination_state_id`, `unstable_direction` and `saddle_positions` rejected; an event id
in no catalog rejected for all nine (current, destination) combinations; and the captured
minimiser request pinned to `saddle ± endpoint_displacement · unstable_direction`. See F9.

## Findings

### F1 — `IO-004` output-collision hole. FIXED.

`ReferenceEngine._validate_output_preflight` applied the resume exemption to all three
output paths at once, so an existing `summary_path` or `trajectory_path` was tolerated
whenever the checkpoint existed and `resume=true`. `IO-004`: "Existing outputs MUST NOT
be overwritten unless `overwrite=true`. Existing compatible checkpoint plus
`resume=true` resumes; every other collision returns `OUTPUT_EXISTS`." The exemption
covers the checkpoint only. Measured before the fix: `["checkpoint.json",
"summary.json"]` with `resume=true, overwrite=false` returned `True` and proceeded;
after the fix it raises `OUTPUT_EXISTS`. Reverting the fix fails three subtests of
`test_resume_tolerates_only_the_checkpoint_collision`. The emitted
`context.requirement_id` is deliberately left at `E2-SCHEMA-010` so the observable
response string does not move relative to the other backend.

### F2 — `E2-CKPT-007` verification order. NOT CHANGED; escalated.

`E2-CKPT-007` orders step 8 as "Philox key/counter/buffer/count relations and complete
substream map" and step 9 as "trajectory sequence/state chain using each historical rate
snapshot". `validate_checkpoint_payload` runs the trajectory replay first. For a
checkpoint damaged in both places the emitted `context.requirement_id` is `E2-KMC-005`
where the spec order implies `E2-RNG-005`. The terminal `status` (`CHECKPOINT_CORRUPT`)
and `exit_code` (74) are unaffected. Not changed here because `context.requirement_id`
is part of the public response that `E2-PAR-003` requires to be byte-identical across
backends, so a unilateral change could create the mismatch it aims to remove. Pinned by
`CheckpointRestoreOrderTests.test_dual_defect_reports_the_trajectory_requirement`.

### F3 — saddle search on the double well: parameter-domain limitation, not a defect.

Analytic answer: minima at x = ±1 (E = 0), first-order saddle at x = 0 (E = 1.0 eV),
forward and reverse barriers exactly 1.000 eV. Hessian eigenvalues at (x, 0, 0) are
`(12x²−4, 2, 2)`.

**Inside its domain the solver is correct.** From x₀ = 0.30 with
`force_tolerance = 1e-4`: saddle at x = 0, E = 1.0 eV, curvature −4, orthogonal
curvature +2, both endpoints at |x| = 1, forward and reverse barriers 1.000 eV to
8 decimal places.

**Accuracy is bounded by `force_tolerance`, not by the algorithm.** Near the saddle
V ≈ 1 − 2x² and |F| ≈ 4|x|, so the force-tolerance ball has radius f/4 and admits an
energy error of at most f²/8:

| `force_tolerance` (eV Å⁻¹) | measured \|E_saddle − 1.0\| (eV) | bound f²/8 (eV) |
|---|---|---|
| 5 × 10⁻² | 2.90 × 10⁻⁴ | 3.13 × 10⁻⁴ |
| 1 × 10⁻² | 1.17 × 10⁻⁵ | 1.25 × 10⁻⁵ |
| 1 × 10⁻³ | 5.99 × 10⁻⁸ | 1.25 × 10⁻⁷ |
| 1 × 10⁻⁴ | 1.73 × 10⁻⁹ | 1.25 × 10⁻⁹ |
| 1 × 10⁻⁵ | 2.43 × 10⁻¹¹ | 1.25 × 10⁻¹¹ |
| 1 × 10⁻⁶ | 2.43 × 10⁻¹³ | 1.25 × 10⁻¹³ |

**Why (0.9, 0, 0) burns 6000 evaluations.** At x = 0.9 the Hessian is
diag(5.72, 2, 2): every eigenvalue is positive and the *lowest* mode is transverse, not
the reaction coordinate. Minimum-mode following therefore locks onto y/z. The reflected
force climbs that transverse mode while the dimer centre slides into the minimum at
x = +1 (centre x reaches 0.99803616 by iteration 50 and 0.99999999 by iteration 200).
Measured transverse growth is 1.020000000 per iteration against the predicted
1 + Δt·λ⊥ = 1 + 0.01 × 2 = 1.02, agreeing to 9 decimal places; the transverse
coordinate reaches 5.2 × 10¹³ Å by iteration 1800 (E ≈ 2.7 × 10²⁷ eV). With
`max_iterations = 2000` this is exactly `SADDLE_NOT_FOUND` after 3 × 2000 = 6000
evaluations.

Two thresholds explain the behaviour: the reaction coordinate has negative curvature
only for |x| < 1/√3 = 0.5774, and stops being the lowest mode above |x| = 1/√2 = 0.7071.
Neither is a threshold in the *outcome*. Measured success rate over eight independent
search substreams per launch point:

| x₀ | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.55 | 0.6 | 0.65 | 0.70 | 0.75 | 0.8 | 0.9 | 1.0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| OK / 8 | 8 | 8 | 8 | 8 | 8 | 7 | 8 | 6 | 6 | 4 | 1 | 0 | 0 | 0 |

The initial dimer mode is drawn from a substream keyed by `state_id`, which hashes the
binary64 energy, so two adjacent doubles for the same launch coordinate draw unrelated
modes: x₀ = `0.58999999999999996891` fails while x₀ = `0.59000000000000007994`
succeeds. A single scan of x₀ therefore *looks* like a sharp ragged threshold and is
not one; only the success rate is reproducible.

**No fix applied.** `SADDLE-002` permits any first-order-saddle solver, `SADDLE-005`
defines the acceptance evidence that the solver does emit, `SADDLE-006`/`SADDLE-007`
require only that a failed attempt report its termination reason and not count as a
redundant discovery, and `E2-STATUS-002` classifies `SADDLE_NOT_FOUND` as a
`candidate reject` that cannot terminate a public operation. Nothing in the
specification obliges convergence from a launch point with no negative eigenvalue, so
inventing an algorithm here would be unmandated behaviour.

**One real consequence worth escalating.** Beyond about 2010 dimer iterations the
transverse runaway overflows binary64: at `max_iterations = 20000` the run stops at
calculator evaluation 6031 with `CALCULATOR_FAILURE` (`CALC-006`). That status is *not*
in the engine's candidate-reject continue set
(`SADDLE_NOT_FOUND, INVALID_SADDLE, SADDLE_WRONG_BASIN, ENDPOINT_COLLAPSED,
RELAX_NOT_CONVERGED`), so a larger iteration budget converts one rejected candidate into
a terminated run with exit 69.

### F4 — the shared `e2_minimal_model.json` cannot execute. Reported, NOT changed.

Two measured, independent problems with the shared golden model. It is in the Rust tree
and shared by both backends, so it is deliberately left untouched; the resolution is a
cross-language decision.

1. `resources.resident_memory_bytes` is `1000000` (1 MB). `ResourceLedger.check_wall_time`
   raises `RESOURCE_LIMIT` (`RES-002`) when the resident set reaches the limit, and it is
   consulted on every `reserve_evaluation`. Measured `ru_maxrss` of a bare CPython 3.12.3
   process running this package is 20.6–26.1 MB, i.e. 21–26× the limit, so the first
   calculator evaluation of any run using this model fails before doing any work.
2. `relaxation.force_tolerance` is `0.05` eV Å⁻¹ while `kinetics.barrier_tolerance` is
   `1e-10` eV. On the corpus double well the force tolerance admits an energy error of
   f²/8 = 3.125 × 10⁻⁴ eV, which is 3.1 × 10⁶ times the barrier tolerance. A measured
   consequence: from x₀ = 0.50 the "+" endpoint is force-converged while still on the
   barrier top and the raw forward barrier comes out **negative**, which `E2-RATE-001`
   turns into `RATE_INVALID` — fatal in strict mode. Pinned by
   `test_loose_tolerance_can_produce_a_negative_forward_barrier`.
   `fixtures.checkpoint_model` raises `resident_memory_bytes` for this reason and
   documents it in place.

### F5 — shared corpus metadata is now stale. Reported, NOT changed.

`../spark-atomistic-rs/tests/corpus/philox_errata1.json` carries `"executed": false` and
`e2_parity_manifest.json` carries `"dynamic_execution": "NOT_RUN"`. Both are now false
for the Python backend: the Philox boundary, zero-key/counter, derivation, carry, and
restore cases are executed by `PhiloxCorpusTests`. Those files are shared with the Rust
backend and were left unmodified.

### F6 — historical pre-Errata-3 cross-language run: 16 of 46 fixtures diverged. RESOLVED.

D-127 adopted Errata 3. The current 2026-08-12 gate is **85/85 cases and 46/46
fixtures**, with zero divergence (core 83, adapter 2). The record below is retained as
the 2026-08-11 before-state.

`E2-PAR-005` was previously untested in both directions: each backend had only been run
against its own suite. `tests/xlang_harness.py` now feeds the shared corpus
`../spark-atomistic-rs/tests/corpus/xlang/` to both backends and compares the emitted
canonical JSON byte-for-byte. Recorded run, 2026-08-11, 85 artifacts per backend:

| | count |
|---|---:|
| cases byte-identical | 51 |
| cases divergent | 34 |
| cases with a missing artifact | 0 |
| mandatory fixtures byte-identical | 30 of 46 |

Everything scientific agrees byte-for-byte: candidate/state identities and geometry
certificates for the translated, periodic-image and same-species-permuted cases; every
Philox block, uniform, raw bit pattern, stream derivation, counter carry and substream
digest; search IDs and the deterministic discovery class sequence; pair/event/reverse ID
algebra including unstable-direction sign canonicalisation; `COMMON_PREFACTOR` rates,
detailed-balance residuals and the negative-barrier tolerance band; rate-table snapshots
including the lost-rate log upper bound and `payload_sha256`; KMC selections and
residence times; the capability value and the schema/config/tolerance/identity digests;
and all 29 status messages and all 29 exit codes.

Every divergence falls into five groups, none of which is a scientific disagreement:

1. `context.component` — this backend names the failing subsystem (`json`, `schema`,
   `adapter`), the Rust backend hardcodes `"api"`. 22 cases differ in this field alone.
   `E2-STATUS-001` fixes the key set but no requirement fixes the vocabulary.
2. `context.details` — this backend populates it (`path`, `missing`, `unknown`,
   `required_keys`, `extension_field`), the Rust backend always emits `{}`.
   `E2-STATUS-001` says only that `details` is an object.
3. `severity` — the Rust enum kebab-cases 12 of the 29 `E2-STATUS-002` severity strings
   (`candidate reject` vs `candidate-reject`, `pause/qualified` vs `pause-qualified`,
   `transaction fail` vs `transaction-fail`, `fatal in strict mode` vs
   `fatal-in-strict-mode`, `terminal-success if requested, else fatal` vs
   `terminal-success-if-requested-else-fatal`). This one **is** decidable: the
   `E2-STATUS-002` table is normative and this backend reproduces it verbatim.
4. Checkpoint restore — this backend restores the shared 28439-byte checkpoint carrying
   three reciprocal pairs and one committed step; the Rust backend refuses it with
   `CHECKPOINT_CORRUPT` (defect `D-E2-02` there). The paired control, the same model with
   an empty catalog, is accepted by both and is byte-identical, so the divergence is
   caused by the catalog events and nothing else.
5. `run` — the two backends bind different adapter transports (`E2-SCOPE-004` puts the
   transport outside the frozen surface), so a `run` request cannot agree: this backend
   requires `extension.calculator_command` and returns `INVALID_INPUT`/`E2-SCOPE-004`,
   the Rust backend requires a caller-supplied `RunAdapter` and returns
   `CALCULATOR_FAILURE`/`E2-API-003`.

Two defects were found by the harness rather than by either suite:

* **D-E2-03 (Rust only).** `serde_json` is pinned without `float_roundtrip`, so its
  number parser is not exactly rounding: canonical `9.999999999999997e-7`
  (`0x3eb0c6f7a0b5ed8c`) reparses as `0x3eb0c6f7a0b5ed8b` and `1.0000000000000001e+21`
  (`0x444b1ae4d6e2ef51`) as `0x444b1ae4d6e2ef50` — both are mandatory `E2-PAR-002`
  item-4 values. The Rust backend cannot re-read 4 of its own 85 artifacts. This backend
  re-reads all 85, asserted by
  `test_the_backend_can_reread_every_artifact_it_emits`.
* **D-E2-04 (both backends).** `E2-CAN-004` emits `9.999999999999999e20` as the digit-form
  token `999999999999999900000`, and `E2-JSON-002` classifies that token as an integer
  outside the portable domain and refuses it on input. Both backends do both things, so
  parity holds and the tension is between the two requirements. Recorded by
  `test_canonical_binary64_output_can_reenter_the_integer_domain`.

The `E2-JSON-002` domain itself is now identical in both backends: the Rust defect
`D-E2-01` (integer literals from 2^64 upward accepted and silently rounded) was fixed on
2026-08-11, and this backend was tested explicitly for the same hole and does not have it
(`PortableIntegerDomain`, five tests, every rejection paired with an accepted baseline).
After the fix both backends return `INVALID_INPUT` / exit 64 / `input invalid` /
`requirement_id` `E2-JSON-002` for all seven out-of-domain literals; only
`context.component` still differs.

### F7 — the commit-time detailed-balance gate was an algebraic identity. FIXED (gate deleted).

`catalog.validate_candidate` compared
`residual = (log_forward - log_reverse) + beta*(E_j - E_i)` against
`kinetics.detailed_balance_tolerance` and raised `DETAILED_BALANCE_VIOLATION`. Both legs
come from one saddle energy by `E2-RATE-001` ("Raw barriers are exactly `b_f=E_s-E_i` and
`b_r=E_s-E_j`") and `E2-EVENT-001` ("Raw same-saddle differences"), and `E2-RATE-002`
fixes the rates to those barriers, so `log_forward - log_reverse = beta*(E_i - E_j)`
identically and the residual is algebraically zero.

Measured through the real `Catalog.validate_candidate`, 20000 randomised physically valid
candidates (T in {100,300,700,1200} K, nu in {1e12,1e13,6.2e12} s^-1, |dE| <= 3 eV,
|E| <= 500 eV):

| tolerance | committed | `DETAILED_BALANCE_VIOLATION` | worst \|residual\| |
|---|---:|---:|---|
| `1e-8` (corpus) | 20000 | 0 | `1.1368683772161603e-13` |
| `1e-18` (below the rounding floor) | 8171 | 11829 | `0.0` on the survivors |

The status can therefore only ever reject valid physics on binary64 rounding. The shared
corpus says the same thing independently: `xlang/probes.json` reaches
`DETAILED_BALANCE_VIOLATION` in exactly one `rate_pair` case, which differs from its
accepted neighbour only by `detailed_balance_tolerance: 0` and whose "violation" is
`-8.881784197001252e-16`.

Repair option (a) — an independently obtained reverse barrier — is not reachable:
`E2-EVENT-001` requires `destination_state_id` to be a COMMITTED state and `E2-RATE-001`
pins `b_r` to that state's energy, which is the quantity already used, so the only
alternative would be a second saddle search or a Hessian prefactor, i.e. new physics no
requirement asks for. Option (b) was taken and the gate is deleted. The
`detailed_balance_residual` FIELD is retained because `E2-EVENT-003` fixes the
`rate_model` key set exactly; it is documented in source as an information-free column.

The verification that does have power is at restore: `E2-CKPT-007`(5) makes
`checkpoint._restore_event` recompute the residual from the stored barriers and energies
and refuse a stored value that disagrees, because there the stored value is an
independent input. Inside that same predicate the
`abs(expected_residual) > detailed_balance_tolerance` term is the same dead identity (the
two lines above it already pinned the barriers to the same-saddle differences); it is
annotated in source and left in place rather than removed unilaterally, since the
tolerance is the specification's own and the restore predicate is shared wire behaviour.

Tests: `test_rate_detailed_balance.py`, 8 tests. Mutation-verified — reinstating the
deleted gate turns `test_a_residual_tolerance_gate_could_only_reject_valid_physics` red
(446 of 750 valid candidates rejected). Every restore tamper re-seals the catalog digest,
and `test_reseal_alone_is_not_what_makes_a_tamper_pass_or_fail` proves the tamper reaches
`E2-EVENT-003` instead of being answered by `E2-CKPT-006`.

### F8 — the saddle-order gate had no statistical power. FIXED (sample replaced by a bound).

`_orthogonal_curvatures` reported one Rayleigh quotient `d^T H d` per configured
direction along a RANDOM `d` perpendicular to the mode. That quantity is a positively
weighted average of the orthogonal spectrum and is dominated by the stiff modes, so
`SADDLE-005`'s "nonnegative sampled orthogonal curvatures" was being certified by a
measurement that cannot distinguish index 1 from index 7, while `E2-EVENT-002` recorded
`evidence_level = "DIRECTIONAL"` on it.

Measured on the real method, 5 movable atoms (15 dof), Householder-rotated eigenbasis,
500 independent restarts per cell:

| spectrum (mode first) | true answer | before | after |
|---|---|---|---|
| `[-1.0, -0.5, +5.0 x 13]` (index 2) | reject | accepted 500/500 at 1, 3 and 10 directions; lowest sample `+1.869633` | rejected 500/500; bound reaches `-0.500000` |
| `[-1.0, -0.5 x 6, +5.0 x 8]` (index 7) | reject | accepted 500/500 at 1 and 3, 498/500 at 10 | rejected 500/500 |
| `[-1.0, +5.0 x 14]` (index 1) | accept | accepted 500/500 | accepted 500/500, bound `+5.000000` |
| `[-1.0, -0.5 x 14]` (local maximum) | reject | rejected 500/500 | rejected 500/500 |

The gate was not vacuous before — it could see a local maximum — it could only see the
case where EVERY orthogonal direction is soft. A negative sample of the index-2 cell
needs the random direction to place `c1^2 > 5/5.5 = 0.909` of its norm on the single soft
mode; 200000 draws produced none.

The replacement minimises the Rayleigh quotient over the complement with a second
dimer-style rotation confined to `d` perpendicular to `mode`: the projected residual
`P(Hd) - lambda d` gives the in-plane ascent direction, `dC/dt(0) = 2|residual|` is known
analytically, one trial angle fixes the Fourier model `a0/2 + a1 cos 2t + b1 sin 2t`, and
its minimising angle is taken. Every reported number is still a central-difference
directional curvature, so nothing is claimed that was not measured; it is now an upper
BOUND on the lowest curvature of the complement instead of one draw from it.

Cost is bounded and paid only while the bound moves: 2 evaluations for an isotropic
complement, 6 for the index-2 cell, at most `2 + 4*40 = 162` per restart. The budget is
what sets tightness on a deliberately clustered complement (`lambda_min = -0.05` against
`+0.3` and `+9.2`), measured over 200 restarts: cap 6 rejects 4, 12 rejects 52, 24 rejects
155, 40 rejects 186, 60 rejects 196. The cap is 40. A too-small cap under-reports a
negative curvature and never invents one, so the failure direction is a false ACCEPT, and
that residual limitation is real and unresolved for near-degenerate soft clusters.

The analytic double well is unaffected: launched from `x0 = 0.30` with
`force_tolerance = 1e-4` the search still converges, saddle energy
`1.0000000017345347` eV, curvature `-3.9999959984857427`, orthogonal curvature
`1.9999999999999998`, forward and reverse barriers `1.0000000014` eV, i.e. `1.00000000`
to 8 decimals, evidence `DIRECTIONAL`.

Tests: `test_saddle_order_gate.py`, 7 tests. Mutation-verified — restoring the
random-direction sampler produces 55 failures in that module.

### F9 — `_verify_application` compared the destination with itself. FIXED.

The verification relaxation was launched from `_state_as_system(destination)`, the
destination's own coordinates. `RELAX-003` had already committed that state at
`max_movable_force <= relaxation.force_tolerance`, so the minimiser returned at step 0
and the assertion was `X == X`.

Measured before the repair, three mutually non-matching committed states of the analytic
double well and the event id `event:fabricated`, which is in no catalog: all 9
(current, destination) combinations ACCEPTED after exactly 1 calculator evaluation, and
the returned state was a pure function of `destination_state_id`.

`EVENT-004` says "Event application MUST use the validated destination minimum, followed
by one verification relaxation. Failure to recover the destination within state
tolerances returns `EVENT_APPLICATION_FAILED`" — the validated destination minimum is
what the relaxation must RECOVER, not where it starts. The geometry is now rebuilt the
way `catalog.validate_candidate` built the endpoints: committed saddle positions
displaced along the committed unstable direction by
`saddle_search.endpoint_displacement`. Both signs are attempted because no requirement
fixes the sign of `saddle.unstable_direction` relative to the destination
(`E2-EVENT-006` canonicalises it only inside the `pair_id` hash) and
`validate_candidate` does not record which endpoint it kept.

After the repair, same fixture: the 9 fabricated-event combinations are refused with
`EVENT_APPLICATION_FAILED` before any calculator work (0 evaluations); a genuine A->B
event is accepted after 48 evaluations; single-field fabrications of
`destination_state_id`, `unstable_direction` and `saddle_positions` are each refused (96,
104 and 30 evaluations). `CALCULATOR_FAILURE`, `RESOURCE_LIMIT` and `CANCELLED` propagate
unchanged instead of being relabelled, matching the Rust twin's `application_failure`.
The emitted `context.requirement_id` stays `E2-KMC-003` so the public response string
does not move relative to the other backend, as with the F1 repair.

Tests: `test_event_application.py`, 8 tests. Mutation-verified — restoring the old body
produces 15 failures in that module.

## Deliberately not done

- No end-to-end `run` operation is exercised. The user testing policy forbids running
  anything under `SPARK/examples/`, and the solver diagnosis was required to be done
  in process rather than through the API. `ProcessCalculator`'s subprocess transport is
  therefore covered only by its response-validation unit assertions.
- No timing or performance assertion exists anywhere in the suite (`VAL-020`,
  user policy).
- The shared corpus under `../spark-atomistic-rs/tests/corpus/` was read and referenced
  but never written.
- The `E2-CKPT-007` order deviation (F2) and the corpus problems (F4, F5) were not
  changed unilaterally.
- `xlang_harness.py`'s `rate_pair` emitter still applies the detailed-balance tolerance
  that F7 deleted from `catalog.validate_candidate`. It is a wire-parity fixture
  generator, not the engine path, and the shared `xlang/probes.json` contains a case
  (`detailed_balance_tolerance: 0`) whose expected artifact is
  `DETAILED_BALANCE_VIOLATION`; removing the emitter's branch would create a new
  divergence against a backend this pass may not modify. `RATE-006` and the
  `E2-STATUS-002` row for the status are likewise left in the taxonomy.
- No `orthogonal_directions`-style configuration key was added for the F8 rotation
  budget. `E2-SCOPE-003` refuses unknown keys in every normative object, so a new key
  would invalidate the shared `e2_minimal_model.json` and move the config digest.
