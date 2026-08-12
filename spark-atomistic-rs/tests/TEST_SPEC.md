# Test specification

Status: executed. `cargo test` runs 57 tests and all 57 pass — 12 in `conformance_spec.rs`, 38 in
`errata2_parity.rs`, 1 in `xlang_emit.rs`, 5 in `run_e2e.rs`, and 1 in `run_parity.rs`. A green run is
evidence for the requirements the tests assert and for nothing else; see `../REQUIREMENTS.md` for the
defects and findings that execution produced.

Errata 2 mandatory coverage is frozen in `ERRATA2_PARITY_SPEC.md` and `corpus/e2_parity_manifest.json`
and is driven by `errata2_parity.rs`. The manifest records a per-fixture execution status and is
checked against the suite by `parity_manifest_matches_this_suite`, so the two cannot drift apart.
Of the 46 mandatory fixtures in this backend's own suite: 37 `PASS`, 8 `PARTIAL`, 0 `FAIL_DEFECT`,
1 `BLOCKED`.

`xlang_emit.rs` is the Rust half of the cross-language experiment `E2-PAR-005` demands. It drives
`corpus/xlang/` and writes 85 canonical artifacts to `$SPARK_XLANG_OUT/rust/`; the Python half and
the byte-for-byte comparator are `../../spark-atomistic/tests/xlang_harness.py`. Recorded run
2026-08-11, re-run after `src/run.rs` gave the corpus an adapter-backed `run`: 52 of 85 cases
byte-identical, 33 divergent, 0 missing; 31 of 46 fixtures identical. It also records, in
`_roundtrip.json`, the artifacts this backend cannot re-read (`D-E2-03`).

`run_e2e.rs` and `run_parity.rs` execute the run adapter against the shared analytic calculator
`../../spark-atomistic/corpus/mock_calculator.py`. They spawn one child process per energy/force
evaluation and therefore dominate the suite's wall time; no timing figure is asserted anywhere
(`VAL-020`, user policy). Neither test can silently skip: the calculator fixture's absence is an
assertion failure, and `run_parity.rs` falls back to its own temporary directory when
`SPARK_RUN_PARITY_DIR` is unset rather than doing nothing.

## Executed

- `E2-API-001..008`, `E2-SCOPE-003/004`, `E2-SCHEMA-001..012`: the three-operation wire surface,
  exact response and value key sets, metadata and extension digest exclusion, unknown-key rejection.
- `E2-JSON-001..003`, `E2-CAN-001..007`: duplicate keys at depth 0 and 2, malformed UTF-8, BOM, lone
  surrogates, nonfinite tokens and overflow, the portable-integer boundary, RFC 8785 numbers across
  the negative-zero, 1e-6 and 1e21 neighbourhoods, escaping rules, and the schema descriptor.
- `E2-ID-001..007`: translation, periodic-image and same-species-permutation invariance with negative
  controls; pair/event ID algebra including unstable-direction sign canonicalisation.
- `E1-DET-002-B/C/D`, `E2-RNG-001..007`: both spec-stated boundary golden-line SHA-256 values, the
  zero-key/zero-counter reference block, trajectory and saddle stream derivation against independently
  computed key/counter words, state restore, and 128-bit counter carry.
- `E2-RATE-001/002`, `RATE-004/005/006`: the negative-barrier tolerance band with a rejected control,
  detailed-balance violation with an accepted baseline, lost-rate log upper bound.
- `E2-DISC-001/002/003`: deterministic class selection with the real `class-selection` substream,
  seed dependence, non-consumption of the trajectory stream, and completion-order-independent commit.
- `E2-KMC-001/002/003/005`: snapshot ordering, `exp`/`log` bit equality, compensated total,
  historical immutability, and the two-uniform propose/rollback/commit boundary.
- `E2-CKPT-001..009`, `E2-BASIN-001/002/003`: zero-step envelope round trip with stable bytes,
  corrupt hash, non-canonical bytes, incompatible digests, recursive corruption at levels 1-4 and
  6, 8, 9, 10, pause flags, resource limits, and the exact `v1-disabled` basin record.
- `CALC-001..006`, `RELAX-001..004`, `SADDLE-001..007`, `DISC-001..006`, `EVENT-001..004`,
  `CAT-001..007`, `KMC-001..005`, `E2-API-003`, `E2-SCOPE-004` (`run_e2e.rs`, `run_parity.rs`): a
  complete run over the shared analytic double well. Measured against the closed form — saddle at
  `x = 0`, `E = 1.0` eV above the minimum — the search launched from `x0 = 0.30` with
  `force_tolerance = 1e-4` returns forward and reverse barriers of `1.000000000000` eV; the paired
  rejection launched from `x0 = 0.90`, identical in every other respect, returns `SADDLE_NOT_FOUND`
  because the Hessian there is `diag(5.72, 2, 2)` and the lowest mode is transverse. The end-to-end
  run commits 2 minima, 1 reciprocal pair with both barriers at 1.0 eV, and 1 KMC step, then is
  refused its own checkpoint by `D-E2-02`.

## Not executed, and why

- `VAL-001..003`, `VAL-005`, `VAL-014`: need multi-exit and multi-well fixtures. The bundled
  adapter (`src/run.rs`) can now drive them, but only the single-exit double well has been run.
- `VAL-004`: executed in the sense that a periodic-free analytic potential's minimum, saddle,
  endpoints and both barriers match the closed form; the periodic-surface and wrapped-identity half
  of the fixture is not run.
- `VAL-006`: the wire-level invariance is executed through `parity::state_ids`; the internal
  `identity::{closest_periodic_vector, discrete_identity, match_states}` path still lacks golden states.
- `VAL-007`: the ambiguity band is unreachable on the Errata 2 wire, where `environment_key` is fixed
  to `disabled` and `environment_version` to `none/1`.
- `VAL-008`: parallel saddles are distinguished at the ID level and `CAT-001` deduplication is
  exercised by the run (1 new event plus 1 duplicate from 2 successful searches); the analytic
  fixture has only one mechanism, so retaining two parallel events in one catalog is untested, and
  persisting either is blocked by defect `D-E2-02`.
- `VAL-009..012`: discovery config-digest mutation is executed by `conformance_spec.rs`; receipt
  tolerance/budget tamper, scoped overrun, and bit-exact multi-step replay need a live run loop.
- `VAL-013`: both backends have been run against the same corpus (`corpus/xlang/`), and the result
  is recorded rather than claimed: 31 of 46 mandatory fixtures agree byte-for-byte and 15 do not, so
  `E2-PAR-005` execution conformance is refused. A separate RUN comparison over one shared model is
  recorded in `../REQUIREMENTS.md`: all four model-derived digests agree, every emitted artifact
  diverges, and the `PAR-004` geometry tolerance of 1e-12 A is missed by 2.7e6 times because no
  requirement pins a minimiser (`D-E2-08`).
- `VAL-015..017`: basin acceleration is always disabled until an internal completeness validator
  exists; there is no exact-basin claim to test.
- `VAL-018..020`: publication, statistical, and performance validation are not run. No benchmark was
  executed and no timing claim is made.

No validation, production, release, performance, parity, or publication claim is made. The
cross-language run is a measurement, not a parity claim.

Additional negative corpus still required: discovery-digest decode mismatch, leading `../..` path
components, active-region search-ID changes, and symlink/hardlink output aliases. These exercise
`config`/`checkpoint` paths that the Errata 2 wire surface does not reach.
