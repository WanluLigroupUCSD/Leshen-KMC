# Provenance

## Design inputs

- Independent implementation date: 2026-08-09. Test and closeout pass: 2026-08-11.
- Base specification read path: `/home/shidi/ai-chemist/catgo-projects/data/ai-kmc/specs/OFFLATTICE_OTF_KMC_SPEC_V1.md`.
- Base self-excluding SHA-256: `8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84`.
- Erratum 1 read path: `.../OFFLATTICE_OTF_KMC_SPEC_V1_ERRATA_1.md`.
- Erratum 1 SHA-256: `52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40`.
- Erratum 2 read path: `.../OFFLATTICE_OTF_KMC_SPEC_V1_ERRATA_2_PARITY.md`.
- Erratum 2 SHA-256 (frozen, normative): `eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995`.
- Erratum 3 read path: `.../OFFLATTICE_OTF_KMC_SPEC_V1_ERRATA_3.md`.
- Erratum 3 SHA-256: `eba384af3694c5f3997caf28829e56d188ef9929f29d99a2116520f0067d8a96`.
- All four self-excluding digests match.
- Clean-room exceptions: zero.
- Other source, history, tests, symbols, layouts, papers, or competitor material read: none.

## Dependencies

- Runtime dependencies: Python standard library only; Python license; dynamic runtime linkage.
- Build dependency: `setuptools>=68`; MIT license; build-time only.
- Package license: PolyForm Strict License 1.0.0, with the package boundary stated in
  `LICENSE_SCOPE.md`. It is source-available, not open source. Dependency compatibility,
  contributor/institutional authority, and patent clearance remain release gates.

## Validation state

**`implemented_unvalidated`.**

Superseded statement: earlier revisions of this file recorded "static inspection only;
no import, execution, test, benchmark, or compile performed". That is **no longer
true** and is corrected here.

What has been performed:

- Package import of all 15 modules under `src/spark_atomistic/`.
- Executed test suite: `python3 -m unittest discover -s tests -p 'test_*.py'` —
  **100 tests: 99 passed, 0 failed, 1 skipped** without `SPARK_XLANG_OUT`.
- Executed tiered cross-language gate: **85/85 cases, 46/46 fixtures**; core 83,
  adapter 2.
- Byte-identity of the capability response against the shared golden
  `../spark-atomistic-rs/tests/corpus/e2_capabilities.response.json`.
- Canonical-number goldens 7/7; Philox zero-key/zero-counter words and both Errata-1
  uniform boundary bit patterns exact, including the two normative line hashes of
  `E2-RNG-006`/`E2-RNG-007`.
- In-process numerical characterisation of `DirectionalDimerSearcher` on the analytic
  double well (`tests/TEST_SPEC.md` finding F3).
- Mutation verification of the `E2-KMC-005` historical-snapshot replay assertion.

What has **not** been performed:

- No benchmark and no timing measurement of any kind (`VAL-020`, and user policy).
- Rust compilation/tests and dynamic cross-language parity were run. This does not
  establish scientific validation.
- No end-to-end `run` operation and nothing under `SPARK/examples/` was executed
  (user policy).
- No scientific validation and no release.
- `E2-PAR-002` categories 8 (deterministic discovery class sequence, out-of-order
  completion commit) and parts of 5, 7, 9, 10, 11 have no executed Python test yet;
  the gaps are enumerated in `tests/TEST_SPEC.md`.

## Change record for this pass (2026-08-11)

Source change, one file:

- `src/spark_atomistic/engine.py` — `_validate_output_preflight` now applies the
  `resume=true` exemption to the checkpoint path only, per `IO-004` ("Existing
  compatible checkpoint plus `resume=true` resumes; every other collision returns
  `OUTPUT_EXISTS`"). Previously a pre-existing `summary_path` or `trajectory_path` was
  silently tolerated and would have been overwritten at run completion without
  `overwrite=true`. The emitted `context.requirement_id` is unchanged
  (`E2-SCHEMA-010`) so the public response string does not move relative to the other
  backend. Regression-confirmed: reverting the predicate fails three subtests.

## Change record for the defect-repair pass (2026-08-11)

Three source files, three new test modules, one fixture module extended. Full numbers,
before and after, in `tests/TEST_SPEC.md` findings F7, F8 and F9.

- `src/spark_atomistic/catalog.py` — the commit-time `DETAILED_BALANCE_VIOLATION` gate is
  deleted. `E2-RATE-001` fixes `b_f=E_s-E_i` and `b_r=E_s-E_j` to one saddle energy and
  `E2-RATE-002` fixes the rates to those barriers, so the residual of `E2-EVENT-003` is an
  algebraic identity (worst |residual| `1.1368683772161603e-13` over 20000 randomised
  physical inputs against a `1e-8` tolerance) and the gate could only ever reject valid
  physics on rounding noise (11829 of the same 20000 rejected at a `1e-18` tolerance). An
  independently obtained reverse barrier is not available without inventing physics,
  because `E2-EVENT-001` requires the destination to be a committed state and
  `E2-RATE-001` pins `b_r` to that state's energy. The `detailed_balance_residual` FIELD
  is retained: `E2-EVENT-003` fixes the `rate_model` key set exactly.
- `src/spark_atomistic/checkpoint.py` — comment only. The `E2-CKPT-007`(5) recomputation
  is kept because there the stored residual is an independent input, and the dead
  tolerance term inside the same predicate is annotated rather than removed, since the
  restore predicate is shared wire behaviour.
- `src/spark_atomistic/solvers.py` — `_orthogonal_curvatures` now minimises the Rayleigh
  quotient over the complement of the unstable mode (a second dimer-style rotation
  confined to `d` ⟂ `mode`) and reports the minimum found per restart instead of one
  random-direction sample. A random-direction sample accepted a genuine index-2 saddle
  500/500 and a genuine index-7 saddle 500/500; the bound rejects both 500/500 and still
  accepts index 1 500/500. No configuration key was added (`E2-SCOPE-003`).
- `src/spark_atomistic/engine.py` — `_verify_application` reconverges the selected event's
  `saddle ± saddle_search.endpoint_displacement · unstable_direction` and matches the
  result against the validated destination, per `EVENT-004`. It previously relaxed the
  destination's own coordinates, which returns at step 0 and asserts `X == X`: all nine
  (current, destination) combinations of three mutually non-matching states were accepted
  after one calculator evaluation, including for an event id in no catalog. The emitted
  `context.requirement_id` stays `E2-KMC-003` so the public response does not move
  relative to the other backend. `CALCULATOR_FAILURE`, `RESOURCE_LIMIT` and `CANCELLED`
  now propagate instead of being relabelled `EVENT_APPLICATION_FAILED`.
- `tests/fixtures.py` — additive only: a scoped `evaluation_scope` on `StubCalculator`
  matching `ProcessCalculator`'s contract, an N-atom double-well state builder, and an
  exact quadratic calculator with a prescribed Hessian spectrum.
- `tests/test_rate_detailed_balance.py`, `tests/test_saddle_order_gate.py`,
  `tests/test_event_application.py` — 23 new tests. Each repair is mutation-verified:
  restoring the old code turns the new module red (1, 55 and 15 failures respectively).

Deliberately **not** changed in the defect-repair pass:

- `tests/xlang_harness.py`'s `rate_pair` emitter, which still applies the deleted
  tolerance. It is a wire-parity fixture generator and the shared `xlang/probes.json`
  contains a case with `detailed_balance_tolerance: 0` whose expected artifact is
  `DETAILED_BALANCE_VIOLATION` (its "violation" is `-8.881784197001252e-16`). Removing the
  branch would create a new divergence against a backend this pass may not modify.
- The `DETAILED_BALANCE_VIOLATION` status itself, which `E2-STATUS-002`, `ERR-002` and
  `RATE-006` all name; it stays in the taxonomy with no reachable raise site on the
  commit path.
- Anything under `../spark-atomistic-rs/`. `src/kmc.rs` was READ as the reference design
  for the `EVENT-004` repair and not modified; nothing in the shared corpus was written.

Historical pre-Errata-3 differences in restore-order citations,
`context.component`, and `context.details` were retained until an owner decision. D-127
adopted Errata 3 and resolved them; the current tiered gate has zero divergence. Finding
F6 in `tests/TEST_SPEC.md` preserves the before-state.
- `DirectionalDimerSearcher`. Its non-convergence from a launch point with no negative
  Hessian eigenvalue is a characterised parameter-domain limitation; no requirement
  mandates convergence there.
- Anything under `../spark-atomistic-rs/src/`. The shared corpus under
  `../spark-atomistic-rs/tests/corpus/` was extended on 2026-08-11 with the
  cross-language directory `xlang/` and its per-fixture verdicts, because `E2-PAR-001`
  requires one corpus consumed by both backends and there was nowhere else it could
  live; nothing pre-existing in that directory was rewritten except
  `e2_parity_manifest.json`, whose Rust-side `json-integer-out-of-domain` verdict changed
  from `FAIL_DEFECT` to `PASS` when that defect was fixed. Two shared-corpus problems were
  measured and reported rather than edited:
  `e2_minimal_model.json` sets `resources.resident_memory_bytes` to 1 MB (21–26× below
  the measured resident set of a CPython process running this package, so any run using
  it fails immediately with `RESOURCE_LIMIT`), and it pairs
  `relaxation.force_tolerance = 0.05` eV Å⁻¹ with `kinetics.barrier_tolerance = 1e-10`
  eV, a 3.1 × 10⁶ mismatch that can produce a negative raw barrier and hence
  `RATE_INVALID`. `philox_errata1.json` (`"executed": false`) is still stale for the
  Python backend and was left unmodified.

## Integrity manifest

Per-file SHA-256 is written to `FILE_HASHES.sha256` in plain `sha256sum -c` format. It
covers every file under `src/spark_atomistic/`, plus `pyproject.toml`, the three
documentation files, `corpus/`, and `tests/`. It excludes itself and excludes
`__pycache__/` (regenerated bytecode is a build artifact, not source). Regenerate it as
the last step of a change:

```text
find src/spark_atomistic pyproject.toml README.md PROVENANCE.md LICENSE LICENSE_SCOPE.md corpus tests \
     -type f -not -path '*/__pycache__/*' | sort | xargs sha256sum > FILE_HASHES.sha256
sha256sum -c FILE_HASHES.sha256
```

**A refreshed hash manifest records what the bytes are; it is not evidence that they are
correct.** It is an integrity record only, and it carries no information about test
results, specification conformance, or cross-language parity.
