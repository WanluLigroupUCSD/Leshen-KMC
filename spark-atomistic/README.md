# spark-atomistic

Independent Python reference for fixed-composition off-lattice/on-the-fly KMC.

Public API:

```python
from spark_atomistic import (
    capabilities,
    run_atomistic_reference_unvalidated_json,
    validate_atomistic_model_json,
)
```

The IR identifier is `spark-atomistic-model/1`. `validate_atomistic_model_json`
returns a canonical JSON response and performs no scientific work.
`run_atomistic_reference_unvalidated_json` refuses to run unless the request object
carries `"allow_unvalidated": true`. Its only calculator transport is a
process-isolated JSON-lines executable with enforced timeout, process-group
kill-and-wait, and separate stdout/stderr byte caps.
The executable schema is the single strict validator in `model.py`; validation and
runtime call that same validator, so no separate drifting schema artifact exists.
Its complete language-neutral nested descriptor is revision 2 and exports a real
canonical SHA-256 as `SCHEMA_SHA256`; the same digest is embedded in capabilities,
validation responses, catalogs, and checkpoints.

## Status

The implementation is reported as `implemented_unvalidated` and `partial`.

It **has** been imported, executed, and tested: `python3 -m unittest discover -s tests
-p 'test_*.py'` runs 100 tests: 99 passed and 1 cross-language comparison skipped when
`SPARK_XLANG_OUT` is absent. See `tests/TEST_SPEC.md` for the run record, the requirement-ID coverage
map against `E2-PAR-002`, and every finding with its measured numbers.

It **has** now been run against the Rust backend at runtime. `tests/xlang_harness.py`
feeds the shared corpus `../spark-atomistic-rs/tests/corpus/xlang/` to both backends and
compares the emitted canonical JSON byte-for-byte:

```text
python3 tests/xlang_harness.py emit    --out DIR
python3 tests/xlang_harness.py compare --out DIR   # nonzero exit on any divergence
```

Recorded Errata 3 run 2026-08-12: **85/85 cases and 46/46 mandatory fixtures pass**.
The tier split is core 83 and adapter 2. `tests/xlang_harness.py compare` gates core
unconditionally and reports adapter separately.

It has **not** been benchmarked, scientifically validated, or released. The status stays
`implemented_unvalidated`: parity establishes a shared contract, not physical accuracy.

Exact finite-basin acceleration is disabled by default and cannot be enabled under this
IR revision (`E2-BASIN-003`). Enabling it in the model is schema-valid but records
`BASIN_DISABLED`, preserves the RNG, and continues serially. The basin capability is
therefore `partial`, not a validated exactness claim.

Retry is not implemented: `retry_count` and `retry_backoff_s` must both be zero.

## Behavioural notes

**No detailed-balance verification happens when an event is committed.** `E2-RATE-001`
fixes both barriers to one saddle energy and `E2-RATE-002` fixes the rates to those
barriers, so the `detailed_balance_residual` of `E2-EVENT-003` is an algebraic identity:
measured worst |residual| over 20000 randomised physical inputs is
`1.1368683772161603e-13`, pure binary64 rounding, against a configured tolerance of
`1e-8`. The field is still emitted because `E2-EVENT-003` fixes the `rate_model` key set
exactly, but it is an information-free column and no gate is applied to it. Detailed
balance IS recomputed at checkpoint restore (`E2-CKPT-007`(5)), where the stored residual
is an independent input and a tampered value is refused. Finding F7 in
`tests/TEST_SPEC.md` has the numbers.

The orthogonal saddle-order evidence of `SADDLE-005` is a BOUND, not a sample: per
configured restart the Rayleigh quotient is minimised over the complement of the unstable
mode by a second dimer-style rotation, and the minimum found is reported. A single random
direction — what was reported before — accepted a genuine index-2 saddle 500/500 times
and a genuine index-7 saddle 500/500 times; the bound rejects both 500/500 while still
accepting index 1 500/500. The bound is honest but not tight for near-degenerate soft
clusters, where the rotation budget binds; it under-reports negative curvature and never
invents it. Finding F8.

Event application reconverges the selected event's own `saddle ± endpoint_displacement`
geometry and matches the result against the validated destination (`EVENT-004`). It does
not relax the destination's own coordinates, which is what made the check `X == X`
before. Finding F9.

Each physical saddle is committed atomically as two immutable directed records whose
`reverse_event_id` values point to each other. State matching first checks a
geometry-free fixed contract, then performs full configured-tolerance triclinic
geometry and energy verification. Candidate hashes are hints and never hard filters.

Every committed step stores a hashed full rate-table snapshot; checkpoint replay uses
that historical snapshot, never the later expanded catalog (`E2-KMC-005`). This is
covered by a mutation-verified test: rebuilding the table from the restored catalog
makes the accepted-baseline test fail.

Output bytes are a reserved resource. A checkpoint write drives the self-referential
`output_bytes` field to a fixed point, reserves exactly that many bytes, and only then
writes through a sibling temporary, content flush, atomic replacement, and
parent-directory flush (`E2-CKPT-009`). A refused reservation leaves nothing on disk;
a failure at either write stage leaves the previous artifact byte-identical and no
temporary behind.

Resource checkpoints retain callback, solver, and endpoint attempt audits with
calculator-evaluation deltas, final statuses, and termination reasons. Calculator
cleanup kills the process group and uses only bounded waits.

## Known limitation: saddle search domain

`DirectionalDimerSearcher` is a minimum-mode follower. On the analytic double well
`V = Σᵢ (xᵢ²−1)² + yᵢ² + zᵢ²` it reproduces the exact answer — saddle at x = 0,
E = 1.0 eV, forward and reverse barriers 1.000 eV — when launched from a point whose
Hessian still has a negative eigenvalue and when `force_tolerance` is tight enough
(achieved accuracy is bounded by `force_tolerance²/8`).

Launched from (0.9, 0, 0) the Hessian is diag(5.72, 2, 2): every eigenvalue is positive
and the lowest mode is transverse, so the search climbs the transverse mode instead and
returns `SADDLE_NOT_FOUND` after exhausting its iteration budget. This is a
characterised parameter-domain limitation, not a specification violation:
`SADDLE-002` permits any first-order-saddle solver and `E2-STATUS-002` classifies
`SADDLE_NOT_FOUND` as a non-terminal candidate reject. Full numbers, the measured
success rates, and the one escalated consequence are in `tests/TEST_SPEC.md` finding F3.

## Tests

```text
python3 -m unittest discover -s tests -p 'test_*.py'
```

`tests/fixtures.py` puts `src/` on `sys.path`, so no install and no `PYTHONPATH` is
needed. The shared Errata-2 golden corpus lives at `../spark-atomistic-rs/tests/corpus/`
and is referenced, never copied (`E2-PAR-001`).

`test_cross_language_parity.py` runs the Python half unconditionally. Its byte-for-byte
comparison against the Rust artifacts runs only when `SPARK_XLANG_OUT` points at a
directory holding them, and SKIPS with the exact producing command otherwise, so an
unavailable comparison is never reported as a clean one:

```text
cd ../spark-atomistic-rs && SPARK_XLANG_OUT=<dir> CARGO_TARGET_DIR=<tmp> \
    cargo test --test xlang_emit
cd ../spark-atomistic && SPARK_XLANG_OUT=<dir> python3 -m unittest discover -s tests -p 'test_*.py'
```

Package license: [PolyForm Strict License 1.0.0](LICENSE), limited by
[LICENSE_SCOPE.md](LICENSE_SCOPE.md). It permits noncommercial use but not
distribution, modification, or derivative works. It is source-available, not
open source. See `PROVENANCE.md`.
