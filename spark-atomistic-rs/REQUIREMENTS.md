# Requirement coverage

Status: the implementation map below is static. Errata 2 observable behavior is executed by
`tests/errata2_parity.rs` (39 tests), `tests/conformance_spec.rs` (13 tests), `tests/xlang_emit.rs`
(1 test), `tests/errata3_candidate.rs` (5 tests), `tests/run_e2e.rs` (10 tests) and
`tests/run_parity.rs` (1 test); all 69 pass. Passing tests are evidence for the requirements they
assert and for nothing else. Cross-language dynamic parity against the Python backend has been RUN
over a shared corpus and the result is recorded, not claimed. A complete `run` operation is
implemented (`src/run.rs`) and executed against the shared analytic calculator.
The `VAL-001..020` scientific matrix and every performance question remain unrun.

**ERRATA 3 NORMATIVE (2026-08-12) — adopted by `D-127`.** Current tiered gate:

| tier | cases | result |
|---|---:|---:|
| core | 83 | 83 PASS |
| adapter | 2 | 2 PASS |
| mandatory fixtures | 46 | 46 PASS |

Adopted items are `E3-EVENT-001` options 1 + 4, `E3-PAR-001`, `E3-PAR-002`, `E3-CAN-001`, and
`E3-STATUS-001`/`E3-STATUS-002` option (a). `parity::CheckpointPolicy::PreErrata3` and
`parity::PayloadHashPolicy::PreErrata3Full` remain compatibility paths, not writer defaults.

Errata 2 overrides observable wire behavior. `src/parity.rs` owns `E2-API/SCHEMA/JSON/CAN/ID/RNG/EVENT/RATE/DISC/KMC/CKPT/BASIN/PAR`; adapter transport remains only in run-request `extension`.

| Contract | Module | Static coverage | Executed |
|---|---|---|---|
| `SCOPE-*`, `ATOM-*`, `STATE-*` | `model`, `identity`, `config` | Private validated state ctor; recomputed movable force/canonical ID/provenance; fixed ID→species/constraint binding. | Wire-level only, through `E2-SCHEMA-*` validation and `E2-ID-006` state records. |
| `ENV-*` | `identity` | Versioned invariant key, ambiguity failure, refinement provenance, reuse counters; specific validation remains mandatory. | No. `environment_key` is fixed to `disabled` on the Errata 2 wire. |
| `CALC-*`, `RELAX-*`, `SADDLE-*` | `callbacks`, `run` | Receipt revalidates request, exact force tolerance, convergence, steps/evaluations, callbacks; search ID includes sorted active region. `run` adds a process-isolated JSON-lines calculator, a backtracking steepest-descent minimiser and a minimum-mode-following saddle solver. | Yes, through `src/run.rs`: `tests/run_e2e.rs` drives the shared analytic double well and reproduces its closed-form 1.000 eV barriers. The `callbacks` receipt layer itself is still unexecuted. |
| `EVENT-*`, `CAT-*` | `catalog` | Request-bound endpoint receipts, two directed reciprocal records, undirected-pair/saddle dedup, recursive snapshot recomputation, epoch/digest rate cache. | ID algebra only (`E2-EVENT-006`). Catalog admission is blocked by `D-E2-02`. |
| `DISC-*` | `discovery` | Externally private statistics; exact config digest checked on mutation, checkpoint validation, and decode. | `E2-DISC-001/002/003` executed; `E2-DISC-004/005` counters exercised only through checkpoint validation. |
| `RATE-*` | `rate` | Common prefactor, shared saddle, log form, detailed balance, fail-closed extreme domain. Harmonic TST disabled. | Yes: `E2-RATE-001/002`, `RATE-004/005/006`. |
| `KMC-*` | `kmc` | Private trajectory, validated construction/resume, canonical event order, compensated sum, two-uniform clone/commit, destination reconvergence. | `E2-KMC-001/002/005` and the two-uniform rollback executed; the crate's own selection replay is blocked by `D-E2-02`. |
| `BASIN-*` | `basin` | Hard-disabled: no internal completeness validator, enabling receipt, checkpoint self-attestation, or exact claim. | Yes: `E2-BASIN-001/002/003`, `E2-SCHEMA-011`, `E2-CKPT-003`. |
| `DET-*` | `rng` | Philox4x32-10, Errata-1 lane pairing/midpoint, checked counter/buffer/consumption relation, deterministic substreams. | Yes: `E1-DET-002-B/C/D`, `E2-RNG-001..007` including both spec-stated golden-line hashes. |
| `CKPT-*`, `IO-*` | `checkpoint`, `config` | Exact expected discovery digest; `ParentDir` only removes `Normal`; canonical existing-parent/target same-file gate. | Zero-step envelopes only: `E2-CKPT-001..009` except the event/trajectory levels blocked by `D-E2-02`. |
| `RES-*` | `resource`, `callbacks` | Scoped request-bound exact actual/counter overrun evidence; checkpoint requires resource-limited flag; retry `0`. | Checkpoint-side counters and flags executed; live budget enforcement is not. |
| `ERR-*` | `status`, `checkpoint` | Exact closed status/severity/context parsing and exit mapping; cleanup cause retained; code `74` only without valid checkpoint. | `E2-STATUS-001/002/003` executed for the statuses the wire surface can reach. |
| `PAR-*`, `VAL-*` | `tests` | Corpus and spec. | `E2-PAR-001..005` partially: 37 of 46 mandatory fixtures pass this backend's own suite, and 31 of 46 agree byte-for-byte with the Python backend over the shared `corpus/xlang/`. A separate RUN comparison over one shared model is recorded below: every model-derived digest agrees, every emitted run artifact diverges. `VAL-001..020` remain unrun. |
| `LIC-*` | all files, `PROVENANCE.json` | Per-file clean-room headers, dependency record, SHA allowlist, zero-source attestation. | Manifest verified by `SPARK/tools/check_provenance.py`. |

| Errata 2 group | Static coverage | Execution |
|---|---|---|
| `E2-API/SCHEMA` | Exact three-operation JSON surface; strict wire model and cross-field gates. | Executed. |
| `E2-JSON/CAN/STATUS` | Portable integers, nonfinite classification, Unicode/duplicate rejection, shortest canonical numbers, exact status response vocabulary. | Executed. `D-E2-01` fixed; `D-E2-03` (parser is not exactly rounding), `D-E2-04` (both backends), `D-E2-05`/`D-E2-06` (unpinned context fields) and `D-E2-07` (severity strings) stand. |
| `E2-ID/RNG` | Anchor-minimum certificate, fixed/constraint/candidate/state/event IDs, revision-2 trajectory/saddle streams. | Executed; one recorded arithmetic limit (`F-E2-01`). |
| `E2-EVENT/RATE/DISC/KMC` | Exact wire records, raw negative barriers, cutoff/lost bound, deterministic class/search order, immutable historical snapshots. | Executed at constructor level; catalog admission blocked (`D-E2-02`). |
| `E2-CKPT/BASIN/PAR` | Revision-2 typed envelope and recursive restore, `v1-disabled`, complete mandatory fixture index. | Executed for zero-step payloads; event-bearing payloads blocked (`D-E2-02`), including the Python-produced checkpoint in `corpus/xlang/` that the Python backend restores. |

## Fixture execution

`tests/corpus/e2_parity_manifest.json` carries the authoritative per-fixture record and is checked
against this suite by `parity_manifest_matches_this_suite`, so the two cannot drift apart.

| Status | Count | Meaning |
|---|---:|---|
| `PASS` | 37 | an executing test asserts the fixture's requirement and passes |
| `PARTIAL` | 8 | the public-constructor half executes; the catalog/trajectory half is blocked |
| `FAIL_DEFECT` | 0 | executes, and the implementation does not meet the cited requirement |
| `BLOCKED` | 1 | no conforming input can be constructed, so it cannot execute at all |

Those counts describe this backend's own suite. `E2-PAR-005` execution conformance is a *cross*
-language property, and the manifest's `cross_language_execution` block carries that record
separately: **31 of 46 mandatory fixtures agree byte-for-byte with the Python backend; 15 diverge.**
See "Cross-language parity" below.

## Defects found by execution

**`D-E2-01` — `E2-JSON-002` portable-integer domain had a hole from 2^64 upward. Severity: high. FIXED 2026-08-11.**
`E2-JSON-002` requires that "a syntactically valid integer outside this domain returns
`INVALID_INPUT` before schema validation". The strict parser enforced the domain in its `i64` and
`u64` visitors only. Any integer literal too large for `u64` was routed by `serde_json` to the
binary64 visitor, which accepts every finite value. Measured before the fix: `9007199254740992`
(2^53) and `18446744073709551615` (2^64-1) were correctly rejected; `18446744073709551616` (2^64)
and `1000000000000000000000000000000` (1e30) were accepted and silently rounded to binary64. The
hole covered every integer literal with magnitude in [2^64, ~1.8e308]. "Silently" is exact: no
status, no warning, and where the literal was not exactly representable the stored value differed
from the input — the digit-form 1e30 became `1000000000000000019884624838656`, while 2^64 happens to
be exactly representable and was accepted unrounded. Free-form `metadata` and run-request
`extension` were the reachable carriers, since `E2-SCOPE-003` requires their values to obey the
portable JSON domain too.

Fix: the visitor guards cannot see the source token, because `serde_json` has already reduced a wide
literal to a binary64 by the time `visit_f64` runs. `checkpoint::portable_integer_violation` decides
the domain on the source text instead, using exactly the classification `E2-JSON-002`/`E2-JSON-003`
word: a number token with no fraction and no exponent is an integer and must lie in
`[-9007199254740991,9007199254740991]`; any other number token is a binary64. It is applied by
`parse_strict_json` (so checkpoint bytes are covered too) and by `parity::parse_public`, which cites
`E2-JSON-002` as the governing requirement rather than the `E2-JSON-001` the generic parse-failure
path used. Measured after the fix, on the public wire: `18446744073709551616`,
`-18446744073709551616`, a digit-form `1e30`, a 401-digit literal and `999999999999999900000` all
return `INVALID_INPUT`, exit `64`, message `input invalid`, `requirement_id` `E2-JSON-002`; the
in-domain endpoints `±9007199254740991` and `0`, the exponent forms `1e30` /
`1.8446744073709552e19` / `9.999999999999999e20`, and the string `"18446744073709551616"` are all
still accepted. Witness: `json_integer_out_of_domain` (replacing the defect witness), plus
`json_integer_boundary`. The Python backend never had this hole: it decodes integer literals as
arbitrary-precision `int` and applies `validate_portable_value` to the parsed tree. Verified
explicitly by `spark-atomistic/tests/test_cross_language_parity.py::PortableIntegerDomain`.

**`D-E2-03` — canonical binary64 output does not survive this crate's own parser. Severity: high.**
Found by the cross-language harness's self-check. `E2-CAN-004` requires "the shortest decimal that
round-trips to the same binary64 value"; the write side satisfies it. The read side does not.
`serde_json` is pinned at `=1.0.140` **without** the `float_roundtrip` feature, so its number parser
is not exactly rounding. Measured: `-11.999999984770021` (`0xc027ffffff7d2cec`) reparses as
`0xc027ffffff7d2ceb`; `9.999999999999997e-7` (`0x3eb0c6f7a0b5ed8c`) as `0x3eb0c6f7a0b5ed8b`;
`1.0000000000000001e+21` (`0x444b1ae4d6e2ef51`) as `0x444b1ae4d6e2ef50`. The last two are mandatory
`E2-PAR-002` item-4 corpus values. Accepted baseline that does round-trip:
`0.0000010000000000000002` (`0x3eb0c6f7a0b5ed8e`). Consequence: `E2-CKPT-007` step 1 re-canonicalises
the parsed payload and compares it to the input bytes, so a valid checkpoint carrying such a value is
refused as `CHECKPOINT_CORRUPT`, and `E2-PAR-003` byte-identity cannot be guaranteed. Measured blast
radius on the shared corpus: 4 of this backend's 85 canonical artifacts cannot be re-read
(`kmc-lost-rate-bound__snapshot`, `kmc-two-event-selection__snapshot`,
`kmc-two-event-selection__selection`, `kmc-historical-snapshot-growth__snapshot`); the Python backend
re-reads all 85. Witnesses:
`defect_e2_can_004_canonical_binary64_does_not_survive_this_crates_own_parser` and the
`_roundtrip.json` report written by `xlang_emit`. Not fixed here: the remedy is a dependency-feature
change (`serde_json` `float_roundtrip`) that alters the pinned dependency contract and `Cargo.lock`,
which is an owner decision, and the affected corpus checkpoints happen to round-trip so the D-E2-02
attribution below is unaffected.

**`D-E2-04` — `E2-CAN-004` can emit a token that `E2-JSON-002` refuses. Severity: medium. Specification tension, both backends.**
`E2-CAN-004` serialises binary64 by RFC 8785, which uses digit form for `1e-6 <= |v| < 1e21`. For
`|v|` in `[2^53, 1e21)` the result is a token with no fraction and no exponent, i.e. syntactically an
integer outside the portable domain, which `E2-JSON-002` requires to be refused on input. Measured:
`9.999999999999999e20` (`0x444b1ae4d6e2ef4f`, a mandatory `E2-PAR-002` item-4 value) serialises to
`999999999999999900000` in **both** backends and is refused on input by **both** backends. The
neighbouring value `1e21` crosses into exponent form (`1e+21`) and round-trips. Parity therefore
holds; the defect is in the requirement pair, not in either implementation. Witness:
`test_canonical_binary64_output_can_reenter_the_integer_domain` in the Python suite.

**`D-E2-05` — `context.component` has no specified vocabulary. Severity: medium. Both backends.**
`E2-STATUS-001` fixes the context key set but no requirement fixes the value space of `component`.
Python names the failing subsystem (`json`, `schema`, `adapter`, `state`, `checkpoint`); this crate
hardcodes `"api"` on every public failure. Because `E2-PAR-001`/`E2-PAR-003` require the public
response to be byte-identical, an unpinned field is a parity blocker: 22 of the 34 divergent
cross-language cases differ in this field alone. No requirement decides which is right.

**`D-E2-06` — `context.details` content is unspecified. Severity: medium. Both backends.**
`E2-STATUS-001` says only that `details` is an object. Python populates it (`path`, `missing`,
`unknown`, `required_keys`, `extension_field`); this crate always emits `{}`. Same parity
consequence as `D-E2-05`, same absence of a deciding requirement.

**`D-E2-07` — severity strings do not match the `E2-STATUS-002` table. Severity: medium. This crate.**
`E2-STATUS-002` gives the exact severity column. The `Severity` enum here is
`#[serde(rename_all = "kebab-case")]`, so 12 of the 29 status rows serialise differently from the
table: `candidate reject` -> `candidate-reject`, `pause/qualified` -> `pause-qualified`,
`transaction fail` -> `transaction-fail`, `fatal in strict mode` -> `fatal-in-strict-mode`,
`terminal-success if requested, else fatal` -> `terminal-success-if-requested-else-fatal`. All 29
exit codes and all 29 messages agree byte-for-byte with Python. This one **is** decided by a
requirement: the `E2-STATUS-002` table is normative and Python reproduces it verbatim, so this crate
is the nonconforming side. Recorded rather than changed here because it alters the public wire for
every nonterminal and strict-mode status and is therefore an owner-visible contract change.

**`D-E2-02` — a reciprocal directed-event pair cannot exist in a revision-2 checkpoint. Severity: critical (blocking).**
`parity::validate_checkpoint_v2` and `parity::validate_event_v2` jointly impose three constraints:

1. every directed record must satisfy `saddle.search_id == discovery_provenance.search_id`;
2. the reciprocal record's `saddle` object must equal the forward record's `saddle`, hence both
   records carry one `search_id`;
3. every record's `discovery_provenance.search_id`, and its substream key/counter, must re-derive
   from **that record's own** `origin_state_id` per `E2-DISC-002` and `E2-RNG-004`.

`E2-EVENT-001` requires the origin and destination to be distinct committed IDs, and `E2-DISC-002`
binds the search ID to the state ID, so constraint 3 forces two different search IDs while
constraints 1 and 2 force one. The set is unsatisfiable. Verified exhaustively: all 64 assignments of
`(saddle.search_id, discovery_provenance.search_id, provenance substream)` to `{origin, destination}`
across both directed records are refused with `CHECKPOINT_CORRUPT`, while the byte-identical payload
with an empty catalog encodes and decodes cleanly. Witness:
`defect_reciprocal_event_pair_is_unrepresentable_in_a_v2_checkpoint`.

The binding in constraint 3 appears to be the error: a reciprocal record is produced by the search
launched from the *other* endpoint, so its provenance legitimately names that endpoint's search. The
correct binding is to the discovering state recorded in the pair, not to each record's own origin.
This is a spec-interpretation decision that changes which checkpoints are accepted on the wire, so it
has been recorded rather than silently changed here.

**Cross-language update (2026-08-11).** The Python backend emits *and* restores
`tests/corpus/xlang/checkpoints/checkpoint-clean.json`: 28439 canonical bytes carrying six directed
records (three reciprocal pairs), four committed states, one committed KMC step and its historical
rate snapshot. This crate refuses the same bytes with `CHECKPOINT_CORRUPT` / `E2-EVENT-001`. A
conforming encoding therefore **exists**, and D-E2-02 is an implementation over-constraint rather
than a demonstrated specification impossibility. What Python does instead:

* it accepts one `search_id` per pair and requires the two records to carry an *identical*
  `discovery_provenance`, then checks that the id re-derives from **exactly one** of
  `{origin_state_id, destination_state_id}` — satisfying `E2-DISC-002` without demanding that each
  record's own origin be the discovering state;
* it does not require `reverse.saddle == event.saddle`. It requires the reciprocal saddle to be the
  forward saddle **mapped through** `validation.destination_match.atom_mapping` with the unstable
  direction negated, which is what a reciprocal record physically is.

Measured ladder against this crate, starting from the Python checkpoint and relaxing one constraint
at a time (each variant resealed):

| variant | this crate's verdict |
|---|---|
| Python bytes as written | `CHECKPOINT_CORRUPT` `E2-EVENT-001` |
| + `validation.method` set to `saddle.evidence_level` | `CHECKPOINT_CORRUPT` `E2-CKPT-007` |
| + reciprocal `saddle` made byte-identical (constraint 2 holds, 3 violated), substream map and `discovery_statistics.config_digest` rebuilt to this crate's expectations | `CHECKPOINT_CORRUPT` `E2-CKPT-007` |
| + provenance re-derived from each record's own origin (constraint 3 holds, 2 violated), same rebuilds | `CHECKPOINT_CORRUPT` `E2-CKPT-007` |
| control: same model, empty catalog, zero steps | accepted, and byte-identical to Python's answer |

The first row also exposes a fourth unspecified constraint: `validate_event_v2` requires
`validation.method == saddle.evidence_level`, while `E2-EVENT-005` lists `method` without
constraining it and `E2-EVENT-002` constrains only `saddle.evidence_level`. Python stores
`full-endpoint-relaxation/1`, which reads as the validation method rather than the saddle evidence
level. A fifth: `validate_checkpoint_v2` compares `discovery_statistics[*].config_digest` against the
digest of the `discovery` object alone, while Python stores the full model `config_digest`;
`E2-DISC-004` names the field `config_digest` and does not say which config.

Consequence: `E2-PAR-002` items 7, 9 and the exact-next-event resume of item 10 cannot execute. Every
`kmc-*` and `event-*` fixture is limited to its public-constructor half, `checkpoint-clean` is limited
to a zero-step payload, `checkpoint-recursive-corruption` cannot reach `E2-CKPT-007` recursion level 5,
and `checkpoint-exact-next-event-resume` cannot execute at all.

## Findings (recorded, not defects)

**`F-E2-01` — the identity certificate is ULP-sensitive to translation.** `E2-ID-002` states the
certificate is invariant to whole-cell translation. That holds in exact arithmetic. In binary64 a
translation whose coordinate differences are not exactly representable shifts the certificate rows
by 1-2 ULP (measured: 4.441e-16 A for a 3.3/4.4/5.5 A offset), and because `E2-CAN-005` hashes those
exact bytes the candidate and state IDs change. The drift is 4.4e-4 of the `PAR-004` geometry
tolerance of 1e-12 A, so the two geometries are the same state scientifically. `E2-ID-005` and
`STATE-006` are the intended mitigation. An exactly representable translation is invariant, so the
anchor/sort/closest-image construction itself is sound. Witness:
`identity_certificate_is_ulp_sensitive_to_translation`.

**`F-E2-02` — the Philox counter wraps in place on the failing increment. Severity: low, unreachable.**
`rng::increment_counter` mutates each lane before it discovers saturation, so a generator that
refuses the terminal block is left with `next_counter = [0,0,0,0]` while `consumed_blocks` is still
`0`, violating `E2-RNG-002`'s `next_counter = initial_counter + consumed_blocks` in the live object.
The broken state cannot be persisted: `rng::validate_state` recomputes the same relation and refuses
it, so the failure stays closed. Reaching it requires consuming 2^128 blocks. Witness:
`rng_counter_carry`.

**`F-E2-03` — requirement IDs in rejection contexts are coarse. Partly closed 2026-08-11.** Every
strict-parse rejection reported `context.requirement_id = "E2-JSON-001"`, including out-of-domain
integers, whose rule is `E2-JSON-002`. The integer case is now cited correctly by
`parity::parse_public` and agrees with Python. The remainder stands, and the cross-language run makes
it concrete: for an unknown key in a normative object this crate reports `E2-API-002` (the request
deserialiser is where `deny_unknown_fields` fires) while Python reports `E2-SCOPE-003`, which is the
requirement that actually governs unknown keys; for a checkpoint digest mismatch this crate reports
`E2-CKPT-007` and Python reports `E2-CKPT-003`; for the recursive-corruption fixture this crate
reports `E2-EVENT-001` and Python reports `E2-RATE-001`. Errata 2 does not prescribe which
requirement ID accompanies which failure, so only the unknown-key case is decidable, and there
`E2-SCOPE-003` is the governing requirement.

**`F-E2-04` — an out-of-domain integer in an internal record is refused as `NONFINITE_RESULT`.**
Setting a checkpoint resource counter to `u64::MAX` is refused by `checkpoint::canonical_json_bytes`
with `NONFINITE_RESULT` / `E2-CAN-005` rather than `INVALID_INPUT` / `E2-JSON-002`. This is an
internal encode path and the value is genuinely unrepresentable, so the encode correctly fails
closed; only the status token differs from what the portable-integer rule would suggest.

## Cross-language parity

`E2-PAR-005`: "Execution conformance additionally requires every fixture to pass in every
implementation." That was previously untested in either direction: each backend had only ever been
run against its own suite. It is now measured.

* Shared corpus: `tests/corpus/xlang/` — 46 raw request byte strings, 30 probe inputs and 9
  checkpoint cases, covering all 46 mandatory `E2-PAR-002` fixtures. Every rejection case is paired
  with an accepted baseline differing only in the property under test.
* Rust driver: `tests/xlang_emit.rs` (`cargo test --test xlang_emit`, output directory from
  `SPARK_XLANG_OUT`). Python driver and comparator:
  `../spark-atomistic/tests/xlang_harness.py`.
* Each backend emits 85 canonical artifacts; the comparator diffs them byte-for-byte and reports the
  first differing offset with both surrounding strings. Recorded run, 2026-08-11, re-run after
  `src/run.rs` was added so the corpus is answered by an adapter-backed backend: **52 cases
  identical, 33 divergent, 0 missing; 31 of 46 fixtures identical, 15 divergent.**
* Everything scientific agrees byte-for-byte: all state/candidate identities and geometry
  certificates (translated, periodic-image and same-species-permuted), all Philox blocks, uniforms,
  raw bit patterns, stream derivations, counter carries and substream digests, all search IDs and
  the deterministic discovery class sequence, all event/pair/reverse ID algebra, all
  `COMMON_PREFACTOR` rates and detailed-balance residuals, all rate-table snapshots including the
  lost-rate bound and `payload_sha256`, all KMC selections and residence times, the capability
  value, the schema/config/tolerance/identity digests, and all 29 status messages and exit codes.
* Every divergence is in the failure-response envelope or the checkpoint validator:
  `D-E2-05` (`context.component`), `D-E2-06` (`context.details`), `D-E2-07` (severity strings),
  `D-E2-02` (every trajectory-bearing checkpoint) and the coarse requirement IDs of `F-E2-03`.
  The two `run` cases no longer differ in status, exit code or requirement ID: with the adapter
  present both backends answer `INVALID_INPUT` / exit `64` / `E2-SCOPE-004` /
  `details.extension_field = "calculator_command"` and differ only in `context.component`
  (`adapter` in Python, `api` here). Errata 3 section 5 called the adapter surface unreachable; the
  reachable part of it is now reached, and what is left is the section-1 vocabulary item.
* The per-fixture verdicts live in `tests/corpus/e2_parity_manifest.json` under
  `cross_language_execution` and are checked against the fixture list by
  `parity_manifest_matches_this_suite`, which also refuses to let `e2_par_005_claimable` be `true`
  while any fixture diverges.

## Run adapter and cross-language run parity

`src/run.rs` supplies the `parity::RunAdapter` the crate previously lacked, so an `E2-API-003` run
request now executes: process-isolated JSON-lines calculator, backtracking steepest-descent
minimiser, minimum-mode-following (`directional-dimer`) saddle solver, endpoint validation,
directed reversible catalog with tolerance deduplication, `COMMON_PREFACTOR` rates, deterministic
serial KMC, canonical checkpoint write and resume, trajectory and summary artifacts.
`E2-SCOPE-004` is enforced structurally: the argv comes from `extension.calculator_command` and is
used only to spawn the child; every recorded identity is model-derived
(`calculator.model_name`/`model_version`) or a constant of the adapter.

### Analytic acceptance (`VAL-004`-shaped, one potential)

Calculator: the shared clean-room fixture `../spark-atomistic/corpus/mock_calculator.py`,
`V = sum_i (x_i^2-1)^2 + y_i^2 + z_i^2` eV/A, model digest
`7c799e3c0c25eb952d433430027d3d73de8d9f8f3d06064b6374f4b6eab4dd47`. Closed form: minima at
`x = +-1` with `E = 0`, first-order saddle at `x = 0` with `E = 1.0` eV, forward and reverse
barriers exactly 1.000 eV, Hessian eigenvalues `(12x^2-4, 2, 2)` at `(x,0,0)`.

Launched from `x0 = 0.30` with `force_tolerance = 1e-4` (`analytic_double_well_barriers_reproduce_the_closed_form_from_x0_030`):
forward barrier `1.000000000000` eV, reverse barrier `1.000000000000` eV, absolute saddle energy
`1.250000002418` eV against the exact `1.25` (the frozen marker atom contributes a constant
`0.25` eV), 34 dimer iterations, 85 calculator evaluations.

Accuracy is set by the force tolerance, not by the algorithm
(`analytic_saddle_accuracy_is_bounded_by_the_force_tolerance_not_by_the_algorithm`). Along the
reaction coordinate `V ~ E_s - 2x^2` with `|F_x| ~ 4|x|` gives `f^2/8`; the transverse coordinates
`V ~ y^2` with `|F_y| = 2|y|` give `f^2/4`, and the residual force may sit anywhere on that sphere,
so `f^2/4` is the attainable bound for this potential:

| `force_tolerance` (eV A^-1) | measured \|E_saddle - 1.25\| (eV) | `f^2/8` | `f^2/4` | evaluations |
|---|---|---|---|---:|
| 5e-2 | 1.8161099414504278e-4 | 3.125e-4 | 6.25e-4 | 29 |
| 1e-2 | 1.801240568211071e-5 | 1.25e-5 | 2.5e-5 | 45 |
| 1e-3 | 2.097489799623986e-7 | 1.25e-7 | 2.5e-7 | 65 |
| 1e-4 | 2.418317102126366e-9 | 1.25e-9 | 2.5e-9 | 85 |
| 1e-5 | 1.7843948540985366e-11 | 1.25e-11 | 2.5e-11 | 107 |

Launched from `x0 = 0.90`, with the identical model, the identical substream-seeded initial mode
and the identical tolerances -- only the launch coordinate differs -- the search returns
`SADDLE_NOT_FOUND` / `iteration-budget-exhausted` after 400 iterations and 1200 evaluations. At
`x = 0.9` the Hessian is `diag(5.72, 2, 2)`: every eigenvalue is positive and the lowest mode is
transverse, so minimum-mode following climbs the wrong mode. `SADDLE-002` permits any first-order
saddle solver and `E2-STATUS-002` classifies `SADDLE_NOT_FOUND` as a `candidate reject` that cannot
terminate a public operation, so this is a characterised parameter-domain limitation, not a
violation.

### End-to-end run

`run_operation_executes_the_whole_pipeline_and_reports_where_it_stops`: two atoms (one movable H,
one frozen He marker), 12-attempt budget, strict mode. Measured: 6 attempts, 2 successes, 1
duplicate, 4 `SADDLE_NOT_FOUND`, `stopping_state = CONVERGED_HEURISTIC` by the `DISC-004` redundant
-success rule; 2 committed minima; exactly 1 reciprocal pair (2 directed records) with forward and
reverse barriers `0.999999999160` and `1.000000000000` eV; 1 committed KMC step with
`k = 1.5875938077951148e-4` s^-1 and `dt = 2.0980472679934087e3` s; 2041 calculator evaluations.

The run then stops at `CHECKPOINT_CORRUPT` / `E2-CKPT-007` / exit `74`, which is `D-E2-02` reached
from the production path instead of from a hand-built fixture: the first checkpoint that carries a
committed event is refused by this crate's own validator. The paired accepted baseline is in the
same run: the pre-discovery checkpoint, whose only difference is an empty catalog, is written and
validated normally.

### Cross-language RUN comparison

One canonical model (`config_digest`
`sha256:fa6ec6f4d1d632a3771ab1233364952122552cb4d70286f6d413e958fe370180`) written once and read
verbatim by both backends, both driving the same `mock_calculator.py`.

| Artifact | Verdict | Python bytes | Rust bytes | First differing offset |
|---|---|---:|---:|---:|
| `checkpoint.json` | DIVERGENT | 16164 | 5438 | 89 (inside `catalog.digest`) |
| `trajectory.json` | MISSING in Python | - | 1214 | - |
| `summary.json` | MISSING in Python | - | 282 | - |
| public response | DIVERGENT | 383 | 282 | 46 (`context.component`) |

What agrees: the checkpoint envelope and payload key sets at every level (`E2-CKPT-002`..`006`
conformance on both sides), and all four `digests` -- `config`, `model`, `schema` and `tolerances`
are byte-identical. Everything derived from the MODEL agrees; everything derived from a SOLVER does
not.

**`D-E2-08` — no requirement pins the minimiser or the saddle solver, so a `run` artifact cannot be
byte-identical across implementations. Severity: high; specification-level.** `E2-PAR-003` requires
checkpoint bytes to be identical. `E2-ID-004` makes `state_id` a SHA-256 over the exact binary64
geometry certificate and energy. `RELAX-001..004` and `SADDLE-002` deliberately leave the algorithm
free ("The interface MUST permit minimum-mode following, dimer-family, ART-family, or another
first-order-saddle solver"). Measured on the identical model and calculator: the two minimisers
commit initial minima 2.661e-6 A apart with a 3.585e-10 eV energy difference (Python
`x0 = [0.9999999999999983, 3.234618539845976e-5, -3.234618539845976e-5]`, 62 evaluations, 61 steps;
Rust `x0 = [1.0000000434990044, 3.500739258819418e-5, -3.500739258819418e-5]`, 17 evaluations, 16
steps). Both satisfy `STATE-004` at `force_tolerance = 1e-4`. Different geometry, therefore
different certificate, therefore different `candidate_identity`, `state_id`, `catalog.digest` and
`payload_sha256`. This is not a defect in either backend: the requirement pair is unsatisfiable
unless the solvers are pinned or `E2-PAR-003` names its run-artifact exclusions. It also makes
`PAR-004`'s "geometry absolute 1e-12 A" unreachable for solver output by 2.7e6 times.

**`D-E2-09` — the trajectory and summary artifacts have no schema in any requirement. Severity:
medium; specification-level.** `E2-SCHEMA-010` names `trajectory_path` and `summary_path` and fixes
nothing about their content, and no other requirement mentions them. This crate emits canonical
`{"records":[<E2-KMC-004 step records>],"schema":"spark-atomistic-trajectory/2"}` and
`{"schema":"spark-atomistic-summary/2","status":<token>,"value":{<E2-API-008 keys>}}`, which are the
only spec-anchored records available. The Python backend writes neither file when the run pauses at
`DISCOVERY_INCOMPLETE`. Even the existence of the artifacts is therefore unpinned, before any
question of their bytes.

**`D-E2-10` — checkpoint bytes are not reproducible even within ONE backend, because
`resources.wall_elapsed_s` is a wall-clock measurement. Severity: high; specification-level.**
`E2-PAR-003` lists "checkpoint bytes" among the artifacts that must be byte-identical, and
`E2-CKPT-008` requires replayed IDs, uniforms, rate arrays, selections and counters to be
byte-identical. `E2-CKPT-005` requires `resources` to contain `wall_elapsed_s`. Measured: the same
model run twice through this backend, same directory, same calculator, produced checkpoints that are
both 5439 bytes and differ at byte offset 4909 -- `"wall_elapsed_s":0.353757169` against
`"wall_elapsed_s":0.355789421` -- which changes `payload_sha256` from
`sha256:296b164bf5390f9c9f0bfa9016e48d14854016d4d1529aecea3b5b3cf20eb7f4` to
`sha256:17a77469f3add5c1da0d0b0a41d9f636cd120278a8d047b827db3de330396f8a`. The trajectory, summary
and public response of the same two runs are bit-identical, so the nondeterminism is confined to
that one field. No wording in Errata 2 exempts it. A byte-identity requirement and a mandatory
wall-clock field cannot both hold; the erratum must either exclude `wall_elapsed_s` from the
identity set, quantise it, or drop it.

**`F-E2-05` — the two backends disagree about `flags.incomplete_catalog` in strict mode, on an
unpinned relation.** This crate's `validate_checkpoint_v2` requires
`flags.incomplete_catalog == any(discovery_statistics[*].permanently_incomplete_catalog)`. The
Python checkpoint measured here carries `flags.incomplete_catalog = true` with
`permanently_incomplete_catalog = false` (strict mode, `DISCOVERY_INCOMPLETE`), so this crate would
refuse it for a second, independent reason beyond `D-E2-02`. `E2-CKPT-004` fixes the key set of
`flags` and `E2-DISC-004` the key set of the statistics, but neither states the relation between the
two fields; `DISC-006` distinguishes a strict-mode pause from an exploratory permanent tag, which is
the reading this crate implements.

**`F-E2-06` — the two backends disagree about `resources.output_bytes` and
`resources.resident_memory_bytes`, also unpinned.** Python drives `output_bytes` to the
self-referential fixed point of the checkpoint it is writing (16164, its own size) and records a
measured resident set (22540288). This crate records the pre-write value (0) and does not instrument
resident memory at all (0). `E2-CKPT-005` fixes the key set and `RES-001` requires the limits to
exist; neither states what the counters must contain at the moment of writing. The missing
resident-memory instrumentation is a real gap in this crate: `RES-001`'s resident-memory limit is
validated but never enforced during a run.

**`F-E2-07` — the initial discovery perturbation is unpinned, and the two backends therefore build
different catalogs from the same model.** `DISC-001` requires "a configured mixture of global,
local, and optional targeted initial perturbations" and `E2-SCHEMA-007` configures the class names,
kinds and probabilities -- but no requirement fixes the perturbation amplitude or geometry. This
crate uses `saddle_search.endpoint_displacement` as the amplitude, which is declared in
`src/run.rs`. Measured on the identical model: this crate finds the saddle in 2 of 6 attempts and
converges heuristically; the Python backend returns `SADDLE_NOT_FOUND` on all 12 attempts and
finishes `DISCOVERY_INCOMPLETE` with an empty catalog. Neither backend is wrong under `DISC-001`.
The consequence is stronger than a byte mismatch: the same model produces different SCIENCE in the
two implementations, and no requirement forbids it.

## Claims deliberately not made

P0 conformance and publishability are not claimed. `src/run.rs` now supplies the documented
calculator adapter, minimiser and saddle solver that `ACCEPT-002` requires, but `VAL-001..014` is
still unrun as a matrix: only `VAL-004`-shaped analytic agreement has been measured, on one
one-degree-of-freedom potential. `PAR-004`/`PAR-005` numeric-tolerance agreement across
implementations is measured for the first time below and is **not met** for a run: the two
minimisers stop 2.661e-6 A apart, which is 2.7e6 times the `PAR-004` geometry tolerance of 1e-12 A,
because no requirement pins a minimiser. No P1 missing-rate, harmonic-prefactor, or local-reuse
efficacy claim is made. No performance, throughput, or timing claim is made, and no benchmark was
run. `E2-PAR-005` execution conformance is explicitly **false**: 15 of 46 mandatory fixtures do not
agree between the two implementations, and `D-E2-02`, `D-E2-03`, `D-E2-04`, `D-E2-05`, `D-E2-06`,
`D-E2-07`, `D-E2-08`, `D-E2-09` and `D-E2-10` all stand. `D-E2-01` is fixed, which is a precondition for the
claim, not the claim itself.
