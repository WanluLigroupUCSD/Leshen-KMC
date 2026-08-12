# spark-atomistic-rs

Independent clean-room Rust core for fixed-composition, serial, off-lattice/on-the-fly KMC. Normative wire revision: Errata 3; IR `spark-atomistic-model/1`; backend `rust`.

Status: `validated=false`, `production=false`, `release=false`. Trajectory correctness remains conditional on Markovian state dynamics, valid rates, and sufficient exit discovery. A finite saddle budget never proves catalog completeness.

## Build and test status

Current gate recorded 2026-08-12:

| Command | Result |
|---|---|
| `cargo check --all-targets` | 0 errors, 0 warnings |
| `cargo test --all-targets` | 69 tests, 69 passed, 0 failed, 0 ignored |
| Cross-language core/adapter | 85/85 cases and 46/46 fixtures pass; core 83, adapter 2 |

Re-recorded 2026-08-11 after the `D-E2-01` fix, the cross-language emitter, and the addition of the
run adapter (`src/run.rs`) with its execution tests. `tests/run_e2e.rs` spawns the shared Python
calculator fixture once per energy/force evaluation, so it dominates the suite wall time; no timing
figure is reported or claimed.

The crate had never been compiled before this session. Reaching a clean `cargo check` required three source repairs, none of which changed observable behavior: `src/resource.rs` used an illegal comma-separated `let` binding, `src/status.rs` needed an explicit type annotation on a `next_value` call (E0282), and `src/catalog.rs` carried two unused imports. Formatting is left alone on purpose: the source is a deliberate compressed single-line style and reformatting would invalidate every provenance hash.

A green test run is **not** a conformance claim. See `tests/corpus/e2_parity_manifest.json` for the per-fixture execution record and `REQUIREMENTS.md` for the defects and blockers found while executing it.

## What has and has not been exercised

- **Executed**: the whole `parity::dispatch_json` wire surface (`capabilities`, `validate`, adapter-backed `run`), portable-JSON and canonical-encoding rules, schema/config/tolerance/identity digests, geometry certificates and state IDs, event/pair ID algebra, `COMMON_PREFACTOR` rates and detailed balance, deterministic discovery class selection and search IDs, Errata-1 Philox with revision-2 stream derivation, rate-table snapshots, and revision-2 checkpoint encode/decode for a zero-step payload.
- **Executed cross-language (2026-08-12)**: the shared corpus `tests/corpus/xlang/` is consumed by both backends. All 85 cases and 46 fixtures pass the Errata 3 tiered gate: core 83, adapter 2.
- **Executed end to end (2026-08-11)**: a complete `run` operation. `src/run.rs` binds a process-isolated JSON-lines calculator from `extension.calculator_command`, relaxes, follows the minimum mode to a first-order saddle, validates both downhill endpoints, builds a directed reversible catalog, runs serial KMC and writes trajectory/summary/checkpoint artifacts. On the shared analytic double well it reproduces the closed-form saddle and both 1.000 eV barriers; see `REQUIREMENTS.md`.
- **Not executed**: the `callbacks` receipt layer (`src/run.rs` produces the Errata-2 wire records directly); the `VAL-001..020` scientific matrix as a matrix; resident-memory enforcement; benchmarks or any timing measurement. Nothing under `SPARK/examples/` was run.
- **Unblocked by normative Errata 3**: a checkpoint containing a catalog event encodes and restores, `run` returns `OK`, and resume commits the same next event as an uninterrupted run. `CheckpointPolicy::PreErrata3` remains only as a compatibility reader/test path.

## Scope

The public parity surface is `parity::dispatch_json`: exact `capabilities`, `validate`, and adapter-backed `run` requests/responses. It enforces the frozen ten-object model, extension/digest separation, portable JSON, RFC 8785-style canonical bytes, exact schema/config/tolerance/identity digests, revision-2 state/event/search IDs, Errata-1 Philox with revision-2 stream derivation, historical rate snapshots, revision-2 checkpoint envelopes, and uniformly disabled basin acceleration.

P1 finite-basin sampling is hard-disabled because no internal catalog/discovery completeness validator exists. No exact-basin claim or enabling receipt exists; requests fail with `BASIN_DISABLED`. Missing-rate estimation, harmonic prefactors, variable composition, GPU, Time Warp, tau-leap, and parallel-time KMC are disabled/not implemented.

The run-request field `allow_unvalidated` is mandatory and explicit, and it lives on the request root, not inside `saddle_search` (`E2-API-003`, `E2-API-004`). Only the JSON boolean `true` authorizes execution of an explicitly unvalidated reference; `false`, numeric `1`, a string, `null`, or a missing field returns `INVALID_INPUT` before any calculator, output, or checkpoint work. It does not authorize unvalidated events: unvalidated candidates cannot construct the crate-private catalog validation receipt and are never selectable.

The package is Rust-only. No PyO3 module, Python extension, automatic backend, or dynamic-language adapter is declared.

## Safety and limitations

- `src/run.rs` bundles a calculator transport, a minimiser and a dimer-family saddle solver. There is still no CLI and no Python orchestration adapter. The solver constants (perturbation amplitude, dimer separation, translation and rotation steps) are declared constants of that file because no requirement pins them; changing them changes the discovered catalog, and the two backends already differ on exactly this point.
- Callback implementations must honor deadline/cancellation and return full termination/evaluation accounting.
- Retry implementation is absent; `resources.retry_count` must be exactly `0`.
- Exact identity search fails closed on pathological cells or closest-vector/permutation budget exhaustion.
- Rate overflow/underflow fails closed; no undeclared cutoff or clipping exists.
- Revision-2 restore is exposed through `parity::decode_checkpoint_v2`; old internal records are not parity artifacts and are incompatible on the public wire.
- `run` requires `allow_unvalidated:true` and a `RunAdapter`; `parity::dispatch_json` still takes the adapter as a parameter, and `run::ProcessRunAdapter` is the one this crate provides. Adapter `extension` never enters scientific digests or records: its only consumer is the child-process spawn, and every recorded identity is model-derived.
- Mandatory fixtures are indexed by `tests/corpus/e2_parity_manifest.json` and tiered by `tests/corpus/xlang/manifest.json`. The current gate passes all 46 fixtures.
- A separate cross-language RUN comparison over one shared model is recorded in `REQUIREMENTS.md`: every model-derived digest (`config`, `model`, `schema`, `tolerances`) is byte-identical and the checkpoint key sets match at every level, but every emitted run artifact diverges, because no requirement pins the minimiser, the saddle solver, or the discovery perturbation (`D-E2-08`, `D-E2-09`, `F-E2-07`), and because `resources.wall_elapsed_s` made checkpoint bytes irreproducible even within one backend (`D-E2-10`, addressed by the `E3-PAR-002` candidate below). A refreshed hash manifest is an integrity record, never evidence of correctness.
- Canonical binary64 output is round-trippable under `E3-CAN-001`, including finite values with magnitude at least `2^53` and below `1e21`.
- Checkpoint byte accounting reserves a conservative fixed 24-byte slot for each canonical `wall_elapsed_s` token. This preserves the real output limit and prevents wall-clock formatting from leaking into the hashed cumulative `resources.output_bytes` counter (`D-126`).
- The candidate/state identity is byte-exact over binary64 coordinates. A whole-cell translation that is not exactly representable perturbs the certificate by 1-2 ULP and therefore changes the state ID, even though the geometries agree to 4.4e-16 A, far inside the `PAR-004` tolerance of 1e-12 A. `E2-ID-005` and `STATE-006` are the intended mitigation: candidate IDs are hints and geometry verification decides equivalence.
- Package license: [PolyForm Strict License 1.0.0](LICENSE), limited by [LICENSE_SCOPE.md](LICENSE_SCOPE.md). It permits noncommercial use but not distribution, modification, or derivative works. It is source-available, not open source. Dependency license/patent review remains a release gate.

## Errata 3 — normative, adopted by D-127

`OFFLATTICE_OTF_KMC_SPEC_V1_ERRATA_3.md` is normative. D-127 adopts the implementation below.
Pre-Errata-3 read/test paths remain for old checkpoints; writers use Errata 3.

| Item | Change | Pre-Errata-3 path |
|---|---|---|
| `E3-EVENT-001` options 1 + 4 | C3 dropped (the search ID must re-derive from exactly one of `{origin, destination}` and both records of a pair carry an identical `discovery_provenance`); C2 relaxed to mapped equality (`validation.destination_match.atom_mapping`, `unstable_direction` negated); `validation.method` no longer forced to equal `saddle.evidence_level`; `discovery_statistics[*].config_digest` bound to `E2-CAN-007`'s model digest | `parity::CheckpointPolicy::PreErrata3` |
| `E3-PAR-002` | `resources.wall_elapsed_s` excluded from the HASHED payload; the field stays in the record; both seals are accepted on read and a re-encode reproduces the seal it read | `parity::PayloadHashPolicy::PreErrata3Full` |
| `E3-STATUS-001` / `E3-STATUS-002` option (a) | `context.component` names the failing subsystem (`api`, `json`, `schema`, `adapter`, `output`) and `context.details` carries `path`/`missing`/`unknown`/`required_keys`/`extension_field`, threaded on `ApiFailure` rather than mapped at the boundary | delete the `component` field and the `.in_component(..)` / `.with_detail(..)` call sites |

Two further items were decidable from the FROZEN text and are not Errata 3 changes: a barrier that
does not reproduce from its own endpoints is now cited to `E2-RATE-001` (the requirement that states
`b_f = E_s - E_i`) rather than to the `E2-EVENT-001` field table, and `E2-EVENT-004`'s
`rng_substream_digest` is now taken from the substream state the checkpoint actually stores, without
which `E2-CKPT-007`(8) can never verify the map.

`tests/errata3_candidate.rs` retains its historical filename and is the accept/reject evidence: the 28439-byte Python checkpoint is
accepted, the same bytes are still refused under `PreErrata3`, and three corruptions are refused
(an unmapped reciprocal saddle, a provenance deriving from neither endpoint, a forged payload hash),
each against a baseline that differs only in the property under test.

## Dead gates repaired

- **Detailed balance (`src/rate.rs`)**: under `COMMON_PREFACTOR` the `RATE-005` residual is an
  algebraic identity of `RATE-002` plus `E2-RATE-001`, so the comparison verifies nothing. Neither
  escape is available: an independently obtained reverse barrier is forbidden by `E2-RATE-001`
  ("Raw barriers are exactly `b_f=E_s-E_i` and `b_r=E_s-E_j`"), and deleting the recorded residual
  is forbidden by `E2-EVENT-003`, `E2-CKPT-007`(5) and `E2-PAR-002`(7). The comparison is therefore
  kept exactly as `RATE-005` words it and documented at the site as a ROUNDING guard, with the
  measurement in `conformance_spec::detailed_balance_gate_measures_rounding_not_physics`: worst
  `|residual| = 7.275957614183426e-12` over 200,000 randomized inputs against the `1e-8` default,
  and every one of the 90059 nonzero residuals rejected once the tolerance is pushed to `0`.
- **Saddle order (`src/run.rs`)**: `SADDLE-005` allows "nonnegative sampled orthogonal curvatures",
  but `E2-EVENT-005` also mandates `unstable_mode_count = 1`, which a fixed sample cannot establish:
  `v^T H v` on any chosen basis is a positively weighted average of the spectrum in that direction.
  This crate sampled DETERMINISTIC Gram-Schmidt Cartesian probes, not random ones, and had the same
  defect. The sample is replaced by a MINIMISATION of the Rayleigh quotient over the orthogonal
  complement, seeded from exactly those probes, so the reported value can only fall. The result is
  an upper bound on the smallest restricted eigenvalue: it can refute index 1 and cannot prove it,
  which is what the `DIRECTIONAL` evidence level of `E2-EVENT-002` means. Witness in
  `run_e2e::saddle_order_evidence_is_minimised_over_the_orthogonal_complement_not_sampled`: a
  transverse block with Cartesian diagonal `(2, 2)` and eigenvalues `5, -1` was accepted by the
  fixed sample and is refuted by the minimisation.

See `REQUIREMENTS.md`, `tests/TEST_SPEC.md`, `tests/ERRATA2_PARITY_SPEC.md`, and `PROVENANCE.json`.
