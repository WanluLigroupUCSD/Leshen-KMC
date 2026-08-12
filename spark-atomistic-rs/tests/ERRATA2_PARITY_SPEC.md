# Errata 2 parity fixture specification

Status: executed by `errata2_parity.rs`. No fixture is waived. `e2_parity_manifest.json` enumerates
every `E2-PAR-002` fixture and now records a per-fixture execution status, checked against the suite
by `parity_manifest_matches_this_suite`. Result in this backend's own suite: 37 `PASS`, 8 `PARTIAL`,
0 `FAIL_DEFECT`, 1 `BLOCKED` of 46. Cross-language against the Python backend over the same corpus:
30 `PASS`, 16 `DIVERGENT` of 46, recorded in the manifest's `cross_language_execution` block.

Concrete shared inputs cover capability, minimal model, RFC 8785 boundaries, portable integers, the
duplicate-key document, and the Errata-1 midpoint records. Remaining scientific fixtures are generated
through the public revision-2 constructors in `parity`: `state_ids`, `event_ids`, `search_id`,
`geometry_certificate`, `fixed_contract_digest`, `constraint_digest`, `choose_discovery_class`,
`ordered_candidate_ids`, `make_rate_snapshot`, `CheckpointEnvelopeV2::new`, `encode_checkpoint_v2`,
and `decode_checkpoint_v2`.

## Assertion discipline

- Every assertion names the Errata 2 requirement that decides it.
- Every fixture asserted to be REJECTED is paired with an ACCEPTED baseline that differs only in the
  property under test, so a rejection can never be credited to an unrelated malformation.
- Golden values that Errata 2 states normatively — the schema digest `sha256:583d580d…`, and the two
  boundary-record hashes `6ce1fb52…` and `a15157604d…` — are recomputed from this crate's own
  canonical encoder rather than copied, and the trajectory/saddle stream key and counter words were
  derived independently from the `E2-RNG-004` byte recipe before being written into the suite.
- Two tests are defect witnesses. They are named `defect_*` and assert observed, nonconforming
  behavior so that the suite stays green without hiding the defect. `../REQUIREMENTS.md` carries the
  analysis; the manifest carries the per-fixture verdict. The `D-E2-01` witness was retired on
  2026-08-11 when the defect was fixed and became the conformance test `json_integer_out_of_domain`;
  a new witness `defect_e2_can_004_canonical_binary64_does_not_survive_this_crates_own_parser`
  records `D-E2-03`.
- The cross-language half lives in `xlang_emit.rs` and `corpus/xlang/`. Its per-case artifacts are
  compared byte-for-byte against the Python backend by
  `../../spark-atomistic/tests/xlang_harness.py compare`; divergences are reported with the exact
  first differing byte offset and both surrounding strings, never normalised away.

## Acceptance

Acceptance requires canonical byte equality for every public response, digest, ID, RNG state, event,
rate snapshot, step, catalog, and checkpoint. Scientific values use `PAR-004`; time uses relative
`5e-15` plus absolute `1e-18 s`.

Execution readiness is **false**. `E2-PAR-005` requires zero Critical, High, or Medium observable
mismatches plus every mandatory fixture passing in every implementation. `D-E2-01` (high) is fixed,
but `D-E2-02` (critical/blocking), `D-E2-03` (high), `D-E2-04` (medium, both backends), `D-E2-05`
and `D-E2-06` (medium, both backends) and `D-E2-07` (medium) stand, nine fixtures are not fully
passing in this backend's own suite, and the cross-language run shows 16 of 46 mandatory fixtures
diverging between the two implementations. No parity claim exists.
