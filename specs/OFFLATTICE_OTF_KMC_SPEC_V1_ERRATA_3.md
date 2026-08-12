# Off-lattice/on-the-fly KMC specification v1 — Errata 3

Date: 2026-08-12  
Status: normative; adopted by owner-delegated decision `D-127`  
Applies to: `OFFLATTICE_OTF_KMC_SPEC_V1.md`, Errata 1, Errata 2  
Base self-excluding SHA-256: `8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84`  
Errata 1 SHA-256: `52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40`  
Errata 2 self-excluding SHA-256: `eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995`

## 1. Scope and precedence

This erratum resolves the six contradictions measured in the 2026-08-11 cross-language run. It overrides conflicting Errata 2 interpretations of `STATUS`, reciprocal checkpoint events, canonical large-real encoding, and `PAR`. Pre-Errata-3 checkpoints remain readable; new checkpoints use this erratum.

## 2. Status context

`E3-STATUS-001` `context.component` names the refusing subsystem. The exhaustive vocabulary is `adapter`, `api`, `calculator`, `catalog`, `checkpoint`, `discovery`, `geometry`, `json`, `kinetics`, `output`, `rates`, `relaxation`, `resources`, `rng`, `saddle`, `schema`, and `state`.

`E3-STATUS-002` `context.details` is normative for every shared-corpus response. Its exact keys and values are the shared corpus outputs. Adding a detail key requires a corpus revision; implementations MUST NOT silently omit a defined key.

## 3. Reciprocal checkpoint events

`E3-EVENT-001` One physical saddle produces one reciprocal pair. Both directed records carry identical discovery provenance. Its search ID MUST re-derive from exactly one of the two endpoint state IDs: the state from which the saddle was discovered. The reciprocal saddle MUST equal the forward saddle after applying `validation.destination_match.atom_mapping` to positions and negating `unstable_direction`. Byte equality before this mapping is not required. `validation.method` and `saddle.evidence_level` are independent fields. `discovery_statistics[*].config_digest` binds the full model config digest.

## 4. Closed canonical real domain

`E3-CAN-001` A binary64 real with absolute value at least `2^53` MUST use shortest round-tripping exponent form, including values below `1e21`. Integer tokens remain restricted to `[-9007199254740991,9007199254740991]`. Thus every canonical real emitted by `E2-CAN-004` is accepted as a real by `E2-JSON-002` and re-encodes byte-identically.

## 5. Parity tiers and run artifacts

`E3-PAR-001` The mandatory corpus has two tiers:

- `core`: capability, validation, error/status, model-derived digests and IDs, RNG, event/catalog/rate/KMC records, and checkpoint validation. Every core artifact is byte-identical.
- `adapter`: executable, IPC, process isolation, and run transport governed by `E2-SCOPE-004`. Adapter cases are conditional on equivalent adapter bindings and do not block core conformance.

The corpus manifest labels every case. A conformance report MUST give counts per tier; a missing or divergent core case fails conformance.

`E3-PAR-002` `resources.wall_elapsed_s` remains mandatory in the checkpoint record but is excluded from `payload_sha256`. No other field is excluded. Readers accept both the Errata-2 full-payload seal and the Errata-3 environment-excluded seal; writers use the latter. Run artifacts are compared by declared field/tolerance rules, not blanket byte identity. Exact-coordinate state IDs and discovery outcomes are not cross-implementation parity claims unless the minimizer and discovery perturbation are also frozen.

## 6. Status boundary

This erratum makes the core parity contract attainable. It does not establish physical accuracy, catalytic validity, production readiness, release rights, or equivalence of discovery chemistry. Backends remain `implemented_unvalidated` until their separate validation and release gates pass.

## Errata 3 digest

Hash rule: SHA-256 of the exact UTF-8 bytes before the line beginning `## Errata 3 digest`; the separator blank line is included and this digest section is excluded.

Errata 3 SHA-256: `eba384af3694c5f3997caf28829e56d188ef9929f29d99a2116520f0067d8a96`
