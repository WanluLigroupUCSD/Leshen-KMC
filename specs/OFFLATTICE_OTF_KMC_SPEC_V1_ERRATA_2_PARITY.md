# Off-lattice/on-the-fly KMC specification v1 — Errata 2: parity wire contract

Date: 2026-08-09  
Status: normative  
Applies to: `OFFLATTICE_OTF_KMC_SPEC_V1.md` and Errata 1  
Base self-excluding SHA-256: `8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84`  
Errata 1 SHA-256: `52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40`

## 1. Scope and precedence

This erratum freezes observable JSON, mathematical, identity, RNG, event, discovery, and checkpoint behavior for `spark-atomistic-model/1`. It defines no source layout, implementation symbol, dependency, solver algorithm, or transport process.

`E2-SCOPE-001` This erratum overrides conflicting v1 or Errata 1 interpretations of `DET-*`, `CKPT-*`, `IO-*`, `ERR-*`, `EVENT-*`, `DISC-*`, `BASIN-*`, and `PAR-*`. Errata 1 remains normative for Philox word pairing and midpoint uniforms.

`E2-SCOPE-002` The frozen target is the smallest common P0 subset: fixed composition, fixed cell, validated directed reversible events, `COMMON_PREFACTOR`, deterministic serial KMC, deterministic discovery scheduling, canonical checkpoint/resume, and basin acceleration disabled.

`E2-SCOPE-003` Unknown keys in every normative object return `INVALID_INPUT`. Only the root model `metadata` object and the run-request `extension` object accept free-form keys. Their values still obey the portable JSON domain.

`E2-SCOPE-004` Adapter-specific transport, executable, IPC, environment, timeout-enforcement mechanism, and process-isolation settings belong only under run-request `extension`. They MUST NOT enter the model, capability value, config digest, schema digest, state/event/search ID, catalog digest, checkpoint payload, trajectory record, or scientific result.

## 2. Common public JSON surface

All operations accept one UTF-8 JSON object and return one UTF-8 canonical JSON object. There are exactly three operations.

### 2.1 Requests

`E2-API-001` Capability request:

```json
{"operation":"capabilities"}
```

`E2-API-002` Validation request:

```json
{"model":{},"operation":"validate"}
```

Here `model` is replaced by one complete object conforming to §3.

`E2-API-003` Run request:

```json
{"allow_unvalidated":true,"extension":{},"model":{},"operation":"run"}
```

Here `model` is replaced by one complete object conforming to §3. `extension` is an adapter-only object governed by `E2-SCOPE-004`.

`E2-API-004` `allow_unvalidated` is required only for `run`, MUST be the JSON boolean `true`, and authorizes execution of an explicitly unvalidated reference. Missing, `false`, numeric `1`, or any other value returns `INVALID_INPUT` before calculator, output, or checkpoint work. It does not authorize unvalidated events: only fully validated events are selectable.

`E2-API-005` Validation and capability operations perform no scientific work, create no output, and do not consume RNG.

### 2.2 Public response

`E2-API-006` Every public response has exactly these keys:

```json
{"causal_status":null,"context":{"component":"api","details":{},"requirement_id":"E2-API-001","retryable":false,"search_or_event_id":null,"state_id":null},"exit_code":0,"message":"transaction committed","operation":"capabilities","severity":"success","status":"OK","value":{}}
```

`operation` equals the request operation. `value` is `null` on failure. `causal_status` is `null` or one exact status token. A nonterminal status never appears as the final public response.

`E2-API-007` Successful capability `value` is exactly:

```json
{"api":"spark-atomistic-json/1","basin_acceleration":"disabled","conformance":"unvalidated","features":{"common_prefactor":true,"deterministic_checkpoint":true,"fixed_composition_off_lattice":true,"harmonic_tst":false,"local_environment_generic_reuse":false,"serial_kmc":true,"variable_composition":false},"ir":"spark-atomistic-model/1","operations":["capabilities","validate","run"],"production":false,"release":false,"validated":false}
```

`E2-API-008` Successful validation `value` has exactly `config_digest`, `ir`, and `schema_digest`. Successful run `value` has exactly `checkpoint_sequence`, `current_state_id`, `incomplete_catalog`, `simulation_time_s`, and `step_index`.

## 3. Exact `spark-atomistic-model/1` schema

`E2-SCHEMA-001` The root has exactly ten required keys and one optional key:

| Key | Type | Rule |
|---|---|---|
| `schema` | object | Exact fields below. |
| `system` | object | Exact fields below. |
| `calculator` | object | Exact fields below. |
| `relaxation` | object | Exact fields below. |
| `saddle_search` | object | Exact fields below. |
| `discovery` | object | Exact fields below. |
| `kinetics` | object | Exact fields below. |
| `resources` | object | Exact fields below. |
| `output` | object | Exact fields below. |
| `basin` | object | Must declare disabled. |
| `metadata` | object, optional | Free-form, nonbehavioral, excluded from every digest. |

`E2-SCHEMA-002` `schema` has exactly `id`; its value is `spark-atomistic-model/1`.

`E2-SCHEMA-003` `system` has exactly:

| Field | Type and constraint |
|---|---|
| `atom_ids` | Nonempty array of unique, nonempty strings. |
| `species` | Nonempty string array with length `N`. |
| `positions` | `N x 3` finite binary64 numbers, Å. |
| `cell` | `3 x 3` finite binary64 numbers, Å; nonsingular. |
| `pbc` | Three booleans. |
| `movable` | `N` booleans; at least one `true`. A value applies to all three Cartesian components of that atom. |
| `constraints` | Exactly `{"kind":"fixed-mask"}`. |
| `charge` | Finite binary64 number, elementary-charge units. |
| `spin` | Finite binary64 number, declared calculator convention. |
| `calculator_model_digest` | Nonempty string equal to `calculator.model_digest`. |

`E2-SCHEMA-004` `calculator` has exactly `deterministic`, `model_digest`, `model_name`, and `model_version`. `deterministic` is boolean; the strings are nonempty. Adapter configuration is forbidden here.

`E2-SCHEMA-005` `relaxation` has exactly `force_tolerance`, `max_evaluations`, and `max_steps`. Tolerance is finite and positive; limits are portable integers at least `1`.

`E2-SCHEMA-006` `saddle_search` has exactly `curvature_tolerance`, `endpoint_displacement`, `force_tolerance`, `max_iterations`, `method`, and `orthogonal_directions`. `method` is `directional-dimer`; numeric tolerances/displacement are finite and positive; integer limits are at least `1`.

`E2-SCHEMA-007` `discovery` has exactly:

| Field | Type and constraint |
|---|---|
| `mode` | `strict` or `exploratory`. |
| `classes` | Nonempty ordered array of class objects. |
| `minimum_successful` | Portable integer at least `1`. |
| `consecutive_redundant` | Portable integer at least `1`. |
| `maximum_attempts` | Portable integer at least `1`. |
| `maximum_evaluations` | Portable integer at least `1`. |
| `relevance_rate_min` | Finite number at least `0`, s⁻¹; `0` makes every positive rate relevant. |
| `alpha` | `null` or finite number in `(0,1]`. |
| `alpha_calibration` | `null` iff `alpha` is `null`; otherwise exact object `{id,source,version}` with nonempty strings. |

Each class object has exactly `kind`, `name`, and `probability`. `kind` is `global`, `local`, or `targeted`; names are nonempty and unique; probabilities are finite and positive, sum to `1` within absolute `1e-12`, and at least one non-targeted class has nonzero mass. Array order is behavioral.

`E2-SCHEMA-008` `kinetics` has exactly:

| Field | Type and constraint |
|---|---|
| `temperature` | Finite positive number, K. |
| `rate_model` | Exact string `COMMON_PREFACTOR`. |
| `prefactor` | Finite positive number, s⁻¹. |
| `barrier_tolerance` | Finite number at least `0`, eV. |
| `detailed_balance_tolerance` | Finite positive number; default `1e-8`. |
| `log_rate_cutoff` | Finite binary64 natural-log rate. |
| `absorbing_ok` | Boolean. |
| `maximum_steps` | Portable integer at least `1`. |
| `run_seed` | Portable integer in `[0,9007199254740991]`. |
| `state_rms_tolerance` | Finite positive number, Å; default `1e-3`. |
| `state_max_tolerance` | Finite positive number, Å; default `5e-3`. |
| `state_energy_tolerance_per_atom` | Finite positive number, eV atom⁻¹; default `1e-6`. |
| `saddle_rms_tolerance` | Finite positive number, Å; default `1e-3`. |
| `saddle_max_tolerance` | Finite positive number, Å; default `5e-3`. |
| `saddle_energy_tolerance` | Finite positive number, eV; default `1e-5`. |

`E2-SCHEMA-009` `resources` has exactly `callback_timeout_s`, `catalog_events`, `evaluations_per_relaxation`, `evaluations_per_saddle_attempt`, `output_bytes`, `resident_memory_bytes`, `retry_backoff_s`, `retry_count`, `saddle_attempts_per_state`, `total_calculator_evaluations`, and `wall_time_s`. Time values are finite; positive except `retry_backoff_s=0`. Count/byte values are portable integers at least `1`, except `retry_count=0`.

`E2-SCHEMA-010` `output` has exactly `checkpoint_every_steps`, `checkpoint_path`, `checkpoint_wall_time_s`, `overwrite`, `resume`, `summary_path`, and `trajectory_path`. Paths are nonempty strings and resolve relative to the request source when a source path exists; otherwise they must be absolute. The three resolved paths must be distinct. Frequency values are positive; booleans are exact JSON booleans.

`E2-SCHEMA-011` `basin` has exactly one boolean field:

```json
{"enabled":false}
```

`true` is schema-valid but cannot enable an accelerator under this revision. Validation still returns `OK`; a run records `BASIN_DISABLED` and continues serially.

`E2-SCHEMA-012` Cross-field limits require `relaxation.max_evaluations <= resources.evaluations_per_relaxation`, `discovery.maximum_attempts <= resources.saddle_attempts_per_state`, and `discovery.maximum_evaluations <= resources.total_calculator_evaluations`.

## 4. Portable JSON and canonical encoding

`E2-JSON-001` Input is UTF-8 without BOM. Strings contain Unicode scalar values only; lone UTF-16 surrogate code points are forbidden. Object keys are strings. Duplicate keys at any nesting depth return `INVALID_INPUT`.

`E2-JSON-002` Every integer is in `[-9007199254740991,9007199254740991]`. A syntactically valid integer outside this domain returns `INVALID_INPUT` before schema validation.

`E2-JSON-003` Every noninteger number is interpreted as IEEE-754 binary64 and must be finite. Overflow to infinity and the nonstandard tokens `NaN`, `Infinity`, and `-Infinity` return `NONFINITE_RESULT`. Other malformed JSON returns `INVALID_INPUT`. Values are never clipped or replaced.

`E2-CAN-001` Canonical JSON has no insignificant whitespace and no trailing newline. Object keys are sorted by increasing Unicode scalar sequence. Array order is preserved.

`E2-CAN-002` Strings are emitted as UTF-8. Only quotation mark, reverse solidus, and U+0000–U+001F are escaped. Quotation mark and reverse solidus use `\"` and `\\`; controls use lowercase four-digit `\u00xx`. Solidus and non-ASCII scalars are not escaped.

`E2-CAN-003` `null`, `true`, and `false` use those lowercase tokens. Portable integers use base-10 without leading zeros; zero is `0`.

`E2-CAN-004` Binary64 values use RFC 8785 JSON-number serialization: the shortest decimal that round-trips to the same binary64 value under round-to-nearest ties-to-even; negative zero serializes as `0`; lowercase `e`; no leading zero in the exponent; positive exponents include `+`. Nonfinite values are forbidden.

`E2-CAN-005` SHA-256 digests use lowercase 64-digit hexadecimal and the prefix `sha256:`. The hash input is exactly the canonical UTF-8 bytes, with no BOM or newline.

`E2-CAN-006` The schema-descriptor payload is exactly:

```json
{"base_spec_sha256":"8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84","errata_1_sha256":"52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40","ir":"spark-atomistic-model/1","revision":2}
```

Its `schema_digest` is `sha256:583d580d54e3847ef92f1b1456dda006161689c0bac27fd7ea896a093f48c02c`.

`E2-CAN-007` `config_digest` hashes the complete validated model after removing root `metadata`. The run-request `extension` is outside the model and excluded. `tolerance_digest` hashes the exact object containing the nine tolerance/cutoff fields listed in `kinetics`, in their canonical key order.

`E2-CAN-008` `identity_digest` hashes exactly `{"reflection_invariant":false,"rotation_invariant":false,"state_energy_tolerance_per_atom":value,"state_max_tolerance":value,"state_rms_tolerance":value,"version":"spark-state-identity/2"}`. `calculator_digest` equals the configured calculator model digest.

## 5. Exact status, severity, exit, and context

`E2-STATUS-001` A stored status record has exactly `causal_status`, `context`, `exit_code`, `message`, `severity`, and `status`. `exit_code` is an integer for a terminal record and `null` for a nonterminal record. Context has exactly `component`, `details`, `requirement_id`, `retryable`, `search_or_event_id`, and `state_id`. IDs are strings or `null`; `details` is an object; other fields have the types shown in `E2-API-006`. A public response adds exactly `operation` and `value` to this record.

`E2-STATUS-002` The exact vocabulary, severity, final exit, and stable message are:

| Status | Severity | Final exit | Exact message |
|---|---|---:|---|
| `OK` | `success` | `0` | `transaction committed` |
| `DISCOVERY_CONVERGED_HEURISTIC` | `success-with-qualification` | `0` | `heuristic discovery criterion passed` |
| `DUPLICATE_EVENT` | `candidate reject` | — | `duplicate event rejected` |
| `SADDLE_NOT_FOUND` | `candidate reject` | — | `saddle not found` |
| `INVALID_SADDLE` | `candidate reject` | — | `saddle validation failed` |
| `SADDLE_WRONG_BASIN` | `candidate reject` | — | `neither endpoint matches origin` |
| `ENDPOINT_COLLAPSED` | `candidate reject` | — | `both endpoints match origin` |
| `ENVIRONMENT_AMBIGUOUS` | `recoverable` | — | `environment identity ambiguous` |
| `BASIN_DISABLED` | `recoverable` | — | `basin acceleration disabled` |
| `DISCOVERY_INCOMPLETE` | `pause/qualified` | `75` strict; `0` exploratory | `discovery budget exhausted` |
| `RELAX_NOT_CONVERGED` | `transaction fail` | `65` after policy exhaustion | `relaxation not converged` |
| `EVENT_APPLICATION_FAILED` | `transaction fail` | `65` | `event application failed` |
| `CALCULATOR_FAILURE` | `transaction fail` | `69` after policy exhaustion | `calculator callback failed` |
| `NONFINITE_RESULT` | `fatal` | `65` | `nonfinite value rejected` |
| `INVALID_INPUT` | `fatal` | `64` | `input invalid` |
| `SCHEMA_UNSUPPORTED` | `fatal` | `64` | `schema unsupported` |
| `INVALID_STATE` | `fatal` | `65` | `state invalid` |
| `RATE_INVALID` | `fatal in strict mode` | `65` when terminal | `rate invalid` |
| `DETAILED_BALANCE_VIOLATION` | `fatal in strict mode` | `65` when terminal | `detailed balance violated` |
| `CATALOG_CONFLICT` | `fatal` | `65` | `catalog conflict` |
| `CATALOG_INCOMPATIBLE` | `fatal` | `65` | `catalog incompatible` |
| `ATOM_COUNT_CHANGE_UNSUPPORTED` | `fatal` | `65` | `atom count change unsupported` |
| `NO_ENABLED_EVENT` | `terminal-success if requested, else fatal` | `0` if requested; else `65` | `no enabled event` |
| `RESOURCE_LIMIT` | `pause` | `75` | `resource limit reached` |
| `OUTPUT_EXISTS` | `fatal` | `64` | `output exists` |
| `CHECKPOINT_CORRUPT` | `fatal` | `74` | `checkpoint corrupt` |
| `CHECKPOINT_INCOMPATIBLE` | `fatal` | `74` | `checkpoint incompatible` |
| `CANCELLED` | `pause` | `75` | `run cancelled` |
| `INTERNAL_ERROR` | `fatal` | `70` | `internal error` |

`E2-STATUS-003` An em dash in the table means the status cannot terminate a public operation and has JSON `exit_code:null` in its stored internal record. Public final responses always contain an integer exit code.

`E2-STATUS-004` The earliest causal failure remains top-level. `causal_status` records the immediately preceding exact status token or `null`. Cleanup/checkpoint failure details are appended under `context.details`; exit `74` replaces the mapped exit only when no valid last-committed checkpoint remains.

## 6. Canonical identities and digests

`E2-ID-001` All identity payloads use canonical JSON. All prefixes are exact: `candidate:sha256:`, `state:sha256:`, `pair:sha256:`, `event:sha256:`, `search:sha256:`, and `sha256:` for non-ID digests.

`E2-ID-002` A canonical geometry certificate is constructed as follows:

1. For each atom as anchor, subtract its position from every position.
2. For each periodic displacement, choose the closest lattice image. A norm tie chooses the lexicographically smallest integer lattice-shift triple.
3. Form rows `[species,movable,dx,dy,dz]`; sort rows lexicographically using canonical scalar encodings.
4. Choose the lexicographically smallest complete row array over all anchors.
5. The certificate is exactly `{"cell":cell,"pbc":pbc,"rows":rows,"version":"anchor-minimum-closest-image/1"}`.

This is invariant to whole-cell translation, periodic image choice, atom ID, and same-species permutation. Rotation/reflection equivalence is disabled.

`E2-ID-003` `fixed_contract_digest` hashes exactly:

```json
{"atom_contracts":[],"calculator_model_digest":"","cell":[],"charge":0,"constraints":{"kind":"fixed-mask"},"pbc":[],"spin":0}
```

`atom_contracts` is the sorted array of `[species,movable]` rows. The placeholders above are replaced by actual values.

`constraint_digest` hashes exactly `{"constraints":{"kind":"fixed-mask"},"movable":movable}` with the actual movable array.

`E2-ID-004` `candidate_identity` hashes exactly `{"fixed_contract_digest":digest,"geometry":certificate,"version":"spark-state-identity/2"}`. `state_id` hashes exactly `{"candidate_identity":id,"energy_ev":energy,"fixed_contract_digest":digest,"version":"spark-state-identity/2"}`.

`E2-ID-005` Candidate IDs are hints only. Geometry and energy verification decides equivalence. A passing candidate takes the already committed state ID; a nonmatching candidate keeps its computed state ID. Atom IDs alone never establish equality.

`E2-ID-006` A committed state record has exactly `atom_ids`, `calculator_model_digest`, `candidate_identity`, `cell`, `charge`, `constraint_digest`, `constraints`, `energy_ev`, `fixed_contract_digest`, `force_tolerance_ev_per_angstrom`, `forces_ev_per_angstrom`, `identity_version`, `max_movable_force_ev_per_angstrom`, `movable`, `pbc`, `positions`, `relaxation_provenance`, `schema`, `species`, `spin`, and `state_id`. `schema` is `spark-atomistic-state/2`; `identity_version` is `spark-state-identity/2`. Relaxation provenance has exactly `calculator_evaluations`, `calculator_identity`, `minimizer_identity`, `steps`, and `termination_reason`.

`E2-ID-007` `saddle_geometry_digest` is `sha256:` plus the SHA-256 of the §6 canonical geometry certificate after replacing minimum positions by saddle positions while retaining the same cell, periodicity, species, and movable mask.

## 7. Philox and deterministic substreams

`E2-RNG-001` The exact algorithm ID is `Philox4x32-10:errata-1-midpoint52`. Philox round constants, lane order, `(w0,w1)/(w2,w3)` pairing, 52-bit `q`, and midpoint mapping are Errata 1.

`E2-RNG-002` A Philox state has exactly `algorithm`, `buffered_block`, `consumed_blocks`, `consumed_uniforms`, `initial_counter`, `key`, `next_counter`, and `next_pair`. `key` is two unsigned 32-bit integers; counters and a present block are four unsigned 32-bit integers; `buffered_block` is `null` or a block; `next_pair` is `0` or `1`. With no buffered block, `next_pair=0`. After consuming the first pair, the generated block is retained and `next_pair=1`; after the second pair, the block becomes `null` and `next_pair=0`. `consumed_blocks=ceil(consumed_uniforms/2)` and `next_counter=initial_counter+consumed_blocks` under §7 counter arithmetic.

`E2-RNG-003` Counter lane `c0` is least significant. Block increment adds one to the unsigned 128-bit value `c0 + 2^32 c1 + 2^64 c2 + 2^96 c3`, propagating carry toward `c3`.

`E2-RNG-004` Let `u64be`, `u32be`, and `len32(s)` mean unsigned big-endian encodings and the UTF-8 byte length encoded by `u32be`. The trajectory digest material is:

```text
UTF8("spark-trajectory-stream/2\0") || u64be(run_seed)
```

The saddle/class substream digest material is:

```text
UTF8("spark-saddle-substream/2\0") || u64be(run_seed) || len32(state_id) || UTF8(state_id) || len32(search_class) || UTF8(search_class) || u64be(search_index)
```

SHA-256 bytes `0..7` become two big-endian key words; bytes `8..23` become four big-endian initial-counter words; remaining bytes are discarded. The next counter initially equals the initial counter.

`E2-RNG-005` Search-class choice uses the first uniform of the substream whose `search_class` is exact string `class-selection`. The chosen saddle search uses a separate substream derived with the chosen configured class name. Neither consumes the trajectory stream.

`E2-RNG-006` The first boundary golden record is this exact canonical JSON line:

```json
{"a":0,"b":0,"q":0,"raw_binary64_bits":"0x3ca0000000000000","uniform_hex":"0x1.0000000000000p-53"}
```

SHA-256 of the line without newline: `6ce1fb5214530ba6b04e4bf75aaeba5d02acf6694cd462004faf7640a665fc03`.

`E2-RNG-007` The second boundary golden record is this exact canonical JSON line. `b=0xfffff000`; its low 12 discarded bits are zero.

```json
{"a":4294967295,"b":4294963200,"q":4503599627370495,"raw_binary64_bits":"0x3fefffffffffffff","uniform_hex":"0x1.fffffffffffffp-1"}
```

SHA-256 of the line without newline: `a15157604d319e5525e3b83eba02259088e317ea2b0b0e1a3bb28060e093cf43`.

## 8. Directed reversible event record and rates

`E2-EVENT-001` A directed event record has exactly:

| Field | Rule |
|---|---|
| `schema` | `spark-atomistic-directed-event/2`. |
| `event_id` | Canonical directed ID. |
| `reverse_event_id` | Reciprocal directed ID. |
| `pair_id` | Shared physical saddle-pair ID. |
| `origin_state_id`, `destination_state_id` | Distinct committed IDs. |
| `saddle` | Exact saddle object below. |
| `barrier_ev`, `reverse_barrier_ev` | Raw same-saddle differences. |
| `rate_model` | Exact common-prefactor object below. |
| `selectable` | Boolean at checkpoint temperature/cutoff. |
| `active_atom_mapping` | Sorted origin-to-destination `[origin_index,destination_index]` pairs; species-preserving partial bijection. |
| `environment_key`, `environment_version` | Exact `disabled` and `none/1`. |
| `discovery_provenance` | Exact object below. |
| `validation` | Exact object below. |
| `calculator_digest`, `identity_digest`, `schema_digest`, `tolerance_digest` | Nonempty canonical digests. |

`E2-EVENT-002` `saddle` has exactly `curvature_ev_per_angstrom2`, `energy_ev`, `evaluation_count`, `evidence_level`, `forces_ev_per_angstrom`, `orthogonal_curvatures_ev_per_angstrom2`, `positions`, `search_id`, `termination_reason`, and `unstable_direction`. Evidence is `HESSIAN` or `DIRECTIONAL`.

`E2-EVENT-003` `rate_model` has exactly `common_prefactor_per_s`, `detailed_balance_residual`, `log_forward_rate_per_s`, `log_reverse_rate_per_s`, `model`, and `temperature_k`. `model` is `COMMON_PREFACTOR`. For the reciprocal record, forward/reverse logs swap and the residual changes sign.

`E2-EVENT-004` `discovery_provenance` has exactly `rng_substream_digest`, `search_class`, `search_id`, and `search_index`.

`E2-EVENT-005` `validation` has exactly `calculator_model_digest`, `constraint_digest`, `destination_match`, `full_endpoint_relaxations`, `method`, `origin_match`, and `unstable_mode_count`. A match has exactly `atom_mapping`, `energy_difference_ev`, `max_displacement_angstrom`, and `rms_displacement_angstrom`. `full_endpoint_relaxations=true` and `unstable_mode_count=1` are mandatory.

`E2-EVENT-006` Let the endpoint IDs be sorted. Orient the pair mapping from the lower ID to the higher ID. Canonicalize unstable-direction sign by choosing the lexicographically smaller canonical flattened array of `v` and `-v`. `pair_id` hashes exactly:

```json
{"active_atom_mapping":[],"endpoint_state_ids":[],"saddle_energy_ev":0,"saddle_geometry_digest":"sha256:","schema":"spark-atomistic-event-pair/2","unstable_direction":[]}
```

Placeholders are replaced by actual values. `event_id` hashes exactly `{"destination_state_id":destination,"origin_state_id":origin,"pair_id":pair,"schema":"spark-atomistic-directed-event/2"}`.

`E2-RATE-001` Raw barriers are exactly `b_f=E_s-E_i` and `b_r=E_s-E_j`. If either is less than `-barrier_tolerance`, return `RATE_INVALID`. If `-barrier_tolerance <= b < 0`, retain the exact negative value in the record and rate formula. Clamping to zero is forbidden.

`E2-RATE-002` Rates use the raw barriers and the base-spec `COMMON_PREFACTOR` equations. A log rate below `log_rate_cutoff` is nonselectable. Every snapshot containing such disabled rates reports their summed lost-rate log upper bound.

## 9. Deterministic discovery

`E2-DISC-001` For state ID and zero-based search index, derive the `class-selection` substream per §7, consume its first uniform `u`, and traverse the configured class array in input order. Choose the first class whose cumulative probability is strictly greater than `u`; if rounding leaves none, choose the last class.

`E2-DISC-002` Search ID hashes exactly `{"run_seed":seed,"search_class":name,"search_index":index,"state_id":state_id}` and uses prefix `search:sha256:`. The saddle substream uses the chosen class name. Retry uses the identical search ID and substream state.

`E2-DISC-003` Parallel completion cannot alter catalog visibility. Candidate transactions are committed in ascending search ID order, then deduplicated transactionally.

`E2-DISC-004` Discovery statistics have exactly `alpha`, `alpha_calibration`, `attempts`, `config_digest`, `consecutive_redundant_successes`, `duplicates`, `evaluations`, `event_log_rates`, `failures_by_status`, `heuristic_confidence`, `permanently_incomplete_catalog`, `relevance_rate_min`, `state_id`, `stopping_state`, and `successes`. `stopping_state` is `RUNNING`, `CONVERGED_HEURISTIC`, or `INCOMPLETE`. Confidence is a finite number or exact string `UNAVAILABLE`.

`E2-DISC-005` A failed/invalid/wrong-basin attempt increments `attempts` and one `failures_by_status` count, resets consecutive redundant successes, and never increments success/duplicate counts. A new valid event increments success and resets redundancy. A relevant duplicate increments success, duplicate, and consecutive redundancy without adding rate twice.

## 10. Serial KMC and historical rate snapshot

`E2-KMC-001` Before a step, selectable directed events from the current state are sorted by event ID. Their linear rates are generated from stored log rates and summed by Neumaier compensated summation in that order.

`E2-KMC-002` A rate-table snapshot envelope has exactly `payload` and `payload_sha256`. Payload has exactly `destination_state_ids`, `event_ids`, `log_rates`, `lost_rate_log_upper_bound`, `origin_state_id`, `rates`, `schema`, and `total_rate_per_s`. Schema is `spark-atomistic-rate-table-snapshot/1`; parallel arrays use event-ID order; lost bound is finite or `null`.

`E2-KMC-003` Each proposed step clones trajectory RNG, consumes selection then time uniforms, selects the first cumulative rate strictly greater than `u_s K`, and computes `-ln(u_t)/K`. Event application and verification occur before committing state, time, step, log, and cloned RNG.

`E2-KMC-004` A committed step record has exactly `checkpoint_sequence`, `log_sequence`, `post_state_id`, `pre_state_id`, `rate_table_snapshot`, `selected_event_id`, `selected_rate_per_s`, `selection_uniform`, `step_index`, `time_increment_s`, `time_uniform`, and `total_rate_per_s`.

`E2-KMC-005` The per-step rate-table snapshot is historical and immutable. Restore/replay MUST use that snapshot, never a later expanded catalog or current rate cache, to re-evaluate the earlier selection.

## 11. Checkpoint wire schema and recursive restore

`E2-CKPT-001` Checkpoint envelope has exactly:

```json
{"payload":{},"payload_sha256":"sha256:"}
```

The hash covers canonical payload bytes.

`E2-CKPT-002` Payload has exactly `basin`, `catalog`, `checkpoint_sequence`, `digests`, `discovery_statistics`, `flags`, `initial_state`, `log_sequence`, `resources`, `rng`, `schema`, `simulation_time_s`, `step_index`, `trajectory`, and `current_state`. Schema is `spark-atomistic-checkpoint/2`.

`E2-CKPT-003` `digests` has exactly `config`, `model`, `schema`, and `tolerances`. `rng` has exactly `run_seed`, `substream_map`, and `trajectory`; each stream uses §7. `basin` is exactly `{"enabled":false,"reason":"v1-disabled"}`.

`E2-CKPT-004` `flags` has exactly `cancelled`, `complete`, `incomplete_catalog`, `last_status`, and `resource_limited`. Booleans are exact; `last_status` is an exact status token.

`E2-CKPT-005` `resources` has exactly `calculator_evaluations`, `catalog_events`, `output_bytes`, `resident_memory_bytes`, `retry_history`, `saddle_attempts_by_state`, and `wall_elapsed_s`. Retry history is empty in v1 because configured retry count is zero.

`E2-CKPT-006` `catalog` has exactly `digest`, `events`, `multiplicity`, `schema`, and `states`. Schema is `spark-atomistic-catalog/2`. States/events/multiplicity are objects keyed by their IDs. Catalog digest hashes the same object without `digest`.

`E2-CKPT-007` Restore is validate-before-mutation and recursively verifies, in order:

1. strict/canonical JSON and envelope shape;
2. payload SHA-256;
3. schema/config/model/tolerance digests;
4. complete state records, state IDs, fixed contracts, forces, tolerances, and provenance;
5. catalog digest, event IDs, reciprocal pairs, saddle/barrier/rate recomputation, detailed balance, mappings, and validation evidence;
6. discovery counter identities, stopping state, rates, and config binding;
7. resource counters/limits and zero retry history;
8. Philox key/counter/buffer/count relations and complete substream map;
9. trajectory sequence/state chain using each historical rate snapshot;
10. final state, simulation time, step/log counts, and trajectory RNG equality.

Any failure before mutation returns `CHECKPOINT_CORRUPT` for damaged content or `CHECKPOINT_INCOMPATIBLE` for a valid but mismatched run contract.

`E2-CKPT-008` Step and log sequences equal trajectory length. Trajectory RNG consumed uniforms equal `2 * step_index`. Replayed IDs, uniforms, rate arrays, selections, and counters are byte-identical. Replayed time uses base `PAR-005` tolerance.

`E2-CKPT-009` Crash-safe write remains sibling temporary creation, content flush, atomic replacement, and parent-directory flush. Canonical bytes are the envelope bytes with no trailing newline.

## 12. Basin capability

`E2-BASIN-001` Basin acceleration is uniformly unavailable in v1. Capability reports `disabled`; model `enabled` is only a request flag; checkpoint stores the disabled record; no basin path appears in a trajectory record.

`E2-BASIN-002` Any attempt to enable or invoke basin acceleration returns nonterminal `BASIN_DISABLED`, preserves state/time/trajectory RNG/step, and continues with serial KMC when a run is active.

`E2-BASIN-003` No exact, partial, approximate, pathwise, mean-rate, or phase-type basin capability may be advertised under this IR revision. A later enabled capability requires a new IR or normative extension.

## 13. Parity golden fixtures

`E2-PAR-001` Every conforming implementation MUST consume the same canonical request corpus and emit the same canonical response, status record, state/event/catalog record, rate snapshot, RNG state, and checkpoint schema. No mandatory fixture may be skipped.

`E2-PAR-002` The mandatory corpus includes:

1. capability request and exact capability value;
2. one valid minimal model and metadata-only variants with identical config digest;
3. duplicate-key, malformed UTF-8, lone-surrogate, nonfinite, integer-boundary, out-of-domain integer, unknown-key, and adapter-extension cases;
4. canonical-number cases including negative zero, `1e-6`, values adjacent to `1e-6`, `1e21`, and values adjacent to `1e21`;
5. translated, periodic-image, and same-species-permuted state identities;
6. both Errata 1 boundary records, Philox zero-key/counter words, trajectory derivation, saddle derivation, state restore, and counter-carry cases;
7. reversible event pair, reciprocal mapping, tolerated negative barrier, rejected negative barrier, detailed-balance failure, and parallel saddle cases;
8. deterministic discovery class sequence and out-of-order completion commit;
9. serial two-event KMC selection, application rollback, lost-rate bound, and historical snapshot replay after catalog growth;
10. clean checkpoint, corrupted hash, incompatible digest, recursive nested corruption, cancellation, resource limit, and exact next-event resume;
11. basin-disabled capability, input rejection, checkpoint record, and serial fallback.

`E2-PAR-003` Canonical JSON, SHA-256, schema/config/tolerance digests, ID strings, Philox words/uniform bits/counters, selected-event sequence, status sequence, event records, rate snapshots, and checkpoint bytes are byte-identical.

`E2-PAR-004` Energies, forces, barriers, rates, and geometries use base `PAR-004`; time uses base `PAR-005`. Fixture expectations store both canonical decimal and raw binary64 bits for RNG uniforms.

`E2-PAR-005` Static parity readiness requires zero Critical, High, or Medium observable mismatches against this erratum, plus all mandatory fixtures present. Execution conformance additionally requires every fixture to pass in every implementation.

## 14. Effect

This erratum replaces divergent language-specific wire contracts with one implementation-neutral revision. Existing artifacts that do not match the exact schema, IDs, RNG derivation, event records, historical snapshots, or checkpoint schema are incompatible and MUST be rejected rather than migrated implicitly.

## Errata 2 digest

Hash rule: SHA-256 of the exact UTF-8 bytes before the line beginning `## Errata 2 digest`; the separator blank line is included and this digest section is excluded.

Errata 2 SHA-256: `eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995`
