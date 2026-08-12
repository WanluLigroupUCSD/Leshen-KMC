# Off-lattice/on-the-fly KMC neutral specification v1

Date: 2026-08-09  
Status: clean-room behavioral/math contract  
Normative words: `MUST`, `MUST NOT`, `SHOULD`, `MAY`  
Evidence map: `../literature_notes/SPARK_OFFLATTICE_OTF_CLEANROOM_REVIEW_2026-08-09.md`

## 1. Scope and priority

This specification defines interfaces and observable behavior. It contains no source code or implementation-derived symbol names.

| Priority | Included |
|---|---|
| P0 | Fixed atom count, fixed species, fixed cell; off-lattice minima; energy/force calculator callback; minimizer and first-order-saddle callbacks; on-the-fly event discovery; event catalog and deduplication; reversible rates; serial residence-time KMC; deterministic RNG/checkpoint; Python/Rust conformance. |
| P1 | Exact finite-superbasin exit sampler; calibrated missing-rate estimator; topology/local-environment reuse; harmonic prefactors when a mode provider exists. |
| Out of v1 | Atom insertion/deletion, adsorption/desorption to a reservoir, variable cell, hybrid MD, Time Warp, GPU execution, tau-leap, chemical-reaction templates, quantum nuclear effects. |

`SCOPE-001` P0 MUST not require a lattice, predefined sites, or a predefined event table.

`SCOPE-002` P0 MUST treat topology changes caused by continuous atomic motion, bond rearrangement, concerted diffusion, and surface reconstruction. Bonds are derived observations, not state variables.

`SCOPE-003` P0 MUST keep atom count, species sequence modulo same-species permutation, simulation cell, periodicity flags, and constraints fixed throughout one trajectory.

`SCOPE-004` Time Warp, GPU kernels, tau-leap, and parallel-discrete-event execution MUST NOT be P0 dependencies. Independent saddle searches MAY run concurrently, but scheduling MUST NOT alter accepted events.

## 2. Mathematical model

A state (i) is a locally minimized atomic configuration. A directed event (i\rightarrow j) is a validated first-order saddle (s) connected by downhill relaxation to (i) and (j). At fixed temperature (T>0), enabled events form a continuous-time Markov chain (CTMC).

For rates (k_{ij}>0), (K_i=\sum_j k_{ij}). Two independent uniforms (u_s,u_t\in(0,1)) select the first event in canonical order whose cumulative rate exceeds (u_sK_i), and advance time by

\[
\Delta t=-\ln(u_t)/K_i.
\]

`MATH-001` The engine MUST state that trajectory correctness is conditional on Markovian state-to-state dynamics, valid supplied rates, and sufficient discovery of relevant exits.

`MATH-002` A finite saddle-search budget MUST NOT be described as proof that the event catalog is complete.

`MATH-003` All energy calculations use eV, length Å, force eV Å(^{-1}), time s, rate s(^{-1}), temperature K, and (k_B=8.617333262145\times10^{-5}\) eV K(^{-1}).

## 3. State contract

`STATE-001` A state request MUST contain: schema version; atom IDs; species; Cartesian positions; (3\times3) cell; three periodicity flags; fixed/movable constraints; charge/spin metadata needed by the calculator; and calculator model digest.

`STATE-002` Atom IDs are trajectory-stable labels, not physical identity claims. Same-species permutation equivalence MUST be handled during structure matching; IDs alone MUST NOT establish state equality.

`STATE-003` A committed state MUST contain finite relaxed energy, finite forces, maximum movable-atom force, minimizer tolerance, constraint digest, and relaxation provenance.

`STATE-004` A committed state MUST satisfy `max_force <= force_tolerance`. Otherwise the status is `RELAX_NOT_CONVERGED`, and the state MUST NOT enter the kinetic catalog.

`STATE-005` State equivalence MUST be invariant to periodic image choice, whole-cell translation, and same-species atom permutation. Optional rotational or space-group equivalence MUST be declared in configuration and versioned.

`STATE-006` State matching MUST use both a discrete candidate identity and a geometry verification. Geometry verification MUST report atom mapping, cell-aware displacement norm, maximum displacement, and energy difference.

`STATE-007` Two minima are equal only when species/cell/constraints match and configured geometry and energy tolerances pass. Default tolerances: RMS displacement `1e-3 Å`, maximum displacement `5e-3 Å`, energy difference `1e-6 eV/atom`. Overrides are part of the run digest.

`STATE-008` Nonfinite coordinates, energy, force, cell, or tolerance MUST return `NONFINITE_RESULT`; silently replacing or clipping values is forbidden.

## 4. Local-environment identity

Local identity is a cache key, never a physical approximation by itself.

`ENV-001` An environment definition MUST declare center selection, radial extent, neighbor rule, element labels, periodic-image rule, and identity-version string.

`ENV-002` The identity MUST be invariant to global translation, periodic image choice, and permutation of indistinguishable atoms. Rotation/reflection invariance MUST be declared and tested if enabled.

`ENV-003` The identity MAY use a labeled neighbor graph, geometric invariant, or their combination. No particular graph library or fingerprint algorithm is normative.

`ENV-004` Neighbor rules MUST define an ambiguity band or hysteresis rule. A geometry inside an unresolved band returns `ENVIRONMENT_AMBIGUOUS` and MUST bypass generic-event reuse.

`ENV-005` Identical local keys MUST NOT imply identical barriers. A reused event MUST be reconstructed in the current full geometry and reconverged before its rate becomes selectable.

`ENV-006` If one key maps to geometries that fail reconstruction or geometry verification, the implementation MUST split the identity, version the refinement, and retain provenance. Silent key collision is forbidden.

`ENV-007` Long-range elastic or electrostatic changes outside the local region MUST be included through full-state calculator evaluation during reconvergence.

`ENV-008` Cache efficacy MUST be reported: key count, reuse attempts, successful reconvergences, rejected mappings, ambiguous identities, and fresh searches. No minimum hit rate is required.

## 5. Calculator, minimizer, and saddle interfaces

### 5.1 Calculator callback

`CALC-001` The calculator request MUST provide the complete state contract and requested properties. P0 properties are total potential energy and all atomic forces.

`CALC-002` A successful calculator response MUST provide status `OK`, energy, forces with shape (N\times3), units, model name/version, model digest, evaluation ID, and deterministic flag.

`CALC-003` The callback MUST NOT change atom count, species, atom IDs, cell, periodicity, or constraints.

`CALC-004` Same canonical request and model digest SHOULD produce identical results. If not, the response MUST set `deterministic=false`; bitwise replay claims then become invalid, while tolerance-based replay remains allowed.

`CALC-005` Optional properties are stress, Hessian, mass-weighted stable-mode logs, and uncertainty. Their absence MUST NOT block constant-prefactor P0.

`CALC-006` Callback timeout, process failure, malformed shape, unit mismatch, or model-digest change returns `CALCULATOR_FAILURE`. Partial energy/force data MUST be discarded.

### 5.2 Minimizer callback

`RELAX-001` Input MUST contain state, force tolerance, step/evaluation limits, constraints, and calculator callback identity.

`RELAX-002` Output MUST contain status, relaxed state, energy, forces, maximum movable-atom force, step count, calculator evaluations, and termination reason.

`RELAX-003` `OK` is legal only if `STATE-004` passes. Budget exhaustion returns `RELAX_NOT_CONVERGED`; calculator errors propagate as `CALCULATOR_FAILURE`.

`RELAX-004` Endpoint relaxations MUST use the same calculator model digest and constraints as the originating minimum and saddle.

### 5.3 Saddle-search callback

`SADDLE-001` Input MUST contain the committed minimum, active-region hint or `none`, independent search ID, RNG substream, convergence tolerances, resource limits, and calculator/minimizer identities.

`SADDLE-002` The interface MUST permit minimum-mode following, dimer-family, ART-family, or another first-order-saddle solver without changing the event contract.

`SADDLE-003` A successful candidate MUST report saddle geometry/energy/forces, unstable direction, curvature evidence, search provenance, and two perturbed downhill endpoint requests.

`SADDLE-004` Validation MUST relax both sides. Exactly one endpoint MUST match the origin state and the other MUST be a different committed minimum. Zero matching endpoints returns `SADDLE_WRONG_BASIN`; two matching endpoints returns `ENDPOINT_COLLAPSED`.

`SADDLE-005` A valid saddle MUST have one unstable mode under the declared validation method. If a Hessian is available, exactly one nonzero imaginary mode is required after excluding constrained and rigid modes. Without a Hessian, negative curvature along the reported mode plus nonnegative sampled orthogonal curvatures is required and the evidence level MUST be `DIRECTIONAL`, not `HESSIAN`.

`SADDLE-006` Failed searches, invalid saddles, and wrong-basin saddles MUST NOT count as redundant discoveries for a convergence claim.

`SADDLE-007` The solver MUST return its termination reason and evaluation counts for every attempt, including failures.

## 6. Event contract, discovery, catalog, and deduplication

### 6.1 Event contract

`EVENT-001` A committed reversible event MUST contain immutable event ID; origin and destination state IDs; saddle record; forward and reverse barriers; rate-model record; active-atom mapping; environment key/version; discovery provenance; validation evidence; and reverse-pair ID.

`EVENT-002` Barriers MUST be derived from the same saddle: (\Delta E^\ddagger_{ij}=E_s-E_i\) and (\Delta E^\ddagger_{ji}=E_s-E_j\). Negative barriers below `-barrier_tolerance` return `RATE_INVALID`.

`EVENT-003` A process with the same endpoints but a geometrically distinct saddle is a distinct event. Parallel mechanisms MUST remain separate because rates add.

`EVENT-004` Event application MUST use the validated destination minimum, followed by one verification relaxation. Failure to recover the destination within state tolerances returns `EVENT_APPLICATION_FAILED`; no kinetic time is committed.

### 6.2 Discovery

`DISC-001` Each new state MUST launch a configured mixture of global, local, and optional targeted initial perturbations. At least one non-targeted class MUST have nonzero probability.

`DISC-002` Every attempt MUST use a scheduling-independent RNG substream derived from run seed, state ID, and search index.

`DISC-003` A relevant event window MUST be defined by rate, not barrier alone, when prefactors differ. Constant-prefactor mode MAY use a barrier window.

`DISC-004` P0 stopping uses all of: a minimum successful-search count; a consecutive redundant-success count; a maximum search/evaluation budget; and an explicit relevance window.

`DISC-005` If the redundant-success criterion passes, status is `DISCOVERY_CONVERGED_HEURISTIC`. The report MUST include the assumed minimum relative discovery probability `alpha`, redundant count (N_r), and heuristic confidence (C=1-1/(\alpha N_r)). `alpha` MUST come from a declared calibration; otherwise `C` is `UNAVAILABLE`.

`DISC-006` Budget exhaustion before the stopping criterion returns `DISCOVERY_INCOMPLETE`. In strict mode no KMC step may occur. In exploratory mode a step MAY occur only with the trajectory and checkpoint permanently tagged `INCOMPLETE_CATALOG`.

`DISC-007` P1 MAY estimate missing total rate (k_U) and catalog validity time (\tau_V=-\ln(1-\delta)/k_U\). It MUST identify the estimator, calibration, target error `delta`, and uncertainty. An estimate is not an exact bound unless its proof assumptions are verified.

`DISC-008` Discovery statistics MUST survive checkpoint/resume: attempts, successes, failures by status, duplicate counts, event rates, relevance window, calibration, and stopping state.

### 6.3 Deduplication and reuse

`CAT-001` Deduplication MUST compare origin equivalence, destination equivalence, saddle geometry, unstable direction up to sign, active-atom mapping, and barrier tolerance.

`CAT-002` Default event tolerances are saddle RMS displacement `1e-3 Å`, saddle maximum displacement `5e-3 Å`, and saddle-energy difference `1e-5 eV`. Overrides are versioned.

`CAT-003` A passing duplicate returns `DUPLICATE_EVENT`, increments multiplicity/search statistics, and MUST NOT add its rate twice.

`CAT-004` Same endpoints with nonmatching saddles MUST be retained as parallel events. Same saddle with inconsistent energies beyond tolerance returns `CATALOG_CONFLICT`.

`CAT-005` Generic/local event reuse MAY propose a candidate. Only a fully reconverged, endpoint-validated specific event may enter the selectable catalog.

`CAT-006` Catalog records MUST include calculator, identity, tolerance, and schema digests. Digest mismatch forbids reuse and returns `CATALOG_INCOMPATIBLE`.

`CAT-007` Catalog mutation is transactional: a record becomes visible only after state, saddle, endpoints, rate, deduplication, and checksum validations pass.

## 7. Rates and detailed balance

`RATE-001` P0 MUST implement `COMMON_PREFACTOR` and MAY implement `HARMONIC_TST` or externally supplied free-energy rates.

`RATE-002` In `COMMON_PREFACTOR`, one positive finite (\nu) is shared by both directions of a reversible saddle pair:

\[
k_{ij}=\nu e^{-\beta(E_s-E_i)},\qquad k_{ji}=\nu e^{-\beta(E_s-E_j)}.
\]

`RATE-003` In `HARMONIC_TST`, the rate is

\[
k_{ij}=\frac{\prod_m \nu_{i,m}}{\prod'_m \nu_{s,m}}
e^{-\beta(E_s-E_i)},
\]

where constrained/rigid zero modes are removed and the saddle unstable mode is excluded. Mode-selection policy and masses MUST be recorded.

`RATE-004` Rates MUST be computed and summed in log-safe form. A selectable rate MUST be finite and strictly positive. Underflowed rates MAY be disabled only if the rate cutoff and lost-rate upper bound are reported.

`RATE-005` For equilibrium reversible pairs, detailed balance MUST satisfy

\[
\left|\ln(k_{ij}/k_{ji})+\beta(F_j-F_i)\right|\le \epsilon_{DB}.
\]

Default `epsilon_DB=1e-8`. (F=E) in common-prefactor energy mode; a free-energy rate model MUST supply consistent (F).

`RATE-006` Detailed-balance failure returns `DETAILED_BALANCE_VIOLATION`; the pair MUST NOT be selectable in strict mode.

`RATE-007` Irreversible reservoir or driven events MUST NOT be encoded as reversible P0 events. They require a later protocol defining activities, chemical potentials, degeneracy/proposal factors, and local detailed balance.

`RATE-008` Temperature changes invalidate all cached numeric rates but not validated barriers/prefactor inputs. The recalculation event MUST be logged.

## 8. Serial KMC selection

`KMC-001` Selectable events MUST be sorted by immutable event ID before every rate sum and cumulative selection.

`KMC-002` Rate summation MUST use a declared deterministic order and compensated or pairwise summation. Total rate `0` returns `NO_ENABLED_EVENT`; nonfinite total returns `RATE_INVALID`.

`KMC-003` Each committed KMC step consumes exactly two uniforms from the trajectory stream: selection, then time. Rejected discovery/callback work MUST NOT consume this stream.

`KMC-004` The selected event, uniforms, total rate, individual selected rate, time increment, pre/post state IDs, and checkpoint sequence MUST be logged before the next step.

`KMC-005` Event and time commit are atomic. If event application fails, state, time, RNG counter, and step index MUST remain at their pre-step values.

`KMC-006` P0 is time-homogeneous during one selection. Time-dependent rates or scheduled external events are out of v1.

## 9. Basin/superbasin acceleration

`BASIN-001` Basin acceleration is disabled by default and is not required for P0 publication.

`BASIN-002` An optional exact accelerator MUST represent a finite set of transient states, all known internal rates, all known exit rates, the entry-state distribution, and the transient CTMC generator (Q).

`BASIN-003` Exact pathwise acceleration MUST sample the joint first-exit time and exit channel from the phase-type law generated by (Q), or use a mathematically equivalent exact CTMC construction. Sampling only a mean residence time is insufficient.

`BASIN-004` A mean-rate construction that preserves mean occupancies or exit probabilities but loses time/exit correlation MUST be labeled `APPROXIMATE_MEAN_RATE`, not exact.

`BASIN-005` Assuming an exponential basin-exit time is permitted only when an explicit test accepts that approximation at a configured level; otherwise the phase-type distribution is required.

`BASIN-006` Exactness gate fails and serial KMC MUST be used when any transient/exit rate is missing or invalid, the basin can be nonabsorbing, detailed balance fails where required, the linear problem is ill-conditioned beyond tolerance, catalog status is incomplete, or concurrent non-basin events cannot be represented.

`BASIN-007` Fallback from accelerated to serial KMC is status `BASIN_DISABLED`, not a run failure. The reason and gate metrics MUST be recorded.

`BASIN-008` Exact-accelerator acceptance requires agreement with unaccelerated CTMC for joint exit-channel probabilities and exit-time distributions, not only mean time or diffusion coefficient.

## 10. Atom insertion, deletion, and topology scope

`ATOM-001` P0 MUST reject any callback or event changing atom count or species with `ATOM_COUNT_CHANGE_UNSUPPORTED`.

`ATOM-002` P0 MAY represent desorption-like motion only while all atoms remain in the fixed simulation state. Removing an atom into an implicit reservoir is forbidden.

`ATOM-003` Continuous surface reconstruction, concerted exchanges, defect migration, molecule dissociation, and rebonding are in scope when composition is fixed and the calculator supports them.

`ATOM-004` A future variable-composition extension MUST define reservoir identity, stoichiometric change, chemical potential/standard state, forward/reverse proposal multiplicity, volume/area convention, and local-detailed-balance test before it is selectable.

## 11. Determinism, RNG, and checkpoint

`DET-001` The normative RNG is counter-based `Philox4x32-10`. Checkpoint stores algorithm ID, key, counter, substream map, and consumed-block counts.

`DET-002` A uniform double is generated from a 53-bit integer (m\in[0,2^{53}-1]) as (u=(m+0.5)2^{-53}). Thus (0<u<1).

`DET-003` Saddle-search substreams MUST depend only on run seed, state ID, search class, and search index. Worker count, completion order, and retry timing MUST NOT change them.

`DET-004` Given deterministic callbacks, identical config/input/checkpoint MUST reproduce state IDs, event IDs, selected-event sequence, RNG counters, and status sequence. Geometry and time follow the parity tolerances in §14.

`CKPT-001` A checkpoint MUST contain schema/config/model digests; full current state; simulation time; step index; catalog; discovery statistics; RNG state; optional basin state; resource counters; log sequence; and completion/incomplete flags.

`CKPT-002` Checkpoint serialization is UTF-8 canonical JSON: sorted object keys, no duplicate keys, finite numbers only, and enough decimal digits for round-trip IEEE-754 binary64.

`CKPT-003` Each checkpoint MUST store SHA-256 for its canonical payload. Write behavior MUST be crash-safe: write a sibling temporary file, flush file content, atomically replace, then flush the parent directory.

`CKPT-004` Resume MUST verify every digest and payload hash before mutation. Failure returns `CHECKPOINT_CORRUPT` or `CHECKPOINT_INCOMPATIBLE`; best-effort partial resume is forbidden.

`CKPT-005` Checkpoint frequency is configurable by committed steps and wall time. Cancellation and resource exhaustion MUST attempt one checkpoint without changing kinetic state.

## 12. Input and resource contract

`IO-001` P0 accepts one UTF-8 JSON configuration with required top-level objects: `schema`, `system`, `calculator`, `relaxation`, `saddle_search`, `discovery`, `kinetics`, `resources`, and `output`.

`IO-002` Unknown keys in normative objects return `INVALID_INPUT`. Free-form user content is allowed only under `metadata` and has no behavioral effect.

`IO-003` Paths are resolved relative to the input file unless absolute. Input and output paths, encoding, and overwrite policy MUST be validated before calculator work starts.

`IO-004` Existing outputs MUST NOT be overwritten unless `overwrite=true`. Existing compatible checkpoint plus `resume=true` resumes; every other collision returns `OUTPUT_EXISTS`.

`RES-001` Limits MUST exist for wall time, total calculator evaluations, evaluations per relaxation, evaluations per saddle attempt, saddle attempts per state, catalog events, resident memory, output bytes, and callback timeout.

`RES-002` A limit hit returns `RESOURCE_LIMIT`, aborts the current uncommitted transaction, preserves last committed state/time/RNG, and attempts checkpoint.

`RES-003` Retries are allowed only for callbacks explicitly marked retryable and only within configured count/backoff. Retry MUST use identical request and RNG substream.

`RES-004` Resource counters and retry history MUST be in checkpoints and final summaries.

## 13. Exact status and error semantics

Every response contains one status, severity, message, context object, and causal status where applicable. Messages are descriptive; automation MUST branch on status.

| Status | Severity | Exact effect |
|---|---|---|
| `OK` | success | Requested transaction committed. |
| `DISCOVERY_CONVERGED_HEURISTIC` | success-with-qualification | Discovery stopping rule passed; completeness is not proven. |
| `DUPLICATE_EVENT` | candidate reject | Existing event retained; multiplicity/statistics updated; discovery continues. |
| `SADDLE_NOT_FOUND` | candidate reject | Attempt ended without candidate; does not count as redundant. |
| `INVALID_SADDLE` | candidate reject | First-order-saddle validation failed; candidate discarded. |
| `SADDLE_WRONG_BASIN` | candidate reject | Neither endpoint matches origin; candidate discarded. |
| `ENDPOINT_COLLAPSED` | candidate reject | Both endpoints match origin; candidate discarded. |
| `ENVIRONMENT_AMBIGUOUS` | recoverable | Generic reuse bypassed; fresh discovery required. |
| `BASIN_DISABLED` | recoverable | Accelerator gate failed; serial KMC used. |
| `DISCOVERY_INCOMPLETE` | pause/qualified | Budget ended before stop rule; strict mode pauses, exploratory mode permanently tags trajectory. |
| `RELAX_NOT_CONVERGED` | transaction fail | Candidate/state not committed; retry only by configured policy. |
| `EVENT_APPLICATION_FAILED` | transaction fail | KMC state/time/RNG unchanged. |
| `CALCULATOR_FAILURE` | transaction fail | Partial callback data discarded; retry only if marked retryable. |
| `NONFINITE_RESULT` | fatal | Run stops after checkpoint attempt. |
| `INVALID_INPUT` | fatal | No scientific work starts. |
| `SCHEMA_UNSUPPORTED` | fatal | No scientific work starts. |
| `INVALID_STATE` | fatal | Run stops; current transaction discarded. |
| `RATE_INVALID` | fatal in strict mode | Affected event disabled; strict run stops. |
| `DETAILED_BALANCE_VIOLATION` | fatal in strict mode | Affected pair disabled; strict run stops. |
| `CATALOG_CONFLICT` | fatal | Catalog mutation rolled back. |
| `CATALOG_INCOMPATIBLE` | fatal | Reuse/resume refused. |
| `ATOM_COUNT_CHANGE_UNSUPPORTED` | fatal | Transaction rolled back. |
| `NO_ENABLED_EVENT` | terminal-success if requested, else fatal | State is absorbing; summary/checkpoint written. |
| `RESOURCE_LIMIT` | pause | Uncommitted work discarded; checkpoint attempted. |
| `OUTPUT_EXISTS` | fatal | No output overwritten. |
| `CHECKPOINT_CORRUPT` | fatal | Resume refused. |
| `CHECKPOINT_INCOMPATIBLE` | fatal | Resume refused. |
| `CANCELLED` | pause | Uncommitted work discarded; checkpoint attempted. |
| `INTERNAL_ERROR` | fatal | Uncommitted work discarded; diagnostic and checkpoint attempt required. |

`ERR-001` Candidate rejects MUST NOT alter committed state, KMC time, trajectory RNG, or catalog except explicit duplicate statistics.

`ERR-002` Process exit codes are exact. Code `0`: `OK`, `DISCOVERY_CONVERGED_HEURISTIC`, exploratory continuation after `DISCOVERY_INCOMPLETE`, or requested absorbing `NO_ENABLED_EVENT`. Code `64`: `INVALID_INPUT`, `SCHEMA_UNSUPPORTED`, `OUTPUT_EXISTS`. Code `65`: exhausted `RELAX_NOT_CONVERGED`, `EVENT_APPLICATION_FAILED`, `NONFINITE_RESULT`, `INVALID_STATE`, `RATE_INVALID`, `DETAILED_BALANCE_VIOLATION`, `CATALOG_CONFLICT`, `CATALOG_INCOMPATIBLE`, `ATOM_COUNT_CHANGE_UNSUPPORTED`, or unrequested `NO_ENABLED_EVENT`. Code `69`: exhausted `CALCULATOR_FAILURE`. Code `70`: `INTERNAL_ERROR`. Code `74`: `CHECKPOINT_CORRUPT`, `CHECKPOINT_INCOMPATIBLE`, or fatal checkpoint/output I/O failure. Code `75`: strict `DISCOVERY_INCOMPLETE`, `RESOURCE_LIMIT`, or `CANCELLED`. Candidate-reject and recoverable statuses do not end the process by themselves.

`ERR-003` Error context MUST include component, requirement ID, state ID if known, search/event ID if known, and retryability.

`ERR-004` If several failures are present, the returned status MUST be the earliest causal failure in transaction order. Cleanup/checkpoint failures are appended as causal context and change the process exit code to `74` only when no valid last-committed checkpoint remains.

## 14. Python/Rust parity

`PAR-001` Python and Rust implementations MUST accept the same canonical schema, units, status set, event/checkpoint records, and requirement IDs.

`PAR-002` They MUST pass common golden files. Internal algorithms and data structures need not match.

`PAR-003` With deterministic mock callbacks, event/state IDs, selected events, RNG words/counters, checkpoint canonical JSON, and status sequence MUST be byte-identical.

`PAR-004` Energies, forces, barriers, rates, and geometries MUST agree to configured scientific tolerances. For the analytic conformance suite defaults are relative `1e-13` plus absolute `1e-14`; geometry absolute `1e-12 Å`.

`PAR-005` Because system `log` implementations can differ, (\Delta t) parity is relative `5e-15` plus absolute `1e-18 s`. A stricter bitwise-time claim requires the same correctly rounded logarithm implementation and MUST be declared separately.

`PAR-006` Parallel discovery MUST yield the same committed catalog as serial discovery after sorting by search ID and transactional deduplication.

## 15. Validation and benchmark matrix

No implementation benchmarks were run for this specification. These are acceptance obligations.

| ID | Fixture | Required assertion | Priority |
|---|---|---|---|
| `VAL-001` | One-state, two-rate analytic CTMC | Event frequencies and waiting-time CDF match analytic values; deterministic golden stream passes. | P0 |
| `VAL-002` | Two-state reversible CTMC | Stationary populations and rate ratio pass `RATE-005`. | P0 |
| `VAL-003` | Three-state network with parallel saddles | Parallel rates add; distinct mechanisms are not deduplicated. | P0 |
| `VAL-004` | Periodic analytic surface potential | Minimum, saddle, endpoints, barriers, and wrapped state identity match analytic values. | P0 |
| `VAL-005` | Two-dimensional multiwell potential | Multiple exits found; wrong-basin and collapsed endpoints rejected. | P0 |
| `VAL-006` | Same structure under translation/image/permutation | State and environment identities remain invariant. | P0 |
| `VAL-007` | Neighbor at ambiguity boundary | `ENVIRONMENT_AMBIGUOUS`; fresh search used. | P0 |
| `VAL-008` | Same endpoints, two saddle geometries | Two events retained. | P0 |
| `VAL-009` | Incomplete search budget | Strict mode stops with `DISCOVERY_INCOMPLETE`; exploratory tag persists after resume. | P0 |
| `VAL-010` | Nonfinite/malformed calculator responses | Exact statuses; no partial commit. | P0 |
| `VAL-011` | Failure during event application | State/time/RNG unchanged. | P0 |
| `VAL-012` | Forced cancellation/resource limit | Last checkpoint resumes to same next selected event. | P0 |
| `VAL-013` | Python/Rust golden corpus | All `PAR-*` requirements pass. | P0 |
| `VAL-014` | Fixed-composition reconstructed surface toy model | Concerted reconstruction is discoverable without lattice sites. | P0 |
| `VAL-015` | Finite absorbing CTMC basin | Joint exit-channel distribution and exit-time CDF match unaccelerated CTMC. | P1 exact basin |
| `VAL-016` | Nonexponential phase-type basin | Mean-only/exponential shortcut fails exact gate; exact sampler passes. | P1 exact basin |
| `VAL-017` | Missing exit or singular basin | `BASIN_DISABLED`; serial result preserved. | P1 exact basin |
| `VAL-018` | Published atomistic case chosen by authors | Mechanism classes and barriers compared with declared calculator/tolerances; no universal numeric target is imposed across potentials. | Publication |

`VAL-019` Statistical tests MUST state sample count, seed set, confidence interval, test statistic, and acceptance threshold before execution.

`VAL-020` Performance claims MUST report calculator cost separately from orchestration cost and must compare identical scientific tolerances and discovery budgets.

## 16. P0 acceptance matrix

| Gate | Requirements | Pass condition |
|---|---|---|
| Schema/state | `SCOPE-001..004`, `MATH-001..003`, `STATE-001..008`, `IO-001..004` | All mandatory fields validated; off-lattice fixture passes. |
| Callbacks | `CALC-001..006`, `RELAX-001..004`, `SADDLE-001..007` | Transaction and validation negative tests pass. |
| Discovery/catalog | `DISC-001..008`, `EVENT-001..004`, `CAT-001..007`, `ENV-001..008` | No unvalidated event becomes selectable; incomplete status preserved. |
| Kinetics | `RATE-001..008`, `KMC-001..006` | Analytic CTMC and detailed-balance tests pass. |
| Scope safety | `ATOM-001..004`, `BASIN-001` | Variable composition rejected; basin disabled by default. |
| Reproducibility | `DET-001..004`, `CKPT-001..005`, `RES-001..004`, `ERR-001..004` | Crash/resume and exact status tests pass. |
| Language parity | `PAR-001..006` | `VAL-013` passes with no skipped mandatory case. |
| Scientific minimum | `VAL-001..014`, `VAL-019..020` | All P0 tests pass; limitations reported. |

`ACCEPT-001` P0 conformance requires every P0 `MUST` and every `VAL-001..014` test. P1 and publication rows MAY be pending and MUST be labeled so.

`ACCEPT-002` A publishable software claim additionally requires public schemas/golden fixtures, at least one documented calculator adapter, one minimizer, one saddle solver, Python/Rust results, provenance, and a limitations section stating catalog incompleteness.

## 17. Evidence-backed and speculative extensions

`EXT-001` Evidence-backed P1: local-environment reuse with mandatory specific reconvergence; calibrated discovery stopping; harmonic prefactors; exact absorbing-CTMC superbasin exits; fixed-composition reconstruction/diffusion applications.

`EXT-002` Evidence-backed but deliberately deferred: hybrid MD for externally scheduled deposition/fast relaxation and predefined-table variable-composition events. They add a different event semantics and are not needed for the OTF core.

`EXT-003` Speculative until separately validated: learned saddle-proposal directions, uncertainty-calibrated fingerprints across calculator families, adaptive environment radii, grand-canonical OTF discovery, and automatic catalyst reconstruction/reaction-network co-learning.

`EXT-004` Time Warp, GPU execution, and tau-leap remain excluded unless a benchmark first proves serial P0 is the dominant bottleneck and an exactness analysis covers altered event ordering/time semantics.

## 18. Clean-room and license boundary

`LIC-001` Implementers MUST work from this specification, cited papers, and official public behavior/API documents. Prohibited implementation source, source history, copied pseudocode, tests, symbols, comments, and file layout MUST NOT be consulted or reproduced.

`LIC-002` Every implementation file MUST carry author/provenance metadata identifying this specification revision and independently authored status.

`LIC-003` Dependency names, versions, licenses, and linking mode MUST be recorded. License compatibility and patent clearance are separate release gates; paper publication does not grant a software or patent license.

`LIC-004` The literature review reports article licenses where verified. It does not infer software licenses from article licenses.

`LIC-005` A permissive project license such as `Apache-2.0 OR MIT` is a recommendation only and requires repository-owner approval. No license is assigned by this specification.

`LIC-006` Clean-room review MUST compare observable requirements and independently authored tests, never competitor source structure.

## 19. Final design decision

The smallest publishable v1 is P0, not the full roadmap. Its scientific claim is: a deterministic, calculator-agnostic, fixed-composition, serial off-lattice/on-the-fly KMC engine with validated saddle-to-minimum events, reversible rates, explicit discovery uncertainty, portable checkpoints, and Python/Rust behavioral parity. Exact basin acceleration is an optional later module. Variable composition and parallel-time acceleration require new specifications.

## 20. Spec digest

Hash rule: SHA-256 of the exact UTF-8 bytes through the newline immediately before this heading; this digest section is excluded.

Spec SHA-256: `8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84`
