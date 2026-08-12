# Cross-language parity corpus

One canonical corpus, consumed byte-identically by both backends. `E2-PAR-001` requires exactly
that; `E2-PAR-003` requires the emitted canonical JSON to be identical; `E2-PAR-005` makes execution
conformance conditional on every mandatory fixture passing in *every* implementation, which no
single-backend suite can establish.

## Layout

| Path | Contents |
|---|---|
| `manifest.json` | case index: request cases with their exact byte length, probe IDs, checkpoint cases |
| `requests/*.req` | 46 raw request byte strings fed to the public wire surface. Some are deliberately not valid JSON or not valid UTF-8, which is why they are files rather than JSON string literals. |
| `probes.json` | 30 probe inputs for artifacts that are not on the public wire: canonical numbers, digests, state identities, Philox blocks and streams, search IDs, discovery class selection, event ID algebra, rates, rate-table snapshots, KMC selection, and the full status table |
| `checkpoints/` | `checkpoint_model.json` plus 7 checkpoint envelopes: clean, zero-step control, corrupt hash, incompatible digest, recursive corruption, cancellation, resource limit |

`checkpoints/checkpoint-clean.json` is 28439 canonical bytes holding six directed records (three
reciprocal pairs), four committed states, one committed KMC step and its historical rate snapshot.
It is the control that separates defect `D-E2-02` from everything else: the paired
`checkpoint-zero-step.json` uses the same model with an empty catalog, and both backends accept it
and answer byte-identically.

## Drivers

```text
cd ../../..                                            # spark-atomistic-rs
SPARK_XLANG_OUT=DIR CARGO_TARGET_DIR=TMP cargo test --test xlang_emit
cd ../spark-atomistic
python3 tests/xlang_harness.py emit    --out DIR
python3 tests/xlang_harness.py compare --out DIR       # nonzero exit on any divergence
```

Each backend writes 85 canonical artifacts to `DIR/<backend>/<case>.out`. The comparator diffs them
byte-for-byte and prints the first differing offset with both surrounding strings. `DIR/rust/
_roundtrip.json` records the artifacts the Rust backend cannot re-read (defect `D-E2-03`).

## Rules the corpus obeys

- Every rejection case is paired with an accepted baseline that differs only in the property under
  test, so a rejection can never be credited to an unrelated malformation.
- Nothing is normalised before comparison. A divergence is reported, quantified, and attributed to a
  requirement ID or recorded as unpinned by the erratum.
- Two probes are marked SHARED-ALGEBRA in `xlang_harness.py` (`event_ids`, `digests`): the Python
  backend exposes no standalone function for them, so the harness composes the spec-quoted payload
  and only the backend's canonical encoder and SHA-256 are under test. Both currently agree, and the
  authoritative event-record comparison is the checkpoint case, which uses real records from both
  sides.

The per-fixture verdicts are recorded in `../e2_parity_manifest.json` under
`cross_language_execution` and are checked against the mandatory fixture list by
`parity_manifest_matches_this_suite`.
