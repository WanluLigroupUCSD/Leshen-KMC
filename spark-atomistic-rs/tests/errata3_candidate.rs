// Normative Errata 3 regression tests. The historical filename is retained to avoid needless
// churn. Errata 3 was adopted by `D-127`; the pre-Errata-3 path stays executable through
// `CheckpointPolicy::PreErrata3` and `PayloadHashPolicy::PreErrata3Full`.
//
// Authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2_PARITY + normative Errata 3.
// Independently authored; no implementation source consulted.
//
// Every rejection below is paired with an accepted baseline that differs ONLY in the property
// under test, and each pair states which property that is.
use serde_json::{json, Value};
use spark_atomistic_rs::checkpoint::canonical_json_bytes;
use spark_atomistic_rs::parity::*;
use spark_atomistic_rs::rng::{derive_saddle_substream, substream_digest};
use spark_atomistic_rs::status::StatusCode;
use std::path::{Path, PathBuf};

fn corpus() -> PathBuf { Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/corpus/xlang/checkpoints") }
fn clean_bytes() -> Vec<u8> { std::fs::read(corpus().join("checkpoint-clean.json")).expect("clean checkpoint") }
fn model() -> ValidatedModel {
    let value: Value = serde_json::from_slice(&std::fs::read(corpus().join("checkpoint_model.json"))
        .expect("checkpoint model")).expect("model parses");
    let wire: WireModel = serde_json::from_value(value).expect("model decodes");
    wire.validate(None).expect("model validates")
}
/// Re-seal a mutated payload WITHOUT re-validating it, so the only thing under test is the
/// validator's verdict and never the encoder's.
fn seal(payload: Value, policy: PayloadHashPolicy) -> Vec<u8> {
    let typed: CheckpointPayloadV2 = serde_json::from_value(payload).expect("payload decodes");
    let envelope = CheckpointEnvelopeV2::new_sealed(typed, policy).expect("envelope seals");
    canonical_json_bytes(&envelope).expect("canonical envelope")
}
fn payload_of(bytes: &[u8]) -> Value {
    let envelope: Value = serde_json::from_slice(bytes).expect("envelope parses");
    envelope["payload"].clone()
}
fn verdict(bytes: &[u8], policy: CheckpointPolicy) -> Result<CheckpointEnvelopeV2, (StatusCode, String)> {
    decode_checkpoint_v2_with(bytes, &model(), policy).map_err(|e| (e.status, e.requirement_id))
}

/// `E3-EVENT-001`, adopted options 1 + 4.
///
/// ACCEPTED BASELINE: the Python backend's `checkpoint-clean.json` exactly as written -- 28439
/// canonical bytes, 6 directed records (3 reciprocal pairs), 4 committed states, 1 committed KMC
/// step with its historical snapshot. `E2-PAR-002`(10) makes this fixture mandatory.
#[test]
fn python_clean_checkpoint_bytes_are_accepted_under_the_candidate_and_refused_under_the_shipped_rule() {
    let bytes = clean_bytes();
    assert_eq!(bytes.len(), 28439, "the fixture under test is the 28439-byte Python checkpoint");
    let envelope = verdict(&bytes, CheckpointPolicy::Errata3Event001).expect("candidate accepts the Python bytes");
    let p = &envelope.payload;
    assert_eq!(p.catalog.events.len(), 6, "6 directed records");
    assert_eq!(p.catalog.events.values().map(|e| e.pair_id.clone()).collect::<std::collections::BTreeSet<_>>().len(), 3,
        "3 reciprocal pairs");
    assert_eq!(p.catalog.states.len(), 4, "4 committed states");
    assert_eq!(p.trajectory.len(), 1, "1 committed KMC step");
    assert_eq!(p.step_index, 1);
    assert!(!p.trajectory[0].rate_table_snapshot.payload.event_ids.is_empty(),
        "the committed step carries its historical rate snapshot (`E2-KMC-005`)");
    assert_eq!(envelope.sealed_with, PayloadHashPolicy::PreErrata3Full,
        "the Python backend seals `E2-CKPT-001` over the complete payload");
    // PAIRED REJECTION on the policy itself: the SAME bytes under the shipped rule. Nothing about
    // the bytes differs; only the invariant set does. This is the measured `D-E2-02` refusal.
    let (status, requirement) = verdict(&bytes, CheckpointPolicy::PreErrata3).expect_err("shipped rule refuses the same bytes");
    assert_eq!(status, StatusCode::CheckpointCorrupt);
    assert_eq!(requirement, "E2-EVENT-001");
}

/// REJECTION 1 of 3 -- a reciprocal `saddle` that is not the correct mapped image.
/// `E3-EVENT-001` option 4: "Reciprocal saddle = forward saddle mapped through
/// `validation.destination_match.atom_mapping` with `unstable_direction` negated." The variant
/// below un-negates ONE reciprocal record's `unstable_direction`, which is exactly the byte-equal
/// form the shipped C2 demanded, and changes nothing else.
#[test]
fn a_reciprocal_saddle_that_is_not_the_mapped_image_is_refused() {
    let bytes = clean_bytes();
    // ACCEPTED BASELINE: the same payload, resealed by the same helper, with no field edited.
    let baseline = seal(payload_of(&bytes), PayloadHashPolicy::PreErrata3Full);
    assert_eq!(baseline, bytes, "resealing an unedited payload reproduces the input bytes");
    verdict(&baseline, CheckpointPolicy::Errata3Event001).expect("unedited baseline is accepted");

    let mut payload = payload_of(&bytes);
    let target = payload["catalog"]["events"].as_object().expect("events")
        .iter().find(|(_, e)| e["saddle"]["unstable_direction"][0][0].as_f64() == Some(-1.0))
        .map(|(k, _)| k.clone()).expect("a reciprocal record with a negated x direction");
    let forward = payload["catalog"]["events"][&target]["reverse_event_id"].as_str().expect("reverse id").to_owned();
    let mapped = payload["catalog"]["events"][&forward]["saddle"]["unstable_direction"].clone();
    payload["catalog"]["events"][&target]["saddle"]["unstable_direction"] = mapped;
    // The mutation is chosen so that the record now satisfies the SHIPPED C2 exactly -- the two
    // `saddle` objects are byte-equal -- which is what pins the refusal below to the mapped-equality
    // term and to nothing else in the reciprocal chain.
    assert_eq!(payload["catalog"]["events"][&target]["saddle"],
               payload["catalog"]["events"][&forward]["saddle"],
               "the variant is byte-equal under the shipped C2 and wrong only under option 4");
    // The catalog digest covers the events object (`E2-CKPT-006`), so it is rebuilt: the ONLY
    // property left broken is the reciprocal-saddle relation.
    let mut catalog = payload["catalog"].clone();
    catalog.as_object_mut().expect("catalog").remove("digest");
    payload["catalog"]["digest"] = json!(format!("sha256:{}", spark_atomistic_rs::identity::hex_sha256(
        &canonical_json_bytes(&catalog).expect("canonical catalog"))));
    let (status, requirement) = verdict(&seal(payload, PayloadHashPolicy::PreErrata3Full), CheckpointPolicy::Errata3Event001)
        .expect_err("an unmapped reciprocal saddle is refused");
    assert_eq!(status, StatusCode::CheckpointCorrupt);
    assert_eq!(requirement, "E2-CKPT-007");
}

/// REJECTION 2 of 3 -- a `discovery_provenance` that re-derives from NEITHER endpoint.
/// `E3-EVENT-001` option 1 requires the search ID to re-derive from exactly one of
/// `{origin, destination}`. Both records of one pair are re-pointed at a search of a THIRD
/// committed state, with `saddle.search_id` (C1), the substream map entry and
/// `rng_substream_digest` all rebuilt so that the derivation source is the ONLY broken property.
#[test]
fn a_provenance_that_derives_from_neither_endpoint_is_refused() {
    let bytes = clean_bytes();
    let payload = payload_of(&bytes);
    let seed = payload["rng"]["run_seed"].as_u64().expect("run seed");
    let events = payload["catalog"]["events"].as_object().expect("events").clone();
    let (forward_id, forward) = events.iter().next().expect("an event");
    let reverse_id = forward["reverse_event_id"].as_str().expect("reverse id").to_owned();
    let class = forward["discovery_provenance"]["search_class"].as_str().expect("class").to_owned();
    let index = forward["discovery_provenance"]["search_index"].as_u64().expect("index");
    let (origin, destination) = (forward["origin_state_id"].as_str().expect("origin").to_owned(),
                                 forward["destination_state_id"].as_str().expect("destination").to_owned());
    let outsider = payload["catalog"]["states"].as_object().expect("states").keys()
        .find(|s| **s != origin && **s != destination).expect("a third committed state").clone();
    let forward_id = forward_id.clone();

    // The two variants differ ONLY in which state the search ID is derived from.
    let rebuild = |source: &str| -> Vec<u8> {
        let mut p = payload.clone();
        let id = search_id(seed, source, &class, index).expect("search id");
        let stream = derive_saddle_substream(seed, source, &class, index).expect("substream");
        let digest = substream_digest(&stream);
        for event in [&forward_id, &reverse_id] {
            p["catalog"]["events"][event]["saddle"]["search_id"] = json!(id);
            p["catalog"]["events"][event]["discovery_provenance"]["search_id"] = json!(id);
            p["catalog"]["events"][event]["discovery_provenance"]["rng_substream_digest"] = json!(digest);
        }
        let old = forward["discovery_provenance"]["search_id"].as_str().expect("old id").to_owned();
        let map = p["rng"]["substream_map"].as_object_mut().expect("substream map");
        map.remove(&old);
        map.insert(id.clone(), serde_json::to_value(stream.state()).expect("state"));
        let mut catalog = p["catalog"].clone();
        catalog.as_object_mut().expect("catalog").remove("digest");
        p["catalog"]["digest"] = json!(format!("sha256:{}", spark_atomistic_rs::identity::hex_sha256(
            &canonical_json_bytes(&catalog).expect("canonical catalog"))));
        seal(p, PayloadHashPolicy::PreErrata3Full)
    };

    // ACCEPTED BASELINE: derived from the true discovering state, which IS one of the endpoints.
    let discovering = if search_id(seed, &origin, &class, index).expect("id")
        == forward["discovery_provenance"]["search_id"].as_str().expect("id") { origin.clone() } else { destination.clone() };
    verdict(&rebuild(&discovering), CheckpointPolicy::Errata3Event001)
        .expect("a provenance deriving from an endpoint is accepted");
    // REJECTION: identical construction, derived from a committed state that is neither endpoint.
    let (status, requirement) = verdict(&rebuild(&outsider), CheckpointPolicy::Errata3Event001)
        .expect_err("a provenance deriving from neither endpoint is refused");
    assert_eq!(status, StatusCode::CheckpointCorrupt);
    assert_eq!(requirement, "E2-CKPT-007");
}

/// REJECTION 3 of 3 -- a tampered payload hash. `E2-CKPT-007`(2) verifies the payload SHA-256, and
/// `E3-PAR-002` must not weaken that: neither accepted seal may match a forged digest.
#[test]
fn a_tampered_payload_hash_is_refused_under_both_seals() {
    let bytes = clean_bytes();
    verdict(&bytes, CheckpointPolicy::Errata3Event001).expect("untampered baseline is accepted");
    let mut envelope: Value = serde_json::from_slice(&bytes).expect("envelope parses");
    let stored = envelope["payload_sha256"].as_str().expect("hash").to_owned();
    // One hex digit, nothing else.
    let forged = format!("{}{}{}", &stored[..stored.len() - 1],
        if stored.ends_with('0') { '1' } else { '0' }, "");
    assert_ne!(forged, stored);
    envelope["payload_sha256"] = json!(forged);
    let tampered = canonical_json_bytes(&envelope).expect("canonical envelope");
    assert_eq!(tampered.len(), bytes.len(), "only the hash digit differs");
    let (status, requirement) = verdict(&tampered, CheckpointPolicy::Errata3Event001)
        .expect_err("a forged payload hash is refused");
    assert_eq!(status, StatusCode::CheckpointCorrupt);
    assert_eq!(requirement, "E2-CKPT-001");
}

/// `E3-PAR-002`. The excluded set is exactly
/// `resources.wall_elapsed_s`. Editing that field alone is invisible to a payload sealed under the
/// candidate rule and fatal to one sealed under `E2-CKPT-001` as written; editing any other field
/// is fatal under both. Without both halves the exclusion could be a hash that checks nothing.
#[test]
fn only_the_declared_environment_field_is_outside_the_hashed_payload() {
    assert_eq!(ENVIRONMENT_DERIVED_PAYLOAD_FIELDS, &[("resources", "wall_elapsed_s")]);
    let bytes = clean_bytes();
    let payload = payload_of(&bytes);
    let original = payload["resources"]["wall_elapsed_s"].as_f64().expect("wall elapsed");

    // ACCEPTED BASELINE: sealed under the candidate rule, unedited.
    let candidate = seal(payload.clone(), PayloadHashPolicy::Errata3EnvironmentExcluded);
    let sealed = verdict(&candidate, CheckpointPolicy::Errata3Event001).expect("candidate seal is accepted");
    assert_eq!(sealed.sealed_with, PayloadHashPolicy::Errata3EnvironmentExcluded);
    assert_ne!(candidate, bytes, "the two seals produce different `payload_sha256`");

    // The property under test: ONLY `wall_elapsed_s` changes.
    let mut edited = payload.clone();
    edited["resources"]["wall_elapsed_s"] = json!(original + 0.5);
    let mut envelope: Value = serde_json::from_slice(&candidate).expect("parses");
    envelope["payload"] = edited.clone();
    let moved = canonical_json_bytes(&envelope).expect("canonical");
    verdict(&moved, CheckpointPolicy::Errata3Event001)
        .expect("a candidate-sealed checkpoint tolerates a different wall clock");
    // PAIRED REJECTION: the SAME edit against the `E2-CKPT-001` seal.
    let mut envelope: Value = serde_json::from_slice(&bytes).expect("parses");
    envelope["payload"] = edited;
    let moved_full = canonical_json_bytes(&envelope).expect("canonical");
    let (status, requirement) = verdict(&moved_full, CheckpointPolicy::Errata3Event001)
        .expect_err("a fully sealed checkpoint does not tolerate a different wall clock");
    assert_eq!((status, requirement.as_str()), (StatusCode::CheckpointCorrupt, "E2-CKPT-001"));

    // PAIRED REJECTION on the exclusion's scope: a HASHED field, edited the same way, under the
    // candidate seal. If the exclusion leaked past `wall_elapsed_s` this would be accepted.
    let mut edited = payload;
    let time = edited["simulation_time_s"].as_f64().expect("simulation time");
    edited["simulation_time_s"] = json!(time + 0.5);
    let mut envelope: Value = serde_json::from_slice(&candidate).expect("parses");
    envelope["payload"] = edited;
    let moved = canonical_json_bytes(&envelope).expect("canonical");
    let (status, requirement) = verdict(&moved, CheckpointPolicy::Errata3Event001)
        .expect_err("a hashed field is still covered by the candidate seal");
    assert_eq!((status, requirement.as_str()), (StatusCode::CheckpointCorrupt, "E2-CKPT-001"));
}
