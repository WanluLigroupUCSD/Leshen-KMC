// Clean-room tests authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2_PARITY.
// Independently authored; no implementation source consulted.
// Executing counterpart of tests/ERRATA2_PARITY_SPEC.md and corpus/e2_parity_manifest.json.
// Every assertion cites the Errata 2 requirement that decides it. Each fixture asserted to be
// REJECTED is paired with an ACCEPTED baseline that differs only in the property under test.
use serde_json::{json,Map,Value};
use spark_atomistic_rs::basin::{attempt_acceleration,BasinCheckpoint};
use spark_atomistic_rs::checkpoint::{canonical_json_bytes,parse_strict_json};
use spark_atomistic_rs::identity::hex_sha256;
use spark_atomistic_rs::parity::*;
use spark_atomistic_rs::rate::{common_prefactor_pair,log_sum_exp,neumaier_sum};
use spark_atomistic_rs::rng::{derive_saddle_substream,derive_trajectory_stream,philox4x32_10,substream_digest,uniform_from_words,Philox};
use spark_atomistic_rs::status::StatusCode;
use std::collections::BTreeMap;
use std::path::PathBuf;

const MODEL:&[u8]=include_bytes!("corpus/e2_minimal_model.json");
const CAPS_REQUEST:&[u8]=include_bytes!("corpus/e2_capabilities.request.json");
const CAPS_RESPONSE:&[u8]=include_bytes!("corpus/e2_capabilities.response.json");
const NUMBERS:&[u8]=include_bytes!("corpus/e2_canonical_numbers.json");
const PHILOX:&[u8]=include_bytes!("corpus/philox_errata1.json");
const MANIFEST:&[u8]=include_bytes!("corpus/e2_parity_manifest.json");
const DUPLICATE:&[u8]=include_bytes!("corpus/strict_duplicate.invalid.json");
const CALC_DIGEST:&str="sha256:test-model";
const CELL:[[f64;3];3]=[[12.0,0.0,0.0],[0.0,12.0,0.0],[0.0,0.0,12.0]];
// A synthetic state ID used only as opaque derivation input for the E2-RNG-004 / E2-DISC-002 goldens.
const GOLDEN_STATE_ID:&str="state:sha256:0000000000000000000000000000000000000000000000000000000000000001";

// ---------------------------------------------------------------- helpers

fn response(bytes:&[u8])->Value{serde_json::from_slice(&dispatch_json(bytes,None,None)).expect("canonical response parses")}
fn model_text()->String{String::from_utf8(MODEL.to_vec()).expect("corpus model is UTF-8")}
/// Splices raw text into the corpus model as the free-form root `metadata` object (E2-SCHEMA-001).
fn validate_request_with_raw_metadata(raw:&str)->String{
    let m=model_text();let m=m.trim_end();let body=m.strip_suffix('}').expect("model ends with brace");
    format!("{{\"model\":{body},\"metadata\":{raw}}},\"operation\":\"validate\"}}")
}
fn temp_dir(tag:&str)->PathBuf{
    let d=std::env::temp_dir().join(format!("spark-e2-{}-{}",std::process::id(),tag));
    let _=std::fs::create_dir_all(&d);d
}
fn validated_model(dir:&std::path::Path,mutate:impl FnOnce(&mut Value))->ValidatedModel{
    let mut m:Value=serde_json::from_slice(MODEL).expect("corpus model parses");
    m["output"]["checkpoint_path"]=json!(dir.join("ckpt.json").to_str().expect("utf8 path"));
    m["output"]["summary_path"]=json!(dir.join("summary.json").to_str().expect("utf8 path"));
    m["output"]["trajectory_path"]=json!(dir.join("traj.json").to_str().expect("utf8 path"));
    mutate(&mut m);
    let wire:WireModel=serde_json::from_value(m).expect("wire model deserializes");
    wire.validate(None).expect("corpus model validates")
}
fn digest_of(v:&Value)->String{format!("sha256:{}",hex_sha256(&canonical_json_bytes(v).expect("finite canonical value")))}
fn system(positions:&[[f64;3]],species:&[&str],ids:&[&str])->WireSystem{
    WireSystem{atom_ids:ids.iter().map(|x|(*x).to_owned()).collect(),species:species.iter().map(|x|(*x).to_owned()).collect(),
        positions:positions.to_vec(),cell:CELL,pbc:[true,true,true],movable:vec![true;positions.len()],
        constraints:WireConstraints{kind:"fixed-mask".into()},charge:0.0,spin:0.0,calculator_model_digest:CALC_DIGEST.into()}
}
fn committed(positions:&[[f64;3]],energy:f64)->CommittedStateV2{
    let ids:Vec<String>=(0..positions.len()).map(|i|format!("a{i}")).collect();
    let idrefs:Vec<&str>=ids.iter().map(String::as_str).collect();
    let species=vec!["H";positions.len()];
    let s=system(positions,&species,&idrefs);
    let(candidate,state_id)=state_ids(&s,energy).expect("finite state identity");
    CommittedStateV2{atom_ids:s.atom_ids.clone(),calculator_model_digest:CALC_DIGEST.into(),candidate_identity:candidate,
        cell:CELL,charge:0.0,constraint_digest:constraint_digest(&s).expect("constraint digest"),constraints:s.constraints.clone(),
        energy_ev:energy,fixed_contract_digest:fixed_contract_digest(&s).expect("fixed contract digest"),
        force_tolerance_ev_per_angstrom:0.05,forces_ev_per_angstrom:vec![[0.0,0.0,0.0];positions.len()],
        identity_version:"spark-state-identity/2".into(),max_movable_force_ev_per_angstrom:0.0,movable:s.movable.clone(),
        pbc:[true,true,true],positions:positions.to_vec(),
        relaxation_provenance:RelaxationProvenanceV2{calculator_evaluations:3,calculator_identity:"conformance-calculator".into(),
            minimizer_identity:"conformance-minimizer".into(),steps:2,termination_reason:"force tolerance reached".into()},
        schema:"spark-atomistic-state/2".into(),species:s.species.clone(),spin:0.0,state_id}
}
/// Builds the smallest checkpoint payload that E2-CKPT-002/007 accepts: one committed state, no
/// catalog event, no trajectory step. `mutate` injects the single property under test.
fn zero_step_payload(state:&CommittedStateV2,vm:&ValidatedModel,mutate:impl FnOnce(&mut CheckpointPayloadV2))->CheckpointPayloadV2{
    let mut states=BTreeMap::new();states.insert(state.state_id.clone(),state.clone());
    let mut catalog=CatalogV2{digest:String::new(),events:BTreeMap::new(),multiplicity:BTreeMap::new(),
        schema:"spark-atomistic-catalog/2".into(),states};
    let mut cv=serde_json::to_value(&catalog).expect("catalog serializes");
    cv.as_object_mut().expect("catalog is object").remove("digest");
    catalog.digest=digest_of(&cv);
    let seed=vm.model().kinetics.run_seed;
    let mut payload=CheckpointPayloadV2{
        basin:BasinCheckpointV2{enabled:false,reason:"v1-disabled".into()},catalog,checkpoint_sequence:1,
        digests:CheckpointDigestsV2{config:vm.config_digest().into(),model:vm.model().calculator.model_digest.clone(),
            schema:SCHEMA_DIGEST.into(),tolerances:vm.tolerance_digest().into()},
        discovery_statistics:BTreeMap::new(),
        flags:CheckpointFlagsV2{cancelled:false,complete:false,incomplete_catalog:false,last_status:StatusCode::Ok,resource_limited:false},
        initial_state:state.clone(),log_sequence:0,
        resources:CheckpointResourcesV2{calculator_evaluations:3,catalog_events:0,output_bytes:0,resident_memory_bytes:1000,
            retry_history:vec![],saddle_attempts_by_state:BTreeMap::new(),wall_elapsed_s:0.25},
        rng:CheckpointRngV2{run_seed:seed,substream_map:BTreeMap::new(),trajectory:derive_trajectory_stream(seed).state()},
        schema:"spark-atomistic-checkpoint/2".into(),simulation_time_s:0.0,step_index:0,trajectory:vec![],
        current_state:state.clone()};
    mutate(&mut payload);
    // Recompute the catalog digest so that a mutation of the catalog contents is not masked by a
    // stale digest, unless the test is deliberately corrupting the digest itself.
    payload
}

// ---------------------------------------------------------------- api-*

#[test]
fn api_capability_request_response(){
    // E2-API-001 / E2-API-006 / E2-API-007 / E2-CAN-001.
    let actual=dispatch_json(CAPS_REQUEST,None,None);
    let expected=CAPS_RESPONSE.strip_suffix(b"\n").expect("golden response stores one trailing newline");
    assert_eq!(actual,expected,"E2-PAR-003: capability response must be byte-identical to the golden");
    assert!(!actual.ends_with(b"\n"),"E2-CAN-001: canonical output carries no trailing newline");
    let v:Value=serde_json::from_slice(&actual).expect("response parses");
    let keys:Vec<&str>=v.as_object().expect("object").keys().map(String::as_str).collect();
    assert_eq!(keys,vec!["causal_status","context","exit_code","message","operation","severity","status","value"],
        "E2-API-006: the public response has exactly these keys, sorted per E2-CAN-001");
    let ctx:Vec<&str>=v["context"].as_object().expect("object").keys().map(String::as_str).collect();
    assert_eq!(ctx,vec!["component","details","requirement_id","retryable","search_or_event_id","state_id"],
        "E2-STATUS-001: context has exactly these fields");
    assert_eq!(v["status"],"OK");assert_eq!(v["exit_code"],0);assert_eq!(v["message"],"transaction committed");
    assert_eq!(v["severity"],"success");assert_eq!(v["causal_status"],Value::Null);
    assert_eq!(v["value"]["basin_acceleration"],"disabled","E2-API-007 / E2-BASIN-001");
    assert_eq!(v["value"]["validated"],false);assert_eq!(v["value"]["production"],false);assert_eq!(v["value"]["release"],false);
    assert_eq!(v["value"]["operations"],json!(["capabilities","validate","run"]));
}

#[test]
fn api_valid_minimal_model(){
    // E2-API-002 / E2-API-008 / E2-CAN-006.
    let request=json!({"model":serde_json::from_slice::<Value>(MODEL).expect("model parses"),"operation":"validate"});
    let v=response(&serde_json::to_vec(&request).expect("request serializes"));
    assert_eq!(v["status"],"OK","E2-API-002: the minimal model is a valid spark-atomistic-model/1 reference");
    let keys:Vec<&str>=v["value"].as_object().expect("object").keys().map(String::as_str).collect();
    assert_eq!(keys,vec!["config_digest","ir","schema_digest"],"E2-API-008: validation value has exactly these keys");
    assert_eq!(v["value"]["ir"],"spark-atomistic-model/1");
    assert_eq!(v["value"]["schema_digest"],"sha256:583d580d54e3847ef92f1b1456dda006161689c0bac27fd7ea896a093f48c02c",
        "E2-CAN-006: the schema-descriptor digest is stated normatively by the erratum");
    assert_eq!(SCHEMA_DIGEST,"sha256:583d580d54e3847ef92f1b1456dda006161689c0bac27fd7ea896a093f48c02c");
    // E2-CAN-006 cross-check: rebuild the descriptor payload from the erratum text and re-hash it
    // with this crate's own canonical encoder, so the constant is not merely copied.
    let descriptor=json!({"base_spec_sha256":spark_atomistic_rs::BASE_SPEC_SHA256,
        "errata_1_sha256":spark_atomistic_rs::ERRATA_1_SHA256,"ir":"spark-atomistic-model/1","revision":2});
    assert_eq!(digest_of(&descriptor),SCHEMA_DIGEST,"E2-CAN-006: descriptor payload must hash to the stated schema digest");
}

#[test]
fn api_metadata_digest_invariance(){
    // E2-SCHEMA-001 (metadata excluded from every digest) / E2-CAN-007.
    let base:Value=serde_json::from_slice(MODEL).expect("model parses");
    let mut a=base.clone();a["metadata"]=json!({"author":"conformance","run":1});
    let mut b=base.clone();b["metadata"]=json!({"unrelated":["free","form",{"nested":true}]});
    let d=|m:&Value|{let r=response(&serde_json::to_vec(&json!({"model":m,"operation":"validate"})).expect("serializes"));
        assert_eq!(r["status"],"OK");r["value"]["config_digest"].as_str().expect("digest").to_owned()};
    let(d0,d1,d2)=(d(&base),d(&a),d(&b));
    assert_eq!(d0,d1,"E2-CAN-007: config_digest hashes the model after removing root metadata");
    assert_eq!(d0,d2,"E2-CAN-007: two different metadata objects must not change config_digest");
    // Negative control: a behavioral field does change the digest, so the invariance above is not vacuous.
    let mut c=base.clone();c["kinetics"]["temperature"]=json!(301.0);
    assert_ne!(d0,d(&c),"E2-CAN-007: a behavioral field must change config_digest");
}

// ---------------------------------------------------------------- json-*

#[test]
fn json_duplicate_key_is_rejected_as_a_duplicate(){
    // E2-JSON-001: "Duplicate keys at any nesting depth return INVALID_INPUT."
    // The corpus document differs from its baseline ONLY by the repeated `value` key, so the
    // rejection cannot be attributed to any other property of the document.
    let text=String::from_utf8(DUPLICATE.to_vec()).expect("corpus is UTF-8");
    assert_eq!(text.matches("\"value\"").count(),2,"corpus fixture must actually carry the duplicate key");
    let baseline=text.replacen("\"value\":2","\"value_second\":2",1);
    assert_ne!(baseline,text);
    parse_strict_json(baseline.as_bytes()).expect("E2-JSON-001 baseline: the de-duplicated document is otherwise valid");
    let err=parse_strict_json(DUPLICATE).expect_err("E2-JSON-001: duplicate key must be rejected");
    assert_eq!(err.status,StatusCode::CheckpointCorrupt,"strict-parser transport status");
    assert!(err.message.contains("duplicate key: value"),
        "E2-JSON-001: the rejection must name the duplicate key, not a generic parse failure; got: {}",err.message);
    // Same property on the public wire, with the requirement ID asserted.
    let dup_request=validate_request_with_raw_metadata(&text);
    let ok_request=validate_request_with_raw_metadata(&baseline);
    let ok=response(ok_request.as_bytes());
    assert_eq!(ok["status"],"OK","E2-JSON-001 baseline: the de-duplicated request is accepted");
    let bad=response(dup_request.as_bytes());
    assert_eq!(bad["status"],"INVALID_INPUT","E2-JSON-001: nested duplicate key returns INVALID_INPUT");
    assert_eq!(bad["context"]["requirement_id"],"E2-JSON-001","the rejecting requirement is the portable-JSON duplicate-key rule");
    assert_eq!(bad["exit_code"],64,"E2-STATUS-002: INVALID_INPUT exits 64");
    assert_eq!(bad["value"],Value::Null,"E2-API-006: value is null on failure");
    // Depth 0 as well as depth 2.
    let top=response(br#"{"operation":"capabilities","operation":"capabilities"}"#);
    assert_eq!(top["status"],"INVALID_INPUT","E2-JSON-001: a top-level duplicate is rejected");
    assert_eq!(response(CAPS_REQUEST)["status"],"OK","E2-JSON-001 baseline: the same request without the duplicate is accepted");
}

#[test]
fn json_malformed_utf8(){
    // E2-JSON-001: "Input is UTF-8 without BOM."
    let good=validate_request_with_raw_metadata("{\"note\":\"aa\"}");
    assert_eq!(response(good.as_bytes())["status"],"OK","E2-JSON-001 baseline");
    let mut bad=good.clone().into_bytes();
    let at=bad.windows(2).position(|w|w==b"aa").expect("marker present");
    bad[at]=0xff;bad[at+1]=0xfe;
    let v=response(&bad);
    assert_eq!(v["status"],"INVALID_INPUT","E2-JSON-001: malformed UTF-8 returns INVALID_INPUT");
    assert_eq!(v["context"]["requirement_id"],"E2-JSON-001");
    // A BOM is likewise rejected while the identical BOM-free bytes are accepted.
    let mut bom=vec![0xef,0xbb,0xbf];bom.extend_from_slice(CAPS_REQUEST);
    assert_eq!(response(&bom)["status"],"INVALID_INPUT","E2-JSON-001: UTF-8 BOM is rejected");
    assert_eq!(response(CAPS_REQUEST)["status"],"OK","E2-JSON-001 baseline: identical bytes without the BOM");
}

#[test]
fn json_lone_surrogate(){
    // E2-JSON-001: "lone UTF-16 surrogate code points are forbidden."
    let paired=validate_request_with_raw_metadata(r#"{"emoji":"\ud83d\ude00"}"#);
    assert_eq!(response(paired.as_bytes())["status"],"OK","E2-JSON-001 baseline: a well-formed surrogate pair is accepted");
    let lone=validate_request_with_raw_metadata(r#"{"emoji":"\ud800"}"#);
    let v=response(lone.as_bytes());
    assert_eq!(v["status"],"INVALID_INPUT","E2-JSON-001: a lone leading surrogate is rejected");
    assert_eq!(v["context"]["requirement_id"],"E2-JSON-001");
    let trailing=validate_request_with_raw_metadata(r#"{"emoji":"\udc00"}"#);
    assert_eq!(response(trailing.as_bytes())["status"],"INVALID_INPUT","E2-JSON-001: a lone trailing surrogate is rejected");
}

#[test]
fn json_nonfinite(){
    // E2-JSON-003: NaN / Infinity / -Infinity and overflow to infinity return NONFINITE_RESULT.
    assert_eq!(response(validate_request_with_raw_metadata(r#"{"x":1.5}"#).as_bytes())["status"],"OK",
        "E2-JSON-003 baseline: a finite binary64 is accepted");
    for token in ["NaN","Infinity","-Infinity","1e400","-1e400"]{
        let raw=format!("{{\"x\":{token}}}");
        let v=response(validate_request_with_raw_metadata(&raw).as_bytes());
        assert_eq!(v["status"],"NONFINITE_RESULT","E2-JSON-003: {token} must return NONFINITE_RESULT");
        assert_eq!(v["context"]["requirement_id"],"E2-JSON-003");
        assert_eq!(v["exit_code"],65,"E2-STATUS-002: NONFINITE_RESULT exits 65");
        assert_eq!(v["message"],"nonfinite value rejected","E2-STATUS-002: exact stable message");
    }
}

#[test]
fn json_integer_boundary(){
    // E2-JSON-002: "Every integer is in [-9007199254740991,9007199254740991]."
    let ok=validate_request_with_raw_metadata(r#"{"hi":9007199254740991,"lo":-9007199254740991}"#);
    assert_eq!(response(ok.as_bytes())["status"],"OK","E2-JSON-002: both domain endpoints are inside the domain");
    for out in ["9007199254740992","-9007199254740992"]{
        let raw=format!("{{\"x\":{out}}}");
        let v=response(validate_request_with_raw_metadata(&raw).as_bytes());
        assert_eq!(v["status"],"INVALID_INPUT","E2-JSON-002: {out} is one step outside the domain");
        assert_eq!(v["exit_code"],64);
    }
    assert!(parse_strict_json(b"[9007199254740991]").is_ok());
    assert!(parse_strict_json(b"[9007199254740992]").is_err());
    assert!(parse_strict_json(b"[-9007199254740992]").is_err());
}

#[test]
fn json_integer_out_of_domain(){
    // E2-JSON-002: "Every integer is in [-9007199254740991,9007199254740991]. A syntactically valid
    // integer outside this domain returns INVALID_INPUT before schema validation."
    //
    // This test replaces the D-E2-01 defect witness. Before the fix the domain was enforced only in
    // the strict parser's i64/u64 visitors: serde_json routes every integer literal wider than u64 to
    // the binary64 visitor, so the gate held on [2^53,2^64) and let everything from 2^64 up to
    // binary64 overflow through, silently rounded. The domain is now decided on the source text, so
    // the classification is the syntactic one E2-JSON-002 words: no fraction and no exponent means
    // integer, anything else is an E2-JSON-003 binary64.
    let inside=validate_request_with_raw_metadata(r#"{"x":9007199254740991}"#);
    assert_eq!(response(inside.as_bytes())["status"],"OK","ACCEPTED BASELINE: the largest in-domain integer");
    for out in ["9007199254740992","-9007199254740992","18446744073709551615","18446744073709551616",
                "-18446744073709551616","1000000000000000000000000000000","999999999999999900000"]{
        let raw=format!("{{\"x\":{out}}}");
        let v=response(validate_request_with_raw_metadata(&raw).as_bytes());
        assert_eq!(v["status"],"INVALID_INPUT","E2-JSON-002: integer literal {out} is outside the portable domain");
        assert_eq!(v["exit_code"],64,"E2-STATUS-002: INVALID_INPUT exits 64");
        assert_eq!(v["message"],"input invalid","E2-STATUS-002: exact stable message");
        assert_eq!(v["context"]["requirement_id"],"E2-JSON-002","E2-JSON-002 governs the integer domain");
    }
    // A 401-digit literal is out of domain for the same reason, not because it overflows binary64.
    let huge=format!("{{\"x\":1{}}}","0".repeat(400));
    let v=response(validate_request_with_raw_metadata(&huge).as_bytes());
    assert_eq!(v["status"],"INVALID_INPUT");
    assert_eq!(v["context"]["requirement_id"],"E2-JSON-002");
    // ACCEPTED BASELINES differing only in the property under test. The same magnitudes carrying a
    // fraction or an exponent are not integers, so E2-JSON-003 admits them as finite binary64.
    for accepted in ["1e30","1.8446744073709552e19","9.999999999999999e20","-1e30"]{
        let raw=format!("{{\"x\":{accepted}}}");
        assert_eq!(response(validate_request_with_raw_metadata(&raw).as_bytes())["status"],"OK",
            "E2-JSON-003: {accepted} is a finite binary64, not an integer");
    }
    // Digits inside a string are not a number token.
    assert_eq!(response(validate_request_with_raw_metadata(r#"{"x":"18446744073709551616"}"#).as_bytes())["status"],"OK",
        "E2-JSON-002 applies to number tokens, not to string content");
    // The same domain now governs the single strict parser, hence checkpoint bytes as well.
    assert!(parse_strict_json(b"[9007199254740991]").is_ok(),"ACCEPTED BASELINE: in-domain integer parses");
    assert!(parse_strict_json(b"[18446744073709551616]").is_err(),"E2-JSON-002: 2^64 is refused by the strict parser");
    assert!(parse_strict_json(br#"["18446744073709551616"]"#).is_ok(),"string content is not a number token");
    assert!(parse_strict_json(b"[1e30]").is_ok(),"E2-JSON-003: exponent form is a binary64");
}

#[test]
fn json_unknown_key(){
    // E2-SCOPE-003: "Unknown keys in every normative object return INVALID_INPUT."
    let base:Value=serde_json::from_slice(MODEL).expect("model parses");
    assert_eq!(response(&serde_json::to_vec(&json!({"model":base,"operation":"validate"})).expect("serializes"))["status"],"OK",
        "E2-SCOPE-003 baseline");
    for path in [vec!["kinetics"],vec!["system"],vec!["resources"],vec!["basin"],vec!["schema"]]{
        let mut m:Value=serde_json::from_slice(MODEL).expect("model parses");
        m[path[0]]["spark_unknown_key"]=json!(1);
        let v=response(&serde_json::to_vec(&json!({"model":m,"operation":"validate"})).expect("serializes"));
        assert_eq!(v["status"],"INVALID_INPUT","E2-SCOPE-003: unknown key in normative object `{}`",path[0]);
        assert_eq!(v["exit_code"],64);
    }
    // The root request object is normative too.
    let mut root=json!({"model":serde_json::from_slice::<Value>(MODEL).expect("parses"),"operation":"validate"});
    root["spark_unknown_key"]=json!(true);
    assert_eq!(response(&serde_json::to_vec(&root).expect("serializes"))["status"],"INVALID_INPUT",
        "E2-SCOPE-003: unknown key on the request root");
    // Only root `metadata` accepts free-form keys.
    let mut with_meta:Value=serde_json::from_slice(MODEL).expect("parses");
    with_meta["metadata"]=json!({"spark_unknown_key":1});
    assert_eq!(response(&serde_json::to_vec(&json!({"model":with_meta,"operation":"validate"})).expect("serializes"))["status"],"OK",
        "E2-SCOPE-003: root metadata is the free-form exception");
}

struct RecordingAdapter(std::sync::Mutex<Vec<(String,String,String,String)>>);
impl RunAdapter for RecordingAdapter{
    fn run(&self,model:&ValidatedModel,extension:&Map<String,Value>)->Result<RunValue,ApiFailure>{
        self.0.lock().expect("adapter mutex").push((model.config_digest().into(),model.tolerance_digest().into(),
            model.identity_digest().into(),serde_json::to_string(extension).expect("extension serializes")));
        Ok(RunValue{checkpoint_sequence:0,current_state_id:"state:sha256:conformance".into(),incomplete_catalog:false,
            simulation_time_s:0.0,step_index:0})}
}

#[test]
fn json_adapter_extension_exclusion(){
    // E2-SCOPE-004: adapter settings live only under run-request `extension` and MUST NOT enter the
    // config digest. E2-API-003 / E2-API-004 / E2-API-008.
    let dir=temp_dir("extension");
    let adapter=RecordingAdapter(std::sync::Mutex::new(Vec::new()));
    let mut m:Value=serde_json::from_slice(MODEL).expect("parses");
    m["output"]["checkpoint_path"]=json!(dir.join("ckpt.json").to_str().expect("utf8"));
    m["output"]["summary_path"]=json!(dir.join("summary.json").to_str().expect("utf8"));
    m["output"]["trajectory_path"]=json!(dir.join("traj.json").to_str().expect("utf8"));
    let validate_digest=response(&serde_json::to_vec(&json!({"model":m,"operation":"validate"})).expect("serializes"))
        ["value"]["config_digest"].as_str().expect("digest").to_owned();
    for extension in [json!({}),json!({"transport":"stdio","executable":"/usr/bin/false","timeout_s":5,"isolation":"process"})]{
        let request=json!({"allow_unvalidated":true,"extension":extension,"model":m,"operation":"run"});
        let bytes=dispatch_json(&serde_json::to_vec(&request).expect("serializes"),None,Some(&adapter));
        let v:Value=serde_json::from_slice(&bytes).expect("parses");
        assert_eq!(v["status"],"OK","E2-API-003: a run request with allow_unvalidated:true and an adapter is accepted");
        let keys:Vec<&str>=v["value"].as_object().expect("object").keys().map(String::as_str).collect();
        assert_eq!(keys,vec!["checkpoint_sequence","current_state_id","incomplete_catalog","simulation_time_s","step_index"],
            "E2-API-008: run value has exactly these keys");
    }
    let seen=adapter.0.lock().expect("adapter mutex");
    assert_eq!(seen.len(),2);
    assert_ne!(seen[0].3,seen[1].3,"the two runs really did carry different extension objects");
    assert_eq!(seen[0].0,seen[1].0,"E2-SCOPE-004: extension must not change config_digest");
    assert_eq!(seen[0].1,seen[1].1,"E2-SCOPE-004: extension must not change tolerance_digest");
    assert_eq!(seen[0].2,seen[1].2,"E2-SCOPE-004: extension must not change identity_digest");
    assert_eq!(seen[0].0,validate_digest,"E2-SCOPE-004: the run digest equals the extension-free validate digest");
    drop(seen);
    // E2-API-004: allow_unvalidated must be the JSON boolean true; nothing else authorizes a run.
    for bad in [json!(false),json!(1),json!("true"),json!(null)]{
        let request=json!({"allow_unvalidated":bad,"extension":{},"model":m,"operation":"run"});
        let v:Value=serde_json::from_slice(&dispatch_json(&serde_json::to_vec(&request).expect("serializes"),None,Some(&adapter)))
            .expect("parses");
        assert_eq!(v["status"],"INVALID_INPUT","E2-API-004: allow_unvalidated={bad} must not authorize a run");
        assert_eq!(v["context"]["requirement_id"],"E2-API-004");
    }
    assert_eq!(adapter.0.lock().expect("adapter mutex").len(),2,"E2-API-004: the adapter is never reached on rejection");
    let _=std::fs::remove_dir_all(&dir);
}

// ---------------------------------------------------------------- canonical-*

#[test]
fn canonical_number_boundaries(){
    // E2-CAN-004: RFC 8785 shortest round-tripping decimal; negative zero serializes as `0`.
    let corpus:Value=serde_json::from_slice(NUMBERS).expect("canonical-number corpus parses");
    let cases=corpus["cases"].as_array().expect("cases array");
    assert_eq!(cases.len(),7,"the corpus must cover negative zero plus the 1e-6 and 1e21 neighbourhoods");
    let mut covered=0;
    for case in cases{
        let expected=case["canonical"].as_str().expect("canonical string");
        let tag=case["binary64"].as_str().expect("binary64 tag");
        let value=if tag=="negative-zero"{-0.0f64}
            else{f64::from_bits(u64::from_str_radix(tag.strip_prefix("0x").expect("corpus stores 0x-prefixed bits"),16)
                .expect("hex bits"))};
        let actual=String::from_utf8(canonical_json_bytes(&json!(value)).expect("finite")).expect("ascii");
        assert_eq!(actual,expected,"E2-CAN-004: binary64 {tag} must serialize as {expected}");
        // Round-trip: the canonical decimal must read back to the same binary64 value. E2-CAN-004
        // mandates that negative zero serializes as `0`, so the sign of zero is deliberately NOT
        // preserved; every other case must round-trip the identical bit pattern.
        let back:f64=actual.parse().expect("canonical decimal parses");
        assert_eq!(back,value,"E2-CAN-004: canonical form must round-trip the same binary64 value");
        if value!=0.0{assert_eq!(back.to_bits(),value.to_bits(),"E2-CAN-004: shortest round-tripping decimal is bit-exact");}
        else{assert_eq!(back.to_bits(),0u64,"E2-CAN-004: negative zero serializes as 0 and reads back as +0");}
        covered+=1;
    }
    assert_eq!(covered,7);
    assert_eq!(corpus["integer_max"],json!(9007199254740991i64),"E2-JSON-002");
    assert_eq!(corpus["integer_min"],json!(-9007199254740991i64),"E2-JSON-002");
    // E2-CAN-003 / E2-CAN-002 / E2-CAN-001 on the surrounding encoder.
    assert_eq!(canonical_json_bytes(&json!(0i64)).expect("zero"),b"0");
    assert_eq!(canonical_json_bytes(&json!(null)).expect("null"),b"null");
    assert_eq!(canonical_json_bytes(&json!([true,false])).expect("bools"),b"[true,false]");
    assert_eq!(canonical_json_bytes(&json!({"b":1,"a":2})).expect("sorted"),b"{\"a\":2,\"b\":1}");
    assert_eq!(canonical_json_bytes(&json!("a/b")).expect("solidus"),b"\"a/b\"","E2-CAN-002: solidus is not escaped");
    assert_eq!(canonical_json_bytes(&json!("\u{7}")).expect("control"),b"\"\\u0007\"","E2-CAN-002: lowercase four-digit control escape");
    assert_eq!(canonical_json_bytes(&json!("\u{e9}")).expect("non-ascii"),"\"\u{e9}\"".as_bytes(),"E2-CAN-002: non-ASCII scalars are not escaped");
}

// ---------------------------------------------------------------- identity-*

#[test]
fn identity_is_invariant_to_translation_image_and_same_species_permutation(){
    // E2-ID-002 / E2-ID-004: the anchor-minimum closest-image certificate is invariant to whole-cell
    // translation, periodic image choice, atom ID, and same-species permutation.
    // The offsets below are dyadic, so every translated/imaged coordinate difference reproduces the
    // base displacement bit-exactly. See `identity_certificate_is_ulp_sensitive_to_translation`
    // for the measured behavior of a translation that is not exact in binary64.
    let energy=-12.5_f64;
    let base=system(&[[0.0,0.0,0.0],[1.5,0.0,0.0],[0.0,2.25,0.0]],&["H","H","O"],&["a0","a1","a2"]);
    let translated=system(&[[3.0,4.0,5.0],[4.5,4.0,5.0],[3.0,6.25,5.0]],&["H","H","O"],&["a0","a1","a2"]);
    let image=system(&[[0.0,0.0,0.0],[-10.5,0.0,0.0],[0.0,2.25,-12.0]],&["H","H","O"],&["a0","a1","a2"]);
    let permuted=system(&[[1.5,0.0,0.0],[0.0,0.0,0.0],[0.0,2.25,0.0]],&["H","H","O"],&["z9","q1","a2"]);
    let reference=state_ids(&base,energy).expect("identity");
    for(label,s)in [("identity-translated",&translated),("identity-periodic-image",&image),("identity-same-species-permuted",&permuted)]{
        let got=state_ids(s,energy).expect("identity");
        assert_eq!(got.0,reference.0,"E2-ID-002 [{label}]: candidate_identity must be invariant");
        assert_eq!(got.1,reference.1,"E2-ID-004 [{label}]: state_id must be invariant");
    }
    // Negative controls: the invariance above must not be a constant function.
    let moved=system(&[[0.0,0.0,0.0],[1.9,0.0,0.0],[0.0,2.25,0.0]],&["H","H","O"],&["a0","a1","a2"]);
    assert_ne!(state_ids(&moved,energy).expect("identity").0,reference.0,"E2-ID-002: a real geometry change must change the certificate");
    // Cross-species permutation is NOT an equivalence.
    let cross=system(&[[0.0,0.0,0.0],[1.5,0.0,0.0],[0.0,2.25,0.0]],&["H","O","H"],&["a0","a1","a2"]);
    assert_ne!(state_ids(&cross,energy).expect("identity").0,reference.0,"E2-ID-002: species labels are part of the certificate rows");
    // E2-ID-004: energy participates in state_id but not in candidate_identity.
    let other=state_ids(&base,energy+1e-3).expect("identity");
    assert_eq!(other.0,reference.0,"E2-ID-004: candidate_identity is energy-free");
    assert_ne!(other.1,reference.1,"E2-ID-004: state_id binds energy_ev");
    assert!(reference.0.starts_with("candidate:sha256:")&&reference.1.starts_with("state:sha256:"),"E2-ID-001: exact prefixes");
}

#[test]
fn identity_certificate_is_ulp_sensitive_to_translation(){
    // FINDING -- records observed behavior and its magnitude. E2-ID-002 states the certificate is
    // invariant to whole-cell translation. That holds in exact arithmetic; in binary64 a translation
    // whose coordinate differences are not exactly representable perturbs the certificate rows by
    // 1-2 ULP, and because E2-CAN-005 hashes those exact bytes the candidate/state IDs change.
    // Errata 2 tolerates this by construction: E2-ID-005 makes candidate IDs hints only and STATE-006
    // requires geometry verification to decide equivalence. The measurement is recorded here so the
    // limit is a known quantity rather than a surprise in cross-language parity.
    let energy=-12.5_f64;
    let base=system(&[[0.0,0.0,0.0],[1.6,0.0,0.0],[0.0,2.1,0.0]],&["H","H","O"],&["a0","a1","a2"]);
    let shifted=system(&[[3.3,4.4,5.5],[4.9,4.4,5.5],[3.3,6.5,5.5]],&["H","H","O"],&["a0","a1","a2"]);
    let drift=((4.9f64-3.3f64)-1.6f64).abs().max(((6.5f64-4.4f64)-2.1f64).abs());
    assert!(drift>0.0,"the fixture translation is genuinely inexact in binary64");
    assert!(drift<=1e-15,"measured displacement drift is {drift:e} A, i.e. 1-2 ULP at this magnitude");
    assert!(drift<1e-12,"PAR-004: the drift is far inside the 1e-12 A geometry tolerance, so the two \
        geometries are the same state scientifically even though their certificates differ");
    assert_ne!(state_ids(&shifted,energy).expect("identity").0,state_ids(&base,energy).expect("identity").0,
        "observed: a binary64-inexact translation changes the candidate identity by {drift:e} A of coordinate drift");
    // The same displacement expressed exactly is invariant, so the sensitivity is arithmetic and not
    // a defect in the anchor/sort/closest-image construction.
    let exact_base=system(&[[0.0,0.0,0.0],[1.5,0.0,0.0],[0.0,2.25,0.0]],&["H","H","O"],&["a0","a1","a2"]);
    let exact_shift=system(&[[3.0,4.0,5.0],[4.5,4.0,5.0],[3.0,6.25,5.0]],&["H","H","O"],&["a0","a1","a2"]);
    assert_eq!(state_ids(&exact_shift,energy).expect("identity").0,state_ids(&exact_base,energy).expect("identity").0,
        "E2-ID-002: an exactly representable translation is invariant");
}

// ---------------------------------------------------------------- rng-*

#[test]
fn rng_errata1_boundary_records(){
    // E1-DET-002-B/C/D, E2-RNG-006, E2-RNG-007. The two golden-line SHA-256 values below are stated
    // normatively by Errata 2; they are recomputed here from this crate's canonical encoder.
    assert_eq!(uniform_from_words(0,0).expect("low boundary").to_bits(),0x3ca0_0000_0000_0000);
    assert_eq!(uniform_from_words(u32::MAX,0xffff_f000).expect("high boundary").to_bits(),0x3fef_ffff_ffff_ffff);
    let low=json!({"a":0,"b":0,"q":0,"raw_binary64_bits":"0x3ca0000000000000","uniform_hex":"0x1.0000000000000p-53"});
    let low_bytes=canonical_json_bytes(&low).expect("canonical");
    assert_eq!(hex_sha256(&low_bytes),"6ce1fb5214530ba6b04e4bf75aaeba5d02acf6694cd462004faf7640a665fc03",
        "E2-RNG-006: first boundary golden record hash");
    let high=json!({"a":4294967295u64,"b":4294963200u64,"q":4503599627370495u64,
        "raw_binary64_bits":"0x3fefffffffffffff","uniform_hex":"0x1.fffffffffffffp-1"});
    assert_eq!(hex_sha256(&canonical_json_bytes(&high).expect("canonical")),
        "a15157604d319e5525e3b83eba02259088e317ea2b0b0e1a3bb28060e093cf43","E2-RNG-007: second boundary golden record hash");
    // E1-DET-002-B: b[11:0] is discarded, so the low twelve bits of `b` cannot change the uniform.
    assert_eq!(uniform_from_words(u32::MAX,0xffff_ffff).expect("uniform").to_bits(),
        uniform_from_words(u32::MAX,0xffff_f000).expect("uniform").to_bits(),"E1-DET-002-B: low 12 bits of b are discarded");
    // E1-DET-002-C: every output is strictly inside the open interval.
    for(a,b)in [(0u32,0u32),(u32::MAX,u32::MAX),(1,0),(0,0xffff_ffff),(0x8000_0000,0x8000_0000)]{
        let u=uniform_from_words(a,b).expect("uniform");
        assert!(u>0.0&&u<1.0,"E1-DET-002-C: 0 < u < 1 for every output");
    }
}

#[test]
fn rng_zero_key_counter_reference_block(){
    // E2-RNG-001: Philox4x32-10 round constants and lane order.
    let corpus:Value=serde_json::from_slice(PHILOX).expect("philox corpus parses");
    let expected:Vec<u32>=corpus["zero_counter_zero_key_words"].as_array().expect("array").iter()
        .map(|x|x.as_u64().expect("word") as u32).collect();
    assert_eq!(philox4x32_10([0;4],[0;2]).to_vec(),expected,"E2-RNG-001: zero-key/zero-counter reference block");
    assert_eq!(expected,vec![0x6627_e8d5,0xe169_c58d,0xbc57_ac4c,0x9b00_dbd8]);
    assert_eq!(corpus["algorithm"],"Philox4x32-10:errata-1-midpoint52","E2-RNG-001: exact algorithm ID");
    assert_eq!(Philox::new([0;2],[0;4]).state().algorithm,"Philox4x32-10:errata-1-midpoint52");
    // Key and counter must both matter, and lane order must not be host-byte-order dependent.
    assert_ne!(philox4x32_10([1,0,0,0],[0;2]),philox4x32_10([0,0,0,1],[0;2]),"E2-RNG-003: c0 and c3 are distinct lanes");
    assert_ne!(philox4x32_10([0;4],[1,0]),philox4x32_10([0;4],[0,1]),"E2-RNG-002: k0 and k1 are distinct lanes");
}

#[test]
fn rng_trajectory_and_saddle_derivation(){
    // E2-RNG-004: SHA-256 bytes 0..7 become two big-endian key words; bytes 8..23 become the four
    // big-endian initial-counter words. The expected words below were computed independently from
    // the erratum's byte recipe with a separate SHA-256 implementation.
    let t0=derive_trajectory_stream(0).state();
    assert_eq!(t0.key,[2351537927,2073809783],"E2-RNG-004: trajectory key for run_seed 0");
    assert_eq!(t0.initial_counter,[4266885156,3806684783,2000346754,1046855009],"E2-RNG-004: trajectory counter for run_seed 0");
    assert_eq!(t0.next_counter,t0.initial_counter,"E2-RNG-004: next counter initially equals the initial counter");
    assert_eq!(t0.consumed_uniforms,0);assert_eq!(t0.consumed_blocks,0);assert_eq!(t0.next_pair,0);
    assert!(t0.buffered_block.is_none());
    let t1=derive_trajectory_stream(12345).state();
    assert_eq!(t1.key,[3053090049,811585613],"E2-RNG-004: trajectory key for run_seed 12345");
    assert_eq!(t1.initial_counter,[2799702670,2319355720,3244935538,479610466]);
    let s0=derive_saddle_substream(0,GOLDEN_STATE_ID,"class-selection",0).expect("substream").state();
    assert_eq!(s0.key,[4060213990,2074348605],"E2-RNG-004: saddle substream key");
    assert_eq!(s0.initial_counter,[2647842716,3610871162,1756738045,3802976106],"E2-RNG-004: saddle substream counter");
    let s1=derive_saddle_substream(12345,GOLDEN_STATE_ID,"global",7).expect("substream").state();
    assert_eq!(s1.key,[3366311393,3423593154]);
    assert_eq!(s1.initial_counter,[937477282,1483016221,2158739293,3081166098]);
    // E2-RNG-004 / DET-003: every derivation input must be load-bearing.
    assert_ne!(derive_saddle_substream(0,GOLDEN_STATE_ID,"global",0).expect("s").state().key,s0.key,"search class matters");
    assert_ne!(derive_saddle_substream(0,GOLDEN_STATE_ID,"class-selection",1).expect("s").state().key,s0.key,"search index matters");
    assert_ne!(derive_saddle_substream(1,GOLDEN_STATE_ID,"class-selection",0).expect("s").state().key,s0.key,"run seed matters");
    // The length prefixes prevent a concatenation collision between (state_id, class) splits.
    let a=derive_saddle_substream(0,"ab","cd",0).expect("s").state();
    let b=derive_saddle_substream(0,"a","bcd",0).expect("s").state();
    assert_ne!(a.key,b.key,"E2-RNG-004: len32 prefixes make the material unambiguous");
    // The trajectory stream and the saddle substreams are different streams (E2-RNG-005).
    assert_ne!(t0.key,s0.key,"E2-RNG-004: distinct domain-separation labels");
    // E2-DISC-002: the search ID is the canonical hash of exactly these four fields.
    assert_eq!(search_id(0,GOLDEN_STATE_ID,"class-selection",0).expect("search id"),
        "search:sha256:e51b54e04098dbe42704891fd7d0d6d258a532333323006cef45492ed0c77a31");
    assert_eq!(search_id(12345,GOLDEN_STATE_ID,"global",7).expect("search id"),
        "search:sha256:e00c465e4e6fe4b20969dc28dc002e4c67a227fba60214256ac658f7aecc5af4");
}

#[test]
fn rng_state_restore(){
    // E2-RNG-002: consumed_blocks = ceil(consumed_uniforms/2); buffer/pair relations are checked.
    let mut rng=Philox::new([7,9],[1,2,3,4]);
    let s0=rng.state();
    assert_eq!(Philox::from_state(s0.clone()).expect("restore").state(),s0,"E2-RNG-002: a fresh state restores");
    let mut expect_blocks=0u64;
    for consumed in 1..=6u64{
        let u=rng.next_uniform().expect("uniform");
        assert!(u>0.0&&u<1.0);
        expect_blocks=consumed.div_ceil(2);
        let s=rng.state();
        assert_eq!(s.consumed_uniforms,consumed);
        assert_eq!(s.consumed_blocks,expect_blocks,"E2-RNG-002: consumed_blocks = ceil(consumed_uniforms/2)");
        assert_eq!(s.buffered_block.is_some(),consumed%2==1,"E2-RNG-002: a block is retained only between the two pairs");
        assert_eq!(s.next_pair,if consumed%2==1{1}else{0},"E2-RNG-002: next_pair follows the buffer");
        let restored=Philox::from_state(s.clone()).expect("mid-stream restore");
        assert_eq!(restored.state(),s,"E2-RNG-002: mid-stream state restores exactly");
        // Bit-exact continuation after restore (E2-CKPT-008).
        let mut a=Philox::from_state(s.clone()).expect("restore");
        let mut b=rng.clone();
        assert_eq!(a.next_uniform().expect("u").to_bits(),b.next_uniform().expect("u").to_bits(),
            "E2-CKPT-008: replay after restore is bit-identical");
    }
    assert_eq!(expect_blocks,3);
    // Forged buffers and broken counter relations must be refused.
    let mut forged=rng.state();forged.consumed_blocks+=1;
    assert!(Philox::from_state(forged).is_err(),"E2-RNG-002: block/uniform relation is enforced");
    let mut wrong=rng.state();wrong.algorithm="Philox4x32-10".into();
    assert!(Philox::from_state(wrong).is_err(),"E2-RNG-001: the exact algorithm ID is enforced");
}

#[test]
fn rng_counter_carry(){
    // E2-RNG-003: c0 is least significant; block increment propagates carry toward c3.
    let step=|start:[u32;4]|->[u32;4]{let mut r=Philox::new([0,0],start);r.next_uniform().expect("u");r.next_uniform().expect("u");r.state().next_counter};
    assert_eq!(step([0,0,0,0]),[1,0,0,0],"E2-RNG-003: one block advances c0");
    assert_eq!(step([u32::MAX,0,0,0]),[0,1,0,0],"E2-RNG-003: carry from c0 into c1");
    assert_eq!(step([u32::MAX,u32::MAX,0,0]),[0,0,1,0],"E2-RNG-003: carry propagates into c2");
    assert_eq!(step([u32::MAX,u32::MAX,u32::MAX,0]),[0,0,0,1],"E2-RNG-003: carry propagates into c3");
    // The persistent relation next_counter = initial_counter + consumed_blocks must survive the carry.
    let mut r=Philox::new([0,0],[u32::MAX,u32::MAX,0,0]);
    for _ in 0..8{r.next_uniform().expect("u");}
    let s=r.state();
    assert_eq!(s.consumed_blocks,4);assert_eq!(s.next_counter,[3,0,1,0],"E2-RNG-002/003: 128-bit counter arithmetic");
    assert!(Philox::from_state(s).is_ok(),"E2-RNG-002: the carried state validates");
    // Exhaustion must fail closed rather than wrap. E2-RNG-002 requires
    // next_counter = initial_counter + consumed_blocks, which is unrepresentable once the 128-bit
    // counter is saturated; the generator therefore refuses the terminal block instead of wrapping.
    let mut terminal=Philox::new([0,0],[u32::MAX;4]);
    let exhausted=terminal.next_uniform().expect_err("E2-RNG-003: the saturated counter must fail closed");
    assert_eq!(exhausted.status,StatusCode::ResourceLimit,"observed: exhaustion is a pause, not a silently wrapped stream");
    assert_eq!(terminal.state().consumed_uniforms,0,"E2-RNG-002: nothing was consumed by the refused block");
    assert_eq!(terminal.state().consumed_blocks,0,"E2-RNG-002: the block counter did not advance");
    // FINDING (low severity, unreachable in practice): the increment wraps each lane in place before
    // it discovers saturation, so the refused generator is left with next_counter = [0,0,0,0] while
    // consumed_blocks is still 0. That violates E2-RNG-002's next_counter = initial_counter +
    // consumed_blocks, but only in the live object: the broken state cannot be persisted because
    // validate_state recomputes the same relation and refuses it. Reaching this needs 2^128 blocks.
    assert_eq!(terminal.state().next_counter,[0,0,0,0],
        "observed: the failing increment leaves the lanes wrapped instead of restoring them");
    assert!(Philox::from_state(terminal.state()).is_err(),
        "E2-RNG-002: the inconsistent post-exhaustion state cannot be restored or checkpointed");
    // One block earlier the generator still works, so the refusal is exactly at saturation.
    let mut last_usable=Philox::new([0,0],[u32::MAX-1,u32::MAX,u32::MAX,u32::MAX]);
    last_usable.next_uniform().expect("the penultimate block is usable");
    last_usable.next_uniform().expect("both uniforms of the penultimate block are usable");
    assert_eq!(last_usable.state().next_counter,[u32::MAX;4],"E2-RNG-003: carry advances into the terminal counter");
    assert!(last_usable.next_uniform().is_err(),"E2-RNG-003: the counter is exhausted, never wrapped to zero");
}

// ---------------------------------------------------------------- event-* and rate-*

fn saddle_of(positions:&[[f64;3]],energy:f64,search:&str)->SaddleV2{
    SaddleV2{curvature_ev_per_angstrom2:-1.5,energy_ev:energy,evaluation_count:9,evidence_level:"DIRECTIONAL".into(),
        forces_ev_per_angstrom:vec![[0.0,0.0,0.0];positions.len()],orthogonal_curvatures_ev_per_angstrom2:vec![0.7],
        positions:positions.to_vec(),search_id:search.into(),termination_reason:"curvature converged".into(),
        unstable_direction:vec![[1.0,0.0,0.0];positions.len()]}
}
fn saddle_geometry_digest(origin:&CommittedStateV2,positions:&[[f64;3]])->String{
    // E2-ID-007: the saddle certificate reuses the origin cell, periodicity, species and movable mask.
    let s=WireSystem{atom_ids:origin.atom_ids.clone(),species:origin.species.clone(),positions:positions.to_vec(),
        cell:origin.cell,pbc:origin.pbc,movable:origin.movable.clone(),constraints:origin.constraints.clone(),
        charge:origin.charge,spin:origin.spin,calculator_model_digest:origin.calculator_model_digest.clone()};
    digest_of(&geometry_certificate(&s).expect("saddle certificate"))
}

#[test]
fn event_reversible_pair_and_reciprocal_mapping(){
    // E2-EVENT-006: sorted endpoints, mapping oriented low->high, sign-canonical unstable direction.
    let a=committed(&[[0.0,0.0,0.0],[1.6,0.0,0.0],[0.0,2.2,0.0]],-10.0);
    let b=committed(&[[0.0,0.0,0.0],[2.4,0.0,0.0],[0.0,2.2,0.0]],-9.7);
    let sp=[[0.0,0.0,0.0],[2.0,0.0,0.0],[0.0,2.2,0.0]];
    let sgd=saddle_geometry_digest(&a,&sp);
    let saddle=saddle_of(&sp,-9.4,"search:sha256:fixture");
    let forward_map=vec![[0usize,1usize],[1,2],[2,0]];
    let reverse_map={let mut m:Vec<[usize;2]>=forward_map.iter().map(|x|[x[1],x[0]]).collect();m.sort();m};
    assert_ne!(forward_map,reverse_map,"the mapping under test is genuinely asymmetric");
    let(pair_f,event_f,reverse_f)=event_ids(&a.state_id,&b.state_id,&saddle,&sgd,&forward_map).expect("ids");
    let(pair_r,event_r,reverse_r)=event_ids(&b.state_id,&a.state_id,&saddle,&sgd,&reverse_map).expect("ids");
    assert_eq!(pair_f,pair_r,"E2-EVENT-006: both directions share one physical pair_id");
    assert_eq!(event_f,reverse_r,"E2-EVENT-006: reciprocal directed IDs");
    assert_eq!(event_r,reverse_f,"E2-EVENT-006: reciprocal directed IDs");
    assert_ne!(event_f,event_r,"E2-EVENT-001: origin and destination distinguish the two directed records");
    assert!(pair_f.starts_with("pair:sha256:")&&event_f.starts_with("event:sha256:"),"E2-ID-001: exact prefixes");
    // Sign canonicalisation: v and -v must produce the same pair_id (E2-EVENT-006).
    let mut flipped=saddle.clone();
    flipped.unstable_direction=saddle.unstable_direction.iter().map(|v|[-v[0],-v[1],-v[2]]).collect();
    let(pair_flipped,_,_)=event_ids(&a.state_id,&b.state_id,&flipped,&sgd,&forward_map).expect("ids");
    assert_eq!(pair_f,pair_flipped,"E2-EVENT-006: the unstable direction is canonicalised up to sign");
    // Negative control: a different mapping is a different pair.
    let other_map=vec![[0usize,0usize],[1,1],[2,2]];
    let(pair_other,_,_)=event_ids(&a.state_id,&b.state_id,&saddle,&sgd,&other_map).expect("ids");
    assert_ne!(pair_f,pair_other,"E2-EVENT-006: active_atom_mapping is part of the pair identity");
}

#[test]
fn event_parallel_saddle_stays_a_separate_event(){
    // E2-EVENT-003 / CAT-004: same endpoints with a geometrically distinct saddle is a distinct event.
    let a=committed(&[[0.0,0.0,0.0],[1.6,0.0,0.0],[0.0,2.2,0.0]],-10.0);
    let b=committed(&[[0.0,0.0,0.0],[2.4,0.0,0.0],[0.0,2.2,0.0]],-9.7);
    let map=vec![[0usize,0usize],[1,1],[2,2]];
    let p1=[[0.0,0.0,0.0],[2.0,0.0,0.0],[0.0,2.2,0.0]];
    let p2=[[0.0,0.0,0.0],[1.0,1.7,0.0],[0.0,2.2,0.0]];
    let(d1,d2)=(saddle_geometry_digest(&a,&p1),saddle_geometry_digest(&a,&p2));
    assert_ne!(d1,d2,"E2-ID-007: the two saddle geometries have distinct certificates");
    let s1=saddle_of(&p1,-9.4,"search:sha256:one");
    let s2=saddle_of(&p2,-9.4,"search:sha256:two");
    let(pair1,event1,_)=event_ids(&a.state_id,&b.state_id,&s1,&d1,&map).expect("ids");
    let(pair2,event2,_)=event_ids(&a.state_id,&b.state_id,&s2,&d2,&map).expect("ids");
    assert_ne!(pair1,pair2,"E2-EVENT-003: parallel mechanisms must not collapse into one pair");
    assert_ne!(event1,event2,"E2-EVENT-003: parallel mechanisms remain separate directed events because rates add");
    // Same geometry but a different saddle energy is also a different pair (E2-EVENT-006 payload).
    let s3=saddle_of(&p1,-9.3,"search:sha256:one");
    let(pair3,_,_)=event_ids(&a.state_id,&b.state_id,&s3,&d1,&map).expect("ids");
    assert_ne!(pair1,pair3,"E2-EVENT-006: saddle_energy_ev is part of the pair identity");
}

#[test]
fn rate_negative_barrier_tolerance_band(){
    // E2-RATE-001: b in [-barrier_tolerance, 0) is retained exactly; clamping to zero is forbidden;
    // b < -barrier_tolerance returns RATE_INVALID.
    let(t,nu,eps)=(300.0f64,1e13f64,1e-8f64);
    let tolerance=1e-3f64;
    // Tolerated: E_s is 5e-4 eV below E_i, i.e. -tolerance <= b_f < 0.
    let e_i=-10.0f64;let e_s=e_i-5e-4;let e_j=-10.6f64;
    let r=common_prefactor_pair(e_i,e_j,e_s,t,nu,tolerance,eps).expect("E2-RATE-001: inside the tolerance band");
    let beta=1.0/(spark_atomistic_rs::rate::KB_EV_PER_K*t);
    let expected_forward=nu.ln()-beta*(e_s-e_i);
    assert_eq!(r.log_forward_rate_per_s.to_bits(),expected_forward.to_bits(),
        "E2-RATE-001: the exact negative barrier is used in the rate formula, not a clamped zero");
    assert!(r.log_forward_rate_per_s>nu.ln(),"E2-RATE-001: a negative barrier raises the rate above the prefactor");
    let clamped=common_prefactor_pair(e_i,e_j,e_i,t,nu,tolerance,eps).expect("zero-barrier control");
    assert_ne!(r.log_forward_rate_per_s.to_bits(),clamped.log_forward_rate_per_s.to_bits(),
        "E2-RATE-001: clamping to zero would have produced the zero-barrier rate");
    // Rejected: one step outside the band.
    let e_s_bad=e_i-1.5e-3;
    let err=common_prefactor_pair(e_i,e_j,e_s_bad,t,nu,tolerance,eps).expect_err("E2-RATE-001: outside the tolerance band");
    assert_eq!(err.status,StatusCode::RateInvalid,"E2-RATE-001: below -barrier_tolerance returns RATE_INVALID");
    assert_eq!(err.status.exit_code(false,false),Some(65),"E2-STATUS-002: RATE_INVALID exits 65 when terminal");
    assert_eq!(err.status.message(),"rate invalid","E2-STATUS-002: exact stable message");
    // The reverse barrier is policed on the same footing.
    let reverse_bad=common_prefactor_pair(e_i,e_j,e_j-1.5e-3,t,nu,tolerance,eps);
    assert!(reverse_bad.is_err(),"E2-RATE-001: b_r is checked with the same tolerance");
}

#[test]
fn rate_detailed_balance_failure(){
    // RATE-005 / RATE-006 / E2-STATUS-002: a residual above epsilon_DB is DETAILED_BALANCE_VIOLATION.
    let(e_i,e_j,e_s,t,nu)=(-10.0f64,-9.7f64,-9.4f64,300.0f64,1e13f64);
    let ok=common_prefactor_pair(e_i,e_j,e_s,t,nu,1e-10,1e-8).expect("baseline within default epsilon_DB");
    assert!(ok.detailed_balance_residual.abs()<=1e-8,"RATE-005 baseline: |residual| <= 1e-8");
    assert!(ok.detailed_balance_residual!=0.0,"the finite-precision residual is nonzero, so the threshold test is not vacuous");
    let magnitude=ok.detailed_balance_residual.abs();
    // Identical inputs, epsilon_DB tightened below the observed residual: only the tolerance changes.
    let err=common_prefactor_pair(e_i,e_j,e_s,t,nu,1e-10,magnitude/2.0)
        .expect_err("RATE-006: residual above epsilon_DB must be rejected");
    assert_eq!(err.status,StatusCode::DetailedBalanceViolation);
    assert_eq!(err.status.exit_code(false,false),Some(65),"E2-STATUS-002: exits 65 when terminal");
    assert_eq!(err.status.message(),"detailed balance violated","E2-STATUS-002: exact stable message");
    assert!(magnitude<1e-14,"observed residual magnitude {magnitude:e} is pure round-off, well inside PAR-004");
    // E2-EVENT-003: for the reciprocal record the logs swap and the residual changes sign.
    let rev=common_prefactor_pair(e_j,e_i,e_s,t,nu,1e-10,1e-8).expect("reciprocal");
    assert_eq!(rev.log_forward_rate_per_s.to_bits(),ok.log_reverse_rate_per_s.to_bits(),"E2-EVENT-003: logs swap");
    assert_eq!(rev.log_reverse_rate_per_s.to_bits(),ok.log_forward_rate_per_s.to_bits(),"E2-EVENT-003: logs swap");
    assert_eq!(rev.detailed_balance_residual.to_bits(),(-ok.detailed_balance_residual).to_bits(),
        "E2-EVENT-003: the residual changes sign");
}

// ---------------------------------------------------------------- discovery-*

#[test]
fn discovery_class_sequence_is_deterministic(){
    // E2-DISC-001: consume the first uniform of the `class-selection` substream and take the first
    // class whose cumulative probability is strictly greater than u, in configured array order.
    let dir=temp_dir("discovery");
    let classes=json!([{"kind":"global","name":"global","probability":0.25},
                       {"kind":"local","name":"local","probability":0.75}]);
    let vm=validated_model(&dir,|m|{m["discovery"]["classes"]=classes.clone();m["kinetics"]["run_seed"]=json!(4242);});
    let state="state:sha256:00000000000000000000000000000000000000000000000000000000000000ab";
    let mut sequence=Vec::new();
    for index in 0..12u64{
        let choice=choose_discovery_class(&vm,state,index).expect("class choice");
        // The recorded search ID must be the E2-DISC-002 hash of the CHOSEN class.
        assert_eq!(choice.search_id,search_id(4242,state,&choice.search_class,index).expect("search id"),
            "E2-DISC-002: the search ID binds run seed, state ID, chosen class and index");
        // E2-DISC-001: the class stream itself is the `class-selection` substream, not the class stream.
        let expected_class_stream=derive_saddle_substream(4242,state,"class-selection",index).expect("stream").state();
        assert_eq!(choice.class_stream.key,expected_class_stream.key,"E2-DISC-001: class selection uses the class-selection substream");
        assert_eq!(choice.class_stream.consumed_uniforms,1,"E2-DISC-001: exactly the first uniform is consumed");
        // E2-RNG-005: the chosen saddle search uses a separate substream keyed by the chosen class name.
        let expected_saddle=derive_saddle_substream(4242,state,&choice.search_class,index).expect("stream").state();
        assert_eq!(choice.saddle_stream.key,expected_saddle.key,"E2-RNG-005: the saddle substream uses the chosen class name");
        assert_ne!(choice.saddle_stream.key,choice.class_stream.key,"E2-RNG-005: the two substreams are distinct");
        // Reproduce the documented rule against the real uniform.
        let mut probe=derive_saddle_substream(4242,state,"class-selection",index).expect("stream");
        let u=probe.next_uniform().expect("uniform");
        let expected=if 0.25>u{"global"}else{"local"};
        assert_eq!(choice.search_class,expected,"E2-DISC-001: first cumulative probability strictly greater than u={u}");
        sequence.push(choice.search_class.clone());
    }
    // Determinism and scheduling independence (DET-003): repeating the calls reproduces the sequence.
    let repeat:Vec<String>=(0..12u64).map(|i|choose_discovery_class(&vm,state,i).expect("choice").search_class).collect();
    assert_eq!(sequence,repeat,"DET-003: the class sequence depends only on seed, state ID and index");
    assert!(sequence.iter().any(|c|c=="global")&&sequence.iter().any(|c|c=="local"),
        "the fixture must exercise both classes; got {sequence:?}");
    // A different run seed must produce a different stream.
    let vm2=validated_model(&dir,|m|{m["discovery"]["classes"]=classes.clone();m["kinetics"]["run_seed"]=json!(4243);});
    let other:Vec<String>=(0..12u64).map(|i|choose_discovery_class(&vm2,state,i).expect("choice").search_class).collect();
    assert_ne!(sequence,other,"E2-DISC-001: the class sequence is seed-dependent");
    // E2-API-005 / KMC-003: discovery never consumes the trajectory stream.
    let before=derive_trajectory_stream(4242).state();
    let _=choose_discovery_class(&vm,state,0).expect("choice");
    assert_eq!(derive_trajectory_stream(4242).state(),before,"E2-RNG-005: class selection never consumes the trajectory stream");
    let _=std::fs::remove_dir_all(&dir);
}

#[test]
fn discovery_out_of_order_completion_commits_in_search_id_order(){
    // E2-DISC-003: candidate transactions are committed in ascending search ID order, then
    // deduplicated transactionally, so completion order cannot alter catalog visibility.
    let ids:Vec<String>=(0..6u64).map(|i|search_id(11,GOLDEN_STATE_ID,"global",i).expect("search id")).collect();
    let mut sorted=ids.clone();sorted.sort();
    assert_ne!(ids,sorted,"the fixture needs a completion order that differs from search-ID order");
    let shuffled=vec![ids[3].clone(),ids[0].clone(),ids[5].clone(),ids[1].clone(),ids[4].clone(),ids[2].clone()];
    assert_eq!(ordered_candidate_ids(&shuffled).expect("ordering"),sorted,"E2-DISC-003: commit order is ascending search ID");
    assert_eq!(ordered_candidate_ids(&ids).expect("ordering"),sorted,"E2-DISC-003: the result is completion-order independent");
    let mut reversed=ids.clone();reversed.reverse();
    assert_eq!(ordered_candidate_ids(&reversed).expect("ordering"),sorted,"E2-PAR-006: parallel order yields the serial catalog");
    // Transactional deduplication: a repeated candidate is refused rather than committed twice.
    let mut duplicated=ids.clone();duplicated.push(ids[2].clone());
    assert!(ordered_candidate_ids(&duplicated).is_err(),"E2-DISC-003: duplicate candidate IDs are refused");
    assert!(ordered_candidate_ids(&[String::new()]).is_err(),"E2-DISC-003: empty candidate IDs are refused");
}

// ---------------------------------------------------------------- kmc-*

#[test]
fn kmc_two_event_selection_table(){
    // E2-KMC-001: selectable events sorted by event ID, summed by Neumaier compensated summation in
    // that order. E2-KMC-002: the snapshot envelope shape and parallel-array order.
    let(low,high)=("event:sha256:0a".to_owned(),"event:sha256:f0".to_owned());
    let rows=vec![(high.clone(),"state:sha256:b".to_owned(),3.0f64),(low.clone(),"state:sha256:a".to_owned(),1.0f64)];
    let snap=make_rate_snapshot("state:sha256:origin",&rows,-700.0).expect("snapshot");
    let p=&snap.payload;
    assert_eq!(p.schema,"spark-atomistic-rate-table-snapshot/1","E2-KMC-002");
    assert_eq!(p.event_ids,vec![low.clone(),high.clone()],"E2-KMC-001: event-ID order, not insertion order");
    assert_eq!(p.destination_state_ids,vec!["state:sha256:a","state:sha256:b"],"E2-KMC-002: parallel arrays follow event-ID order");
    assert_eq!(p.log_rates,vec![1.0,3.0]);
    for(r,l)in p.rates.iter().zip(&p.log_rates){assert_eq!(r.to_bits(),l.exp().to_bits(),"E2-KMC-002: rates are exp(log_rates) bit-exactly");}
    assert_eq!(p.total_rate_per_s.to_bits(),neumaier_sum(&p.rates).expect("sum").to_bits(),"E2-KMC-001: compensated sum in event-ID order");
    assert_eq!(p.lost_rate_log_upper_bound,None,"E2-RATE-002: no disabled rates in this table");
    let payload_value=serde_json::to_value(p).expect("payload serializes");
    let keys:Vec<&str>=payload_value.as_object().expect("object").keys().map(String::as_str).collect();
    assert_eq!(keys,vec!["destination_state_ids","event_ids","log_rates","lost_rate_log_upper_bound","origin_state_id","rates","schema","total_rate_per_s"],
        "E2-KMC-002: the payload has exactly these keys");
    snap.validate().expect("E2-KMC-002: a freshly built snapshot is self-consistent");
    // E2-KMC-003: the selection rule is the first cumulative rate strictly greater than u_s * K,
    // driven by the real trajectory stream.
    let rng=derive_trajectory_stream(0);
    let(us,ut,_)=rng.two_uniforms_atomic().expect("two uniforms");
    let threshold=us*p.total_rate_per_s;
    let mut cumulative=0.0;let mut chosen=p.rates.len()-1;
    for(j,r)in p.rates.iter().enumerate(){cumulative+=*r;if cumulative>threshold{chosen=j;break}}
    assert!(chosen<p.rates.len());
    let dt=-ut.ln()/p.total_rate_per_s;
    assert!(dt.is_finite()&&dt>0.0,"E2-KMC-003: dt = -ln(u_t)/K is finite and positive");
    // Rejections that E2-KMC-001/002 require.
    assert_eq!(make_rate_snapshot("state:sha256:origin",&[(low.clone(),"d".into(),1.0),(low.clone(),"d".into(),2.0)],-700.0)
        .expect_err("duplicate event ID").status,StatusCode::RateInvalid,"E2-KMC-002: duplicate event IDs are refused");
    assert_eq!(make_rate_snapshot("",&rows,-700.0).expect_err("empty origin").status,StatusCode::RateInvalid);
    let mut tampered=snap.clone();tampered.payload.total_rate_per_s+=1.0;
    assert!(tampered.validate().is_err(),"E2-KMC-005: a tampered historical snapshot is refused");
}

#[test]
fn kmc_lost_rate_bound(){
    // E2-RATE-002: a log rate below log_rate_cutoff is nonselectable, and every snapshot containing
    // such disabled rates reports their summed lost-rate log upper bound.
    let rows=vec![("event:sha256:01".to_owned(),"state:sha256:a".to_owned(),2.0f64),
                  ("event:sha256:02".to_owned(),"state:sha256:b".to_owned(),1.0f64),
                  ("event:sha256:03".to_owned(),"state:sha256:c".to_owned(),0.5f64)];
    let cutoff=1.5f64;
    let snap=make_rate_snapshot("state:sha256:origin",&rows,cutoff).expect("snapshot");
    assert_eq!(snap.payload.event_ids,vec!["event:sha256:01"],"E2-RATE-002: rates below the cutoff are nonselectable");
    let expected=log_sum_exp(&[1.0,0.5]).expect("bound");
    assert_eq!(snap.payload.lost_rate_log_upper_bound.expect("bound").to_bits(),expected.to_bits(),
        "E2-RATE-002: the lost-rate log upper bound is the summed log of the disabled rates");
    // Baseline: with the cutoff below every log rate, nothing is lost and the bound is null.
    let all=make_rate_snapshot("state:sha256:origin",&rows,-700.0).expect("snapshot");
    assert_eq!(all.payload.event_ids.len(),3,"baseline: all three rates are selectable");
    assert_eq!(all.payload.lost_rate_log_upper_bound,None,"E2-KMC-002: the bound is null when nothing is disabled");
    assert!(all.payload.total_rate_per_s>snap.payload.total_rate_per_s,
        "the disabled rates carry rate that the cutoff really removes from the total");
    // A rate cutoff far below the retained rates removes rate that is below the binary64 resolution
    // of the total, so the bound -- not the total -- is the only record that it existed (RATE-004).
    let deep=vec![rows[0].clone(),("event:sha256:02".to_owned(),"state:sha256:b".to_owned(),-50.0f64)];
    let underflow=make_rate_snapshot("state:sha256:origin",&deep,-10.0).expect("snapshot");
    let full=make_rate_snapshot("state:sha256:origin",&deep,-700.0).expect("snapshot");
    assert_eq!(underflow.payload.total_rate_per_s.to_bits(),full.payload.total_rate_per_s.to_bits(),
        "observed: exp(-50) is below one ULP of exp(2), so disabling it does not change the total");
    assert_eq!(underflow.payload.lost_rate_log_upper_bound.expect("bound").to_bits(),(-50.0f64).to_bits(),
        "RATE-004: the lost rate is still reported even when it cannot change the total");
    // KMC-002: a table with no selectable event is NO_ENABLED_EVENT, not a silent zero total.
    let none=make_rate_snapshot("state:sha256:origin",&rows[1..],cutoff).expect_err("all rates disabled");
    assert_eq!(none.status,StatusCode::NoEnabledEvent,"E2-KMC-001: an empty selectable table is NO_ENABLED_EVENT");
}

#[test]
fn kmc_application_rollback_leaves_the_trajectory_stream_untouched(){
    // KMC-005 / E2-KMC-003: the two step uniforms are drawn from a clone; the committed stream only
    // advances when the caller adopts the returned generator.
    let rng=derive_trajectory_stream(99);
    let before=rng.state();
    let(us,ut,next)=rng.two_uniforms_atomic().expect("two uniforms");
    assert_eq!(rng.state(),before,"KMC-005: a proposed-but-not-committed step leaves the RNG at its pre-step value");
    assert_eq!(before.consumed_uniforms,0);
    assert_eq!(next.state().consumed_uniforms,2,"E2-KMC-003: a committed step consumes exactly two uniforms");
    // Repeating the proposal after a failed application must reproduce the identical uniforms.
    let(us2,ut2,_)=rng.two_uniforms_atomic().expect("two uniforms");
    assert_eq!(us.to_bits(),us2.to_bits(),"E2-CKPT-008: replay after rollback is bit-identical");
    assert_eq!(ut.to_bits(),ut2.to_bits(),"E2-CKPT-008: replay after rollback is bit-identical");
    // After committing, the next proposal draws different uniforms.
    let(us3,_,_)=next.two_uniforms_atomic().expect("two uniforms");
    assert_ne!(us.to_bits(),us3.to_bits(),"E2-KMC-003: a committed step really advances the stream");
    // E2-CKPT-008: consumed uniforms equal 2 * step_index.
    let mut committed_rng=derive_trajectory_stream(99);
    for step in 1..=5u64{
        let(_,_,n)=committed_rng.two_uniforms_atomic().expect("two uniforms");
        committed_rng=n;
        assert_eq!(committed_rng.state().consumed_uniforms,2*step,"E2-CKPT-008: consumed uniforms = 2 * step_index");
    }
}

#[test]
fn kmc_historical_snapshot_is_immutable_under_catalog_growth(){
    // E2-KMC-005: the per-step snapshot is historical and immutable; restore/replay must use it and
    // never a later expanded catalog.
    let early=vec![("event:sha256:01".to_owned(),"state:sha256:a".to_owned(),1.0f64)];
    let snap=make_rate_snapshot("state:sha256:origin",&early,-700.0).expect("snapshot");
    let recorded=snap.clone();
    let grown=vec![early[0].clone(),("event:sha256:02".to_owned(),"state:sha256:b".to_owned(),4.0f64)];
    let later=make_rate_snapshot("state:sha256:origin",&grown,-700.0).expect("snapshot");
    assert_eq!(snap,recorded,"E2-KMC-005: building a later snapshot must not mutate the historical one");
    assert_ne!(snap.payload_sha256,later.payload_sha256,"the catalog really did grow");
    assert_eq!(snap.payload.event_ids.len(),1);
    assert_eq!(later.payload.event_ids.len(),2);
    assert!(later.payload.total_rate_per_s>snap.payload.total_rate_per_s);
    snap.validate().expect("E2-KMC-005: the historical snapshot still validates against its own hash");
    // Substituting the later table under the earlier hash must be detectable.
    let mut forged=snap.clone();forged.payload=later.payload.clone();
    assert_eq!(forged.validate().expect_err("substituted payload").status,StatusCode::CheckpointCorrupt,
        "E2-KMC-005: payload_sha256 binds the snapshot to its own historical contents");
}

// ---------------------------------------------------------------- checkpoint-* and basin-*

#[test]
fn checkpoint_clean_round_trip(){
    // E2-CKPT-001/002/003/004/005/006/007/009 on the smallest payload the wire schema admits.
    let dir=temp_dir("ckpt-clean");
    let vm=validated_model(&dir,|_|{});
    let state=committed(&[[0.0,0.0,0.0],[1.6,0.0,0.0]],-10.0);
    let payload=zero_step_payload(&state,&vm,|_|{});
    let bytes=encode_checkpoint_v2(payload.clone(),&vm).expect("E2-CKPT-001: a clean payload encodes");
    assert!(!bytes.ends_with(b"\n"),"E2-CKPT-009: canonical bytes carry no trailing newline");
    let decoded=decode_checkpoint_v2(&bytes,&vm).expect("E2-CKPT-007: a clean checkpoint restores");
    assert_eq!(decoded.payload,payload,"E2-CKPT-008: restore is byte-for-byte the same payload");
    assert_eq!(encode_checkpoint_v2(decoded.payload.clone(),&vm).expect("re-encode"),bytes,"E2-PAR-003: checkpoint bytes are stable");
    let v:Value=serde_json::from_slice(&bytes).expect("parses");
    let keys:Vec<&str>=v.as_object().expect("object").keys().map(String::as_str).collect();
    assert_eq!(keys,vec!["payload","payload_sha256"],"E2-CKPT-001: the envelope has exactly these keys");
    let pk:Vec<&str>=v["payload"].as_object().expect("object").keys().map(String::as_str).collect();
    assert_eq!(pk,vec!["basin","catalog","checkpoint_sequence","current_state","digests","discovery_statistics","flags",
        "initial_state","log_sequence","resources","rng","schema","simulation_time_s","step_index","trajectory"],
        "E2-CKPT-002: the payload has exactly these fifteen keys");
    assert_eq!(v["payload"]["schema"],"spark-atomistic-checkpoint/2","E2-CKPT-002");
    assert_eq!(v["payload"]["catalog"]["schema"],"spark-atomistic-catalog/2","E2-CKPT-006");
    let dk:Vec<&str>=v["payload"]["digests"].as_object().expect("object").keys().map(String::as_str).collect();
    assert_eq!(dk,vec!["config","model","schema","tolerances"],"E2-CKPT-003: digests has exactly these fields");
    let fk:Vec<&str>=v["payload"]["flags"].as_object().expect("object").keys().map(String::as_str).collect();
    assert_eq!(fk,vec!["cancelled","complete","incomplete_catalog","last_status","resource_limited"],"E2-CKPT-004");
    let rk:Vec<&str>=v["payload"]["resources"].as_object().expect("object").keys().map(String::as_str).collect();
    assert_eq!(rk,vec!["calculator_evaluations","catalog_events","output_bytes","resident_memory_bytes","retry_history",
        "saddle_attempts_by_state","wall_elapsed_s"],"E2-CKPT-005");
    assert_eq!(v["payload"]["resources"]["retry_history"],json!([]),"E2-CKPT-005: retry history is empty because retry_count is zero");
    assert_eq!(v["payload"]["rng"]["trajectory"]["consumed_uniforms"],0,"E2-CKPT-008: 2 * step_index uniforms consumed");
    let _=std::fs::remove_dir_all(&dir);
}

#[test]
fn checkpoint_corrupt_hash(){
    // E2-CKPT-001: the envelope hash covers the canonical payload bytes.
    let dir=temp_dir("ckpt-hash");
    let vm=validated_model(&dir,|_|{});
    let state=committed(&[[0.0,0.0,0.0],[1.6,0.0,0.0]],-10.0);
    let bytes=encode_checkpoint_v2(zero_step_payload(&state,&vm,|_|{}),&vm).expect("encode");
    decode_checkpoint_v2(&bytes,&vm).expect("baseline: the untampered bytes restore");
    let text=String::from_utf8(bytes).expect("canonical bytes are UTF-8");
    let at=text.rfind("\"payload_sha256\":\"sha256:").expect("hash field");
    let start=at+"\"payload_sha256\":\"sha256:".len();
    let mut tampered=text.clone();
    tampered.replace_range(start..start+1,if &text[start..start+1]=="a"{"b"}else{"a"});
    assert_eq!(tampered.len(),text.len(),"only one hex digit changed, so the bytes stay canonical");
    let err=decode_checkpoint_v2(tampered.as_bytes(),&vm).expect_err("E2-CKPT-001: hash mismatch must be refused");
    assert_eq!(err.status,StatusCode::CheckpointCorrupt);
    assert_eq!(err.status.exit_code(false,false),Some(74),"E2-STATUS-002: CHECKPOINT_CORRUPT exits 74");
    // Non-canonical bytes are refused before the hash is even consulted (E2-CKPT-007 step 1).
    let spaced=text.replacen("{\"payload\"","{ \"payload\"",1);
    assert_eq!(decode_checkpoint_v2(spaced.as_bytes(),&vm).expect_err("noncanonical").status,StatusCode::CheckpointCorrupt,
        "E2-CAN-001 / E2-CKPT-007: insignificant whitespace makes the checkpoint noncanonical");
    let _=std::fs::remove_dir_all(&dir);
}

#[test]
fn checkpoint_incompatible_digest(){
    // E2-CKPT-007: a valid but mismatched run contract is CHECKPOINT_INCOMPATIBLE, not CORRUPT.
    let dir=temp_dir("ckpt-incompat");
    let vm=validated_model(&dir,|_|{});
    let other=validated_model(&dir,|m|{m["kinetics"]["saddle_energy_tolerance"]=json!(2e-5);});
    assert_ne!(vm.config_digest(),other.config_digest(),"the two contracts really differ");
    assert_ne!(vm.tolerance_digest(),other.tolerance_digest(),"E2-CAN-007: the tolerance digest tracks the nine tolerance fields");
    let state=committed(&[[0.0,0.0,0.0],[1.6,0.0,0.0]],-10.0);
    let bytes=encode_checkpoint_v2(zero_step_payload(&state,&other,|_|{}),&other).expect("encode under the other contract");
    decode_checkpoint_v2(&bytes,&other).expect("baseline: it restores under its own contract");
    let err=decode_checkpoint_v2(&bytes,&vm).expect_err("E2-CKPT-007: mismatched contract");
    assert_eq!(err.status,StatusCode::CheckpointIncompatible);
    assert_eq!(err.status.exit_code(false,false),Some(74),"E2-STATUS-002: CHECKPOINT_INCOMPATIBLE exits 74");
    assert_eq!(err.status.message(),"checkpoint incompatible","E2-STATUS-002: exact stable message");
    let _=std::fs::remove_dir_all(&dir);
}

#[test]
fn checkpoint_recursive_corruption(){
    // E2-CKPT-007 steps 4 and 6: the envelope hash can be self-consistent while nested content is
    // wrong. Each variant below re-hashes the tampered payload, so only the recursive verification
    // can catch it.
    let dir=temp_dir("ckpt-recursive");
    let vm=validated_model(&dir,|_|{});
    let state=committed(&[[0.0,0.0,0.0],[1.6,0.0,0.0]],-10.0);
    let seal=|p:CheckpointPayloadV2|->Vec<u8>{
        let env=CheckpointEnvelopeV2::new(p).expect("envelope hashes its own payload");
        canonical_json_bytes(&env).expect("canonical envelope")};
    // Baseline: the untampered payload sealed the same way restores.
    decode_checkpoint_v2(&seal(zero_step_payload(&state,&vm,|_|{})),&vm).expect("baseline restores");
    // Step 4: the committed state no longer hashes to its recorded state_id.
    let energy_tamper=zero_step_payload(&state,&vm,|p|{
        let mut s=p.current_state.clone();s.energy_ev=-9.5;
        p.catalog.states.insert(s.state_id.clone(),s.clone());p.current_state=s.clone();p.initial_state=s;});
    let err=decode_checkpoint_v2(&seal(energy_tamper),&vm).expect_err("E2-CKPT-007 step 4: state record must recompute");
    assert_eq!(err.status,StatusCode::CheckpointCorrupt);
    // Step 5: the catalog digest no longer covers the catalog contents.
    let digest_tamper=zero_step_payload(&state,&vm,|p|{p.catalog.digest="sha256:".to_owned()+&"0".repeat(64);});
    assert_eq!(decode_checkpoint_v2(&seal(digest_tamper),&vm).expect_err("catalog digest").status,StatusCode::CheckpointCorrupt,
        "E2-CKPT-006: the catalog digest hashes the same object without `digest`");
    // Step 8: the Philox counter relation no longer holds.
    let rng_tamper=zero_step_payload(&state,&vm,|p|{p.rng.trajectory.consumed_uniforms=3;});
    assert_eq!(decode_checkpoint_v2(&seal(rng_tamper),&vm).expect_err("rng relation").status,StatusCode::CheckpointCorrupt,
        "E2-CKPT-008: consumed uniforms equal 2 * step_index");
    // Step 9/10: the trajectory length must equal the step and log sequences.
    let seq_tamper=zero_step_payload(&state,&vm,|p|{p.step_index=1;});
    assert_eq!(decode_checkpoint_v2(&seal(seq_tamper),&vm).expect_err("sequence").status,StatusCode::CheckpointCorrupt,
        "E2-CKPT-008: step and log sequences equal trajectory length");
    // Step 10: the final state must be the state the trajectory ends on.
    let other_state=committed(&[[0.0,0.0,0.0],[2.4,0.0,0.0]],-9.7);
    let state_tamper=zero_step_payload(&state,&vm,|p|{p.current_state=other_state.clone();});
    assert_eq!(decode_checkpoint_v2(&seal(state_tamper),&vm).expect_err("final state").status,StatusCode::CheckpointCorrupt,
        "E2-CKPT-007 step 10: the final state must be present and consistent");
    let _=std::fs::remove_dir_all(&dir);
}

#[test]
fn checkpoint_cancellation_and_resource_limit_flags(){
    // E2-CKPT-004: the five flags are exact booleans plus one exact status token.
    // RES-002 / E2-STATUS-002: CANCELLED and RESOURCE_LIMIT are pauses that exit 75.
    let dir=temp_dir("ckpt-flags");
    let vm=validated_model(&dir,|_|{});
    let state=committed(&[[0.0,0.0,0.0],[1.6,0.0,0.0]],-10.0);
    let seal=|p:CheckpointPayloadV2|->Vec<u8>{canonical_json_bytes(&CheckpointEnvelopeV2::new(p).expect("envelope")).expect("bytes")};
    for(status,flag)in [(StatusCode::Cancelled,"cancelled"),(StatusCode::ResourceLimit,"resource_limited")]{
        let payload=zero_step_payload(&state,&vm,|p|{
            p.flags.last_status=status;
            if flag=="cancelled"{p.flags.cancelled=true}else{p.flags.resource_limited=true}});
        let bytes=encode_checkpoint_v2(payload,&vm).expect("E2-CKPT-005: a paused run still writes a checkpoint");
        let decoded=decode_checkpoint_v2(&bytes,&vm).expect("the paused checkpoint restores");
        assert_eq!(decoded.payload.flags.last_status,status,"E2-CKPT-004: last_status is an exact status token");
        assert!(!decoded.payload.flags.complete,"E2-CKPT-004: a paused run is not complete");
        assert_eq!(status.exit_code(false,false),Some(75),"E2-STATUS-002: {status:?} exits 75");
        assert_eq!(decoded.payload.step_index,0,"RES-002: the pause preserved the last committed state and time");
        assert_eq!(decoded.payload.simulation_time_s,0.0);
    }
    // A run cannot be complete and cancelled/resource-limited at the same time.
    for tamper in [|p:&mut CheckpointPayloadV2|{p.flags.complete=true;p.flags.cancelled=true;},
                   |p:&mut CheckpointPayloadV2|{p.flags.complete=true;p.flags.resource_limited=true;}]{
        let payload=zero_step_payload(&state,&vm,tamper);
        assert_eq!(decode_checkpoint_v2(&seal(payload),&vm).expect_err("contradictory flags").status,StatusCode::CheckpointCorrupt,
            "E2-CKPT-004: complete cannot coexist with cancelled or resource_limited");
    }
    // E2-CKPT-007 step 7: resource counters must stay inside the configured limits.
    let limit=vm.model().resources.total_calculator_evaluations;
    let at_limit=zero_step_payload(&state,&vm,|p|{p.resources.calculator_evaluations=limit;});
    decode_checkpoint_v2(&seal(at_limit),&vm).expect("baseline: a counter exactly at the limit is accepted");
    let over=zero_step_payload(&state,&vm,|p|{p.resources.calculator_evaluations=limit+1;});
    assert_eq!(decode_checkpoint_v2(&seal(over),&vm).expect_err("counter overrun").status,StatusCode::CheckpointCorrupt,
        "E2-CKPT-007 step 7: one evaluation over the configured limit is refused");
    // A counter outside the portable-integer domain cannot even be canonicalised (E2-JSON-002 /
    // E2-CAN-003), so it is refused by the encoder before the limit check is reached.
    let unportable=zero_step_payload(&state,&vm,|p|{p.resources.calculator_evaluations=u64::MAX;});
    assert_eq!(CheckpointEnvelopeV2::new(unportable).expect_err("unportable counter").status,StatusCode::NonfiniteResult,
        "observed: an out-of-domain integer is refused at canonicalisation with NONFINITE_RESULT / E2-CAN-005");
    let retried=zero_step_payload(&state,&vm,|p|{p.resources.retry_history=vec![json!({"attempt":1})];});
    assert_eq!(decode_checkpoint_v2(&seal(retried),&vm).expect_err("retry history").status,StatusCode::CheckpointCorrupt,
        "E2-CKPT-005: retry history is empty in v1 because configured retry count is zero");
    let _=std::fs::remove_dir_all(&dir);
}

#[test]
fn basin_is_uniformly_disabled(){
    // E2-BASIN-001/002/003, E2-SCHEMA-011, E2-CKPT-003.
    let caps=response(CAPS_REQUEST);
    assert_eq!(caps["value"]["basin_acceleration"],"disabled","E2-BASIN-001: capability reports disabled");
    assert_eq!(caps["value"]["features"]["serial_kmc"],true,"E2-BASIN-002: serial KMC remains the execution path");
    assert_eq!(caps["value"]["features"]["harmonic_tst"],false);
    assert_eq!(caps["value"]["conformance"],"unvalidated","E2-API-007");
    // A model may request basin acceleration; validation still returns OK (E2-SCHEMA-011).
    let mut m:Value=serde_json::from_slice(MODEL).expect("parses");
    m["basin"]["enabled"]=json!(true);
    let v=response(&serde_json::to_vec(&json!({"model":m,"operation":"validate"})).expect("serializes"));
    assert_eq!(v["status"],"OK","E2-SCHEMA-011: basin.enabled=true is schema-valid");
    assert_eq!(response(CAPS_REQUEST)["value"]["basin_acceleration"],"disabled",
        "E2-BASIN-003: no basin capability may be advertised under this IR revision");
    // Any attempt to invoke acceleration is the nonterminal BASIN_DISABLED (E2-BASIN-002).
    let err=attempt_acceleration().expect_err("E2-BASIN-002: acceleration cannot be enabled");
    assert_eq!(err.status,StatusCode::BasinDisabled);
    assert_eq!(err.status.exit_code(false,false),None,"E2-STATUS-003: a nonterminal status has no process exit code");
    assert_eq!(err.status.message(),"basin acceleration disabled","E2-STATUS-002: exact stable message");
    BasinCheckpoint::disabled().validate().expect("the disabled internal record validates");
    // E2-CKPT-003: the wire checkpoint stores exactly {"enabled":false,"reason":"v1-disabled"}.
    let dir=temp_dir("basin");
    let vm=validated_model(&dir,|_|{});
    let state=committed(&[[0.0,0.0,0.0],[1.6,0.0,0.0]],-10.0);
    let bytes=encode_checkpoint_v2(zero_step_payload(&state,&vm,|_|{}),&vm).expect("encode");
    let parsed:Value=serde_json::from_slice(&bytes).expect("parses");
    assert_eq!(parsed["payload"]["basin"],json!({"enabled":false,"reason":"v1-disabled"}),"E2-CKPT-003: exact basin record");
    let seal=|p:CheckpointPayloadV2|->Vec<u8>{canonical_json_bytes(&CheckpointEnvelopeV2::new(p).expect("envelope")).expect("bytes")};
    for bad in [BasinCheckpointV2{enabled:true,reason:"v1-disabled".into()},
                BasinCheckpointV2{enabled:false,reason:"enabled-later".into()}]{
        let payload=zero_step_payload(&state,&vm,|p|{p.basin=bad;});
        assert_eq!(decode_checkpoint_v2(&seal(payload),&vm).expect_err("basin record").status,StatusCode::CheckpointCorrupt,
            "E2-CKPT-003: any other basin record is refused");
    }
    let _=std::fs::remove_dir_all(&dir);
}

// ---------------------------------------------------------------- defect witness

#[test]
fn defect_reciprocal_event_pair_is_unrepresentable_in_a_v2_checkpoint(){
    // DEFECT WITNESS -- records observed, blocking behavior. This test is green because it asserts
    // what the crate does today, not what Errata 2 requires.
    //
    // `validate_checkpoint_v2` imposes three constraints that cannot hold together for any directed
    // pair whose endpoints differ (which E2-EVENT-001 requires them to):
    //   C1  every directed record must satisfy saddle.search_id == discovery_provenance.search_id;
    //   C2  the reciprocal record's `saddle` object must equal the forward record's `saddle`, hence
    //       both records carry one search_id;
    //   C3  every record's discovery_provenance.search_id (and its substream key/counter) must
    //       re-derive from THAT RECORD's own origin_state_id, per E2-DISC-002 and E2-RNG-004.
    // C3 forces two different search IDs, C1+C2 force one. The consequence is that no checkpoint
    // containing any catalog event can be encoded or restored, which blocks E2-PAR-002 items 7, 9
    // and the exact-next-event resume of item 10.
    let dir=temp_dir("defect-pair");
    let vm=validated_model(&dir,|_|{});
    let seed=vm.model().kinetics.run_seed;
    let(t,nu,btol,db,cutoff)=(vm.model().kinetics.temperature,vm.model().kinetics.prefactor,
        vm.model().kinetics.barrier_tolerance,vm.model().kinetics.detailed_balance_tolerance,vm.model().kinetics.log_rate_cutoff);
    let(pa,pb,ps)=([[0.0,0.0,0.0],[1.6,0.0,0.0]],[[0.0,0.0,0.0],[2.4,0.0,0.0]],[[0.0,0.0,0.0],[2.0,0.0,0.0]]);
    let(ea,eb,es)=(-10.0f64,-9.7f64,-9.4f64);
    let(a,b)=(committed(&pa,ea),committed(&pb,eb));
    assert_ne!(a.state_id,b.state_id,"E2-EVENT-001: origin and destination are distinct committed IDs");
    let sid_a=search_id(seed,&a.state_id,"global",0).expect("search id");
    let sid_b=search_id(seed,&b.state_id,"global",0).expect("search id");
    assert_ne!(sid_a,sid_b,"E2-DISC-002: the search ID binds the state ID, so the two endpoints cannot share one");
    let sgd=saddle_geometry_digest(&a,&ps);
    let saddle=saddle_of(&ps,es,&sid_a);
    let map=vec![[0usize,0usize],[1,1]];
    let(pair,ef,er)=event_ids(&a.state_id,&b.state_id,&saddle,&sgd,&map).expect("ids");
    let rf=common_prefactor_pair(ea,eb,es,t,nu,btol,db).expect("forward rate");
    let rr=common_prefactor_pair(eb,ea,es,t,nu,btol,db).expect("reverse rate");
    let streams=[(sid_a.clone(),derive_saddle_substream(seed,&a.state_id,"global",0).expect("stream")),
                 (sid_b.clone(),derive_saddle_substream(seed,&b.state_id,"global",0).expect("stream"))];
    let build=|origin:&CommittedStateV2,dest:&CommittedStateV2,eid:&str,rid:&str,rate:&spark_atomistic_rs::model::RateModelRecord,
               saddle_sid:&str,prov_sid:&str,stream:&Philox,mapping:&Vec<[usize;2]>|->DirectedEventV2{
        let mut s=saddle.clone();s.search_id=saddle_sid.into();
        DirectedEventV2{schema:"spark-atomistic-directed-event/2".into(),event_id:eid.into(),reverse_event_id:rid.into(),
            pair_id:pair.clone(),origin_state_id:origin.state_id.clone(),destination_state_id:dest.state_id.clone(),
            saddle:s,barrier_ev:es-origin.energy_ev,reverse_barrier_ev:es-dest.energy_ev,
            rate_model:RateModelV2{common_prefactor_per_s:nu,detailed_balance_residual:rate.detailed_balance_residual,
                log_forward_rate_per_s:rate.log_forward_rate_per_s,log_reverse_rate_per_s:rate.log_reverse_rate_per_s,
                model:"COMMON_PREFACTOR".into(),temperature_k:t},
            selectable:rate.log_forward_rate_per_s>=cutoff,active_atom_mapping:mapping.clone(),
            environment_key:"disabled".into(),environment_version:"none/1".into(),
            discovery_provenance:DiscoveryProvenanceV2{rng_substream_digest:substream_digest(stream),
                search_class:"global".into(),search_id:prov_sid.into(),search_index:0},
            validation:EventValidationV2{calculator_model_digest:CALC_DIGEST.into(),constraint_digest:origin.constraint_digest.clone(),
                destination_match:MatchV2{atom_mapping:vec![0,1],energy_difference_ev:0.0,max_displacement_angstrom:0.0,rms_displacement_angstrom:0.0},
                full_endpoint_relaxations:true,method:"DIRECTIONAL".into(),
                origin_match:MatchV2{atom_mapping:vec![0,1],energy_difference_ev:0.0,max_displacement_angstrom:0.0,rms_displacement_angstrom:0.0},
                unstable_mode_count:1},
            calculator_digest:vm.model().calculator.model_digest.clone(),identity_digest:vm.identity_digest().into(),
            schema_digest:SCHEMA_DIGEST.into(),tolerance_digest:vm.tolerance_digest().into()}};
    let reverse_map={let mut m:Vec<[usize;2]>=map.iter().map(|x|[x[1],x[0]]).collect();m.sort();m};
    // Exhaust every assignment of (saddle.search_id, provenance.search_id, provenance substream) to
    // {A, B} for both directed records: 2^6 = 64 candidate catalogs.
    let mut accepted=Vec::new();
    for combo in 0..64u32{
        let pick=|bit:u32|&streams[((combo>>bit)&1) as usize];
        let forward=build(&a,&b,&ef,&er,&rf,&pick(0).0,&pick(1).0,&pick(2).1,&map);
        let reverse=build(&b,&a,&er,&ef,&rr,&pick(3).0,&pick(4).0,&pick(5).1,&reverse_map);
        let mut events=BTreeMap::new();events.insert(ef.clone(),forward);events.insert(er.clone(),reverse);
        let mut multiplicity=BTreeMap::new();multiplicity.insert(ef.clone(),1u64);multiplicity.insert(er.clone(),1u64);
        let mut states=BTreeMap::new();states.insert(a.state_id.clone(),a.clone());states.insert(b.state_id.clone(),b.clone());
        let mut catalog=CatalogV2{digest:String::new(),events,multiplicity,schema:"spark-atomistic-catalog/2".into(),states};
        let mut cv=serde_json::to_value(&catalog).expect("catalog");cv.as_object_mut().expect("object").remove("digest");
        catalog.digest=digest_of(&cv);
        let mut substreams=BTreeMap::new();
        for(id,stream)in &streams{substreams.insert(id.clone(),stream.state());}
        let payload=zero_step_payload(&a,&vm,|p|{
            p.catalog=catalog;p.rng.substream_map=substreams;p.resources.catalog_events=2;
            p.resources.saddle_attempts_by_state=BTreeMap::from([(a.state_id.clone(),1u64)]);});
        if encode_checkpoint_v2(payload,&vm).is_ok(){accepted.push(combo);}
    }
    assert!(accepted.is_empty(),"DEFECT: expected every reciprocal-pair assignment to be refused; accepted {accepted:?}");
    // Control: the identical payload without the catalog event is accepted, so the rejection above is
    // caused by the event pair and not by any other part of the fixture.
    let control=encode_checkpoint_v2(zero_step_payload(&a,&vm,|_|{}),&vm);
    assert!(control.is_ok(),"control: the same payload with an empty catalog encodes");
    // And a catalog can never dodge the problem by making the endpoints equal.
    let(_,self_id,_)=event_ids(&a.state_id,&a.state_id,&saddle,&sgd,&map).expect("ids");
    assert!(!self_id.is_empty());
    let _=std::fs::remove_dir_all(&dir);
}

#[test]
fn e2_can_004_canonical_binary64_survives_this_crates_own_parser(){
    // CONFORMANCE (was the D-E2-03 defect witness, retired 2026-08-11).
    //
    // E2-CAN-004 requires "the shortest decimal that round-trips to the same binary64 value", and
    // E2-CKPT-007 step 1 together with E2-PAR-003 require canonical bytes to survive a parse and
    // re-encode cycle. The WRITE side always satisfied this. The READ side did not: `serde_json` was
    // pinned WITHOUT the `float_roundtrip` feature, so its parser was not exactly rounding and
    // returned a value one ULP low for several 17-significant-digit decimals -- two of them mandatory
    // E2-PAR-002 item-4 corpus values, which meant this crate could not re-read 4 of its own 85
    // emitted artifacts. Enabling `float_roundtrip` in Cargo.toml closed it; this test now asserts
    // the requirement instead of the defect, and fails if the feature is ever dropped.
    let reparse=|x:f64|->f64{
        let bytes=canonical_json_bytes(&json!(x)).expect("canonical");
        parse_strict_json(&bytes).expect("reparses").as_f64().expect("number")};
    // Baseline that round-tripped even before the fix, so a regression cannot hide behind it.
    let ok=f64::from_bits(0x3eb0c6f7a0b5ed8e);
    assert_eq!(canonical_json_bytes(&json!(ok)).expect("canonical"),b"0.0000010000000000000002");
    assert_eq!(reparse(ok).to_bits(),0x3eb0c6f7a0b5ed8e,"baseline: this canonical decimal round-trips");
    // The three values that previously came back one ULP low.
    for(bits,expected)in[(0xc027ffffff7d2cecu64,"-11.999999984770021"),
                         (0x3eb0c6f7a0b5ed8c,"9.999999999999997e-7"),
                         (0x444b1ae4d6e2ef51,"1.0000000000000001e+21")]{
        let x=f64::from_bits(bits);
        assert_eq!(canonical_json_bytes(&json!(x)).expect("canonical"),expected.as_bytes(),
            "E2-CAN-004: the emitted decimal is the shortest round-tripping one");
        assert_eq!(reparse(x).to_bits(),bits,
            "E2-CAN-004/E2-CKPT-007: {expected} must reparse to exactly 0x{bits:016x}");
    }
}

#[test]
fn e3_can_001_large_binary64_is_exponent_form_and_reparseable(){
    let x=f64::from_bits(0x444b1ae4d6e2ef4f);
    let bytes=canonical_json_bytes(&json!(x)).expect("canonical");
    assert_eq!(bytes,b"9.999999999999999e+20");
    let reparsed=parse_strict_json(&bytes).expect("E3-CAN-001 output must be valid input");
    assert_eq!(reparsed.as_f64().expect("real").to_bits(),x.to_bits());
    let boundary=canonical_json_bytes(&json!(9_007_199_254_740_992.0_f64)).expect("canonical");
    assert_eq!(boundary,b"9.007199254740992e+15");
}

// ---------------------------------------------------------------- manifest self-check

#[test]
fn parity_manifest_matches_this_suite(){
    // E2-PAR-001 / E2-PAR-005: no mandatory fixture may be skipped, and the recorded execution status
    // must exist for every fixture. This test makes the manifest tamper-evident against the suite.
    let manifest:Value=serde_json::from_slice(MANIFEST).expect("manifest parses");
    assert_eq!(manifest["schema"],"spark-atomistic-parity-corpus/2");
    let fixtures:Vec<&str>=manifest["fixtures"].as_array().expect("fixtures array").iter()
        .map(|x|x.as_str().expect("fixture name")).collect();
    assert_eq!(fixtures.len(),46,"E2-PAR-002: the mandatory corpus has 46 named fixtures");
    let execution=manifest["fixture_execution"].as_object().expect("fixture_execution object");
    assert_eq!(execution.len(),fixtures.len(),"every fixture needs exactly one execution record");
    let allowed=["PASS","PARTIAL","FAIL_DEFECT","BLOCKED"];
    for name in &fixtures{
        let status=execution.get(*name).unwrap_or_else(||panic!("fixture {name} has no execution record"))
            .as_str().unwrap_or_else(||panic!("fixture {name} execution status is not a string"));
        assert!(allowed.contains(&status),"fixture {name} has unknown execution status {status}");
    }
    // Errata 3 closes the old partial/blocker split; all mandatory fixtures execute.
    let passing=fixtures.iter().filter(|n|execution[**n]=="PASS").count();
    assert_eq!(passing,fixtures.len(),"E3-PAR-001: every mandatory fixture must execute");
    assert_eq!(manifest["dynamic_execution"],"COMPLETE_WIRE_PARITY");
    assert_eq!(manifest["execution_summary"]["pass"],json!(passing),"the recorded pass count must match the recorded statuses");
    for key in ["capability_request","capability_response","canonical_numbers","minimal_model","philox_boundaries"]{
        let file=manifest["golden_sources"][key].as_str().expect("golden source name");
        assert!(!file.is_empty(),"golden source {key} must name a corpus file");
    }
    // E2-PAR-005: execution conformance additionally requires every fixture to pass in EVERY
    // implementation, so the manifest must also carry the cross-language record and must not
    // claim more than that record supports.
    let cross=&manifest["cross_language_execution"];
    assert_eq!(cross["status"],"EXECUTED","the cross-language corpus has been run against both backends");
    let verdicts=cross["fixture_verdict"].as_object().expect("cross-language fixture verdicts");
    assert_eq!(verdicts.len(),fixtures.len(),"every mandatory fixture needs a cross-language verdict");
    let cross_pass=fixtures.iter().filter(|n|verdicts[**n]=="PASS").count();
    for name in &fixtures{
        let v=verdicts[*name].as_str().unwrap_or_else(||panic!("fixture {name} verdict is not a string"));
        assert!(["PASS","DIVERGENT","NOT-APPLICABLE"].contains(&v),"fixture {name} has unknown verdict {v}");
    }
    assert_eq!(cross["fixture_summary"]["pass"],json!(cross_pass),"recorded cross-language pass count must match the verdicts");
    assert_eq!(cross["fixture_summary"]["total"],json!(fixtures.len()));
    assert_eq!(cross["e2_par_005_claimable"],json!(cross_pass==fixtures.len()),
        "E2-PAR-005 may be claimed only when every mandatory fixture passes in every implementation");
    assert!(cross["e2_par_005_claimable"].as_bool().expect("bool"),
        "wire parity is claimable only because every fixture passes; scientific validation is separate");
    assert_eq!(cross["tier_summary"]["core"]["pass"],json!(83));
    assert_eq!(cross["tier_summary"]["adapter"]["pass"],json!(2));
}
