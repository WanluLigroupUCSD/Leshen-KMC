// Clean-room cross-language parity emitter, authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 +
// ERRATA_2_PARITY. Independently authored; no implementation source consulted.
//
// `E2-PAR-001` requires every conforming implementation to consume the SAME canonical request corpus
// and emit the same canonical response, status record, state/event/catalog record, rate snapshot,
// RNG state and checkpoint schema; `E2-PAR-003` requires those bytes to be identical. This test is
// the Rust half of that experiment: it drives `tests/corpus/xlang/` and writes one canonical byte
// string per case. `spark-atomistic/tests/xlang_harness.py` is the Python half and the comparator.
//
// It is a real test, not only an emitter: every case must produce canonical bytes that re-parse and
// re-canonicalise to themselves (E2-CAN-001), and the corpus manifest must name every case it emits.
use serde_json::{json, Map, Value};
use spark_atomistic_rs::checkpoint::{canonical_json_bytes, parse_strict_json};
use spark_atomistic_rs::parity::*;
use spark_atomistic_rs::rate::common_prefactor_pair;
use spark_atomistic_rs::rng::{
    derive_saddle_substream, derive_trajectory_stream, philox4x32_10, substream_digest,
    uniform_from_words, Philox,
};
use spark_atomistic_rs::run::ProcessRunAdapter;
use spark_atomistic_rs::status::StatusCode;
use std::path::{Path, PathBuf};

const ALL_STATUS: [StatusCode; 29] = [
    StatusCode::Ok, StatusCode::DiscoveryConvergedHeuristic, StatusCode::DuplicateEvent,
    StatusCode::SaddleNotFound, StatusCode::InvalidSaddle, StatusCode::SaddleWrongBasin,
    StatusCode::EndpointCollapsed, StatusCode::EnvironmentAmbiguous, StatusCode::BasinDisabled,
    StatusCode::DiscoveryIncomplete, StatusCode::RelaxNotConverged,
    StatusCode::EventApplicationFailed, StatusCode::CalculatorFailure, StatusCode::NonfiniteResult,
    StatusCode::InvalidInput, StatusCode::SchemaUnsupported, StatusCode::InvalidState,
    StatusCode::RateInvalid, StatusCode::DetailedBalanceViolation, StatusCode::CatalogConflict,
    StatusCode::CatalogIncompatible, StatusCode::AtomCountChangeUnsupported,
    StatusCode::NoEnabledEvent, StatusCode::ResourceLimit, StatusCode::OutputExists,
    StatusCode::CheckpointCorrupt, StatusCode::CheckpointIncompatible, StatusCode::Cancelled,
    StatusCode::InternalError,
];

fn corpus() -> PathBuf { Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/corpus/xlang") }
fn out_dir() -> PathBuf {
    let base = std::env::var_os("SPARK_XLANG_OUT")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("spark-xlang-out"));
    base.join("rust")
}
fn canon(v: &Value) -> Vec<u8> { canonical_json_bytes(v).expect("emitted value is canonical") }
fn num(x: f64) -> String { String::from_utf8(canon(&json!(x))).expect("canonical number is UTF-8") }
fn bits(x: f64) -> String { format!("0x{:016x}", x.to_bits()) }
fn scalar(x: f64) -> Value { json!({"bits": bits(x), "canonical": num(x)}) }
fn failure(status: StatusCode, requirement: &str) -> Value {
    json!({"outcome": "failure", "requirement_id": requirement, "status": status})
}
fn from_api(e: &ApiFailure) -> Value { failure(e.status, &e.requirement_id) }
fn from_status(s: &spark_atomistic_rs::status::Status) -> Value {
    failure(s.status, &s.context.requirement_id)
}

fn validated(model: &Value) -> Result<ValidatedModel, ApiFailure> {
    let wire: WireModel = serde_json::from_value(model.clone())
        .map_err(|_| ApiFailure::new(StatusCode::InvalidInput, "E2-SCOPE-003"))?;
    wire.validate(None)
}

fn philox_state_value(p: &Philox) -> Value {
    serde_json::to_value(p.state()).expect("Philox state serialises")
}

fn draw(mut stream: Philox, draws: u64) -> Value {
    let mut uniforms = Vec::new();
    for _ in 0..draws {
        match stream.next_uniform() {
            Ok(u) => uniforms.push(scalar(u)),
            Err(e) => { uniforms.push(from_status(&e)); break }
        }
    }
    json!({"final_state": philox_state_value(&stream), "substream_digest": substream_digest(&stream),
           "uniforms": uniforms})
}

fn saddle_stub(energy_ev: f64, mode: &Value) -> SaddleV2 {
    SaddleV2 {
        curvature_ev_per_angstrom2: -1.0, energy_ev, evaluation_count: 0,
        evidence_level: "DIRECTIONAL".into(), forces_ev_per_angstrom: Vec::new(),
        orthogonal_curvatures_ev_per_angstrom2: Vec::new(), positions: Vec::new(),
        search_id: String::new(), termination_reason: "stub".into(),
        unstable_direction: serde_json::from_value(mode.clone()).expect("mode vectors"),
    }
}

fn rows_of(input: &Value) -> Vec<(String, String, f64)> {
    input["rows"].as_array().expect("rows array").iter().map(|r| {
        (r[0].as_str().expect("event id").to_owned(),
         r[1].as_str().expect("destination id").to_owned(),
         r[2].as_f64().expect("log rate"))
    }).collect()
}

fn snapshot_value(s: &RateTableSnapshot) -> Value {
    serde_json::to_value(s).expect("rate-table snapshot serialises")
}

fn probe(input: &Value, model_template: &Value) -> Value {
    match input["probe"].as_str().expect("probe name") {
        "canonical_number" => {
            let out: Vec<Value> = input["bits"].as_array().expect("bits array").iter().map(|b| {
                let raw = b.as_str().expect("hex bits");
                let x = f64::from_bits(u64::from_str_radix(&raw[2..], 16).expect("hex"));
                json!({"canonical": num(x), "input_bits": raw, "round_trip_bits": bits(x)})
            }).collect();
            json!({"cases": out})
        }
        "digests" => {
            let out: Vec<Value> = input["models"].as_array().expect("models").iter().map(|m| {
                match validated(m) {
                    Ok(v) => json!({"config_digest": v.config_digest(), "identity_digest": v.identity_digest(),
                                    "outcome": "validated", "schema_digest": SCHEMA_DIGEST,
                                    "tolerance_digest": v.tolerance_digest()}),
                    Err(e) => from_api(&e),
                }
            }).collect();
            json!({"cases": out})
        }
        "state_identity" => {
            let out: Vec<Value> = input["cases"].as_array().expect("cases").iter().map(|c| {
                let system: WireSystem = serde_json::from_value(c["system"].clone()).expect("system");
                let energy = c["energy_ev"].as_f64().expect("energy");
                match (fixed_contract_digest(&system), constraint_digest(&system),
                       geometry_certificate(&system), state_ids(&system, energy)) {
                    (Ok(fixed), Ok(constraint), Ok(certificate), Ok((candidate, state))) =>
                        json!({"candidate_identity": candidate, "constraint_digest": constraint,
                               "energy_ev": scalar(energy), "fixed_contract_digest": fixed,
                               "geometry_certificate": certificate, "state_id": state}),
                    _ => failure(StatusCode::InvalidState, "E2-ID-002"),
                }
            }).collect();
            json!({"cases": out})
        }
        "philox_block" => {
            let out: Vec<Value> = input["blocks"].as_array().expect("blocks").iter().map(|b| {
                let key: [u32; 2] = serde_json::from_value(b["key"].clone()).expect("key words");
                let counter: [u32; 4] = serde_json::from_value(b["counter"].clone()).expect("counter words");
                json!({"counter": counter, "key": key, "words": philox4x32_10(counter, key)})
            }).collect();
            json!({"cases": out})
        }
        "uniform_words" => {
            let out: Vec<Value> = input["pairs"].as_array().expect("pairs").iter().map(|p| {
                let a = p[0].as_u64().expect("a") as u32;
                let b = p[1].as_u64().expect("b") as u32;
                let q = ((a as u64) << 20) | ((b as u64) >> 12);
                match uniform_from_words(a, b) {
                    Ok(u) => json!({"a": a, "b": b, "q": q, "raw_binary64_bits": bits(u),
                                    "uniform": num(u)}),
                    Err(e) => from_status(&e),
                }
            }).collect();
            json!({"cases": out})
        }
        "rng_derivation" => {
            let draws = input["draws"].as_u64().expect("draws");
            let trajectory: Vec<Value> = input["trajectory_seeds"].as_array().expect("seeds").iter()
                .map(|s| {
                    let seed = s.as_u64().expect("seed");
                    let stream = derive_trajectory_stream(seed);
                    json!({"initial_state": philox_state_value(&stream), "run_seed": seed,
                           "stream": draw(stream, draws)})
                }).collect();
            let saddle: Vec<Value> = input["saddle"].as_array().expect("saddle entries").iter()
                .map(|e| {
                    let seed = e["run_seed"].as_u64().expect("run_seed");
                    let state = e["state_id"].as_str().expect("state_id");
                    let class = e["search_class"].as_str().expect("search_class");
                    let index = e["search_index"].as_u64().expect("search_index");
                    match derive_saddle_substream(seed, state, class, index) {
                        Ok(stream) => json!({"entry": e.clone(),
                                             "initial_state": philox_state_value(&stream),
                                             "stream": draw(stream, draws)}),
                        Err(err) => from_status(&err),
                    }
                }).collect();
            json!({"saddle": saddle, "trajectory": trajectory})
        }
        "rng_state_sequence" => {
            let key: [u32; 2] = serde_json::from_value(input["key"].clone()).expect("key");
            let counter: [u32; 4] = serde_json::from_value(input["initial_counter"].clone()).expect("counter");
            let mut stream = Philox::new(key, counter);
            let mut steps = Vec::new();
            for _ in 0..input["draws"].as_u64().expect("draws") {
                let before = philox_state_value(&stream);
                match stream.next_uniform() {
                    Ok(u) => steps.push(json!({"state_after": philox_state_value(&stream),
                                               "state_before": before, "uniform": scalar(u)})),
                    Err(e) => { steps.push(from_status(&e)); break }
                }
            }
            json!({"final_substream_digest": substream_digest(&stream), "steps": steps})
        }
        "search_ids" => {
            let out: Vec<Value> = input["entries"].as_array().expect("entries").iter().map(|e| {
                let seed = e["run_seed"].as_u64().expect("run_seed");
                let state = e["state_id"].as_str().expect("state_id");
                let class = e["search_class"].as_str().expect("search_class");
                let index = e["search_index"].as_u64().expect("search_index");
                match search_id(seed, state, class, index) {
                    Ok(id) => json!({"entry": e.clone(), "search_id": id}),
                    Err(err) => from_api(&err),
                }
            }).collect();
            // E2-DISC-003: commit order is ascending search ID, independent of completion order.
            let mut ordered: Vec<String> = out.iter().filter_map(|x| x["search_id"].as_str().map(str::to_owned)).collect();
            ordered.sort();
            json!({"ascending_commit_order": ordered, "cases": out})
        }
        "discovery_class" => {
            let mut model = model_template.clone();
            model["discovery"]["classes"] = input["classes"].clone();
            model["kinetics"]["run_seed"] = input["run_seed"].clone();
            let vm = match validated(&model) { Ok(v) => v, Err(e) => return from_api(&e) };
            let state = input["state_id"].as_str().expect("state_id");
            let out: Vec<Value> = input["indices"].as_array().expect("indices").iter().map(|i| {
                let index = i.as_u64().expect("index");
                match choose_discovery_class(&vm, state, index) {
                    Ok(c) => json!({"class_stream": serde_json::to_value(&c.class_stream).expect("state"),
                                    "saddle_stream": serde_json::to_value(&c.saddle_stream).expect("state"),
                                    "search_class": c.search_class, "search_id": c.search_id,
                                    "search_index": index}),
                    Err(e) => from_api(&e),
                }
            }).collect();
            json!({"cases": out})
        }
        "event_ids" => {
            let out: Vec<Value> = input["cases"].as_array().expect("cases").iter().map(|c| {
                let saddle = saddle_stub(c["saddle_energy_ev"].as_f64().expect("energy"),
                                         &c["unstable_direction"]);
                let mapping: Vec<[usize; 2]> =
                    serde_json::from_value(c["active_atom_mapping"].clone()).expect("mapping");
                match event_ids(c["origin_state_id"].as_str().expect("origin"),
                                c["destination_state_id"].as_str().expect("destination"), &saddle,
                                c["saddle_geometry_digest"].as_str().expect("saddle digest"), &mapping) {
                    Ok((pair, event, reverse)) =>
                        json!({"event_id": event, "input": c.clone(), "pair_id": pair,
                               "reverse_event_id": reverse}),
                    Err(e) => from_api(&e),
                }
            }).collect();
            json!({"cases": out})
        }
        "rate_pair" => {
            let out: Vec<Value> = input["cases"].as_array().expect("cases").iter().map(|c| {
                let g = |k: &str| c[k].as_f64().unwrap_or_else(|| panic!("rate field {k}"));
                match common_prefactor_pair(g("origin_ev"), g("destination_ev"), g("saddle_ev"),
                                            g("temperature"), g("prefactor"), g("barrier_tolerance"),
                                            g("detailed_balance_tolerance")) {
                    Ok(r) => json!({"barrier_ev": scalar(g("saddle_ev") - g("origin_ev")),
                                    "detailed_balance_residual": scalar(r.detailed_balance_residual),
                                    "log_forward_rate_per_s": scalar(r.log_forward_rate_per_s),
                                    "log_reverse_rate_per_s": scalar(r.log_reverse_rate_per_s),
                                    "outcome": "rated",
                                    "reverse_barrier_ev": scalar(g("saddle_ev") - g("destination_ev"))}),
                    Err(e) => from_status(&e),
                }
            }).collect();
            json!({"cases": out})
        }
        "rate_snapshot" => {
            match make_rate_snapshot(input["origin_state_id"].as_str().expect("origin"),
                                     &rows_of(input),
                                     input["log_rate_cutoff"].as_f64().expect("cutoff")) {
                Ok(s) => json!({"outcome": "snapshot", "snapshot": snapshot_value(&s)}),
                Err(e) => from_api(&e),
            }
        }
        "kmc_selection" => {
            let snapshot = match make_rate_snapshot(input["origin_state_id"].as_str().expect("origin"),
                                                    &rows_of(input),
                                                    input["log_rate_cutoff"].as_f64().expect("cutoff")) {
                Ok(s) => s, Err(e) => return from_api(&e),
            };
            let table = &snapshot.payload;
            let mut rng = derive_trajectory_stream(input["run_seed"].as_u64().expect("run_seed"));
            let mut steps = Vec::new();
            for index in 1..=input["steps"].as_u64().expect("steps") {
                let (us, ut, next) = match rng.two_uniforms_atomic() {
                    Ok(x) => x, Err(e) => { steps.push(from_status(&e)); break }
                };
                let threshold = us * table.total_rate_per_s;
                let mut sum = 0.0;
                let mut chosen = table.rates.len() - 1;
                for (j, r) in table.rates.iter().enumerate() {
                    sum += *r;
                    if sum > threshold { chosen = j; break }
                }
                let dt = -ut.ln() / table.total_rate_per_s;
                steps.push(json!({"post_state_id": table.destination_state_ids[chosen],
                                  "rng_after": serde_json::to_value(next.state()).expect("state"),
                                  "selected_event_id": table.event_ids[chosen],
                                  "selected_rate_per_s": scalar(table.rates[chosen]),
                                  "selection_uniform": scalar(us), "step_index": index,
                                  "time_increment_s": scalar(dt), "time_uniform": scalar(ut),
                                  "total_rate_per_s": scalar(table.total_rate_per_s)}));
                rng = next;
            }
            json!({"snapshot": snapshot_value(&snapshot), "steps": steps})
        }
        "status_records" => {
            let out: Vec<Value> = ALL_STATUS.iter().map(|s| json!({
                "exit_code_absorbing_requested": s.exit_code(true, false),
                "exit_code_default": s.exit_code(false, false),
                "exit_code_exploratory": s.exit_code(false, true),
                "message": stable_message(*s), "severity": s.severity(), "status": s})).collect();
            json!({"cases": out})
        }
        other => json!({"outcome": "unsupported_probe", "probe": other}),
    }
}

fn checkpoint_case(bytes: &[u8], model: &Value) -> Value {
    let vm = match validated(model) { Ok(v) => v, Err(e) => return from_api(&e) };
    match decode_checkpoint_v2(bytes, &vm) {
        Err(e) => from_api(&e),
        Ok(envelope) => {
            let p = &envelope.payload;
            let steps: Vec<Value> = p.trajectory.iter().map(|s| json!({
                "post_state_id": s.post_state_id, "pre_state_id": s.pre_state_id,
                "rate_table_snapshot": snapshot_value(&s.rate_table_snapshot),
                "selected_event_id": s.selected_event_id,
                "selected_rate_per_s": scalar(s.selected_rate_per_s),
                "selection_uniform": scalar(s.selection_uniform), "step_index": s.step_index,
                "time_increment_s": scalar(s.time_increment_s),
                "time_uniform": scalar(s.time_uniform),
                "total_rate_per_s": scalar(s.total_rate_per_s)})).collect();
            // `E2-CKPT-007`(1)/(2): the question this field answers is whether the backend
            // re-canonicalises and re-hashes the input to exactly the input. With two accepted
            // sealing rules (`E3-PAR-002`) the correct re-encode is under
            // the rule the input itself declared, which `decode_checkpoint_v2` recorded.
            let reencoded = encode_checkpoint_v2_sealed(
                p.clone(), &vm, CheckpointPolicy::default(), envelope.sealed_with);
            json!({"basin": serde_json::to_value(&p.basin).expect("basin"),
                   "catalog_digest": p.catalog.digest,
                   "catalog_event_ids": p.catalog.events.keys().collect::<Vec<_>>(),
                   "catalog_state_ids": p.catalog.states.keys().collect::<Vec<_>>(),
                   "checkpoint_sequence": p.checkpoint_sequence,
                   "current_state_id": p.current_state.state_id,
                   "flags": serde_json::to_value(&p.flags).expect("flags"),
                   "initial_state_id": p.initial_state.state_id,
                   "log_sequence": p.log_sequence, "outcome": "restored",
                   "reencoded_equals_input": reencoded.map(|b| b == bytes).unwrap_or(false),
                   "rng_trajectory": serde_json::to_value(&p.rng.trajectory).expect("rng"),
                   "simulation_time_s": scalar(p.simulation_time_s),
                   "step_index": p.step_index,
                   "substream_ids": p.rng.substream_map.keys().collect::<Vec<_>>(),
                   "trajectory": steps})
        }
    }
}

#[test]
fn emit_cross_language_parity_artifacts() {
    let dir = corpus();
    let out = out_dir();
    std::fs::create_dir_all(&out).expect("output directory");
    let manifest: Value = serde_json::from_slice(&std::fs::read(dir.join("manifest.json"))
        .expect("xlang manifest")).expect("manifest parses");
    let mut emitted: Vec<String> = Vec::new();
    let mut not_reread: Vec<Value> = Vec::new();
    let mut write = |id: &str, bytes: Vec<u8>| {
        // E2-CAN-001 / E2-CKPT-007(1): a canonical artifact must re-parse and re-canonicalise to
        // itself. Any case that does not is recorded, never swallowed; see D-E2-03 below.
        let parsed = parse_strict_json(&bytes)
            .unwrap_or_else(|e| panic!("{id} emitted JSON its own strict parser rejects: {e:?}"));
        let again = canonical_json_bytes(&parsed).expect("recanonicalises");
        if again != bytes {
            let offset = again.iter().zip(bytes.iter()).position(|(a, b)| a != b).unwrap_or(0);
            let window = |b: &[u8]| String::from_utf8_lossy(
                &b[offset.saturating_sub(40)..(offset + 40).min(b.len())]).into_owned();
            not_reread.push(json!({"emitted": window(&bytes), "first_diff_offset": offset,
                                   "id": id, "reparsed": window(&again)}));
        }
        std::fs::write(out.join(format!("{id}.out")), &bytes).expect("write case output");
        emitted.push(id.to_owned());
    };

    for entry in manifest["requests"].as_array().expect("requests array") {
        let id = entry["id"].as_str().expect("case id");
        assert!(matches!(entry["tier"].as_str(),Some("core"|"adapter")),
                "E3-PAR-001: {id} must declare core or adapter tier");
        let bytes = std::fs::read(dir.join(entry["file"].as_str().expect("file"))).expect("request bytes");
        assert_eq!(bytes.len() as u64, entry["bytes"].as_u64().expect("byte count"),
                   "{id}: corpus request length changed");
        // The crate now supplies a `RunAdapter` (`src/run.rs`), so the corpus is answered by the
        // backend as it actually is. `E2-SCOPE-004` puts the executable under the run-request
        // `extension`, so both `run` fixtures -- which carry no `calculator_command` -- are decided
        // by that rule rather than by the absence of an adapter.
        write(id, dispatch_json(&bytes, None, Some(&ProcessRunAdapter)));
    }

    let probes: Value = serde_json::from_slice(&std::fs::read(dir.join("probes.json"))
        .expect("probes")).expect("probes parse");
    let template: Value = serde_json::from_slice(include_bytes!("corpus/e2_minimal_model.json"))
        .expect("minimal model parses");
    for case in probes["cases"].as_array().expect("probe cases") {
        let id = case["id"].as_str().expect("probe id");
        assert_eq!(case["tier"],json!("core"),"E3-PAR-001: probes are core");
        write(id, canon(&probe(&case["input"], &template)));
    }

    let model: Value = serde_json::from_slice(&std::fs::read(dir.join("checkpoints/checkpoint_model.json"))
        .expect("checkpoint model")).expect("checkpoint model parses");
    for entry in manifest["checkpoints"].as_array().expect("checkpoint array") {
        let id = entry["id"].as_str().expect("case id");
        assert_eq!(entry["tier"],json!("core"),"E3-PAR-001: checkpoints are core");
        let bytes = std::fs::read(dir.join(entry["file"].as_str().expect("file"))).expect("checkpoint bytes");
        write(id, canon(&checkpoint_case(&bytes, &model)));
    }

    let expected = manifest["requests"].as_array().expect("requests").len()
        + probes["cases"].as_array().expect("probes").len()
        + manifest["checkpoints"].as_array().expect("checkpoints").len();
    assert_eq!(emitted.len(), expected, "every manifest case must emit exactly one artifact");
    let mut sorted = emitted.clone();
    sorted.sort();
    sorted.dedup();
    assert_eq!(sorted.len(), emitted.len(), "case IDs must be unique");
    let index: Map<String, Value> = [("backend".to_owned(), json!("rust")),
                                     ("cases".to_owned(), json!(sorted))].into_iter().collect();
    std::fs::write(out.join("_index.json"), canon(&Value::Object(index))).expect("write index");
    // D-E2-03 CLOSED (2026-08-11). `serde_json` was pinned without the `float_roundtrip` feature, so
    // its number parser was not exactly rounding and several canonical binary64 decimals came back one
    // ULP low: this crate could not re-read 4 of its own emitted artifacts, two of which carry
    // mandatory E2-PAR-002 item-4 corpus values. Enabling the feature closed it. E2-CKPT-007 step 1
    // and E2-PAR-003 both require canonical bytes to survive a parse and re-encode cycle, so the
    // requirement is now asserted directly: EVERY emitted artifact must re-read to identical bytes.
    // `_roundtrip.json` is still written so a regression names the exact cases rather than a count.
    let witness: Vec<&str> = not_reread.iter().map(|x| x["id"].as_str().expect("id")).collect();
    std::fs::write(out.join("_roundtrip.json"),
                   canon(&json!({"cases": not_reread, "defect": "D-E2-03", "status": "closed"})))
        .expect("write report");
    assert!(witness.is_empty(),
            "E2-CAN-004/E2-CKPT-007: every emitted artifact must re-read to identical bytes; \
             these did not: {witness:?}");
    eprintln!("xlang: wrote {} Rust artifacts to {}", emitted.len(), out.display());
}
