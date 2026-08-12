// Cross-language RUN parity: the Rust half. Authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 +
// ERRATA_2_PARITY. Independently authored; no implementation source consulted.
//
// `E2-PAR-001` requires both implementations to consume the SAME canonical corpus and
// `E2-PAR-003` requires the emitted checkpoint bytes to be identical. `tests/corpus/xlang/` covers
// every static artifact; nothing covered a `run`, because this crate had no adapter. This test
// emits the run half: it writes ONE canonical model (so the config digest is shared verbatim),
// executes it, and files the three artifacts under `<dir>/rust/` for byte comparison against the
// Python backend driven from the same `<dir>/model.json`.
//
// It is a real test and never a silent skip: with no `SPARK_RUN_PARITY_DIR` it still runs the model
// in a temporary directory and asserts the same invariants.
use serde_json::{json, Map, Value};
use spark_atomistic_rs::checkpoint::{canonical_json_bytes, parse_strict_json};
use spark_atomistic_rs::parity::*;
use spark_atomistic_rs::run::ProcessRunAdapter;
use std::path::{Path, PathBuf};

const MODEL_DIGEST: &str = "7c799e3c0c25eb952d433430027d3d73de8d9f8f3d06064b6374f4b6eab4dd47";

fn shared_model(dir: &Path) -> Value {
    json!({
        "basin": {"enabled": false},
        "calculator": {"deterministic": true, "model_digest": MODEL_DIGEST,
            "model_name": "analytic-double-well", "model_version": "1"},
        "discovery": {"alpha": null, "alpha_calibration": null,
            "classes": [{"kind": "global", "name": "global", "probability": 1.0}],
            "consecutive_redundant": 1, "maximum_attempts": 12, "maximum_evaluations": 100000,
            "minimum_successful": 1, "mode": "strict", "relevance_rate_min": 0.0},
        "kinetics": {"absorbing_ok": false, "barrier_tolerance": 1e-10,
            "detailed_balance_tolerance": 1e-08, "log_rate_cutoff": -700.0, "maximum_steps": 2,
            "prefactor": 1e13, "rate_model": "COMMON_PREFACTOR", "run_seed": 0,
            "saddle_energy_tolerance": 1e-05, "saddle_max_tolerance": 0.005,
            "saddle_rms_tolerance": 0.001, "state_energy_tolerance_per_atom": 1e-06,
            "state_max_tolerance": 0.005, "state_rms_tolerance": 0.001, "temperature": 300.0},
        "output": {"checkpoint_every_steps": 1,
            "checkpoint_path": dir.join("checkpoint.json").to_str().expect("path"),
            "checkpoint_wall_time_s": 600.0, "overwrite": false, "resume": false,
            "summary_path": dir.join("summary.json").to_str().expect("path"),
            "trajectory_path": dir.join("trajectory.json").to_str().expect("path")},
        "relaxation": {"force_tolerance": 0.0001, "max_evaluations": 2000, "max_steps": 2000},
        "resources": {"callback_timeout_s": 60.0, "catalog_events": 100,
            "evaluations_per_relaxation": 2000, "evaluations_per_saddle_attempt": 8000,
            "output_bytes": 100000000, "resident_memory_bytes": 8589934592_u64,
            "retry_backoff_s": 0.0, "retry_count": 0, "saddle_attempts_per_state": 12,
            "total_calculator_evaluations": 200000, "wall_time_s": 1200.0},
        "saddle_search": {"curvature_tolerance": 1e-06, "endpoint_displacement": 1.0,
            "force_tolerance": 0.0001, "max_iterations": 150, "method": "directional-dimer",
            "orthogonal_directions": 2},
        "schema": {"id": "spark-atomistic-model/1"},
        "system": {"atom_ids": ["a0", "a1"], "calculator_model_digest": MODEL_DIGEST,
            "cell": [[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]], "charge": 0.0,
            "constraints": {"kind": "fixed-mask"}, "movable": [true, false],
            "pbc": [false, false, false], "positions": [[0.95, 0.02, -0.02], [1.0, 0.5, 0.0]],
            "species": ["H", "He"], "spin": 0.0},
    })
}

#[test]
fn emit_cross_language_run_artifacts() {
    let dir = match std::env::var_os("SPARK_RUN_PARITY_DIR") {
        Some(d) => PathBuf::from(d),
        None => std::env::temp_dir().join(format!("spark-run-parity-{}", std::process::id())),
    };
    std::fs::create_dir_all(&dir).expect("parity directory");
    let model_path = dir.join("model.json");
    // ONE canonical model, written once and read verbatim by the other backend (`E2-PAR-001`).
    let model_value = shared_model(&dir);
    if !model_path.exists() {
        std::fs::write(&model_path, canonical_json_bytes(&model_value).expect("canonical model"))
            .expect("model written");
    }
    let stored: Value = serde_json::from_slice(&std::fs::read(&model_path).expect("model readable"))
        .expect("model parses");
    assert_eq!(canonical_json_bytes(&stored).expect("canonical"),
        canonical_json_bytes(&model_value).expect("canonical"),
        "the stored model must be the model this test drives");

    let calculator = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../spark-atomistic/corpus/mock_calculator.py");
    assert!(calculator.exists(), "shared calculator fixture missing: {}", calculator.display());
    let mut extension = Map::new();
    extension.insert("calculator_command".into(),
        json!(["python3", calculator.to_str().expect("UTF-8 fixture path")]));

    // ONE execution, through the public wire, so the response envelope and the three artifacts are
    // all products of the same run.
    let request = json!({"allow_unvalidated": true, "extension": extension,
        "model": stored, "operation": "run"});
    let bytes = serde_json::to_vec(&request).expect("request encodes");
    let wire_response = dispatch_json(&bytes, None, Some(&ProcessRunAdapter));

    let out = dir.join("rust");
    std::fs::create_dir_all(&out).expect("rust artifact directory");
    for name in ["checkpoint.json", "trajectory.json", "summary.json"] {
        let src = dir.join(name);
        if src.exists() { std::fs::rename(&src, out.join(name)).expect("artifact filed"); }
    }
    std::fs::write(out.join("_response.json"), &wire_response).expect("response written");
    println!("rust run response: {}", String::from_utf8_lossy(&wire_response));
    for name in ["checkpoint.json", "trajectory.json", "summary.json", "_response.json"] {
        let path = out.join(name);
        println!("  rust/{name}: exists={} bytes={}", path.exists(),
            std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0));
    }
    assert!(out.join("checkpoint.json").exists(), "a checkpoint was emitted for comparison");
    assert!(out.join("trajectory.json").exists(), "a trajectory was emitted for comparison");
    assert!(out.join("summary.json").exists(), "a summary was emitted for comparison");
    // `E2-CAN-001`/`E2-JSON-002`: every emitted artifact must re-parse under the strict reader and
    // re-canonicalise to exactly its own bytes. This is the `D-E2-03` round-trip check applied to
    // the run artifacts, which the static corpus could not reach.
    for name in ["checkpoint.json", "trajectory.json", "summary.json", "_response.json"] {
        let bytes = std::fs::read(out.join(name)).expect("artifact readable");
        let again = canonical_json_bytes(&parse_strict_json(&bytes).expect("artifact re-parses"))
            .expect("artifact re-canonicalises");
        assert_eq!(again, bytes, "{name} does not re-read to its own canonical bytes");
    }
}
