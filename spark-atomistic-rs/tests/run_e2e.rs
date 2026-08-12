// Clean-room execution tests for the `E2-API-003` run adapter, authored from
// OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2_PARITY. Independently authored; no
// implementation source consulted.
//
// The calculator is the shared clean-room fixture `../spark-atomistic/corpus/mock_calculator.py`,
// a deterministic stdlib-only JSON-lines process implementing
// `V = sum_i (x_i^2 - 1)^2 + y_i^2 + z_i^2` in eV and Angstrom. Its closed form is the acceptance
// oracle: minima at x = +-1 with E = 0, a first-order saddle at x = 0 with E = 1.0 eV, forward and
// reverse barriers exactly 1.000 eV, and Hessian eigenvalues (12x^2-4, 2, 2) at (x,0,0).
//
// Every rejection assertion below is paired with an accepted baseline that differs ONLY in the
// property under test.
use serde_json::{json, Map, Value};
use spark_atomistic_rs::checkpoint::{canonical_json_bytes, parse_strict_json};
use spark_atomistic_rs::model::Vec3;
use spark_atomistic_rs::parity::*;
use spark_atomistic_rs::resource::{CancelToken, ResourceLedger, ResourceLimits};
use spark_atomistic_rs::rng::derive_saddle_substream;
use spark_atomistic_rs::run::*;
use spark_atomistic_rs::status::StatusCode;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

const MODEL_DIGEST: &str = "7c799e3c0c25eb952d433430027d3d73de8d9f8f3d06064b6374f4b6eab4dd47";
/// The frozen marker atom is immovable and sits at a zero-force-free offset that keeps every state
/// certificate distinct: the movable atom at x = +1, x = 0 and x = -1 gives three different
/// anchor-relative geometries, so `E2-EVENT-001` "distinct committed IDs" is reachable.
const MARKER: Vec3 = [1.0, 0.5, 0.0];

fn calculator_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../spark-atomistic/corpus/mock_calculator.py")
}
fn calculator_argv() -> Value {
    let path = calculator_path();
    assert!(path.exists(), "shared calculator fixture missing: {}", path.display());
    json!(["python3", path.to_str().expect("UTF-8 fixture path")])
}
fn temp_dir(tag: &str) -> PathBuf {
    static SEQUENCE: AtomicU64 = AtomicU64::new(0);
    let dir = std::env::temp_dir().join(format!(
        "spark-run-{tag}-{}-{}", std::process::id(), SEQUENCE.fetch_add(1, Ordering::SeqCst)));
    std::fs::create_dir_all(&dir).expect("temporary directory");
    dir
}

fn model_json(dir: &PathBuf, start: Vec3) -> Value {
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
            "pbc": [false, false, false], "positions": [start, MARKER],
            "species": ["H", "He"], "spin": 0.0},
    })
}

fn validated(model: &Value) -> ValidatedModel {
    let wire: WireModel = serde_json::from_value(model.clone()).expect("model decodes");
    wire.validate(None).expect("model validates")
}
fn ledger_for(model: &ValidatedModel) -> ResourceLedger {
    let r = &model.model().resources;
    ResourceLedger::new(ResourceLimits {
        wall_time_s: r.wall_time_s, total_calculator_evaluations: r.total_calculator_evaluations,
        evaluations_per_relaxation: r.evaluations_per_relaxation,
        evaluations_per_saddle_attempt: r.evaluations_per_saddle_attempt,
        saddle_attempts_per_state: r.saddle_attempts_per_state, catalog_events: r.catalog_events,
        resident_memory_bytes: r.resident_memory_bytes, output_bytes: r.output_bytes,
        callback_timeout_s: r.callback_timeout_s, retry_count: 0, retry_backoff_s: 0.0,
    }, Default::default(), CancelToken::new()).expect("ledger")
}
fn extension() -> Map<String, Value> {
    let mut m = Map::new();
    m.insert("calculator_command".into(), calculator_argv());
    m
}

/// A fixed substream-derived unit mode, so the launch coordinate is the ONLY difference between
/// the accepted and the rejected search below (`DISC-002`, `E2-RNG-004`).
fn seeded_mode() -> Vec<Vec3> {
    let mut stream = derive_saddle_substream(0, "analytic-launch", "global", 0).expect("substream");
    let mut v = [0.0_f64; 3];
    for slot in v.iter_mut() { *slot = 2.0 * stream.next_uniform().expect("uniform") - 1.0; }
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    vec![[v[0] / n, v[1] / n, v[2] / n], [0.0, 0.0, 0.0]]
}

/// A deterministic JSON-lines calculator for `V = sum_i (x_i^2-1)^2 + y_i^2 + z_i^2 + 2c y_i z_i`,
/// written to a temporary directory at test time. It is the double-well fixture plus ONE extra
/// term, so the Hessian at the `x = 0` saddle is block `[[12x^2-4, 0, 0], [0, 2, 2c], [0, 2c, 2]]`:
/// the reaction coordinate is unchanged and the transverse block has eigenvalues `2 +- 2c` while
/// its Cartesian DIAGONAL stays `(2, 2)` for every `c`. That is exactly the configuration in which a
/// fixed set of sampled Rayleigh quotients certifies index 1 and is wrong.
fn coupled_calculator(dir: &PathBuf, coupling: f64) -> Value {
    let path = dir.join("coupled_calculator.py");
    let source = format!(r#"import hashlib, json, sys
C = {coupling:?}
request = json.loads(sys.stdin.buffer.readline().decode("utf-8"))
positions = request["state"]["positions"]
energy = 0.0
forces = []
for x, y, z in positions:
    energy += (x * x - 1.0) ** 2 + y * y + z * z + 2.0 * C * y * z
    forces.append([-4.0 * x * (x * x - 1.0), -(2.0 * y + 2.0 * C * z), -(2.0 * z + 2.0 * C * y)])
response = {{"status": "OK", "energy": energy, "forces": forces,
            "units": {{"energy": "eV", "forces": "eV/angstrom"}},
            "model_name": "analytic-coupled-well", "model_version": "1",
            "model_digest": "{MODEL_DIGEST}",
            "evaluation_id": "eval-" + hashlib.sha256(request["request_digest"].encode("ascii")).hexdigest(),
            "deterministic": True, "request_digest": request["request_digest"]}}
sys.stdout.write(json.dumps(response, sort_keys=True, separators=(",", ":")) + "\n")
"#);
    std::fs::write(&path, source).expect("calculator fixture written");
    json!(["python3", path.to_str().expect("UTF-8 path")])
}

/// Runs one `directional_dimer` from `(0.30, 0, 0)` with the reaction coordinate as the initial
/// mode, so the mode-following stage is not the variable under test and the only thing that changes
/// between arms is the coupling `c` and the curvature tolerance.
fn coupled_search(coupling: f64, curvature_tolerance: f64) -> Result<SaddlePoint, (StatusCode, String)> {
    let dir = temp_dir("coupled");
    let mut json_model = model_json(&dir, [0.30, 0.0, 0.0]);
    json_model["saddle_search"]["curvature_tolerance"] = json!(curvature_tolerance);
    json_model["saddle_search"]["max_iterations"] = json!(400);
    let model = validated(&json_model);
    let mut ext = Map::new();
    ext.insert("calculator_command".into(), coupled_calculator(&dir, coupling));
    let calc = CalculatorProcess::from_extension(&model, &ext).expect("adapter binding");
    let mut ledger = ledger_for(&model);
    let system = model.model().system.clone();
    let launch = vec![[0.30, 0.0, 0.0], MARKER];
    let mode = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    match directional_dimer(&calc, &system, &launch, &mode, &model.model().saddle_search,
        &mut ledger, "search:sha256:coupled").expect("no hard failure") {
        SearchOutcome::Found(p) => Ok(*p),
        SearchOutcome::Reject { status, termination_reason, .. } => Err((status, termination_reason)),
    }
}

/// `SADDLE-005` without a Hessian requires "negative curvature along the reported mode plus
/// nonnegative sampled orthogonal curvatures", and `E2-EVENT-005` makes `unstable_mode_count = 1`
/// mandatory. A FIXED sample cannot establish the second claim, because `v^T H v` on any chosen
/// basis is a positively weighted average of the spectrum in that direction. This crate's sampler
/// was deterministic (Gram-Schmidt Cartesian probes), not random, but it had the same defect. It is
/// now a MINIMISATION of the Rayleigh quotient over the orthogonal complement; see
/// `rayleigh_minimum_over_complement` in `src/run.rs`.
#[test]
fn saddle_order_evidence_is_minimised_over_the_orthogonal_complement_not_sampled() {
    // ACCEPTED BASELINE: c = 0. The transverse block is diag(2, 2); the fixed samples and the
    // minimum agree because the basis IS the eigenbasis.
    let baseline = coupled_search(0.0, 1e-6).expect("the uncoupled saddle is found");
    println!("c=0.0: curvature={:.9} orthogonal={:?}", baseline.curvature, baseline.orthogonal_curvatures);
    assert!((baseline.curvature + 4.0).abs() < 1e-3, "reaction-coordinate curvature is -4");
    assert!(baseline.orthogonal_curvatures.iter().all(|x| (*x - 2.0).abs() < 1e-2),
        "every recorded quotient is the true transverse eigenvalue 2: {:?}", baseline.orthogonal_curvatures);

    // WITNESS: c = 1.5 makes the transverse block [[2, 3], [3, 2]] -- eigenvalues 5 and -1 -- while
    // its Cartesian diagonal is unchanged at (2, 2). `curvature_tolerance = 2.0` keeps the record
    // acceptable so its contents can be read: the seeds are the old fixed sample and are BOTH
    // nonnegative, and the minimised quotient is the true -1 that the sample never saw.
    let coupled = coupled_search(1.5, 2.0).expect("the coupled saddle is found at a loose tolerance");
    println!("c=1.5: curvature={:.9} orthogonal={:?}", coupled.curvature, coupled.orthogonal_curvatures);
    assert_eq!(coupled.orthogonal_curvatures.len(), 3, "two seeds plus the minimised quotient");
    assert!(coupled.orthogonal_curvatures[0] >= 0.0 && coupled.orthogonal_curvatures[1] >= 0.0,
        "the fixed Cartesian samples are nonnegative and would have certified index 1: {:?}",
        coupled.orthogonal_curvatures);
    assert!((coupled.orthogonal_curvatures[2] + 1.0).abs() < 5e-2,
        "the minimisation reaches the true smallest transverse eigenvalue -1: {:?}",
        coupled.orthogonal_curvatures);

    // REJECTION, paired with the c = 0 baseline: identical launch, identical initial mode,
    // identical tolerances; ONLY the coupling differs, and the minimisation refutes index 1.
    let (status, reason) = coupled_search(1.5, 1e-6).expect_err("the coupled saddle is not first order");
    println!("c=1.5 at the default tolerance: status={status:?} termination_reason={reason}");
    assert_eq!(status, StatusCode::InvalidSaddle);
    assert_eq!(reason, "additional-negative-curvature");
}

struct Analytic { saddle_energy: f64, forward: f64, reverse: f64, evaluations: u64, iterations: u64 }

fn analytic_from(x0: f64, force_tolerance: f64, max_iterations: u64) -> Result<Analytic, (StatusCode, String, u64, u64)> {
    let dir = temp_dir("analytic");
    let mut json_model = model_json(&dir, [x0, 0.0, 0.0]);
    json_model["relaxation"]["force_tolerance"] = json!(force_tolerance);
    json_model["saddle_search"]["force_tolerance"] = json!(force_tolerance);
    json_model["saddle_search"]["max_iterations"] = json!(max_iterations);
    let model = validated(&json_model);
    let calc = CalculatorProcess::from_extension(&model, &extension()).expect("adapter binding");
    let mut ledger = ledger_for(&model);
    let system = model.model().system.clone();
    let launch = vec![[x0, 0.0, 0.0], MARKER];
    let outcome = directional_dimer(&calc, &system, &launch, &seeded_mode(),
        &model.model().saddle_search, &mut ledger, "search:sha256:analytic").expect("no hard failure");
    let point = match outcome {
        SearchOutcome::Found(p) => *p,
        SearchOutcome::Reject { status, termination_reason, evaluations, iterations } =>
            return Err((status, termination_reason, evaluations, iterations)),
    };
    let mut endpoints = Vec::new();
    for sign in [-1.0_f64, 1.0] {
        let mut positions = point.positions.clone();
        for k in 0..3 { positions[0][k] += sign * model.model().saddle_search.endpoint_displacement * point.direction[0][k]; }
        let mut moved = system.clone();
        moved.positions = positions;
        let relaxed = relax(&calc, &moved, force_tolerance, model.model().relaxation.max_steps,
            model.model().relaxation.max_evaluations, &mut ledger, &format!("relax:{sign}")).expect("endpoint relaxes");
        endpoints.push(relaxed);
    }
    Ok(Analytic {
        saddle_energy: point.energy_ev,
        forward: point.energy_ev - endpoints[0].energy_ev,
        reverse: point.energy_ev - endpoints[1].energy_ev,
        evaluations: point.evaluations,
        iterations: point.iterations,
    })
}

#[test]
fn analytic_double_well_barriers_reproduce_the_closed_form_from_x0_030() {
    // ACCEPTED BASELINE. At (0.30, 0, 0) the Hessian is diag(12*0.09-4, 2, 2) = (-2.92, 2, 2):
    // the reaction coordinate is the lowest mode and it is negative, which is the documented
    // domain of a minimum-mode follower (`SADDLE-002`, `SADDLE-005`).
    let a = analytic_from(0.30, 1e-4, 400).expect("search converges inside its domain");
    // The marker atom contributes a constant (1-1)^2 + 0.5^2 + 0 = 0.25 eV to every energy, so the
    // absolute saddle energy is 1.25 eV and both BARRIERS are the closed-form 1.000 eV.
    println!("x0=0.30 ft=1e-4: E_saddle={:.12} forward={:.12} reverse={:.12} iterations={} evaluations={}",
        a.saddle_energy, a.forward, a.reverse, a.iterations, a.evaluations);
    assert!((a.saddle_energy - 1.25).abs() < 1e-8, "absolute saddle energy {}", a.saddle_energy);
    assert!((a.forward - 1.0).abs() < 5e-9, "forward barrier {} deviates by {:e}", a.forward, (a.forward - 1.0).abs());
    assert!((a.reverse - 1.0).abs() < 5e-9, "reverse barrier {} deviates by {:e}", a.reverse, (a.reverse - 1.0).abs());
}

#[test]
fn analytic_saddle_accuracy_is_bounded_by_the_force_tolerance_not_by_the_algorithm() {
    // The saddle-energy error is set by where inside the force-tolerance ball the search stops,
    // not by the algorithm. Along the reaction coordinate V ~ E_s - 2x^2 with |F_x| ~ 4|x|, which
    // alone gives f^2/8; the transverse coordinates contribute V ~ y^2 with |F_y| = 2|y|, which
    // alone gives f^2/4. The residual force can sit anywhere on that sphere, so f^2/4 is the
    // attainable bound for this potential and f^2/8 is the reaction-coordinate-only reference.
    for tolerance in [5e-2_f64, 1e-2, 1e-3, 1e-4, 1e-5] {
        let a = analytic_from(0.30, tolerance, 400).expect("search converges inside its domain");
        let error = (a.saddle_energy - 1.25).abs();
        println!("force_tolerance={:e} |E_saddle-1.25|={:e} f^2/8={:e} f^2/4={:e} forward={:.12} reverse={:.12} evaluations={}",
            tolerance, error, tolerance * tolerance / 8.0, tolerance * tolerance / 4.0, a.forward, a.reverse, a.evaluations);
        assert!(error <= tolerance * tolerance / 4.0 + 1e-15, "error {error:e} exceeds the f^2/4 bound");
    }
}

#[test]
fn analytic_double_well_search_from_x0_090_leaves_the_solver_domain() {
    // REJECTION, paired with the x0 = 0.30 baseline above: identical model, identical seeded
    // initial mode, identical tolerances; ONLY the launch coordinate differs. At (0.90, 0, 0) the
    // Hessian is diag(5.72, 2, 2) -- every eigenvalue positive and the lowest mode transverse --
    // so minimum-mode following climbs the wrong mode. `E2-STATUS-002` classifies the answer as a
    // `candidate reject` that cannot terminate a public operation.
    match analytic_from(0.90, 1e-4, 400) {
        Ok(a) => panic!("unexpected convergence from x0=0.90: E_saddle={}", a.saddle_energy),
        Err((status, reason, evaluations, iterations)) => {
            println!("x0=0.90 ft=1e-4: status={status:?} termination_reason={reason} iterations={iterations} evaluations={evaluations}");
            assert_eq!(status, StatusCode::SaddleNotFound);
        }
    }
}

#[test]
fn run_request_without_a_calculator_command_is_invalid_input() {
    // `E2-SCOPE-004`: the executable belongs only to the run-request `extension`, so a request
    // that carries none has no adapter binding at all.
    let dir = temp_dir("no-command");
    let model = model_json(&dir, [0.95, 0.02, -0.02]);
    let request = json!({"allow_unvalidated": true, "extension": {}, "model": model, "operation": "run"});
    let bytes = serde_json::to_vec(&request).expect("request encodes");
    let response = parse_strict_json(&dispatch_json(&bytes, None, Some(&ProcessRunAdapter))).expect("canonical response");
    assert_eq!(response["status"], json!("INVALID_INPUT"));
    assert_eq!(response["exit_code"], json!(64));
    assert_eq!(response["context"]["requirement_id"], json!("E2-SCOPE-004"));
    // ACCEPTED BASELINE: the same request with the executable present gets past the binding and
    // reaches scientific work, which is what proves the rejection above is the binding rule and
    // not an unrelated malformation.
    let dir = temp_dir("with-command");
    let mut model = model_json(&dir, [0.95, 0.02, -0.02]);
    model["kinetics"]["maximum_steps"] = json!(1);
    model["discovery"]["maximum_attempts"] = json!(1);
    model["resources"]["saddle_attempts_per_state"] = json!(1);
    model["saddle_search"]["max_iterations"] = json!(5);
    let request = json!({"allow_unvalidated": true, "extension": extension(), "model": model, "operation": "run"});
    let bytes = serde_json::to_vec(&request).expect("request encodes");
    let response = parse_strict_json(&dispatch_json(&bytes, None, Some(&ProcessRunAdapter))).expect("canonical response");
    assert_ne!(response["context"]["requirement_id"], json!("E2-SCOPE-004"),
        "the accepted baseline must not stop at the adapter-binding rule: {response}");
    assert!(dir.join("checkpoint.json").exists(), "the initial minimum was committed and checkpointed");
}

#[test]
fn calculator_stderr_is_drained_instead_of_deadlocking_the_child() {
    let dir = temp_dir("stderr-drain");
    let mut model = model_json(&dir, [0.30, 0.0, 0.0]);
    model["resources"]["callback_timeout_s"] = json!(1.0);
    let model = validated(&model);
    let extension = json!({"calculator_command": ["/usr/bin/python3", "-c",
        "import sys; sys.stdin.buffer.readline(); sys.stderr.write('x'*1000000); sys.stderr.flush(); sys.stdout.write('{}\\n')"]})
        .as_object().expect("object").clone();
    let calculator = CalculatorProcess::from_extension(&model, &extension).expect("adapter");
    let error = calculator.evaluate(&model.model().system).expect_err("empty response object is invalid");
    assert_eq!(error.status, StatusCode::CalculatorFailure);
    assert_ne!(error.details.get("transport"), Some(&json!("callback timeout")));
}

#[test]
fn callback_timeout_covers_process_exit_after_stdout_closes() {
    let dir = temp_dir("process-timeout");
    let mut model = model_json(&dir, [0.30, 0.0, 0.0]);
    model["resources"]["callback_timeout_s"] = json!(0.05);
    let model = validated(&model);
    let extension = json!({"calculator_command": ["/usr/bin/python3", "-c",
        "import sys,time; sys.stdin.buffer.readline(); sys.stdout.close(); time.sleep(1.0)"]})
        .as_object().expect("object").clone();
    let calculator = CalculatorProcess::from_extension(&model, &extension).expect("adapter");
    let error = calculator.evaluate(&model.model().system).expect_err("child must time out");
    assert_eq!(error.status, StatusCode::CalculatorFailure);
    assert_eq!(error.details.get("transport"), Some(&json!("callback timeout")));
}

#[test]
fn run_operation_executes_the_whole_pipeline_and_reports_where_it_stops() {
    let dir = temp_dir("e2e");
    let json_model = model_json(&dir, [0.95, 0.02, -0.02]);
    let model = validated(&json_model);
    let (outcome, report) = execute_run(&model, &extension());
    println!("outcome: {:?}", outcome.as_ref().map_err(|e| (e.status, e.requirement_id.clone())));
    println!("calculator evaluations: {} checkpoints written: {}", report.calculator_evaluations, report.checkpoints_written);
    println!("committed states: {} directed events: {}", report.states.len(), report.events.len());
    for (id, state) in &report.states {
        println!("  state {id} E={:.12} x={:?}", state.energy_ev, state.positions[0]);
    }
    for (id, e) in &report.events {
        println!("  event {id}\n    {} -> {} barrier={:.12} reverse={:.12} E_saddle={:.12} x_saddle={:?} log_k={:.9} selectable={} pair={}",
            e.origin_state_id, e.destination_state_id, e.barrier_ev, e.reverse_barrier_ev,
            e.saddle.energy_ev, e.saddle.positions[0], e.rate_model.log_forward_rate_per_s, e.selectable, e.pair_id);
    }
    for (id, d) in &report.discovery {
        println!("  discovery {id}: attempts={} successes={} duplicates={} consecutive={} failures={} stopping={}",
            d.attempts, d.successes, d.duplicates, d.consecutive_redundant_successes,
            serde_json::to_string(&d.failures_by_status).unwrap_or_default(), d.stopping_state);
    }
    for record in &report.trajectory {
        println!("  step {} {} -> {} rate={} dt={:e}", record.step_index, record.pre_state_id,
            record.post_state_id, record.selected_rate_per_s, record.time_increment_s);
    }
    for name in ["trajectory.json", "summary.json", "checkpoint.json"] {
        let path = dir.join(name);
        println!("artifact {name}: exists={} bytes={}", path.exists(),
            std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0));
    }
    // The physical catalog of this fixture is exactly one reciprocal pair: the movable atom
    // crossing x = 0 between the two wells.
    assert_eq!(report.states.len(), 2, "two distinct committed minima");
    assert_eq!(report.events.len(), 2, "one reciprocal pair = two directed records");
    for e in report.events.values() {
        assert!((e.barrier_ev - 1.0).abs() < 1e-6, "forward barrier {} is the closed-form 1.0 eV", e.barrier_ev);
        assert!((e.reverse_barrier_ev - 1.0).abs() < 1e-6, "reverse barrier {} is the closed-form 1.0 eV", e.reverse_barrier_ev);
        assert!(e.selectable && e.validation.full_endpoint_relaxations && e.validation.unstable_mode_count == 1);
    }
    let stats = report.discovery.values().next().expect("discovery statistics for the origin state");
    assert_eq!(stats.stopping_state, "CONVERGED_HEURISTIC", "DISC-004 stopping rule fired");
    assert!(stats.duplicates >= 1 && stats.consecutive_redundant_successes >= 1,
        "the redundant-success criterion is what stopped discovery (`DISC-005`)");
    assert_eq!(report.trajectory.len(), 2, "both KMC steps of `kinetics.maximum_steps` committed");

    // `D-E2-02` resolved by normative `E3-EVENT-001` (adopted by `D-127`).
    // Before it, `validate_checkpoint_v2` demanded simultaneously that
    //   C1 saddle.search_id == discovery_provenance.search_id,
    //   C2 the reciprocal record's `saddle` equal the forward record's `saddle`, and
    //   C3 each record's provenance re-derive from THAT record's own `origin_state_id`,
    // which is unsatisfiable once `E2-EVENT-001` forces the endpoints to be distinct, so the run
    // committed real science and was then refused its own checkpoint (`CHECKPOINT_CORRUPT` /
    // `E2-CKPT-007`, exit 74, exactly one accepted checkpoint -- the empty-catalog one).
    // Options 1 + 4 drop C3 and relax C2 to mapped equality, and `run` now completes.
    // PAIRED CONTROL below: the SAME committed payload under `CheckpointPolicy::PreErrata3` is
    // still refused, so the difference is the invariant set and not the science.
    let value = outcome.expect("run completes end to end under the E3-EVENT-001 candidate");
    assert_eq!(value.step_index, 2);
    assert_eq!(value.checkpoint_sequence, 4, "bootstrap + one per step + the final completion write");
    assert!(!value.incomplete_catalog);
    assert!(value.simulation_time_s.is_finite() && value.simulation_time_s > 0.0);
    for name in ["trajectory.json", "summary.json", "checkpoint.json"] {
        assert!(dir.join(name).exists(), "{name} was written");
    }
    let bytes = std::fs::read(dir.join("checkpoint.json")).expect("checkpoint readable");
    let envelope = decode_checkpoint_v2_with(&bytes, &model, CheckpointPolicy::Errata3Event001)
        .expect("the crate re-reads the catalog-bearing checkpoint it just wrote");
    assert_eq!(envelope.payload.catalog.events.len(), 2, "one reciprocal pair survives the round trip");
    assert_eq!(envelope.sealed_with, PayloadHashPolicy::Errata3EnvironmentExcluded);
    let (status, requirement) = decode_checkpoint_v2_with(&bytes, &model, CheckpointPolicy::PreErrata3)
        .map_err(|e| (e.status, e.requirement_id)).expect_err("the shipped invariant set still refuses it");
    // The shipped set refuses at its FIRST unstated constraint, `validation.method ==
    // saddle.evidence_level`, which is inside `validate_event_v2` and therefore cited to
    // `E2-EVENT-001`; C2 and C3 are never reached. That is the same first rung the proposal
    // measured for the Python bytes ("Python bytes as written -> CHECKPOINT_CORRUPT /
    // `E2-EVENT-001`"), so this backend's own records now fail the shipped rule in the same way the
    // other backend's did.
    assert_eq!(status, StatusCode::CheckpointCorrupt);
    assert_eq!(requirement, "E2-EVENT-001");
}

/// `E2-CKPT-008`: "Replayed IDs, uniforms, rate arrays, selections, and counters are byte-identical."
/// This drives it FORWARD: a run that resumes from a catalog-bearing checkpoint must commit the same
/// NEXT event as an uninterrupted run of the same science. Before the `E3-EVENT-001` candidate no
/// such checkpoint could exist, so this path was unreachable.
///
/// One substitution is unavoidable and is stated rather than hidden: `E2-CAN-007` hashes the
/// COMPLETE validated model into `config_digest`, and `output.resume` is part of that model, so a
/// checkpoint written by a non-resuming run can never satisfy `E2-CKPT-007`(3) for the resuming run.
/// The priming checkpoint therefore has `digests.config` and `discovery_statistics[*].config_digest`
/// -- the two places that binding appears -- rewritten to the resuming model's digest and is
/// resealed. Nothing else is touched; every other field is what the priming run actually committed.
#[test]
fn resume_from_a_catalog_bearing_checkpoint_commits_the_same_next_event() {
    // REFERENCE: one uninterrupted run of two steps.
    let reference_dir = temp_dir("resume-reference");
    let reference_model = validated(&model_json(&reference_dir, [0.95, 0.02, -0.02]));
    let (outcome, reference) = execute_run(&reference_model, &extension());
    outcome.expect("the uninterrupted reference run completes");
    assert_eq!(reference.trajectory.len(), 2, "the reference committed two steps");

    // PRIMING: the same science stopped after one step, which is the only model difference.
    let priming_dir = temp_dir("resume-priming");
    let mut priming_json = model_json(&priming_dir, [0.95, 0.02, -0.02]);
    priming_json["kinetics"]["maximum_steps"] = json!(1);
    let (outcome, priming) = execute_run(&validated(&priming_json), &extension());
    outcome.expect("the priming run completes its single step");
    assert_eq!(priming.trajectory.len(), 1);
    assert_eq!(normalise_sequence(&priming.trajectory[0]), normalise_sequence(&reference.trajectory[0]),
        "the priming run and the reference agree on step 1 before any resume happens");

    // RESUME: the reference model plus `output.resume`, fed the priming checkpoint.
    let resume_dir = temp_dir("resume-run");
    let mut resume_json = model_json(&resume_dir, [0.95, 0.02, -0.02]);
    resume_json["output"]["resume"] = json!(true);
    resume_json["output"]["overwrite"] = json!(true);
    let resume_model = validated(&resume_json);
    let primed: Value = parse_strict_json(&std::fs::read(priming_dir.join("checkpoint.json"))
        .expect("priming checkpoint")).expect("checkpoint parses");
    let mut payload = primed["payload"].clone();
    payload["digests"]["config"] = json!(resume_model.config_digest());
    for (_, stats) in payload["discovery_statistics"].as_object_mut().expect("discovery statistics") {
        stats["config_digest"] = json!(resume_model.config_digest());
    }
    let typed: CheckpointPayloadV2 = serde_json::from_value(payload).expect("payload decodes");
    let envelope = CheckpointEnvelopeV2::new_sealed(typed, PayloadHashPolicy::Errata3EnvironmentExcluded)
        .expect("resealed");
    std::fs::write(resume_dir.join("checkpoint.json"),
        canonical_json_bytes(&envelope).expect("canonical envelope")).expect("primed checkpoint written");

    let (outcome, resumed) = execute_run(&resume_model, &extension());
    outcome.expect("the resumed run completes");
    assert_eq!(resumed.trajectory.len(), 2, "one restored step plus one newly committed step");
    // `checkpoint_sequence` counts how many checkpoints THIS process wrote before the step, so it
    // is bookkeeping and not science; every other field of the committed step must be identical.
    assert_eq!(normalise_sequence(&resumed.trajectory[1]), normalise_sequence(&reference.trajectory[1]),
        "resume reproduces the next committed event exactly");
    println!("next event after resume: {} -> {} (uniform {:?}, dt {:?})",
        resumed.trajectory[1].selected_event_id, resumed.trajectory[1].post_state_id,
        resumed.trajectory[1].selection_uniform, resumed.trajectory[1].time_increment_s);
}

fn normalise_sequence(step: &KmcStepV2) -> KmcStepV2 {
    let mut out = step.clone();
    out.checkpoint_sequence = 0;
    out
}

/// `run` end to end through the PUBLIC wire on a model that discovers an event, plus the
/// `E3-PAR-002` reproducibility claim.
///
/// `E2-CKPT-005` mandates `resources.wall_elapsed_s` and `E2-PAR-003` demands byte identity, which
/// no implementation can satisfy against ITSELF. Excluding the declared environment-derived set
/// from the hashed payload makes `payload_sha256` reproducible while the field stays in the record.
/// ACCEPTED BASELINE and PAIRED MEASUREMENT are the same two runs: identical model, identical
/// directory, identical everything except the wall clock.
#[test]
fn run_returns_ok_on_the_wire_and_two_identical_runs_agree_on_the_payload_hash() {
    let dir = temp_dir("reproducible");
    let mut model = model_json(&dir, [0.95, 0.02, -0.02]);
    model["output"]["overwrite"] = json!(true);
    let request = json!({"allow_unvalidated": true, "extension": extension(),
        "model": model, "operation": "run"});
    let bytes = serde_json::to_vec(&request).expect("request encodes");
    let mut sealed = Vec::new();
    let mut elapsed = Vec::new();
    let mut payloads = Vec::new();
    for pass in 0..2 {
        let response = parse_strict_json(&dispatch_json(&bytes, None, Some(&ProcessRunAdapter)))
            .expect("canonical response");
        assert_eq!(response["status"], json!("OK"), "pass {pass}: {response}");
        assert_eq!(response["exit_code"], json!(0));
        assert_eq!(response["message"], json!("transaction committed"));
        assert_eq!(response["value"]["step_index"], json!(2));
        let checkpoint = parse_strict_json(&std::fs::read(dir.join("checkpoint.json"))
            .expect("checkpoint readable")).expect("checkpoint parses");
        sealed.push(checkpoint["payload_sha256"].as_str().expect("hash").to_owned());
        elapsed.push(checkpoint["payload"]["resources"]["wall_elapsed_s"].as_f64().expect("elapsed"));
        let mut payload = checkpoint["payload"].clone();
        payload["resources"].as_object_mut().expect("resources").remove("wall_elapsed_s");
        payloads.push(payload);
    }
    println!("payload_sha256: {} / {}", sealed[0], sealed[1]);
    println!("wall_elapsed_s: {} / {}", elapsed[0], elapsed[1]);
    // Assert the CONTENT before the hash: if the two hashed payloads differ, this names the
    // offending field, whereas a hash comparison alone only reports two opaque digests. A
    // failure here means the declared environment-derived set is incomplete; a failure of the
    // hash assertion below with this one passing would instead mean the sealing is wrong.
    if payloads[0] != payloads[1] {
        let (a, b) = (payloads[0].as_object().expect("obj"), payloads[1].as_object().expect("obj"));
        for (k, va) in a {
            if b.get(k) != Some(va) {
                println!("DIFFERING TOP-LEVEL KEY: {k}");
                if let (Some(oa), Some(ob)) = (va.as_object(), b.get(k).and_then(|x| x.as_object())) {
                    for (k2, v2) in oa {
                        if ob.get(k2) != Some(v2) {
                            println!("  {k}.{k2}:\n    A = {v2}\n    B = {}", ob.get(k2).map(|x| x.to_string()).unwrap_or_else(|| "<absent>".into()));
                        }
                    }
                }
            }
        }
    }
    assert_eq!(payloads[0], payloads[1], "everything the hash covers is identical");
    assert_eq!(sealed[0], sealed[1], "two identical runs must now agree on `payload_sha256`");
    assert_ne!(elapsed[0], elapsed[1], "`E2-CKPT-005`'s field is still recorded and still differs");
}
