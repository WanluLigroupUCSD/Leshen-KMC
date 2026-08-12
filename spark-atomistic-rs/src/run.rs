// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2_PARITY.
// Independently authored; no implementation source consulted.
//
// This module supplies the `parity::RunAdapter` the crate previously lacked, so that an
// `E2-API-003` run request can execute end to end: process-isolated JSON-lines calculator,
// minimiser, first-order-saddle solver, endpoint validation, directed reversible catalog,
// `COMMON_PREFACTOR` rates, deterministic serial KMC, and canonical checkpoint/resume.
//
// `E2-SCOPE-004` is the hard boundary of this file: the transport, executable, IPC, environment
// and timeout-enforcement mechanism are read from the run-request `extension` and are used for
// NOTHING except spawning the child process. No value taken from `extension` is ever written into
// the model, a capability value, a digest, an ID, a checkpoint payload, a trajectory record, or a
// scientific result. Every provenance identity recorded below is derived from the MODEL
// (`calculator.model_name`/`model_version`) or is a fixed constant of this adapter.
//
// Every numeric policy constant of the solvers is a declared constant of this file, not a hidden
// configuration channel: `E2-SCHEMA-006` fixes the tolerances the erratum pins and the erratum
// pins no perturbation amplitude, dimer separation, or step size, so those are stated here.
use crate::checkpoint::{canonical_json_bytes, parse_strict_json};
use crate::identity::{closest_periodic_vector, hex_sha256};
use crate::model::{Cell, Vec3};
use crate::parity::*;
use crate::rate::common_prefactor_pair;
use crate::resource::{CancelToken, ResourceLedger, ResourceLimits};
use crate::rng::{derive_saddle_substream, derive_trajectory_stream, substream_digest, Philox};
use crate::status::StatusCode;
use serde_json::{json, Map, Value};
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

/// Finite-difference dimer half-separation, A. Not pinned by any requirement.
pub const DIMER_SEPARATION_ANGSTROM: f64 = 1.0e-3;
/// Minimum-mode-following translation step, A per eV A^-1. Not pinned by any requirement.
pub const DIMER_TRANSLATION_STEP: f64 = 0.1;
/// Cap on one dimer translation, A. Not pinned by any requirement.
pub const DIMER_MAX_TRANSLATION_ANGSTROM: f64 = 0.2;
/// Dimer rotations evaluated per translation. Not pinned by any requirement.
pub const DIMER_ROTATIONS_PER_ITERATION: u32 = 2;
/// Cap on one dimer rotation, rad (the step itself is the damped Newton angle). Not pinned.
pub const DIMER_ROTATION_ANGLE: f64 = 0.3;
/// Budget for the orthogonal-complement Rayleigh-quotient MINIMISATION (`SADDLE-005`). Each
/// iteration costs exactly one calculator evaluation. Not pinned by any requirement.
pub const RAYLEIGH_MAX_ITERATIONS: u32 = 24;
/// Initial and smallest arc step of that minimisation, in units of the unit sphere. Not pinned.
pub const RAYLEIGH_STEP_INITIAL: f64 = 0.25;
pub const RAYLEIGH_STEP_MIN: f64 = 1.0e-3;
/// Residual `|H v - rho v|` at which `v` is treated as stationary, eV A^-2. Not pinned.
pub const RAYLEIGH_GRADIENT_TOLERANCE: f64 = 1.0e-9;
/// Minimiser initial and maximum steepest-descent step, A per eV A^-1.
pub const RELAX_STEP_INITIAL: f64 = 0.05;
pub const RELAX_STEP_MAX: f64 = 0.2;
/// Largest atom count for which same-species permutation matching is exhaustive (`STATE-005`).
pub const EXHAUSTIVE_PERMUTATION_LIMIT: usize = 8;
/// Identity of this minimiser, recorded in `E2-ID-006` relaxation provenance.
pub const MINIMIZER_IDENTITY: &str = "spark-backtracking-steepest-descent/1";
/// Identity of this saddle solver.
pub const SEARCHER_IDENTITY: &str = "spark-directional-dimer/1";
pub const TRAJECTORY_SCHEMA: &str = "spark-atomistic-trajectory/2";
pub const SUMMARY_SCHEMA: &str = "spark-atomistic-summary/2";

fn fail(status: StatusCode, requirement: &str) -> ApiFailure { ApiFailure::new(status, requirement) }
fn detail(status: StatusCode, requirement: &str, key: &str, value: Value) -> ApiFailure {
    let mut e = ApiFailure::new(status, requirement);
    e.details.insert(key.into(), value);
    e
}
/// `E3-STATUS-001` (normative; adopted by `D-127`): `E2-SCOPE-004` puts transport, executable
/// and IPC under the run-request `extension`, so a failure to bind them is the ADAPTER refusing the
/// request, not the public API surface. Revert by deleting the `.in_component("adapter")` calls.
fn adapter_field(field: &str) -> ApiFailure {
    detail(StatusCode::InvalidInput, "E2-SCOPE-004", "extension_field", json!(field)).in_component("adapter")
}

// ---------------------------------------------------------------------------------------------
// Calculator transport (`CALC-001`..`CALC-006`, `E2-SCOPE-004`)
// ---------------------------------------------------------------------------------------------

/// One energy/force evaluation result.
#[derive(Clone, Debug, PartialEq)]
pub struct Evaluation { pub energy_ev: f64, pub forces: Vec<Vec3> }

/// Process-isolated JSON-lines calculator. The argv comes from run-request
/// `extension.calculator_command` and is never recorded anywhere (`E2-SCOPE-004`).
pub struct CalculatorProcess {
    argv: Vec<String>,
    timeout: Duration,
    model_digest: String,
    deterministic: bool,
    identity: String,
}

impl CalculatorProcess {
    /// `E2-SCOPE-004`: "Adapter-specific transport, executable, IPC, environment,
    /// timeout-enforcement mechanism, and process-isolation settings belong only under
    /// run-request `extension`." A run request without an executable therefore has no adapter
    /// binding at all, which is `INVALID_INPUT` under `E2-API-003`/`E2-SCOPE-004`.
    pub fn from_extension(model: &ValidatedModel, extension: &Map<String, Value>) -> Result<Self, ApiFailure> {
        let command = extension.get("calculator_command").ok_or_else(|| adapter_field("calculator_command"))?;
        let items = command.as_array().ok_or_else(|| adapter_field("calculator_command"))?;
        let mut argv = Vec::new();
        for item in items {
            match item.as_str() {
                Some(s) if !s.is_empty() => argv.push(s.to_owned()),
                _ => return Err(adapter_field("calculator_command")),
            }
        }
        if argv.is_empty() {
            return Err(adapter_field("calculator_command"));
        }
        let m = model.model();
        let seconds = m.resources.callback_timeout_s;
        if !(seconds.is_finite() && seconds > 0.0) {
            return Err(fail(StatusCode::InvalidInput, "E2-SCHEMA-009").in_component("schema"));
        }
        Ok(Self {
            argv,
            timeout: Duration::from_secs_f64(seconds),
            model_digest: m.calculator.model_digest.clone(),
            deterministic: m.calculator.deterministic,
            // Model-derived, never extension-derived (`E2-SCOPE-004`).
            identity: format!("{}/{}", m.calculator.model_name, m.calculator.model_version),
        })
    }

    pub fn identity(&self) -> &str { &self.identity }

    /// `CALC-001`: the request carries the complete state contract and the requested properties.
    fn request_bytes(&self, system: &WireSystem) -> Result<Vec<u8>, ApiFailure> {
        let state = serde_json::to_value(system).map_err(|_| fail(StatusCode::InternalError, "CALC-001"))?;
        let core = json!({"requested_properties": ["energy", "forces"], "state": state});
        let digest = format!("sha256:{}", hex_sha256(&canonical_json_bytes(&core).map_err(|_| fail(StatusCode::NonfiniteResult, "STATE-008"))?));
        let full = json!({
            "request_digest": digest,
            "requested_properties": ["energy", "forces"],
            "state": core.get("state").cloned().unwrap_or(Value::Null),
        });
        let mut bytes = canonical_json_bytes(&full).map_err(|_| fail(StatusCode::NonfiniteResult, "STATE-008"))?;
        bytes.push(b'\n');
        Ok(bytes)
    }

    /// One process-isolated evaluation. `CALC-006`: timeout, process failure, malformed shape,
    /// unit mismatch, or model-digest change returns `CALCULATOR_FAILURE`; partial data is discarded.
    pub fn evaluate(&self, system: &WireSystem) -> Result<Evaluation, ApiFailure> {
        let request = self.request_bytes(system)?;
        let mut child = Command::new(&self.argv[0])
            .args(&self.argv[1..])
            .stdin(Stdio::piped()).stdout(Stdio::piped()).stderr(Stdio::piped())
            .spawn()
            .map_err(|e| detail(StatusCode::CalculatorFailure, "CALC-006", "spawn_error", json!(e.to_string())))?;
        {
            let mut stdin = child.stdin.take().ok_or_else(|| fail(StatusCode::CalculatorFailure, "CALC-006"))?;
            if stdin.write_all(&request).and_then(|()| stdin.flush()).is_err() {
                let _ = child.kill();
                let _ = child.wait();
                return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "transport", json!("stdin write failed")));
            }
        }
        let mut stdout = match child.stdout.take() {
            Some(stdout) => stdout,
            None => {
                let _ = child.kill();
                let _ = child.wait();
                return Err(fail(StatusCode::CalculatorFailure, "CALC-006"));
            }
        };
        let mut stderr = match child.stderr.take() {
            Some(stderr) => stderr,
            None => {
                let _ = child.kill();
                let _ = child.wait();
                return Err(fail(StatusCode::CalculatorFailure, "CALC-006"));
            }
        };
        // Drain both pipes concurrently. Ignoring stderr can fill its OS pipe and deadlock an
        // otherwise valid calculator before stdout closes. The process deadline is polled here,
        // rather than applied only to stdout, because a child may close stdout and keep running.
        let stdout_reader = thread::spawn(move || {
            let mut buffer = Vec::new();
            stdout.read_to_end(&mut buffer).map(|_| buffer)
        });
        let stderr_reader = thread::spawn(move || std::io::copy(&mut stderr, &mut std::io::sink()));
        let deadline = Instant::now() + self.timeout;
        let exit = loop {
            match child.try_wait() {
                Ok(Some(exit)) => break exit,
                Ok(None) if Instant::now() < deadline => {
                    thread::sleep(Duration::from_millis(1).min(deadline.saturating_duration_since(Instant::now())));
                }
                Ok(None) => {
                    let _ = child.kill();
                    let _ = child.wait();
                    let _ = stdout_reader.join();
                    let _ = stderr_reader.join();
                    return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "transport", json!("callback timeout")));
                }
                Err(e) => {
                    let _ = child.kill();
                    let _ = child.wait();
                    let _ = stdout_reader.join();
                    let _ = stderr_reader.join();
                    return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "transport", json!(e.to_string())));
                }
            }
        };
        let received = match stdout_reader.join() {
            Ok(Ok(buffer)) => buffer,
            Ok(Err(e)) => {
                let _ = stderr_reader.join();
                let _ = child.kill();
                let _ = child.wait();
                return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "transport", json!(e.to_string())));
            }
            Err(_) => {
                let _ = stderr_reader.join();
                return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "transport", json!("stdout reader panicked")));
            }
        };
        match stderr_reader.join() {
            Ok(Ok(_)) => {}
            Ok(Err(e)) => return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "transport", json!(e.to_string()))),
            Err(_) => return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "transport", json!("stderr reader panicked"))),
        }
        if !exit.success() {
            return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "exit_status", json!(exit.code())));
        }
        self.decode(system, &received)
    }

    /// `CALC-002`: status, energy, forces N x 3, units, model name/version/digest, evaluation ID,
    /// deterministic flag. `STATE-008`: nonfinite values are `NONFINITE_RESULT`, never clipped.
    fn decode(&self, system: &WireSystem, bytes: &[u8]) -> Result<Evaluation, ApiFailure> {
        let line: &[u8] = match bytes.iter().position(|b| *b == b'\n') { Some(i) => &bytes[..i], None => bytes };
        let value = parse_strict_json(line).map_err(|_| detail(StatusCode::CalculatorFailure, "CALC-006", "transport", json!("response is not strict JSON")))?;
        let object = value.as_object().ok_or_else(|| fail(StatusCode::CalculatorFailure, "CALC-006"))?;
        let text = |k: &str| object.get(k).and_then(Value::as_str).unwrap_or_default().to_owned();
        if text("status") != "OK" { return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "calculator_status", json!(text("status")))); }
        if text("model_digest") != self.model_digest { return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "reason", json!("model digest changed"))); }
        if text("model_name").is_empty() || text("model_version").is_empty() || text("evaluation_id").is_empty() {
            return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "reason", json!("incomplete calculator metadata")));
        }
        let units = object.get("units").and_then(Value::as_object).ok_or_else(|| fail(StatusCode::CalculatorFailure, "CALC-006"))?;
        if units.get("energy").and_then(Value::as_str) != Some("eV")
            || units.get("forces").and_then(Value::as_str) != Some("eV/angstrom") {
            return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "reason", json!("unit mismatch")));
        }
        // `CALC-004`/`DET-004`: a model that declares determinism must not be answered by a
        // response that disclaims it, because bitwise replay claims would silently become invalid.
        let deterministic = object.get("deterministic").and_then(Value::as_bool)
            .ok_or_else(|| fail(StatusCode::CalculatorFailure, "CALC-006"))?;
        if self.deterministic && !deterministic {
            return Err(detail(StatusCode::CalculatorFailure, "CALC-004", "reason", json!("calculator disclaims determinism")));
        }
        let energy = object.get("energy").and_then(Value::as_f64).ok_or_else(|| fail(StatusCode::CalculatorFailure, "CALC-006"))?;
        let rows = object.get("forces").and_then(Value::as_array).ok_or_else(|| fail(StatusCode::CalculatorFailure, "CALC-006"))?;
        if rows.len() != system.positions.len() {
            return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "reason", json!("malformed force shape")));
        }
        let mut forces = Vec::with_capacity(rows.len());
        for row in rows {
            let triple = row.as_array().ok_or_else(|| fail(StatusCode::CalculatorFailure, "CALC-006"))?;
            if triple.len() != 3 { return Err(detail(StatusCode::CalculatorFailure, "CALC-006", "reason", json!("malformed force shape"))); }
            let mut v = [0.0_f64; 3];
            for (k, x) in triple.iter().enumerate() {
                v[k] = x.as_f64().ok_or_else(|| fail(StatusCode::CalculatorFailure, "CALC-006"))?;
            }
            forces.push(v);
        }
        if !energy.is_finite() || !forces.iter().flatten().all(|x| x.is_finite()) {
            return Err(fail(StatusCode::NonfiniteResult, "STATE-008"));
        }
        Ok(Evaluation { energy_ev: energy, forces })
    }
}

/// One budgeted calculator evaluation. `RES-001`/`RES-002`: the reservation is taken before the
/// callback and a limit hit aborts the uncommitted transaction.
fn charged_evaluate(
    calc: &CalculatorProcess, system: &WireSystem, ledger: &mut ResourceLedger,
    relaxation: Option<&str>, saddle: Option<&str>,
) -> Result<Evaluation, ApiFailure> {
    ledger.reserve_calculator(1, relaxation, saddle)
        .map_err(|s| fail(s.status, &s.context.requirement_id))?;
    calc.evaluate(system)
}

// ---------------------------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------------------------

fn max_movable_force(forces: &[Vec3], movable: &[bool]) -> f64 {
    let mut max = 0.0_f64;
    for (i, f) in forces.iter().enumerate() {
        if movable[i] { max = max.max((f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt()); }
    }
    max
}
fn norm(v: &[Vec3]) -> f64 { v.iter().flatten().map(|x| x * x).sum::<f64>().sqrt() }
fn scaled_add(base: &[Vec3], direction: &[Vec3], factor: f64, movable: &[bool]) -> Vec<Vec3> {
    base.iter().enumerate().map(|(i, p)| {
        if movable[i] { [p[0] + factor * direction[i][0], p[1] + factor * direction[i][1], p[2] + factor * direction[i][2]] } else { *p }
    }).collect()
}
fn masked(v: &[Vec3], movable: &[bool]) -> Vec<Vec3> {
    v.iter().enumerate().map(|(i, x)| if movable[i] { *x } else { [0.0, 0.0, 0.0] }).collect()
}
fn normalized(v: &[Vec3], movable: &[bool]) -> Option<Vec<Vec3>> {
    let m = masked(v, movable);
    let n = norm(&m);
    if !(n.is_finite() && n > 0.0) { return None; }
    Some(m.iter().map(|x| [x[0] / n, x[1] / n, x[2] / n]).collect())
}
fn dot(a: &[Vec3], b: &[Vec3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x[0] * y[0] + x[1] * y[1] + x[2] * y[2]).sum()
}
fn with_positions(system: &WireSystem, positions: Vec<Vec3>) -> WireSystem {
    let mut out = system.clone();
    out.positions = positions;
    out
}

// ---------------------------------------------------------------------------------------------
// Minimiser (`RELAX-001`..`RELAX-004`, `STATE-004`)
// ---------------------------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct Relaxed {
    pub positions: Vec<Vec3>,
    pub energy_ev: f64,
    pub forces: Vec<Vec3>,
    pub max_movable_force: f64,
    pub steps: u64,
    pub evaluations: u64,
    pub termination_reason: String,
}

/// Backtracking steepest descent on movable atoms only. `RELAX-003`: `OK` is legal only when
/// `STATE-004` passes; budget exhaustion is `RELAX_NOT_CONVERGED`.
pub fn relax(
    calc: &CalculatorProcess, system: &WireSystem, force_tolerance: f64, max_steps: u64,
    max_evaluations: u64, ledger: &mut ResourceLedger, request_id: &str,
) -> Result<Relaxed, ApiFailure> {
    let movable = system.movable.clone();
    let mut positions = system.positions.clone();
    let mut evaluations = 0_u64;
    let mut current = charged_evaluate(calc, system, ledger, Some(request_id), None)?;
    evaluations += 1;
    let mut step_size = RELAX_STEP_INITIAL;
    let mut steps = 0_u64;
    loop {
        let fmax = max_movable_force(&current.forces, &movable);
        if fmax <= force_tolerance {
            return Ok(Relaxed { positions, energy_ev: current.energy_ev, forces: current.forces,
                max_movable_force: fmax, steps, evaluations, termination_reason: "force-tolerance".into() });
        }
        if steps >= max_steps || evaluations >= max_evaluations {
            return Err(detail(StatusCode::RelaxNotConverged, "RELAX-003", "max_movable_force_ev_per_angstrom", json!(fmax)));
        }
        let direction = masked(&current.forces, &movable);
        let length = norm(&direction);
        let mut factor = step_size;
        if length * factor > RELAX_STEP_MAX { factor = RELAX_STEP_MAX / length; }
        let mut accepted = None;
        for _ in 0..8 {
            if evaluations >= max_evaluations { break; }
            let trial = scaled_add(&positions, &direction, factor, &movable);
            let probe = charged_evaluate(calc, &with_positions(system, trial.clone()), ledger, Some(request_id), None)?;
            evaluations += 1;
            if probe.energy_ev <= current.energy_ev { accepted = Some((trial, probe, factor)); break; }
            factor *= 0.5;
        }
        match accepted {
            Some((trial, probe, used)) => {
                positions = trial;
                current = probe;
                step_size = (used * 1.25).min(RELAX_STEP_MAX);
                steps += 1;
            }
            None => {
                let fmax = max_movable_force(&current.forces, &movable);
                return Err(detail(StatusCode::RelaxNotConverged, "RELAX-003", "max_movable_force_ev_per_angstrom", json!(fmax)));
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// First-order saddle solver (`SADDLE-001`..`SADDLE-007`)
// ---------------------------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct SaddlePoint {
    pub positions: Vec<Vec3>,
    pub energy_ev: f64,
    pub forces: Vec<Vec3>,
    pub curvature: f64,
    pub orthogonal_curvatures: Vec<f64>,
    pub direction: Vec<Vec3>,
    pub iterations: u64,
    pub evaluations: u64,
    pub termination_reason: String,
}

#[derive(Clone, Debug)]
pub enum SearchOutcome {
    Found(Box<SaddlePoint>),
    /// `SADDLE-007`: every attempt, including a failure, reports its termination reason and counts.
    Reject { status: StatusCode, termination_reason: String, evaluations: u64, iterations: u64 },
}

/// Minimum-mode following, `saddle_search.method = directional-dimer` (`E2-SCHEMA-006`).
/// `SADDLE-005` without a Hessian: negative curvature along the reported mode plus nonnegative
/// sampled orthogonal curvatures, and the evidence level is `DIRECTIONAL`.
pub fn directional_dimer(
    calc: &CalculatorProcess, system: &WireSystem, launch: &[Vec3], initial_mode: &[Vec3],
    search: &WireSaddleSearch, ledger: &mut ResourceLedger, search_id: &str,
) -> Result<SearchOutcome, ApiFailure> {
    let movable = system.movable.clone();
    let mut positions = launch.to_vec();
    let mut mode = match normalized(initial_mode, &movable) {
        Some(m) => m,
        None => return Ok(SearchOutcome::Reject { status: StatusCode::InvalidSaddle, termination_reason: "degenerate-initial-mode".into(), evaluations: 0, iterations: 0 }),
    };
    let mut evaluations = 0_u64;
    let mut iterations = 0_u64;
    let reject = |status, reason: &str, evaluations, iterations| {
        Ok(SearchOutcome::Reject { status, termination_reason: reason.into(), evaluations, iterations })
    };
    while iterations < search.max_iterations {
        let centre = match charged_evaluate(calc, &with_positions(system, positions.clone()), ledger, None, Some(search_id)) {
            Ok(v) => v,
            Err(e) if e.status == StatusCode::CalculatorFailure || e.status == StatusCode::NonfiniteResult =>
                return reject(StatusCode::SaddleNotFound, "calculator-domain-exceeded", evaluations, iterations),
            Err(e) => return Err(e),
        };
        evaluations += 1;
        let mut curvature = f64::NAN;
        for _ in 0..DIMER_ROTATIONS_PER_ITERATION {
            let probe_positions = scaled_add(&positions, &mode, DIMER_SEPARATION_ANGSTROM, &movable);
            let probe = match charged_evaluate(calc, &with_positions(system, probe_positions), ledger, None, Some(search_id)) {
                Ok(v) => v,
                Err(e) if e.status == StatusCode::CalculatorFailure || e.status == StatusCode::NonfiniteResult =>
                    return reject(StatusCode::SaddleNotFound, "calculator-domain-exceeded", evaluations, iterations),
                Err(e) => return Err(e),
            };
            evaluations += 1;
            // delta ~ -H n from the force difference; curvature is its projection on n.
            let delta: Vec<Vec3> = (0..positions.len()).map(|i| {
                let (a, b) = (probe.forces[i], centre.forces[i]);
                [(a[0] - b[0]) / DIMER_SEPARATION_ANGSTROM, (a[1] - b[1]) / DIMER_SEPARATION_ANGSTROM, (a[2] - b[2]) / DIMER_SEPARATION_ANGSTROM]
            }).collect();
            let delta = masked(&delta, &movable);
            curvature = -dot(&delta, &mode);
            if !curvature.is_finite() { return reject(StatusCode::InvalidSaddle, "nonfinite-curvature", evaluations, iterations); }
            let parallel = dot(&delta, &mode);
            let perpendicular: Vec<Vec3> = (0..positions.len()).map(|i| {
                [delta[i][0] - parallel * mode[i][0], delta[i][1] - parallel * mode[i][1], delta[i][2] - parallel * mode[i][2]]
            }).collect();
            let magnitude = norm(&perpendicular);
            if !(magnitude.is_finite() && magnitude > 1e-14) { break; }
            // The rotation must CONTRACT, not jitter: a constant-angle rotation leaves a residual
            // mode error of that angle, which makes the reported curvature and unstable direction
            // depend on the launch point and defeats `CAT-001` deduplication. The damped Newton
            // angle |g| / (2|C|) drives the perpendicular curvature gradient to zero.
            let angle = (0.5 * magnitude / curvature.abs().max(1e-300)).min(DIMER_ROTATION_ANGLE);
            if !(angle.is_finite() && angle > 0.0) { break; }
            let rotated: Vec<Vec3> = (0..positions.len()).map(|i| {
                [mode[i][0] + angle * perpendicular[i][0] / magnitude,
                 mode[i][1] + angle * perpendicular[i][1] / magnitude,
                 mode[i][2] + angle * perpendicular[i][2] / magnitude]
            }).collect();
            match normalized(&rotated, &movable) { Some(m) => mode = m, None => break }
        }
        let fmax = max_movable_force(&centre.forces, &movable);
        if curvature <= -search.curvature_tolerance && fmax <= search.force_tolerance {
            let orthogonal = match rayleigh_minimum_over_complement(
                calc, system, &positions, &centre.forces, &mode, &movable, search, ledger, search_id,
                &mut evaluations)? {
                Ok(x) => x,
                Err(reason) => return reject(reason.0, reason.1, evaluations, iterations),
            };
            if orthogonal.is_empty() {
                return reject(StatusCode::InvalidSaddle, "no-orthogonal-evidence", evaluations, iterations);
            }
            // `SADDLE-005`: more than one unstable direction is not a first-order saddle. The last
            // entry is the MINIMISED Rayleigh quotient over the whole orthogonal complement, so it
            // is the binding one; the earlier entries are the deterministic seeds it started from.
            if orthogonal.iter().any(|x| *x < -search.curvature_tolerance) {
                return reject(StatusCode::InvalidSaddle, "additional-negative-curvature", evaluations, iterations);
            }
            return Ok(SearchOutcome::Found(Box::new(SaddlePoint {
                positions, energy_ev: centre.energy_ev, forces: centre.forces, curvature,
                orthogonal_curvatures: orthogonal, direction: mode, iterations, evaluations,
                termination_reason: "curvature-and-force-tolerance".into(),
            })));
        }
        let parallel = dot(&centre.forces, &mode);
        let effective: Vec<Vec3> = if curvature < 0.0 {
            (0..positions.len()).map(|i| {
                [centre.forces[i][0] - 2.0 * parallel * mode[i][0],
                 centre.forces[i][1] - 2.0 * parallel * mode[i][1],
                 centre.forces[i][2] - 2.0 * parallel * mode[i][2]]
            }).collect()
        } else {
            (0..positions.len()).map(|i| [-parallel * mode[i][0], -parallel * mode[i][1], -parallel * mode[i][2]]).collect()
        };
        let effective = masked(&effective, &movable);
        let length = norm(&effective);
        if !length.is_finite() { return reject(StatusCode::InvalidSaddle, "nonfinite-effective-force", evaluations, iterations); }
        if length == 0.0 { return reject(StatusCode::SaddleNotFound, "stationary-effective-force", evaluations, iterations); }
        let mut factor = DIMER_TRANSLATION_STEP;
        if length * factor > DIMER_MAX_TRANSLATION_ANGSTROM { factor = DIMER_MAX_TRANSLATION_ANGSTROM / length; }
        positions = scaled_add(&positions, &effective, factor, &movable);
        if !positions.iter().flatten().all(|x| x.is_finite()) {
            return reject(StatusCode::SaddleNotFound, "translation-diverged", evaluations, iterations);
        }
        iterations += 1;
    }
    reject(StatusCode::SaddleNotFound, "iteration-budget-exhausted", evaluations, iterations)
}

/// `SADDLE-005` without a Hessian: "negative curvature along the reported mode plus nonnegative
/// sampled orthogonal curvatures". A FIXED set of sampled Rayleigh quotients cannot support
/// `E2-EVENT-005`'s mandatory `unstable_mode_count = 1`: each quotient `v^T H v` is a positively
/// weighted average of the orthogonal spectrum in that direction, so a Hessian whose restriction to
/// the complement has negative eigenvalues can still return only nonnegative values on any chosen
/// basis (e.g. `[[1,3],[3,1]]` has diagonal `1, 1` and eigenvalues `4, -2`). This function replaces
/// the fixed sample with a MINIMISATION of the Rayleigh quotient over the orthogonal complement:
/// backtracking descent along the projected Rayleigh gradient `H v - rho v`, seeded from exactly
/// the deterministic probes the fixed sampler used, so the reported value can only ever be lower
/// than the pre-change one. Each iteration costs one calculator evaluation, and the same evaluation
/// supplies both `H v` and `rho`.
///
/// What it establishes and what it does not: the returned minimum is an UPPER bound on the smallest
/// eigenvalue of the Hessian restricted to the complement (Courant-Fischer), so a negative value
/// refutes index 1, while a nonnegative value is evidence and not a proof. That is exactly the
/// weaker `DIRECTIONAL` evidence level `SADDLE-005` names, and `E2-EVENT-002` keeps the record
/// labelled `DIRECTIONAL` rather than `HESSIAN`.
///
/// Returns the seed quotients followed by the minimised quotient, or the rejection to report.
#[allow(clippy::too_many_arguments)]
fn rayleigh_minimum_over_complement(
    calc: &CalculatorProcess, system: &WireSystem, positions: &[Vec3], centre_forces: &[Vec3],
    mode: &[Vec3], movable: &[bool], search: &WireSaddleSearch, ledger: &mut ResourceLedger,
    search_id: &str, evaluations: &mut u64,
) -> Result<Result<Vec<f64>, (StatusCode, &'static str)>, ApiFailure> {
    // One finite-difference Hessian-vector product: `delta = (F(x + eps v) - F(x)) / eps = -H v`.
    let mut hessian_product = |direction: &[Vec3], evaluations: &mut u64|
        -> Result<Result<(Vec<Vec3>, f64), (StatusCode, &'static str)>, ApiFailure> {
        let probe_positions = scaled_add(positions, direction, DIMER_SEPARATION_ANGSTROM, movable);
        let probe = match charged_evaluate(calc, &with_positions(system, probe_positions), ledger, None, Some(search_id)) {
            Ok(v) => v,
            Err(e) if e.status == StatusCode::CalculatorFailure || e.status == StatusCode::NonfiniteResult =>
                return Ok(Err((StatusCode::SaddleNotFound, "calculator-domain-exceeded"))),
            Err(e) => return Err(e),
        };
        *evaluations += 1;
        let hv: Vec<Vec3> = (0..positions.len()).map(|i| {
            let (a, b) = (probe.forces[i], centre_forces[i]);
            [-(a[0] - b[0]) / DIMER_SEPARATION_ANGSTROM,
             -(a[1] - b[1]) / DIMER_SEPARATION_ANGSTROM,
             -(a[2] - b[2]) / DIMER_SEPARATION_ANGSTROM]
        }).collect();
        let hv = masked(&hv, movable);
        let rho = dot(&hv, direction);
        if !rho.is_finite() || !hv.iter().flatten().all(|x| x.is_finite()) {
            return Ok(Err((StatusCode::InvalidSaddle, "nonfinite-orthogonal-curvature")));
        }
        Ok(Ok((hv, rho)))
    };

    let seeds = orthogonal_probes(mode, movable, search.orthogonal_directions as usize);
    if seeds.is_empty() { return Ok(Ok(Vec::new())); }
    let mut recorded = Vec::new();
    let mut best: Option<(Vec<Vec3>, Vec<Vec3>, f64)> = None;
    for direction in &seeds {
        let (hv, rho) = match hessian_product(direction, evaluations)? { Ok(x) => x, Err(e) => return Ok(Err(e)) };
        recorded.push(rho);
        if best.as_ref().map_or(true, |(_, _, b)| rho < *b) { best = Some((direction.clone(), hv, rho)); }
    }
    let (mut v, mut hv, mut rho) = best.expect("at least one seed");
    let mut step = RAYLEIGH_STEP_INITIAL;
    for _ in 0..RAYLEIGH_MAX_ITERATIONS {
        // Projected Rayleigh gradient (half of it): `H v - rho v`, with the reported unstable mode
        // removed so the search never leaves the orthogonal complement.
        let mut gradient: Vec<Vec3> = (0..v.len())
            .map(|i| [hv[i][0] - rho * v[i][0], hv[i][1] - rho * v[i][1], hv[i][2] - rho * v[i][2]])
            .collect();
        let along = dot(&gradient, mode);
        for i in 0..gradient.len() { for k in 0..3 { gradient[i][k] -= along * mode[i][k]; } }
        let gradient = masked(&gradient, movable);
        let magnitude = norm(&gradient);
        if !(magnitude.is_finite() && magnitude > RAYLEIGH_GRADIENT_TOLERANCE) { break; }
        let trial: Vec<Vec3> = (0..v.len())
            .map(|i| [v[i][0] - step * gradient[i][0] / magnitude,
                      v[i][1] - step * gradient[i][1] / magnitude,
                      v[i][2] - step * gradient[i][2] / magnitude]).collect();
        let along = dot(&trial, mode);
        let trial: Vec<Vec3> = (0..trial.len())
            .map(|i| [trial[i][0] - along * mode[i][0], trial[i][1] - along * mode[i][1],
                      trial[i][2] - along * mode[i][2]]).collect();
        let trial = match normalized(&trial, movable) { Some(t) => t, None => break };
        let (trial_hv, trial_rho) = match hessian_product(&trial, evaluations)? { Ok(x) => x, Err(e) => return Ok(Err(e)) };
        if trial_rho < rho {
            v = trial; hv = trial_hv; rho = trial_rho;
        } else {
            step *= 0.5;
            if step < RAYLEIGH_STEP_MIN { break; }
        }
    }
    recorded.push(rho);
    Ok(Ok(recorded))
}

/// Deterministic orthonormal probes for the `SADDLE-005` orthogonal-curvature evidence: Cartesian
/// basis directions on movable atoms, Gram-Schmidt orthogonalised against the mode and each other,
/// taken in descending residual norm so the samples are the least degenerate available.
fn orthogonal_probes(mode: &[Vec3], movable: &[bool], count: usize) -> Vec<Vec<Vec3>> {
    let mut basis: Vec<Vec<Vec3>> = vec![mode.to_vec()];
    let mut out = Vec::new();
    let mut candidates: Vec<Vec<Vec3>> = Vec::new();
    for (i, m) in movable.iter().enumerate() {
        if !m { continue; }
        for k in 0..3 {
            let mut v = vec![[0.0_f64; 3]; mode.len()];
            v[i][k] = 1.0;
            candidates.push(v);
        }
    }
    while out.len() < count {
        let mut best: Option<(f64, Vec<Vec3>)> = None;
        for candidate in &candidates {
            let mut v = candidate.clone();
            for b in &basis {
                let p = dot(&v, b);
                for i in 0..v.len() { for k in 0..3 { v[i][k] -= p * b[i][k]; } }
            }
            let n = norm(&v);
            if !(n.is_finite() && n > 1e-8) { continue; }
            let unit: Vec<Vec3> = v.iter().map(|x| [x[0] / n, x[1] / n, x[2] / n]).collect();
            if best.as_ref().map_or(true, |(bn, _)| n > *bn) { best = Some((n, unit)); }
        }
        match best { Some((_, unit)) => { basis.push(unit.clone()); out.push(unit); } None => break }
    }
    out
}

// ---------------------------------------------------------------------------------------------
// State matching (`STATE-005`..`STATE-007`, `E2-ID-005`)
// ---------------------------------------------------------------------------------------------

fn permutations(groups: &BTreeMap<(String, bool), Vec<usize>>, n: usize) -> Vec<Vec<usize>> {
    // Cartesian product of the per-(species, movable) permutations, so every mapping is
    // species-preserving and constraint-preserving (`STATE-005`).
    let mut out: Vec<Vec<usize>> = vec![vec![usize::MAX; n]];
    for indices in groups.values() {
        let mut orders: Vec<Vec<usize>> = Vec::new();
        permute(indices, &mut vec![], &mut vec![false; indices.len()], &mut orders);
        let mut next = Vec::new();
        for base in &out {
            for order in &orders {
                let mut candidate = base.clone();
                for (slot, target) in indices.iter().zip(order) { candidate[*slot] = *target; }
                next.push(candidate);
            }
        }
        out = next;
    }
    out
}
fn permute(items: &[usize], current: &mut Vec<usize>, used: &mut Vec<bool>, out: &mut Vec<Vec<usize>>) {
    if current.len() == items.len() { out.push(current.clone()); return; }
    for i in 0..items.len() {
        if used[i] { continue; }
        used[i] = true; current.push(items[i]);
        permute(items, current, used, out);
        current.pop(); used[i] = false;
    }
}

/// Geometry verification for state equivalence. `STATE-006` requires the report to carry the atom
/// mapping, the cell-aware displacement norm, the maximum displacement and the energy difference;
/// `E2-ID-005` makes this verification, not the candidate hash, the decision.
pub fn geometry_match(
    reference_positions: &[Vec3], reference_energy: f64, candidate_positions: &[Vec3],
    candidate_energy: f64, species: &[String], movable: &[bool], cell: Cell, pbc: [bool; 3],
    kinetics: &WireKinetics,
) -> Result<Option<MatchV2>, ApiFailure> {
    let n = reference_positions.len();
    if n == 0 || candidate_positions.len() != n { return Ok(None); }
    let mut groups: BTreeMap<(String, bool), Vec<usize>> = BTreeMap::new();
    for i in 0..n { groups.entry((species[i].clone(), movable[i])).or_default().push(i); }
    let orders = if n <= EXHAUSTIVE_PERMUTATION_LIMIT { permutations(&groups, n) } else { vec![(0..n).collect()] };
    let mut best: Option<MatchV2> = None;
    for order in orders {
        // `STATE-005`: remove the whole-cell translation before comparing.
        let mut shift = [0.0_f64; 3];
        for i in 0..n {
            let delta = [candidate_positions[order[i]][0] - reference_positions[i][0],
                         candidate_positions[order[i]][1] - reference_positions[i][1],
                         candidate_positions[order[i]][2] - reference_positions[i][2]];
            let v = closest_periodic_vector(delta, cell, pbc, 1_000_000, 1e-12)
                .map_err(|_| fail(StatusCode::InvalidState, "STATE-005"))?.displacement;
            for k in 0..3 { shift[k] += v[k] / n as f64; }
        }
        let mut sum2 = 0.0_f64;
        let mut max = 0.0_f64;
        for i in 0..n {
            let delta = [candidate_positions[order[i]][0] - reference_positions[i][0] - shift[0],
                         candidate_positions[order[i]][1] - reference_positions[i][1] - shift[1],
                         candidate_positions[order[i]][2] - reference_positions[i][2] - shift[2]];
            let v = closest_periodic_vector(delta, cell, pbc, 1_000_000, 1e-12)
                .map_err(|_| fail(StatusCode::InvalidState, "STATE-005"))?.displacement;
            let d2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            sum2 += d2;
            max = max.max(d2.sqrt());
        }
        let rms = (sum2 / n as f64).sqrt();
        let energy_difference = candidate_energy - reference_energy;
        let candidate = MatchV2 { atom_mapping: order, energy_difference_ev: energy_difference,
            max_displacement_angstrom: max, rms_displacement_angstrom: rms };
        if best.as_ref().map_or(true, |b| candidate.rms_displacement_angstrom < b.rms_displacement_angstrom) {
            best = Some(candidate);
        }
    }
    let chosen = match best { Some(x) => x, None => return Ok(None) };
    let per_atom = chosen.energy_difference_ev.abs() / n as f64;
    if chosen.rms_displacement_angstrom <= kinetics.state_rms_tolerance
        && chosen.max_displacement_angstrom <= kinetics.state_max_tolerance
        && per_atom <= kinetics.state_energy_tolerance_per_atom {
        Ok(Some(chosen))
    } else {
        Ok(None)
    }
}

// ---------------------------------------------------------------------------------------------
// Engine
// ---------------------------------------------------------------------------------------------

struct Engine<'a> {
    model: &'a ValidatedModel,
    calc: CalculatorProcess,
    ledger: ResourceLedger,
    states: BTreeMap<String, CommittedStateV2>,
    events: BTreeMap<String, DirectedEventV2>,
    multiplicity: BTreeMap<String, u64>,
    discovery: BTreeMap<String, DiscoveryStatsV2>,
    trajectory: Vec<KmcStepV2>,
    rng: Philox,
    substreams: BTreeMap<String, crate::rng::PhiloxState>,
    current: String,
    initial: String,
    simulation_time_s: f64,
    step_index: u64,
    checkpoint_sequence: u64,
    last_status: StatusCode,
    cancelled: bool,
    resource_limited: bool,
    complete: bool,
}

/// Everything the run committed, for inspection by an in-process caller. It is NOT part of the
/// public wire surface and never reaches a response, a digest, or an artifact.
#[derive(Clone, Debug, Default)]
pub struct RunReport {
    pub states: BTreeMap<String, CommittedStateV2>,
    pub events: BTreeMap<String, DirectedEventV2>,
    pub multiplicity: BTreeMap<String, u64>,
    pub discovery: BTreeMap<String, DiscoveryStatsV2>,
    pub trajectory: Vec<KmcStepV2>,
    pub calculator_evaluations: u64,
    pub checkpoints_written: u64,
}

/// The complete `E2-API-003` run adapter of this crate.
#[derive(Clone, Copy, Debug, Default)]
pub struct ProcessRunAdapter;

impl RunAdapter for ProcessRunAdapter {
    fn run(&self, model: &ValidatedModel, extension: &Map<String, Value>) -> Result<RunValue, ApiFailure> {
        execute_run(model, extension).0
    }
}

/// Same execution as `RunAdapter::run`, additionally returning the committed catalog for tests.
pub fn execute_run(model: &ValidatedModel, extension: &Map<String, Value>) -> (Result<RunValue, ApiFailure>, RunReport) {
    let calc = match CalculatorProcess::from_extension(model, extension) {
        Ok(c) => c,
        Err(e) => return (Err(e), RunReport::default()),
    };
    {
        let m = model.model();
        let limits = ResourceLimits {
            wall_time_s: m.resources.wall_time_s,
            total_calculator_evaluations: m.resources.total_calculator_evaluations,
            evaluations_per_relaxation: m.resources.evaluations_per_relaxation,
            evaluations_per_saddle_attempt: m.resources.evaluations_per_saddle_attempt,
            saddle_attempts_per_state: m.resources.saddle_attempts_per_state,
            catalog_events: m.resources.catalog_events,
            resident_memory_bytes: m.resources.resident_memory_bytes,
            output_bytes: m.resources.output_bytes,
            callback_timeout_s: m.resources.callback_timeout_s,
            retry_count: 0,
            retry_backoff_s: 0.0,
        };
        let ledger = match ResourceLedger::new(limits, Default::default(), CancelToken::new()) {
            Ok(l) => l,
            Err(s) => return (Err(fail(s.status, &s.context.requirement_id)), RunReport::default()),
        };
        let mut engine = Engine {
            model, calc, ledger,
            states: BTreeMap::new(), events: BTreeMap::new(), multiplicity: BTreeMap::new(),
            discovery: BTreeMap::new(), trajectory: Vec::new(),
            rng: derive_trajectory_stream(m.kinetics.run_seed),
            substreams: BTreeMap::new(),
            current: String::new(), initial: String::new(),
            simulation_time_s: 0.0, step_index: 0, checkpoint_sequence: 0,
            last_status: StatusCode::Ok, cancelled: false, resource_limited: false, complete: false,
        };
        let outcome = engine.execute();
        // `RES-002`/`CKPT-005`: the artifacts of committed work are written on both paths, so a
        // refused checkpoint never destroys the trajectory that was already committed.
        let value = engine.run_value();
        let write = engine.write_artifacts(&outcome);
        let report = RunReport {
            states: engine.states.clone(), events: engine.events.clone(),
            multiplicity: engine.multiplicity.clone(), discovery: engine.discovery.clone(),
            trajectory: engine.trajectory.clone(),
            calculator_evaluations: engine.ledger.counters.calculator_evaluations,
            checkpoints_written: engine.checkpoint_sequence,
        };
        (match outcome { Ok(()) => write.map(|()| value), Err(e) => Err(e) }, report)
    }
}

impl Engine<'_> {
    fn kinetics(&self) -> &WireKinetics { &self.model.model().kinetics }
    fn system(&self) -> &WireSystem { &self.model.model().system }

    fn run_value(&self) -> RunValue {
        RunValue {
            checkpoint_sequence: self.checkpoint_sequence,
            current_state_id: self.current.clone(),
            incomplete_catalog: self.discovery.values().any(|d| d.permanently_incomplete_catalog),
            simulation_time_s: self.simulation_time_s,
            step_index: self.step_index,
        }
    }

    fn execute(&mut self) -> Result<(), ApiFailure> {
        let output = &self.model.model().output;
        if output.resume { self.restore()?; } else { self.bootstrap()?; }
        while self.step_index < self.kinetics().maximum_steps {
            self.discover()?;
            let stepped = self.step()?;
            if !stepped { break; }
            if self.step_index % self.model.model().output.checkpoint_every_steps == 0 {
                self.write_checkpoint()?;
            }
        }
        self.complete = true;
        self.last_status = StatusCode::Ok;
        self.write_checkpoint()
    }

    /// `STATE-003`/`STATE-004`: the initial minimum is relaxed and committed before any discovery.
    fn bootstrap(&mut self) -> Result<(), ApiFailure> {
        let system = self.system().clone();
        let relaxation = self.model.model().relaxation.clone();
        let relaxed = relax(&self.calc, &system, relaxation.force_tolerance, relaxation.max_steps,
            relaxation.max_evaluations, &mut self.ledger, "relax:initial")?;
        let state = self.commit_state(&relaxed)?;
        self.initial = state.clone();
        self.current = state;
        self.write_checkpoint()
    }

    fn commit_state(&mut self, relaxed: &Relaxed) -> Result<String, ApiFailure> {
        let system = self.system().clone();
        let kinetics = self.kinetics().clone();
        // `E2-ID-005`: a passing candidate takes the already committed state ID.
        let mut existing: Vec<&CommittedStateV2> = self.states.values().collect();
        existing.sort_by(|a, b| a.state_id.cmp(&b.state_id));
        for state in existing {
            if geometry_match(&state.positions, state.energy_ev, &relaxed.positions, relaxed.energy_ev,
                &system.species, &system.movable, system.cell, system.pbc, &kinetics)?.is_some() {
                return Ok(state.state_id.clone());
            }
        }
        let placed = with_positions(&system, relaxed.positions.clone());
        let (candidate, state_id) = state_ids(&placed, relaxed.energy_ev)?;
        let record = CommittedStateV2 {
            atom_ids: placed.atom_ids.clone(),
            calculator_model_digest: placed.calculator_model_digest.clone(),
            candidate_identity: candidate,
            cell: placed.cell,
            charge: placed.charge,
            constraint_digest: constraint_digest(&placed)?,
            constraints: placed.constraints.clone(),
            energy_ev: relaxed.energy_ev,
            fixed_contract_digest: fixed_contract_digest(&placed)?,
            force_tolerance_ev_per_angstrom: self.model.model().relaxation.force_tolerance,
            forces_ev_per_angstrom: relaxed.forces.clone(),
            identity_version: IDENTITY_VERSION_V2.into(),
            max_movable_force_ev_per_angstrom: relaxed.max_movable_force,
            movable: placed.movable.clone(),
            pbc: placed.pbc,
            positions: relaxed.positions.clone(),
            relaxation_provenance: RelaxationProvenanceV2 {
                calculator_evaluations: relaxed.evaluations,
                calculator_identity: self.calc.identity().to_owned(),
                minimizer_identity: MINIMIZER_IDENTITY.into(),
                steps: relaxed.steps,
                termination_reason: relaxed.termination_reason.clone(),
            },
            schema: "spark-atomistic-state/2".into(),
            species: placed.species.clone(),
            spin: placed.spin,
            state_id: state_id.clone(),
        };
        self.states.insert(state_id.clone(), record);
        Ok(state_id)
    }

    /// `DISC-001`..`DISC-006` for the current state.
    fn discover(&mut self) -> Result<(), ApiFailure> {
        let state_id = self.current.clone();
        let discovery = self.model.model().discovery.clone();
        let kinetics = self.kinetics().clone();
        let mut stats = self.discovery.remove(&state_id).unwrap_or_else(|| DiscoveryStatsV2 {
            alpha: discovery.alpha, alpha_calibration: discovery.alpha_calibration.clone(),
            attempts: 0,
            // `E2-DISC-004` names `config_digest` but not which config; `E2-CAN-007` is the only
            // definition of the term ("hashes the complete validated model after removing root
            // `metadata`"), so that is what is stored. Pre-Errata-3 this was the digest of the
            // `discovery` object alone; see `E3-EVENT-001`.
            config_digest: self.model.config_digest().to_owned(),
            consecutive_redundant_successes: 0, duplicates: 0, evaluations: 0,
            event_log_rates: BTreeMap::new(), failures_by_status: BTreeMap::new(),
            heuristic_confidence: ConfidenceV2::Unavailable("UNAVAILABLE".into()),
            permanently_incomplete_catalog: false, relevance_rate_min: discovery.relevance_rate_min,
            state_id: state_id.clone(), stopping_state: "RUNNING".into(), successes: 0,
        });
        if stats.stopping_state == "CONVERGED_HEURISTIC" { self.discovery.insert(state_id, stats); return Ok(()); }
        let before = self.ledger.counters.calculator_evaluations;
        let mut result = Ok(());
        while stats.attempts < discovery.maximum_attempts {
            if stats.successes >= discovery.minimum_successful
                && stats.consecutive_redundant_successes >= discovery.consecutive_redundant { break; }
            if self.ledger.counters.calculator_evaluations.saturating_sub(before) >= discovery.maximum_evaluations { break; }
            let index = stats.attempts;
            match self.attempt(&state_id, index, &mut stats) {
                Ok(()) => {}
                Err(e) => { result = Err(e); break; }
            }
        }
        stats.evaluations = self.ledger.counters.calculator_evaluations.saturating_sub(before);
        let converged = stats.successes >= discovery.minimum_successful
            && stats.consecutive_redundant_successes >= discovery.consecutive_redundant;
        stats.stopping_state = if converged { "CONVERGED_HEURISTIC".into() } else { "INCOMPLETE".into() };
        stats.heuristic_confidence = match (converged, discovery.alpha) {
            // `DISC-005`: C = 1 - 1/(alpha * N_r), and without a declared calibration it is
            // the exact string `UNAVAILABLE`.
            (true, Some(alpha)) if alpha > 0.0 && stats.consecutive_redundant_successes > 0 =>
                ConfidenceV2::Number(1.0 - 1.0 / (alpha * stats.consecutive_redundant_successes as f64)),
            _ => ConfidenceV2::Unavailable("UNAVAILABLE".into()),
        };
        // `DISC-006`: strict mode pauses; only exploratory mode tags the trajectory permanently.
        stats.permanently_incomplete_catalog = !converged && discovery.mode == "exploratory";
        let strict_incomplete = !converged && discovery.mode == "strict";
        self.discovery.insert(state_id.clone(), stats);
        result?;
        let _ = kinetics;
        if strict_incomplete {
            self.last_status = StatusCode::DiscoveryIncomplete;
            let mut e = fail(StatusCode::DiscoveryIncomplete, "E2-DISC-004");
            e.details.insert("state_id".into(), json!(state_id));
            // `DISC-006`/`RES-002`: the pause attempts one checkpoint. `E2-STATUS-004`: a failure
            // of that attempt is appended as context, never swallowed.
            if let Err(c) = self.write_checkpoint() {
                e.details.insert("checkpoint_attempt".into(), json!({
                    "requirement_id": c.requirement_id, "status": c.status, "details": c.details}));
            }
            if let Some(d) = self.discovery.get(&state_id) {
                e.details.insert("failures_by_status".into(),
                    serde_json::to_value(&d.failures_by_status).unwrap_or(Value::Null));
                e.details.insert("attempts".into(), json!(d.attempts));
            }
            return Err(e);
        }
        Ok(())
    }

    /// One saddle-search attempt with its scheduling-independent substream (`DISC-002`,
    /// `E2-DISC-001`, `E2-DISC-002`, `E2-RNG-004`).
    fn attempt(&mut self, state_id: &str, index: u64, stats: &mut DiscoveryStatsV2) -> Result<(), ApiFailure> {
        let choice = choose_discovery_class(self.model, state_id, index)?;
        let seed = self.kinetics().run_seed;
        let mut stream = derive_saddle_substream(seed, state_id, &choice.search_class, index)
            .map_err(|s| fail(s.status, &s.context.requirement_id))?;
        if let Err(s) = self.ledger.reserve_saddle_attempt(state_id) {
            return Err(fail(s.status, &s.context.requirement_id));
        }
        stats.attempts += 1;
        let system = self.system().clone();
        let origin = self.states.get(state_id).cloned().ok_or_else(|| fail(StatusCode::InvalidState, "E2-ID-006"))?;
        let search = self.model.model().saddle_search.clone();
        let launch = self.perturb(&origin, &choice.search_class, &mut stream, search.endpoint_displacement)?;
        let mode = self.random_direction(&mut stream, system.movable.len(), &system.movable)?;
        // `E2-CKPT-007`(8) verifies the "complete substream map", and `E2-EVENT-004`'s
        // `rng_substream_digest` is the only per-event handle on it, so the digest must identify
        // the substream state the checkpoint actually stores. Taking it BEFORE the perturbation and
        // mode draws made every catalog-bearing checkpoint of this backend unverifiable against its
        // own map; the digest is now taken from exactly the stored state.
        let substream_digest_value = substream_digest(&stream);
        self.substreams.insert(choice.search_id.clone(), stream.state());
        let outcome = directional_dimer(&self.calc, &system, &launch, &mode, &search, &mut self.ledger, &choice.search_id)?;
        let point = match outcome {
            SearchOutcome::Found(p) => *p,
            SearchOutcome::Reject { status, .. } => {
                // `E2-DISC-005`: a failed attempt increments one `failures_by_status` count and
                // resets the consecutive redundant counter.
                *stats.failures_by_status.entry(status).or_insert(0) += 1;
                stats.consecutive_redundant_successes = 0;
                return Ok(());
            }
        };
        match self.validate_and_commit(&origin, &point, &choice, &substream_digest_value, index)? {
            Commit::Committed(event_id) => {
                stats.successes += 1;
                stats.consecutive_redundant_successes = 0;
                if let Some(e) = self.events.get(&event_id) {
                    stats.event_log_rates.insert(event_id.clone(), e.rate_model.log_forward_rate_per_s);
                }
            }
            Commit::Duplicate(event_id) => {
                stats.successes += 1;
                stats.duplicates += 1;
                stats.consecutive_redundant_successes += 1;
                if let Some(e) = self.events.get(&event_id) {
                    stats.event_log_rates.insert(event_id.clone(), e.rate_model.log_forward_rate_per_s);
                }
            }
            Commit::Rejected(status) => {
                *stats.failures_by_status.entry(status).or_insert(0) += 1;
                stats.consecutive_redundant_successes = 0;
            }
        }
        Ok(())
    }

    /// `DISC-001` initial perturbation. The erratum pins the class mixture but no amplitude, so the
    /// amplitude is `saddle_search.endpoint_displacement` and that choice is declared here.
    fn perturb(&self, origin: &CommittedStateV2, class: &str, stream: &mut Philox, amplitude: f64) -> Result<Vec<Vec3>, ApiFailure> {
        let movable = origin.movable.clone();
        let kind = self.model.model().discovery.classes.iter()
            .find(|c| c.name == class).map(|c| c.kind.clone()).unwrap_or_else(|| "global".into());
        let mut mask = movable.clone();
        match kind.as_str() {
            "global" => {}
            "local" => {
                let movable_indices: Vec<usize> = movable.iter().enumerate().filter(|(_, m)| **m).map(|(i, _)| i).collect();
                let u = stream.next_uniform().map_err(|s| fail(s.status, &s.context.requirement_id))?;
                let pick = movable_indices[((u * movable_indices.len() as f64) as usize).min(movable_indices.len() - 1)];
                mask = movable.iter().enumerate().map(|(i, m)| *m && i == pick).collect();
            }
            _ => {
                // `targeted`: the movable atom carrying the largest committed force, chosen with no
                // RNG so a targeted class never perturbs the substream ordering of the others.
                let mut pick = 0;
                let mut best = f64::NEG_INFINITY;
                for (i, f) in origin.forces_ev_per_angstrom.iter().enumerate() {
                    if !movable[i] { continue; }
                    let n = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
                    if n > best { best = n; pick = i; }
                }
                mask = movable.iter().enumerate().map(|(i, m)| *m && i == pick).collect();
            }
        }
        let direction = self.random_direction(stream, movable.len(), &mask)?;
        Ok(scaled_add(&origin.positions, &direction, amplitude, &mask))
    }

    fn random_direction(&self, stream: &mut Philox, n: usize, mask: &[bool]) -> Result<Vec<Vec3>, ApiFailure> {
        for _ in 0..64 {
            let mut v = vec![[0.0_f64; 3]; n];
            for (i, slot) in v.iter_mut().enumerate() {
                if !mask[i] { continue; }
                for k in 0..3 {
                    let u = stream.next_uniform().map_err(|s| fail(s.status, &s.context.requirement_id))?;
                    slot[k] = 2.0 * u - 1.0;
                }
            }
            if let Some(unit) = normalized(&v, mask) { return Ok(unit); }
        }
        Err(fail(StatusCode::InternalError, "DISC-002"))
    }

    fn validate_and_commit(
        &mut self, origin: &CommittedStateV2, point: &SaddlePoint, choice: &DiscoveryChoice,
        substream: &str, index: u64,
    ) -> Result<Commit, ApiFailure> {
        let system = self.system().clone();
        let kinetics = self.kinetics().clone();
        let relaxation = self.model.model().relaxation.clone();
        let displacement = self.model.model().saddle_search.endpoint_displacement;
        // `SADDLE-004`: relax BOTH downhill endpoints; exactly one must match the origin.
        let mut endpoints = Vec::new();
        for sign in [-1.0_f64, 1.0] {
            let launch = scaled_add(&point.positions, &point.direction, sign * displacement, &system.movable);
            let relaxed = match relax(&self.calc, &with_positions(&system, launch), relaxation.force_tolerance,
                relaxation.max_steps, relaxation.max_evaluations, &mut self.ledger, &format!("relax:{}:{}", choice.search_id, sign)) {
                Ok(r) => r,
                Err(e) if e.status == StatusCode::RelaxNotConverged => return Ok(Commit::Rejected(StatusCode::RelaxNotConverged)),
                Err(e) => return Err(e),
            };
            endpoints.push(relaxed);
        }
        let mut matches_origin = Vec::new();
        for relaxed in &endpoints {
            matches_origin.push(geometry_match(&origin.positions, origin.energy_ev, &relaxed.positions,
                relaxed.energy_ev, &system.species, &system.movable, system.cell, system.pbc, &kinetics)?);
        }
        let origin_side = match (matches_origin[0].is_some(), matches_origin[1].is_some()) {
            (false, false) => return Ok(Commit::Rejected(StatusCode::SaddleWrongBasin)),
            (true, true) => return Ok(Commit::Rejected(StatusCode::EndpointCollapsed)),
            (true, false) => 0,
            (false, true) => 1,
        };
        let origin_match = matches_origin[origin_side].clone().ok_or_else(|| fail(StatusCode::InternalError, "SADDLE-004"))?;
        let far = endpoints.remove(1 - origin_side);
        let destination_id = self.commit_state(&far)?;
        if destination_id == origin.state_id { return Ok(Commit::Rejected(StatusCode::EndpointCollapsed)); }
        let destination = self.states.get(&destination_id).cloned().ok_or_else(|| fail(StatusCode::InternalError, "E2-ID-006"))?;
        let destination_match = geometry_match(&destination.positions, destination.energy_ev, &far.positions,
            far.energy_ev, &system.species, &system.movable, system.cell, system.pbc, &kinetics)?
            .ok_or_else(|| fail(StatusCode::InternalError, "E2-ID-005"))?;
        let saddle = SaddleV2 {
            curvature_ev_per_angstrom2: point.curvature,
            energy_ev: point.energy_ev,
            evaluation_count: point.evaluations,
            evidence_level: "DIRECTIONAL".into(),
            forces_ev_per_angstrom: point.forces.clone(),
            orthogonal_curvatures_ev_per_angstrom2: point.orthogonal_curvatures.clone(),
            positions: point.positions.clone(),
            search_id: choice.search_id.clone(),
            termination_reason: point.termination_reason.clone(),
            unstable_direction: point.direction.clone(),
        };
        // `CAT-001`/`CAT-003`: a passing duplicate must not add its rate twice.
        if let Some(existing) = self.duplicate_of(&origin.state_id, &destination_id, &saddle, &kinetics)? {
            let reverse = self.events.get(&existing).map(|e| e.reverse_event_id.clone());
            *self.multiplicity.entry(existing.clone()).or_insert(0) += 1;
            if let Some(r) = reverse { *self.multiplicity.entry(r).or_insert(0) += 1; }
            return Ok(Commit::Duplicate(existing));
        }
        let mapping: Vec<[usize; 2]> = {
            let mut v: Vec<[usize; 2]> = (0..system.positions.len()).map(|i| [i, destination_match.atom_mapping[i]]).collect();
            v.sort();
            v
        };
        let saddle_system = with_positions(&system, point.positions.clone());
        let saddle_geometry = geometry_certificate(&saddle_system)?;
        let saddle_digest = format!("sha256:{}", hex_sha256(&canonical_json_bytes(&saddle_geometry)
            .map_err(|_| fail(StatusCode::NonfiniteResult, "E2-ID-007"))?));
        let (pair_id, forward_id, reverse_id) = event_ids(&origin.state_id, &destination_id, &saddle, &saddle_digest, &mapping)?;
        let rates = common_prefactor_pair(origin.energy_ev, destination.energy_ev, saddle.energy_ev,
            kinetics.temperature, kinetics.prefactor, kinetics.barrier_tolerance, kinetics.detailed_balance_tolerance);
        let rates = match rates {
            Ok(r) => r,
            Err(s) if s.status == StatusCode::RateInvalid || s.status == StatusCode::DetailedBalanceViolation =>
                return Ok(Commit::Rejected(s.status)),
            Err(s) => return Err(fail(s.status, &s.context.requirement_id)),
        };
        let provenance = DiscoveryProvenanceV2 {
            rng_substream_digest: substream.to_owned(),
            search_class: choice.search_class.clone(),
            search_id: choice.search_id.clone(),
            search_index: index,
        };
        let barrier = saddle.energy_ev - origin.energy_ev;
        let reverse_barrier = saddle.energy_ev - destination.energy_ev;
        let cutoff = kinetics.log_rate_cutoff;
        let validation = |o: MatchV2, d: MatchV2| EventValidationV2 {
            calculator_model_digest: origin.calculator_model_digest.clone(),
            constraint_digest: origin.constraint_digest.clone(),
            destination_match: d,
            full_endpoint_relaxations: true,
            // `E2-EVENT-005` enumerates the `validation` fields but constrains only
            // `full_endpoint_relaxations` and `unstable_mode_count`, so `method` is unstated. It
            // names the VALIDATION method, which `SADDLE-004` fixes ("Validation MUST relax both
            // sides"), and not the saddle curvature evidence, which is a different requirement
            // (`SADDLE-005`) and a different field (`saddle.evidence_level` = `DIRECTIONAL`).
            // Conflating the two was this crate's own unstated constraint; see `E3-EVENT-001`
            // (normative; adopted by `D-127`). The token matches the Python backend so that
            // `E2-PAR-003` byte identity is attainable for a catalog-bearing checkpoint.
            method: "full-endpoint-relaxation/1".into(),
            origin_match: o,
            unstable_mode_count: 1,
        };
        let reverse_mapping: Vec<[usize; 2]> = {
            let mut v: Vec<[usize; 2]> = mapping.iter().map(|x| [x[1], x[0]]).collect();
            v.sort();
            v
        };
        let forward = DirectedEventV2 {
            schema: "spark-atomistic-directed-event/2".into(),
            event_id: forward_id.clone(), reverse_event_id: reverse_id.clone(), pair_id: pair_id.clone(),
            origin_state_id: origin.state_id.clone(), destination_state_id: destination_id.clone(),
            saddle: saddle.clone(), barrier_ev: barrier, reverse_barrier_ev: reverse_barrier,
            rate_model: RateModelV2 {
                common_prefactor_per_s: kinetics.prefactor,
                detailed_balance_residual: rates.detailed_balance_residual,
                log_forward_rate_per_s: rates.log_forward_rate_per_s,
                log_reverse_rate_per_s: rates.log_reverse_rate_per_s,
                model: "COMMON_PREFACTOR".into(), temperature_k: kinetics.temperature,
            },
            selectable: rates.log_forward_rate_per_s >= cutoff,
            active_atom_mapping: mapping.clone(),
            environment_key: "disabled".into(), environment_version: "none/1".into(),
            discovery_provenance: provenance.clone(),
            validation: validation(origin_match.clone(), destination_match.clone()),
            calculator_digest: self.model.model().calculator.model_digest.clone(),
            identity_digest: self.model.identity_digest().to_owned(),
            schema_digest: SCHEMA_DIGEST.into(),
            tolerance_digest: self.model.tolerance_digest().to_owned(),
        };
        let mut reverse = forward.clone();
        // `E3-EVENT-001` option 4: the reciprocal record carries
        // the SAME physical saddle mapped through `validation.destination_match.atom_mapping` with
        // `unstable_direction` negated, because the reciprocal transition leaves the saddle along
        // the opposite branch. `E2-EVENT-006` canonicalises the unstable-direction sign before
        // hashing, so `pair_id` is unchanged by the negation. Pre-Errata-3 this record carried a
        // byte-identical `saddle`; revert by deleting this block.
        reverse.saddle = mapped_saddle(&saddle, &destination_match.atom_mapping)
            .ok_or_else(|| fail(StatusCode::InternalError, "E2-EVENT-005"))?;
        reverse.event_id = reverse_id.clone();
        reverse.reverse_event_id = forward_id.clone();
        reverse.origin_state_id = destination_id.clone();
        reverse.destination_state_id = origin.state_id.clone();
        reverse.barrier_ev = reverse_barrier;
        reverse.reverse_barrier_ev = barrier;
        reverse.rate_model.log_forward_rate_per_s = rates.log_reverse_rate_per_s;
        reverse.rate_model.log_reverse_rate_per_s = rates.log_forward_rate_per_s;
        reverse.rate_model.detailed_balance_residual = -rates.detailed_balance_residual;
        reverse.selectable = rates.log_reverse_rate_per_s >= cutoff;
        reverse.active_atom_mapping = reverse_mapping;
        reverse.validation = validation(destination_match, origin_match);
        // `CAT-007`: the pair becomes visible only after every validation above has passed.
        if let Err(s) = self.ledger.reserve_catalog_event(2) { return Err(fail(s.status, &s.context.requirement_id)); }
        self.events.insert(forward_id.clone(), forward);
        self.events.insert(reverse_id.clone(), reverse);
        self.multiplicity.insert(forward_id.clone(), 1);
        self.multiplicity.insert(reverse_id, 1);
        Ok(Commit::Committed(forward_id))
    }

    /// `CAT-001`: origin equivalence, destination equivalence, saddle geometry, unstable direction
    /// up to sign, active-atom mapping and barrier tolerance.
    fn duplicate_of(&self, origin: &str, destination: &str, saddle: &SaddleV2, kinetics: &WireKinetics) -> Result<Option<String>, ApiFailure> {
        let system = self.system();
        for (id, event) in &self.events {
            if event.origin_state_id != origin || event.destination_state_id != destination { continue; }
            if (event.saddle.energy_ev - saddle.energy_ev).abs() > kinetics.saddle_energy_tolerance { continue; }
            let mut sum2 = 0.0_f64;
            let mut max = 0.0_f64;
            for i in 0..saddle.positions.len() {
                let delta = [saddle.positions[i][0] - event.saddle.positions[i][0],
                             saddle.positions[i][1] - event.saddle.positions[i][1],
                             saddle.positions[i][2] - event.saddle.positions[i][2]];
                let v = closest_periodic_vector(delta, system.cell, system.pbc, 1_000_000, 1e-12)
                    .map_err(|_| fail(StatusCode::InvalidState, "CAT-001"))?.displacement;
                let d2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                sum2 += d2;
                max = max.max(d2.sqrt());
            }
            let rms = (sum2 / saddle.positions.len() as f64).sqrt();
            // `CAT-001` compares the unstable direction "up to sign". `CAT-002` supplies no
            // separate tolerance for it, so the sign-canonicalised unit vectors are compared with
            // the same saddle RMS/maximum tolerances; requiring near-exact equality instead would
            // be an undeclared criterion and would retain one parallel event per attempt.
            let sign = if dot(&saddle.unstable_direction, &event.saddle.unstable_direction) < 0.0 { -1.0 } else { 1.0 };
            let mut direction_sum2 = 0.0_f64;
            let mut direction_max = 0.0_f64;
            for i in 0..saddle.unstable_direction.len() {
                let mut d2 = 0.0;
                for k in 0..3 {
                    let delta = sign * saddle.unstable_direction[i][k] - event.saddle.unstable_direction[i][k];
                    d2 += delta * delta;
                }
                direction_sum2 += d2;
                direction_max = direction_max.max(d2.sqrt());
            }
            let direction_rms = (direction_sum2 / saddle.unstable_direction.len() as f64).sqrt();
            if rms <= kinetics.saddle_rms_tolerance && max <= kinetics.saddle_max_tolerance
                && direction_rms <= kinetics.saddle_rms_tolerance && direction_max <= kinetics.saddle_max_tolerance {
                return Ok(Some(id.clone()));
            }
        }
        Ok(None)
    }

    /// One serial KMC step (`E2-KMC-001`..`E2-KMC-005`). Returns `false` when the current state is
    /// absorbing and `kinetics.absorbing_ok` allows the run to end there.
    fn step(&mut self) -> Result<bool, ApiFailure> {
        let kinetics = self.kinetics().clone();
        let rows: Vec<(String, String, f64)> = self.events.values()
            .filter(|e| e.origin_state_id == self.current && e.selectable)
            .map(|e| (e.event_id.clone(), e.destination_state_id.clone(), e.rate_model.log_forward_rate_per_s))
            .collect();
        if rows.is_empty() {
            if kinetics.absorbing_ok { self.last_status = StatusCode::NoEnabledEvent; return Ok(false); }
            return Err(fail(StatusCode::NoEnabledEvent, "E2-KMC-001"));
        }
        let snapshot = make_rate_snapshot(&self.current, &rows, kinetics.log_rate_cutoff)?;
        let table = snapshot.payload.clone();
        let (selection, time, next) = self.rng.two_uniforms_atomic().map_err(|s| fail(s.status, &s.context.requirement_id))?;
        let threshold = selection * table.total_rate_per_s;
        let mut cumulative = 0.0;
        let mut chosen = table.rates.len() - 1;
        for (j, rate) in table.rates.iter().enumerate() {
            cumulative += *rate;
            if cumulative > threshold { chosen = j; break; }
        }
        let increment = -time.ln() / table.total_rate_per_s;
        if !(increment.is_finite() && increment > 0.0) { return Err(fail(StatusCode::RateInvalid, "E2-KMC-003")); }
        // `EVENT-004`: apply the validated destination minimum, then one verification relaxation.
        let destination = self.states.get(&table.destination_state_ids[chosen]).cloned()
            .ok_or_else(|| fail(StatusCode::EventApplicationFailed, "EVENT-004"))?;
        let system = with_positions(self.system(), destination.positions.clone());
        let relaxation = self.model.model().relaxation.clone();
        let verified = match relax(&self.calc, &system, relaxation.force_tolerance, relaxation.max_steps,
            relaxation.max_evaluations, &mut self.ledger, &format!("verify:{}", table.event_ids[chosen])) {
            Ok(r) => r,
            Err(e) if e.status == StatusCode::RelaxNotConverged => return Err(fail(StatusCode::EventApplicationFailed, "EVENT-004")),
            Err(e) => return Err(e),
        };
        if geometry_match(&destination.positions, destination.energy_ev, &verified.positions, verified.energy_ev,
            &system.species, &system.movable, system.cell, system.pbc, &kinetics)?.is_none() {
            // `KMC-005`: state, time, RNG counter and step index stay at their pre-step values.
            return Err(fail(StatusCode::EventApplicationFailed, "EVENT-004"));
        }
        self.step_index += 1;
        self.trajectory.push(KmcStepV2 {
            checkpoint_sequence: self.checkpoint_sequence,
            log_sequence: self.step_index,
            post_state_id: table.destination_state_ids[chosen].clone(),
            pre_state_id: self.current.clone(),
            rate_table_snapshot: snapshot,
            selected_event_id: table.event_ids[chosen].clone(),
            selected_rate_per_s: table.rates[chosen],
            selection_uniform: selection,
            step_index: self.step_index,
            time_increment_s: increment,
            time_uniform: time,
            total_rate_per_s: table.total_rate_per_s,
        });
        self.simulation_time_s += increment;
        self.current = table.destination_state_ids[chosen].clone();
        self.rng = next;
        self.last_status = StatusCode::Ok;
        Ok(true)
    }

    fn checkpoint_payload(&self) -> Result<CheckpointPayloadV2, ApiFailure> {
        let counters = self.ledger.checkpoint_counters();
        let mut catalog = CatalogV2 {
            digest: String::new(), events: self.events.clone(), multiplicity: self.multiplicity.clone(),
            schema: "spark-atomistic-catalog/2".into(), states: self.states.clone(),
        };
        let mut value = serde_json::to_value(&catalog).map_err(|_| fail(StatusCode::InternalError, "E2-CKPT-006"))?;
        value.as_object_mut().ok_or_else(|| fail(StatusCode::InternalError, "E2-CKPT-006"))?.remove("digest");
        catalog.digest = format!("sha256:{}", hex_sha256(&canonical_json_bytes(&value)
            .map_err(|_| fail(StatusCode::NonfiniteResult, "E2-CKPT-006"))?));
        let initial = self.states.get(&self.initial).cloned().ok_or_else(|| fail(StatusCode::InternalError, "E2-CKPT-002"))?;
        let current = self.states.get(&self.current).cloned().ok_or_else(|| fail(StatusCode::InternalError, "E2-CKPT-002"))?;
        Ok(CheckpointPayloadV2 {
            basin: BasinCheckpointV2 { enabled: false, reason: "v1-disabled".into() },
            catalog,
            checkpoint_sequence: self.checkpoint_sequence + 1,
            digests: CheckpointDigestsV2 {
                config: self.model.config_digest().to_owned(),
                model: self.model.model().calculator.model_digest.clone(),
                schema: SCHEMA_DIGEST.into(),
                tolerances: self.model.tolerance_digest().to_owned(),
            },
            discovery_statistics: self.discovery.clone(),
            flags: CheckpointFlagsV2 {
                cancelled: self.cancelled, complete: self.complete,
                incomplete_catalog: self.discovery.values().any(|d| d.permanently_incomplete_catalog),
                last_status: self.last_status, resource_limited: self.resource_limited,
            },
            initial_state: initial,
            log_sequence: self.step_index,
            resources: CheckpointResourcesV2 {
                calculator_evaluations: counters.calculator_evaluations,
                catalog_events: counters.catalog_events,
                output_bytes: counters.output_bytes,
                resident_memory_bytes: counters.resident_memory_bytes,
                retry_history: Vec::new(),
                saddle_attempts_by_state: counters.saddle_attempts_by_state.clone(),
                wall_elapsed_s: counters.baseline_wall_elapsed_s,
            },
            rng: CheckpointRngV2 {
                run_seed: self.kinetics().run_seed,
                substream_map: self.substreams.clone(),
                trajectory: self.rng.state(),
            },
            schema: "spark-atomistic-checkpoint/2".into(),
            simulation_time_s: self.simulation_time_s,
            step_index: self.step_index,
            trajectory: self.trajectory.clone(),
            current_state: current,
        })
    }

    /// `E2-CKPT-009`/`CKPT-003`: sibling temporary, content flush, atomic replace, parent flush.
    fn write_checkpoint(&mut self) -> Result<(), ApiFailure> {
        let payload = self.checkpoint_payload()?;
        let path = self.model.resolved_paths()[0].clone();
        let has_valid = self.checkpoint_sequence > 0;
        write_checkpoint_v2(&path, payload, self.model, has_valid, &mut self.ledger)?;
        self.checkpoint_sequence += 1;
        Ok(())
    }

    /// `CKPT-004`: resume verifies every digest and payload hash before mutation.
    fn restore(&mut self) -> Result<(), ApiFailure> {
        let path = self.model.resolved_paths()[0].clone();
        let bytes = std::fs::read(&path).map_err(|e| detail(StatusCode::CheckpointCorrupt, "E2-CKPT-007", "io_error", json!(e.to_string())))?;
        let envelope = decode_checkpoint_v2(&bytes, self.model)?;
        let p = envelope.payload;
        self.states = p.catalog.states;
        self.events = p.catalog.events;
        self.multiplicity = p.catalog.multiplicity;
        self.discovery = p.discovery_statistics;
        self.trajectory = p.trajectory;
        self.rng = Philox::from_state(p.rng.trajectory).map_err(|s| fail(s.status, &s.context.requirement_id))?;
        self.substreams = p.rng.substream_map;
        self.initial = p.initial_state.state_id;
        self.current = p.current_state.state_id;
        self.simulation_time_s = p.simulation_time_s;
        self.step_index = p.step_index;
        self.checkpoint_sequence = p.checkpoint_sequence;
        self.last_status = p.flags.last_status;
        Ok(())
    }

    /// Trajectory and summary artifacts. `E2-SCHEMA-010` names both paths; NO requirement in the
    /// base specification or in either erratum fixes their content, so the records emitted here
    /// are exactly the `E2-KMC-004` committed-step records and the `E2-API-008` run value.
    fn write_artifacts(&mut self, outcome: &Result<(), ApiFailure>) -> Result<(), ApiFailure> {
        let value = self.run_value();
        let status = match outcome { Ok(()) => self.last_status, Err(e) => e.status };
        let trajectory = json!({
            "records": self.trajectory,
            "schema": TRAJECTORY_SCHEMA,
        });
        let summary = json!({
            "schema": SUMMARY_SCHEMA,
            "status": status,
            "value": serde_json::to_value(&value).map_err(|_| fail(StatusCode::InternalError, "E2-API-008"))?,
        });
        let paths = self.model.resolved_paths().clone();
        self.write_artifact(&paths[1], &trajectory)?;
        self.write_artifact(&paths[2], &summary)
    }

    fn write_artifact(&mut self, path: &Path, value: &Value) -> Result<(), ApiFailure> {
        let bytes = canonical_json_bytes(value).map_err(|_| fail(StatusCode::NonfiniteResult, "E2-CAN-004"))?;
        let reserved: u64 = bytes.len().try_into().map_err(|_| fail(StatusCode::ResourceLimit, "RES-001"))?;
        self.ledger.reserve_output(reserved).map_err(|s| fail(s.status, &s.context.requirement_id))?;
        let parent = path.parent().ok_or_else(|| fail(StatusCode::InvalidInput, "IO-003"))?;
        let name = path.file_name().and_then(|x| x.to_str()).ok_or_else(|| fail(StatusCode::InvalidInput, "IO-003"))?;
        let temporary = parent.join(format!(".{name}.tmp.{}", std::process::id()));
        let write = (|| -> std::io::Result<()> {
            let mut file = std::fs::File::create(&temporary)?;
            file.write_all(&bytes)?;
            file.sync_all()?;
            std::fs::rename(&temporary, path)?;
            std::fs::File::open(parent)?.sync_all()
        })();
        write.map_err(|e| detail(StatusCode::CheckpointCorrupt, "E2-CKPT-009", "io_error", json!(e.to_string())))
    }
}

enum Commit { Committed(String), Duplicate(String), Rejected(StatusCode) }
