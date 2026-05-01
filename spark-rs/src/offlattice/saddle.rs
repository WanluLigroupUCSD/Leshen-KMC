//! Dimer-method saddle-point search — pure-Rust port of `spark/offlattice/saddle.py`.
//!
//! Algorithm (Henkelman & Jónsson 1999, JCP 111, 7010; Heyden et al. 2005):
//!   1. Initialize a dimer (two configurations separated by `dimer_sep` along
//!      a direction `axis`).
//!   2. Rotate `axis` toward the eigenvector of lowest curvature using a
//!      finite-difference torque (`_rotate_dimer`).
//!   3. Translate along the *effective gradient* — invert the parallel
//!      component when curvature < 0 so we climb only the soft mode.
//!   4. Trust-radius adaptation, convex-region exit, collision detection
//!      against a history of known saddle points.
//!
//! Each dimer iteration costs **3 force evaluations**:
//!   - 2 inside the rotor loop (`pos ± h·axis`)
//!   - 1 at the dimer center for the effective gradient
//! Plus N additional force evals if rotor does N rotation passes per step.
//!
//! Force eval is the **only** place this module crosses the Calculator
//! boundary — everything else is `f64` arithmetic. When the Phase C
//! `PyCalculator` plugs in, the Python boundary is hit only inside `Calculator::eval`.

use crate::offlattice::calc::Calculator;

/// Tunable knobs for the dimer search. Defaults match `spark/offlattice/saddle.py`.
#[derive(Clone, Debug)]
pub struct DimerParams {
    pub dimer_sep: f64,        // half-separation h (Å)
    pub f_tol: f64,            // |g_eff| < f_tol → converged
    pub max_steps: usize,      // outer translation steps
    pub max_rotor_steps: usize,// inner rotor iterations per outer step
    pub rotor_tol: f64,        // torque norm threshold for rotor break
    pub trust: f64,            // initial trust radius (Å)
    pub trust_grow: f64,       // multiplier on a "good" step
    pub trust_shrink: f64,     // multiplier on a "bad" step
    pub trust_max: f64,        // upper cap on trust radius
    pub convex_max: usize,     // consecutive convex steps before bailout
    pub debug: bool,
}

impl Default for DimerParams {
    fn default() -> Self {
        Self {
            dimer_sep: 0.005,
            f_tol: 0.05,
            max_steps: 300,
            max_rotor_steps: 10,
            rotor_tol: 0.01,
            trust: 0.1,
            trust_grow: 1.1,
            trust_shrink: 0.5,
            trust_max: 0.5,
            convex_max: 5,
            debug: false,
        }
    }
}

/// Why the search exited.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DimerStatus {
    /// |g_eff| converged below f_tol.
    Success,
    /// Hit `convex_max` consecutive convex (curv > 0) iterations.
    Convex,
    /// Axis aligned with a previously found saddle within theta_tol.
    Collision,
    /// Outer loop exhausted `max_steps`.
    MaxIter,
}

/// Result of `DimerSearch::find_saddle`. Mirrors Python's `result` dict.
pub struct DimerResult {
    pub positions: Vec<[f64; 3]>,
    pub axis: Vec<[f64; 3]>,
    pub energy: f64,
    pub curvature: f64,
    pub status: DimerStatus,
}

/// Rotation phase output (just the rotated axis + final curvature estimate).
struct RotorOut {
    axis: Vec<[f64; 3]>,
    curvature: f64,
}

/// Run a single dimer saddle search.
///
/// Mirrors `DimerSearch.find_saddle` from saddle.py exactly except that the
/// random-axis init is delegated to the caller (caller supplies `axis`); this
/// keeps the function deterministic given inputs.
///
/// # Arguments
/// * `calc`         — Calculator (Müller-Brown, PyCalculator, etc.)
/// * `start_pos`    — initial atomic positions, length N (Å), shape (N, 3)
/// * `start_axis`   — initial dimer axis, same shape; assumed normalized
/// * `params`       — tunable knobs
/// * `history_sp`   — optional previously-found SP positions for collision check
/// * `theta_tol`    — collision angular tolerance (rad); ignored if history is None
pub fn find_saddle<C: Calculator>(
    calc: &mut C,
    start_pos: &[[f64; 3]],
    start_axis: &[[f64; 3]],
    params: &DimerParams,
    history_sp: Option<&[Vec<[f64; 3]>]>,
    theta_tol: Option<f64>,
) -> DimerResult {
    let n = start_pos.len();
    assert_eq!(start_axis.len(), n, "axis/pos length mismatch");
    let frozen: Option<Vec<bool>> = calc.frozen_mask().map(|m| m.to_vec());

    let mut pos: Vec<[f64; 3]> = start_pos.to_vec();
    let mut axis: Vec<[f64; 3]> = start_axis.to_vec();

    let mut trust = params.trust;
    let mut convex_count: usize = 0;
    let mut status = DimerStatus::MaxIter;
    let mut last_energy = f64::NAN;
    let mut last_curvature = f64::NAN;

    for step in 0..params.max_steps {
        // ---- Rotation phase: align axis with lowest-curvature mode ----
        let rot = rotate_dimer(calc, &pos, &axis, frozen.as_deref(), params);
        axis = rot.axis;
        let curvature = rot.curvature;
        last_curvature = curvature;

        // ---- Effective gradient at center ----
        let (e_center, grad_center) = calc.eval(&pos);
        last_energy = e_center;
        let mut g_eff = effective_gradient(&grad_center, &axis, curvature);
        if let Some(ref m) = frozen {
            for i in 0..n {
                if m[i] { g_eff[i] = [0.0; 3]; }
            }
        }

        // ---- Convergence check ----
        let g_norm = norm_flat(&g_eff);
        if g_norm < params.f_tol {
            status = DimerStatus::Success;
            if params.debug {
                println!("DIMER: converged step={} E={:.6} |g|={:.4e} curv={:.4}",
                         step, e_center, g_norm, curvature);
            }
            break;
        }

        // ---- Convex-region check ----
        if curvature > 0.0 {
            convex_count += 1;
            if convex_count >= params.convex_max {
                status = DimerStatus::Convex;
                if params.debug {
                    println!("DIMER: convex exit step={}", step);
                }
                break;
            }
        } else {
            convex_count = 0;
        }

        // ---- Collision detection against past SPs ----
        if let (Some(hist), Some(theta_tol)) = (history_sp, theta_tol) {
            for sp_prev in hist.iter() {
                if collision_check(&pos, sp_prev, &axis, theta_tol) {
                    status = DimerStatus::Collision;
                    if params.debug {
                        println!("DIMER: collision step={}", step);
                    }
                    return DimerResult {
                        positions: pos,
                        axis,
                        energy: e_center,
                        curvature,
                        status,
                    };
                }
            }
        }

        // ---- Translation step: descend g_eff with trust-region clamp ----
        let mut step_vec: Vec<[f64; 3]> = g_eff.iter()
            .map(|g| [-g[0], -g[1], -g[2]])
            .collect();
        let step_norm = norm_flat(&step_vec);
        if step_norm > 0.0 {
            let scale = (1.0_f64).min(trust / step_norm);
            for v in step_vec.iter_mut() {
                v[0] *= scale; v[1] *= scale; v[2] *= scale;
            }
        }
        for i in 0..n {
            for d in 0..3 {
                pos[i][d] += step_vec[i][d];
            }
        }

        // ---- Trust-radius adaptation (skip first step) ----
        if step > 0 {
            let proj = dot_flat(&step_vec, &g_eff);
            if proj < -0.5 * trust {
                trust *= params.trust_shrink;
            } else if proj > 0.5 * trust {
                trust = (trust * params.trust_grow).min(params.trust_max);
            }
        }
    }

    // Final energy snapshot (may overwrite last_energy if we broke from
    // success branch — same value).
    let energy = if last_energy.is_finite() {
        last_energy
    } else {
        calc.eval(&pos).0
    };

    DimerResult {
        positions: pos,
        axis,
        energy,
        curvature: last_curvature,
        status,
    }
}

/// Inner rotor loop. Returns the rotated axis (unit-normalized) and the final
/// curvature estimate at the rotated axis.
fn rotate_dimer<C: Calculator>(
    calc: &mut C,
    pos: &[[f64; 3]],
    initial_axis: &[[f64; 3]],
    frozen: Option<&[bool]>,
    params: &DimerParams,
) -> RotorOut {
    let n = pos.len();
    let h = params.dimer_sep;
    let mut axis: Vec<[f64; 3]> = initial_axis.to_vec();
    let mut curvature = 0.0;

    for _ in 0..params.max_rotor_steps {
        // Evaluate gradients at pos ± h·axis.
        let pos_plus: Vec<[f64; 3]> = (0..n)
            .map(|i| [
                pos[i][0] + h * axis[i][0],
                pos[i][1] + h * axis[i][1],
                pos[i][2] + h * axis[i][2],
            ])
            .collect();
        let pos_minus: Vec<[f64; 3]> = (0..n)
            .map(|i| [
                pos[i][0] - h * axis[i][0],
                pos[i][1] - h * axis[i][1],
                pos[i][2] - h * axis[i][2],
            ])
            .collect();
        let (_, g_plus) = calc.eval(&pos_plus);
        let (_, g_minus) = calc.eval(&pos_minus);

        // dg = g_plus - g_minus; curvature = (dg · axis) / (2h)
        let mut dg = vec![[0.0; 3]; n];
        for i in 0..n {
            for d in 0..3 {
                dg[i][d] = g_plus[i][d] - g_minus[i][d];
            }
        }
        curvature = dot_flat(&dg, &axis) / (2.0 * h);

        // Rotor "torque" — we want to MINIMIZE curvature κ(N̂) = dg·N̂/(2h).
        // ∇_{N̂} κ = dg/(2h), so to decrease κ we step in direction -dg/(2h),
        // projected perpendicular to N̂ (to keep |N̂|=1).
        //
        //   rotor_dir = -[dg/(2h) - (dg/(2h) · N̂) N̂]
        //
        // Equivalent: just negate the perpendicular projection of dg/(2h).
        // NB: This sign matters! Without it the rotor walks toward the
        // *highest* curvature direction (porting bug carried from openFLY
        // path through spark/offlattice/saddle.py — fixed here, fix in
        // Python landed at d5267c6+).
        let mut dg_per_2h = dg;
        for i in 0..n {
            for d in 0..3 {
                dg_per_2h[i][d] /= 2.0 * h;
            }
        }
        let proj = dot_flat(&dg_per_2h, &axis);
        let mut torque: Vec<[f64; 3]> = dg_per_2h.iter().enumerate()
            .map(|(i, v)| [
                -(v[0] - proj * axis[i][0]),
                -(v[1] - proj * axis[i][1]),
                -(v[2] - proj * axis[i][2]),
            ])
            .collect();
        if let Some(m) = frozen {
            for i in 0..n {
                if m[i] { torque[i] = [0.0; 3]; }
            }
        }
        let torque_norm = norm_flat(&torque);
        if torque_norm < params.rotor_tol {
            break;
        }

        // Rodrigues rotation toward torque direction
        let theta_dir: Vec<[f64; 3]> = torque.iter()
            .map(|t| [t[0]/torque_norm, t[1]/torque_norm, t[2]/torque_norm])
            .collect();
        let angle = (0.1_f64).min(torque_norm * 0.1);
        let (c, s) = (angle.cos(), angle.sin());

        let mut new_axis: Vec<[f64; 3]> = (0..n)
            .map(|i| [
                axis[i][0] * c + theta_dir[i][0] * s,
                axis[i][1] * c + theta_dir[i][1] * s,
                axis[i][2] * c + theta_dir[i][2] * s,
            ])
            .collect();
        if let Some(m) = frozen {
            for i in 0..n {
                if m[i] { new_axis[i] = [0.0; 3]; }
            }
        }
        let nrm = norm_flat(&new_axis);
        if nrm > 0.0 {
            for i in 0..n {
                for d in 0..3 {
                    new_axis[i][d] /= nrm;
                }
            }
            axis = new_axis;
        }
    }
    RotorOut { axis, curvature }
}

/// Effective gradient for the translation step.
///   curv < 0  → invert parallel component: g_eff = g - 2 (g·axis) axis
///   curv > 0  → climb only along axis:     g_eff = -(g·axis) axis
fn effective_gradient(
    grad: &[[f64; 3]],
    axis: &[[f64; 3]],
    curvature: f64,
) -> Vec<[f64; 3]> {
    let n = grad.len();
    let g_dot_a = dot_flat(grad, axis);
    let mut out = vec![[0.0; 3]; n];
    if curvature < 0.0 {
        for i in 0..n {
            for d in 0..3 {
                out[i][d] = grad[i][d] - 2.0 * g_dot_a * axis[i][d];
            }
        }
    } else {
        for i in 0..n {
            for d in 0..3 {
                out[i][d] = -g_dot_a * axis[i][d];
            }
        }
    }
    out
}

/// Has the dimer drifted so close to a previous SP that we should bail?
/// Returns true if the displacement (pos - sp_prev) is within `theta_tol`
/// radians of the dimer axis.
fn collision_check(
    pos: &[[f64; 3]],
    sp_prev: &[[f64; 3]],
    axis: &[[f64; 3]],
    theta_tol: f64,
) -> bool {
    let n = pos.len();
    let mut disp = vec![[0.0; 3]; n];
    for i in 0..n {
        for d in 0..3 {
            disp[i][d] = pos[i][d] - sp_prev[i][d];
        }
    }
    let d_norm = norm_flat(&disp);
    if d_norm <= 0.0 { return false; }
    let a_norm = norm_flat(axis);
    if a_norm <= 0.0 { return false; }
    let cos_theta = dot_flat(&disp, axis) / (d_norm * a_norm);
    cos_theta.abs() > theta_tol.cos()
}

// ---------------------------------------------------------------------------
//                              SCRATCH HELPERS
// ---------------------------------------------------------------------------

fn dot_flat(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let mut s = 0.0;
    for (av, bv) in a.iter().zip(b.iter()) {
        for d in 0..3 {
            s += av[d] * bv[d];
        }
    }
    s
}

fn norm_flat(a: &[[f64; 3]]) -> f64 {
    let mut s = 0.0;
    for v in a {
        for d in 0..3 {
            s += v[d] * v[d];
        }
    }
    s.sqrt()
}

// ---------------------------------------------------------------------------
//                                TESTS
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::offlattice::calc::MullerBrown;

    /// Quadratic saddle V = -x²/2 + y²/2. The saddle is at origin, x is
    /// the unstable mode. Dimer started at any point with x ≠ 0 should
    /// converge to (0, 0) cleanly. This is the simplest possible saddle
    /// test — purely analytic, exercises every code path.
    #[test]
    fn dimer_refines_quadratic_saddle() {
        use crate::offlattice::calc::QuadraticSaddle;
        let mut q = QuadraticSaddle::new(1.0, 1.0);
        let start_pos = vec![[0.3, 0.2, 0.0]];
        // Initial axis aligned ROUGHLY with unstable mode (rotor will refine)
        let start_axis = vec![[0.9, 0.4, 0.0]];
        let an: f64 = (0.81_f64 + 0.16).sqrt();
        let start_axis = vec![[start_axis[0][0]/an, start_axis[0][1]/an, 0.0]];

        let mut params = DimerParams::default();
        params.dimer_sep = 0.001;
        params.f_tol = 1e-4;
        params.max_steps = 200;
        params.trust = 0.05;
        params.trust_max = 0.2;
        params.max_rotor_steps = 30;
        params.rotor_tol = 1e-6;
        params.convex_max = 100;     // shouldn't hit on a quadratic saddle

        let result = find_saddle(
            &mut q,
            &start_pos,
            &start_axis,
            &params,
            None,
            None,
        );

        assert_eq!(
            result.status, DimerStatus::Success,
            "quadratic dimer didn't converge: status={:?}, curv={:.4}, pos={:?}",
            result.status, result.curvature, result.positions[0]
        );
        // Saddle is at origin
        let pos = result.positions[0];
        let dist: f64 = (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]).sqrt();
        assert!(
            dist < 0.01,
            "did not converge to origin: got {:?}, dist={:.4e}",
            pos, dist
        );
        // Curvature at saddle is the unstable eigenvalue: -k_x = -1.0
        assert!(
            result.curvature < -0.5,
            "curvature at saddle should be ~-1.0, got {:.4}", result.curvature
        );
    }

    /// Anisotropic quadratic saddle. Stable mode 10x stiffer than unstable.
    /// Verifies the rotor finds the correct unstable mode starting from an
    /// axis tilted MOSTLY toward the wrong (stable) direction — the rotor
    /// must rotate it ~80° to find the unstable x-mode.
    ///
    /// NB: Starting axis EXACTLY on an eigenvector (e.g. pure y) gives
    /// zero torque mathematically — production code uses random axis init
    /// to avoid this degenerate case.
    #[test]
    fn dimer_finds_correct_unstable_mode() {
        use crate::offlattice::calc::QuadraticSaddle;
        let mut q = QuadraticSaddle::new(1.0, 10.0);
        let start_pos = vec![[0.05, 0.05, 0.0]];
        // Start with axis 80% along WRONG (stable) y-direction, 20% along
        // unstable x — the rotor has to flip it.
        let raw: [f64; 3] = [0.2, 0.95, 0.0];
        let n: f64 = (raw[0] * raw[0] + raw[1] * raw[1]).sqrt();
        let start_axis = vec![[raw[0]/n, raw[1]/n, 0.0]];

        let mut params = DimerParams::default();
        params.dimer_sep = 0.001;
        params.f_tol = 1e-4;
        params.max_steps = 300;
        params.trust = 0.02;
        params.trust_max = 0.1;
        params.max_rotor_steps = 50;
        params.rotor_tol = 1e-6;
        params.convex_max = 100;

        let result = find_saddle(&mut q, &start_pos, &start_axis, &params, None, None);
        assert_eq!(
            result.status, DimerStatus::Success,
            "anisotropic dimer didn't converge: status={:?}, pos={:?}",
            result.status, result.positions[0]
        );
        // After convergence, axis should be ~aligned with x (the unstable mode)
        let ax = result.axis[0];
        let along_x: f64 = ax[0].abs();
        assert!(
            along_x > 0.95,
            "rotor failed to align with x-mode: axis={:?}",
            ax
        );
    }

    /// Effective-gradient invariants: when curv < 0, g_eff has *opposite*
    /// component along axis; when curv > 0, g_eff is purely along -axis.
    #[test]
    fn effective_gradient_branch_check() {
        let grad = vec![[1.0, 2.0, 3.0]];
        let axis = vec![[1.0, 0.0, 0.0]]; // unit x
        let g_eff_neg = effective_gradient(&grad, &axis, -1.0);
        // g·a = 1.0 → g_eff = g - 2*1.0*[1,0,0] = [-1, 2, 3]
        assert!((g_eff_neg[0][0] - (-1.0)).abs() < 1e-12);
        assert!((g_eff_neg[0][1] - 2.0).abs() < 1e-12);
        assert!((g_eff_neg[0][2] - 3.0).abs() < 1e-12);
        let g_eff_pos = effective_gradient(&grad, &axis, 1.0);
        // -(g·a) * a = [-1, 0, 0]
        assert!((g_eff_pos[0][0] - (-1.0)).abs() < 1e-12);
        assert_eq!(g_eff_pos[0][1], 0.0);
        assert_eq!(g_eff_pos[0][2], 0.0);
    }

    /// Collision check: a position offset by 1.0 along axis should trigger.
    #[test]
    fn collision_detected() {
        let pos = vec![[1.0, 0.0, 0.0]];
        let sp_prev = vec![[0.0, 0.0, 0.0]];
        let axis = vec![[1.0, 0.0, 0.0]];
        assert!(collision_check(&pos, &sp_prev, &axis, 0.1));
    }

    /// Collision check: if displacement is perpendicular to axis, no collision.
    #[test]
    fn no_collision_when_perpendicular() {
        let pos = vec![[0.0, 1.0, 0.0]];
        let sp_prev = vec![[0.0, 0.0, 0.0]];
        let axis = vec![[1.0, 0.0, 0.0]];
        assert!(!collision_check(&pos, &sp_prev, &axis, 0.1));
    }
}
