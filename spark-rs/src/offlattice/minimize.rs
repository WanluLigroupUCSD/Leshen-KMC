//! Geometry minimizer — FIRE (Fast Inertial Relaxation Engine).
//!
//! Replaces the Python scipy L-BFGS-B minimizer for the OTF KMC hot loop.
//! FIRE was chosen over L-BFGS for the Rust port because:
//!   1. Pure Rust, no external solver crate needed.
//!   2. Comparable convergence to L-BFGS for KMC structural relaxation
//!      (Bitzek et al., PRL 2006, DOI 10.1103/PhysRevLett.97.170201).
//!   3. ~150 LOC vs >500 for an L-BFGS-B reimplementation.
//!   4. No history matrix, so memory is O(N) not O(N·m).
//!
//! Algorithm (Bitzek et al. 2006, with standard parameters):
//!     v_new = v_old + dt * F                          # MD-like step
//!     P     = F · v
//!     if P > 0 (downhill power):
//!         v = (1-α) v + α |v| F̂
//!         if step_count > N_min: dt = min(dt * f_inc, dt_max); α = α * f_α
//!     else:
//!         v = 0
//!         dt = dt * f_dec; α = α_start
//!         step_count = 0
//!     x_new = x_old + dt * v_new
//!
//! Convergence: max ‖F_i‖ < f_tol, where F_i is the per-atom force.

use crate::offlattice::calc::Calculator;

/// Outcome of a minimization run.
pub struct MinResult {
    pub positions: Vec<[f64; 3]>,
    pub energy: f64,
    pub converged: bool,
    pub n_steps: usize,
    pub final_max_force: f64,
}

/// FIRE minimizer parameters. Defaults match Bitzek 2006.
#[derive(Clone, Debug)]
pub struct FireParams {
    pub f_tol: f64,        // force convergence threshold (eV/Å)
    pub max_steps: usize,  // upper bound on iterations
    pub max_step: f64,     // hard clamp on per-atom displacement (Å)
    pub dt_start: f64,     // initial timestep (fs-like, in MD units the user picks)
    pub dt_max: f64,       // upper bound on dt
    pub dt_min: f64,       // lower bound on dt (after consecutive uphill events)
    pub alpha_start: f64,
    pub f_inc: f64,        // dt growth factor on downhill steps
    pub f_dec: f64,        // dt shrink factor on uphill steps
    pub f_alpha: f64,      // alpha decay factor on downhill steps
    pub n_min: usize,      // require this many downhill steps before growing dt
    pub debug: bool,
}

impl Default for FireParams {
    fn default() -> Self {
        Self {
            f_tol: 0.01,
            max_steps: 500,
            max_step: 0.2,
            dt_start: 0.05,
            dt_max:   0.5,
            dt_min:   0.001,
            alpha_start: 0.1,
            f_inc:   1.1,
            f_dec:   0.5,
            f_alpha: 0.99,
            n_min:   5,
            debug:   false,
        }
    }
}

/// Run FIRE minimization in pure Rust.
///
/// `calc` evaluates energy + gradient at each step. Frozen atoms (per
/// `calc.frozen_mask()`) are held fixed via zero force + zero velocity.
pub fn minimize<C: Calculator>(
    calc: &mut C,
    initial_positions: &[[f64; 3]],
    params: &FireParams,
) -> MinResult {
    let n = initial_positions.len();
    let mut positions: Vec<[f64; 3]> = initial_positions.to_vec();
    let mut velocities: Vec<[f64; 3]> = vec![[0.0; 3]; n];

    // Snapshot the frozen mask once. Cheap to clone since it's a Vec<bool>.
    let frozen: Option<Vec<bool>> = calc.frozen_mask().map(|m| m.to_vec());

    let mut dt = params.dt_start;
    let mut alpha = params.alpha_start;
    let mut downhill_streak: usize = 0;
    let mut last_energy = f64::NAN;

    for step in 0..params.max_steps {
        // 1. Force evaluation (the expensive part — single boundary crossing
        //    if calc is a Python callback).
        let (energy, gradient) = calc.eval(&positions);
        last_energy = energy;

        // Force = -gradient. Zero out frozen atoms.
        let mut forces: Vec<[f64; 3]> = gradient.iter().map(|g| [-g[0], -g[1], -g[2]]).collect();
        if let Some(ref m) = frozen {
            for i in 0..n {
                if m[i] {
                    forces[i] = [0.0; 3];
                    velocities[i] = [0.0; 3];
                }
            }
        }

        // 2. Convergence check — max ‖F_i‖ over free atoms.
        let max_force = forces.iter().fold(0.0_f64, |acc, f| {
            let norm = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
            if norm > acc { norm } else { acc }
        });
        if max_force < params.f_tol {
            if params.debug {
                println!("FIRE: converged in {} steps, |F|max = {:.4e}", step, max_force);
            }
            return MinResult {
                positions,
                energy,
                converged: true,
                n_steps: step,
                final_max_force: max_force,
            };
        }

        // 3. v <- v + dt * F  (semi-implicit Euler velocity update)
        for i in 0..n {
            for d in 0..3 {
                velocities[i][d] += dt * forces[i][d];
            }
        }

        // 4. FIRE mixing: P = F·v, ‖F‖, ‖v‖, then v ← (1-α) v + α |v| F̂
        let p = dot_flat(&forces, &velocities);
        let f_norm = norm_flat(&forces);
        let v_norm = norm_flat(&velocities);

        if p > 0.0 {
            // Downhill — mix v toward F direction.
            if f_norm > 0.0 {
                let scale = alpha * v_norm / f_norm;
                for i in 0..n {
                    for d in 0..3 {
                        velocities[i][d] = (1.0 - alpha) * velocities[i][d]
                                          + scale * forces[i][d];
                    }
                }
            }
            downhill_streak += 1;
            if downhill_streak > params.n_min {
                dt = (dt * params.f_inc).min(params.dt_max);
                alpha *= params.f_alpha;
            }
        } else {
            // Uphill — kill velocity, shrink dt, reset alpha.
            for i in 0..n {
                velocities[i] = [0.0; 3];
            }
            dt = (dt * params.f_dec).max(params.dt_min);
            alpha = params.alpha_start;
            downhill_streak = 0;
        }

        // 5. Position update: x ← x + dt * v, with per-atom step clamp.
        for i in 0..n {
            // Skip frozen atoms (already have v=0 and F=0).
            if let Some(ref m) = frozen {
                if m[i] { continue; }
            }
            let dxv = [
                dt * velocities[i][0],
                dt * velocities[i][1],
                dt * velocities[i][2],
            ];
            let dxv_norm = (dxv[0]*dxv[0] + dxv[1]*dxv[1] + dxv[2]*dxv[2]).sqrt();
            let scale = if dxv_norm > params.max_step {
                params.max_step / dxv_norm
            } else {
                1.0
            };
            for d in 0..3 {
                positions[i][d] += dxv[d] * scale;
            }
        }

        if params.debug && step % 20 == 0 {
            println!(
                "FIRE step {:4}: E = {:.6}, |F|max = {:.4e}, dt = {:.4e}, alpha = {:.4e}",
                step, energy, max_force, dt, alpha,
            );
        }
    }

    // Did not converge within max_steps.
    let (energy, gradient) = calc.eval(&positions);
    let max_force_final = gradient.iter().fold(0.0_f64, |acc, g| {
        let f_norm = (g[0]*g[0] + g[1]*g[1] + g[2]*g[2]).sqrt();
        if f_norm > acc { f_norm } else { acc }
    });
    MinResult {
        positions,
        energy: if energy.is_finite() { energy } else { last_energy },
        converged: false,
        n_steps: params.max_steps,
        final_max_force: max_force_final,
    }
}

// ---------------------------------------------------------------------------
//                              SCRATCH HELPERS
// ---------------------------------------------------------------------------

/// Inner product of two flat (n, 3) arrays.
fn dot_flat(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let mut s = 0.0;
    for (av, bv) in a.iter().zip(b.iter()) {
        for d in 0..3 {
            s += av[d] * bv[d];
        }
    }
    s
}

/// Euclidean norm of a flat (n, 3) array.
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

    /// Starting near the (-0.558, 1.442) Müller-Brown minimum, FIRE should
    /// converge to within tight gradient tolerance of an actual stationary
    /// point.
    #[test]
    fn fire_converges_on_muller_brown() {
        let mut mb = MullerBrown::standard();
        let start = vec![[-0.50, 1.40, 0.0]];
        let mut params = FireParams::default();
        params.f_tol = 1e-3;
        params.max_steps = 5000;
        params.dt_start = 0.001;
        params.dt_max = 0.01;
        params.max_step = 0.05;
        let result = minimize(&mut mb, &start, &params);
        assert!(
            result.converged,
            "FIRE failed to converge: max_force={:.4e}, energy={:.4}, steps={}",
            result.final_max_force, result.energy, result.n_steps
        );
        // Should land near (-0.558, 1.442)
        let dx = result.positions[0][0] - (-0.558);
        let dy = result.positions[0][1] - 1.442;
        let dist = (dx * dx + dy * dy).sqrt();
        assert!(
            dist < 0.05,
            "did not converge to the (-0.558, 1.442) minimum: got {:?}, dist={:.4}",
            result.positions[0], dist
        );
    }

    /// Multi-atom: free atom 0 (carries MB potential), frozen atom 1.
    /// Atom 1 should stay put.
    #[test]
    fn fire_respects_frozen_mask() {
        struct FrozenWrap(MullerBrown, Vec<bool>);
        impl Calculator for FrozenWrap {
            fn eval(&mut self, p: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
                self.0.eval(p)
            }
            fn frozen_mask(&self) -> Option<&[bool]> {
                Some(&self.1)
            }
        }
        let mut wrap = FrozenWrap(
            MullerBrown::standard(),
            vec![false, true],
        );
        let start = vec![[-0.50, 1.40, 0.0], [3.0, 3.0, 0.0]];
        let mut params = FireParams::default();
        params.f_tol = 1e-3;
        params.max_steps = 5000;
        params.dt_start = 0.001;
        params.dt_max = 0.01;
        params.max_step = 0.05;
        let result = minimize(&mut wrap, &start, &params);
        assert!(result.converged);
        assert_eq!(result.positions[1], [3.0, 3.0, 0.0], "frozen atom moved!");
    }
}
