/// Mean-field microkinetic ODE solver.
///
/// Solves dθ_i/dt = sum_j(ν_ij * r_j) using explicit RK4 integration
/// and Newton's method for steady-state (dθ/dt = 0).

use std::collections::HashMap;

/// A reaction in the mean-field model.
pub struct Reaction {
    pub name: String,
    /// (species_index_or_EMPTY, stoichiometry)
    /// species_index = usize::MAX means "empty" (implicit 1 - sum(θ))
    pub reactants: Vec<(usize, f64)>,
    pub products: Vec<(usize, f64)>,
    /// Forward rate constant (evaluated from parameters)
    pub rate_fwd_fn: Box<dyn Fn(&HashMap<String, f64>) -> f64>,
    /// Reverse rate constant (None = irreversible)
    pub rate_rev_fn: Option<Box<dyn Fn(&HashMap<String, f64>) -> f64>>,
    pub tof_count: HashMap<String, f64>,
}

pub struct MicroKineticModel {
    pub species: Vec<String>,
    species_map: HashMap<String, usize>,
    pub reactions: Vec<Reaction>,
    pub parameters: HashMap<String, f64>,
}

const EMPTY: usize = usize::MAX;

impl MicroKineticModel {
    pub fn new() -> Self {
        Self {
            species: Vec::new(),
            species_map: HashMap::new(),
            reactions: Vec::new(),
            parameters: HashMap::new(),
        }
    }

    pub fn add_species(&mut self, name: &str) {
        if name != "empty" && name != "*" {
            let id = self.species.len();
            self.species_map.insert(name.to_string(), id);
            self.species.push(name.to_string());
        }
    }

    pub fn add_reaction(
        &mut self,
        name: &str,
        reactants: &[(&str, f64)],
        products: &[(&str, f64)],
        rate_fwd: Box<dyn Fn(&HashMap<String, f64>) -> f64>,
        rate_rev: Option<Box<dyn Fn(&HashMap<String, f64>) -> f64>>,
        tof_count: HashMap<String, f64>,
    ) {
        let to_indices = |pairs: &[(&str, f64)]| -> Vec<(usize, f64)> {
            pairs.iter().map(|(name, stoich)| {
                if *name == "empty" || *name == "*" {
                    (EMPTY, *stoich)
                } else {
                    (self.species_map[*name], *stoich)
                }
            }).collect()
        };

        self.reactions.push(Reaction {
            name: name.to_string(),
            reactants: to_indices(reactants),
            products: to_indices(products),
            rate_fwd_fn: rate_fwd,
            rate_rev_fn: rate_rev,
            tof_count,
        });
    }

    /// Build dθ/dt vector.
    fn build_odes(&self, theta: &[f64]) -> Vec<f64> {
        let n = self.species.len();
        let mut dtheta = vec![0.0; n];

        let theta_empty = (1.0 - theta.iter().sum::<f64>()).max(0.0);

        let get_cov = |idx: usize| -> f64 {
            if idx == EMPTY { theta_empty } else { theta[idx].max(0.0) }
        };

        for rxn in &self.reactions {
            // Forward
            let kf = (rxn.rate_fwd_fn)(&self.parameters);
            let mut rf = kf;
            for &(sp, stoich) in &rxn.reactants {
                rf *= get_cov(sp).powf(stoich);
            }

            // Reverse
            let mut rr = 0.0;
            if let Some(ref rev_fn) = rxn.rate_rev_fn {
                let kr = rev_fn(&self.parameters);
                rr = kr;
                for &(sp, stoich) in &rxn.products {
                    rr *= get_cov(sp).powf(stoich);
                }
            }

            let net = rf - rr;

            for &(sp, stoich) in &rxn.reactants {
                if sp != EMPTY {
                    dtheta[sp] -= stoich * net;
                }
            }
            for &(sp, stoich) in &rxn.products {
                if sp != EMPTY {
                    dtheta[sp] += stoich * net;
                }
            }
        }

        dtheta
    }

    /// Solve ODE using RK4, return (times, coverages_matrix).
    pub fn solve_ode(
        &self,
        t_end: f64,
        dt: f64,
        theta0: Option<&[f64]>,
    ) -> (Vec<f64>, Vec<Vec<f64>>) {
        let n = self.species.len();
        let mut theta: Vec<f64> = match theta0 {
            Some(t0) => t0.to_vec(),
            None => vec![0.0; n],
        };

        let mut times = Vec::new();
        let mut trajectory = Vec::new();
        let mut t = 0.0;

        times.push(t);
        trajectory.push(theta.clone());

        while t < t_end {
            let h = dt.min(t_end - t);

            // RK4 step
            let k1 = self.build_odes(&theta);
            let y2: Vec<f64> = (0..n).map(|i| (theta[i] + 0.5 * h * k1[i]).clamp(0.0, 1.0)).collect();
            let k2 = self.build_odes(&y2);
            let y3: Vec<f64> = (0..n).map(|i| (theta[i] + 0.5 * h * k2[i]).clamp(0.0, 1.0)).collect();
            let k3 = self.build_odes(&y3);
            let y4: Vec<f64> = (0..n).map(|i| (theta[i] + h * k3[i]).clamp(0.0, 1.0)).collect();
            let k4 = self.build_odes(&y4);

            for i in 0..n {
                theta[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
                theta[i] = theta[i].clamp(0.0, 1.0);
            }

            // Renormalize if total > 1
            let total: f64 = theta.iter().sum();
            if total > 1.0 {
                for v in theta.iter_mut() {
                    *v /= total;
                }
            }

            t += h;
            times.push(t);
            trajectory.push(theta.clone());
        }

        (times, trajectory)
    }

    /// Solve for steady state using simple iteration (run ODE to large time).
    ///
    /// If `initial_theta` is provided, uses it as the starting point.
    /// Otherwise, uses a default where the first species gets most coverage.
    pub fn solve_steady_state(&self) -> HashMap<String, f64> {
        self.solve_steady_state_with(None)
    }

    pub fn solve_steady_state_with(&self, initial: Option<&[f64]>) -> HashMap<String, f64> {
        let n = self.species.len();

        let mut theta = match initial {
            Some(t0) => t0.to_vec(),
            None => {
                let mut t = vec![0.01; n];
                t[0] = 1.0 - 0.01 * (n as f64 - 1.0);
                t
            }
        };

        // Phase 1: Pre-equilibrate with adaptive Euler to get near basin
        for &(dt, steps) in &[
            (1e-15, 200), (1e-12, 200), (1e-9, 200),
            (1e-6, 200), (1e-3, 200), (1.0, 200),
        ] {
            for _ in 0..steps {
                let k1 = self.build_odes(&theta);
                for i in 0..n {
                    theta[i] += dt * k1[i];
                    theta[i] = theta[i].clamp(0.0, 1.0);
                }
                let total: f64 = theta.iter().sum();
                if total > 1.0 {
                    for v in theta.iter_mut() { *v /= total; }
                }
            }
            let dtheta = self.build_odes(&theta);
            let max_d: f64 = dtheta.iter().map(|&d| d.abs()).fold(0.0, f64::max);
            if max_d < 1e-12 { break; }
        }

        // Phase 2: Newton-Raphson to polish solution — f(θ) = dθ/dt = 0
        // J[i][j] = ∂f_i/∂θ_j  (finite-difference Jacobian)
        let eps = 1e-8;
        for _newton_iter in 0..50 {
            let f0 = self.build_odes(&theta);
            let max_f: f64 = f0.iter().map(|&v| v.abs()).fold(0.0, f64::max);
            if max_f < 1e-15 { break; }

            // Build Jacobian by finite differences
            let mut jac = vec![vec![0.0; n]; n];
            for j in 0..n {
                let h = (theta[j].abs() * eps).max(eps);
                let mut theta_p = theta.clone();
                theta_p[j] += h;
                // Clamp and renormalize perturbed state
                theta_p[j] = theta_p[j].min(1.0);
                let fp = self.build_odes(&theta_p);
                for i in 0..n {
                    jac[i][j] = (fp[i] - f0[i]) / h;
                }
            }

            // Solve J * delta = -f0  via Gaussian elimination with partial pivoting
            let delta = solve_linear_system(&jac, &f0.iter().map(|v| -v).collect::<Vec<_>>(), n);

            if let Some(delta) = delta {
                // Damped Newton step: line search with backtracking
                let mut alpha = 1.0;
                for _ in 0..10 {
                    let mut theta_new = theta.clone();
                    for i in 0..n {
                        theta_new[i] += alpha * delta[i];
                        theta_new[i] = theta_new[i].clamp(0.0, 1.0);
                    }
                    let total: f64 = theta_new.iter().sum();
                    if total > 1.0 {
                        for v in theta_new.iter_mut() { *v /= total; }
                    }
                    let f_new = self.build_odes(&theta_new);
                    let new_norm: f64 = f_new.iter().map(|v| v * v).sum::<f64>();
                    let old_norm: f64 = f0.iter().map(|v| v * v).sum::<f64>();
                    if new_norm < old_norm {
                        theta = theta_new;
                        break;
                    }
                    alpha *= 0.5;
                }
                if alpha < 1.0 / 1024.0 {
                    // Newton stalled, fall back to Euler
                    for _ in 0..500 {
                        let k = self.build_odes(&theta);
                        for i in 0..n {
                            theta[i] += 1e3 * k[i];
                            theta[i] = theta[i].clamp(0.0, 1.0);
                        }
                        let total: f64 = theta.iter().sum();
                        if total > 1.0 {
                            for v in theta.iter_mut() { *v /= total; }
                        }
                    }
                }
            } else {
                // Singular Jacobian — continue with Euler
                for _ in 0..500 {
                    let k = self.build_odes(&theta);
                    for i in 0..n {
                        theta[i] += 1e3 * k[i];
                        theta[i] = theta[i].clamp(0.0, 1.0);
                    }
                    let total: f64 = theta.iter().sum();
                    if total > 1.0 {
                        for v in theta.iter_mut() { *v /= total; }
                    }
                }
            }
        }

        let mut result = HashMap::new();
        for (i, name) in self.species.iter().enumerate() {
            result.insert(name.clone(), theta[i]);
        }
        result
    }

    /// Get TOF at given coverages.
    pub fn get_tof(&self, theta: &HashMap<String, f64>) -> HashMap<String, f64> {
        let theta_arr: Vec<f64> = self.species.iter()
            .map(|sp| *theta.get(sp).unwrap_or(&0.0))
            .collect();
        let theta_empty = (1.0 - theta_arr.iter().sum::<f64>()).max(0.0);

        let get_cov = |idx: usize| -> f64 {
            if idx == EMPTY { theta_empty } else { theta_arr[idx].max(0.0) }
        };

        let mut tof = HashMap::new();
        for rxn in &self.reactions {
            let kf = (rxn.rate_fwd_fn)(&self.parameters);
            let mut rf = kf;
            for &(sp, stoich) in &rxn.reactants {
                rf *= get_cov(sp).powf(stoich);
            }

            let mut rr = 0.0;
            if let Some(ref rev_fn) = rxn.rate_rev_fn {
                let kr = rev_fn(&self.parameters);
                rr = kr;
                for &(sp, stoich) in &rxn.products {
                    rr *= get_cov(sp).powf(stoich);
                }
            }

            let net = rf - rr;
            for (obs, &coeff) in &rxn.tof_count {
                *tof.entry(obs.clone()).or_insert(0.0) += coeff * net;
            }
        }
        tof
    }

    pub fn print_summary(&self, theta: &HashMap<String, f64>) {
        println!("{}", "=".repeat(65));
        println!("  Microkinetic Model Summary");
        println!("{}", "=".repeat(65));

        println!("\n  Parameters:");
        for (k, v) in &self.parameters {
            println!("    {} = {}", k, v);
        }

        let theta_empty = 1.0 - theta.values().sum::<f64>();
        println!("\n  Steady-state coverages:");
        println!("    {:<20} {:.6e}", "empty", theta_empty);
        let mut sorted: Vec<_> = theta.iter().collect();
        sorted.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
        for (name, val) in sorted {
            if *val > 1e-12 {
                println!("    {:<20} {:.6e}", name, val);
            }
        }

        let tof = self.get_tof(theta);
        if !tof.is_empty() {
            println!("\n  Turn-over frequencies [s^-1]:");
            for (obs, val) in &tof {
                println!("    {:<20} {:.6e}", obs, val);
            }
        }
        println!("{}", "=".repeat(65));
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  Linear solver: Gaussian elimination with partial pivoting
// ═══════════════════════════════════════════════════════════════════════

/// Solve A·x = b via Gaussian elimination with partial pivoting.
/// Returns None if the matrix is singular.
fn solve_linear_system(a: &[Vec<f64>], b: &[f64], n: usize) -> Option<Vec<f64>> {
    // Augmented matrix [A | b]
    let mut aug: Vec<Vec<f64>> = (0..n).map(|i| {
        let mut row = a[i].clone();
        row.push(b[i]);
        row
    }).collect();

    // Forward elimination
    for col in 0..n {
        // Partial pivoting: find max absolute value in column
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..n {
            let v = aug[row][col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-30 {
            return None;  // Singular
        }
        if max_row != col {
            aug.swap(col, max_row);
        }

        let pivot = aug[col][col];
        for row in (col + 1)..n {
            let factor = aug[row][col] / pivot;
            for j in col..=n {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i][n];
        for j in (i + 1)..n {
            sum -= aug[i][j] * x[j];
        }
        x[i] = sum / aug[i][i];
    }

    Some(x)
}
