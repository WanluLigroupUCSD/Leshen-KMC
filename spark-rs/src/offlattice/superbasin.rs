//! SuperBasin — group of low-barrier basins with mean-rate-method acceleration.

use rand::Rng;
use super::basin::Basin;
use super::mechanism::Mechanism;

/// SuperBasin KMC choice result.
pub struct SBChoice {
    pub mechanism: Mechanism,
    pub atom_index: usize,
    pub dt: f64,
    pub source_basin: usize,
    pub switched: bool,
}

/// A group of basins connected by low barriers.
pub struct SuperBasin {
    pub basins: Vec<Basin>,
    pub occupied: usize,
    prob: Vec<Vec<f64>>, // Transition probability matrix
}

impl SuperBasin {
    pub fn new(initial_basin: Basin) -> Self {
        SuperBasin {
            basins: vec![initial_basin],
            occupied: 0,
            prob: vec![vec![0.0]],
        }
    }

    pub fn size(&self) -> usize {
        self.basins.len()
    }

    /// Add a new basin and occupy it.
    pub fn expand_occupy(&mut self, new_basin: Basin) {
        let n = self.size();
        self.basins.push(new_basin);

        // Expand probability matrix
        for row in &mut self.prob {
            row.push(0.0);
        }
        self.prob.push(vec![0.0; n + 1]);

        self.occupied = n;
    }

    /// Find and occupy a basin matching the state hash.
    pub fn find_occupy(
        &mut self, state_hash: u64, positions: &[[f64; 3]], tol: f64,
    ) -> Option<usize> {
        let n = positions.len();
        let centroid_x = centroid(positions);

        for (i, basin) in self.basins.iter().enumerate() {
            if state_hash == basin.state_hash {
                let delta = sub3(&centroid_x, &centroid(&basin.positions));
                let mut sum_sq = 0.0;
                for j in 0..n.min(basin.positions.len()) {
                    let diff = [
                        positions[j][0] - basin.positions[j][0] - delta[0],
                        positions[j][1] - basin.positions[j][1] - delta[1],
                        positions[j][2] - basin.positions[j][2] - delta[2],
                    ];
                    sum_sq += diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
                }
                if sum_sq < tol * tol {
                    let prev = self.occupied;
                    self.occupied = i;
                    return Some(prev);
                }
            }
        }
        None
    }

    /// Mark a mechanism as internal connection.
    pub fn connect_from(
        &mut self, basin_idx: usize, atom_index: usize, mechanism: &Mechanism,
    ) {
        let to = self.occupied;
        let from = basin_idx;

        // Search and capture the matching mechanism's rate without holding
        // a mutable borrow of self.basins past the loop.
        let mut hit_rate: Option<f64> = None;
        for lm in &mut self.basins[basin_idx].mechs {
            if lm.atom_index == atom_index
                && (lm.mechanism.barrier - mechanism.barrier).abs() < 1e-10
                && (lm.mechanism.delta - mechanism.delta).abs() < 1e-10
            {
                lm.exit_mech = false;
                hit_rate = Some(lm.rate);
                break;
            }
        }
        if let Some(rate) = hit_rate {
            let rate_sum = self.basins[basin_idx].rate_sum;
            self.basins[basin_idx].connected = true;
            self.prob[to][from] = rate / rate_sum;
        }
    }

    /// KMC selection with mean-rate method.
    pub fn kmc_choice<R: Rng>(&self, rng: &mut R) -> SBChoice {
        if !self.basins[self.occupied].connected {
            let ch = self.basins[self.occupied].kmc_choice(rng);
            return SBChoice {
                mechanism: ch.mechanism,
                atom_index: ch.atom_index,
                dt: ch.dt,
                source_basin: self.occupied,
                switched: false,
            };
        }

        let tau = self.compute_tau();
        let tau_sum: f64 = tau.iter().sum();
        if tau_sum <= 0.0 {
            let ch = self.basins[self.occupied].kmc_choice(rng);
            return SBChoice {
                mechanism: ch.mechanism,
                atom_index: ch.atom_index,
                dt: ch.dt,
                source_basin: self.occupied,
                switched: false,
            };
        }

        let tau_norm: Vec<f64> = tau.iter().map(|t| t / tau_sum).collect();

        // Collect exit mechanisms with effective rates
        let mut exits: Vec<(usize, usize, f64, &Mechanism)> = Vec::new();
        let mut r_sum = 0.0;

        for (i, basin) in self.basins.iter().enumerate() {
            for lm in &basin.mechs {
                if lm.exit_mech {
                    let eff = tau_norm[i] * lm.rate;
                    exits.push((i, lm.atom_index, eff, &lm.mechanism));
                    r_sum += eff;
                }
            }
        }

        if exits.is_empty() || r_sum <= 0.0 {
            let ch = self.basins[self.occupied].kmc_choice(rng);
            return SBChoice {
                mechanism: ch.mechanism,
                atom_index: ch.atom_index,
                dt: ch.dt,
                source_basin: self.occupied,
                switched: false,
            };
        }

        let target = rng.gen::<f64>() * r_sum;
        let mut cumul = 0.0;
        let mut sel = &exits[exits.len() - 1];

        for e in &exits {
            cumul += e.2;
            if cumul >= target {
                sel = e;
                break;
            }
        }

        let dt = -(rng.gen::<f64>().ln()) / r_sum;

        SBChoice {
            mechanism: sel.3.clone(),
            atom_index: sel.1,
            dt,
            source_basin: sel.0,
            switched: self.occupied != sel.0,
        }
    }

    /// Compute mean residence times via (I - P)^-1 @ theta.
    fn compute_tau(&self) -> Vec<f64> {
        let n = self.size();
        // Simple approach: Gaussian elimination
        // A = I - P, solve A * tau = theta
        let mut a = vec![vec![0.0; n + 1]; n]; // augmented matrix

        for i in 0..n {
            for j in 0..n {
                a[i][j] = if i == j { 1.0 } else { 0.0 };
                a[i][j] -= self.prob[i][j];
            }
            a[i][n] = if i == self.occupied { 1.0 } else { 0.0 };
        }

        // Gaussian elimination with partial pivoting
        for col in 0..n {
            let mut max_row = col;
            let mut max_val = a[col][col].abs();
            for row in (col + 1)..n {
                if a[row][col].abs() > max_val {
                    max_val = a[row][col].abs();
                    max_row = row;
                }
            }
            a.swap(col, max_row);

            if a[col][col].abs() < 1e-15 { continue; }

            for row in (col + 1)..n {
                let factor = a[row][col] / a[col][col];
                for j in col..=n {
                    a[row][j] -= factor * a[col][j];
                }
            }
        }

        // Back substitution
        let mut tau = vec![0.0; n];
        for i in (0..n).rev() {
            let mut sum = a[i][n];
            for j in (i + 1)..n {
                sum -= a[i][j] * tau[j];
            }
            tau[i] = if a[i][i].abs() > 1e-15 { sum / a[i][i] } else { 0.0 };
        }

        // Convert: tau_i / rate_sum_i
        for i in 0..n {
            if self.basins[i].rate_sum > 0.0 {
                tau[i] /= self.basins[i].rate_sum;
            }
        }

        tau
    }
}

fn centroid(positions: &[[f64; 3]]) -> [f64; 3] {
    let n = positions.len() as f64;
    let mut c = [0.0; 3];
    for p in positions {
        for k in 0..3 { c[k] += p[k]; }
    }
    for k in 0..3 { c[k] /= n; }
    c
}

fn sub3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
