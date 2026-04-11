//! Local environment detection, fingerprinting, and matching.
//!
//! Three-step matching: graph hash → fingerprint → Kabsch + greedy permutation.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// 3D vector operations.
fn norm(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn sub(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// SVD of a 3x3 matrix using Jacobi iterations.
/// Returns (U, S, Vt) where M = U * diag(S) * Vt.
fn svd_3x3(m: &[[f64; 3]; 3]) -> ([[f64; 3]; 3], [f64; 3], [[f64; 3]; 3]) {
    // Compute M^T M
    let mut mtm = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                mtm[i][j] += m[k][i] * m[k][j];
            }
        }
    }

    // Eigendecompose M^T M via Jacobi rotations
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut a = mtm;

    for _ in 0..50 {
        // Find largest off-diagonal
        let pairs = [(0, 1), (0, 2), (1, 2)];
        for &(p, q) in &pairs {
            if a[p][q].abs() < 1e-15 {
                continue;
            }
            let tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = if tau >= 0.0 {
                1.0 / (tau + (1.0 + tau * tau).sqrt())
            } else {
                -1.0 / (-tau + (1.0 + tau * tau).sqrt())
            };
            let c = 1.0 / (1.0 + t * t).sqrt();
            let s = t * c;

            // Apply Jacobi rotation to a
            let mut new_a = a;
            new_a[p][p] = c * c * a[p][p] - 2.0 * s * c * a[p][q] + s * s * a[q][q];
            new_a[q][q] = s * s * a[p][p] + 2.0 * s * c * a[p][q] + c * c * a[q][q];
            new_a[p][q] = 0.0;
            new_a[q][p] = 0.0;
            for r in 0..3 {
                if r != p && r != q {
                    new_a[p][r] = c * a[p][r] - s * a[q][r];
                    new_a[r][p] = new_a[p][r];
                    new_a[q][r] = s * a[p][r] + c * a[q][r];
                    new_a[r][q] = new_a[q][r];
                }
            }
            a = new_a;

            // Update V
            let mut new_v = v;
            for r in 0..3 {
                new_v[r][p] = c * v[r][p] - s * v[r][q];
                new_v[r][q] = s * v[r][p] + c * v[r][q];
            }
            v = new_v;
        }
    }

    // Singular values
    let s = [a[0][0].sqrt().max(0.0), a[1][1].sqrt().max(0.0), a[2][2].sqrt().max(0.0)];

    // U = M V S^-1
    let mut u = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut val = 0.0;
            for k in 0..3 {
                val += m[i][k] * v[k][j];
            }
            u[i][j] = if s[j] > 1e-12 { val / s[j] } else { 0.0 };
        }
    }

    // Vt = V^T
    let mut vt = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            vt[i][j] = v[j][i];
        }
    }

    (u, s, vt)
}

/// Kabsch algorithm: find optimal rotation minimizing RMSD between X and Y.
/// Returns (rotation_3x3, rmsd).
pub fn kabsch(x: &[[f64; 3]], y: &[[f64; 3]]) -> ([[f64; 3]; 3], f64) {
    let n = x.len();
    assert_eq!(n, y.len());

    // Covariance matrix H = X^T Y
    let mut h = [[0.0f64; 3]; 3];
    for k in 0..n {
        for i in 0..3 {
            for j in 0..3 {
                h[i][j] += x[k][i] * y[k][j];
            }
        }
    }

    let (u, _s, vt) = svd_3x3(&h);

    // det(V * U^T)
    let mut vut = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                vut[i][j] += vt[k][i] * u[j][k]; // V * U^T
            }
        }
    }
    let det = vut[0][0] * (vut[1][1] * vut[2][2] - vut[1][2] * vut[2][1])
            - vut[0][1] * (vut[1][0] * vut[2][2] - vut[1][2] * vut[2][0])
            + vut[0][2] * (vut[1][0] * vut[2][1] - vut[1][1] * vut[2][0]);

    let sign = if det < 0.0 { -1.0 } else { 1.0 };

    // O = V * sign_matrix * U^T
    let mut o = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut val = 0.0;
            for k in 0..3 {
                let s_k = if k == 2 { sign } else { 1.0 };
                val += vt[k][i] * s_k * u[j][k];
            }
            o[i][j] = val;
        }
    }

    // RMSD
    let mut sum_sq = 0.0;
    for k in 0..n {
        for i in 0..3 {
            let mut rotated = 0.0;
            for j in 0..3 {
                rotated += o[i][j] * x[k][j];
            }
            let d = rotated - y[k][i];
            sum_sq += d * d;
        }
    }
    let rmsd = (sum_sq / n as f64).sqrt();

    (o, rmsd)
}

/// Fast fingerprint for local environment equivalence.
#[derive(Debug, Clone)]
pub struct Fingerprint {
    /// Sorted center-to-neighbor distances.
    pub r_0j: Vec<f64>,
    /// Sorted inter-neighbor distances.
    pub r_ij: Vec<f64>,
    /// Minimum center-neighbor distance.
    pub r_min: f64,
}

impl Fingerprint {
    /// Build fingerprint from local positions (center at origin = index 0).
    pub fn new(positions: &[[f64; 3]]) -> Self {
        let n = positions.len();

        let mut r_0j: Vec<f64> = if n > 1 {
            positions[1..].iter().map(|p| norm(p)).collect()
        } else {
            vec![]
        };
        r_0j.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mut r_ij = Vec::new();
        if n > 2 {
            for i in 1..n {
                for j in (i + 1)..n {
                    r_ij.push(norm(&sub(&positions[i], &positions[j])));
                }
            }
            r_ij.sort_by(|a, b| a.partial_cmp(b).unwrap());
        }

        let r_min = if r_0j.is_empty() { 1.0 } else { r_0j[0] };

        Fingerprint { r_0j, r_ij, r_min }
    }

    /// Fast equivalence check within tolerance delta * sqrt(2).
    pub fn equiv(&self, other: &Fingerprint, delta: f64) -> bool {
        let tol = delta * std::f64::consts::SQRT_2;

        if self.r_0j.len() != other.r_0j.len() {
            return false;
        }
        for (a, b) in self.r_0j.iter().zip(&other.r_0j) {
            if (a - b).abs() > tol {
                return false;
            }
        }
        if self.r_ij.len() != other.r_ij.len() {
            return false;
        }
        for (a, b) in self.r_ij.iter().zip(&other.r_ij) {
            if (a - b).abs() > tol {
                return false;
            }
        }
        true
    }
}

/// Local environment around a central atom.
#[derive(Debug, Clone)]
pub struct Geometry {
    /// Positions in local frame (center at origin).
    pub positions: Vec<[f64; 3]>,
    /// Colour encoding: 2 * type_id + is_frozen.
    pub colours: Vec<i32>,
    /// Global atom indices.
    pub indices: Vec<usize>,
    /// Cached fingerprint.
    fingerprint: Option<Fingerprint>,
    /// Cached graph hash.
    graph_hash: Option<u64>,
}

impl Geometry {
    pub fn new(positions: Vec<[f64; 3]>, colours: Vec<i32>, indices: Vec<usize>) -> Self {
        Geometry {
            positions,
            colours,
            indices,
            fingerprint: None,
            graph_hash: None,
        }
    }

    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }

    pub fn center_index(&self) -> usize {
        self.indices[0]
    }

    /// Compute and cache fingerprint.
    pub fn fingerprint(&mut self) -> &Fingerprint {
        if self.fingerprint.is_none() {
            self.fingerprint = Some(Fingerprint::new(&self.positions));
        }
        self.fingerprint.as_ref().unwrap()
    }

    /// Compute graph hash using sorted adjacency invariants.
    pub fn graph_hash(&mut self, r_edge: f64) -> u64 {
        if let Some(h) = self.graph_hash {
            return h;
        }

        let n = self.n_atoms();
        let mut invariants: Vec<(i32, usize, Vec<i32>)> = Vec::with_capacity(n);

        for i in 0..n {
            let mut neighbors = Vec::new();
            for j in 0..n {
                if i == j { continue; }
                let d = norm(&sub(&self.positions[i], &self.positions[j]));
                if d < r_edge {
                    neighbors.push(self.colours[j]);
                }
            }
            neighbors.sort();
            invariants.push((self.colours[i], neighbors.len(), neighbors));
        }
        invariants.sort();

        // Hash
        let mut hasher = DefaultHasher::new();
        self.colours[0].hash(&mut hasher);
        for (c, deg, neigh) in &invariants {
            c.hash(&mut hasher);
            deg.hash(&mut hasher);
            for n in neigh {
                n.hash(&mut hasher);
            }
        }
        let h = hasher.finish();
        self.graph_hash = Some(h);
        h
    }

    /// Greedy permutation matching onto reference.
    /// Returns Some((rotation, rmsd, permutation)) or None.
    pub fn permute_onto(&self, reference: &Geometry, delta: f64)
        -> Option<([[f64; 3]; 3], f64, Vec<usize>)>
    {
        let n = self.n_atoms();
        if n != reference.n_atoms() {
            return None;
        }

        let mut best: Option<([[f64; 3]; 3], f64, Vec<usize>)> = None;
        let mut best_rmsd = delta;

        fn greedy(
            slf: &Geometry, reference: &Geometry, n: usize,
            assigned: &mut Vec<bool>, perm: &mut Vec<usize>, depth: usize,
            best: &mut Option<([[f64; 3]; 3], f64, Vec<usize>)>,
            best_rmsd: &mut f64, delta: f64,
        ) {
            if depth == n {
                let x: Vec<[f64; 3]> = perm.iter().map(|&i| slf.positions[i]).collect();
                let y: Vec<[f64; 3]> = reference.positions.iter().copied().collect();
                let (o, rmsd) = kabsch(&x, &y);
                if rmsd < *best_rmsd {
                    *best_rmsd = rmsd;
                    *best = Some((o, rmsd, perm.clone()));
                }
                return;
            }

            let ref_col = reference.colours[depth];
            let tol = delta * std::f64::consts::SQRT_2;

            for self_idx in 0..n {
                if assigned[self_idx] { continue; }
                if slf.colours[self_idx] != ref_col { continue; }

                // Distance pruning
                let mut ok = true;
                for k in 0..depth {
                    let d_self = norm(&sub(&slf.positions[perm[k]], &slf.positions[self_idx]));
                    let d_ref = norm(&sub(&reference.positions[k], &reference.positions[depth]));
                    if (d_self - d_ref).abs() > tol {
                        ok = false;
                        break;
                    }
                }
                if !ok { continue; }

                perm.push(self_idx);
                assigned[self_idx] = true;
                greedy(slf, reference, n, assigned, perm, depth + 1, best, best_rmsd, delta);
                assigned[self_idx] = false;
                perm.pop();

                // Early exit
                if best.is_some() { return; }
            }
        }

        let mut assigned = vec![false; n];
        let mut perm = Vec::with_capacity(n);
        greedy(self, reference, n, &mut assigned, &mut perm, 0,
               &mut best, &mut best_rmsd, delta);

        best
    }
}

/// Build local environment from atomic data.
pub fn build_geometry(
    center_idx: usize,
    positions: &[[f64; 3]],
    types: &[i32],
    frozen: &[bool],
    neighbor_indices: &[usize],
    r_env: f64,
) -> Geometry {
    let center_pos = positions[center_idx];

    let mut local_pos = vec![[0.0, 0.0, 0.0]]; // center at origin
    let mut local_col = vec![2 * types[center_idx] + frozen[center_idx] as i32];
    let mut local_idx = vec![center_idx];

    for &j in neighbor_indices {
        if j == center_idx { continue; }
        let dr = sub(&positions[j], &center_pos);
        let dist = norm(&dr);
        if dist < r_env {
            local_pos.push(dr);
            local_col.push(2 * types[j] + frozen[j] as i32);
            local_idx.push(j);
        }
    }

    // Center (shift centroid to origin)
    let n = local_pos.len();
    if n > 0 {
        let mut centroid = [0.0; 3];
        for p in &local_pos {
            for k in 0..3 { centroid[k] += p[k]; }
        }
        for k in 0..3 { centroid[k] /= n as f64; }
        for p in &mut local_pos {
            for k in 0..3 { p[k] -= centroid[k]; }
        }
    }

    Geometry::new(local_pos, local_col, local_idx)
}
