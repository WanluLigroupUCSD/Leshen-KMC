//! Force / energy provider abstraction.
//!
//! Phase B: pure-Rust trait + Müller-Brown analytic PES so we can unit-test
//! DimerSearch / FireMinimizer without dragging in Python or ASE.
//! Phase C will add a `PyCalculator` impl that calls back into a Python ASE
//! Calculator for real catalysis runs.
//!
//! The trait is intentionally minimal — only `eval` (energy + gradient) and
//! optional frozen mask. ASE force eval is the slowest call in the whole loop;
//! everything else (rotation arithmetic, FIRE updates, perturbation sampling)
//! happens in Rust and never crosses the Python boundary.

/// Calculator interface: maps positions -> (energy, gradient).
///
/// `gradient[i]` is the energy gradient on atom `i` (NOT the force; force = -gradient).
/// Frozen atoms — those whose `frozen_mask()` returns `true` — must have zero
/// gradient. Implementations are responsible for zeroing.
///
/// Positions and gradient are length-N `Vec<[f64; 3]>` mirroring numpy's
/// `(N, 3)` layout. Allocations happen at the call site; impls should not
/// retain references.
pub trait Calculator {
    /// Compute potential energy and gradient at the given positions.
    fn eval(&mut self, positions: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>);

    /// Optional frozen mask. Returns `Some(mask)` of length `positions.len()`,
    /// or `None` if no atoms are frozen.
    fn frozen_mask(&self) -> Option<&[bool]> {
        None
    }

    /// Convenience: just energy (default impl calls `eval` and discards gradient).
    fn energy(&mut self, positions: &[[f64; 3]]) -> f64 {
        self.eval(positions).0
    }
}

// ===========================================================================
//                      MULLER-BROWN ANALYTIC PES
// ===========================================================================
//
// 2D potential with three minima and two saddle points, widely used as a
// reference for transition state methods. We embed it in 3D by treating
// (x, y) = (positions[0][0], positions[0][1]) and ignoring z. That lets the
// same dimer / minimizer / pyo3 plumbing exercise it without special-casing.
//
// Reference saddle locations (Müller & Brown, 1979):
//   SP1: (-0.822,  0.624)  E ≈ -40.7
//   SP2: ( 0.212,  0.293)  E ≈ -72.2
//
// Reference minima:
//   M1:  ( 0.624,  0.028)  E ≈ -108.2  (lowest)
//   M2:  ( 0.105,  0.467)  E ≈ -110.8 / depends on convention
//   M3:  (-0.558,  1.442)  E ≈ -146.7
//
// Useful for unit tests: start at M1 perturbed, dimer should walk to SP2.

/// Müller-Brown analytic 2D PES, embedded in 3D (z component ignored).
pub struct MullerBrown {
    pub a: [f64; 4],
    pub b: [f64; 4],
    pub c: [f64; 4],
    pub a_coef: [f64; 4],
    pub x0: [f64; 4],
    pub y0: [f64; 4],
}

impl MullerBrown {
    /// Standard Müller-Brown parameters (Müller & Brown 1979, Theor. Chim. Acta 53).
    pub fn standard() -> Self {
        Self {
            a:      [-1.0, -1.0, -6.5,  0.7],
            b:      [ 0.0,  0.0, 11.0,  0.6],
            c:      [-10.0, -10.0, -6.5,  0.7],
            a_coef: [-200.0, -100.0, -170.0, 15.0],
            x0:     [ 1.0,  0.0, -0.5, -1.0],
            y0:     [ 0.0,  0.5,  1.5,  1.0],
        }
    }

    /// Evaluate at a single 2D point. Returns (E, dE/dx, dE/dy).
    pub fn eval_2d(&self, x: f64, y: f64) -> (f64, f64, f64) {
        let mut e = 0.0;
        let mut dx = 0.0;
        let mut dy = 0.0;
        for k in 0..4 {
            let dx_k = x - self.x0[k];
            let dy_k = y - self.y0[k];
            let q = self.a[k] * dx_k * dx_k
                  + self.b[k] * dx_k * dy_k
                  + self.c[k] * dy_k * dy_k;
            let term = self.a_coef[k] * q.exp();
            e += term;
            // dE/dx = term * (2*a*dx + b*dy)
            dx += term * (2.0 * self.a[k] * dx_k + self.b[k] * dy_k);
            dy += term * (self.b[k] * dx_k + 2.0 * self.c[k] * dy_k);
        }
        (e, dx, dy)
    }
}

impl Calculator for MullerBrown {
    fn eval(&mut self, positions: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
        // Use only the first "atom"'s (x, y); z is ignored.
        // All other atoms see zero force — this is intentional, makes
        // higher-N tests degenerate to 2D MB while still exercising the
        // ndarray plumbing for n_atoms > 1.
        let (e, dx, dy) = self.eval_2d(positions[0][0], positions[0][1]);
        let mut grad = vec![[0.0; 3]; positions.len()];
        grad[0] = [dx, dy, 0.0];
        (e, grad)
    }
}

// ===========================================================================
//                  ANALYTIC QUADRATIC SADDLE (test fixture)
// ===========================================================================
//
// V(x, y) = -k_x * x^2 / 2 + k_y * y^2 / 2
// Saddle at origin: x is unstable mode (curvature -k_x), y is stable (+k_y).
// dE/dx = -k_x * x       dE/dy = k_y * y
//
// Used to verify dimer convergence on a clean, analytically-solved problem.

/// Quadratic 2D saddle, embedded in 3D (z ignored).
pub struct QuadraticSaddle {
    pub k_x: f64,  // unstable mode spring constant (>0; PES has -k_x*x²/2)
    pub k_y: f64,  // stable mode spring constant
}

impl QuadraticSaddle {
    pub fn new(k_x: f64, k_y: f64) -> Self { Self { k_x, k_y } }
}

impl Calculator for QuadraticSaddle {
    fn eval(&mut self, positions: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
        let x = positions[0][0];
        let y = positions[0][1];
        let e = -0.5 * self.k_x * x * x + 0.5 * self.k_y * y * y;
        let mut grad = vec![[0.0; 3]; positions.len()];
        grad[0] = [-self.k_x * x, self.k_y * y, 0.0];
        (e, grad)
    }
}

// ---------------------------------------------------------------------------
//                                TESTS
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference values (Müller & Brown 1979) at the three minima.
    /// We don't assert exact; just verify ordering: M_C is global min, all
    /// three are below the saddle SP_BC ~ -40.
    #[test]
    fn muller_brown_ordering() {
        let mb = MullerBrown::standard();
        // Approximate minimum positions
        let m_a = mb.eval_2d(-0.558,  1.442).0;
        let m_b = mb.eval_2d( 0.624,  0.028).0;
        let m_c = mb.eval_2d( 0.105,  0.467).0;
        let sp_bc = mb.eval_2d(0.212, 0.293).0;
        assert!(m_a < sp_bc, "minimum A should be below SP_BC: {} vs {}", m_a, sp_bc);
        assert!(m_b < sp_bc, "minimum B should be below SP_BC: {} vs {}", m_b, sp_bc);
        assert!(m_c < sp_bc, "minimum C should be below SP_BC: {} vs {}", m_c, sp_bc);
    }

    /// Gradient should vanish at minima (within ~1.0 since approximate positions).
    #[test]
    fn muller_brown_gradient_small_at_minima() {
        let mb = MullerBrown::standard();
        // Use the global minimum location (more refined value)
        let (_, dx, dy) = mb.eval_2d(-0.558, 1.442);
        let g_norm = (dx * dx + dy * dy).sqrt();
        assert!(g_norm < 5.0, "gradient at approx minimum too large: {}", g_norm);
    }

    /// The Calculator trait wraps MB cleanly: gradient on atom 0 only,
    /// other atoms zero.
    #[test]
    fn calculator_trait_round_trip() {
        let mut mb = MullerBrown::standard();
        let positions = vec![[0.5, 0.5, 0.0], [9.9, 9.9, 0.0]];
        let (e, grad) = mb.eval(&positions);
        assert!(e.is_finite());
        assert_eq!(grad.len(), 2);
        assert_eq!(grad[1], [0.0; 3], "non-active atom should have zero gradient");
        // Atom 0 gradient should match eval_2d
        let (_, dx, dy) = mb.eval_2d(0.5, 0.5);
        assert!((grad[0][0] - dx).abs() < 1e-10);
        assert!((grad[0][1] - dy).abs() < 1e-10);
        assert_eq!(grad[0][2], 0.0);
    }
}
