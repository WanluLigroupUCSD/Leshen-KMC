/// Polarization curve calculation from constant-potential DFT data.
///
/// Workflow:
///   1. Load DFT energies at multiple potentials from JSON
///   2. Build energy landscape with Ea(U) via cubic spline or quadratic fit
///   3. Create microkinetic model with potential-dependent rates
///   4. Sweep potential → steady state → TOF → current density
///
/// Three fitting modes for Ea(U):
///   (a) Cubic spline interpolation through tabulated Ea(U) points  [default]
///   (b) Quadratic regression:  Ω(U) = aU² + bU + c  fitted to DFT data
///   (c) Direct quadratic coefficients (user provides a, b, c per state)

use std::collections::HashMap;
use std::fs;
use std::io::Write;
use crate::units::*;
use crate::rates::{tst_rate, hertz_knudsen};
use crate::microkinetic::MicroKineticModel;

// ═══════════════════════════════════════════════════════════════════
//  Natural cubic spline interpolation
// ═══════════════════════════════════════════════════════════════════

/// Natural cubic spline with extrapolation via boundary segments.
#[derive(Clone)]
struct CubicSpline {
    x: Vec<f64>,
    y: Vec<f64>,
    /// Second derivatives at each knot (M₀ = Mₙ = 0 for natural BC)
    m: Vec<f64>,
}

impl CubicSpline {
    fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        let n = x.len();
        assert!(n >= 2, "Need at least 2 data points for interpolation");
        assert_eq!(x.len(), y.len(), "x and y must have the same length");

        // Linear case: second derivatives are all zero
        if n == 2 {
            return Self { x, y, m: vec![0.0; 2] };
        }

        // Compute intervals h_i = x_{i+1} - x_i
        let h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();

        // Tridiagonal system for M[1]..M[n-2]  (n_inner = n-2 unknowns)
        let ni = n - 2;
        let mut diag = vec![0.0; ni];
        let mut sub = vec![0.0; ni];
        let mut sup = vec![0.0; ni];
        let mut rhs = vec![0.0; ni];

        for i in 0..ni {
            let j = i + 1; // actual knot index
            diag[i] = 2.0 * (h[j - 1] + h[j]);
            rhs[i] = 6.0 * ((y[j + 1] - y[j]) / h[j] - (y[j] - y[j - 1]) / h[j - 1]);
            if i > 0 {
                sub[i] = h[j - 1];
            }
            if i < ni - 1 {
                sup[i] = h[j];
            }
        }

        // Thomas algorithm (tridiagonal solve)
        for i in 1..ni {
            let factor = sub[i] / diag[i - 1];
            diag[i] -= factor * sup[i - 1];
            rhs[i] -= factor * rhs[i - 1];
        }

        // Back substitution → m[1..n-1]
        let mut m = vec![0.0; n]; // m[0] = m[n-1] = 0 (natural BC)
        if ni > 0 {
            m[ni] = rhs[ni - 1] / diag[ni - 1];
            for i in (0..ni - 1).rev() {
                m[i + 1] = (rhs[i] - sup[i] * m[i + 2]) / diag[i];
            }
        }

        Self { x, y, m }
    }

    fn evaluate(&self, xval: f64) -> f64 {
        let n = self.x.len();

        // Linear fallback for 2 points
        if n == 2 {
            let h = self.x[1] - self.x[0];
            return self.y[0] + (xval - self.x[0]) / h * (self.y[1] - self.y[0]);
        }

        // Find interval index (boundary segments used for extrapolation)
        let i = if xval <= self.x[0] {
            0
        } else if xval >= self.x[n - 1] {
            n - 2
        } else {
            let mut lo = 0usize;
            let mut hi = n - 1;
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                if self.x[mid] <= xval {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        };

        let h = self.x[i + 1] - self.x[i];
        let a = self.x[i + 1] - xval; // distance to right knot
        let b = xval - self.x[i]; // distance to left knot

        // S_i(x) = M_i*(x_{i+1}-x)^3/(6h) + M_{i+1}*(x-x_i)^3/(6h)
        //        + (y_i/h - M_i*h/6)*(x_{i+1}-x) + (y_{i+1}/h - M_{i+1}*h/6)*(x-x_i)
        self.m[i] * a * a * a / (6.0 * h)
            + self.m[i + 1] * b * b * b / (6.0 * h)
            + (self.y[i] / h - self.m[i] * h / 6.0) * a
            + (self.y[i + 1] / h - self.m[i + 1] * h / 6.0) * b
    }
}

// ═══════════════════════════════════════════════════════════════════
//  Grand canonical quadratic energy: Ω(U) = aU² + bU + c
// ═══════════════════════════════════════════════════════════════════

/// Grand canonical energy as a quadratic function of applied potential.
///
/// Ω(U) = a·U² + b·U + c
///
/// Physical origin: in constant-potential DFT, the grand canonical energy
/// varies parabolically with potential due to the roughly constant
/// interfacial capacitance:  a ≈ -C/2,  b relates to the equilibrium
/// charge,  c = Ω(0).
#[derive(Clone, Debug)]
pub struct QuadraticEnergy {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl QuadraticEnergy {
    pub fn new(a: f64, b: f64, c: f64) -> Self {
        Self { a, b, c }
    }

    /// Evaluate Ω(U).
    pub fn evaluate(&self, u: f64) -> f64 {
        self.a * u * u + self.b * u + self.c
    }

    /// Fit Ω(U) = aU² + bU + c to tabulated (U, E) data via least squares.
    ///
    /// Returns (QuadraticEnergy, r_squared).
    pub fn fit(u_data: &[f64], e_data: &[f64]) -> (Self, f64) {
        assert!(u_data.len() >= 3, "Need >= 3 points for quadratic fit");
        assert_eq!(u_data.len(), e_data.len());

        let n = u_data.len() as f64;

        // Normal equations:  X^T X β = X^T y
        //   X = [[u², u, 1], ...],  y = [e, ...],  β = [a, b, c]
        let su4: f64 = u_data.iter().map(|u| u.powi(4)).sum();
        let su3: f64 = u_data.iter().map(|u| u.powi(3)).sum();
        let su2: f64 = u_data.iter().map(|u| u.powi(2)).sum();
        let su1: f64 = u_data.iter().sum();

        let su2e: f64 = u_data.iter().zip(e_data).map(|(u, e)| u * u * e).sum();
        let su1e: f64 = u_data.iter().zip(e_data).map(|(u, e)| u * e).sum();
        let se: f64 = e_data.iter().sum();

        let (a, b, c) = solve_3x3(
            [su4, su3, su2, su3, su2, su1, su2, su1, n],
            [su2e, su1e, se],
        );

        // R² = 1 - SS_res / SS_tot
        let e_mean = se / n;
        let ss_tot: f64 = e_data.iter().map(|e| (e - e_mean).powi(2)).sum();
        let ss_res: f64 = u_data
            .iter()
            .zip(e_data)
            .map(|(u, e)| {
                let predicted = a * u * u + b * u + c;
                (e - predicted).powi(2)
            })
            .sum();
        let r2 = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { 1.0 };

        (Self { a, b, c }, r2)
    }
}

/// Solve 3×3 linear system Ax = b via Gaussian elimination with partial pivoting.
fn solve_3x3(mat: [f64; 9], rhs: [f64; 3]) -> (f64, f64, f64) {
    let mut a = [
        [mat[0], mat[1], mat[2]],
        [mat[3], mat[4], mat[5]],
        [mat[6], mat[7], mat[8]],
    ];
    let mut b = [rhs[0], rhs[1], rhs[2]];

    for col in 0..3 {
        // Partial pivoting
        let mut max_row = col;
        let mut max_val = a[col][col].abs();
        for row in col + 1..3 {
            if a[row][col].abs() > max_val {
                max_val = a[row][col].abs();
                max_row = row;
            }
        }
        a.swap(col, max_row);
        b.swap(col, max_row);

        // Eliminate below
        for row in col + 1..3 {
            let factor = a[row][col] / a[col][col];
            for j in col..3 {
                a[row][j] -= factor * a[col][j];
            }
            b[row] -= factor * b[col];
        }
    }

    let c = b[2] / a[2][2];
    let bv = (b[1] - a[1][2] * c) / a[1][1];
    let av = (b[0] - a[0][1] * bv - a[0][2] * c) / a[0][0];
    (av, bv, c)
}

/// Fitting mode for constructing Ea(U).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FitMode {
    /// Natural cubic spline interpolation (default).
    Spline,
    /// Quadratic regression: Ω(U) = aU² + bU + c.
    Quadratic,
}

// ═══════════════════════════════════════════════════════════════════
//  InterpolatedBarrier: Ea(U) from tabulated or fitted data
// ═══════════════════════════════════════════════════════════════════

#[derive(Clone)]
enum BarrierKind {
    Spline(CubicSpline),
    Quadratic(QuadraticEnergy),
}

/// Activation barrier Ea(U) supporting multiple representations:
///  - Cubic spline through tabulated data
///  - Quadratic fit (aU² + bU + c) to tabulated data
///  - Direct quadratic coefficients
#[derive(Clone)]
pub struct InterpolatedBarrier {
    kind: BarrierKind,
}

impl InterpolatedBarrier {
    /// Create from tabulated data using cubic spline.
    pub fn new(u_data: Vec<f64>, ea_data: Vec<f64>) -> Self {
        Self {
            kind: BarrierKind::Spline(CubicSpline::new(u_data, ea_data)),
        }
    }

    /// Create from tabulated data using quadratic fit. Returns (barrier, R²).
    pub fn new_quadratic_fit(u_data: &[f64], ea_data: &[f64]) -> (Self, f64) {
        let (qe, r2) = QuadraticEnergy::fit(u_data, ea_data);
        (Self { kind: BarrierKind::Quadratic(qe) }, r2)
    }

    /// Create directly from quadratic coefficients: Ea(U) = aU² + bU + c.
    pub fn from_quadratic(a: f64, b: f64, c: f64) -> Self {
        Self {
            kind: BarrierKind::Quadratic(QuadraticEnergy::new(a, b, c)),
        }
    }

    /// Return Ea [eV] at potential U [V]. Clamped to >= 0.
    pub fn evaluate(&self, u: f64) -> f64 {
        match &self.kind {
            BarrierKind::Spline(s) => s.evaluate(u).max(0.0),
            BarrierKind::Quadratic(q) => q.evaluate(u).max(0.0),
        }
    }

    /// Create a TST rate function for use in MicroKineticModel.
    ///
    /// k(U, T) = (kB*T/h) * exp(-Ea(U) / (kB*T))
    pub fn make_rate_fn(&self) -> Box<dyn Fn(&HashMap<String, f64>) -> f64> {
        let barrier = self.clone();
        Box::new(move |params: &HashMap<String, f64>| {
            let ea = barrier.evaluate(params["U"]);
            tst_rate(ea, params["T"])
        })
    }
}

// ═══════════════════════════════════════════════════════════════════
//  EnergyLandscape: compute Ea(U) from DFT energy surfaces
// ═══════════════════════════════════════════════════════════════════

/// Energy landscape from constant-potential DFT calculations.
///
/// Given E_IS(U), E_TS(U), E_FS(U) at several potentials, computes:
///   Ea_fwd(U) = max(E_TS(U) - E_IS(U), 0)
///   Ea_rev(U) = max(E_TS(U) - E_FS(U), 0)
pub struct EnergyLandscape {
    pub potentials: Vec<f64>,
    barriers_fwd: HashMap<(String, String), InterpolatedBarrier>,
    barriers_rev: HashMap<(String, String), InterpolatedBarrier>,
}

impl EnergyLandscape {
    pub fn new(
        potentials: Vec<f64>,
        state_energies: &HashMap<String, Vec<f64>>,
        ts_energies: &HashMap<(String, String), Vec<f64>>,
    ) -> Self {
        let n = potentials.len();
        let mut barriers_fwd = HashMap::new();
        let mut barriers_rev = HashMap::new();

        for ((is_name, fs_name), e_ts) in ts_energies {
            let e_is = state_energies
                .get(is_name)
                .unwrap_or_else(|| panic!("Missing state energy for '{}'", is_name));
            let e_fs = state_energies
                .get(fs_name)
                .unwrap_or_else(|| panic!("Missing state energy for '{}'", fs_name));

            assert_eq!(e_ts.len(), n, "TS energy length mismatch for {}->{}",
                       is_name, fs_name);
            assert_eq!(e_is.len(), n, "IS energy length mismatch for {}", is_name);
            assert_eq!(e_fs.len(), n, "FS energy length mismatch for {}", fs_name);

            let ea_fwd: Vec<f64> = (0..n)
                .map(|i| (e_ts[i] - e_is[i]).max(0.0))
                .collect();
            let ea_rev: Vec<f64> = (0..n)
                .map(|i| (e_ts[i] - e_fs[i]).max(0.0))
                .collect();

            let key = (is_name.clone(), fs_name.clone());
            barriers_fwd.insert(
                key.clone(),
                InterpolatedBarrier::new(potentials.clone(), ea_fwd),
            );
            barriers_rev.insert(
                key,
                InterpolatedBarrier::new(potentials.clone(), ea_rev),
            );
        }

        Self {
            potentials,
            barriers_fwd,
            barriers_rev,
        }
    }

    /// Create from tabulated data using quadratic regression on each state.
    ///
    /// Fits Ω(U) = aU² + bU + c to each intermediate state and TS,
    /// then computes barriers as differences of the fitted parabolas:
    ///   Ea_fwd(U) = Ω_TS(U) - Ω_IS(U)   (a quadratic itself)
    ///
    /// Prints fit quality (R²) for each state.
    pub fn new_quadratic_fit(
        potentials: Vec<f64>,
        state_energies: &HashMap<String, Vec<f64>>,
        ts_energies: &HashMap<(String, String), Vec<f64>>,
    ) -> Self {
        let mut barriers_fwd = HashMap::new();
        let mut barriers_rev = HashMap::new();

        // Fit each state energy to a quadratic
        let mut state_fits: HashMap<String, QuadraticEnergy> = HashMap::new();
        eprintln!("\nQuadratic fits for state energies:");
        eprintln!("  {:>12}  {:>10}  {:>10}  {:>10}  {:>6}",
                  "State", "a", "b", "c", "R²");
        eprintln!("  {}  {}  {}  {}  {}",
                  "-".repeat(12), "-".repeat(10), "-".repeat(10),
                  "-".repeat(10), "-".repeat(6));

        for (name, energies) in state_energies {
            let (qe, r2) = QuadraticEnergy::fit(&potentials, energies);
            eprintln!("  {:>12}  {:>10.5}  {:>10.5}  {:>10.5}  {:>6.4}",
                      name, qe.a, qe.b, qe.c, r2);
            state_fits.insert(name.clone(), qe);
        }

        // Fit each TS energy and compute barrier quadratics
        eprintln!("\nQuadratic fits for TS energies:");
        eprintln!("  {:>20}  {:>10}  {:>10}  {:>10}  {:>6}",
                  "Transition", "a", "b", "c", "R²");
        eprintln!("  {}  {}  {}  {}  {}",
                  "-".repeat(20), "-".repeat(10), "-".repeat(10),
                  "-".repeat(10), "-".repeat(6));

        for ((is_name, fs_name), e_ts) in ts_energies {
            let (ts_fit, r2) = QuadraticEnergy::fit(&potentials, e_ts);
            eprintln!("  {:>20}  {:>10.5}  {:>10.5}  {:>10.5}  {:>6.4}",
                      format!("{}->{}",is_name, fs_name), ts_fit.a, ts_fit.b, ts_fit.c, r2);

            let is_fit = state_fits.get(is_name)
                .unwrap_or_else(|| panic!("Missing state fit for '{}'", is_name));
            let fs_fit = state_fits.get(fs_name)
                .unwrap_or_else(|| panic!("Missing state fit for '{}'", fs_name));

            // Barrier = Ω_TS - Ω_IS  (difference of quadratics is a quadratic)
            let key = (is_name.clone(), fs_name.clone());
            barriers_fwd.insert(
                key.clone(),
                InterpolatedBarrier::from_quadratic(
                    ts_fit.a - is_fit.a,
                    ts_fit.b - is_fit.b,
                    ts_fit.c - is_fit.c,
                ),
            );
            barriers_rev.insert(
                key,
                InterpolatedBarrier::from_quadratic(
                    ts_fit.a - fs_fit.a,
                    ts_fit.b - fs_fit.b,
                    ts_fit.c - fs_fit.c,
                ),
            );
        }

        Self { potentials, barriers_fwd, barriers_rev }
    }

    /// Create directly from user-provided quadratic coefficients.
    ///
    /// Each state/TS is described by Ω(U) = aU² + bU + c.
    /// Barrier coefficients are computed as differences.
    pub fn from_quadratic_coefficients(
        state_coefficients: &HashMap<String, [f64; 3]>,
        ts_coefficients: &HashMap<(String, String), [f64; 3]>,
    ) -> Self {
        let mut barriers_fwd = HashMap::new();
        let mut barriers_rev = HashMap::new();

        for ((is_name, fs_name), ts_abc) in ts_coefficients {
            let is_abc = state_coefficients.get(is_name)
                .unwrap_or_else(|| panic!("Missing coefficients for '{}'", is_name));
            let fs_abc = state_coefficients.get(fs_name)
                .unwrap_or_else(|| panic!("Missing coefficients for '{}'", fs_name));

            let key = (is_name.clone(), fs_name.clone());
            barriers_fwd.insert(
                key.clone(),
                InterpolatedBarrier::from_quadratic(
                    ts_abc[0] - is_abc[0],
                    ts_abc[1] - is_abc[1],
                    ts_abc[2] - is_abc[2],
                ),
            );
            barriers_rev.insert(
                key,
                InterpolatedBarrier::from_quadratic(
                    ts_abc[0] - fs_abc[0],
                    ts_abc[1] - fs_abc[1],
                    ts_abc[2] - fs_abc[2],
                ),
            );
        }

        Self {
            potentials: Vec::new(),
            barriers_fwd,
            barriers_rev,
        }
    }

    /// Check if a transition IS→FS exists in the landscape.
    pub fn has_transition(&self, is: &str, fs: &str) -> bool {
        self.barriers_fwd
            .contains_key(&(is.to_string(), fs.to_string()))
    }

    /// Get forward rate function for IS → FS.
    pub fn get_fwd_rate_fn(
        &self,
        is: &str,
        fs: &str,
    ) -> Option<Box<dyn Fn(&HashMap<String, f64>) -> f64>> {
        self.barriers_fwd
            .get(&(is.to_string(), fs.to_string()))
            .map(|b| b.make_rate_fn())
    }

    /// Get reverse rate function for FS → IS.
    pub fn get_rev_rate_fn(
        &self,
        is: &str,
        fs: &str,
    ) -> Option<Box<dyn Fn(&HashMap<String, f64>) -> f64>> {
        self.barriers_rev
            .get(&(is.to_string(), fs.to_string()))
            .map(|b| b.make_rate_fn())
    }

    /// Print barrier summary. Uses tabulated potentials if available,
    /// otherwise generates a default range [-2.0, 0.0] with 5 points.
    pub fn summary(&self) {
        let display_u: Vec<f64> = if self.potentials.is_empty() {
            vec![-2.0, -1.5, -1.0, -0.5, 0.0]
        } else {
            self.potentials.clone()
        };

        println!("Energy Landscape Summary");
        println!("{}", "=".repeat(70));
        for ((is, fs), fwd) in &self.barriers_fwd {
            let rev = &self.barriers_rev[&(is.clone(), fs.clone())];
            println!("\n  {} -> {}:", is, fs);
            println!(
                "    {:>8}  {:>12}  {:>12}",
                "U [V]", "Ea_fwd [eV]", "Ea_rev [eV]"
            );
            println!(
                "    {}  {}  {}",
                "-".repeat(8),
                "-".repeat(12),
                "-".repeat(12)
            );
            for &u in &display_u {
                println!(
                    "    {:>8.3}  {:>12.4}  {:>12.4}",
                    u,
                    fwd.evaluate(u),
                    rev.evaluate(u)
                );
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════
//  DFT data loading from JSON
// ═══════════════════════════════════════════════════════════════════

/// Parsed DFT data from JSON, supporting two formats.
pub enum DftData {
    /// Tabulated energies at discrete potentials.
    Tabulated {
        potentials: Vec<f64>,
        state_energies: HashMap<String, Vec<f64>>,
        ts_energies: HashMap<(String, String), Vec<f64>>,
    },
    /// Direct quadratic coefficients:  Ω(U) = aU² + bU + c.
    QuadraticCoefficients {
        state_coefficients: HashMap<String, [f64; 3]>,
        ts_coefficients: HashMap<(String, String), [f64; 3]>,
    },
}

/// Load DFT energy data from a JSON file.
///
/// Auto-detects two formats:
///
/// **Format 1 — Tabulated** (default):
/// ```json
/// {
///   "potentials": [-2.0, -1.5, -1.0, -0.5, 0.0],
///   "state_energies": { "N2*": [...], "NNH*": [...] },
///   "ts_energies": { "N2*->NNH*": [...] }
/// }
/// ```
///
/// **Format 2 — Quadratic coefficients**:
/// ```json
/// {
///   "mode": "quadratic_coefficients",
///   "state_coefficients": { "N2*": [a, b, c], ... },
///   "ts_coefficients": { "N2*->NNH*": [a, b, c], ... }
/// }
/// ```
/// where Ω(U) = a·U² + b·U + c  in eV.
pub fn load_dft_data(path: &str) -> DftData {
    let text =
        fs::read_to_string(path).unwrap_or_else(|e| panic!("Cannot read '{}': {}", path, e));
    let data: serde_json::Value =
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("Invalid JSON in '{}': {}", path, e));

    // Detect format
    let mode = data.get("mode").and_then(|v| v.as_str()).unwrap_or("tabulated");

    if mode == "quadratic_coefficients" {
        // Format 2: direct coefficients
        let mut state_coefficients = HashMap::new();
        for (key, vals) in data["state_coefficients"]
            .as_object()
            .expect("'state_coefficients' must be an object")
        {
            let abc: Vec<f64> = vals.as_array().unwrap().iter()
                .map(|v| v.as_f64().unwrap()).collect();
            assert_eq!(abc.len(), 3, "Coefficients for '{}' must have 3 elements [a, b, c]", key);
            state_coefficients.insert(key.clone(), [abc[0], abc[1], abc[2]]);
        }

        let mut ts_coefficients = HashMap::new();
        for (key, vals) in data["ts_coefficients"]
            .as_object()
            .expect("'ts_coefficients' must be an object")
        {
            let abc: Vec<f64> = vals.as_array().unwrap().iter()
                .map(|v| v.as_f64().unwrap()).collect();
            assert_eq!(abc.len(), 3, "Coefficients for '{}' must have 3 elements [a, b, c]", key);
            if let Some((is, fs)) = key.split_once("->") {
                ts_coefficients.insert(
                    (is.trim().to_string(), fs.trim().to_string()),
                    [abc[0], abc[1], abc[2]],
                );
            }
        }

        DftData::QuadraticCoefficients { state_coefficients, ts_coefficients }
    } else {
        // Format 1: tabulated
        let potentials: Vec<f64> = data["potentials"]
            .as_array()
            .expect("'potentials' must be an array")
            .iter()
            .map(|v| v.as_f64().expect("potential values must be numbers"))
            .collect();

        let mut state_energies = HashMap::new();
        for (key, vals) in data["state_energies"]
            .as_object()
            .expect("'state_energies' must be an object")
        {
            let energies: Vec<f64> = vals
                .as_array()
                .unwrap_or_else(|| panic!("state_energies['{}'] must be an array", key))
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect();
            state_energies.insert(key.clone(), energies);
        }

        let mut ts_energies: HashMap<(String, String), Vec<f64>> = HashMap::new();
        for (key, vals) in data["ts_energies"]
            .as_object()
            .expect("'ts_energies' must be an object")
        {
            let energies: Vec<f64> = vals
                .as_array()
                .unwrap_or_else(|| panic!("ts_energies['{}'] must be an array", key))
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect();
            if let Some((is, fs)) = key.split_once("->") {
                ts_energies.insert((is.trim().to_string(), fs.trim().to_string()), energies);
            } else {
                eprintln!("Warning: skipping TS key '{}' (expected 'IS->FS' format)", key);
            }
        }

        DftData::Tabulated { potentials, state_energies, ts_energies }
    }
}

// ═══════════════════════════════════════════════════════════════════
//  Current density conversion
// ═══════════════════════════════════════════════════════════════════

/// Convert TOF [s⁻¹ site⁻¹] to current density [mA/cm²].
///
/// j [A/m²] = n_e × e × TOF / A_site
/// j [mA/cm²] = j [A/m²] × 0.1
pub fn tof_to_current_density(tof: f64, n_electrons: f64, a_site: f64) -> f64 {
    n_electrons * E_CHARGE * tof / a_site * 0.1
}

// ═══════════════════════════════════════════════════════════════════
//  Polarization curve computation
// ═══════════════════════════════════════════════════════════════════

/// Results from a polarization curve computation.
pub struct PolarizationResult {
    pub u_values: Vec<f64>,
    /// TOF arrays per observable
    pub tof: HashMap<String, Vec<f64>>,
    /// Current density arrays per observable [mA/cm²]
    pub j_per_obs: HashMap<String, Vec<f64>>,
    /// Total current density [mA/cm²]
    pub j_total: Vec<f64>,
    /// Coverage at each potential
    pub coverages: Vec<HashMap<String, f64>>,
    pub temp: f64,
    pub a_site: f64,
}

impl PolarizationResult {
    /// Save results to a tab-separated file.
    pub fn save(&self, filename: &str, species: &[String]) {
        let mut f = fs::File::create(filename)
            .unwrap_or_else(|e| panic!("Cannot create '{}': {}", filename, e));

        // Header
        let obs_names: Vec<&String> = self.tof.keys().collect();
        let mut cols = vec!["U_V".to_string()];
        for obs in &obs_names {
            cols.push(format!("TOF_{}", obs));
            cols.push(format!("j_{}_mA_cm2", obs));
        }
        cols.push("j_total_mA_cm2".to_string());
        cols.push("log10_abs_j".to_string());
        for sp in species {
            cols.push(format!("theta_{}", sp));
        }
        cols.push("theta_empty".to_string());

        writeln!(f, "# {}", cols.join("\t")).unwrap();
        writeln!(f, "# T = {} K, A_site = {:.4e} m^2", self.temp, self.a_site).unwrap();

        for i in 0..self.u_values.len() {
            let j = self.j_total[i];
            let log_j = if j.abs() > 0.0 {
                j.abs().log10()
            } else {
                -30.0
            };

            let mut vals = vec![format!("{:.4}", self.u_values[i])];
            for obs in &obs_names {
                vals.push(format!("{:.8e}", self.tof[*obs][i]));
                vals.push(format!("{:.8e}", self.j_per_obs[*obs][i]));
            }
            vals.push(format!("{:.8e}", j));
            vals.push(format!("{:.4}", log_j));

            let cov = &self.coverages[i];
            let mut theta_sum = 0.0;
            for sp in species {
                let v = *cov.get(sp).unwrap_or(&0.0);
                vals.push(format!("{:.8e}", v));
                theta_sum += v;
            }
            vals.push(format!("{:.8e}", (1.0 - theta_sum).max(0.0)));

            writeln!(f, "{}", vals.join("\t")).unwrap();
        }
        println!("Results saved to {}", filename);
    }

    /// Print summary table to stdout.
    pub fn print_table(&self) {
        let obs_names: Vec<&String> = self.tof.keys().collect();

        println!(
            "\n{:>8} {:>15} {:>15} {:>10}",
            "U [V]", "TOF [s⁻¹]", "j [mA/cm²]", "log₁₀|j|"
        );
        println!("{}", "-".repeat(52));

        for i in 0..self.u_values.len() {
            let j = self.j_total[i];
            let log_j = if j.abs() > 0.0 {
                j.abs().log10()
            } else {
                -30.0
            };
            let tof_main = obs_names
                .first()
                .map(|obs| self.tof[*obs][i])
                .unwrap_or(0.0);
            println!(
                "{:>8.3} {:>15.4e} {:>15.4e} {:>10.3}",
                self.u_values[i], tof_main, j, log_j
            );
        }
    }
}

/// Compute a polarization curve by sweeping potential.
///
/// Uses continuation (previous steady-state as next initial guess)
/// for robust convergence across the potential range.
pub fn compute_polarization_curve(
    mkm: &mut MicroKineticModel,
    u_range: &[f64],
    temp: f64,
    n_electrons: &HashMap<String, f64>,
    a_site: f64,
    verbose: bool,
) -> PolarizationResult {
    mkm.parameters.insert("T".into(), temp);

    let mut tof_arrays: HashMap<String, Vec<f64>> = n_electrons
        .keys()
        .map(|k| (k.clone(), Vec::new()))
        .collect();
    let mut j_arrays: HashMap<String, Vec<f64>> = n_electrons
        .keys()
        .map(|k| (k.clone(), Vec::new()))
        .collect();
    let mut j_total = Vec::new();
    let mut u_values = Vec::new();
    let mut coverages = Vec::new();

    let mut theta_guess: Option<Vec<f64>> = None;
    let report_interval = (u_range.len() / 10).max(1);

    if verbose {
        eprintln!(
            "Computing polarization curve: {} points",
            u_range.len()
        );
        eprintln!(
            "  T = {} K, U = [{:.2}, {:.2}] V vs RHE",
            temp,
            u_range.first().unwrap_or(&0.0),
            u_range.last().unwrap_or(&0.0)
        );
        eprintln!("  Site area = {:.4e} m²\n", a_site);
    }

    for (i, &u) in u_range.iter().enumerate() {
        mkm.parameters.insert("U".into(), u);

        let ss = mkm.solve_steady_state_with(theta_guess.as_deref());

        // Continuation: use current solution as next initial guess
        theta_guess = Some(
            mkm.species
                .iter()
                .map(|sp| *ss.get(sp).unwrap_or(&0.0))
                .collect(),
        );

        let tof = mkm.get_tof(&ss);

        let mut j_sum = 0.0;
        for (obs, &n_e) in n_electrons {
            let tof_val = *tof.get(obs).unwrap_or(&0.0);
            let j_val = tof_to_current_density(tof_val, n_e, a_site);
            tof_arrays.get_mut(obs).unwrap().push(tof_val);
            j_arrays.get_mut(obs).unwrap().push(j_val);
            j_sum += j_val;
        }

        u_values.push(u);
        j_total.push(j_sum);
        coverages.push(ss);

        if verbose && (i + 1) % report_interval == 0 {
            let tof_str: String = n_electrons
                .keys()
                .map(|obs| {
                    format!(
                        "{}={:.3e}",
                        obs,
                        tof.get(obs).unwrap_or(&0.0)
                    )
                })
                .collect::<Vec<_>>()
                .join(", ");
            eprintln!(
                "  [{:5.1}%] U = {:+.3} V, j = {:.4e} mA/cm², TOF: {}",
                100.0 * (i + 1) as f64 / u_range.len() as f64,
                u,
                j_sum,
                tof_str
            );
        }
    }

    if verbose {
        let j_max = j_total.iter().map(|j| j.abs()).fold(0.0_f64, f64::max);
        eprintln!("\nDone. Max |j| = {:.4e} mA/cm²", j_max);
    }

    PolarizationResult {
        u_values,
        tof: tof_arrays,
        j_per_obs: j_arrays,
        j_total,
        coverages,
        temp,
        a_site,
    }
}

// ═══════════════════════════════════════════════════════════════════
//  N₂ reduction DFT model builder
// ═══════════════════════════════════════════════════════════════════

/// Build a microkinetic model for N₂ reduction using DFT energy landscape.
///
/// Uses interpolated potential-dependent barriers from constant-potential DFT.
/// Adsorption/desorption steps use gas-phase kinetics (Hertz-Knudsen / TST).
pub fn build_n2_mkm_dft(landscape: &EnergyLandscape) -> MicroKineticModel {
    let mut mkm = MicroKineticModel::new();

    // Surface species (simplified model matching DFT data)
    for sp in &["N2", "NNH", "HNNH", "NNH2", "N", "NH", "NH2", "NH3"] {
        mkm.add_species(sp);
    }

    mkm.parameters.insert("T".into(), 300.0);
    mkm.parameters.insert("U".into(), -0.5);
    mkm.parameters.insert("p_N2".into(), 1.0);
    mkm.parameters.insert("E_bind_N2".into(), 0.8);
    mkm.parameters.insert("E_bind_NH3".into(), 0.5);

    // N₂ adsorption/desorption (gas-phase kinetics, not from DFT landscape)
    mkm.add_reaction(
        "N2_adsorption",
        &[("empty", 1.0)],
        &[("N2", 1.0)],
        Box::new(|p: &HashMap<String, f64>| {
            hertz_knudsen(
                p["p_N2"] * 1e5,
                p["T"],
                28.014 * UMASS,
                (3.2e-10_f64).powi(2),
            )
        }),
        Some(Box::new(|p: &HashMap<String, f64>| {
            tst_rate(p["E_bind_N2"], p["T"])
        })),
        HashMap::new(),
    );

    // Electrochemical steps with interpolated DFT barriers
    // (is_dft, fs_dft, reaction_name, reactant_species, product_species, tof_observable)
    let steps: &[(&str, &str, &str, &str, &str, Option<(&str, f64)>)] = &[
        ("N2*", "NNH*", "N2_to_NNH", "N2", "NNH",
         Some(("N2_consumption", 1.0))),
        ("NNH*", "HNNH*", "NNH_to_HNNH", "NNH", "HNNH", None),
        ("NNH*", "NNH2*", "NNH_to_NNH2", "NNH", "NNH2", None),
        ("NNH2*", "N*", "NNH2_to_N", "NNH2", "N",
         Some(("NH3_production", 1.0))),
        ("N*", "NH*", "N_to_NH", "N", "NH", None),
        ("NH*", "NH2*", "NH_to_NH2", "NH", "NH2", None),
        ("NH2*", "NH3*", "NH2_to_NH3", "NH2", "NH3", None),
    ];

    for &(is, fs, name, reactant, product, ref tof_opt) in steps {
        if landscape.has_transition(is, fs) {
            let fwd = landscape.get_fwd_rate_fn(is, fs).unwrap();
            let rev = landscape.get_rev_rate_fn(is, fs);
            let tof_map: HashMap<String, f64> = tof_opt
                .map(|(k, v)| [(k.to_string(), v)].into_iter().collect())
                .unwrap_or_default();
            mkm.add_reaction(name, &[(reactant, 1.0)], &[(product, 1.0)], fwd, rev, tof_map);
        }
    }

    // NH₃ desorption
    mkm.add_reaction(
        "NH3_desorption",
        &[("NH3", 1.0)],
        &[("empty", 1.0)],
        Box::new(|p: &HashMap<String, f64>| tst_rate(p["E_bind_NH3"], p["T"])),
        None,
        [("NH3_production".into(), 1.0)].into_iter().collect(),
    );

    mkm
}

// ═══════════════════════════════════════════════════════════════════
//  Tests
// ═══════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_spline_linear_data() {
        // Spline through linear data should reproduce exactly
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![0.0, 2.0, 4.0, 6.0, 8.0]; // y = 2x
        let spline = CubicSpline::new(x, y);

        assert!((spline.evaluate(0.5) - 1.0).abs() < 1e-10);
        assert!((spline.evaluate(2.5) - 5.0).abs() < 1e-10);
        assert!((spline.evaluate(3.7) - 7.4).abs() < 1e-10);
    }

    #[test]
    fn test_cubic_spline_passes_through_knots() {
        let x = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        let y = vec![4.0, 1.0, 0.0, 1.0, 4.0]; // y = x²
        let spline = CubicSpline::new(x.clone(), y.clone());

        for (xi, yi) in x.iter().zip(y.iter()) {
            assert!(
                (spline.evaluate(*xi) - yi).abs() < 1e-10,
                "Spline({}) = {} != {}",
                xi,
                spline.evaluate(*xi),
                yi
            );
        }
    }

    #[test]
    fn test_cubic_spline_two_points() {
        let spline = CubicSpline::new(vec![0.0, 1.0], vec![0.0, 1.0]);
        assert!((spline.evaluate(0.5) - 0.5).abs() < 1e-10);
        // Extrapolation
        assert!((spline.evaluate(2.0) - 2.0).abs() < 1e-10);
        assert!((spline.evaluate(-1.0) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_interpolated_barrier_clamps_to_zero() {
        // Barrier that goes negative by extrapolation
        let barrier = InterpolatedBarrier::new(
            vec![0.0, 1.0, 2.0],
            vec![2.0, 1.0, 0.0],
        );
        // At x=3.0 (extrapolated), raw value would be negative
        assert!(barrier.evaluate(3.0) >= 0.0);
    }

    #[test]
    fn test_tof_to_current_density() {
        // 1 s⁻¹ TOF, 1 electron, 1e-20 m² site
        let j = tof_to_current_density(1.0, 1.0, 1e-20);
        let expected = E_CHARGE / 1e-20 * 0.1; // e / A_site * 0.1
        assert!((j - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_quadratic_fit_exact() {
        // Data generated from y = 0.5*x² - 1.0*x + 2.0
        let u = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        let e: Vec<f64> = u.iter().map(|x| 0.5 * x * x - 1.0 * x + 2.0).collect();
        let (qe, r2) = QuadraticEnergy::fit(&u, &e);

        assert!((qe.a - 0.5).abs() < 1e-10, "a = {} != 0.5", qe.a);
        assert!((qe.b - (-1.0)).abs() < 1e-10, "b = {} != -1.0", qe.b);
        assert!((qe.c - 2.0).abs() < 1e-10, "c = {} != 2.0", qe.c);
        assert!((r2 - 1.0).abs() < 1e-10, "R² = {} != 1.0", r2);
    }

    #[test]
    fn test_quadratic_fit_noisy() {
        // Slightly noisy quadratic data — R² should still be high
        let u = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        let e = vec![4.01, 1.02, -0.01, 0.98, 3.99]; // ≈ x²
        let (qe, r2) = QuadraticEnergy::fit(&u, &e);

        assert!(r2 > 0.999, "R² = {} too low for near-perfect data", r2);
        assert!((qe.a - 1.0).abs() < 0.02);
    }

    #[test]
    fn test_barrier_from_quadratic() {
        let barrier = InterpolatedBarrier::from_quadratic(0.1, -0.5, 2.0);
        // At U=0: Ea = 2.0
        assert!((barrier.evaluate(0.0) - 2.0).abs() < 1e-10);
        // At U=2.5: Ea = 0.1*6.25 - 1.25 + 2.0 = 1.375
        assert!((barrier.evaluate(2.5) - 1.375).abs() < 1e-10);
        // Clamped to >= 0
        assert!(barrier.evaluate(100.0) >= 0.0);
    }

    #[test]
    fn test_landscape_quadratic_fit() {
        let potentials = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        let mut state_energies = HashMap::new();
        // IS: constant at 0
        state_energies.insert("A*".to_string(), vec![0.0, 0.0, 0.0, 0.0, 0.0]);
        // FS: linear in U
        state_energies.insert("B*".to_string(), vec![-2.0, -1.0, 0.0, 1.0, 2.0]);

        let mut ts_energies = HashMap::new();
        // TS: quadratic, Ea_fwd = TS - IS = 0.1*U² + 2.0
        let ts: Vec<f64> = potentials.iter().map(|u| 0.1 * u * u + 2.0).collect();
        ts_energies.insert(("A*".to_string(), "B*".to_string()), ts);

        let landscape = EnergyLandscape::new_quadratic_fit(
            potentials, &state_energies, &ts_energies);

        assert!(landscape.has_transition("A*", "B*"));
        // At U=0: Ea_fwd = TS(0) - IS(0) = 2.0 - 0.0 = 2.0
        let ea = landscape.barriers_fwd[&("A*".to_string(), "B*".to_string())].evaluate(0.0);
        assert!((ea - 2.0).abs() < 0.01, "Ea_fwd(0) = {} != 2.0", ea);
    }

    #[test]
    fn test_landscape_from_coefficients() {
        let mut state_coeff = HashMap::new();
        state_coeff.insert("A*".to_string(), [0.0, 0.0, 0.0]); // E=0
        state_coeff.insert("B*".to_string(), [0.0, 0.0, -1.0]); // E=-1

        let mut ts_coeff = HashMap::new();
        ts_coeff.insert(("A*".to_string(), "B*".to_string()), [0.0, 0.0, 1.5]);

        let landscape = EnergyLandscape::from_quadratic_coefficients(
            &state_coeff, &ts_coeff);

        // Ea_fwd = 1.5 - 0.0 = 1.5, Ea_rev = 1.5 - (-1.0) = 2.5
        let fwd = landscape.barriers_fwd[&("A*".to_string(), "B*".to_string())].evaluate(0.0);
        let rev = landscape.barriers_rev[&("A*".to_string(), "B*".to_string())].evaluate(0.0);
        assert!((fwd - 1.5).abs() < 1e-10);
        assert!((rev - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_energy_landscape() {
        let potentials = vec![-2.0, -1.0, 0.0];
        let mut state_energies = HashMap::new();
        state_energies.insert("A*".to_string(), vec![0.0, 0.0, 0.0]);
        state_energies.insert("B*".to_string(), vec![-1.0, -0.5, 0.0]);

        let mut ts_energies = HashMap::new();
        ts_energies.insert(
            ("A*".to_string(), "B*".to_string()),
            vec![1.0, 1.5, 2.0],
        );

        let landscape = EnergyLandscape::new(potentials, &state_energies, &ts_energies);

        assert!(landscape.has_transition("A*", "B*"));
        assert!(!landscape.has_transition("B*", "A*"));

        // At U=-2: Ea_fwd = 1.0 - 0.0 = 1.0, Ea_rev = 1.0 - (-1.0) = 2.0
        let fwd = landscape.barriers_fwd[&("A*".to_string(), "B*".to_string())].evaluate(-2.0);
        let rev = landscape.barriers_rev[&("A*".to_string(), "B*".to_string())].evaluate(-2.0);
        assert!((fwd - 1.0).abs() < 1e-10);
        assert!((rev - 2.0).abs() < 1e-10);
    }
}
