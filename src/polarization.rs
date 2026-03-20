/// Polarization curve calculation from constant-potential DFT data.
///
/// Workflow:
///   1. Load DFT energies at multiple potentials from JSON
///   2. Build energy landscape with interpolated Ea(U) via cubic spline
///   3. Create microkinetic model with potential-dependent rates
///   4. Sweep potential → steady state → TOF → current density
///
/// Two input modes:
///   (a) Direct Ea(U) tables:  InterpolatedBarrier
///   (b) Energy surfaces E_IS(U), E_TS(U), E_FS(U):  EnergyLandscape

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
//  InterpolatedBarrier: Ea(U) from tabulated DFT data
// ═══════════════════════════════════════════════════════════════════

/// Interpolated activation barrier Ea(U) from tabulated DFT data.
///
/// Uses natural cubic spline (falls back to linear for < 3 points).
/// Extrapolates beyond data range using boundary polynomial segments.
#[derive(Clone)]
pub struct InterpolatedBarrier {
    spline: CubicSpline,
}

impl InterpolatedBarrier {
    pub fn new(u_data: Vec<f64>, ea_data: Vec<f64>) -> Self {
        Self {
            spline: CubicSpline::new(u_data, ea_data),
        }
    }

    /// Return Ea [eV] at potential U [V]. Clamped to >= 0.
    pub fn evaluate(&self, u: f64) -> f64 {
        self.spline.evaluate(u).max(0.0)
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

    /// Print barrier summary at tabulated potentials.
    pub fn summary(&self) {
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
            for &u in &self.potentials {
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

/// Load DFT energy data from a JSON file.
///
/// Expected format:
/// ```json
/// {
///   "potentials": [-2.0, -1.5, -1.0, -0.5, 0.0],
///   "state_energies": { "N2*": [...], "NNH*": [...], ... },
///   "ts_energies": { "N2*->NNH*": [...], ... }
/// }
/// ```
///
/// Returns (potentials, state_energies, ts_energies) with TS keys as (IS, FS) tuples.
pub fn load_dft_data(
    path: &str,
) -> (
    Vec<f64>,
    HashMap<String, Vec<f64>>,
    HashMap<(String, String), Vec<f64>>,
) {
    let text =
        fs::read_to_string(path).unwrap_or_else(|e| panic!("Cannot read '{}': {}", path, e));
    let data: serde_json::Value =
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("Invalid JSON in '{}': {}", path, e));

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
        // Parse "IS->FS" format
        if let Some((is, fs)) = key.split_once("->") {
            ts_energies.insert((is.trim().to_string(), fs.trim().to_string()), energies);
        } else {
            eprintln!("Warning: skipping TS key '{}' (expected 'IS->FS' format)", key);
        }
    }

    (potentials, state_energies, ts_energies)
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
