/// mykmc-rs: Kinetic Monte Carlo & Microkinetic Modeling
///
/// Two modes:
///   default  — full model (adsorption + desorption + reactions)
///   --cycle  — pure catalytic cycle (only chemical steps, Σθᵢ=1)
///
/// Usage:
///   mykmc kmc                          # full model KMC
///   mykmc kmc --cycle --u=-1.0         # cycle-only KMC
///   mykmc mkm --u=-1.0                 # full mean-field
///   mykmc mkm --cycle --u=-1.0         # cycle-only mean-field
///   mykmc scan-u --cycle               # potential scan (cycle mode)
///   mykmc validate                     # CO/Pd100 Langmuir test

mod units;
mod model;
mod rates;
mod engine;
mod microkinetic;
mod analysis;
mod models;
mod polarization;

use clap::{Parser, Subcommand};
use crate::engine::KMCEngine;
use crate::analysis::TrajectoryRecorder;
use crate::models::*;

#[derive(Parser)]
#[command(name = "mykmc", about = "Kinetic Monte Carlo & Microkinetic Modeling for N2 Reduction")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run lattice KMC simulation
    Kmc {
        #[arg(long, default_value_t = 300.0)]
        t: f64,
        #[arg(long, default_value_t = 1.0)]
        p_n2: f64,
        #[arg(long, default_value_t = -1.0, allow_hyphen_values = true)]
        u: f64,
        #[arg(long, default_value_t = 20)]
        size: usize,
        #[arg(long, default_value_t = 1000000)]
        steps: u64,
        /// Use pure catalytic cycle mode (no adsorption/desorption/migration)
        #[arg(long, default_value_t = false)]
        cycle: bool,
    },
    /// Run mean-field microkinetic model
    Mkm {
        #[arg(long, default_value_t = 300.0)]
        t: f64,
        #[arg(long, default_value_t = 1.0)]
        p_n2: f64,
        #[arg(long, default_value_t = -1.0, allow_hyphen_values = true)]
        u: f64,
        /// Use pure catalytic cycle mode
        #[arg(long, default_value_t = false)]
        cycle: bool,
    },
    /// Scan applied potential
    ScanU {
        #[arg(long, default_value_t = 300.0)]
        t: f64,
        #[arg(long, default_value_t = -0.5, allow_hyphen_values = true)]
        u_min: f64,
        #[arg(long, default_value_t = -2.0, allow_hyphen_values = true)]
        u_max: f64,
        #[arg(long, default_value_t = 31)]
        u_steps: usize,
        #[arg(long, default_value_t = false)]
        cycle: bool,
    },
    /// Scan temperature
    ScanT {
        #[arg(long, default_value_t = -1.0, allow_hyphen_values = true)]
        u: f64,
        #[arg(long, default_value_t = 250.0)]
        t_min: f64,
        #[arg(long, default_value_t = 500.0)]
        t_max: f64,
        #[arg(long, default_value_t = 11)]
        t_steps: usize,
        #[arg(long, default_value_t = false)]
        cycle: bool,
    },
    /// Validate: CO/Pd(100) vs Langmuir isotherm
    Validate,
    /// Compute polarization curve from DFT data
    Polarization {
        /// Path to DFT energy JSON file
        #[arg(long)]
        dft_data: String,
        #[arg(long, default_value_t = 300.0)]
        t: f64,
        #[arg(long, default_value_t = 1.0)]
        p_n2: f64,
        #[arg(long, default_value_t = -2.0, allow_hyphen_values = true)]
        u_min: f64,
        #[arg(long, default_value_t = 0.0)]
        u_max: f64,
        #[arg(long, default_value_t = 50)]
        u_steps: usize,
        /// Fitting mode: "spline" (cubic spline, default) or "quadratic" (Ω=aU²+bU+c)
        #[arg(long, default_value = "spline")]
        fit: String,
        /// Output filename
        #[arg(long, default_value = "polarization_curve.dat")]
        output: String,
    },
}

// ══════════════════════════════════════════════════════════════════

fn run_kmc(temp: f64, p_n2: f64, u: f64, size: usize, steps: u64, cycle: bool) {
    let mode = if cycle { "cycle (pure chemistry)" } else { "full (with ads/des)" };
    println!("{}", "=".repeat(65));
    println!("  Lattice KMC: N2 Reduction on Mo [{}]", mode);
    println!("{}", "=".repeat(65));

    let pt = if cycle { build_n2_kmc_cycle() } else { build_n2_kmc_full() };
    let mut engine = KMCEngine::new(pt, &[size, size]);
    engine.set_parameter("T", temp);
    engine.set_parameter("U", u);
    if !cycle { engine.set_parameter("p_N2", p_n2); }
    engine.print_rates();

    println!("T={} K, U={} V vs RHE, Lattice={}x{}", temp, u, size, size);
    println!("Running {} steps ...\n", steps);

    let mut recorder = TrajectoryRecorder::new();
    let chunk = (steps / 100).max(1);

    for i in 0..100u64 {
        engine.do_steps(chunk, false);
        recorder.record(&mut engine);
        if (i + 1) % 10 == 0 {
            let pct = (i + 1) as f64;
            let cov = engine.get_coverage();
            let dominant = cov.iter()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .map(|(k, v)| format!("{}({:.4})", k, v))
                .unwrap_or_default();
            eprintln!("  [{:3.0}%] time={:.4e} s  {}", pct, engine.kmc_time, dominant);
        }
    }

    println!("\n--- Results ---");
    println!("Steps: {}  Time: {:.6e} s", engine.kmc_step, engine.kmc_time);
    engine.print_coverages();

    let tof = engine.get_tof();
    if !tof.is_empty() {
        println!("\nTOF:");
        for (obs, val) in &tof { println!("  {}: {:.6e} site^-1 s^-1", obs, val); }
    }
    println!("\nProcess counts:");
    for (name, count) in engine.get_process_stats() {
        if count > 0 { println!("  {}: {}", name, count); }
    }
}

fn run_mkm(temp: f64, p_n2: f64, u: f64, cycle: bool) {
    let mode = if cycle { "cycle" } else { "full" };
    println!("{}", "=".repeat(65));
    println!("  Mean-Field Microkinetic: N2 Reduction on Mo [{}]", mode);
    println!("{}", "=".repeat(65));

    let mut mkm = if cycle { build_n2_mkm_cycle() } else { build_n2_mkm_full() };
    mkm.parameters.insert("T".into(), temp);
    mkm.parameters.insert("U".into(), u);
    if !cycle { mkm.parameters.insert("p_N2".into(), p_n2); }

    println!("T={} K, U={} V vs RHE  mode={}\n", temp, u, mode);

    let ss = mkm.solve_steady_state();
    mkm.print_summary(&ss);
}

fn scan_potential(temp: f64, u_min: f64, u_max: f64, u_steps: usize, cycle: bool) {
    let mode = if cycle { "cycle" } else { "full" };
    println!("{}", "=".repeat(65));
    println!("  Potential Scan [{}] T={} K", mode, temp);
    println!("{}", "=".repeat(65));

    let mut mkm = if cycle { build_n2_mkm_cycle() } else { build_n2_mkm_full() };
    mkm.parameters.insert("T".into(), temp);

    println!("\n{:>8} {:>15} {:>12} {:>8}",
             "U [V]", "TOF_NH3 [s-1]", "dominant", "theta");
    println!("{}", "-".repeat(48));

    for i in 0..u_steps {
        let u = u_min + (u_max - u_min) * i as f64 / (u_steps - 1).max(1) as f64;
        mkm.parameters.insert("U".into(), u);
        let ss = mkm.solve_steady_state();
        let tof = mkm.get_tof(&ss);
        let tof_nh3 = tof.get("NH3_production").copied().unwrap_or(0.0);
        let (dom, dv) = ss.iter()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(k, v)| (k.as_str(), *v))
            .unwrap_or(("?", 0.0));
        println!("{:>8.2} {:>15.4e} {:>12} {:>8.4}", u, tof_nh3, dom, dv);
    }
}

fn scan_temperature(u: f64, t_min: f64, t_max: f64, t_steps: usize, cycle: bool) {
    let mode = if cycle { "cycle" } else { "full" };
    println!("{}", "=".repeat(65));
    println!("  Temperature Scan [{}] U={} V", mode, u);
    println!("{}", "=".repeat(65));

    let mut mkm = if cycle { build_n2_mkm_cycle() } else { build_n2_mkm_full() };
    mkm.parameters.insert("U".into(), u);

    println!("\n{:>8} {:>15} {:>12} {:>8}",
             "T [K]", "TOF_NH3 [s-1]", "dominant", "theta");
    println!("{}", "-".repeat(48));

    for i in 0..t_steps {
        let t = t_min + (t_max - t_min) * i as f64 / (t_steps - 1).max(1) as f64;
        mkm.parameters.insert("T".into(), t);
        let ss = mkm.solve_steady_state();
        let tof = mkm.get_tof(&ss);
        let tof_nh3 = tof.get("NH3_production").copied().unwrap_or(0.0);
        let (dom, dv) = ss.iter()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(k, v)| (k.as_str(), *v))
            .unwrap_or(("?", 0.0));
        println!("{:>8.1} {:>15.4e} {:>12} {:>8.4}", t, tof_nh3, dom, dv);
    }
}

fn run_polarization(dft_path: &str, temp: f64, p_n2: f64,
                     u_min: f64, u_max: f64, u_steps: usize,
                     fit_mode: &str, output: &str) {
    println!("{}", "=".repeat(65));
    println!("  Polarization Curve from DFT Data");
    println!("{}", "=".repeat(65));
    println!("  DFT data: {}", dft_path);
    println!("  T = {} K, p_N2 = {} bar", temp, p_n2);
    println!("  U range: [{:.2}, {:.2}] V ({} points)", u_min, u_max, u_steps);
    println!("  Fit mode: {}", fit_mode);
    println!();

    // Load DFT data (auto-detects format)
    let dft_data = polarization::load_dft_data(dft_path);

    // Build energy landscape based on data format and fit mode
    let landscape = match dft_data {
        polarization::DftData::Tabulated { potentials, state_energies, ts_energies } => {
            match fit_mode {
                "quadratic" => {
                    println!("Using quadratic regression: Ω(U) = aU² + bU + c\n");
                    polarization::EnergyLandscape::new_quadratic_fit(
                        potentials, &state_energies, &ts_energies)
                }
                _ => {
                    println!("Using cubic spline interpolation\n");
                    polarization::EnergyLandscape::new(
                        potentials, &state_energies, &ts_energies)
                }
            }
        }
        polarization::DftData::QuadraticCoefficients { state_coefficients, ts_coefficients } => {
            println!("Using direct quadratic coefficients\n");
            polarization::EnergyLandscape::from_quadratic_coefficients(
                &state_coefficients, &ts_coefficients)
        }
    };
    landscape.summary();

    // Build MKM with interpolated barriers
    let mut mkm = polarization::build_n2_mkm_dft(&landscape);
    mkm.parameters.insert("p_N2".into(), p_n2);

    // Generate U range
    let u_range: Vec<f64> = (0..u_steps)
        .map(|i| u_min + (u_max - u_min) * i as f64 / (u_steps - 1).max(1) as f64)
        .collect();

    // N₂ + 6H⁺ + 6e⁻ → 2NH₃  →  3 e⁻ per NH₃ molecule
    let n_electrons: std::collections::HashMap<String, f64> =
        [("NH3_production".into(), 3.0)].into_iter().collect();
    let a_site = (3.2e-10_f64).powi(2);

    let result = polarization::compute_polarization_curve(
        &mut mkm, &u_range, temp, &n_electrons, a_site, true);

    result.print_table();
    result.save(output, &mkm.species);
}

fn validate() {
    println!("{}", "=".repeat(60));
    println!("  Validation: CO/Pd(100) vs Langmuir Isotherm");
    println!("{}", "=".repeat(60));

    let langmuir = |t: f64, dg: f64| -> f64 {
        let k = (-dg * units::EV / (units::KB * t)).exp();
        k / (1.0 + k)
    };

    let t = 600.0;
    println!("\nT={} K, p_CO=1 bar, 30x30 lattice\n", t);
    println!("{:>10} {:>10} {:>15} {:>8}", "dG[eV]", "KMC", "Langmuir", "Status");
    println!("{}", "-".repeat(48));

    let mut all_pass = true;
    for dg in &[-0.5, -0.3, -0.1, 0.0, 0.1_f64] {
        let exact = langmuir(t, *dg);
        let pt = build_co_test();
        let mut eng = KMCEngine::new(pt, &[30, 30]);
        eng.set_parameter("T", t);
        eng.set_parameter("deltaG", *dg);
        eng.do_steps(500000, false);
        let mut samples = Vec::new();
        for _ in 0..20 {
            eng.do_steps(10000, false);
            samples.push(*eng.get_coverage().get("CO").unwrap_or(&0.0));
        }
        let kmc: f64 = samples.iter().sum::<f64>() / samples.len() as f64;
        let err = (kmc - exact).abs() / exact.max(1e-10);
        let st = if err < 0.05 { "PASS" } else { all_pass = false; "FAIL" };
        println!("{:>10.1} {:>10.4} {:>15.4} {:>8}", dg, kmc, exact, st);
    }
    println!("\n{}", if all_pass { "All PASSED!" } else { "Some FAILED." });
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Kmc { t, p_n2, u, size, steps, cycle } =>
            run_kmc(t, p_n2, u, size, steps, cycle),
        Commands::Mkm { t, p_n2, u, cycle } =>
            run_mkm(t, p_n2, u, cycle),
        Commands::ScanU { t, u_min, u_max, u_steps, cycle } =>
            scan_potential(t, u_min, u_max, u_steps, cycle),
        Commands::ScanT { u, t_min, t_max, t_steps, cycle } =>
            scan_temperature(u, t_min, t_max, t_steps, cycle),
        Commands::Validate => validate(),
        Commands::Polarization { dft_data, t, p_n2, u_min, u_max, u_steps, fit, output } =>
            run_polarization(&dft_data, t, p_n2, u_min, u_max, u_steps, &fit, &output),
    }
}
