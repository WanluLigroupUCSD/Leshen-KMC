/// Built-in models for N₂ reduction on Mo surface.
///
/// Two modes:
///   "full"  — complete model with adsorption, desorption, all processes
///   "cycle" — pure catalytic cycle, only chemical reaction steps
///             no adsorption/desorption/migration/defects
///             all sites always occupied, coverage conserved: Σθᵢ = 1

use std::collections::HashMap;
use crate::model::*;
use crate::microkinetic::MicroKineticModel;
use crate::rates::*;
use crate::units::UMASS;

// ══════════════════════════════════════════════════════════════════════
//  DFT barriers shared by both modes
// ══════════════════════════════════════════════════════════════════════

/// DFT activation barriers [eV] from Mo surface calculations.
pub struct DftBarriers {
    pub n2_to_nnh: f64,
    pub nnh_to_hnnh: f64,
    pub nnh_to_nnh2: f64,
    pub hnnh_to_hnnh2: f64,
    pub nnh2_to_n: f64,          // N-N cleavage (thermal)
    pub hnnh2_to_h2nnh2: f64,
    pub h2nnh2_to_nh2: f64,     // N-N cleavage (thermal)
    pub n_to_nh: f64,
    pub nh_to_nh2: f64,
    pub nh2_to_nh3: f64,
}

impl Default for DftBarriers {
    fn default() -> Self {
        Self {
            n2_to_nnh: 2.82,
            nnh_to_hnnh: 1.14,
            nnh_to_nnh2: 3.22,
            hnnh_to_hnnh2: 3.32,
            nnh2_to_n: 2.94,
            hnnh2_to_h2nnh2: 3.36,
            h2nnh2_to_nh2: 5.14,
            n_to_nh: 2.01,
            nh_to_nh2: 2.68,
            nh2_to_nh3: 4.35,
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
//  MODE 1: Full model (adsorption + desorption + reactions)
// ══════════════════════════════════════════════════════════════════════

/// Full KMC model: includes adsorption, desorption, and all reaction steps.
pub fn build_n2_kmc_full() -> Project {
    let mut pt = Project::new("N2_reduction_Mo_full", 2);

    for name in &["empty", "N2", "NNH", "HNNH", "NNH2", "HNNH2",
                   "H2NNH2", "N", "NH", "NH2", "NH3"] {
        pt.add_species(name);
    }
    pt.cell_size = vec![3.2, 3.2];

    // Parameters
    pt.add_parameter_adjustable("T", 300.0, 200.0, 500.0);
    pt.add_parameter_adjustable("p_N2", 1.0, 1e-3, 10.0);
    pt.add_parameter_adjustable("U", -0.5, -5.0, 0.5);
    pt.add_parameter("A_site", (3.2e-10_f64).powi(2));
    pt.add_parameter("E_bind_N2", 0.8);
    pt.add_parameter("E_bind_NH3", 0.5);
    pt.add_parameter("beta_BV", 0.5);

    let b = DftBarriers::default();
    pt.add_parameter("Ea_N2_to_NNH", b.n2_to_nnh);
    pt.add_parameter("Ea_NNH_to_HNNH", b.nnh_to_hnnh);
    pt.add_parameter("Ea_NNH_to_NNH2", b.nnh_to_nnh2);
    pt.add_parameter("Ea_HNNH_to_HNNH2", b.hnnh_to_hnnh2);
    pt.add_parameter("Ea_NNH2_to_N", b.nnh2_to_n);
    pt.add_parameter("Ea_HNNH2_to_H2NNH2", b.hnnh2_to_h2nnh2);
    pt.add_parameter("Ea_H2NNH2_to_NH2", b.h2nnh2_to_nh2);
    pt.add_parameter("Ea_N_to_NH", b.n_to_nh);
    pt.add_parameter("Ea_NH_to_NH2", b.nh_to_nh2);
    pt.add_parameter("Ea_NH2_to_NH3", b.nh2_to_nh3);

    let o = vec![0, 0];

    // Adsorption / desorption
    pt.add_process("N2_adsorption",
        vec![Condition::new(&o, "empty")], vec![Action::new(&o, "N2")],
        RateExpr::Expression("p_N2*bar*A_site/sqrt(2*pi*m_N2*umass/beta)".into()),
        HashMap::new());
    pt.add_process("N2_desorption",
        vec![Condition::new(&o, "N2")], vec![Action::new(&o, "empty")],
        RateExpr::Expression("p_N2*bar*A_site/sqrt(2*pi*m_N2*umass/beta)*exp(-beta*E_bind_N2*eV)".into()),
        HashMap::new());
    pt.add_process("NH3_desorption",
        vec![Condition::new(&o, "NH3")], vec![Action::new(&o, "empty")],
        RateExpr::Expression("kB*T/h*exp(-E_bind_NH3*eV/(kB*T))".into()),
        [("NH3_production".to_string(), 1.0)].into_iter().collect());

    // PCET hydrogenation steps
    for (name, from, to, ea) in &[
        ("N2_to_NNH",          "N2",     "NNH",    "Ea_N2_to_NNH"),
        ("NNH_to_HNNH",        "NNH",    "HNNH",   "Ea_NNH_to_HNNH"),
        ("NNH_to_NNH2",        "NNH",    "NNH2",   "Ea_NNH_to_NNH2"),
        ("HNNH_to_HNNH2",      "HNNH",   "HNNH2",  "Ea_HNNH_to_HNNH2"),
        ("HNNH2_to_H2NNH2",    "HNNH2",  "H2NNH2", "Ea_HNNH2_to_H2NNH2"),
        ("N_to_NH",             "N",      "NH",     "Ea_N_to_NH"),
        ("NH_to_NH2",           "NH",     "NH2",    "Ea_NH_to_NH2"),
        ("NH2_to_NH3",          "NH2",    "NH3",    "Ea_NH2_to_NH3"),
    ] {
        let tof = if *name == "N2_to_NNH" {
            [("N2_consumption".to_string(), 1.0)].into_iter().collect()
        } else { HashMap::new() };
        pt.add_process(name,
            vec![Condition::new(&o, from)], vec![Action::new(&o, to)],
            RateExpr::Expression(format!("kB*T/h*exp(-({}+beta_BV*U)*eV/(kB*T))", ea)),
            tof);
    }

    // Thermal N-N cleavage
    pt.add_process("NNH2_to_N_NH3",
        vec![Condition::new(&o, "NNH2")], vec![Action::new(&o, "N")],
        RateExpr::Expression("kB*T/h*exp(-Ea_NNH2_to_N*eV/(kB*T))".into()),
        [("NH3_production".to_string(), 1.0)].into_iter().collect());
    pt.add_process("H2NNH2_to_NH2",
        vec![Condition::new(&o, "H2NNH2")], vec![Action::new(&o, "NH2")],
        RateExpr::Expression("kB*T/h*exp(-Ea_H2NNH2_to_NH2*eV/(kB*T))".into()),
        HashMap::new());

    pt
}

/// Full mean-field model (with adsorption/desorption + empty sites).
pub fn build_n2_mkm_full() -> MicroKineticModel {
    let mut mkm = MicroKineticModel::new();
    for sp in &["N2","NNH","HNNH","NNH2","HNNH2","H2NNH2","N","NH","NH2","NH3"] {
        mkm.add_species(sp);
    }
    mkm.parameters.insert("T".into(), 300.0);
    mkm.parameters.insert("p_N2".into(), 1.0);
    mkm.parameters.insert("U".into(), -0.5);
    mkm.parameters.insert("beta_BV".into(), 0.5);
    mkm.parameters.insert("E_bind_N2".into(), 0.8);
    mkm.parameters.insert("E_bind_NH3".into(), 0.5);

    let b = DftBarriers::default();

    fn echem(ea: f64) -> Box<dyn Fn(&HashMap<String, f64>) -> f64> {
        Box::new(move |p: &HashMap<String, f64>|
            electrochemical_rate(ea, p["T"], p["U"], 0.0, p["beta_BV"]))
    }
    fn therm(ea: f64) -> Box<dyn Fn(&HashMap<String, f64>) -> f64> {
        Box::new(move |p: &HashMap<String, f64>| tst_rate(ea, p["T"]))
    }

    // Adsorption/desorption
    mkm.add_reaction("N2_adsorption",
        &[("empty",1.0)], &[("N2",1.0)],
        Box::new(|p: &HashMap<String,f64>| hertz_knudsen(p["p_N2"]*1e5, p["T"], 28.014*UMASS, (3.2e-10_f64).powi(2))),
        Some(Box::new(|p: &HashMap<String,f64>| tst_rate(p["E_bind_N2"], p["T"]))),
        HashMap::new());

    // Chemical steps
    mkm.add_reaction("N2_to_NNH",   &[("N2",1.0)],  &[("NNH",1.0)],  echem(b.n2_to_nnh), None,
        [("N2_consumption".into(),1.0)].into_iter().collect());
    mkm.add_reaction("NNH_to_HNNH", &[("NNH",1.0)], &[("HNNH",1.0)], echem(b.nnh_to_hnnh), None, HashMap::new());
    mkm.add_reaction("NNH_to_NNH2", &[("NNH",1.0)], &[("NNH2",1.0)], echem(b.nnh_to_nnh2), None, HashMap::new());
    mkm.add_reaction("HNNH_to_HNNH2", &[("HNNH",1.0)], &[("HNNH2",1.0)], echem(b.hnnh_to_hnnh2), None, HashMap::new());
    mkm.add_reaction("NNH2_to_N",   &[("NNH2",1.0)], &[("N",1.0)], therm(b.nnh2_to_n), None,
        [("NH3_production".into(),1.0)].into_iter().collect());
    mkm.add_reaction("HNNH2_to_H2NNH2", &[("HNNH2",1.0)], &[("H2NNH2",1.0)], echem(b.hnnh2_to_h2nnh2), None, HashMap::new());
    mkm.add_reaction("H2NNH2_to_NH2", &[("H2NNH2",1.0)], &[("NH2",1.0)], therm(b.h2nnh2_to_nh2), None, HashMap::new());
    mkm.add_reaction("N_to_NH",     &[("N",1.0)],  &[("NH",1.0)],  echem(b.n_to_nh), None, HashMap::new());
    mkm.add_reaction("NH_to_NH2",   &[("NH",1.0)], &[("NH2",1.0)], echem(b.nh_to_nh2), None, HashMap::new());
    mkm.add_reaction("NH2_to_NH3",  &[("NH2",1.0)], &[("NH3",1.0)], echem(b.nh2_to_nh3), None, HashMap::new());

    // NH3 desorption
    mkm.add_reaction("NH3_desorption", &[("NH3",1.0)], &[("empty",1.0)],
        Box::new(|p: &HashMap<String,f64>| tst_rate(p["E_bind_NH3"], p["T"])), None,
        [("NH3_production".into(),1.0)].into_iter().collect());

    mkm
}

// ══════════════════════════════════════════════════════════════════════
//  MODE 2: Catalytic cycle (pure chemistry, no ads/des/migration)
//
//  Σθᵢ = 1 always. No empty sites. NH₃ released instantly, site → N₂.
//
//  Distal:      N₂→NNH→NNH₂→(N + NH₃↑)→NH→NH₂→(N₂ + NH₃↑)
//  Alternating: N₂→NNH→HNNH→HNNH₂→H₂NNH₂→(NH₂ + NH₃↑)→(N₂ + NH₃↑)
//
//  9 surface species: N2, NNH, HNNH, NNH2, HNNH2, H2NNH2, N, NH, NH2
//  NH₃ is NOT a surface species — it's released immediately.
//  The step "NH₂ + H⁺/e⁻ → NH₃↑" regenerates N₂: NH₂* → N₂* + NH₃↑
// ══════════════════════════════════════════════════════════════════════

/// Cycle-mode KMC model: all sites start as N₂, only chemical steps.
pub fn build_n2_kmc_cycle() -> Project {
    let mut pt = Project::new("N2_reduction_Mo_cycle", 2);

    // No "empty" — first species is N₂ (default occupancy)
    for name in &["N2", "NNH", "HNNH", "NNH2", "HNNH2", "H2NNH2",
                   "N", "NH", "NH2"] {
        pt.add_species(name);
    }
    pt.cell_size = vec![3.2, 3.2];

    pt.add_parameter_adjustable("T", 300.0, 200.0, 500.0);
    pt.add_parameter_adjustable("U", -0.5, -5.0, 0.5);
    pt.add_parameter("beta_BV", 0.5);

    let b = DftBarriers::default();
    pt.add_parameter("Ea_N2_to_NNH", b.n2_to_nnh);
    pt.add_parameter("Ea_NNH_to_HNNH", b.nnh_to_hnnh);
    pt.add_parameter("Ea_NNH_to_NNH2", b.nnh_to_nnh2);
    pt.add_parameter("Ea_HNNH_to_HNNH2", b.hnnh_to_hnnh2);
    pt.add_parameter("Ea_NNH2_to_N", b.nnh2_to_n);
    pt.add_parameter("Ea_HNNH2_to_H2NNH2", b.hnnh2_to_h2nnh2);
    pt.add_parameter("Ea_H2NNH2_to_NH2", b.h2nnh2_to_nh2);
    pt.add_parameter("Ea_N_to_NH", b.n_to_nh);
    pt.add_parameter("Ea_NH_to_NH2", b.nh_to_nh2);
    pt.add_parameter("Ea_NH2_to_NH3", b.nh2_to_nh3);

    let o = vec![0, 0];

    // ── PCET hydrogenation steps ──
    for (name, from, to, ea, tof) in [
        ("N2_to_NNH",       "N2",    "NNH",   "Ea_N2_to_NNH",
         vec![("N2_consumption", 1.0)]),
        ("NNH_to_HNNH",     "NNH",   "HNNH",  "Ea_NNH_to_HNNH",   vec![]),
        ("NNH_to_NNH2",     "NNH",   "NNH2",  "Ea_NNH_to_NNH2",   vec![]),
        ("HNNH_to_HNNH2",   "HNNH",  "HNNH2", "Ea_HNNH_to_HNNH2", vec![]),
        ("HNNH2_to_H2NNH2", "HNNH2", "H2NNH2","Ea_HNNH2_to_H2NNH2", vec![]),
        ("N_to_NH",          "N",     "NH",    "Ea_N_to_NH",        vec![]),
        ("NH_to_NH2",        "NH",    "NH2",   "Ea_NH_to_NH2",      vec![]),
    ] {
        let tof_map: HashMap<String,f64> = tof.into_iter().map(|(k,v)|(k.to_string(),v)).collect();
        pt.add_process(name,
            vec![Condition::new(&o, from)], vec![Action::new(&o, to)],
            RateExpr::Expression(format!("kB*T/h*exp(-({}+beta_BV*U)*eV/(kB*T))", ea)),
            tof_map);
    }

    // ── Thermal N-N cleavage steps ──
    // NNH₂* → N* (releases NH₃)
    pt.add_process("NNH2_to_N_NH3",
        vec![Condition::new(&o, "NNH2")], vec![Action::new(&o, "N")],
        RateExpr::Expression("kB*T/h*exp(-Ea_NNH2_to_N*eV/(kB*T))".into()),
        [("NH3_production".to_string(), 1.0)].into_iter().collect());

    // H₂NNH₂* → NH₂* (releases NH₃ fragment)
    pt.add_process("H2NNH2_to_NH2",
        vec![Condition::new(&o, "H2NNH2")], vec![Action::new(&o, "NH2")],
        RateExpr::Expression("kB*T/h*exp(-Ea_H2NNH2_to_NH2*eV/(kB*T))".into()),
        HashMap::new());

    // ── Cycle closure: NH₂* + H⁺/e⁻ → N₂* + NH₃↑ ──
    // NH₃ released instantly, site regenerates to N₂
    pt.add_process("NH2_to_N2_NH3",
        vec![Condition::new(&o, "NH2")], vec![Action::new(&o, "N2")],
        RateExpr::Expression("kB*T/h*exp(-(Ea_NH2_to_NH3+beta_BV*U)*eV/(kB*T))".into()),
        [("NH3_production".to_string(), 1.0)].into_iter().collect());

    pt
}

/// Cycle-mode mean-field model: pure catalytic cycle.
///
/// No empty sites. Σθᵢ = 1. NH₃ not a surface species.
/// NH₂ + H⁺/e⁻ → N₂ + NH₃↑ closes the cycle.
pub fn build_n2_mkm_cycle() -> MicroKineticModel {
    let mut mkm = MicroKineticModel::new();

    // 9 surface species (no empty, no NH3)
    for sp in &["N2","NNH","HNNH","NNH2","HNNH2","H2NNH2","N","NH","NH2"] {
        mkm.add_species(sp);
    }

    mkm.parameters.insert("T".into(), 300.0);
    mkm.parameters.insert("U".into(), -0.5);
    mkm.parameters.insert("beta_BV".into(), 0.5);

    let b = DftBarriers::default();

    fn echem(ea: f64) -> Box<dyn Fn(&HashMap<String, f64>) -> f64> {
        Box::new(move |p: &HashMap<String, f64>|
            electrochemical_rate(ea, p["T"], p["U"], 0.0, p["beta_BV"]))
    }
    fn therm(ea: f64) -> Box<dyn Fn(&HashMap<String, f64>) -> f64> {
        Box::new(move |p: &HashMap<String, f64>| tst_rate(ea, p["T"]))
    }

    // PCET steps
    mkm.add_reaction("N2_to_NNH",   &[("N2",1.0)],  &[("NNH",1.0)],  echem(b.n2_to_nnh), None,
        [("N2_consumption".into(),1.0)].into_iter().collect());
    mkm.add_reaction("NNH_to_HNNH", &[("NNH",1.0)], &[("HNNH",1.0)], echem(b.nnh_to_hnnh), None, HashMap::new());
    mkm.add_reaction("NNH_to_NNH2", &[("NNH",1.0)], &[("NNH2",1.0)], echem(b.nnh_to_nnh2), None, HashMap::new());
    mkm.add_reaction("HNNH_to_HNNH2", &[("HNNH",1.0)], &[("HNNH2",1.0)], echem(b.hnnh_to_hnnh2), None, HashMap::new());
    mkm.add_reaction("HNNH2_to_H2NNH2", &[("HNNH2",1.0)], &[("H2NNH2",1.0)], echem(b.hnnh2_to_h2nnh2), None, HashMap::new());
    mkm.add_reaction("N_to_NH",     &[("N",1.0)],  &[("NH",1.0)],  echem(b.n_to_nh), None, HashMap::new());
    mkm.add_reaction("NH_to_NH2",   &[("NH",1.0)], &[("NH2",1.0)], echem(b.nh_to_nh2), None, HashMap::new());

    // Thermal N-N cleavage
    mkm.add_reaction("NNH2_to_N",   &[("NNH2",1.0)], &[("N",1.0)], therm(b.nnh2_to_n), None,
        [("NH3_production".into(),1.0)].into_iter().collect());
    mkm.add_reaction("H2NNH2_to_NH2", &[("H2NNH2",1.0)], &[("NH2",1.0)], therm(b.h2nnh2_to_nh2), None, HashMap::new());

    // Cycle closure: NH₂ + H⁺/e⁻ → N₂ + NH₃↑
    mkm.add_reaction("NH2_to_N2_NH3", &[("NH2",1.0)], &[("N2",1.0)], echem(b.nh2_to_nh3), None,
        [("NH3_production".into(),1.0)].into_iter().collect());

    mkm
}

// ══════════════════════════════════════════════════════════════════════
//  CO validation model (unchanged)
// ══════════════════════════════════════════════════════════════════════

pub fn build_co_test() -> Project {
    let mut pt = Project::new("CO_on_Pd100", 2);
    pt.add_species("empty");
    pt.add_species("CO");
    pt.cell_size = vec![3.5, 3.5];
    pt.add_parameter_adjustable("T", 600.0, 400.0, 800.0);
    pt.add_parameter_adjustable("p_CO", 1.0, 1e-10, 100.0);
    pt.add_parameter("A", (3.5e-10_f64).powi(2));
    pt.add_parameter_adjustable("deltaG", -0.5, -1.3, 0.3);
    let o = vec![0, 0];
    pt.add_process("CO_adsorption",
        vec![Condition::new(&o, "empty")], vec![Action::new(&o, "CO")],
        RateExpr::Expression("p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)".into()),
        HashMap::new());
    pt.add_process("CO_desorption",
        vec![Condition::new(&o, "CO")], vec![Action::new(&o, "empty")],
        RateExpr::Expression("p_CO*bar*A/sqrt(2*pi*umass*m_CO/beta)*exp(beta*deltaG*eV)".into()),
        HashMap::new());
    pt
}
