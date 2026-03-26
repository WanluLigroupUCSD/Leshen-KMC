/// Physical constants and unit conversions (SI base).

// Fundamental constants
pub const KB: f64 = 1.380649e-23;          // Boltzmann constant [J/K]
pub const H_PLANCK: f64 = 6.62607015e-34;  // Planck constant [J·s]
pub const NA: f64 = 6.02214076e23;         // Avogadro number [1/mol]
pub const E_CHARGE: f64 = 1.602176634e-19; // Elementary charge [C]

// Energy
pub const EV: f64 = E_CHARGE;             // 1 eV in Joules

// Length
pub const ANGSTROM: f64 = 1e-10;          // 1 Å in meters

// Pressure
pub const BAR: f64 = 1e5;                 // 1 bar in Pa

// Mass
pub const UMASS: f64 = 1.66053906660e-27; // 1 amu in kg

pub const PI: f64 = std::f64::consts::PI;

/// Molecular masses in amu
pub fn molecular_mass(name: &str) -> f64 {
    match name {
        "H" => 1.008,
        "H2" => 2.016,
        "C" => 12.011,
        "N" => 14.007,
        "N2" => 28.014,
        "O" => 15.999,
        "O2" => 31.998,
        "CO" => 28.010,
        "CO2" => 44.009,
        "H2O" => 18.015,
        "NH3" => 17.031,
        "N2H4" => 32.045,
        "NNH" => 29.022,
        "Mo" => 95.95,
        "W" => 183.84,
        _ => 28.0, // default
    }
}
