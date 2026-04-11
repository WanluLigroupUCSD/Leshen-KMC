//! Mechanism — a min→saddle→min transition.

use serde::{Deserialize, Serialize};

/// Inverse Boltzmann constant in eV^-1 (matches openFLY).
pub const INV_BOLTZ: f64 = 16021766340.0 / 1380649.0;

/// A transition mechanism between two local minima via a saddle point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mechanism {
    /// Activation energy E_sp - E_0 (eV).
    pub barrier: f64,
    /// Energy change E_final - E_0 (eV).
    pub delta: f64,
    /// Arrhenius pre-factor (Hz).
    pub kinetic_pre: f64,
    /// Displacement initial→saddle for each local atom, shape [n_local][3].
    pub delta_sp: Vec<[f64; 3]>,
    /// Displacement initial→final for each local atom, shape [n_local][3].
    pub delta_fwd: Vec<[f64; 3]>,
    /// Dimer axis at saddle point (optional).
    pub axis: Option<Vec<[f64; 3]>>,
    /// Reconstruction error for forward state (Å).
    pub err_fwd: f64,
    /// Reconstruction error for saddle point (Å).
    pub err_sp: f64,
    /// Forward state unreconstructable.
    pub poison_fwd: bool,
    /// Saddle point unreconstructable.
    pub poison_sp: bool,
}

impl Mechanism {
    /// Number of atoms in the local environment.
    pub fn n_atoms(&self) -> usize {
        self.delta_fwd.len()
    }

    /// Barrier for the reverse transition.
    pub fn reverse_barrier(&self) -> f64 {
        self.barrier - self.delta
    }

    /// Arrhenius rate at given temperature (K).
    pub fn rate(&self, temperature: f64) -> f64 {
        self.kinetic_pre * (-self.barrier * INV_BOLTZ / temperature).exp()
    }
}
