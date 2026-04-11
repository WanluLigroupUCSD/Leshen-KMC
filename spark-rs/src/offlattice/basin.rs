//! Basin — a local minimum with mechanisms and KMC selection.

use rand::Rng;
use super::mechanism::{Mechanism, INV_BOLTZ};
use super::catalogue::Catalogue;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// A mechanism localized to a specific atom.
pub struct LocalisedMech {
    pub atom_index: usize,
    pub rate: f64,
    pub mechanism: Mechanism,
    pub exit_mech: bool,
}

/// A local energy minimum with all accessible mechanisms.
pub struct Basin {
    pub positions: Vec<[f64; 3]>,
    pub mechs: Vec<LocalisedMech>,
    pub rate_sum: f64,
    pub state_hash: u64,
    pub connected: bool,
    temperature: f64,
}

/// KMC choice result.
pub struct Choice {
    pub mechanism: Mechanism,
    pub atom_index: usize,
    pub dt: f64,
}

impl Basin {
    /// Build basin from positions and catalogue.
    pub fn new(
        positions: Vec<[f64; 3]>,
        catalogue: &Catalogue,
        temperature: f64,
        max_barrier: f64,
    ) -> Self {
        let n = positions.len();
        let mut mechs = Vec::new();
        let mut rate_sum = 0.0;

        // Collect mechanisms
        for i in 0..n {
            for mech in catalogue.get_mechanisms(i) {
                if mech.barrier > 0.0 && mech.barrier < max_barrier {
                    let rate = mech.kinetic_pre
                        * (-mech.barrier * INV_BOLTZ / temperature).exp();
                    mechs.push(LocalisedMech {
                        atom_index: i,
                        rate,
                        mechanism: mech.clone(),
                        exit_mech: true,
                    });
                    rate_sum += rate;
                }
            }
        }

        // State hash
        let mut hasher = DefaultHasher::new();
        for i in 0..n {
            if let Some(entry) = catalogue.get_entry(i) {
                entry.cat_index.hash(&mut hasher);
            } else {
                (usize::MAX).hash(&mut hasher);
            }
        }
        let state_hash = hasher.finish();

        Basin {
            positions,
            mechs,
            rate_sum,
            state_hash,
            connected: false,
            temperature,
        }
    }

    /// KMC event selection (n-fold way).
    pub fn kmc_choice<R: Rng>(&self, rng: &mut R) -> Choice {
        assert!(self.rate_sum > 0.0 && !self.mechs.is_empty());

        let target = rng.gen::<f64>() * self.rate_sum;
        let mut cumul = 0.0;
        let mut selected = &self.mechs[self.mechs.len() - 1];

        for lm in &self.mechs {
            cumul += lm.rate;
            if cumul >= target {
                selected = lm;
                break;
            }
        }

        let dt = -(rng.gen::<f64>().ln()) / self.rate_sum;

        Choice {
            mechanism: selected.mechanism.clone(),
            atom_index: selected.atom_index,
            dt,
        }
    }
}
