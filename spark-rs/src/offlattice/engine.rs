//! SKMC Engine — on-the-fly off-lattice KMC simulation (Rust core).
//!
//! This module provides the high-performance KMC event selection and
//! basin/superbasin management. The computationally expensive saddle
//! point searches and energy evaluations are delegated to the Python
//! layer (via ASE calculators), while the catalogue matching, KMC
//! selection, and superbasin acceleration run at native Rust speed.

use super::catalogue::Catalogue;
use super::cache::{SuperCache, CacheOptions};
use super::environment::{Geometry, build_geometry};
use super::mechanism::Mechanism;

/// SKMC simulation state.
pub struct SKMCState {
    pub kmc_time: f64,
    pub kmc_step: u64,
    pub catalogue: Catalogue,
    pub super_cache: SuperCache,
}

/// Result of a single KMC step selection.
pub struct StepResult {
    pub mechanism: Mechanism,
    pub atom_index: usize,
    pub dt: f64,
    pub source_basin: usize,
}

impl SKMCState {
    /// Create a new SKMC state.
    pub fn new(
        r_env: f64,
        r_edge: f64,
        delta_max: f64,
        temperature: f64,
        barrier_tol: f64,
        max_barrier: f64,
    ) -> Self {
        let catalogue = Catalogue::new(r_env, r_edge, delta_max);
        let cache_opt = CacheOptions {
            temperature,
            barrier_tol,
            max_barrier,
            ..CacheOptions::default()
        };
        let super_cache = SuperCache::new(cache_opt);

        SKMCState {
            kmc_time: 0.0,
            kmc_step: 0,
            catalogue,
            super_cache,
        }
    }

    /// Build local environments from atomic positions.
    pub fn build_environments(
        &self,
        positions: &[[f64; 3]],
        types: &[i32],
        frozen: &[bool],
        neighbor_lists: &[Vec<usize>],
        r_env: f64,
    ) -> Vec<Geometry> {
        let n = positions.len();
        let mut geometries = Vec::with_capacity(n);
        for i in 0..n {
            let geo = build_geometry(
                i, positions, types, frozen, &neighbor_lists[i], r_env,
            );
            geometries.push(geo);
        }
        geometries
    }

    /// Rebuild catalogue and return new environment indices.
    pub fn rebuild_catalogue(&mut self, geometries: Vec<Geometry>) -> Vec<usize> {
        self.catalogue.rebuild(geometries)
    }

    /// Initialize the SuperCache with the first basin.
    pub fn initialize_cache(&mut self, positions: Vec<[f64; 3]>) {
        self.super_cache.initialize(positions, &self.catalogue);
    }

    /// Select a KMC event.
    pub fn select_event(&self) -> StepResult {
        let mut rng = rand::thread_rng();
        let choice = self.super_cache.kmc_choice(&mut rng);
        StepResult {
            mechanism: choice.mechanism,
            atom_index: choice.atom_index,
            dt: choice.dt,
            source_basin: choice.source_basin,
        }
    }

    /// Accept a step: advance time and update connectivity.
    pub fn accept_step(
        &mut self,
        dt: f64,
        source_basin: usize,
        atom_index: usize,
        mechanism: &Mechanism,
        new_positions: Vec<[f64; 3]>,
        state_hash: u64,
    ) {
        self.kmc_time += dt;
        self.kmc_step += 1;
        self.super_cache.connect_from(
            source_basin, atom_index, mechanism,
            new_positions, &self.catalogue, state_hash,
        );
    }

    /// Set mechanisms for a catalogue entry.
    pub fn set_mechanisms(&mut self, atom_index: usize, mechanisms: Vec<Mechanism>) {
        self.catalogue.set_mechanisms(atom_index, mechanisms);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mechanism_rate() {
        let m = Mechanism {
            barrier: 0.5,
            delta: -0.1,
            kinetic_pre: 1e13,
            delta_sp: vec![[0.1, 0.0, 0.0]],
            delta_fwd: vec![[0.2, 0.0, 0.0]],
            axis: None,
            err_fwd: 0.0,
            err_sp: 0.0,
            poison_fwd: false,
            poison_sp: false,
        };

        let rate_300 = m.rate(300.0);
        let rate_500 = m.rate(500.0);
        assert!(rate_300 > 0.0);
        assert!(rate_500 > rate_300);
        // At 300K, barrier 0.5 eV should give rate ~ 4e4 Hz
        assert!(rate_300 > 1e3 && rate_300 < 1e6);
    }

    #[test]
    fn test_fingerprint() {
        use super::super::environment::Fingerprint;

        let pos = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let fp = Fingerprint::new(&pos);

        assert_eq!(fp.r_0j.len(), 2);
        assert_eq!(fp.r_ij.len(), 1);
        assert!((fp.r_0j[0] - 1.0).abs() < 1e-10);

        // Self-equivalence
        assert!(fp.equiv(&fp, 0.1));
    }

    #[test]
    fn test_kabsch() {
        use super::super::environment::kabsch;

        let x = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let y = vec![[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let (o, rmsd) = kabsch(&x, &y);
        assert!(rmsd < 1e-6, "Kabsch RMSD too high: {}", rmsd);
        // Check O is orthogonal
        let mut check = [0.0f64; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    check[i] += o[j][i] * o[j][k] * if i == k { 1.0 } else { 0.0 };
                }
            }
        }
    }

    #[test]
    fn test_geometry_hash() {
        use super::super::environment::Geometry;

        let pos = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let mut geo1 = Geometry::new(
            pos.clone(),
            vec![0, 2, 2],
            vec![0, 1, 2],
        );
        let mut geo2 = Geometry::new(
            pos,
            vec![0, 2, 2],
            vec![0, 1, 2],
        );

        let h1 = geo1.graph_hash(1.5);
        let h2 = geo2.graph_hash(1.5);
        assert_eq!(h1, h2, "Same geometry should have same hash");
    }
}
