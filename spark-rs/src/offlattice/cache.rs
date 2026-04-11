//! SuperCache — superbasin manager with dynamic barrier tolerance.

use super::basin::Basin;
use super::superbasin::SuperBasin;
use super::mechanism::Mechanism;
use super::catalogue::Catalogue;
use rand::Rng;

/// Configuration for the SuperCache.
pub struct CacheOptions {
    pub temperature: f64,
    pub barrier_tol: f64,
    pub state_tol: f64,
    pub cache_size: usize,
    pub max_superbasin_size: usize,
    pub dynamic_tol: bool,
    pub tol_grow: f64,
    pub tol_shrink: f64,
    pub max_barrier: f64,
}

impl Default for CacheOptions {
    fn default() -> Self {
        CacheOptions {
            temperature: 300.0,
            barrier_tol: 0.3,
            state_tol: 0.5,
            cache_size: 10,
            max_superbasin_size: 20,
            dynamic_tol: true,
            tol_grow: 1.5,
            tol_shrink: 0.7,
            max_barrier: 5.0,
        }
    }
}

/// SuperBasin cache with dynamic barrier tolerance.
pub struct SuperCache {
    pub opt: CacheOptions,
    sb: Option<SuperBasin>,
    cache: Vec<SuperBasin>,
    in_cache_count: usize,
}

impl SuperCache {
    pub fn new(opt: CacheOptions) -> Self {
        SuperCache {
            opt,
            sb: None,
            cache: Vec::new(),
            in_cache_count: 0,
        }
    }

    pub fn initialize(&mut self, positions: Vec<[f64; 3]>, catalogue: &Catalogue) {
        let basin = Basin::new(
            positions, catalogue,
            self.opt.temperature, self.opt.max_barrier,
        );
        self.sb = Some(SuperBasin::new(basin));
    }

    pub fn kmc_choice<R: Rng>(&self, rng: &mut R) -> super::superbasin::SBChoice {
        self.sb.as_ref().expect("Not initialized").kmc_choice(rng)
    }

    /// Process a completed transition.
    pub fn connect_from(
        &mut self,
        source_basin: usize,
        atom_index: usize,
        mechanism: &Mechanism,
        positions: Vec<[f64; 3]>,
        catalogue: &Catalogue,
        state_hash: u64,
    ) {
        let sb = self.sb.as_mut().expect("Not initialized");

        // Case 1: existing basin
        if let Some(_prev) = sb.find_occupy(state_hash, &positions, self.opt.state_tol) {
            sb.connect_from(source_basin, atom_index, mechanism);
            return;
        }

        let max_barrier = mechanism.barrier.max(mechanism.barrier - mechanism.delta);

        if max_barrier < self.opt.barrier_tol {
            // Case 2: low barrier → expand
            if self.opt.dynamic_tol && sb.size() >= self.opt.max_superbasin_size {
                self.opt.barrier_tol = (self.opt.barrier_tol * self.opt.tol_shrink).max(0.0);
                let basin = Basin::new(
                    positions, catalogue,
                    self.opt.temperature, self.opt.max_barrier,
                );
                self.sb = Some(SuperBasin::new(basin));
            } else {
                let new_basin = Basin::new(
                    positions, catalogue,
                    self.opt.temperature, self.opt.max_barrier,
                );
                sb.expand_occupy(new_basin);
                sb.connect_from(source_basin, atom_index, mechanism);
            }
            return;
        }

        // Case 3: high barrier → cache current, try restore
        let cached = self.try_restore(state_hash, &positions);

        if let Some(restored) = cached {
            let old_sb = std::mem::replace(&mut self.sb, Some(restored));
            if let Some(sb) = old_sb {
                self.cache_sb(sb);
            }
            self.in_cache_count += 1;
        } else {
            let old_sb = self.sb.take();
            if let Some(sb) = old_sb {
                self.cache_sb(sb);
            }
            let basin = Basin::new(
                positions, catalogue,
                self.opt.temperature, self.opt.max_barrier,
            );
            self.sb = Some(SuperBasin::new(basin));
            self.in_cache_count = 0;
        }

        if self.opt.dynamic_tol && self.in_cache_count > self.opt.cache_size {
            self.opt.barrier_tol *= self.opt.tol_grow;
        }
    }

    fn cache_sb(&mut self, sb: SuperBasin) {
        self.cache.push(sb);
        while self.cache.len() > self.opt.cache_size {
            self.cache.remove(0);
        }
    }

    fn try_restore(
        &mut self, state_hash: u64, positions: &[[f64; 3]],
    ) -> Option<SuperBasin> {
        for i in 0..self.cache.len() {
            if self.cache[i]
                .find_occupy(state_hash, positions, self.opt.state_tol)
                .is_some()
            {
                return Some(self.cache.remove(i));
            }
        }
        None
    }
}
