/// Core KMC engine: BKL rejection-free algorithm.
///
/// Every step advances the clock and changes configuration (no rejected moves).
/// Uses O(1) available-sites bookkeeping via swap-with-last.

use std::collections::HashMap;
use rand::Rng;

use crate::model::Project;
use crate::rates::evaluate_rate_expression;

/// Compiled process data for fast lookup (no String comparisons in hot path).
struct CompiledProcess {
    /// (offset_flat_delta, species_id) for each condition
    conditions: Vec<(Vec<i32>, usize)>,
    /// (offset_flat_delta, new_species_id) for each action
    actions: Vec<(Vec<i32>, usize)>,
    /// TOF counters
    tof_count: HashMap<String, f64>,
}

pub struct KMCEngine {
    // Model
    project: Project,
    ndim: usize,
    lattice_size: Vec<usize>,
    nsites: usize,
    nspecies: usize,
    nproc: usize,

    // State
    lattice: Vec<u16>,
    pub kmc_time: f64,
    pub kmc_step: u64,
    pub procstat: Vec<u64>,

    // For TOF calculation
    prev_procstat: Vec<u64>,
    prev_time: f64,

    // Rate constants
    pub rates: Vec<f64>,
    accum_rates: Vec<f64>,

    // Available sites bookkeeping (per process)
    avail_sites: Vec<Vec<usize>>,
    site_in_avail: Vec<HashMap<usize, usize>>,

    // Compiled processes
    compiled: Vec<CompiledProcess>,
    max_offset: usize,

    // RNG
    rng: rand::rngs::ThreadRng,
}

impl KMCEngine {
    pub fn new(mut project: Project, size: &[usize]) -> Self {
        let ndim = project.model_dimension;
        let lattice_size: Vec<usize> = size[..ndim].to_vec();
        let nsites: usize = lattice_size.iter().product();
        let nspecies = project.species_list.len();
        let nproc = project.process_list.len();

        project.rebuild_maps();

        // Compile processes to numeric form
        let compiled: Vec<CompiledProcess> = project.process_list.iter().map(|proc| {
            let conditions = proc.conditions.iter().map(|c| {
                let sp_id = project.species_map[&c.species];
                (c.offset.clone(), sp_id)
            }).collect();
            let actions = proc.actions.iter().map(|a| {
                let sp_id = project.species_map[&a.species];
                (a.offset.clone(), sp_id)
            }).collect();
            CompiledProcess {
                conditions,
                actions,
                tof_count: proc.tof_count.clone(),
            }
        }).collect();

        // Compute max offset range
        let mut max_off: usize = 1;
        for cp in &compiled {
            for (off, _) in cp.conditions.iter().chain(cp.actions.iter()) {
                for &o in off {
                    max_off = max_off.max(o.unsigned_abs() as usize + 1);
                }
            }
        }

        // Initialize lattice with default species (id=0)
        let lattice = vec![0u16; nsites];

        // Initialize bookkeeping
        let avail_sites = vec![Vec::new(); nproc];
        let site_in_avail = vec![HashMap::new(); nproc];

        let mut engine = Self {
            project,
            ndim,
            lattice_size,
            nsites,
            nspecies,
            nproc,
            lattice,
            kmc_time: 0.0,
            kmc_step: 0,
            procstat: vec![0u64; nproc],
            prev_procstat: vec![0u64; nproc],
            prev_time: 0.0,
            rates: vec![0.0; nproc],
            accum_rates: vec![0.0; nproc],
            avail_sites,
            site_in_avail,
            compiled,
            max_offset: max_off,
            rng: rand::thread_rng(),
        };

        engine.update_rate_constants();
        engine.rebuild_avail_sites();
        engine
    }

    // ── Coordinate conversion (with PBC) ──────────────────────────────

    #[inline]
    fn site_to_coord(&self, site: usize) -> Vec<i32> {
        match self.ndim {
            1 => vec![site as i32],
            2 => vec![
                (site / self.lattice_size[1]) as i32,
                (site % self.lattice_size[1]) as i32,
            ],
            3 => {
                let ly = self.lattice_size[1];
                let lz = self.lattice_size[2];
                vec![
                    (site / (ly * lz)) as i32,
                    ((site % (ly * lz)) / lz) as i32,
                    (site % lz) as i32,
                ]
            }
            _ => unreachable!(),
        }
    }

    #[inline]
    fn coord_to_site(&self, coord: &[i32]) -> usize {
        match self.ndim {
            1 => {
                coord[0].rem_euclid(self.lattice_size[0] as i32) as usize
            }
            2 => {
                let x = coord[0].rem_euclid(self.lattice_size[0] as i32) as usize;
                let y = coord[1].rem_euclid(self.lattice_size[1] as i32) as usize;
                x * self.lattice_size[1] + y
            }
            3 => {
                let x = coord[0].rem_euclid(self.lattice_size[0] as i32) as usize;
                let y = coord[1].rem_euclid(self.lattice_size[1] as i32) as usize;
                let z = coord[2].rem_euclid(self.lattice_size[2] as i32) as usize;
                x * self.lattice_size[1] * self.lattice_size[2]
                    + y * self.lattice_size[2] + z
            }
            _ => unreachable!(),
        }
    }

    // ── Available sites bookkeeping ───────────────────────────────────

    #[inline]
    fn check_process_at_site(&self, proc_id: usize, site: usize) -> bool {
        let coord = self.site_to_coord(site);
        let cp = &self.compiled[proc_id];
        for (offset, sp_id) in &cp.conditions {
            let nc: Vec<i32> = coord.iter().zip(offset.iter())
                .map(|(&c, &o)| c + o).collect();
            if self.lattice[self.coord_to_site(&nc)] as usize != *sp_id {
                return false;
            }
        }
        true
    }

    fn add_to_avail(&mut self, proc_id: usize, site: usize) {
        if !self.site_in_avail[proc_id].contains_key(&site) {
            let idx = self.avail_sites[proc_id].len();
            self.site_in_avail[proc_id].insert(site, idx);
            self.avail_sites[proc_id].push(site);
        }
    }

    fn remove_from_avail(&mut self, proc_id: usize, site: usize) {
        if let Some(idx) = self.site_in_avail[proc_id].remove(&site) {
            let last = *self.avail_sites[proc_id].last().unwrap();
            self.avail_sites[proc_id][idx] = last;
            if last != site {
                self.site_in_avail[proc_id].insert(last, idx);
            }
            self.avail_sites[proc_id].pop();
        }
    }

    fn rebuild_avail_sites(&mut self) {
        for p in 0..self.nproc {
            self.avail_sites[p].clear();
            self.site_in_avail[p].clear();
        }
        for s in 0..self.nsites {
            for p in 0..self.nproc {
                if self.check_process_at_site(p, s) {
                    self.add_to_avail(p, s);
                }
            }
        }
    }

    fn get_affected_sites(&self, site: usize, proc_id: usize) -> Vec<usize> {
        let coord = self.site_to_coord(site);
        let cp = &self.compiled[proc_id];
        let r = self.max_offset as i32;
        let mut affected = Vec::new();

        // Get coordinates of changed sites
        let mut changed: Vec<Vec<i32>> = Vec::new();
        for (offset, _) in &cp.actions {
            let nc: Vec<i32> = coord.iter().zip(offset.iter())
                .map(|(&c, &o)| c + o).collect();
            changed.push(nc);
        }

        // Expand to neighbors within max_offset range
        match self.ndim {
            1 => {
                for cc in &changed {
                    for dx in -r..=r {
                        affected.push(self.coord_to_site(&[cc[0] + dx]));
                    }
                }
            }
            2 => {
                for cc in &changed {
                    for dx in -r..=r {
                        for dy in -r..=r {
                            affected.push(self.coord_to_site(&[cc[0] + dx, cc[1] + dy]));
                        }
                    }
                }
            }
            3 => {
                for cc in &changed {
                    for dx in -r..=r {
                        for dy in -r..=r {
                            for dz in -r..=r {
                                affected.push(self.coord_to_site(
                                    &[cc[0] + dx, cc[1] + dy, cc[2] + dz]));
                            }
                        }
                    }
                }
            }
            _ => unreachable!(),
        }

        affected.sort_unstable();
        affected.dedup();
        affected
    }

    fn update_avail_after_exec(&mut self, affected: &[usize]) {
        for &site in affected {
            for p in 0..self.nproc {
                if self.check_process_at_site(p, site) {
                    self.add_to_avail(p, site);
                } else {
                    self.remove_from_avail(p, site);
                }
            }
        }
    }

    // ── Rates ─────────────────────────────────────────────────────────

    pub fn update_rate_constants(&mut self) {
        for i in 0..self.nproc {
            self.rates[i] = evaluate_rate_expression(
                &self.project.process_list[i].rate_constant,
                &self.project.parameter_list,
            );
        }
    }

    fn update_accum_rates(&mut self) -> f64 {
        let mut total = 0.0;
        for p in 0..self.nproc {
            total += self.rates[p] * self.avail_sites[p].len() as f64;
            self.accum_rates[p] = total;
        }
        total
    }

    // ── KMC Step (BKL) ───────────────────────────────────────────────

    /// Execute one BKL KMC step. Returns false if system is frozen.
    pub fn do_kmc_step(&mut self) -> bool {
        let total_rate = self.update_accum_rates();
        if total_rate <= 0.0 {
            return false;
        }

        let r_time: f64 = self.rng.gen();
        let r_proc: f64 = self.rng.gen();
        let r_site: f64 = self.rng.gen();

        // Time advancement
        self.kmc_time += -r_time.ln() / total_rate;

        // Process selection (binary search)
        let target = r_proc * total_rate;
        let proc_id = match self.accum_rates[..self.nproc]
            .binary_search_by(|v| v.partial_cmp(&target).unwrap())
        {
            Ok(i) => i,
            Err(i) => i.min(self.nproc - 1),
        };

        // Site selection
        let n_avail = self.avail_sites[proc_id].len();
        if n_avail == 0 {
            return false;
        }
        let site_idx = (r_site * n_avail as f64) as usize;
        let site_idx = site_idx.min(n_avail - 1);
        let site = self.avail_sites[proc_id][site_idx];

        // Execute: update lattice
        let coord = self.site_to_coord(site);
        let cp = &self.compiled[proc_id];
        for (offset, new_sp) in &cp.actions {
            let nc: Vec<i32> = coord.iter().zip(offset.iter())
                .map(|(&c, &o)| c + o).collect();
            let ns = self.coord_to_site(&nc);
            self.lattice[ns] = *new_sp as u16;
        }

        // Update bookkeeping
        let affected = self.get_affected_sites(site, proc_id);
        self.update_avail_after_exec(&affected);

        self.kmc_step += 1;
        self.procstat[proc_id] += 1;
        true
    }

    /// Execute n KMC steps.
    pub fn do_steps(&mut self, n: u64, progress: bool) {
        let report = if progress { (n / 10).max(1) } else { u64::MAX };
        for i in 0..n {
            if !self.do_kmc_step() {
                eprintln!("System frozen at step {}, time={:.6e} s",
                          self.kmc_step, self.kmc_time);
                break;
            }
            if progress && (i + 1) % report == 0 {
                let pct = 100.0 * (i + 1) as f64 / n as f64;
                eprintln!("  [{:5.1}%] step={}, time={:.4e} s",
                          pct, self.kmc_step, self.kmc_time);
            }
        }
    }

    // ── Observables ──────────────────────────────────────────────────

    /// Get fractional coverage for each species.
    pub fn get_coverage(&self) -> HashMap<String, f64> {
        let mut counts = vec![0usize; self.nspecies];
        for &sp in &self.lattice {
            counts[sp as usize] += 1;
        }
        let mut cov = HashMap::new();
        for sp in &self.project.species_list {
            cov.insert(sp.name.clone(), counts[sp.id] as f64 / self.nsites as f64);
        }
        cov
    }

    /// Get TOFs since last call.
    pub fn get_tof(&mut self) -> HashMap<String, f64> {
        let dt = self.kmc_time - self.prev_time;
        if dt <= 0.0 {
            return HashMap::new();
        }
        let mut tof = HashMap::new();
        for p in 0..self.nproc {
            let delta = self.procstat[p] - self.prev_procstat[p];
            for (obs, &coeff) in &self.compiled[p].tof_count {
                *tof.entry(obs.clone()).or_insert(0.0) +=
                    coeff * delta as f64 / (dt * self.nsites as f64);
            }
        }
        self.prev_procstat = self.procstat.clone();
        self.prev_time = self.kmc_time;
        tof
    }

    /// Get process execution counts.
    pub fn get_process_stats(&self) -> Vec<(&str, u64)> {
        self.project.process_list.iter()
            .map(|p| (p.name.as_str(), self.procstat[p.id]))
            .collect()
    }

    // ── Parameter modification ───────────────────────────────────────

    pub fn set_parameter(&mut self, name: &str, value: f64) {
        self.project.set_parameter(name, value);
        self.update_rate_constants();
    }

    pub fn get_parameter(&self, name: &str) -> f64 {
        self.project.get_parameter(name)
    }

    // ── State management ─────────────────────────────────────────────

    pub fn reset(&mut self) {
        self.lattice.fill(0);
        self.kmc_time = 0.0;
        self.kmc_step = 0;
        self.procstat.fill(0);
        self.prev_procstat.fill(0);
        self.prev_time = 0.0;
        self.rebuild_avail_sites();
    }

    // ── Printing ─────────────────────────────────────────────────────

    pub fn print_rates(&self) {
        println!("\nRate constants:");
        println!("  {:<35} {:>12}  {:>8}", "Process", "k [s^-1]", "N_avail");
        println!("  {} {}  {}", "-".repeat(35), "-".repeat(12), "-".repeat(8));
        for i in 0..self.nproc {
            println!("  {:<35} {:>12.4e}  {:>8}",
                     self.project.process_list[i].name,
                     self.rates[i],
                     self.avail_sites[i].len());
        }
        println!();
    }

    pub fn print_coverages(&self) {
        let cov = self.get_coverage();
        println!("Coverages:");
        for sp in &self.project.species_list {
            let v = cov[&sp.name];
            if v > 1e-6 {
                println!("  {}: {:.6}", sp.name, v);
            }
        }
    }

    pub fn nsites(&self) -> usize { self.nsites }
}
