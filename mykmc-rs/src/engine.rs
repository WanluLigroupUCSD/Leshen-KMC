/// Core KMC engine: BKL rejection-free algorithm.
///
/// Features:
///   - BKL/VSSM rejection-free algorithm
///   - Neighbor list for spatial events (4-NN 2D, 6-NN 3D)
///   - Pairwise lateral interactions with per-site rates
///   - Surface diffusion support
///   - BEP (Bronsted-Evans-Polanyi) relations
///   - Site type support
///
/// Uses O(1) available-sites bookkeeping via swap-with-last.

use std::collections::HashMap;
use rand::Rng;

use crate::model::Project;
use crate::rates::evaluate_rate_expression;
use crate::units::{KB, EV};

/// Compiled process data for fast lookup (no String comparisons in hot path).
struct CompiledProcess {
    conditions: Vec<(Vec<i32>, usize)>,
    actions: Vec<(Vec<i32>, usize)>,
    tof_count: HashMap<String, f64>,
    site_type: Option<i32>,
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
    site_types: Vec<i32>,
    pub kmc_time: f64,
    pub kmc_step: u64,
    pub procstat: Vec<u64>,

    // For TOF calculation
    prev_procstat: Vec<u64>,
    prev_time: f64,

    // Rate constants (base rates, before lateral correction)
    pub rates: Vec<f64>,
    accum_rates: Vec<f64>,

    // Available sites bookkeeping (per process)
    avail_sites: Vec<Vec<usize>>,
    site_in_avail: Vec<HashMap<usize, usize>>,
    /// Parallel rate array: avail_rates[p][i] = rate for avail_sites[p][i]
    avail_rates: Vec<Vec<f64>>,
    /// Sum of per-site rates per process (only used with lateral)
    proc_total_rates: Vec<f64>,

    // Neighbor list
    neighbors: Vec<Vec<usize>>,

    // Lateral interactions: (species_a, species_b) -> energy [eV]
    lateral_dict: HashMap<(usize, usize), f64>,
    has_lateral: bool,

    // BEP: proc_id -> alpha
    bep_proc: HashMap<usize, f64>,

    // Empty species id
    empty_species: usize,

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

        // Detect empty species
        let empty_species = project.species_list.iter()
            .find(|s| matches!(s.name.to_lowercase().as_str(), "empty" | "vacant" | "*"))
            .map(|s| s.id)
            .unwrap_or(0);

        // Build lateral interaction dict
        let mut lateral_dict: HashMap<(usize, usize), f64> = HashMap::new();
        for li in &project.lateral_interactions {
            if let (Some(&sp1), Some(&sp2)) = (
                project.species_map.get(&li.species1),
                project.species_map.get(&li.species2),
            ) {
                lateral_dict.insert((sp1, sp2), li.energy);
                lateral_dict.insert((sp2, sp1), li.energy);
            }
        }
        let has_lateral = !lateral_dict.is_empty();

        // Build BEP lookup: proc_id -> alpha
        let bep_by_name: HashMap<String, f64> = project.bep_relations.iter()
            .map(|b| (b.process_name.clone(), b.alpha))
            .collect();

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
                site_type: proc.site_type,
            }
        }).collect();

        // Resolve BEP proc IDs
        let mut bep_proc: HashMap<usize, f64> = HashMap::new();
        for (pid, proc) in project.process_list.iter().enumerate() {
            if let Some(&alpha) = bep_by_name.get(&proc.name) {
                bep_proc.insert(pid, alpha);
            }
        }

        // Compute max offset range
        let mut max_off: usize = 1;
        for cp in &compiled {
            for (off, _) in cp.conditions.iter().chain(cp.actions.iter()) {
                for &o in off {
                    max_off = max_off.max(o.unsigned_abs() as usize + 1);
                }
            }
        }

        // Initialize lattice
        let lattice = vec![0u16; nsites];
        let site_types = vec![0i32; nsites];

        // Build neighbor list
        let neighbors = build_neighbor_list(ndim, &lattice_size, nsites);

        // Initialize bookkeeping
        let avail_sites = vec![Vec::new(); nproc];
        let site_in_avail = vec![HashMap::new(); nproc];
        let avail_rates = vec![Vec::new(); nproc];
        let proc_total_rates = vec![0.0; nproc];

        let mut engine = Self {
            project,
            ndim,
            lattice_size,
            nsites,
            nspecies,
            nproc,
            lattice,
            site_types,
            kmc_time: 0.0,
            kmc_step: 0,
            procstat: vec![0u64; nproc],
            prev_procstat: vec![0u64; nproc],
            prev_time: 0.0,
            rates: vec![0.0; nproc],
            accum_rates: vec![0.0; nproc],
            avail_sites,
            site_in_avail,
            avail_rates,
            proc_total_rates,
            neighbors,
            lateral_dict,
            has_lateral,
            bep_proc,
            empty_species,
            compiled,
            max_offset: max_off,
            rng: rand::thread_rng(),
        };

        engine.update_rate_constants();
        engine.rebuild_avail_sites();

        if has_lateral {
            engine.rebuild_per_site_rates();
        }

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
            1 => coord[0].rem_euclid(self.lattice_size[0] as i32) as usize,
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
        let cp = &self.compiled[proc_id];

        // Check site type requirement
        if let Some(req) = cp.site_type {
            if self.site_types[site] != req {
                return false;
            }
        }

        // Check species conditions
        let coord = self.site_to_coord(site);
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

            if self.has_lateral {
                let rate = self.compute_site_rate(proc_id, site);
                self.avail_rates[proc_id].push(rate);
                self.proc_total_rates[proc_id] += rate;
            } else {
                self.avail_rates[proc_id].push(0.0);
            }
        }
    }

    fn remove_from_avail(&mut self, proc_id: usize, site: usize) {
        if let Some(idx) = self.site_in_avail[proc_id].remove(&site) {
            if self.has_lateral {
                self.proc_total_rates[proc_id] -= self.avail_rates[proc_id][idx];
            }

            let last = *self.avail_sites[proc_id].last().unwrap();
            self.avail_sites[proc_id][idx] = last;
            self.avail_rates[proc_id][idx] = *self.avail_rates[proc_id].last().unwrap();
            if last != site {
                self.site_in_avail[proc_id].insert(last, idx);
            }
            self.avail_sites[proc_id].pop();
            self.avail_rates[proc_id].pop();
        }
    }

    fn rebuild_avail_sites(&mut self) {
        for p in 0..self.nproc {
            self.avail_sites[p].clear();
            self.site_in_avail[p].clear();
            self.avail_rates[p].clear();
        }
        for s in 0..self.nsites {
            for p in 0..self.nproc {
                if self.check_process_at_site(p, s) {
                    self.add_to_avail(p, s);
                }
            }
        }
    }

    fn rebuild_per_site_rates(&mut self) {
        for p in 0..self.nproc {
            let mut total = 0.0;
            for i in 0..self.avail_sites[p].len() {
                let site = self.avail_sites[p][i];
                let rate = self.compute_site_rate(p, site);
                self.avail_rates[p][i] = rate;
                total += rate;
            }
            self.proc_total_rates[p] = total;
        }
    }

    fn get_affected_sites(&self, site: usize, proc_id: usize) -> Vec<usize> {
        let coord = self.site_to_coord(site);
        let cp = &self.compiled[proc_id];
        let r = self.max_offset as i32;
        let mut affected = Vec::new();

        let mut changed: Vec<Vec<i32>> = Vec::new();
        for (offset, _) in &cp.actions {
            let nc: Vec<i32> = coord.iter().zip(offset.iter())
                .map(|(&c, &o)| c + o).collect();
            changed.push(nc);
        }

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
                let was = self.site_in_avail[p].contains_key(&site);
                let is = self.check_process_at_site(p, site);

                if is && !was {
                    self.add_to_avail(p, site);
                } else if !is && was {
                    self.remove_from_avail(p, site);
                } else if is && was && self.has_lateral {
                    let idx = self.site_in_avail[p][&site];
                    let old = self.avail_rates[p][idx];
                    let new = self.compute_site_rate(p, site);
                    self.avail_rates[p][idx] = new;
                    self.proc_total_rates[p] += new - old;
                }
            }
        }
    }

    // ── Lateral interactions & BEP ────────────────────────────────────

    fn get_temperature(&self) -> f64 {
        self.project.get_parameter("T")
    }

    /// Compute rate for (process, site) including lateral + BEP corrections.
    fn compute_site_rate(&self, proc_id: usize, site: usize) -> f64 {
        let base_rate = self.rates[proc_id];
        if !self.has_lateral {
            return base_rate;
        }

        let t = self.get_temperature();
        if t <= 0.0 {
            return base_rate;
        }
        let beta_th = EV / (KB * t);

        let e_react = self.compute_interaction_energy(proc_id, site, false);

        if let Some(&alpha) = self.bep_proc.get(&proc_id) {
            let e_prod = self.compute_interaction_energy(proc_id, site, true);
            let delta_ea = alpha * (e_prod - e_react);
            return base_rate * (-delta_ea * beta_th).exp();
        }

        // Default: reactant-state interaction modifies barrier
        base_rate * (e_react * beta_th).exp()
    }

    /// Sum pairwise interaction energy for condition (reactant) or action (product) sites.
    fn compute_interaction_energy(&self, proc_id: usize, site: usize, is_product: bool) -> f64 {
        let coord = self.site_to_coord(site);
        let cp = &self.compiled[proc_id];
        let entries = if is_product { &cp.actions } else { &cp.conditions };

        // Collect sites involved in this process
        let proc_sites: Vec<usize> = entries.iter()
            .map(|(off, _)| {
                let nc: Vec<i32> = coord.iter().zip(off.iter())
                    .map(|(&c, &o)| c + o).collect();
                self.coord_to_site(&nc)
            })
            .collect();

        let mut e_total = 0.0;
        for (i, (_off, sp_id)) in entries.iter().enumerate() {
            if *sp_id == self.empty_species {
                continue;
            }
            let entry_site = proc_sites[i];
            for &nn in &self.neighbors[entry_site] {
                if proc_sites.contains(&nn) {
                    continue;
                }
                let nn_sp = self.lattice[nn] as usize;
                if let Some(&energy) = self.lateral_dict.get(&(*sp_id, nn_sp)) {
                    e_total += energy;
                }
            }
        }
        e_total
    }

    // ── Rates ─────────────────────────────────────────────────────────

    pub fn update_rate_constants(&mut self) {
        for i in 0..self.nproc {
            self.rates[i] = evaluate_rate_expression(
                &self.project.process_list[i].rate_constant,
                &self.project.parameter_list,
            );
        }
        if self.has_lateral && !self.avail_sites[0].is_empty() {
            self.rebuild_per_site_rates();
        }
    }

    fn update_accum_rates(&mut self) -> f64 {
        let mut total = 0.0;
        for p in 0..self.nproc {
            if self.has_lateral {
                total += self.proc_total_rates[p];
            } else {
                total += self.rates[p] * self.avail_sites[p].len() as f64;
            }
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

        let site_idx = if self.has_lateral {
            // Select site proportional to per-site rate
            let total_proc = self.proc_total_rates[proc_id];
            if total_proc <= 0.0 {
                return false;
            }
            let target = r_site * total_proc;
            let mut cumul = 0.0;
            let mut chosen = n_avail - 1;
            for k in 0..n_avail {
                cumul += self.avail_rates[proc_id][k];
                if cumul >= target {
                    chosen = k;
                    break;
                }
            }
            chosen
        } else {
            // Uniform selection
            ((r_site * n_avail as f64) as usize).min(n_avail - 1)
        };

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

    pub fn get_process_stats(&self) -> Vec<(&str, u64)> {
        self.project.process_list.iter()
            .map(|p| (p.name.as_str(), self.procstat[p.id]))
            .collect()
    }

    pub fn get_neighbor_coverages(&self, site: usize) -> HashMap<String, usize> {
        let mut counts: HashMap<String, usize> = HashMap::new();
        for &nn in &self.neighbors[site] {
            let sp = &self.project.species_list[self.lattice[nn] as usize].name;
            *counts.entry(sp.clone()).or_insert(0) += 1;
        }
        counts
    }

    // ── Parameter modification ───────────────────────────────────────

    pub fn set_parameter(&mut self, name: &str, value: f64) {
        self.project.set_parameter(name, value);
        self.update_rate_constants();
    }

    pub fn get_parameter(&self, name: &str) -> f64 {
        self.project.get_parameter(name)
    }

    // ── Site type manipulation ────────────────────────────────────────

    pub fn set_site_type(&mut self, site: usize, stype: i32) {
        self.site_types[site] = stype;
    }

    pub fn set_site_types_region<F: Fn(&[i32]) -> bool>(&mut self, f: F, stype: i32) {
        for s in 0..self.nsites {
            let coord = self.site_to_coord(s);
            if f(&coord) {
                self.site_types[s] = stype;
            }
        }
        self.rebuild_avail_sites();
        if self.has_lateral {
            self.rebuild_per_site_rates();
        }
    }

    // ── State management ─────────────────────────────────────────────

    pub fn reset(&mut self) {
        self.lattice.fill(0);
        self.site_types.fill(0);
        self.kmc_time = 0.0;
        self.kmc_step = 0;
        self.procstat.fill(0);
        self.prev_procstat.fill(0);
        self.prev_time = 0.0;
        self.rebuild_avail_sites();
        if self.has_lateral {
            self.rebuild_per_site_rates();
        }
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
        if self.has_lateral {
            println!("  Lateral interactions: {}, BEP: {}",
                     self.lateral_dict.len() / 2, self.bep_proc.len());
            if !self.neighbors.is_empty() {
                println!("  Neighbors per site: {}", self.neighbors[0].len());
            }
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
    pub fn has_lateral(&self) -> bool { self.has_lateral }
    pub fn neighbors_of(&self, site: usize) -> &[usize] { &self.neighbors[site] }
}

// ── Neighbor list construction ──────────────────────────────────────────

fn build_neighbor_list(ndim: usize, lattice_size: &[usize], nsites: usize) -> Vec<Vec<usize>> {
    let nn_offsets: Vec<Vec<i32>> = match ndim {
        1 => vec![vec![-1], vec![1]],
        2 => vec![vec![-1, 0], vec![1, 0], vec![0, -1], vec![0, 1]],
        3 => vec![
            vec![-1, 0, 0], vec![1, 0, 0],
            vec![0, -1, 0], vec![0, 1, 0],
            vec![0, 0, -1], vec![0, 0, 1],
        ],
        _ => vec![],
    };

    let mut neighbors = vec![Vec::new(); nsites];
    for s in 0..nsites {
        let coord = match ndim {
            1 => vec![s as i32],
            2 => vec![(s / lattice_size[1]) as i32, (s % lattice_size[1]) as i32],
            3 => {
                let ly = lattice_size[1];
                let lz = lattice_size[2];
                vec![(s / (ly * lz)) as i32, ((s % (ly * lz)) / lz) as i32, (s % lz) as i32]
            }
            _ => vec![],
        };
        for off in &nn_offsets {
            let nc: Vec<i32> = coord.iter().zip(off.iter()).map(|(&c, &o)| c + o).collect();
            let ns = match ndim {
                1 => nc[0].rem_euclid(lattice_size[0] as i32) as usize,
                2 => {
                    let x = nc[0].rem_euclid(lattice_size[0] as i32) as usize;
                    let y = nc[1].rem_euclid(lattice_size[1] as i32) as usize;
                    x * lattice_size[1] + y
                }
                3 => {
                    let x = nc[0].rem_euclid(lattice_size[0] as i32) as usize;
                    let y = nc[1].rem_euclid(lattice_size[1] as i32) as usize;
                    let z = nc[2].rem_euclid(lattice_size[2] as i32) as usize;
                    x * lattice_size[1] * lattice_size[2] + y * lattice_size[2] + z
                }
                _ => 0,
            };
            if ns != s {
                neighbors[s].push(ns);
            }
        }
    }
    neighbors
}
