/// Core KMC engine: BKL rejection-free algorithm with optimized data structures.
///
/// Algorithmic features:
///   - BKL/VSSM rejection-free algorithm
///   - Fenwick tree (BIT) for O(log N) rate-weighted site selection
///   - Neighbor list for spatial events (4-NN 2D, 6-NN 3D)
///   - Pairwise lateral interactions with per-site rates
///   - BEP (Bronsted-Evans-Polanyi) relations
///   - Surface diffusion and site type support
///
/// Performance vs naive implementation:
///   - O(log N) site selection via Fenwick tree (was O(N) linear scan)
///   - Flat array lateral energy lookup (was HashMap)
///   - Direct-indexed availability arrays (was HashMap)
///   - Fixed [i32;3] coordinates (zero heap allocation in hot path)
///   - Neighbor-list-based affected-site detection (fast path)

use std::collections::HashMap;
use rand::Rng;

use crate::model::Project;
use crate::rates::evaluate_rate_expression;
use crate::units::{KB, EV};

// ═══════════════════════════════════════════════════════════════════════
//  Fenwick Tree (Binary Indexed Tree)
// ═══════════════════════════════════════════════════════════════════════

/// Fenwick tree supporting O(log N) point-update, prefix-sum, and
/// rate-weighted random selection.  Used for efficient site selection
/// when lateral interactions make per-site rates heterogeneous.
struct FenwickTree {
    n: usize,
    tree: Vec<f64>,
    running_total: f64,
}

impl FenwickTree {
    fn new(n: usize) -> Self {
        Self { n, tree: vec![0.0; n + 1], running_total: 0.0 }
    }

    /// Add `delta` to position `i` (0-indexed). O(log N).
    #[inline]
    fn update(&mut self, i: usize, delta: f64) {
        self.running_total += delta;
        let mut j = i + 1;
        while j <= self.n {
            self.tree[j] += delta;
            j += j & j.wrapping_neg();
        }
    }

    /// Total sum of all elements. O(1).
    #[inline]
    fn total(&self) -> f64 {
        self.running_total
    }

    /// Find smallest 0-indexed position whose prefix sum >= target.
    /// O(log N) via binary lifting.
    #[inline]
    fn find(&self, mut target: f64) -> usize {
        let mut pos = 0usize;
        let mut bit = 1usize;
        while bit * 2 <= self.n { bit <<= 1; }
        while bit > 0 {
            let next = pos + bit;
            if next <= self.n && self.tree[next] < target {
                target -= self.tree[next];
                pos = next;
            }
            bit >>= 1;
        }
        pos.min(self.n.saturating_sub(1))
    }

    fn clear(&mut self) {
        self.tree.fill(0.0);
        self.running_total = 0.0;
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  Compiled Process
// ═══════════════════════════════════════════════════════════════════════

struct CompiledProcess {
    conditions: Vec<([i32; 3], usize)>,
    actions: Vec<([i32; 3], usize)>,
    tof_count: HashMap<String, f64>,
    site_type: Option<i32>,
}

/// Convert a variable-length offset Vec to fixed [i32; 3].
#[inline(always)]
fn to_offset3(v: &[i32]) -> [i32; 3] {
    match v.len() {
        0 => [0, 0, 0],
        1 => [v[0], 0, 0],
        2 => [v[0], v[1], 0],
        _ => [v[0], v[1], v[2]],
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  KMC Engine
// ═══════════════════════════════════════════════════════════════════════

pub struct KMCEngine {
    // Model
    project: Project,
    ndim: usize,
    lattice_size: [usize; 3],
    nsites: usize,
    nspecies: usize,
    nproc: usize,

    // Lattice state
    lattice: Vec<u16>,
    site_types: Vec<i32>,
    pub kmc_time: f64,
    pub kmc_step: u64,
    pub procstat: Vec<u64>,
    prev_procstat: Vec<u64>,
    prev_time: f64,

    // Base rate constants
    pub rates: Vec<f64>,
    accum_rates: Vec<f64>,

    // ── No-lateral: swap-with-last bookkeeping ──
    avail_sites: Vec<Vec<usize>>,
    // Direct-indexed: site_in_avail[proc][site] = index in avail_sites, or -1
    site_in_avail: Vec<Vec<i32>>,

    // ── Lateral: Fenwick tree bookkeeping ──
    rate_trees: Vec<FenwickTree>,
    site_rates: Vec<Vec<f64>>,  // [nproc][nsites], 0.0 = not available

    // Neighbor list
    neighbors: Vec<Vec<usize>>,

    // Flat lateral energy: lateral_energy[sp1 * nspecies + sp2]
    lateral_energy: Vec<f64>,
    has_lateral: bool,

    // BEP alpha per process (NAN = no BEP for this process)
    bep_alpha: Vec<f64>,

    // Cached thermal factor: eV / (kB * T).  Updated on set_parameter("T", ..).
    cached_beta_th: f64,

    // Periodic Fenwick tree rebuild counter (prevents floating-point drift)
    steps_since_rebuild: u64,

    empty_species: usize,
    compiled: Vec<CompiledProcess>,
    max_offset: usize,

    rng: rand::rngs::ThreadRng,
}

impl KMCEngine {
    pub fn new(mut project: Project, size: &[usize]) -> Self {
        let ndim = project.model_dimension;
        let mut lattice_size = [1usize; 3];
        for i in 0..ndim { lattice_size[i] = size[i]; }
        let nsites: usize = lattice_size[..ndim].iter().product();
        let nspecies = project.species_list.len();
        let nproc = project.process_list.len();

        project.rebuild_maps();

        let empty_species = project.species_list.iter()
            .find(|s| matches!(s.name.to_lowercase().as_str(), "empty" | "vacant" | "*"))
            .map(|s| s.id)
            .unwrap_or(0);

        // ── Flat lateral energy array ──
        let mut lateral_energy = vec![0.0f64; nspecies * nspecies];
        for li in &project.lateral_interactions {
            if let (Some(&sp1), Some(&sp2)) = (
                project.species_map.get(&li.species1),
                project.species_map.get(&li.species2),
            ) {
                lateral_energy[sp1 * nspecies + sp2] = li.energy;
                lateral_energy[sp2 * nspecies + sp1] = li.energy;
            }
        }
        let has_lateral = !project.lateral_interactions.is_empty();

        // ── BEP alpha array ──
        let bep_by_name: HashMap<String, f64> = project.bep_relations.iter()
            .map(|b| (b.process_name.clone(), b.alpha))
            .collect();
        let mut bep_alpha = vec![f64::NAN; nproc];
        for (pid, proc) in project.process_list.iter().enumerate() {
            if let Some(&alpha) = bep_by_name.get(&proc.name) {
                bep_alpha[pid] = alpha;
            }
        }

        // ── Compile processes to fixed-size offsets ──
        let compiled: Vec<CompiledProcess> = project.process_list.iter().map(|proc| {
            CompiledProcess {
                conditions: proc.conditions.iter().map(|c| {
                    (to_offset3(&c.offset), project.species_map[&c.species])
                }).collect(),
                actions: proc.actions.iter().map(|a| {
                    (to_offset3(&a.offset), project.species_map[&a.species])
                }).collect(),
                tof_count: proc.tof_count.clone(),
                site_type: proc.site_type,
            }
        }).collect();

        // ── Max offset ──
        let mut max_off: usize = 1;
        for cp in &compiled {
            for (off, _) in cp.conditions.iter().chain(cp.actions.iter()) {
                for d in 0..ndim {
                    max_off = max_off.max(off[d].unsigned_abs() as usize + 1);
                }
            }
        }

        // ── Neighbor list ──
        let neighbors = build_neighbor_list(ndim, &lattice_size, nsites);

        // ── Bookkeeping arrays ──
        let site_in_avail = vec![vec![-1i32; nsites]; nproc];
        let avail_sites = vec![Vec::new(); nproc];

        let rate_trees = if has_lateral {
            (0..nproc).map(|_| FenwickTree::new(nsites)).collect()
        } else {
            Vec::new()
        };
        let site_rates = if has_lateral {
            vec![vec![0.0f64; nsites]; nproc]
        } else {
            Vec::new()
        };

        let mut engine = Self {
            project,
            ndim,
            lattice_size,
            nsites,
            nspecies,
            nproc,
            lattice: vec![0u16; nsites],
            site_types: vec![0i32; nsites],
            kmc_time: 0.0,
            kmc_step: 0,
            procstat: vec![0u64; nproc],
            prev_procstat: vec![0u64; nproc],
            prev_time: 0.0,
            rates: vec![0.0; nproc],
            accum_rates: vec![0.0; nproc],
            avail_sites,
            site_in_avail,
            rate_trees,
            site_rates,
            neighbors,
            lateral_energy,
            has_lateral,
            bep_alpha,
            cached_beta_th: 0.0,
            steps_since_rebuild: 0,
            empty_species,
            compiled,
            max_offset: max_off,
            rng: rand::thread_rng(),
        };

        engine.update_cached_beta();
        engine.update_rate_constants();
        engine.rebuild_avail();
        engine
    }

    // ── Coordinate conversion (zero heap allocation) ────────────────

    #[inline(always)]
    fn site_to_coord(&self, site: usize) -> [i32; 3] {
        match self.ndim {
            1 => [site as i32, 0, 0],
            2 => [
                (site / self.lattice_size[1]) as i32,
                (site % self.lattice_size[1]) as i32,
                0,
            ],
            3 => {
                let lylz = self.lattice_size[1] * self.lattice_size[2];
                [
                    (site / lylz) as i32,
                    ((site % lylz) / self.lattice_size[2]) as i32,
                    (site % self.lattice_size[2]) as i32,
                ]
            }
            _ => unreachable!(),
        }
    }

    #[inline(always)]
    fn coord_to_site(&self, c: [i32; 3]) -> usize {
        match self.ndim {
            1 => c[0].rem_euclid(self.lattice_size[0] as i32) as usize,
            2 => {
                let x = c[0].rem_euclid(self.lattice_size[0] as i32) as usize;
                let y = c[1].rem_euclid(self.lattice_size[1] as i32) as usize;
                x * self.lattice_size[1] + y
            }
            3 => {
                let x = c[0].rem_euclid(self.lattice_size[0] as i32) as usize;
                let y = c[1].rem_euclid(self.lattice_size[1] as i32) as usize;
                let z = c[2].rem_euclid(self.lattice_size[2] as i32) as usize;
                x * self.lattice_size[1] * self.lattice_size[2]
                    + y * self.lattice_size[2] + z
            }
            _ => unreachable!(),
        }
    }

    /// Apply offset to coordinate and convert back to site index.
    #[inline(always)]
    fn offset_site(&self, coord: [i32; 3], offset: [i32; 3]) -> usize {
        self.coord_to_site([
            coord[0] + offset[0],
            coord[1] + offset[1],
            coord[2] + offset[2],
        ])
    }

    // ── Process availability check ──────────────────────────────────

    #[inline]
    fn check_process_at_site(&self, proc_id: usize, site: usize) -> bool {
        let cp = &self.compiled[proc_id];
        if let Some(req) = cp.site_type {
            if self.site_types[site] != req {
                return false;
            }
        }
        let coord = self.site_to_coord(site);
        for &(offset, sp_id) in &cp.conditions {
            if self.lattice[self.offset_site(coord, offset)] as usize != sp_id {
                return false;
            }
        }
        true
    }

    // ── Available sites bookkeeping ─────────────────────────────────

    fn add_to_avail(&mut self, p: usize, site: usize) {
        if self.site_in_avail[p][site] >= 0 { return; }

        if self.has_lateral {
            let rate = self.compute_site_rate(p, site);
            self.site_rates[p][site] = rate;
            self.rate_trees[p].update(site, rate);
            self.site_in_avail[p][site] = 0;  // flag: available
        } else {
            let idx = self.avail_sites[p].len();
            self.site_in_avail[p][site] = idx as i32;
            self.avail_sites[p].push(site);
        }
    }

    fn remove_from_avail(&mut self, p: usize, site: usize) {
        if self.site_in_avail[p][site] < 0 { return; }

        if self.has_lateral {
            let old_rate = self.site_rates[p][site];
            self.site_rates[p][site] = 0.0;
            self.rate_trees[p].update(site, -old_rate);
        } else {
            let idx = self.site_in_avail[p][site] as usize;
            let last_site = *self.avail_sites[p].last().unwrap();
            self.avail_sites[p][idx] = last_site;
            if last_site != site {
                self.site_in_avail[p][last_site] = idx as i32;
            }
            self.avail_sites[p].pop();
        }
        self.site_in_avail[p][site] = -1;
    }

    fn rebuild_avail(&mut self) {
        for p in 0..self.nproc {
            self.avail_sites[p].clear();
            self.site_in_avail[p].fill(-1);
        }
        if self.has_lateral {
            for p in 0..self.nproc {
                self.site_rates[p].fill(0.0);
                self.rate_trees[p].clear();
            }
        }
        for s in 0..self.nsites {
            for p in 0..self.nproc {
                if self.check_process_at_site(p, s) {
                    self.add_to_avail(p, s);
                }
            }
        }
    }

    // ── Affected sites (with neighbor-based fast path) ──────────────

    fn get_affected_sites(&self, site: usize, proc_id: usize) -> Vec<usize> {
        let cp = &self.compiled[proc_id];

        // Fast path: single-site process, max_offset <= 1
        // Use neighbor list directly instead of coordinate grid
        if cp.actions.len() == 1 && self.max_offset <= 1 {
            let action_site = if cp.actions[0].0 == [0, 0, 0] {
                site
            } else {
                let coord = self.site_to_coord(site);
                self.offset_site(coord, cp.actions[0].0)
            };
            let nbs = &self.neighbors[action_site];
            let mut affected = Vec::with_capacity(1 + nbs.len());
            affected.push(action_site);
            affected.extend_from_slice(nbs);
            return affected;
        }

        // General path: coordinate grid around all changed sites
        let coord = self.site_to_coord(site);
        let r = self.max_offset as i32;
        let mut affected = Vec::new();

        let mut changed: Vec<[i32; 3]> = Vec::new();
        for &(offset, _) in &cp.actions {
            changed.push([
                coord[0] + offset[0],
                coord[1] + offset[1],
                coord[2] + offset[2],
            ]);
        }

        match self.ndim {
            1 => {
                for cc in &changed {
                    for dx in -r..=r {
                        affected.push(self.coord_to_site([cc[0] + dx, 0, 0]));
                    }
                }
            }
            2 => {
                for cc in &changed {
                    for dx in -r..=r {
                        for dy in -r..=r {
                            affected.push(self.coord_to_site(
                                [cc[0] + dx, cc[1] + dy, 0]));
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
                                    [cc[0]+dx, cc[1]+dy, cc[2]+dz]));
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
                let was = self.site_in_avail[p][site] >= 0;
                let is = self.check_process_at_site(p, site);

                match (was, is) {
                    (false, true) => self.add_to_avail(p, site),
                    (true, false) => self.remove_from_avail(p, site),
                    (true, true) if self.has_lateral => {
                        // Rate may have changed due to neighbor state change
                        let old_rate = self.site_rates[p][site];
                        let new_rate = self.compute_site_rate(p, site);
                        if (new_rate - old_rate).abs() > 1e-30 {
                            self.site_rates[p][site] = new_rate;
                            self.rate_trees[p].update(site, new_rate - old_rate);
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    // ── Lateral interactions & BEP ──────────────────────────────────

    /// Update cached beta_thermal = eV / (kB * T).
    fn update_cached_beta(&mut self) {
        let t = self.project.get_parameter("T");
        self.cached_beta_th = if t > 0.0 { EV / (KB * t) } else { 0.0 };
    }

    /// Compute rate for (process, site) including lateral + BEP corrections.
    /// Uses cached beta_thermal to avoid repeated HashMap lookups.
    #[inline]
    fn compute_site_rate(&self, proc_id: usize, site: usize) -> f64 {
        let base_rate = self.rates[proc_id];
        if !self.has_lateral || self.cached_beta_th == 0.0 { return base_rate; }

        let beta_th = self.cached_beta_th;
        let e_react = self.interaction_energy(proc_id, site, false);

        let alpha = self.bep_alpha[proc_id];
        if !alpha.is_nan() {
            let e_prod = self.interaction_energy(proc_id, site, true);
            let delta_ea = alpha * (e_prod - e_react);
            base_rate * (-delta_ea * beta_th).exp()
        } else {
            base_rate * (e_react * beta_th).exp()
        }
    }

    /// Sum pairwise lateral interaction energy for condition/action sites.
    /// Uses flat array lookup instead of HashMap for speed.
    fn interaction_energy(&self, proc_id: usize, site: usize, is_product: bool) -> f64 {
        let coord = self.site_to_coord(site);
        let cp = &self.compiled[proc_id];
        let entries = if is_product { &cp.actions } else { &cp.conditions };
        let ns = self.nspecies;

        // Collect process sites (typically 1-3, allocated on stack via SmallVec pattern)
        let n_entries = entries.len();
        let mut proc_site_buf = [0usize; 8];  // stack buffer, covers up to 8-site processes
        for (i, &(off, _)) in entries.iter().enumerate() {
            proc_site_buf[i] = self.offset_site(coord, off);
        }
        let proc_sites = &proc_site_buf[..n_entries];

        let mut e_total = 0.0;
        for (i, &(_off, sp_id)) in entries.iter().enumerate() {
            if sp_id == self.empty_species { continue; }
            let entry_site = proc_sites[i];
            for &nn in &self.neighbors[entry_site] {
                // Skip intra-process pairs (only check if multi-site process)
                if n_entries > 1 && proc_sites.contains(&nn) { continue; }
                let nn_sp = self.lattice[nn] as usize;
                // Flat array lookup: O(1), no hashing
                e_total += self.lateral_energy[sp_id * ns + nn_sp];
            }
        }
        e_total
    }

    // ── Rates ───────────────────────────────────────────────────────

    pub fn update_rate_constants(&mut self) {
        for i in 0..self.nproc {
            self.rates[i] = evaluate_rate_expression(
                &self.project.process_list[i].rate_constant,
                &self.project.parameter_list,
            );
        }
        if self.has_lateral && !self.rate_trees.is_empty() {
            self.rebuild_lateral_rates();
        }
    }

    fn rebuild_lateral_rates(&mut self) {
        for p in 0..self.nproc {
            self.rate_trees[p].clear();
            for s in 0..self.nsites {
                if self.site_in_avail[p][s] >= 0 {
                    let rate = self.compute_site_rate(p, s);
                    self.site_rates[p][s] = rate;
                    self.rate_trees[p].update(s, rate);
                } else {
                    self.site_rates[p][s] = 0.0;
                }
            }
        }
    }

    fn update_accum_rates(&mut self) -> f64 {
        let mut total = 0.0;
        for p in 0..self.nproc {
            total += if self.has_lateral {
                self.rate_trees[p].total()
            } else {
                self.rates[p] * self.avail_sites[p].len() as f64
            };
            self.accum_rates[p] = total;
        }
        total
    }

    // ── KMC Step (BKL algorithm) ────────────────────────────────────

    /// Execute one BKL KMC step. Returns false if system is frozen.
    pub fn do_kmc_step(&mut self) -> bool {
        let total_rate = self.update_accum_rates();
        if total_rate <= 0.0 { return false; }

        // Clamp to avoid ln(0) and degenerate edge cases
        let r_time: f64 = self.rng.gen::<f64>().max(f64::MIN_POSITIVE);
        let r_proc: f64 = self.rng.gen::<f64>().max(f64::MIN_POSITIVE);
        let r_site: f64 = self.rng.gen::<f64>().max(f64::MIN_POSITIVE);

        self.kmc_time += -r_time.ln() / total_rate;

        // Process selection: binary search O(log N_proc)
        let target = r_proc * total_rate;
        let proc_id = match self.accum_rates[..self.nproc]
            .binary_search_by(|v| v.partial_cmp(&target).unwrap())
        {
            Ok(i) => i,
            Err(i) => i.min(self.nproc - 1),
        };

        // Site selection
        let site = if self.has_lateral {
            // Fenwick tree: O(log N_sites)
            let proc_total = self.rate_trees[proc_id].total();
            if proc_total <= 0.0 { return false; }
            self.rate_trees[proc_id].find(r_site * proc_total)
        } else {
            // Uniform random: O(1)
            let n = self.avail_sites[proc_id].len();
            if n == 0 { return false; }
            let idx = ((r_site * n as f64) as usize).min(n - 1);
            self.avail_sites[proc_id][idx]
        };

        // Execute: update lattice
        let coord = self.site_to_coord(site);
        let cp = &self.compiled[proc_id];
        for &(offset, new_sp) in &cp.actions {
            let ns = self.offset_site(coord, offset);
            self.lattice[ns] = new_sp as u16;
        }

        // Update bookkeeping
        let affected = self.get_affected_sites(site, proc_id);
        self.update_avail_after_exec(&affected);

        self.kmc_step += 1;
        self.procstat[proc_id] += 1;

        // Periodic Fenwick tree rebuild to prevent floating-point drift
        if self.has_lateral {
            self.steps_since_rebuild += 1;
            if self.steps_since_rebuild >= 100_000 {
                self.rebuild_lateral_rates();
                self.steps_since_rebuild = 0;
            }
        }

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

    // ── Observables ─────────────────────────────────────────────────

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
        if dt <= 0.0 { return HashMap::new(); }
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

    // ── Parameters ──────────────────────────────────────────────────

    pub fn set_parameter(&mut self, name: &str, value: f64) {
        self.project.set_parameter(name, value);
        if name == "T" {
            self.update_cached_beta();
        }
        self.update_rate_constants();
    }

    pub fn get_parameter(&self, name: &str) -> f64 {
        self.project.get_parameter(name)
    }

    // ── Site types ──────────────────────────────────────────────────

    pub fn set_site_type(&mut self, site: usize, stype: i32) {
        self.site_types[site] = stype;
    }

    pub fn set_site_types_region<F: Fn(&[i32]) -> bool>(&mut self, f: F, stype: i32) {
        for s in 0..self.nsites {
            let coord = self.site_to_coord(s);
            if f(&coord[..self.ndim]) {
                self.site_types[s] = stype;
            }
        }
        self.rebuild_avail();
    }

    // ── State management ────────────────────────────────────────────

    pub fn reset(&mut self) {
        self.lattice.fill(0);
        self.site_types.fill(0);
        self.kmc_time = 0.0;
        self.kmc_step = 0;
        self.procstat.fill(0);
        self.prev_procstat.fill(0);
        self.prev_time = 0.0;
        self.rebuild_avail();
    }

    // ── Printing ────────────────────────────────────────────────────

    pub fn print_rates(&self) {
        println!("\nRate constants:");
        println!("  {:<35} {:>12}  {:>8}", "Process", "k [s^-1]", "N_avail");
        println!("  {} {}  {}", "-".repeat(35), "-".repeat(12), "-".repeat(8));
        for i in 0..self.nproc {
            let n_avail = if self.has_lateral {
                self.site_in_avail[i].iter().filter(|&&v| v >= 0).count()
            } else {
                self.avail_sites[i].len()
            };
            println!("  {:<35} {:>12.4e}  {:>8}",
                     self.project.process_list[i].name,
                     self.rates[i], n_avail);
        }
        if self.has_lateral {
            let n_lat = self.lateral_energy.iter().filter(|&&v| v != 0.0).count() / 2;
            let n_bep = self.bep_alpha.iter().filter(|v| !v.is_nan()).count();
            println!("  Lateral interactions: {}, BEP: {}", n_lat, n_bep);
            if !self.neighbors.is_empty() {
                println!("  Neighbors per site: {}", self.neighbors[0].len());
            }
            println!("  Site selection: Fenwick tree O(log N)");
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

// ═══════════════════════════════════════════════════════════════════════
//  Neighbor list construction
// ═══════════════════════════════════════════════════════════════════════

fn build_neighbor_list(ndim: usize, ls: &[usize; 3], nsites: usize) -> Vec<Vec<usize>> {
    let nn_offsets: &[[i32; 3]] = match ndim {
        1 => &[[-1,0,0], [1,0,0]],
        2 => &[[-1,0,0], [1,0,0], [0,-1,0], [0,1,0]],
        3 => &[[-1,0,0], [1,0,0], [0,-1,0], [0,1,0], [0,0,-1], [0,0,1]],
        _ => &[],
    };

    let mut neighbors = vec![Vec::with_capacity(nn_offsets.len()); nsites];
    for s in 0..nsites {
        let coord: [i32; 3] = match ndim {
            1 => [s as i32, 0, 0],
            2 => [(s / ls[1]) as i32, (s % ls[1]) as i32, 0],
            3 => {
                let lylz = ls[1] * ls[2];
                [(s / lylz) as i32, ((s % lylz) / ls[2]) as i32, (s % ls[2]) as i32]
            }
            _ => [0, 0, 0],
        };
        for off in nn_offsets {
            let nc = [coord[0]+off[0], coord[1]+off[1], coord[2]+off[2]];
            let ns = match ndim {
                1 => nc[0].rem_euclid(ls[0] as i32) as usize,
                2 => {
                    let x = nc[0].rem_euclid(ls[0] as i32) as usize;
                    let y = nc[1].rem_euclid(ls[1] as i32) as usize;
                    x * ls[1] + y
                }
                3 => {
                    let x = nc[0].rem_euclid(ls[0] as i32) as usize;
                    let y = nc[1].rem_euclid(ls[1] as i32) as usize;
                    let z = nc[2].rem_euclid(ls[2] as i32) as usize;
                    x * ls[1] * ls[2] + y * ls[2] + z
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

// ═══════════════════════════════════════════════════════════════════════
//  Tests
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fenwick_basic() {
        let mut ft = FenwickTree::new(4);
        ft.update(0, 3.0);
        ft.update(1, 4.0);
        ft.update(2, 1.0);
        ft.update(3, 7.0);
        assert!((ft.total() - 15.0).abs() < 1e-10);
        // find: smallest i such that prefix_sum(0..=i) >= target
        assert_eq!(ft.find(3.0), 0);   // prefix[0]=3 >= 3
        assert_eq!(ft.find(3.5), 1);   // prefix[0]=3 < 3.5, prefix[1]=7 >= 3.5
        assert_eq!(ft.find(7.0), 1);   // prefix[1]=7 >= 7
        assert_eq!(ft.find(7.5), 2);   // prefix[2]=8 >= 7.5
        assert_eq!(ft.find(8.0), 2);   // prefix[2]=8 >= 8
        assert_eq!(ft.find(8.5), 3);   // prefix[3]=15 >= 8.5
    }

    #[test]
    fn test_fenwick_update_remove() {
        let mut ft = FenwickTree::new(3);
        ft.update(0, 5.0);
        ft.update(1, 3.0);
        ft.update(2, 2.0);
        assert!((ft.total() - 10.0).abs() < 1e-10);

        // Remove element 1 (set to 0)
        ft.update(1, -3.0);
        assert!((ft.total() - 7.0).abs() < 1e-10);
        assert_eq!(ft.find(5.0), 0);   // prefix[0]=5 >= 5
        assert_eq!(ft.find(5.5), 2);   // prefix[0]=5 < 5.5, prefix[1]=5 < 5.5, prefix[2]=7 >= 5.5
    }

    #[test]
    fn test_fenwick_clear() {
        let mut ft = FenwickTree::new(4);
        ft.update(0, 10.0);
        ft.update(3, 5.0);
        assert!((ft.total() - 15.0).abs() < 1e-10);
        ft.clear();
        assert!((ft.total()).abs() < 1e-10);
    }

    #[test]
    fn test_to_offset3() {
        assert_eq!(to_offset3(&[]), [0, 0, 0]);
        assert_eq!(to_offset3(&[1]), [1, 0, 0]);
        assert_eq!(to_offset3(&[1, -1]), [1, -1, 0]);
        assert_eq!(to_offset3(&[1, -1, 2]), [1, -1, 2]);
    }
}
