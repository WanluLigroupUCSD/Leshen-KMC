"""
Core Kinetic Monte Carlo engine implementing the BKL rejection-free algorithm.

Features:
  - BKL/VSSM rejection-free algorithm
  - Neighbor list for spatial events
  - Pairwise lateral interactions with per-site rates
  - Surface diffusion support
  - BEP (Bronsted-Evans-Polanyi) relations
  - Site type support

Algorithm per step:
  1. Compute cumulative rates (per-site when lateral interactions exist)
  2. Draw random numbers r1, r2, r3 ~ U(0,1)
  3. Advance time: dt = -ln(r1) / R_total
  4. Select process by binary search on cumulative rates using r2
  5. Select site (uniform or rate-weighted) using r3
  6. Execute process (update lattice + bookkeeping + affected rates)
"""

import numpy as np
import time as _time
from collections import defaultdict
from .rates import evaluate_rate_expression
from .units import kB, eV


class ParameterProxy:
    """
    Allows attribute-style access to model parameters.

    Usage:
        model.parameters.T = 550
        print(model.parameters.T)
    """

    def __init__(self, param_list, engine=None):
        object.__setattr__(self, '_params', {p.name: p for p in param_list})
        object.__setattr__(self, '_engine', engine)

    def __getattr__(self, name):
        if name.startswith('_'):
            return object.__getattribute__(self, name)
        params = object.__getattribute__(self, '_params')
        if name in params:
            return params[name].value
        raise AttributeError(f"No parameter '{name}'")

    def __setattr__(self, name, value):
        if name.startswith('_'):
            object.__setattr__(self, name, value)
            return
        params = object.__getattribute__(self, '_params')
        if name in params:
            params[name].value = value
            engine = object.__getattribute__(self, '_engine')
            if engine is not None:
                engine._update_rate_constants()
        else:
            raise AttributeError(f"No parameter '{name}'")

    def __repr__(self):
        params = object.__getattribute__(self, '_params')
        lines = ["Parameters:"]
        for name, p in params.items():
            lines.append(f"  {name} = {p.value}")
        return '\n'.join(lines)

    def as_dict(self):
        params = object.__getattribute__(self, '_params')
        return {name: p.value for name, p in params.items()}


class KMCEngine:
    """
    Lattice Kinetic Monte Carlo simulation engine.

    Parameters
    ----------
    project : Project
        The model definition.
    size : list of int
        Lattice dimensions, e.g. [50, 50] for a 50x50 2D lattice.
    print_rates : bool
        Print rate constants on initialization.
    banner : bool
        Print welcome banner.
    """

    def __init__(self, project, size=None, print_rates=True, banner=True):
        self.project = project
        self.size = size or [20, 20]
        self.ndim = project.meta.get('model_dimension', 2)

        if banner:
            print(f"mykmc KMC Engine - Model: "
                  f"{project.meta.get('model_name', 'unnamed')}")
            print(f"  Lattice: {self.size}, Dimension: {self.ndim}")

        # Species encoding
        self.nspecies = len(project.species_list)
        self.species_names = [sp.name for sp in project.species_list]
        self.species_id = {sp.name: sp.id for sp in project.species_list}

        # Detect empty species
        self._empty_species = 0
        for sp in project.species_list:
            if sp.name.lower() in ('empty', 'vacant', '*'):
                self._empty_species = sp.id
                break

        # Parameters (with auto-update callback)
        self.parameters = ParameterProxy(project.parameter_list, engine=self)

        # Initialize lattice (includes site types)
        self._init_lattice()

        # Build neighbor list
        self._build_neighbor_list()

        # Setup lateral interactions
        self._setup_lateral_interactions()

        # Setup BEP relations
        self._setup_bep_relations()

        # Initialize processes
        self._init_processes()

        # Compute max offset range for neighbor updates
        self._compute_max_offset()

        # Rate constants
        self.rates = np.zeros(self.nproc)
        self._update_rate_constants()

        # Per-site rate tracking (only when lateral interactions exist)
        if self._has_lateral:
            self._proc_total_rates = np.zeros(self.nproc)
            self._rebuild_per_site_rates()

        # Cumulative rates array
        self.accum_rates = np.zeros(self.nproc)

        # Simulation state
        self.kmc_time = 0.0
        self.kmc_step = 0
        self.procstat = np.zeros(self.nproc, dtype=np.int64)
        self._prev_procstat = np.zeros(self.nproc, dtype=np.int64)
        self._prev_time = 0.0

        # Wall time
        self._start_walltime = _time.time()

        if banner and self._has_lateral:
            n_lat = len(project.lateral_interactions)
            n_bep = len(self._bep_by_name)
            print(f"  Lateral interactions: {n_lat}, BEP relations: {n_bep}")
            print(f"  Neighbors per site: {len(self.neighbors[0])}")

        if print_rates:
            self.print_rates()

    # ----------------------------------------------------------------
    # Initialization
    # ----------------------------------------------------------------

    def _init_lattice(self):
        """Initialize the lattice array with default species and site types."""
        self.lattice_size = tuple(self.size[:self.ndim])
        self.nsites = int(np.prod(self.lattice_size))

        # Determine default species
        default_sp = 0
        if self.project.lattice.layers:
            layer = self.project.lattice.layers[0]
            if layer.sites:
                ds_name = layer.sites[0].default_species
                if ds_name in self.species_id:
                    default_sp = self.species_id[ds_name]

        self.lattice = np.full(self.nsites, default_sp, dtype=np.int32)

        # Site types array
        self.site_types = np.zeros(self.nsites, dtype=np.int32)
        if self.project.lattice.layers:
            layer = self.project.lattice.layers[0]
            if layer.sites and layer.sites[0].site_type is not None:
                stype = layer.sites[0].site_type
                if isinstance(stype, int):
                    self.site_types[:] = stype

    def _build_neighbor_list(self):
        """Build nearest-neighbor list for the lattice."""
        self.neighbors = [[] for _ in range(self.nsites)]
        nn_offsets = self._get_nn_offsets()

        for s in range(self.nsites):
            coord = self._site_to_coord(s)
            for offset in nn_offsets:
                nc = tuple(c + o for c, o in zip(coord, offset))
                ns = self._coord_to_site(nc)
                if ns != s:
                    self.neighbors[s].append(ns)

    def _get_nn_offsets(self):
        """Get nearest-neighbor offsets for the lattice geometry."""
        if self.ndim == 1:
            return [(-1,), (1,)]
        elif self.ndim == 2:
            return [(-1, 0), (1, 0), (0, -1), (0, 1)]
        else:
            return [(-1, 0, 0), (1, 0, 0), (0, -1, 0),
                    (0, 1, 0), (0, 0, -1), (0, 0, 1)]

    def _setup_lateral_interactions(self):
        """Build lateral interaction lookup dict from project definition."""
        self._lateral_dict = {}
        if hasattr(self.project, 'lateral_interactions'):
            for li in self.project.lateral_interactions:
                sp1 = self.species_id.get(li.species1)
                sp2 = self.species_id.get(li.species2)
                if sp1 is not None and sp2 is not None:
                    self._lateral_dict[(sp1, sp2)] = li.energy
                    # Symmetric: A-B same as B-A
                    self._lateral_dict[(sp2, sp1)] = li.energy
        self._has_lateral = bool(self._lateral_dict)

    def _setup_bep_relations(self):
        """Build BEP relation lookup from project definition."""
        self._bep_by_name = {}
        if hasattr(self.project, 'bep_relations'):
            self._bep_by_name = dict(self.project.bep_relations)

    def _init_processes(self):
        """Convert process definitions to numeric arrays for fast lookup."""
        self.nproc = len(self.project.process_list)
        self.process_names = [p.name for p in self.project.process_list]

        self._proc_conditions = []
        self._proc_actions = []
        self._proc_rate_exprs = []
        self._proc_tof_count = []
        self._proc_site_types = []

        for proc in self.project.process_list:
            conds = []
            for c in proc.conditions:
                offset = c.coord.offset[:self.ndim]
                sp_id = self.species_id[c.species]
                conds.append((offset, sp_id))
            self._proc_conditions.append(conds)

            acts = []
            for a in proc.actions:
                offset = a.coord.offset[:self.ndim]
                sp_id = self.species_id[a.species]
                acts.append((offset, sp_id))
            self._proc_actions.append(acts)

            self._proc_rate_exprs.append(proc.rate_constant)
            self._proc_tof_count.append(proc.tof_count)
            self._proc_site_types.append(getattr(proc, 'site_type', None))

        # Resolve BEP process IDs
        self._bep_proc_ids = {}
        for pid, name in enumerate(self.process_names):
            if name in self._bep_by_name:
                self._bep_proc_ids[pid] = self._bep_by_name[name]

        # Available sites bookkeeping: O(1) add/remove via swap-with-last
        self._avail_sites = [[] for _ in range(self.nproc)]
        self._site_in_avail = [dict() for _ in range(self.nproc)]
        # Parallel rate array for per-site rates (only used with lateral)
        self._avail_rates = [[] for _ in range(self.nproc)]

        self._rebuild_avail_sites()

    def _compute_max_offset(self):
        """Find the maximum offset range across all process conditions."""
        max_r = 1
        for conds in self._proc_conditions:
            for offset, _ in conds:
                for o in offset:
                    max_r = max(max_r, abs(o) + 1)
        for acts in self._proc_actions:
            for offset, _ in acts:
                for o in offset:
                    max_r = max(max_r, abs(o) + 1)
        self._max_offset = max_r

    # ----------------------------------------------------------------
    # Coordinate conversion (with periodic boundary conditions)
    # ----------------------------------------------------------------

    def _site_to_coord(self, site):
        if self.ndim == 1:
            return (site,)
        elif self.ndim == 2:
            return (site // self.lattice_size[1],
                    site % self.lattice_size[1])
        else:
            Ly = self.lattice_size[1]
            Lz = self.lattice_size[2]
            return (site // (Ly * Lz),
                    (site % (Ly * Lz)) // Lz,
                    site % Lz)

    def _coord_to_site(self, coord):
        if self.ndim == 1:
            return coord[0] % self.lattice_size[0]
        elif self.ndim == 2:
            x = coord[0] % self.lattice_size[0]
            y = coord[1] % self.lattice_size[1]
            return x * self.lattice_size[1] + y
        else:
            x = coord[0] % self.lattice_size[0]
            y = coord[1] % self.lattice_size[1]
            z = coord[2] % self.lattice_size[2]
            return (x * self.lattice_size[1] * self.lattice_size[2]
                    + y * self.lattice_size[2] + z)

    # ----------------------------------------------------------------
    # Available sites bookkeeping
    # ----------------------------------------------------------------

    def _check_process_at_site(self, proc_id, site):
        """Check if process can occur at site (species + site type)."""
        # Check site type requirement
        site_type_req = self._proc_site_types[proc_id]
        if site_type_req is not None:
            if self.site_types[site] != site_type_req:
                return False

        # Check species conditions
        coord = self._site_to_coord(site)
        for offset, sp_id in self._proc_conditions[proc_id]:
            neighbor = tuple(c + o for c, o in zip(coord, offset))
            if self.lattice[self._coord_to_site(neighbor)] != sp_id:
                return False
        return True

    def _add_to_avail(self, proc_id, site):
        if site not in self._site_in_avail[proc_id]:
            idx = len(self._avail_sites[proc_id])
            self._site_in_avail[proc_id][site] = idx
            self._avail_sites[proc_id].append(site)
            if self._has_lateral and hasattr(self, '_proc_total_rates'):
                rate = self._compute_site_rate(proc_id, site)
                self._avail_rates[proc_id].append(rate)
                self._proc_total_rates[proc_id] += rate
            else:
                self._avail_rates[proc_id].append(0.0)

    def _remove_from_avail(self, proc_id, site):
        idx_map = self._site_in_avail[proc_id]
        if site in idx_map:
            idx = idx_map[site]
            avail = self._avail_sites[proc_id]
            rates = self._avail_rates[proc_id]

            if self._has_lateral and hasattr(self, '_proc_total_rates'):
                self._proc_total_rates[proc_id] -= rates[idx]

            # Swap with last
            last_site = avail[-1]
            avail[idx] = last_site
            rates[idx] = rates[-1]
            idx_map[last_site] = idx
            avail.pop()
            rates.pop()
            del idx_map[site]

    def _rebuild_avail_sites(self):
        """Full rebuild of available sites (used at initialization)."""
        for p in range(self.nproc):
            self._avail_sites[p] = []
            self._site_in_avail[p] = {}
            self._avail_rates[p] = []
        for s in range(self.nsites):
            for p in range(self.nproc):
                if self._check_process_at_site(p, s):
                    self._add_to_avail(p, s)

    def _rebuild_per_site_rates(self):
        """Recompute all per-site rates from scratch."""
        for p in range(self.nproc):
            total = 0.0
            for i, site in enumerate(self._avail_sites[p]):
                rate = self._compute_site_rate(p, site)
                self._avail_rates[p][i] = rate
                total += rate
            self._proc_total_rates[p] = total

    def _get_affected_sites(self, site, actions):
        """
        Get all sites whose process availability might change
        after executing actions rooted at site.
        """
        coord = self._site_to_coord(site)
        affected = set()
        r = self._max_offset

        # Sites directly changed by actions
        changed_coords = []
        for offset, _ in actions:
            nc = tuple(c + o for c, o in zip(coord, offset))
            changed_coords.append(nc)

        # Expand to include neighbors within max_offset range
        if self.ndim == 1:
            for cc in changed_coords:
                for dx in range(-r, r + 1):
                    affected.add(self._coord_to_site((cc[0] + dx,)))
        elif self.ndim == 2:
            for cc in changed_coords:
                for dx in range(-r, r + 1):
                    for dy in range(-r, r + 1):
                        affected.add(self._coord_to_site(
                            (cc[0] + dx, cc[1] + dy)))
        else:
            for cc in changed_coords:
                for dx in range(-r, r + 1):
                    for dy in range(-r, r + 1):
                        for dz in range(-r, r + 1):
                            affected.add(self._coord_to_site(
                                (cc[0] + dx, cc[1] + dy, cc[2] + dz)))

        return affected

    def _update_avail_after_execution(self, affected_sites):
        """Update available sites and per-site rates for all affected sites."""
        for site in affected_sites:
            for p in range(self.nproc):
                was_avail = site in self._site_in_avail[p]
                is_avail = self._check_process_at_site(p, site)

                if is_avail and not was_avail:
                    self._add_to_avail(p, site)
                elif not is_avail and was_avail:
                    self._remove_from_avail(p, site)
                elif is_avail and was_avail and self._has_lateral:
                    idx = self._site_in_avail[p][site]
                    old_rate = self._avail_rates[p][idx]
                    new_rate = self._compute_site_rate(p, site)
                    self._avail_rates[p][idx] = new_rate
                    self._proc_total_rates[p] += (new_rate - old_rate)

    # ----------------------------------------------------------------
    # Lateral interactions & BEP
    # ----------------------------------------------------------------

    def _get_temperature(self):
        """Get current temperature from parameters."""
        try:
            return float(self.parameters.T)
        except (AttributeError, TypeError):
            return 300.0

    def _compute_site_rate(self, proc_id, site):
        """
        Compute rate for a specific (process, site) pair including
        lateral interactions and BEP corrections.

        Without lateral interactions, returns the base rate.
        With lateral interactions:
          - Default: k = k_base * exp(+E_lat_react / (kB*T))
            where E_lat_react is the sum of pairwise interactions
            for adsorbates in the reactant state.
          - With BEP: k = k_base * exp(-alpha * delta_delta_H / (kB*T))
            where delta_delta_H = E_lat_product - E_lat_reactant.
        """
        base_rate = self.rates[proc_id]

        if not self._has_lateral:
            return base_rate

        T = self._get_temperature()
        if T <= 0:
            return base_rate

        beta_th = eV / (kB * T)

        # Compute interaction energy for reactant state
        E_lat_react = self._compute_interaction_energy(proc_id, site,
                                                       is_product=False)

        # BEP correction
        if proc_id in self._bep_proc_ids:
            bep = self._bep_proc_ids[proc_id]
            E_lat_prod = self._compute_interaction_energy(proc_id, site,
                                                          is_product=True)
            delta_delta_H = E_lat_prod - E_lat_react
            delta_Ea = bep.alpha * delta_delta_H
            return base_rate * np.exp(-delta_Ea * beta_th)

        # Default: reactant-state interaction modifies barrier
        # Positive E_lat (repulsion) → destabilized adsorbate → higher rate
        return base_rate * np.exp(E_lat_react * beta_th)

    def _compute_interaction_energy(self, proc_id, site, is_product=False):
        """
        Compute total pairwise lateral interaction energy for a process.

        Sums epsilon(species_i, species_neighbor) for all condition/action
        sites and their nearest neighbors (excluding intra-process sites).
        """
        coord = self._site_to_coord(site)
        entries = (self._proc_actions[proc_id] if is_product
                   else self._proc_conditions[proc_id])

        # Collect sites involved in this process (to exclude from counting)
        proc_sites = set()
        for offset, _ in entries:
            ps = self._coord_to_site(
                tuple(c + o for c, o in zip(coord, offset)))
            proc_sites.add(ps)

        E_total = 0.0
        for offset, sp_id in entries:
            if sp_id == self._empty_species:
                continue  # Empty sites have no lateral interactions

            entry_site = self._coord_to_site(
                tuple(c + o for c, o in zip(coord, offset)))

            for nn in self.neighbors[entry_site]:
                if nn in proc_sites:
                    continue  # Skip intra-process neighbor

                nn_sp = self.lattice[nn]
                pair = (sp_id, nn_sp)
                if pair in self._lateral_dict:
                    E_total += self._lateral_dict[pair]

        return E_total

    # ----------------------------------------------------------------
    # Rate constants
    # ----------------------------------------------------------------

    def _update_rate_constants(self):
        """Evaluate all rate expressions with current parameters."""
        params = {p.name: p.value for p in self.project.parameter_list}
        for i, expr in enumerate(self._proc_rate_exprs):
            self.rates[i] = evaluate_rate_expression(expr, params)
        # Recompute per-site rates if lateral interactions are active
        if self._has_lateral and hasattr(self, '_proc_total_rates'):
            self._rebuild_per_site_rates()

    def _update_accum_rates(self):
        """Build cumulative rate array for process selection."""
        total = 0.0
        for p in range(self.nproc):
            if self._has_lateral:
                total += self._proc_total_rates[p]
            else:
                total += self.rates[p] * len(self._avail_sites[p])
            self.accum_rates[p] = total
        return total

    # ----------------------------------------------------------------
    # KMC step (BKL algorithm)
    # ----------------------------------------------------------------

    def do_kmc_step(self):
        """
        Execute one KMC step using the BKL rejection-free algorithm.

        Returns True if a step was executed, False if the system is frozen.
        """
        total_rate = self._update_accum_rates()

        if total_rate <= 0.0:
            return False

        # Three independent random numbers
        r_time = np.random.random()
        r_proc = np.random.random()
        r_site = np.random.random()

        # Time advancement (Poisson process)
        self.kmc_time += -np.log(r_time) / total_rate

        # Process selection (binary search on cumulative rates)
        proc_id = int(np.searchsorted(self.accum_rates,
                                       r_proc * total_rate))
        if proc_id >= self.nproc:
            proc_id = self.nproc - 1

        # Site selection
        n_avail = len(self._avail_sites[proc_id])
        if n_avail == 0:
            return False

        if self._has_lateral:
            # Select site proportional to per-site rate using cumulative sum
            rates_list = self._avail_rates[proc_id]
            total_proc = self._proc_total_rates[proc_id]
            if total_proc <= 0:
                return False
            target = r_site * total_proc
            cumul = 0.0
            site_idx = n_avail - 1
            for k in range(n_avail):
                cumul += rates_list[k]
                if cumul >= target:
                    site_idx = k
                    break
        else:
            # Uniform selection among available sites
            site_idx = min(int(r_site * n_avail), n_avail - 1)

        site = self._avail_sites[proc_id][site_idx]

        # Execute: update lattice
        coord = self._site_to_coord(site)
        for offset, new_sp in self._proc_actions[proc_id]:
            neighbor = tuple(c + o for c, o in zip(coord, offset))
            self.lattice[self._coord_to_site(neighbor)] = new_sp

        # Update bookkeeping
        affected = self._get_affected_sites(
            site, self._proc_actions[proc_id])
        self._update_avail_after_execution(affected)

        # Statistics
        self.kmc_step += 1
        self.procstat[proc_id] += 1

        return True

    def do_steps(self, n, progress=False):
        """
        Execute n KMC steps.

        Parameters
        ----------
        n : int
            Number of steps to execute.
        progress : bool
            Print progress at 10% intervals.
        """
        n = int(n)
        report_interval = max(1, n // 10)
        for i in range(n):
            if not self.do_kmc_step():
                print(f"System frozen at step {self.kmc_step}, "
                      f"time={self.kmc_time:.6e} s")
                break
            if progress and (i + 1) % report_interval == 0:
                pct = 100 * (i + 1) / n
                print(f"  [{pct:5.1f}%] step={self.kmc_step}, "
                      f"time={self.kmc_time:.6e} s")

    # ----------------------------------------------------------------
    # Observables
    # ----------------------------------------------------------------

    def get_coverage(self):
        """Get fractional coverage for each species."""
        coverage = {}
        for sp in self.project.species_list:
            coverage[sp.name] = float(np.sum(self.lattice == sp.id)) / \
                self.nsites
        return coverage

    def get_occupation(self):
        """Get occupation matrix [nspecies, nsites]."""
        occ = np.zeros((self.nspecies, self.nsites))
        for s in range(self.nsites):
            occ[self.lattice[s], s] = 1.0
        return occ

    def get_tof(self):
        """
        Get turn-over frequencies since last call.

        Returns dict of {observable_name: TOF in s^-1 per site}.
        """
        dt = self.kmc_time - self._prev_time
        if dt <= 0:
            return {}

        tof = defaultdict(float)
        for p in range(self.nproc):
            delta = self.procstat[p] - self._prev_procstat[p]
            for obs, coeff in self._proc_tof_count[p].items():
                tof[obs] += coeff * delta / (dt * self.nsites)

        self._prev_procstat = self.procstat.copy()
        self._prev_time = self.kmc_time
        return dict(tof)

    def get_process_stats(self):
        """Get execution count per process."""
        return {self.process_names[i]: int(self.procstat[i])
                for i in range(self.nproc)}

    def get_neighbor_coverages(self, site):
        """
        Get species counts among nearest neighbors of a site.

        Returns dict of {species_name: count}.
        """
        counts = defaultdict(int)
        for nn in self.neighbors[site]:
            sp = self.species_names[self.lattice[nn]]
            counts[sp] += 1
        return dict(counts)

    def get_lattice_2d(self):
        """
        Return 2D array representation of the lattice (for visualization).
        Only works for 2D lattices.
        """
        if self.ndim != 2:
            raise ValueError("get_lattice_2d only works for 2D lattices")
        return self.lattice.reshape(self.lattice_size)

    # ----------------------------------------------------------------
    # Site manipulation
    # ----------------------------------------------------------------

    def put(self, site_coord, species_name):
        """Set species at a specific lattice coordinate."""
        site = self._coord_to_site(tuple(site_coord[:self.ndim]))
        sp_id = self.species_id[species_name]
        self.lattice[site] = sp_id
        affected = self._get_affected_sites(
            site, [((0,) * self.ndim, sp_id)])
        self._update_avail_after_execution(affected)

    def get(self, site_coord):
        """Get species name at a specific lattice coordinate."""
        site = self._coord_to_site(tuple(site_coord[:self.ndim]))
        return self.species_names[self.lattice[site]]

    def set_site_type(self, site_coord, site_type):
        """Set the site type for a specific lattice site."""
        site = self._coord_to_site(tuple(site_coord[:self.ndim]))
        self.site_types[site] = site_type

    def set_site_types_region(self, region_func, site_type):
        """
        Set site type for all sites where region_func(coord) returns True.

        Parameters
        ----------
        region_func : callable
            Function taking a coordinate tuple and returning bool.
        site_type : int
            Site type to assign.
        """
        for s in range(self.nsites):
            coord = self._site_to_coord(s)
            if region_func(coord):
                self.site_types[s] = site_type
        self._rebuild_avail_sites()
        if self._has_lateral and hasattr(self, '_proc_total_rates'):
            self._rebuild_per_site_rates()

    # ----------------------------------------------------------------
    # State management
    # ----------------------------------------------------------------

    def reset(self):
        """Reset simulation to initial state (all sites = default species)."""
        default_sp = 0
        if self.project.lattice.layers:
            layer = self.project.lattice.layers[0]
            if layer.sites:
                ds_name = layer.sites[0].default_species
                if ds_name in self.species_id:
                    default_sp = self.species_id[ds_name]

        self.lattice[:] = default_sp
        self.kmc_time = 0.0
        self.kmc_step = 0
        self.procstat[:] = 0
        self._prev_procstat[:] = 0
        self._prev_time = 0.0
        self._rebuild_avail_sites()
        if self._has_lateral and hasattr(self, '_proc_total_rates'):
            self._rebuild_per_site_rates()

    def get_configuration(self):
        """Return a copy of the current lattice configuration."""
        return self.lattice.copy()

    def set_configuration(self, config):
        """Set lattice configuration and rebuild bookkeeping."""
        self.lattice[:] = config
        self._rebuild_avail_sites()
        if self._has_lateral and hasattr(self, '_proc_total_rates'):
            self._rebuild_per_site_rates()

    # ----------------------------------------------------------------
    # Printing
    # ----------------------------------------------------------------

    def print_rates(self):
        """Print all rate constants and available site counts."""
        print("\nRate constants:")
        print(f"  {'Process':<35s} {'k [s^-1]':>12s}  {'N_avail':>8s}")
        print(f"  {'-'*35} {'-'*12}  {'-'*8}")
        for i in range(self.nproc):
            n_avail = len(self._avail_sites[i])
            print(f"  {self.process_names[i]:<35s} "
                  f"{self.rates[i]:>12.4e}  {n_avail:>8d}")
        print()

    def print_coverages(self):
        """Print current surface coverages."""
        cov = self.get_coverage()
        print("Coverages:")
        for name, val in cov.items():
            if val > 1e-6:
                print(f"  {name}: {val:.6f}")

    # ----------------------------------------------------------------
    # Context manager
    # ----------------------------------------------------------------

    def deallocate(self):
        """Clean up resources (for kmos API compatibility)."""
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.deallocate()

    def __repr__(self):
        lat_str = ''
        if self._has_lateral:
            lat_str = f', lateral={len(self._lateral_dict)//2}'
        return (f"KMCEngine(model='{self.project.meta.get('model_name')}', "
                f"size={list(self.lattice_size)}, "
                f"step={self.kmc_step}, time={self.kmc_time:.6e}"
                f"{lat_str})")
