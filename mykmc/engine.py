"""
Core Kinetic Monte Carlo engine implementing the BKL rejection-free algorithm.

The Variable Step Size Method (VSSM), also known as the BKL algorithm or
n-fold way, ensures every step advances the simulation clock and changes
the system state (no rejected moves).

Algorithm per step:
  1. Compute cumulative rates: R_i = sum_{j<=i} k_j * n_j
  2. Draw random numbers r1, r2, r3 ~ U(0,1)
  3. Advance time: dt = -ln(r1) / R_total
  4. Select process by binary search on cumulative rates using r2
  5. Select site uniformly from available sites for chosen process using r3
  6. Execute process (update lattice + bookkeeping)
"""

import numpy as np
import time as _time
from collections import defaultdict
from .rates import evaluate_rate_expression


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

        # Parameters (with auto-update callback)
        self.parameters = ParameterProxy(project.parameter_list, engine=self)

        # Initialize lattice
        self._init_lattice()

        # Initialize processes
        self._init_processes()

        # Compute max offset range for neighbor updates
        self._compute_max_offset()

        # Rate constants
        self.rates = np.zeros(self.nproc)
        self._update_rate_constants()

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

        if print_rates:
            self.print_rates()

    # ----------------------------------------------------------------
    # Initialization
    # ----------------------------------------------------------------

    def _init_lattice(self):
        """Initialize the lattice array with default species."""
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

    def _init_processes(self):
        """Convert process definitions to numeric arrays for fast lookup."""
        self.nproc = len(self.project.process_list)
        self.process_names = [p.name for p in self.project.process_list]

        self._proc_conditions = []
        self._proc_actions = []
        self._proc_rate_exprs = []
        self._proc_tof_count = []

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

        # Available sites bookkeeping: O(1) add/remove via swap-with-last
        self._avail_sites = [[] for _ in range(self.nproc)]
        self._site_in_avail = [dict() for _ in range(self.nproc)]

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
        """Check if process can occur at site."""
        coord = self._site_to_coord(site)
        for offset, sp_id in self._proc_conditions[proc_id]:
            neighbor = tuple(c + o for c, o in zip(coord, offset))
            if self.lattice[self._coord_to_site(neighbor)] != sp_id:
                return False
        return True

    def _add_to_avail(self, proc_id, site):
        if site not in self._site_in_avail[proc_id]:
            self._site_in_avail[proc_id][site] = len(
                self._avail_sites[proc_id])
            self._avail_sites[proc_id].append(site)

    def _remove_from_avail(self, proc_id, site):
        idx_map = self._site_in_avail[proc_id]
        if site in idx_map:
            idx = idx_map[site]
            avail = self._avail_sites[proc_id]
            last_site = avail[-1]
            avail[idx] = last_site
            idx_map[last_site] = idx
            avail.pop()
            del idx_map[site]

    def _rebuild_avail_sites(self):
        """Full rebuild of available sites (used at initialization)."""
        for p in range(self.nproc):
            self._avail_sites[p] = []
            self._site_in_avail[p] = {}
        for s in range(self.nsites):
            for p in range(self.nproc):
                if self._check_process_at_site(p, s):
                    self._add_to_avail(p, s)

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
        """Update available sites for all affected sites."""
        for site in affected_sites:
            for p in range(self.nproc):
                if self._check_process_at_site(p, site):
                    self._add_to_avail(p, site)
                else:
                    self._remove_from_avail(p, site)

    # ----------------------------------------------------------------
    # Rate constants
    # ----------------------------------------------------------------

    def _update_rate_constants(self):
        """Evaluate all rate expressions with current parameters."""
        params = {p.name: p.value for p in self.project.parameter_list}
        for i, expr in enumerate(self._proc_rate_exprs):
            self.rates[i] = evaluate_rate_expression(expr, params)

    def _update_accum_rates(self):
        """Build cumulative rate array for process selection."""
        total = 0.0
        for p in range(self.nproc):
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

        # Site selection (uniform among available)
        n_avail = len(self._avail_sites[proc_id])
        if n_avail == 0:
            return False
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

    def get_configuration(self):
        """Return a copy of the current lattice configuration."""
        return self.lattice.copy()

    def set_configuration(self, config):
        """Set lattice configuration and rebuild bookkeeping."""
        self.lattice[:] = config
        self._rebuild_avail_sites()

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
        return (f"KMCEngine(model='{self.project.meta.get('model_name')}', "
                f"size={list(self.lattice_size)}, "
                f"step={self.kmc_step}, time={self.kmc_time:.6e})")
