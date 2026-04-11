"""
Dynamic KMC Engine — BKL rejection-free algorithm with local-only updates.

Core algorithm per step (design doc Section 13):
  1. Select event μ with probability r_μ / R_tot
  2. Execute event μ (update adsorbate/site state)
  3. Identify affected local region N(μ)
  4. Recompute descriptors only in N(μ)
  5. For each changed environment: cache hit → reuse; miss → estimate
  6. If uncertainty high → queue refinement
  7. Update local rates, R_tot, advance t += -ln(u)/R_tot
"""

import numpy as np
import time as _time
from collections import defaultdict

from .surface import DynamicSurface
from .descriptor import EnvHash
from .events import EventGenerator, EventType
from .cache import EventCache


class DynamicKMCEngine:
    """
    Environment-dependent dynamic catalytic KMC engine.

    Supports both catalytic events (adsorption, desorption, reaction,
    diffusion) and structural events (site conversion, segregation)
    in a unified BKL rejection-free algorithm with local-only updates
    and environment-aware event caching.

    Parameters
    ----------
    surface : DynamicSurface
        The dynamic surface to simulate on.
    rate_estimator : RateEstimator
        Provides rates for events based on local environment.
    event_generator : EventGenerator
        Generates possible events for each site.
    temperature : float
        Simulation temperature in Kelvin.
    cache_size : int
        Maximum event cache entries (0=unlimited).
    refinement_threshold : float
        Uncertainty threshold for queuing expensive refinement.
    debug : bool
    """

    def __init__(self, surface, rate_estimator, event_generator=None,
                 temperature=300.0, cache_size=0,
                 refinement_threshold=0.5, debug=False):
        self.surface = surface
        self.estimator = rate_estimator
        self.generator = event_generator or EventGenerator()
        self.temperature = temperature
        self.debug = debug
        self.refinement_threshold = refinement_threshold

        # Event cache
        self.cache = EventCache(max_size=cache_size, debug=debug)

        # Simulation state
        self.kmc_time = 0.0
        self.kmc_step = 0

        # Per-site event lists and rates
        self._site_events = [[] for _ in range(surface.n_sites)]
        self._site_rate_sum = np.zeros(surface.n_sites)
        self._total_rate = 0.0

        # TOF tracking
        self._tof_counts = defaultdict(int)
        self._prev_tof_counts = defaultdict(int)
        self._prev_time = 0.0

        # Process statistics
        self.event_stats = defaultdict(int)

        # Refinement queue
        self._refinement_queue = []

        # Trajectory
        self._trajectory = []

        # Build initial event lists for all sites
        self._rebuild_all()

    # ── Initialization ────────────────────────────────────────────────

    def _rebuild_all(self):
        """Build event lists and rates for all sites from scratch."""
        self._total_rate = 0.0
        for i in range(self.surface.n_sites):
            self._rebuild_site(i)

        if self.debug:
            n_events = sum(len(el) for el in self._site_events)
            print(f"ENGINE: Built {n_events} events across "
                  f"{self.surface.n_sites} sites, "
                  f"R_tot={self._total_rate:.4e} s⁻¹")
            self.cache.summary()

    def _rebuild_site(self, site_idx):
        """Rebuild event list and rates for a single site."""
        # Remove old contribution
        self._total_rate -= self._site_rate_sum[site_idx]

        # Get local environment
        env = EnvHash.from_surface(self.surface, site_idx)

        # Try cache first
        cached = self.cache.lookup(env)
        if cached is not None:
            # Instantiate cached events with actual site index
            events = self._instantiate_cached(cached, site_idx)
        else:
            # Generate fresh events
            events = self.generator.generate(self.surface, site_idx)

            # Assign rates from estimator
            for event in events:
                rate, uncertainty = self.estimator.estimate(event, env)
                event.rate = rate
                if uncertainty > self.refinement_threshold:
                    self._refinement_queue.append((site_idx, event, env))

            # Cache the event set (as templates)
            templates = self._make_templates(events, site_idx)
            self.cache.store(env, templates)

        # Filter out zero-rate and impossible events
        active_events = [e for e in events
                         if e.rate > 0 and e.is_possible(self.surface)]

        self._site_events[site_idx] = active_events
        site_rate = sum(e.rate for e in active_events)
        self._site_rate_sum[site_idx] = site_rate
        self._total_rate += site_rate

    def _make_templates(self, events, site_idx):
        """Convert concrete events into cacheable templates."""
        templates = []
        for e in events:
            templates.append({
                'class': type(e),
                'event_type': e.event_type,
                'rate': e.rate,
                'name': e.name,
                'attrs': self._extract_event_attrs(e, site_idx),
            })
        return templates

    def _extract_event_attrs(self, event, site_idx):
        """Extract event-specific attributes relative to the site."""
        attrs = {}
        if hasattr(event, 'species_id'):
            attrs['species_id'] = event.species_id
        if hasattr(event, 'target_site'):
            # Store as NN offset index, not absolute site index
            nn = self.surface.neighbors[site_idx]
            if event.target_site in nn:
                attrs['nn_offset'] = nn.index(event.target_site)
        if hasattr(event, 'partner_site'):
            nn = self.surface.neighbors[site_idx]
            if event.partner_site is not None and event.partner_site in nn:
                attrs['partner_nn_offset'] = nn.index(event.partner_site)
        if hasattr(event, 'reactant_species'):
            attrs['reactant_species'] = event.reactant_species
        if hasattr(event, 'product_species'):
            attrs['product_species'] = event.product_species
        if hasattr(event, 'from_type'):
            attrs['from_type'] = event.from_type
        if hasattr(event, 'to_type'):
            attrs['to_type'] = event.to_type
        if hasattr(event, 'tof_count'):
            attrs['tof_count'] = event.tof_count
        return attrs

    def _instantiate_cached(self, cached, site_idx):
        """Instantiate cached event templates for a specific site."""
        from .events import (Adsorption, Desorption, SurfaceReaction,
                             Diffusion, SiteConversion, Segregation)

        events = []
        nn = self.surface.neighbors[site_idx]

        for tmpl in cached.event_templates:
            attrs = tmpl['attrs']
            rate = tmpl['rate']
            name = tmpl['name']
            et = tmpl['event_type']

            if et == EventType.ADSORPTION:
                events.append(Adsorption(site_idx, attrs['species_id'],
                                         rate, name))
            elif et == EventType.DESORPTION:
                events.append(Desorption(site_idx, attrs['species_id'],
                                         rate, name))
            elif et == EventType.DIFFUSION:
                offset = attrs.get('nn_offset', 0)
                if offset < len(nn):
                    events.append(Diffusion(site_idx, nn[offset], rate, name))
            elif et == EventType.SURFACE_REACTION:
                partner_offset = attrs.get('partner_nn_offset')
                partner = nn[partner_offset] if (partner_offset is not None
                                                  and partner_offset < len(nn)) else None
                events.append(SurfaceReaction(
                    site_idx, partner,
                    attrs.get('reactant_species', (0,)),
                    attrs.get('product_species', (0,)),
                    rate, name,
                    attrs.get('tof_count', {}),
                ))
            elif et == EventType.SITE_CONVERSION:
                events.append(SiteConversion(
                    site_idx, attrs['from_type'], attrs['to_type'],
                    rate, name))
            elif et == EventType.SEGREGATION:
                offset = attrs.get('nn_offset', attrs.get('partner_nn_offset', 0))
                if offset < len(nn):
                    events.append(Segregation(site_idx, nn[offset], rate, name))

        return events

    # ── KMC step ──────────────────────────────────────────────────────

    def do_step(self):
        """
        Execute one KMC step (BKL rejection-free algorithm).

        Returns
        -------
        bool
            True if a step was executed, False if system is frozen.
        """
        if self._total_rate <= 0.0:
            return False

        # 1. Draw random numbers
        r_time = np.random.random()
        r_event = np.random.random()

        # 2. Time advancement
        dt = -np.log(r_time) / self._total_rate
        self.kmc_time += dt

        # 3. Event selection: linear scan on per-site cumulative rates
        target = r_event * self._total_rate
        cumul = 0.0
        selected_event = None

        for i in range(self.surface.n_sites):
            cumul += self._site_rate_sum[i]
            if cumul >= target:
                # Select event within this site
                site_events = self._site_events[i]
                target_in_site = target - (cumul - self._site_rate_sum[i])
                sub_cumul = 0.0
                for event in site_events:
                    sub_cumul += event.rate
                    if sub_cumul >= target_in_site:
                        selected_event = event
                        break
                if selected_event is None and site_events:
                    selected_event = site_events[-1]
                break

        if selected_event is None:
            return False

        # 4. Execute event
        affected_sites = selected_event.execute(self.surface)

        # 5. Track statistics
        self.kmc_step += 1
        self.event_stats[selected_event.event_type] += 1

        # TOF tracking
        if hasattr(selected_event, 'tof_count'):
            for obs, coeff in selected_event.tof_count.items():
                self._tof_counts[obs] += coeff

        # Record trajectory
        self._trajectory.append({
            'step': self.kmc_step,
            'time': self.kmc_time,
            'dt': dt,
            'event_type': selected_event.event_type.name,
            'site': selected_event.site,
            'name': selected_event.name,
        })

        # 6. Local update: rebuild only affected sites
        for site_idx in affected_sites:
            if site_idx < self.surface.n_sites:
                self._rebuild_site(site_idx)

        return True

    def run(self, max_steps=10000, max_time=None, report_interval=None,
            callback=None):
        """
        Run the dynamic KMC simulation.

        Parameters
        ----------
        max_steps : int
        max_time : float, optional
            Maximum simulated KMC time in seconds.
        report_interval : int, optional
            Print progress every N steps.
        callback : callable, optional
            Called after each step: callback(engine) → bool (True to stop).

        Returns
        -------
        trajectory : list of dict
        """
        t_wall_start = _time.time()

        if report_interval is None:
            report_interval = max(1, max_steps // 10)

        for step in range(max_steps):
            if not self.do_step():
                print(f"System frozen at step {self.kmc_step}, "
                      f"t={self.kmc_time:.6e} s")
                break

            if max_time is not None and self.kmc_time >= max_time:
                if self.debug:
                    print(f"Reached max KMC time: {self.kmc_time:.6e} s")
                break

            if self.debug and (step + 1) % report_interval == 0:
                pct = 100 * (step + 1) / max_steps
                print(f"  [{pct:5.1f}%] step={self.kmc_step}, "
                      f"t={self.kmc_time:.6e} s, "
                      f"R_tot={self._total_rate:.4e}")

            if callback is not None:
                if callback(self):
                    break

        t_wall = _time.time() - t_wall_start
        if self.debug:
            print(f"\nSimulation done: {self.kmc_step} steps in "
                  f"{t_wall:.1f}s wall time")

        return self._trajectory

    # ── Observables ───────────────────────────────────────────────────

    def get_coverage(self, species_name=None):
        """Get fractional coverage."""
        if species_name is not None:
            sp_id = self.surface.species_id(species_name)
            return self.surface.get_coverage(sp_id)
        return self.surface.get_coverage()

    def get_composition(self, layer=None):
        """Get surface composition."""
        return self.surface.get_composition(layer)

    def get_tof(self):
        """
        Get turn-over frequencies since last call.

        Returns dict of {observable: TOF in s^-1 per site}.
        """
        dt = self.kmc_time - self._prev_time
        if dt <= 0:
            return {}
        tof = {}
        for obs in set(list(self._tof_counts.keys()) +
                       list(self._prev_tof_counts.keys())):
            delta = self._tof_counts[obs] - self._prev_tof_counts.get(obs, 0)
            tof[obs] = delta / (dt * self.surface.n_sites)

        self._prev_tof_counts = dict(self._tof_counts)
        self._prev_time = self.kmc_time
        return tof

    def get_event_stats(self):
        """Get execution count per event type."""
        return {k.name: v for k, v in self.event_stats.items()}

    def get_site_type_distribution(self):
        """Get current distribution of (atom_type, adsorbate) combos."""
        return self.surface.get_site_type_distribution()

    def get_unique_environments(self):
        """Count unique local environments on the surface."""
        return EnvHash.count_unique(self.surface)

    # ── Printing ──────────────────────────────────────────────────────

    def summary(self):
        """Print simulation summary."""
        print(f"\n{'='*50}")
        print(f"Dynamic KMC Simulation Summary")
        print(f"{'='*50}")
        print(f"  Steps:       {self.kmc_step}")
        print(f"  Time:        {self.kmc_time:.6e} s")
        print(f"  Temperature: {self.temperature} K")
        print(f"  Sites:       {self.surface.n_sites}")
        print(f"  Composition: {self.get_composition()}")
        print(f"  Coverage:    {self.get_coverage():.4f}")
        print(f"\n  Event statistics:")
        for name, count in sorted(self.get_event_stats().items()):
            pct = 100 * count / max(self.kmc_step, 1)
            print(f"    {name:<25s} {count:>8d} ({pct:.1f}%)")
        n_envs = len(self.get_unique_environments())
        print(f"\n  Unique environments: {n_envs}")
        self.cache.summary()
        if self._refinement_queue:
            print(f"  Refinement queue: {len(self._refinement_queue)} pending")

    def __repr__(self):
        return (f"DynamicKMCEngine(sites={self.surface.n_sites}, "
                f"step={self.kmc_step}, t={self.kmc_time:.6e}s)")
