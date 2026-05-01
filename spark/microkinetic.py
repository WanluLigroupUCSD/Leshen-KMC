"""
Mean-field microkinetic modeling using ODE integration.

Solves the coupled ODEs for surface coverages:
    dθ_i/dt = sum_j (ν_ij * r_j)
where ν_ij is the stoichiometric coefficient and r_j is the net rate
of reaction j.

This provides a faster (deterministic) alternative to KMC for
parameter scanning and steady-state analysis.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


class MicroKineticModel:
    """
    Mean-field microkinetic model for surface catalysis.

    Usage:
        mkm = MicroKineticModel()
        mkm.add_species('N2')
        mkm.add_species('NNH')
        mkm.add_species('NH3')

        mkm.add_reaction(
            name='N2_hydrogenation',
            reactants={'N2': 1, 'empty': 0},
            products={'NNH': 1},
            rate_fwd=lambda p: tst_rate(2.82, p['T']),
            rate_rev=lambda p: tst_rate(1.50, p['T']),
        )

        mkm.parameters = {'T': 300, 'U': -0.5}
        ss = mkm.solve_steady_state()
        mkm.print_summary(ss)
    """

    def __init__(self):
        self.species = []
        self.species_id = {}
        self.reactions = []
        self.parameters = {}

    def add_species(self, name):
        """Add a surface species (do NOT add 'empty' - it is implicit)."""
        if name not in ('empty', '*'):
            self.species_id[name] = len(self.species)
            self.species.append(name)

    def add_reaction(self, name, reactants, products,
                     rate_fwd, rate_rev=None, tof_count=None):
        """
        Add an elementary reaction step.

        Parameters
        ----------
        name : str
            Reaction name.
        reactants : dict
            {species_name: stoichiometric_coefficient}.
            Use 'empty' or '*' for vacant sites.
        products : dict
            Same format as reactants.
        rate_fwd : float or callable(params) -> float
            Forward rate constant.
        rate_rev : float, callable, or None
            Reverse rate constant. None = irreversible.
        tof_count : dict, optional
            {observable: coefficient} for TOF tracking.
        """
        self.reactions.append({
            'name': name,
            'reactants': reactants,
            'products': products,
            'rate_fwd': rate_fwd,
            'rate_rev': rate_rev,
            'tof_count': tof_count or {},
        })

    def _get_rate(self, rate_expr):
        """Evaluate a rate expression (constant or callable)."""
        if callable(rate_expr):
            return rate_expr(self.parameters)
        return float(rate_expr)

    def _build_odes(self, theta, t=0):
        """
        Build dθ/dt for all surface species.

        theta[i] = coverage of species i
        theta_empty = 1 - sum(theta)
        """
        n = len(self.species)
        dtheta = np.zeros(n)

        # Clamp to valid range
        theta = np.clip(theta, 0.0, 1.0)
        theta_empty = max(1.0 - np.sum(theta), 0.0)

        cov = {'empty': theta_empty, '*': theta_empty}
        for i, name in enumerate(self.species):
            cov[name] = theta[i]

        for rxn in self.reactions:
            # Forward rate
            kf = self._get_rate(rxn['rate_fwd'])
            rate_f = kf
            for sp, stoich in rxn['reactants'].items():
                rate_f *= cov.get(sp, 0.0) ** stoich

            # Reverse rate
            rate_r = 0.0
            if rxn['rate_rev'] is not None:
                kr = self._get_rate(rxn['rate_rev'])
                rate_r = kr
                for sp, stoich in rxn['products'].items():
                    rate_r *= cov.get(sp, 0.0) ** stoich

            net = rate_f - rate_r

            # Update dθ for each surface species
            for sp, stoich in rxn['reactants'].items():
                if sp in self.species_id:
                    dtheta[self.species_id[sp]] -= stoich * net
            for sp, stoich in rxn['products'].items():
                if sp in self.species_id:
                    dtheta[self.species_id[sp]] += stoich * net

        return dtheta

    def solve_ode(self, t_span, theta0=None, method='BDF',
                  max_step=None, t_eval=None, **kwargs):
        """
        Solve the transient ODE system.

        Parameters
        ----------
        t_span : tuple
            (t_start, t_end) in seconds.
        theta0 : array-like, optional
            Initial coverages (default: all zero = clean surface).
        method : str
            ODE solver method ('BDF', 'Radau', 'LSODA').
        t_eval : array-like, optional
            Times at which to store the solution.

        Returns
        -------
        scipy.integrate.OdeSolution
        """
        n = len(self.species)
        if theta0 is None:
            theta0 = np.zeros(n)

        sol = solve_ivp(
            lambda t, y: self._build_odes(y, t),
            t_span, theta0,
            method=method,
            dense_output=True,
            max_step=max_step or (t_span[1] - t_span[0]) / 1000,
            rtol=1e-8, atol=1e-12,
            t_eval=t_eval,
            **kwargs
        )
        return sol

    def solve_steady_state(self, theta0=None, method='hybr'):
        """
        Find steady-state coverages by solving dθ/dt = 0.

        Parameters
        ----------
        theta0 : array-like, optional
            Initial guess (default: small uniform coverage).
        method : str
            Solver method for fsolve.

        Returns
        -------
        dict
            {species_name: steady-state coverage}
        """
        n = len(self.species)
        if theta0 is None:
            theta0 = np.ones(n) * 0.01

        def equations(theta):
            theta = np.clip(theta, 0.0, 1.0)
            return self._build_odes(theta)

        # Tight tolerance + residual verification. If continuation gives a
        # high-residual solution, fall back to ODE integration from the
        # continuation state to smoothly relax to a proper SS, THEN polish
        # with fsolve.
        #
        # Rationale (2026-04-23 cross-validation on Shaheen kMCOS):
        # - Default fsolve xtol=1.49e-8 causes false convergence in near-onset
        #   regimes where the Jacobian is near-singular (rate spans 10^20).
        # - Continuation (theta0 from previous SS) is the KEY to staying in
        #   the physical (LIVE) basin across a polarization scan. Multi-guess
        #   strategies can flip into DEAD basins (species-saturated, TOF=0).
        # - When continuation+fsolve fails to reach SS, ODE integration from
        #   the current state gives the physically reachable basin reliably.
        def _ode_relax_then_fsolve(start_theta, t_end=1.0, timeout_sec=3):
            """Integrate ODE from start_theta to approach SS, then fsolve.

            2026-04-28 v3 patch: hard signal.alarm timeout 3s per call.
            Even LSODA with t_end=10 deadlocked W Distal near onset (U=+0.030V
            in GS convention, where rates span 10^20 and Jacobian is near
            singular). The fix is a non-negotiable wall-clock cap on the ODE
            attempt; if it fails we fall through to the continuation fsolve
            result rather than perfect SS.
            """
            import signal
            class _OdeTimeout(Exception): pass
            def _on_alarm(signum, frame): raise _OdeTimeout()
            old_handler = signal.signal(signal.SIGALRM, _on_alarm)
            signal.alarm(timeout_sec)
            try:
                ode_sol = solve_ivp(
                    lambda t, y: self._build_odes(y),
                    (0.0, t_end), np.clip(start_theta, 0.0, 1.0),
                    method='LSODA', rtol=1e-5, atol=1e-10,
                    max_step=t_end * 0.1, first_step=1e-9,
                )
                signal.alarm(0)
                theta_dyn = ode_sol.y[:, -1]
                r2, _, _, _ = fsolve(
                    equations, theta_dyn, full_output=True,
                    xtol=1e-12, maxfev=2000)
                return r2, float(np.max(np.abs(equations(r2))))
            except (_OdeTimeout, Exception):
                return None, float('inf')
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)

        def _tof_score(x):
            """Total TOF flux (coverage-weighted) from tagged reactions."""
            if x is None:
                return 0.0
            x_clip = np.clip(x, 0.0, 1.0)
            te = max(1.0 - np.sum(x_clip), 0.0)
            cov = {"empty": te, "*": te}
            for i, name in enumerate(self.species):
                cov[name] = float(x_clip[i])
            tot = 0.0
            for rxn in self.reactions:
                if not rxn['tof_count']:
                    continue
                kf = self._get_rate(rxn['rate_fwd'])
                rf = kf
                for sp, st in rxn['reactants'].items():
                    rf *= cov.get(sp, 0.0) ** st
                rr = 0.0
                if rxn['rate_rev'] is not None:
                    kr = self._get_rate(rxn['rate_rev'])
                    rr = kr
                    for sp, st in rxn['products'].items():
                        rr *= cov.get(sp, 0.0) ** st
                net = rf - rr
                for c in rxn['tof_count'].values():
                    tot += abs(c * net)
            return tot

        # Try two independent approaches and pick the live-basin one.
        # (1) fsolve from continuation theta0
        r_cont, info, ier, msg = fsolve(
            equations, theta0, full_output=True, xtol=1e-12, maxfev=2000)
        resid_cont = float(np.max(np.abs(equations(r_cont))))
        tof_cont = _tof_score(r_cont) if resid_cont < 1e-3 else 0.0

        # (2) ODE from clean surface + fsolve polish (reaches live basin)
        r_ode, resid_ode = _ode_relax_then_fsolve(np.zeros(n))
        tof_ode = _tof_score(r_ode) if (r_ode is not None and resid_ode < 1e-3) else 0.0

        # Prefer whichever has higher TOF (live basin), with physical filter.
        def _is_physical(x):
            if x is None:
                return False
            return np.all(x >= -1e-4) and float(np.sum(np.clip(x, 0, 1))) < 1 + 1e-3

        if _is_physical(r_ode) and tof_ode > tof_cont * 1.1:
            result, resid_max = r_ode, resid_ode
        elif _is_physical(r_cont) and resid_cont < 1e-3:
            result, resid_max = r_cont, resid_cont
        elif _is_physical(r_ode):
            result, resid_max = r_ode, resid_ode
        else:
            result, resid_max = r_cont, resid_cont

        theta_ss = np.clip(result, 0.0, 1.0)
        # Renormalize if total > 1 (defensive — should not trigger at converged SS)
        total = np.sum(theta_ss)
        if total > 1.0:
            theta_ss /= total

        return {name: float(theta_ss[i])
                for i, name in enumerate(self.species)}

    def get_tof(self, theta):
        """
        Calculate turn-over frequencies at given coverages.

        Parameters
        ----------
        theta : dict or array-like
            Surface coverages.

        Returns
        -------
        dict
            {observable: TOF in s^-1}
        """
        if isinstance(theta, dict):
            theta_arr = np.array([theta.get(sp, 0.0)
                                  for sp in self.species])
        else:
            theta_arr = np.array(theta)

        theta_arr = np.clip(theta_arr, 0.0, 1.0)
        theta_empty = max(1.0 - np.sum(theta_arr), 0.0)

        cov = {'empty': theta_empty, '*': theta_empty}
        for i, name in enumerate(self.species):
            cov[name] = theta_arr[i]

        tof = {}
        for rxn in self.reactions:
            kf = self._get_rate(rxn['rate_fwd'])
            rate_f = kf
            for sp, stoich in rxn['reactants'].items():
                rate_f *= cov.get(sp, 0.0) ** stoich

            rate_r = 0.0
            if rxn['rate_rev'] is not None:
                kr = self._get_rate(rxn['rate_rev'])
                rate_r = kr
                for sp, stoich in rxn['products'].items():
                    rate_r *= cov.get(sp, 0.0) ** stoich

            net = rate_f - rate_r

            for obs, coeff in rxn['tof_count'].items():
                tof[obs] = tof.get(obs, 0.0) + coeff * net

        return tof

    def get_reaction_rates(self, theta):
        """
        Get forward and reverse rates for each reaction.

        Returns list of dicts with 'name', 'rate_fwd', 'rate_rev', 'net'.
        """
        if isinstance(theta, dict):
            theta_arr = np.array([theta.get(sp, 0.0)
                                  for sp in self.species])
        else:
            theta_arr = np.array(theta)

        theta_arr = np.clip(theta_arr, 0.0, 1.0)
        theta_empty = max(1.0 - np.sum(theta_arr), 0.0)

        cov = {'empty': theta_empty, '*': theta_empty}
        for i, name in enumerate(self.species):
            cov[name] = theta_arr[i]

        rates = []
        for rxn in self.reactions:
            kf = self._get_rate(rxn['rate_fwd'])
            rf = kf
            for sp, stoich in rxn['reactants'].items():
                rf *= cov.get(sp, 0.0) ** stoich

            rr = 0.0
            if rxn['rate_rev'] is not None:
                kr = self._get_rate(rxn['rate_rev'])
                rr = kr
                for sp, stoich in rxn['products'].items():
                    rr *= cov.get(sp, 0.0) ** stoich

            rates.append({
                'name': rxn['name'],
                'rate_fwd': rf,
                'rate_rev': rr,
                'net': rf - rr,
            })
        return rates

    def degree_of_rate_control(self, theta, observable, dk=0.01):
        """
        Calculate Campbell's degree of rate control (DRC) for each reaction.

        X_RC,i = (k_i / r) * (dr / dk_i)

        Parameters
        ----------
        theta : dict
            Steady-state coverages.
        observable : str
            TOF observable name.
        dk : float
            Fractional perturbation of rate constant.

        Returns
        -------
        dict
            {reaction_name: DRC value}
        """
        tof_base = self.get_tof(theta)
        r_base = tof_base.get(observable, 0.0)
        if abs(r_base) < 1e-30:
            return {rxn['name']: 0.0 for rxn in self.reactions}

        drc = {}
        for i, rxn in enumerate(self.reactions):
            # Perturb forward rate constant
            original_fwd = rxn['rate_fwd']

            if callable(original_fwd):
                k_orig = original_fwd(self.parameters)
                rxn['rate_fwd'] = lambda p, k=k_orig, d=dk: k * (1 + d)
            else:
                k_orig = float(original_fwd)
                rxn['rate_fwd'] = k_orig * (1 + dk)

            # Also perturb reverse if it exists
            original_rev = rxn['rate_rev']
            if original_rev is not None:
                if callable(original_rev):
                    kr_orig = original_rev(self.parameters)
                    rxn['rate_rev'] = lambda p, k=kr_orig, d=dk: k * (1 + d)
                else:
                    kr_orig = float(original_rev)
                    rxn['rate_rev'] = kr_orig * (1 + dk)

            # Recalculate steady state and TOF
            ss_pert = self.solve_steady_state(
                theta0=[theta.get(sp, 0.0) for sp in self.species])
            tof_pert = self.get_tof(ss_pert)
            r_pert = tof_pert.get(observable, 0.0)

            drc[rxn['name']] = (r_pert - r_base) / (dk * r_base)

            # Restore
            rxn['rate_fwd'] = original_fwd
            rxn['rate_rev'] = original_rev

        return drc

    def print_summary(self, theta=None):
        """Print comprehensive model summary."""
        if theta is None:
            theta = self.solve_steady_state()

        print("=" * 65)
        print("  Microkinetic Model Summary")
        print("=" * 65)

        print(f"\n  Parameters:")
        for name, val in self.parameters.items():
            print(f"    {name} = {val}")

        theta_empty = 1.0 - sum(theta.values())
        print(f"\n  Steady-state coverages:")
        print(f"    {'empty':<20s} {theta_empty:.6e}")
        for name, val in sorted(theta.items(), key=lambda x: -x[1]):
            if val > 1e-12:
                print(f"    {name:<20s} {val:.6e}")

        tof = self.get_tof(theta)
        if tof:
            print(f"\n  Turn-over frequencies [s^-1]:")
            for obs, val in tof.items():
                print(f"    {obs:<20s} {val:.6e}")

        rates = self.get_reaction_rates(theta)
        print(f"\n  Reaction rates [s^-1]:")
        print(f"    {'Reaction':<30s} {'Forward':>12s} {'Reverse':>12s} "
              f"{'Net':>12s}")
        print(f"    {'-'*30} {'-'*12} {'-'*12} {'-'*12}")
        for r in rates:
            print(f"    {r['name']:<30s} {r['rate_fwd']:>12.4e} "
                  f"{r['rate_rev']:>12.4e} {r['net']:>12.4e}")

        print("=" * 65)

    def scan_parameter(self, param_name, values, observable=None):
        """
        Scan a parameter and compute steady-state coverages and TOFs.

        Parameters
        ----------
        param_name : str
            Parameter to scan.
        values : array-like
            Values to scan over.
        observable : str, optional
            TOF observable to track.

        Returns
        -------
        dict
            Results with 'values', 'coverages', 'tofs'.
        """
        results = {
            'values': [],
            'coverages': [],
            'tofs': [],
        }

        for val in values:
            self.parameters[param_name] = val
            ss = self.solve_steady_state()
            tof = self.get_tof(ss)

            results['values'].append(val)
            results['coverages'].append(ss)
            if observable:
                results['tofs'].append(tof.get(observable, 0.0))
            else:
                results['tofs'].append(tof)

        return results
