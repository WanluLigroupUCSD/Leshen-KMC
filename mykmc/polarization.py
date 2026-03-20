"""
Polarization curve calculation from constant-potential DFT data.

Workflow:
  1. Provide DFT energies of adsorbates and transition states at multiple potentials
  2. Build microkinetic model with interpolated potential-dependent barriers
  3. Sweep potential -> solve steady state -> compute TOF -> convert to current density

Two input modes:
  (a) Direct Ea(U) tables:  InterpolatedBarrier
  (b) Energy surfaces E_IS(U), E_TS(U), E_FS(U):  EnergyLandscape

Example usage (direct barriers):
    from mykmc.polarization import InterpolatedBarrier, PolarizationCurve

    # DFT data: Ea at several potentials
    barrier = InterpolatedBarrier(
        U_data=[-2.0, -1.5, -1.0, -0.5, 0.0],
        Ea_data=[0.50, 0.82, 1.14, 1.58, 2.10],
    )

    # Use as rate function in MicroKineticModel
    mkm.add_reaction(
        name='N2_to_NNH',
        reactants={'N2': 1}, products={'NNH': 1},
        rate_fwd=barrier.make_rate_func(),
        tof_count={'N2_consumption': 1},
    )

    # Compute polarization curve
    pc = PolarizationCurve(mkm, n_electrons={'NH3_production': 6},
                           A_site=(3.2e-10)**2)
    results = pc.compute(U_range=np.linspace(-2, 0, 50), T=300)

Example usage (energy landscape from DFT):
    from mykmc.polarization import EnergyLandscape, PolarizationCurve

    potentials = [-2.0, -1.5, -1.0, -0.5, 0.0]
    state_energies = {
        'clean':  [0.0,  0.0,  0.0,  0.0,  0.0],
        'N2*':    [-0.80, -0.72, -0.60, -0.50, -0.38],
        'NNH*':   [1.50,  1.65,  1.82,  2.00,  2.20],
    }
    ts_energies = {
        ('N2*', 'NNH*'): [1.80, 2.00, 2.30, 2.60, 3.00],
    }

    landscape = EnergyLandscape(potentials, state_energies, ts_energies)
    # Ea_fwd = E_TS - E_IS, Ea_rev = E_TS - E_FS  are auto-computed
"""

import numpy as np
from scipy.interpolate import interp1d
from .rates import tst_rate, hertz_knudsen
from .units import kB, h, eV, e_charge
from .microkinetic import MicroKineticModel


class InterpolatedBarrier:
    """
    Interpolated activation barrier Ea(U) from tabulated DFT data.

    Parameters
    ----------
    U_data : array-like
        Potentials [V vs RHE] at which barriers were computed.
    Ea_data : array-like
        Activation energies [eV] at each potential.
    kind : str
        Interpolation method ('linear', 'cubic', 'quadratic').
    fill_value : str or float
        Extrapolation handling. 'extrapolate' or a constant value.
    """

    def __init__(self, U_data, Ea_data, kind='cubic', fill_value='extrapolate'):
        self.U_data = np.asarray(U_data, dtype=float)
        self.Ea_data = np.asarray(Ea_data, dtype=float)

        if len(self.U_data) < 2:
            raise ValueError("Need at least 2 data points for interpolation")
        if len(self.U_data) != len(self.Ea_data):
            raise ValueError("U_data and Ea_data must have the same length")

        # Fall back to linear if too few points for cubic
        if len(self.U_data) <= 3 and kind == 'cubic':
            kind = 'quadratic' if len(self.U_data) == 3 else 'linear'

        self._interp = interp1d(
            self.U_data, self.Ea_data,
            kind=kind, fill_value=fill_value, bounds_error=False,
        )
        self.kind = kind

    def __call__(self, U):
        """Return Ea [eV] at potential U [V]. Clamps to >= 0."""
        Ea = float(self._interp(U))
        return max(Ea, 0.0)

    def make_rate_func(self):
        """
        Create a TST rate function for use in MicroKineticModel.

        Returns callable(params) -> float  where params has 'T' and 'U'.
        k(U, T) = (kB*T/h) * exp(-Ea(U) / (kB*T))
        """
        barrier = self

        def rate_func(params):
            Ea = barrier(params['U'])
            return tst_rate(Ea, params['T'])

        return rate_func


class EnergyLandscape:
    """
    Compute potential-dependent barriers from a DFT energy landscape.

    Given E_IS(U), E_TS(U), E_FS(U) at several potentials, computes:
        Ea_fwd(U) = max(E_TS(U) - E_IS(U), 0)
        Ea_rev(U) = max(E_TS(U) - E_FS(U), 0)

    Parameters
    ----------
    potentials : array-like
        Potentials [V vs RHE] from constant-potential DFT.
    state_energies : dict
        {state_name: [E(U1), E(U2), ...]}  energies [eV] per state.
    ts_energies : dict
        {(IS_name, FS_name): [E_TS(U1), ...]}  TS energies [eV].
    kind : str
        Interpolation method.
    """

    def __init__(self, potentials, state_energies, ts_energies, kind='cubic'):
        self.potentials = np.asarray(potentials, dtype=float)
        self.state_energies = {
            k: np.asarray(v, dtype=float) for k, v in state_energies.items()
        }
        self.ts_energies = {
            k: np.asarray(v, dtype=float) for k, v in ts_energies.items()
        }

        # Validate lengths
        n = len(self.potentials)
        for name, arr in self.state_energies.items():
            if len(arr) != n:
                raise ValueError(
                    f"state_energies['{name}'] has {len(arr)} points, "
                    f"expected {n}")
        for key, arr in self.ts_energies.items():
            if len(arr) != n:
                raise ValueError(
                    f"ts_energies[{key}] has {len(arr)} points, expected {n}")

        # Pre-compute interpolated barriers
        self.barriers_fwd = {}
        self.barriers_rev = {}

        for (is_name, fs_name), E_TS in self.ts_energies.items():
            if is_name not in self.state_energies:
                raise KeyError(
                    f"Initial state '{is_name}' not in state_energies")
            if fs_name not in self.state_energies:
                raise KeyError(
                    f"Final state '{fs_name}' not in state_energies")

            E_IS = self.state_energies[is_name]
            E_FS = self.state_energies[fs_name]

            Ea_fwd = np.maximum(E_TS - E_IS, 0.0)
            Ea_rev = np.maximum(E_TS - E_FS, 0.0)

            key = (is_name, fs_name)
            self.barriers_fwd[key] = InterpolatedBarrier(
                self.potentials, Ea_fwd, kind=kind)
            self.barriers_rev[key] = InterpolatedBarrier(
                self.potentials, Ea_rev, kind=kind)

    def get_barrier_fwd(self, is_name, fs_name):
        """Get forward barrier interpolator for IS -> FS."""
        return self.barriers_fwd[(is_name, fs_name)]

    def get_barrier_rev(self, is_name, fs_name):
        """Get reverse barrier interpolator for FS -> IS."""
        return self.barriers_rev[(is_name, fs_name)]

    def get_rate_funcs(self, is_name, fs_name):
        """
        Get (rate_fwd, rate_rev) callables for a reaction step.

        Both are compatible with MicroKineticModel.add_reaction().
        """
        key = (is_name, fs_name)
        return (self.barriers_fwd[key].make_rate_func(),
                self.barriers_rev[key].make_rate_func())

    def summary(self, U=None):
        """Print barrier summary at given potential(s)."""
        potentials = [U] if U is not None else self.potentials
        print("Energy Landscape Summary")
        print("=" * 70)
        for (is_name, fs_name) in self.ts_energies:
            print(f"\n  {is_name} -> {fs_name}:")
            print(f"    {'U [V]':>8s}  {'Ea_fwd [eV]':>12s}  "
                  f"{'Ea_rev [eV]':>12s}")
            print(f"    {'-' * 8}  {'-' * 12}  {'-' * 12}")
            key = (is_name, fs_name)
            for u in potentials:
                ea_f = self.barriers_fwd[key](u)
                ea_r = self.barriers_rev[key](u)
                print(f"    {u:>8.3f}  {ea_f:>12.4f}  {ea_r:>12.4f}")


def tof_to_current_density(tof, n_electrons, A_site):
    """
    Convert turnover frequency to current density.

    j [mA/cm^2] = n_e * e * TOF / A_site * 1e-4 * 1e3
                = n_e * e * TOF / A_site * 0.1

    Parameters
    ----------
    tof : float
        Turnover frequency [s^-1 per site].
    n_electrons : int or float
        Number of electrons transferred per turnover event.
    A_site : float
        Area per active site [m^2].

    Returns
    -------
    float
        Current density [mA/cm^2].
    """
    # j [A/m^2] = n_e * e * TOF / A_site
    j_SI = n_electrons * e_charge * tof / A_site
    # A/m^2 -> mA/cm^2:  1 A/m^2 = 0.1 mA/cm^2
    return j_SI * 0.1


class PolarizationCurve:
    """
    Compute electrochemical polarization curve (j-U) from a microkinetic
    model with potential-dependent rate constants.

    Parameters
    ----------
    mkm : MicroKineticModel
        Model whose rate functions depend on params['U'].
    n_electrons : dict
        {tof_observable: n_electrons_per_event}.
        For NRR producing NH3:  {'NH3_production': 6}  (6 e- per NH3).
    A_site : float
        Active site area [m^2].
    """

    def __init__(self, mkm, n_electrons, A_site):
        self.mkm = mkm
        self.n_electrons = n_electrons
        self.A_site = A_site
        self.results = None

    def compute(self, U_range, T=300.0, verbose=True):
        """
        Sweep potential and compute TOF + current at each point.

        Uses continuation (previous steady-state as next initial guess)
        for robust convergence across the potential range.

        Parameters
        ----------
        U_range : array-like
            Potentials [V vs RHE].
        T : float
            Temperature [K].
        verbose : bool
            Print progress.

        Returns
        -------
        dict with keys:
            'U'             : potential array [V]
            'tof'           : {observable: TOF array [s^-1]}
            'current'       : {observable: j array [mA/cm^2]}
            'total_current' : total j array [mA/cm^2]
            'coverages'     : list of {species: theta} dicts
            'log10_j'       : log10(|j_total|) for Tafel plot
        """
        U_range = np.asarray(U_range, dtype=float)
        self.mkm.parameters['T'] = T

        tof_arrays = {obs: [] for obs in self.n_electrons}
        j_arrays = {obs: [] for obs in self.n_electrons}
        j_total = []
        coverages = []

        if verbose:
            print(f"Computing polarization curve: {len(U_range)} points")
            print(f"  T = {T} K, U = [{U_range[0]:.2f}, "
                  f"{U_range[-1]:.2f}] V vs RHE")
            print(f"  Site area = {self.A_site:.4e} m^2")
            print()

        theta_guess = None
        report_interval = max(1, len(U_range) // 10)

        for i, U in enumerate(U_range):
            self.mkm.parameters['U'] = U

            if theta_guess is not None:
                ss = self.mkm.solve_steady_state(theta0=theta_guess)
            else:
                ss = self.mkm.solve_steady_state()

            theta_guess = [ss.get(sp, 0.0) for sp in self.mkm.species]

            tof = self.mkm.get_tof(ss)

            j_sum = 0.0
            for obs, n_e in self.n_electrons.items():
                tof_val = tof.get(obs, 0.0)
                tof_arrays[obs].append(tof_val)
                j_val = tof_to_current_density(tof_val, n_e, self.A_site)
                j_arrays[obs].append(j_val)
                j_sum += j_val

            j_total.append(j_sum)
            coverages.append(ss)

            if verbose and (i + 1) % report_interval == 0:
                tof_str = ', '.join(
                    f"{obs}={tof.get(obs, 0):.3e}"
                    for obs in self.n_electrons)
                print(f"  [{100*(i+1)/len(U_range):5.1f}%] "
                      f"U = {U:+.3f} V, j = {j_sum:.4e} mA/cm^2, "
                      f"TOF: {tof_str}")

        for obs in tof_arrays:
            tof_arrays[obs] = np.array(tof_arrays[obs])
            j_arrays[obs] = np.array(j_arrays[obs])
        j_total = np.array(j_total)

        with np.errstate(divide='ignore'):
            log10_j = np.where(
                np.abs(j_total) > 0,
                np.log10(np.abs(j_total)),
                -30.0)

        self.results = {
            'U': U_range,
            'tof': tof_arrays,
            'current': j_arrays,
            'total_current': j_total,
            'coverages': coverages,
            'log10_j': log10_j,
            'T': T,
            'A_site': self.A_site,
        }

        if verbose:
            j_max = np.max(np.abs(j_total))
            print(f"\nDone. Max |j| = {j_max:.4e} mA/cm^2")

        return self.results

    def print_results(self):
        """Print formatted results table."""
        if self.results is None:
            print("No results yet. Call compute() first.")
            return

        r = self.results
        obs_names = list(self.n_electrons.keys())

        # Build header
        header = f"{'U [V]':>8s}"
        for obs in obs_names:
            header += f"  {'TOF_' + obs:>14s}"
        header += f"  {'j [mA/cm2]':>14s}  {'log10|j|':>10s}"
        for sp in self.mkm.species:
            header += f"  {'theta_' + sp:>12s}"

        print("\nPolarization Curve Results")
        print("=" * len(header))
        print(header)
        print("-" * len(header))

        for i, U in enumerate(r['U']):
            line = f"{U:>8.3f}"
            for obs in obs_names:
                line += f"  {r['tof'][obs][i]:>14.6e}"
            line += f"  {r['total_current'][i]:>14.6e}"
            line += f"  {r['log10_j'][i]:>10.3f}"
            cov = r['coverages'][i]
            for sp in self.mkm.species:
                line += f"  {cov.get(sp, 0.0):>12.6e}"
            print(line)

    def save(self, filename):
        """Save results to tab-separated file."""
        if self.results is None:
            print("No results yet. Call compute() first.")
            return

        r = self.results
        obs_names = list(self.n_electrons.keys())

        with open(filename, 'w') as f:
            cols = ['U_V']
            for obs in obs_names:
                cols.extend([f'TOF_{obs}', f'j_{obs}_mA_cm2'])
            cols.extend(['j_total_mA_cm2', 'log10_abs_j'])
            for sp in self.mkm.species:
                cols.append(f'theta_{sp}')
            cols.append('theta_empty')

            f.write('# ' + '\t'.join(cols) + '\n')
            f.write(f'# T = {r["T"]} K, '
                    f'A_site = {r["A_site"]:.4e} m^2\n')

            for i, U in enumerate(r['U']):
                vals = [f"{U:.4f}"]
                for obs in obs_names:
                    vals.append(f"{r['tof'][obs][i]:.8e}")
                    vals.append(f"{r['current'][obs][i]:.8e}")
                vals.append(f"{r['total_current'][i]:.8e}")
                vals.append(f"{r['log10_j'][i]:.4f}")

                cov = r['coverages'][i]
                theta_sum = 0.0
                for sp in self.mkm.species:
                    v = cov.get(sp, 0.0)
                    vals.append(f"{v:.8e}")
                    theta_sum += v
                vals.append(f"{max(1.0 - theta_sum, 0.0):.8e}")

                f.write('\t'.join(vals) + '\n')

        print(f"Results saved to {filename}")

    @staticmethod
    def from_energy_landscape(landscape, species_list, reactions,
                              adsorption_reactions=None,
                              n_electrons=None, A_site=None,
                              extra_params=None):
        """
        Build a PolarizationCurve directly from an EnergyLandscape.

        Parameters
        ----------
        landscape : EnergyLandscape
            Energy landscape with interpolated barriers.
        species_list : list of str
            Surface species names (excluding 'empty').
        reactions : list of dict
            Each dict has keys:
              'name'       : str, reaction name
              'is'         : str, initial-state name in landscape
              'fs'         : str, final-state name in landscape
              'reactants'  : dict {species: stoich}
              'products'   : dict {species: stoich}
              'tof_count'  : dict, optional
              'reversible' : bool (default True)
        adsorption_reactions : list of dict, optional
            Keyword arguments for mkm.add_reaction() for adsorption/desorption
            steps that are NOT in the energy landscape (e.g. gas-phase kinetics).
        n_electrons : dict
            {observable: n_e} for current conversion.
        A_site : float
            Site area [m^2].
        extra_params : dict, optional
            Additional parameters (p_N2, etc.) merged into mkm.parameters.

        Returns
        -------
        PolarizationCurve
        """
        mkm = MicroKineticModel()
        for sp in species_list:
            mkm.add_species(sp)

        mkm.parameters = {'T': 300.0, 'U': 0.0}
        if extra_params:
            mkm.parameters.update(extra_params)

        for rxn in reactions:
            key = (rxn['is'], rxn['fs'])
            rate_fwd = landscape.barriers_fwd[key].make_rate_func()
            rate_rev = (landscape.barriers_rev[key].make_rate_func()
                        if rxn.get('reversible', True) else None)
            mkm.add_reaction(
                name=rxn['name'],
                reactants=rxn['reactants'],
                products=rxn['products'],
                rate_fwd=rate_fwd,
                rate_rev=rate_rev,
                tof_count=rxn.get('tof_count', {}),
            )

        if adsorption_reactions:
            for ads_kw in adsorption_reactions:
                mkm.add_reaction(**ads_kw)

        A = A_site if A_site is not None else (3.0e-10) ** 2
        ne = n_electrons if n_electrons is not None else {}
        return PolarizationCurve(mkm, ne, A)


def compute_polarization_curve(mkm, U_range, T, n_electrons, A_site,
                               verbose=True):
    """
    Convenience function: compute polarization curve in one call.

    Parameters
    ----------
    mkm : MicroKineticModel
        Model with U-dependent rate functions.
    U_range : array-like
        Potentials [V vs RHE].
    T : float
        Temperature [K].
    n_electrons : dict
        {observable: n_e_per_event}.
    A_site : float
        Site area [m^2].

    Returns
    -------
    dict
        Same as PolarizationCurve.compute().
    """
    pc = PolarizationCurve(mkm, n_electrons, A_site)
    return pc.compute(U_range, T=T, verbose=verbose)


def load_energy_data(filename):
    """
    Load DFT energy data from a JSON file.

    Expected JSON format::

        {
            "potentials": [-2.0, -1.5, -1.0, -0.5, 0.0],
            "state_energies": {
                "clean": [0.0, 0.0, 0.0, 0.0, 0.0],
                "N2*": [-0.8, -0.7, -0.6, -0.5, -0.4],
                ...
            },
            "ts_energies": {
                "N2*->NNH*": [1.8, 2.0, 2.3, 2.6, 3.0],
                ...
            }
        }

    TS keys use "IS->FS" format and are converted to (IS, FS) tuples.

    Returns
    -------
    dict with 'potentials', 'state_energies', 'ts_energies' (tuple keys).
    """
    import json
    with open(filename) as f:
        data = json.load(f)

    # Convert TS keys from "IS->FS" strings to tuples
    ts = {}
    for key, vals in data['ts_energies'].items():
        if '->' in key:
            is_name, fs_name = key.split('->')
            ts[(is_name.strip(), fs_name.strip())] = vals
        else:
            ts[tuple(key.split(','))] = vals

    return {
        'potentials': data['potentials'],
        'state_energies': data['state_energies'],
        'ts_energies': ts,
    }
