"""
Rate expression evaluation for KMC models.

Supports kmos-style rate expressions with automatic variable resolution:
  - Physical constants: kB, h, hbar, eV, angstrom, bar, umass, pi
  - Temperature: beta = 1/(kB*T)
  - Molecular masses: m_<species> (e.g., m_N2, m_CO)
  - Math functions: exp, log, sqrt, sin, cos, pow, abs
  - All user-defined parameters
"""

import numpy as np
from .units import (
    kB, h, hbar, NA, R_gas, e_charge,
    eV, angstrom, bar, atm, umass, pi,
    ATOMIC_MASSES, MOLECULAR_MASSES,
)


def _build_namespace(params):
    """Build evaluation namespace from parameters and physical constants."""
    ns = {
        # Physical constants
        'kB': kB, 'h': h, 'hbar': hbar, 'pi': pi,
        'eV': eV, 'angstrom': angstrom, 'bar': bar, 'atm': atm,
        'umass': umass, 'NA': NA, 'R': R_gas, 'e': e_charge,
        # Math functions
        'exp': np.exp, 'log': np.log, 'ln': np.log,
        'log10': np.log10, 'sqrt': np.sqrt,
        'sin': np.sin, 'cos': np.cos, 'pow': pow, 'abs': abs,
    }

    # Add molecular masses as m_<species>
    for sp, mass in ATOMIC_MASSES.items():
        ns[f'm_{sp}'] = mass
    for sp, mass in MOLECULAR_MASSES.items():
        ns[f'm_{sp}'] = mass

    # Resolve parameters (some may be expressions themselves)
    resolved = {}
    for name, value in params.items():
        if isinstance(value, str):
            try:
                resolved[name] = eval(value, {"__builtins__": {}}, ns)
            except Exception:
                resolved[name] = value
        else:
            resolved[name] = value
    ns.update(resolved)

    # Compute beta = 1/(kB*T) if T is available
    if 'T' in ns and isinstance(ns['T'], (int, float)):
        ns['beta'] = 1.0 / (kB * float(ns['T']))

    return ns


def evaluate_rate_expression(expr, params):
    """
    Evaluate a rate constant expression.

    Parameters
    ----------
    expr : str, int, or float
        Rate expression (e.g., 'p_CO*bar*A/sqrt(2*pi*m_CO*umass/beta)')
    params : dict
        Parameter name -> value mapping

    Returns
    -------
    float
        Evaluated rate constant in s^-1
    """
    if isinstance(expr, (int, float)):
        return float(expr)

    ns = _build_namespace(params)

    try:
        result = float(eval(str(expr), {"__builtins__": {}}, ns))
        if np.isnan(result) or np.isinf(result):
            raise ValueError(f"Rate expression '{expr}' evaluated to {result}")
        return max(result, 0.0)
    except Exception as e:
        raise ValueError(f"Cannot evaluate rate expression '{expr}': {e}")


def arrhenius(A, Ea, T):
    """
    Arrhenius rate constant: k = A * exp(-Ea / (kB * T))

    Parameters
    ----------
    A : float
        Pre-exponential factor [s^-1]
    Ea : float
        Activation energy [eV]
    T : float
        Temperature [K]

    Returns
    -------
    float
        Rate constant [s^-1]
    """
    return A * np.exp(-Ea * eV / (kB * T))


def tst_rate(Ea, T):
    """
    Transition state theory rate: k = (kB*T/h) * exp(-Ea / (kB*T))

    Parameters
    ----------
    Ea : float
        Activation energy [eV]
    T : float
        Temperature [K]

    Returns
    -------
    float
        Rate constant [s^-1]
    """
    return (kB * T / h) * np.exp(-Ea * eV / (kB * T))


def hertz_knudsen(p, T, m, A_site):
    """
    Hertz-Knudsen impingement rate for gas adsorption.

    k_ads = p * A_site / sqrt(2 * pi * m * kB * T)

    Parameters
    ----------
    p : float
        Partial pressure [Pa]
    T : float
        Temperature [K]
    m : float
        Molecular mass [kg] (use m_amu * umass)
    A_site : float
        Site area [m^2]

    Returns
    -------
    float
        Adsorption rate constant [s^-1]
    """
    return p * A_site / np.sqrt(2 * pi * m * kB * T)


def interpolated_tst_rate(U_data, Ea_data, kind='cubic'):
    """
    Create a TST rate function with Ea(U) interpolated from DFT data.

    Returns a callable rate_func(params) compatible with MicroKineticModel,
    where params must contain 'T' and 'U'.

    Parameters
    ----------
    U_data : array-like
        Potentials [V] at which barriers were computed by DFT.
    Ea_data : array-like
        Activation energies [eV] at each potential.
    kind : str
        Interpolation method ('linear', 'cubic').

    Returns
    -------
    callable
        rate_func(params) -> float  [s^-1]
    """
    from scipy.interpolate import interp1d

    U_arr = np.asarray(U_data, dtype=float)
    Ea_arr = np.asarray(Ea_data, dtype=float)
    if len(U_arr) <= 3 and kind == 'cubic':
        kind = 'quadratic' if len(U_arr) == 3 else 'linear'

    interp = interp1d(U_arr, Ea_arr, kind=kind,
                       fill_value='extrapolate', bounds_error=False)

    def rate_func(params):
        Ea = max(float(interp(params['U'])), 0.0)
        return tst_rate(Ea, params['T'])

    return rate_func


def electrochemical_rate(Ea, T, U, U0=0.0, beta_bv=0.5):
    """
    Electrochemical (Butler-Volmer) rate for proton-coupled electron transfer.

    k = (kB*T/h) * exp(-(Ea + beta_bv * e * (U - U0)) / (kB*T))

    Parameters
    ----------
    Ea : float
        Activation energy at equilibrium potential [eV]
    T : float
        Temperature [K]
    U : float
        Applied potential [V vs RHE]
    U0 : float
        Equilibrium potential [V vs RHE]
    beta_bv : float
        Symmetry factor (default 0.5)

    Returns
    -------
    float
        Rate constant [s^-1]
    """
    effective_Ea = Ea + beta_bv * (U - U0)  # in eV
    if effective_Ea < 0:
        effective_Ea = 0.0
    return (kB * T / h) * np.exp(-effective_Ea * eV / (kB * T))
