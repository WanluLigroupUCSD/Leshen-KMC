"""
Rate estimator protocol and implementations.

The RateEstimator is the key interface that decouples the KMC engine from
the rate computation method. The engine never knows how rates are computed.

Swap V1 (lookup table) → V2 (neural network) → V3 (GNN) with zero
engine changes — just provide a different RateEstimator implementation.

Protocol:
    estimator.estimate(event, local_env) → (rate, uncertainty)
"""

import numpy as np
from typing import Protocol, runtime_checkable
from .descriptor import LocalEnvironment
from .events import Event, EventType


@runtime_checkable
class RateEstimator(Protocol):
    """
    Protocol for rate estimation.

    Any object implementing estimate() can be used as a rate estimator.
    This is the strict decoupling point between the KMC engine and the
    barrier model (PPT slide 15).
    """

    def estimate(self, event: Event, env: LocalEnvironment) -> tuple:
        """
        Estimate the rate for an event in a given local environment.

        Parameters
        ----------
        event : Event
            The event to estimate the rate for.
        env : LocalEnvironment
            The local environment at the event's primary site.

        Returns
        -------
        rate : float
            Rate constant in s^-1.
        uncertainty : float
            Estimated uncertainty in the rate (0 = exact).
        """
        ...


class LookupTableEstimator:
    """
    V1 rate estimator: deterministic lookup table.

    Maps (event_type, center_type, adsorbate, nn_composition) → rate.
    For binary alloys with 6 NN, this is at most ~50 entries — exact,
    no ML needed.

    Parameters
    ----------
    temperature : float
        Temperature in Kelvin.
    """

    # Boltzmann constant in eV/K
    KB_EV = 8.617333e-5

    def __init__(self, temperature=300.0):
        self.temperature = temperature
        self._table = {}       # (event_key) → barrier_eV
        self._prefactors = {}  # (event_key) → prefactor_Hz
        self._default_prefactor = 1e13  # TST order of magnitude

    def add_entry(self, event_type, center_type, adsorbate, nn_comp,
                  barrier, prefactor=None):
        """
        Add a rate table entry.

        Parameters
        ----------
        event_type : EventType
            Type of event.
        center_type : str
            Atom type at the site (e.g., 'Pd').
        adsorbate : int
            Adsorbate species ID at center.
        nn_comp : tuple
            Sorted NN composition, e.g. (('Au', 2), ('Pd', 4)).
        barrier : float
            Activation energy in eV.
        prefactor : float, optional
            Pre-exponential factor in Hz. Default: 1e13.
        """
        key = (event_type, center_type, adsorbate, tuple(sorted(nn_comp)))
        self._table[key] = barrier
        if prefactor is not None:
            self._prefactors[key] = prefactor

    def add_entries_from_dict(self, entries):
        """
        Bulk add from a list of dicts.

        Each dict: {
            'event_type': EventType or str,
            'center_type': str,
            'adsorbate': int,
            'nn_comp': tuple,
            'barrier': float,
            'prefactor': float (optional),
        }
        """
        for e in entries:
            et = e['event_type']
            if isinstance(et, str):
                et = EventType[et.upper()]
            self.add_entry(
                event_type=et,
                center_type=e['center_type'],
                adsorbate=e.get('adsorbate', 0),
                nn_comp=e['nn_comp'],
                barrier=e['barrier'],
                prefactor=e.get('prefactor'),
            )

    def estimate(self, event, env):
        """
        Look up rate from table.

        Falls back to a default barrier if the exact (event, env)
        combination is not in the table.

        Returns (rate, uncertainty).
        """
        key = (event.event_type, env.center_type,
               env.center_adsorbate, env.nn_composition)

        if key in self._table:
            barrier = self._table[key]
            prefactor = self._prefactors.get(key, self._default_prefactor)
            rate = prefactor * np.exp(
                -barrier / (self.KB_EV * self.temperature))
            return max(rate, 0.0), 0.0  # exact → zero uncertainty

        # Try without adsorbate (for structural events)
        key_no_ads = (event.event_type, env.center_type, 0,
                      env.nn_composition)
        if key_no_ads in self._table:
            barrier = self._table[key_no_ads]
            prefactor = self._prefactors.get(key_no_ads, self._default_prefactor)
            rate = prefactor * np.exp(
                -barrier / (self.KB_EV * self.temperature))
            return max(rate, 0.0), 0.0

        # Fuzzy match: find entry with same event_type + center_type
        # but closest NN composition (handles variable coordination)
        best_rate = 0.0
        best_dist = float('inf')
        for tkey, barrier in self._table.items():
            if tkey[0] != event.event_type or tkey[1] != env.center_type:
                continue
            # Compare NN composition as dicts
            env_comp = dict(env.nn_composition)
            tbl_comp = dict(tkey[3])
            all_types = set(env_comp.keys()) | set(tbl_comp.keys())
            dist = sum(abs(env_comp.get(t, 0) - tbl_comp.get(t, 0))
                       for t in all_types)
            if dist < best_dist:
                best_dist = dist
                prefactor = self._prefactors.get(tkey, self._default_prefactor)
                best_rate = prefactor * np.exp(
                    -barrier / (self.KB_EV * self.temperature))

        if best_rate > 0:
            return best_rate, 0.3  # moderate uncertainty for fuzzy match

        # Not found — return zero rate (event disabled)
        return 0.0, 1.0  # uncertainty = 1.0 signals "unknown"

    @classmethod
    def from_dict(cls, data, temperature=300.0):
        """
        Create estimator from a structured dict.

        Format:
        {
            'temperature': 600,
            'entries': [
                {'event_type': 'ADSORPTION', 'center_type': 'Pd',
                 'nn_comp': (('Pd', 6),), 'barrier': 0.0},
                ...
            ]
        }
        """
        T = data.get('temperature', temperature)
        est = cls(temperature=T)
        est.add_entries_from_dict(data.get('entries', []))
        return est

    def summary(self):
        """Print table contents."""
        print(f"LookupTableEstimator: {len(self._table)} entries, "
              f"T={self.temperature} K")
        for key, barrier in sorted(self._table.items(),
                                    key=lambda x: x[1]):
            et, ct, ads, nn = key
            nn_str = ','.join(f"{t}:{n}" for t, n in nn)
            ads_str = f"+sp{ads}" if ads else ""
            print(f"  {et.name:<20s} {ct}{ads_str:<6s} NN=[{nn_str}] "
                  f"Ea={barrier:.3f} eV")

    def __repr__(self):
        return (f"LookupTableEstimator(entries={len(self._table)}, "
                f"T={self.temperature}K)")


class SurrogateEstimator:
    """
    V2 rate estimator: ML surrogate model.

    Wraps a trained model (sklearn, PyTorch, etc.) that predicts
    activation barriers from local environment descriptors.

    Parameters
    ----------
    model : object
        Trained model with .predict(X) → barriers.
    feature_builder : callable
        Converts LocalEnvironment → feature vector.
    temperature : float
    prefactor : float
    uncertainty_model : object, optional
        Model for uncertainty estimation (e.g., RF out-of-bag).
    """

    KB_EV = 8.617333e-5

    def __init__(self, model, feature_builder, temperature=300.0,
                 prefactor=1e13, uncertainty_model=None):
        self.model = model
        self.feature_builder = feature_builder
        self.temperature = temperature
        self.prefactor = prefactor
        self.uncertainty_model = uncertainty_model

    def estimate(self, event, env):
        """Predict rate from ML model."""
        features = self.feature_builder(event, env)
        X = np.array(features).reshape(1, -1)

        barrier = float(self.model.predict(X)[0])
        barrier = max(barrier, 0.0)

        rate = self.prefactor * np.exp(
            -barrier / (self.KB_EV * self.temperature))

        uncertainty = 0.0
        if self.uncertainty_model is not None:
            uncertainty = float(self.uncertainty_model.predict(X)[0])

        return max(rate, 0.0), uncertainty

    def __repr__(self):
        return (f"SurrogateEstimator(model={type(self.model).__name__}, "
                f"T={self.temperature}K)")
