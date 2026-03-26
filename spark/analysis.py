"""
Analysis tools for KMC simulations.

Provides trajectory recording, TOF analysis, coverage tracking,
and steady-state detection.
"""

import numpy as np


class TrajectoryRecorder:
    """
    Records simulation trajectory at regular intervals for post-analysis.

    Usage:
        recorder = TrajectoryRecorder(engine)
        for _ in range(100):
            engine.do_steps(10000)
            recorder.record()

        times = recorder.get_times()
        cov_N = recorder.get_coverage_array('N')
        tof_NH3 = recorder.get_tof_array('NH3_production')
    """

    def __init__(self, engine):
        self.engine = engine
        self.times = []
        self.steps = []
        self.coverages = []
        self.tofs = []
        self.process_stats = []

    def record(self):
        """Record current state snapshot."""
        self.times.append(self.engine.kmc_time)
        self.steps.append(self.engine.kmc_step)
        self.coverages.append(self.engine.get_coverage())
        tof = self.engine.get_tof()
        if tof:
            self.tofs.append(tof)
        self.process_stats.append(self.engine.procstat.copy())

    def get_times(self):
        return np.array(self.times)

    def get_steps(self):
        return np.array(self.steps)

    def get_coverage_array(self, species_name):
        """Get time series of coverage for a species."""
        return np.array([c.get(species_name, 0.0) for c in self.coverages])

    def get_tof_array(self, observable):
        """Get time series of TOF for an observable."""
        return np.array([t.get(observable, 0.0) for t in self.tofs])

    def get_all_coverages(self):
        """Get coverage time series as dict of arrays."""
        if not self.coverages:
            return {}
        species = list(self.coverages[0].keys())
        return {sp: self.get_coverage_array(sp) for sp in species}

    def get_all_tofs(self):
        """Get TOF time series as dict of arrays."""
        if not self.tofs:
            return {}
        observables = set()
        for t in self.tofs:
            observables.update(t.keys())
        return {obs: self.get_tof_array(obs) for obs in observables}


def check_steady_state(values, window=50, rtol=0.05):
    """
    Check if a time series has reached steady state.

    Compares means of two consecutive windows. Returns True if
    the relative change is below rtol.
    """
    if len(values) < 2 * window:
        return False
    recent = values[-window:]
    earlier = values[-2 * window:-window]
    mean_r = np.mean(recent)
    mean_e = np.mean(earlier)
    if mean_e == 0 and mean_r == 0:
        return True
    if abs(mean_e) < 1e-30:
        return abs(mean_r) < 1e-30
    return abs(mean_r - mean_e) / abs(mean_e) < rtol


def run_to_steady_state(engine, max_steps=1e8, sample_interval=10000,
                        rtol=0.05, min_samples=100, verbose=True):
    """
    Run simulation until steady state is reached.

    Parameters
    ----------
    engine : KMCEngine
        The KMC simulation engine.
    max_steps : int
        Maximum number of KMC steps.
    sample_interval : int
        Steps between samples.
    rtol : float
        Relative tolerance for steady-state detection.
    min_samples : int
        Minimum number of samples before checking.
    verbose : bool
        Print progress.

    Returns
    -------
    TrajectoryRecorder
        Complete trajectory recording.
    """
    recorder = TrajectoryRecorder(engine)
    total_steps = 0
    sample_interval = int(sample_interval)

    while total_steps < max_steps:
        engine.do_steps(sample_interval)
        total_steps += sample_interval
        recorder.record()

        if len(recorder.coverages) >= min_samples:
            cov = engine.get_coverage()
            all_steady = True
            for sp_name, val in cov.items():
                if val > 0.01:
                    series = recorder.get_coverage_array(sp_name)
                    if not check_steady_state(series, rtol=rtol):
                        all_steady = False
                        break
            if all_steady:
                if verbose:
                    print(f"Steady state reached: step={engine.kmc_step}, "
                          f"time={engine.kmc_time:.6e} s")
                return recorder

    if verbose:
        print(f"Warning: steady state not reached after "
              f"{int(max_steps)} steps")
    return recorder


def compute_selectivity(tof_dict, target, total_products=None):
    """
    Compute selectivity toward a target product.

    Parameters
    ----------
    tof_dict : dict
        TOF values for each product.
    target : str
        Name of the target product.
    total_products : list of str, optional
        List of product names. If None, uses all keys in tof_dict.

    Returns
    -------
    float
        Selectivity (0 to 1).
    """
    if total_products is None:
        total_products = list(tof_dict.keys())

    total_tof = sum(abs(tof_dict.get(p, 0.0)) for p in total_products)
    if total_tof == 0:
        return 0.0
    return abs(tof_dict.get(target, 0.0)) / total_tof


def apparent_activation_energy(temperatures, tofs, R_eV=8.617333e-5):
    """
    Calculate apparent activation energy from Arrhenius plot.

    Parameters
    ----------
    temperatures : array-like
        Temperatures [K].
    tofs : array-like
        TOF values at each temperature.
    R_eV : float
        Boltzmann constant in eV/K.

    Returns
    -------
    float
        Apparent activation energy [eV].
    """
    T = np.array(temperatures, dtype=float)
    r = np.array(tofs, dtype=float)
    mask = r > 0
    if np.sum(mask) < 2:
        return 0.0
    x = 1.0 / T[mask]
    y = np.log(r[mask])
    slope, _ = np.polyfit(x, y, 1)
    return -slope * R_eV
