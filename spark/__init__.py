"""
SPARK - SPatial Atomistic Reaction Kinetics

A pure Python implementation of lattice KMC (BKL/rejection-free algorithm)
and mean-field microkinetic ODE solver for heterogeneous catalysis.

Features:
  - BKL/VSSM rejection-free KMC algorithm
  - Neighbor list with spatial event matching
  - Pairwise lateral interactions (coverage-dependent rates)
  - BEP (Bronsted-Evans-Polanyi) relations
  - Surface diffusion helper
  - Site type support (top, bridge, hollow, etc.)
  - Butler-Volmer electrochemical rates
  - Mean-field microkinetic ODE solver
  - Polarization curve computation

API style follows kmos conventions.
"""

from .types import (
    Project, Species, Site, Layer, Lattice, Coord,
    Condition, Action, Parameter, Process,
    LateralInteraction, BEPRelation,
)
from .engine import KMCEngine
from .microkinetic import MicroKineticModel
from .analysis import TrajectoryRecorder, run_to_steady_state
from .polarization import (
    InterpolatedBarrier, EnergyLandscape,
    PolarizationCurve, compute_polarization_curve,
    tof_to_current_density, load_energy_data,
)
from .rates import (
    bep_modified_rate, lateral_modified_rate,
    arrhenius, tst_rate, hertz_knudsen, electrochemical_rate,
)
from .io import (
    load_model, load_sparkin, load_yaml, loads_sparkin, loads_yaml,
    project_to_sparkin, project_to_yaml,
    load_spark, loads_spark, project_to_spark,  # aliases
)

__version__ = '0.4.0'
__all__ = [
    'Project', 'Species', 'Site', 'Layer', 'Lattice', 'Coord',
    'Condition', 'Action', 'Parameter', 'Process',
    'LateralInteraction', 'BEPRelation',
    'KMCEngine', 'MicroKineticModel',
    'TrajectoryRecorder', 'run_to_steady_state',
    'InterpolatedBarrier', 'EnergyLandscape',
    'PolarizationCurve', 'compute_polarization_curve',
    'tof_to_current_density', 'load_energy_data',
    'bep_modified_rate', 'lateral_modified_rate',
    'arrhenius', 'tst_rate', 'hertz_knudsen', 'electrochemical_rate',
    'load_model', 'load_sparkin', 'load_yaml',
    'loads_sparkin', 'loads_yaml',
    'project_to_sparkin', 'project_to_yaml',
]
