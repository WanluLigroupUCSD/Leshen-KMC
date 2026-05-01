"""
SPARK - SPatial Atomistic Reaction Kinetics

A unified KMC framework for heterogeneous catalysis combining:

  Lattice KMC (spark.engine):
    - BKL/VSSM rejection-free algorithm
    - Neighbor list with spatial event matching
    - Pairwise lateral interactions (coverage-dependent rates)
    - BEP (Bronsted-Evans-Polanyi) relations
    - Surface diffusion helper
    - Site type support (top, bridge, hollow, etc.)
    - Butler-Volmer electrochemical rates
    - Mean-field microkinetic ODE solver
    - Polarization curve computation
    - Multi-lattice kMC (Hoffmann-Reuter-Scheffler 2015):
        super-lattice with N coexisting commensurate Layers,
        cross-layer "lattice-swap" elementary processes for
        morphological transitions (oxide<->metal, surface
        reconstruction, phase boundary propagation),
        spuck-based 1D site indexing folding layer membership
        into site-within-cell index. See docs/multi_lattice_design.md.

  Off-Lattice KMC (spark.offlattice):
    - On-the-fly saddle point search (dimer method)
    - Local environment detection with 3-step matching
    - Mechanism catalogue with caching and symmetry exploitation
    - Basin / SuperBasin / SuperCache acceleration (bac-MRM)
    - ASE calculator integration (VASP, MACE, EMT, etc.)
    - Full SKMC simulation engine

  Dynamic Catalytic KMC (spark.dynamic):
    - Dynamic surface with mutable site identity (graph/weak-lattice)
    - Environment-dependent rates via RateEstimator protocol
    - Unified catalytic + structural event system
    - Site conversion and segregation events
    - Local-only updates with environment-aware event caching
    - Swappable rate backends: lookup table → ML surrogate → GNN

API style follows kmos conventions for lattice KMC,
ASE conventions for off-lattice KMC,
and graph conventions for dynamic catalytic KMC.
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

__version__ = '0.6.0'

# Off-lattice KMC (lazy import to avoid hard ASE dependency)
from . import offlattice
# Dynamic catalytic KMC
from . import dynamic

__all__ = [
    # Lattice KMC
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
    # Off-lattice KMC
    'offlattice',
    # Dynamic catalytic KMC
    'dynamic',
]
