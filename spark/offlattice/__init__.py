"""
SPARK Off-Lattice KMC — On-the-fly off-lattice kinetic Monte Carlo simulation.

Ported from openFLY (C++) into the SPARK framework. Provides:
  - Continuous (off-lattice) atomic system representation via ASE Atoms
  - Local environment detection with 3-step matching (hash → fingerprint → geometry)
  - On-the-fly saddle point search (dimer method)
  - Mechanism catalogue with caching and symmetry exploitation
  - Basin / SuperBasin / SuperCache acceleration (bac-MRM)
  - Full SKMC simulation engine

Usage:
    from spark.offlattice import SKMCEngine, Catalogue, Mechanism
    from ase.build import bulk

    atoms = bulk('Fe', 'bcc', a=2.87) * (4, 4, 4)
    atoms.calc = ...  # Any ASE calculator

    engine = SKMCEngine(atoms, temperature=500.0)
    engine.run(max_steps=1000, callback=my_callback)
"""

from .mechanism import Mechanism
from .environment import Geometry, Fingerprint
from .catalogue import Catalogue
from .minimize import Minimizer
from .saddle import DimerSearch, SaddleMaster
from .basin import Basin
from .superbasin import SuperBasin
from .cache import SuperCache
from .engine import SKMCEngine
from .potential import CalculatorAdapter

__all__ = [
    'Mechanism',
    'Geometry', 'Fingerprint',
    'Catalogue',
    'Minimizer',
    'DimerSearch', 'SaddleMaster',
    'Basin',
    'SuperBasin',
    'SuperCache',
    'SKMCEngine',
    'CalculatorAdapter',
]
