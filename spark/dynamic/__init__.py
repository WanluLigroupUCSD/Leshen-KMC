"""
SPARK Dynamic Catalytic KMC — environment-dependent KMC for evolving surfaces.

Implements the V1 catalytic method layer from dynamic-catalytic-KMC:
  - Dynamic surface representation (graph with mutable site types)
  - Local environment descriptors and hashing
  - Unified catalytic + structural event system
  - RateEstimator protocol with lookup/surrogate backends
  - Environment-aware event caching
  - BKL engine with local-only updates

Usage:
    from spark.dynamic import DynamicSurface, DynamicKMCEngine
    from spark.dynamic import LookupTableEstimator

    surface = DynamicSurface.fcc111(composition={'Pd': 0.5, 'Au': 0.5}, size=(20, 20))
    estimator = LookupTableEstimator.from_dict({...})
    engine = DynamicKMCEngine(surface, estimator, temperature=600)
    engine.run(max_steps=100000)
"""

from .surface import DynamicSurface, SiteNode
from .descriptor import LocalEnvironment, EnvHash
from .events import (
    Event, CatalyticEvent, StructuralEvent,
    Adsorption, Desorption, SurfaceReaction, Diffusion,
    SiteConversion, Segregation,
    EventGenerator, EventType,
)
from .rates import RateEstimator, LookupTableEstimator, SurrogateEstimator
from .cache import EventCache
from .engine import DynamicKMCEngine

__all__ = [
    'DynamicSurface', 'SiteNode',
    'LocalEnvironment', 'EnvHash',
    'Event', 'CatalyticEvent', 'StructuralEvent',
    'Adsorption', 'Desorption', 'SurfaceReaction', 'Diffusion',
    'SiteConversion', 'Segregation',
    'EventGenerator', 'EventType',
    'RateEstimator', 'LookupTableEstimator', 'SurrogateEstimator',
    'EventCache',
    'DynamicKMCEngine',
]
