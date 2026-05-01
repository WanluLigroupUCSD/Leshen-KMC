"""SPARK Rust acceleration kernels (PyO3 bindings).

Exposes high-performance Rust implementations of the off-lattice OTF KMC
hot-loop:
  * :func:`dimer_find_saddle` — Henkelman-Jónsson dimer saddle search.
  * :func:`fire_minimize`     — Bitzek FIRE minimizer.

Both take a ``force_callback`` Python callable
``(positions: np.ndarray) -> (energy: float, forces: np.ndarray)`` so any
ASE-compatible calculator (EMT, MACE, VASP, ...) can drive Rust without
the inner loop crossing the Python boundary on arithmetic.

Example — dimer saddle from a quadratic test PES:

    >>> import numpy as np
    >>> from spark_rs import dimer_find_saddle
    >>> def f(p):
    ...     # V = -x^2/2 + y^2/2 → forces = (x, -y)
    ...     return -0.5*p[0,0]**2 + 0.5*p[0,1]**2, np.array([[p[0,0], -p[0,1], 0.0]])
    >>> r = dimer_find_saddle(np.array([[0.3, 0.2, 0.0]]),
    ...                       np.array([[0.9, 0.4, 0.0]]) / np.linalg.norm([0.9, 0.4, 0.0]),
    ...                       f, f_tol=1e-4)
    >>> r['status']
    'success'
"""

from spark_rs._native import (
    hello,
    version,
    dimer_find_saddle,
    fire_minimize,
)

__all__ = [
    "hello",
    "version",
    "dimer_find_saddle",
    "fire_minimize",
]
