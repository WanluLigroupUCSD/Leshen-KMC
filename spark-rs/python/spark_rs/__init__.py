"""SPARK Rust acceleration kernels.

This package wraps the spark-rs Rust crate via PyO3. It exposes high-performance
kernels for the off-lattice OTF KMC hot loop (dimer saddle search, local
minimization, environment matching), with Python ASE calculators preserved as
the force-evaluation backend via callbacks.

Usage:
    >>> from spark_rs import hello, version
    >>> version()
    '0.4.0'
    >>> hello()
    'spark-rs is alive'
"""

from spark_rs._native import (
    hello,
    version,
)

__all__ = ["hello", "version"]
