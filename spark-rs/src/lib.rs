//! SPARK Rust acceleration crate.
//!
//! Provides high-performance kernels callable from Python via PyO3 (when the
//! `python` feature is enabled), plus a standalone CLI binary (`src/main.rs`).
//! Off-lattice OTF KMC hot-loop components live in `offlattice/`.
//!
//! Build modes:
//!   - `cargo build --bin spark`     — pure CLI binary, no Python dep.
//!   - `maturin build` (in spark-rs/) — Python wheel exposing `spark_rs._native`.
//!     `pyproject.toml` selects `--features python` automatically.

pub mod units;
pub mod model;
pub mod rates;
pub mod engine;
pub mod microkinetic;
pub mod analysis;
pub mod models;
pub mod polarization;
pub mod offlattice;

#[cfg(feature = "python")]
pub mod python_bindings;

// ---------------------------------------------------------------------------
// PyO3 entry — exposes a Python module named `spark_rs._native`.
// All Python-facing code is gated behind the `python` feature so that
// `cargo build --bin spark` (which does not need pyo3) keeps working.
// ---------------------------------------------------------------------------

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Smoke-test entry: returns a string so Python can verify the .so is loaded.
#[cfg(feature = "python")]
#[pyfunction]
fn hello() -> &'static str {
    "spark-rs is alive"
}

/// Crate version, parsed from Cargo.toml at compile time.
#[cfg(feature = "python")]
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(feature = "python")]
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hello, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    // Phase C — OTF KMC hot-loop kernels with Python force callback.
    m.add_function(wrap_pyfunction!(python_bindings::dimer_find_saddle, m)?)?;
    m.add_function(wrap_pyfunction!(python_bindings::fire_minimize, m)?)?;
    Ok(())
}
