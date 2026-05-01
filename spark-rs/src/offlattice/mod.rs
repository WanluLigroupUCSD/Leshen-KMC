//! Off-lattice on-the-fly kinetic Monte Carlo simulation module.
//!
//! Ported from openFLY (C++) and the Python spark.offlattice package.
//! Provides high-performance implementations of:
//!   - Local environment detection and fingerprinting
//!   - Mechanism storage and catalogue matching
//!   - Basin / SuperBasin / SuperCache acceleration
//!   - SKMC simulation engine

pub mod mechanism;
pub mod environment;
pub mod catalogue;
pub mod basin;
pub mod superbasin;
pub mod cache;
pub mod engine;

// Phase B (2026-05-01): Rust port of the OTF KMC hot loop.
// `calc` defines the Calculator trait + Müller-Brown reference PES.
// `minimize` provides FIRE; `saddle` provides the dimer search.
// All three are pure Rust (no Python dep) so they can be unit-tested via
// `cargo test --lib` without a Python interpreter.
pub mod calc;
pub mod minimize;
pub mod saddle;
