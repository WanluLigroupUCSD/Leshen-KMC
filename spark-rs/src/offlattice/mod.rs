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
