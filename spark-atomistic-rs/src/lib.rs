// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
#![forbid(unsafe_code)]

pub mod basin;
pub mod callbacks;
pub mod catalog;
pub mod checkpoint;
pub mod config;
pub mod discovery;
pub mod identity;
pub mod kmc;
pub mod model;
pub mod parity;
pub mod rate;
pub mod resource;
pub mod rng;
pub mod run;
pub mod status;

pub const IR_SCHEMA: &str = "spark-atomistic-model/1";
pub const BACKEND: &str = "rust";
pub const BASE_SPEC_SHA256: &str =
    "8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84";
pub const ERRATA_1_SHA256: &str =
    "52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40";
pub const ERRATA_2_SHA256: &str =
    "eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReleaseClaims {
    pub validated: bool,
    pub production: bool,
    pub release: bool,
}

pub const RELEASE_CLAIMS: ReleaseClaims = ReleaseClaims {
    validated: false,
    production: false,
    release: false,
};
