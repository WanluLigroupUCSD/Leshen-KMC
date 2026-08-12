// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::status::{Status,StatusCode};
use serde::{Deserialize,Serialize};

pub const BASIN_IMPLEMENTATION_ENABLED:bool=false;

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasinCheckpoint{pub enabled:bool,pub reason:String}

impl BasinCheckpoint{
    pub fn disabled()->Self{Self{enabled:false,reason:"no internal catalog/discovery completeness validator exists".into()}}
    pub fn validate(&self)->Result<(),Status>{if self.enabled||self.reason!="no internal catalog/discovery completeness validator exists"{return Err(Status::simple(StatusCode::CheckpointCorrupt,"basin","BASIN-006","checkpoint cannot self-attest basin exactness/completeness"));}Ok(())}
}

pub fn attempt_acceleration()->Result<(),Status>{Err(Status::simple(StatusCode::BasinDisabled,"basin","BASIN-006","basin acceleration is hard-disabled until an internal catalog/discovery completeness validator exists"))}
