// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::status::{Status, StatusCode};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use crate::identity::{canonical_state_id,IdentityConfig};

pub type Vec3 = [f64; 3];
pub type Cell = [[f64; 3]; 3];

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElectronicState {
    pub charge_e: i32,
    pub spin_2s: i32,
    pub multiplicity: u32,
    #[serde(default)]
    pub metadata: BTreeMap<String, serde_json::Value>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtomicSystem {
    pub schema: String,
    pub atom_ids: Vec<String>,
    pub species: Vec<String>,
    pub positions_angstrom: Vec<Vec3>,
    pub cell_angstrom: Cell,
    pub periodic: [bool; 3],
    pub movable: Vec<[bool; 3]>,
    pub electronic: ElectronicState,
    pub calculator_model_digest: String,
}

impl AtomicSystem {
    pub fn validate(&self) -> Result<(), Status> {
        let n = self.atom_ids.len();
        if self.schema != crate::IR_SCHEMA {
            return Err(Status::simple(StatusCode::SchemaUnsupported, "state", "STATE-001", "unsupported state schema"));
        }
        if n == 0 || self.species.len() != n || self.positions_angstrom.len() != n || self.movable.len() != n {
            return Err(Status::simple(StatusCode::InvalidState, "state", "STATE-001", "atom arrays have inconsistent lengths"));
        }
        let mut ids = self.atom_ids.clone();
        ids.sort();
        if ids.windows(2).any(|w| w[0] == w[1]) {
            return Err(Status::simple(StatusCode::InvalidState, "state", "STATE-001", "atom IDs are not unique"));
        }
        if self.species.iter().any(|x| x.is_empty()) || self.calculator_model_digest.is_empty() {
            return Err(Status::simple(StatusCode::InvalidState, "state", "STATE-001", "empty species or model digest"));
        }
        if self.electronic.multiplicity == 0
            || self.electronic.multiplicity != self.electronic.spin_2s.unsigned_abs() + 1
        {
            return Err(Status::simple(StatusCode::InvalidState, "state", "STATE-001", "spin and multiplicity are inconsistent"));
        }
        if !self.positions_angstrom.iter().flatten().chain(self.cell_angstrom.iter().flatten()).all(|x| x.is_finite()) {
            return Err(Status::simple(StatusCode::NonfiniteResult, "state", "STATE-008", "nonfinite position or cell"));
        }
        let det = determinant(self.cell_angstrom);
        if !det.is_finite() || det.abs() <= 1e-18 {
            return Err(Status::simple(StatusCode::InvalidState, "state", "STATE-001", "singular simulation cell"));
        }
        Ok(())
    }

    pub fn validate_fixed_against(&self, reference: &Self) -> Result<(), Status> {
        self.validate()?;
        if self.atom_ids.len() != reference.atom_ids.len() {
            return Err(Status::simple(StatusCode::AtomCountChangeUnsupported, "state", "ATOM-001", "atom count changed"));
        }
        let binding=|s:&AtomicSystem|{let mut m=BTreeMap::new();for i in 0..s.atom_ids.len(){m.insert(s.atom_ids[i].clone(),(s.species[i].clone(),s.movable[i]));}m};
        if binding(self)!=binding(reference){return Err(Status::simple(StatusCode::AtomCountChangeUnsupported,"state","CALC-003","atom ID set or ID-to-species/constraint binding changed"));}
        if self.cell_angstrom != reference.cell_angstrom || self.periodic != reference.periodic
            || self.electronic != reference.electronic || self.calculator_model_digest != reference.calculator_model_digest
        {
            return Err(Status::simple(StatusCode::InvalidState, "state", "SCOPE-003", "fixed cell, periodicity, electronic state, or model digest changed"));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RelaxationProvenance {
    pub minimizer_name: String,
    pub minimizer_version: String,
    pub minimizer_callback_identity:String,
    pub calculator_callback_identity:String,
    pub calculator_evaluations: u64,
    pub steps: u64,
    pub termination_reason: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CommittedState {
    pub(crate) state_id: String,
    pub(crate) system: AtomicSystem,
    pub(crate) energy_ev: f64,
    pub(crate) forces_ev_per_angstrom: Vec<Vec3>,
    pub(crate) max_movable_force_ev_per_angstrom: f64,
    pub(crate) force_tolerance_ev_per_angstrom: f64,
    pub(crate) constraint_digest: String,
    pub(crate) relaxation: RelaxationProvenance,
    pub(crate) identity_config:IdentityConfig,
    pub(crate) identity_config_digest:String,
}

impl CommittedState {
    pub(crate) fn try_new(system:AtomicSystem,energy_ev:f64,forces_ev_per_angstrom:Vec<Vec3>,force_tolerance_ev_per_angstrom:f64,constraint_digest:String,relaxation:RelaxationProvenance,identity_config:IdentityConfig)->Result<Self,Status>{system.validate()?;identity_config.validate()?;let max_movable_force_ev_per_angstrom=recompute_max_force(&system,&forces_ev_per_angstrom)?;let identity_config_digest=identity_config.digest()?;let state_id=canonical_state_id(&system,&identity_config,&constraint_digest)?;let s=Self{state_id,system,energy_ev,forces_ev_per_angstrom,max_movable_force_ev_per_angstrom,force_tolerance_ev_per_angstrom,constraint_digest,relaxation,identity_config,identity_config_digest};s.validate()?;Ok(s)}
    pub fn state_id(&self)->&str{&self.state_id}pub fn system(&self)->&AtomicSystem{&self.system}pub fn energy_ev(&self)->f64{self.energy_ev}pub fn forces(&self)->&[Vec3]{&self.forces_ev_per_angstrom}pub fn max_movable_force(&self)->f64{self.max_movable_force_ev_per_angstrom}pub fn identity_config(&self)->&IdentityConfig{&self.identity_config}
    pub fn validate(&self) -> Result<(), Status> {
        self.system.validate()?;
        if self.state_id.is_empty() || self.constraint_digest.is_empty()
            || self.forces_ev_per_angstrom.len() != self.system.atom_ids.len()
            || !self.energy_ev.is_finite() || !self.max_movable_force_ev_per_angstrom.is_finite()
            || !self.force_tolerance_ev_per_angstrom.is_finite() || self.force_tolerance_ev_per_angstrom < 0.0
            || !self.forces_ev_per_angstrom.iter().flatten().all(|x| x.is_finite())
        {
            return Err(Status::simple(StatusCode::NonfiniteResult, "state", "STATE-003", "invalid committed-state numeric data"));
        }
        if self.max_movable_force_ev_per_angstrom > self.force_tolerance_ev_per_angstrom {
            return Err(Status::simple(StatusCode::RelaxNotConverged, "state", "STATE-004", "maximum movable force exceeds tolerance"));
        }
        if recompute_max_force(&self.system,&self.forces_ev_per_angstrom)?.to_bits()!=self.max_movable_force_ev_per_angstrom.to_bits()||self.identity_config.digest()?!=self.identity_config_digest||canonical_state_id(&self.system,&self.identity_config,&self.constraint_digest)?!=self.state_id||self.relaxation.minimizer_name.is_empty()||self.relaxation.minimizer_version.is_empty()||self.relaxation.minimizer_callback_identity.is_empty()||self.relaxation.calculator_callback_identity.is_empty()||self.relaxation.termination_reason.is_empty(){return Err(Status::simple(StatusCode::InvalidState,"state","STATE-003","committed-state force/provenance/identity evidence mismatch"));}
        Ok(())
    }
}

fn recompute_max_force(system:&AtomicSystem,forces:&[Vec3])->Result<f64,Status>{if forces.len()!=system.atom_ids.len()||!forces.iter().flatten().all(|x|x.is_finite()){return Err(Status::simple(StatusCode::NonfiniteResult,"state","STATE-008","invalid full-force array"));}let mut max=0.0_f64;for(i,f)in forces.iter().enumerate(){let norm=(0..3).filter(|&k|system.movable[i][k]).map(|k|f[k]*f[k]).sum::<f64>().sqrt();max=max.max(norm);}Ok(max)}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SaddleRecord {
    pub geometry: AtomicSystem,
    pub energy_ev: f64,
    pub forces_ev_per_angstrom: Vec<Vec3>,
    pub unstable_direction: Vec<Vec3>,
    pub curvature_ev_per_angstrom2: f64,
    pub evidence_level: SaddleEvidence,
    pub imaginary_mode_count_after_exclusions: Option<u32>,
    pub rigid_constrained_modes_excluded: bool,
    pub orthogonal_curvatures_ev_per_angstrom2: Vec<f64>,
    pub search_id: String,
    pub evaluation_count: u64,
    pub termination_reason: String,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum SaddleEvidence { Hessian, Directional }

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ValidationEvidence {
    pub method: SaddleEvidence,
    pub origin_match: GeometryMatch,
    pub destination_match: GeometryMatch,
    pub unstable_mode_count: u32,
    pub full_endpoint_relaxations: bool,
    pub calculator_model_digest: String,
    pub constraint_digest: String,
    pub search_id:String,
    pub rng_substream_digest:String,
    pub calculator_callback_identity:String,
    pub minimizer_callback_identity:String,
    pub endpoint_states:[CommittedState;2],
    pub endpoint_receipt_digests:[String;2],
    pub(crate) saddle_search:crate::callbacks::ValidatedSaddle,
    pub(crate) endpoint_relaxations:[crate::callbacks::ValidatedRelaxation;2],
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GeometryMatch {
    pub mapping_a_to_b: Vec<usize>,
    pub translation_angstrom: Vec3,
    pub rms_displacement_angstrom: f64,
    pub max_displacement_angstrom: f64,
    pub energy_difference_ev: f64,
    pub identity_version: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RateModelRecord {
    pub model: String,
    pub temperature_k: f64,
    pub common_prefactor_per_s: Option<f64>,
    pub log_forward_rate_per_s: f64,
    pub log_reverse_rate_per_s: f64,
    pub detailed_balance_residual: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DirectedEvent {
    pub event_id: String,
    pub reverse_pair_id: String,
    pub origin_state_id: String,
    pub destination_state_id: String,
    pub saddle: SaddleRecord,
    pub barrier_ev: f64,
    pub reverse_barrier_ev: f64,
    pub rate_model: RateModelRecord,
    pub active_atom_mapping: Vec<[usize; 2]>,
    pub environment_key: String,
    pub environment_version: String,
    pub discovery_provenance: BTreeMap<String, serde_json::Value>,
    pub validation: ValidationEvidence,
    pub calculator_digest: String,
    pub identity_digest: String,
    pub tolerance_digest: String,
    pub schema_digest: String,
}

pub fn determinant(m: Cell) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
