// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::model::{AtomicSystem,CommittedState,SaddleEvidence,SaddleRecord,Vec3};
use crate::resource::{CancelToken,ResourceLedger};
use crate::rng::{validate_state,substream_digest,Philox,PhiloxState};
use crate::status::{Status,StatusCode};
use serde::{Deserialize,Serialize};
use std::time::Instant;
use crate::identity::IdentityConfig;

#[derive(Clone,Debug,Eq,PartialEq,Serialize,Deserialize)]
#[serde(rename_all="SCREAMING_SNAKE_CASE")]
pub enum Property{Energy,Forces,Stress,Hessian,StableModeLogs,Uncertainty}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CalculatorRequest{pub state:AtomicSystem,pub requested_properties:Vec<Property>,pub request_digest:String}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CalculatorResponse{pub energy_ev:f64,pub forces_ev_per_angstrom:Vec<Vec3>,pub units_energy:String,pub units_force:String,pub model_name:String,pub model_version:String,pub model_digest:String,pub evaluation_id:String,pub deterministic:bool}

pub trait Calculator:Send+Sync{fn identity(&self)->&str;fn evaluate(&self,request:&CalculatorRequest,deadline:Instant,cancel:&CancelToken)->Result<CalculatorResponse,Status>;}

pub fn evaluate_checked(calc:&dyn Calculator,request:&CalculatorRequest,ledger:&mut ResourceLedger,cancel:&CancelToken,
    relaxation_id:Option<&str>,saddle_id:Option<&str>)->Result<CalculatorResponse,Status>{
    request.state.validate()?;
    if !request.requested_properties.contains(&Property::Energy)||!request.requested_properties.contains(&Property::Forces){return Err(calc_fail("energy and full forces were not requested"));}
    let deadline=ledger.callback_deadline()?;ledger.reserve_calculator(1,relaxation_id,saddle_id)?;
    if cancel.is_cancelled(){return Err(Status::simple(StatusCode::Cancelled,"calculator","CKPT-005","cancelled before callback"));}
    let response=calc.evaluate(request,deadline,cancel).map_err(normalize_calc_error)?;
    if Instant::now()>deadline{return Err(calc_fail("calculator deadline exceeded"));}
    validate_calculator_response(request,&response)?;Ok(response)
}

pub fn validate_calculator_response(req:&CalculatorRequest,r:&CalculatorResponse)->Result<(),Status>{
    if r.units_energy!="eV"||r.units_force!="eV/angstrom"||r.model_digest!=req.state.calculator_model_digest{return Err(calc_fail("unit or calculator model digest mismatch"));}
    if r.model_name.is_empty()||r.model_version.is_empty()||r.evaluation_id.is_empty()||r.forces_ev_per_angstrom.len()!=req.state.atom_ids.len(){return Err(calc_fail("missing metadata or malformed full-force shape"));}
    if !r.energy_ev.is_finite()||!r.forces_ev_per_angstrom.iter().flatten().all(|x|x.is_finite()){return Err(Status::simple(StatusCode::NonfiniteResult,"calculator","STATE-008","nonfinite calculator result"));}Ok(())
}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RelaxRequest{pub state:AtomicSystem,pub identity_config:IdentityConfig,pub force_tolerance_ev_per_angstrom:f64,pub step_limit:u64,pub evaluation_limit:u64,pub constraint_digest:String,pub calculator_identity:String,pub request_id:String}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RelaxResponse{pub state:AtomicSystem,pub energy_ev:f64,pub forces_ev_per_angstrom:Vec<Vec3>,pub max_movable_force_ev_per_angstrom:f64,pub steps:u64,pub calculator_evaluations:u64,pub termination_reason:String,pub minimizer_name:String,pub minimizer_version:String}
pub struct CallbackAttempt<T>{pub result:Result<T,Status>,pub calculator_evaluations:u64,pub termination_reason:String}
pub trait Minimizer:Send+Sync{fn identity(&self)->&str;fn relax(&self,request:&RelaxRequest,calculator:&dyn Calculator,deadline:Instant,cancel:&CancelToken)->CallbackAttempt<RelaxResponse>;}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ValidatedRelaxation{pub(crate)request:RelaxRequest,pub(crate)state:CommittedState,pub(crate)minimizer_identity:String,pub(crate)calculator_identity:String,pub(crate)receipt_digest:String}
impl ValidatedRelaxation{pub fn state(&self)->&CommittedState{&self.state}pub(crate)fn validate_receipt(&self)->Result<(),Status>{validate_relax_request(&self.request)?;self.state.validate()?;if self.minimizer_identity!=self.state.relaxation.minimizer_callback_identity||self.calculator_identity!=self.request.calculator_identity||self.calculator_identity!=self.state.relaxation.calculator_callback_identity||self.state.constraint_digest!=self.request.constraint_digest||self.state.identity_config!=self.request.identity_config||self.state.force_tolerance_ev_per_angstrom.to_bits()!=self.request.force_tolerance_ev_per_angstrom.to_bits()||self.state.max_movable_force_ev_per_angstrom>self.request.force_tolerance_ev_per_angstrom||self.state.relaxation.steps>self.request.step_limit||self.state.relaxation.calculator_evaluations>self.request.evaluation_limit||self.state.system.validate_fixed_against(&self.request.state).is_err()||relax_receipt_digest(&self.request,&self.state,&self.minimizer_identity,&self.calculator_identity)?!=self.receipt_digest{return Err(Status::simple(StatusCode::InvalidState,"minimizer","RELAX-002","relaxation receipt/request/output/callback/tolerance/budget binding mismatch"));}Ok(())}}
fn relax_receipt_digest(request:&RelaxRequest,state:&CommittedState,minimizer_identity:&str,calculator_identity:&str)->Result<String,Status>{let bytes=crate::checkpoint::canonical_json_bytes(&("spark-relax-receipt/1",request,state,minimizer_identity,calculator_identity))?;Ok(format!("sha256:{}",crate::identity::hex_sha256(&bytes)))}

pub fn relax_checked(minimizer:&dyn Minimizer,request:&RelaxRequest,calculator:&dyn Calculator,ledger:&mut ResourceLedger,cancel:&CancelToken)->Result<ValidatedRelaxation,Status>{
    validate_relax_request(request)?;if request.calculator_identity!=calculator.identity(){return Err(Status::simple(StatusCode::InvalidInput,"minimizer","RELAX-001","relaxation calculator identity mismatch"));}
    let deadline=ledger.callback_deadline()?;ledger.reserve_calculator(request.evaluation_limit,Some(&request.request_id),None)?;
    let attempt=minimizer.relax(request,calculator,deadline,cancel);ledger.finish_evaluation_reservation(request.evaluation_limit,attempt.calculator_evaluations,"relaxation",Some(&request.request_id),None)?;if attempt.termination_reason.is_empty(){return Err(Status::simple(StatusCode::RelaxNotConverged,"minimizer","RELAX-002","missing termination evidence"));}
    if cancel.is_cancelled(){return Err(Status::simple(StatusCode::Cancelled,"minimizer","CKPT-005","cancelled during relaxation"));}if Instant::now()>deadline{return Err(Status::simple(StatusCode::RelaxNotConverged,"minimizer","RELAX-003","relaxation deadline exceeded"));}
    let response=attempt.result?;if response.calculator_evaluations!=attempt.calculator_evaluations||response.termination_reason!=attempt.termination_reason{return Err(Status::simple(StatusCode::RelaxNotConverged,"minimizer","RELAX-002","attempt and response accounting mismatch"));}validate_relax_response(request,&response)?;let provenance=crate::model::RelaxationProvenance{minimizer_name:response.minimizer_name.clone(),minimizer_version:response.minimizer_version.clone(),minimizer_callback_identity:minimizer.identity().into(),calculator_callback_identity:calculator.identity().into(),calculator_evaluations:response.calculator_evaluations,steps:response.steps,termination_reason:response.termination_reason.clone()};let state=CommittedState::try_new(response.state,response.energy_ev,response.forces_ev_per_angstrom,request.force_tolerance_ev_per_angstrom,request.constraint_digest.clone(),provenance,request.identity_config.clone())?;if state.max_movable_force().to_bits()!=response.max_movable_force_ev_per_angstrom.to_bits(){return Err(Status::simple(StatusCode::RelaxNotConverged,"minimizer","RELAX-002","reported movable max force does not equal internal recomputation"));}let receipt_digest=relax_receipt_digest(request,&state,minimizer.identity(),calculator.identity())?;let validated=ValidatedRelaxation{request:request.clone(),state,minimizer_identity:minimizer.identity().into(),calculator_identity:calculator.identity().into(),receipt_digest};validated.validate_receipt()?;Ok(validated)
}
fn validate_relax_request(request:&RelaxRequest)->Result<(),Status>{request.state.validate()?;request.identity_config.validate()?;if request.calculator_identity.is_empty()||request.evaluation_limit==0||request.step_limit==0||request.force_tolerance_ev_per_angstrom<0.0||!request.force_tolerance_ev_per_angstrom.is_finite()||request.constraint_digest.is_empty()||request.request_id.is_empty(){return Err(Status::simple(StatusCode::InvalidInput,"minimizer","RELAX-001","invalid relaxation request"));}Ok(())}

pub fn validate_relax_response(req:&RelaxRequest,r:&RelaxResponse)->Result<(),Status>{
    r.state.validate_fixed_against(&req.state)?;
    if r.steps>req.step_limit||r.calculator_evaluations>req.evaluation_limit{return Err(Status::simple(StatusCode::RelaxNotConverged,"minimizer","RELAX-003","minimizer budget exceeded"));}
    if r.forces_ev_per_angstrom.len()!=r.state.atom_ids.len()||!r.energy_ev.is_finite()||!r.max_movable_force_ev_per_angstrom.is_finite()||!r.force_tolerance(req.force_tolerance_ev_per_angstrom)||!r.forces_ev_per_angstrom.iter().flatten().all(|x|x.is_finite()){
        return Err(Status::simple(StatusCode::RelaxNotConverged,"minimizer","RELAX-003","illegal successful relaxation"));}Ok(())
}
trait ForceTolerance{fn force_tolerance(&self,t:f64)->bool;}
impl ForceTolerance for RelaxResponse{fn force_tolerance(&self,t:f64)->bool{t.is_finite()&&t>=0.0&&self.max_movable_force_ev_per_angstrom<=t}}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SaddleSearchRequest{pub origin:CommittedState,pub search_family:String,pub searcher_identity:String,pub searcher_version:String,pub search_index:u64,pub active_region_hint:Option<Vec<usize>>,pub search_id:String,pub rng_substream:PhiloxState,pub rng_substream_digest:String,pub force_tolerance_ev_per_angstrom:f64,pub curvature_tolerance_ev_per_angstrom2:f64,pub evaluation_limit:u64,pub calculator_identity:String,pub minimizer_identity:String}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EndpointRequest{pub geometry:AtomicSystem,pub displacement_sign:i8}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SaddleCandidate{pub saddle:SaddleRecord,pub downhill_endpoints:[EndpointRequest;2]}
pub trait SaddleSearcher:Send+Sync{fn identity(&self)->&str;fn version(&self)->&str;fn search(&self,request:&SaddleSearchRequest,calculator:&dyn Calculator,minimizer:&dyn Minimizer,deadline:Instant,cancel:&CancelToken)->CallbackAttempt<SaddleCandidate>;}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ValidatedSaddle{pub(crate)request:SaddleSearchRequest,pub(crate)candidate:SaddleCandidate,pub(crate)receipt_digest:String}
impl ValidatedSaddle{pub fn request(&self)->&SaddleSearchRequest{&self.request}pub fn candidate(&self)->&SaddleCandidate{&self.candidate}pub(crate)fn validate_receipt(&self)->Result<(),Status>{validate_saddle_candidate(&self.request,&self.candidate)?;if saddle_receipt_digest(&self.request,&self.candidate)?!=self.receipt_digest{return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-007","saddle search receipt mismatch"));}Ok(())}}
fn saddle_receipt_digest(request:&SaddleSearchRequest,candidate:&SaddleCandidate)->Result<String,Status>{let bytes=crate::checkpoint::canonical_json_bytes(&("spark-saddle-receipt/1",request,candidate,request.searcher_identity.as_str(),request.searcher_version.as_str()))?;Ok(format!("sha256:{}",crate::identity::hex_sha256(&bytes)))}

pub fn search_checked(searcher:&dyn SaddleSearcher,request:&SaddleSearchRequest,calculator:&dyn Calculator,minimizer:&dyn Minimizer,ledger:&mut ResourceLedger,cancel:&CancelToken)->Result<ValidatedSaddle,Status>{
    validate_search_request(request)?;if request.searcher_identity!=searcher.identity()||request.searcher_version!=searcher.version()||request.calculator_identity!=calculator.identity()||request.minimizer_identity!=minimizer.identity(){return Err(Status::simple(StatusCode::InvalidInput,"saddle","SADDLE-001","searcher/version/callback identity mismatch"));}
    let deadline=ledger.callback_deadline()?;ledger.reserve_saddle_attempt(&request.origin.state_id)?;ledger.reserve_calculator(request.evaluation_limit,None,Some(&request.search_id))?;
    let attempt=searcher.search(request,calculator,minimizer,deadline,cancel);ledger.finish_evaluation_reservation(request.evaluation_limit,attempt.calculator_evaluations,"saddle-search",None,Some(&request.search_id))?;if attempt.termination_reason.is_empty(){return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-007","missing termination evidence"));}
    if cancel.is_cancelled(){return Err(Status::simple(StatusCode::Cancelled,"saddle","CKPT-005","cancelled during saddle search"));}if Instant::now()>deadline{return Err(Status::simple(StatusCode::SaddleNotFound,"saddle","SADDLE-007","saddle-search deadline exceeded"));}
    let candidate=attempt.result?;if candidate.saddle.evaluation_count!=attempt.calculator_evaluations||candidate.saddle.termination_reason!=attempt.termination_reason{return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-007","attempt and candidate accounting mismatch"));}validate_saddle_candidate(request,&candidate)?;let receipt_digest=saddle_receipt_digest(request,&candidate)?;let validated=ValidatedSaddle{request:request.clone(),candidate,receipt_digest};validated.validate_receipt()?;Ok(validated)
}

pub fn validate_saddle_candidate(req:&SaddleSearchRequest,c:&SaddleCandidate)->Result<(),Status>{
    validate_search_request(req)?;
    c.saddle.geometry.validate_fixed_against(&req.origin.system)?;
    let n=req.origin.system.atom_ids.len();
    if c.saddle.search_id!=req.search_id||c.saddle.forces_ev_per_angstrom.len()!=n||c.saddle.unstable_direction.len()!=n||c.saddle.evaluation_count>req.evaluation_limit
        ||!c.saddle.energy_ev.is_finite()||!c.saddle.curvature_ev_per_angstrom2.is_finite()
        ||!c.saddle.forces_ev_per_angstrom.iter().flatten().chain(c.saddle.unstable_direction.iter().flatten()).chain(c.saddle.orthogonal_curvatures_ev_per_angstrom2.iter()).all(|x|x.is_finite()){
        return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-003","malformed saddle response"));}
    if c.saddle.curvature_ev_per_angstrom2>=-req.curvature_tolerance_ev_per_angstrom2{return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-005","reported direction lacks negative curvature"));}
    let mut max_force=0.0_f64;let mut direction_norm2=0.0;for i in 0..n{let mut f2=0.0;for k in 0..3{if req.origin.system.movable[i][k]{f2+=c.saddle.forces_ev_per_angstrom[i][k].powi(2);direction_norm2+=c.saddle.unstable_direction[i][k].powi(2);}else if c.saddle.unstable_direction[i][k].abs()>1e-14{return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-005","unstable direction moves constrained coordinate"));}}max_force=max_force.max(f2.sqrt());}if max_force>req.force_tolerance_ev_per_angstrom||(direction_norm2-1.0).abs()>1e-8{return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-005","saddle force or normalized unstable-direction gate failed"));}
    match c.saddle.evidence_level{
        SaddleEvidence::Hessian if c.saddle.rigid_constrained_modes_excluded&&c.saddle.imaginary_mode_count_after_exclusions==Some(1)=>{},
        SaddleEvidence::Hessian=>return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-005","Hessian evidence does not contain exactly one post-exclusion imaginary mode")),
        SaddleEvidence::Directional if c.saddle.imaginary_mode_count_after_exclusions.is_some()||c.saddle.orthogonal_curvatures_ev_per_angstrom2.is_empty()=>return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-005","directional evidence requires at least one orthogonal sample and no Hessian count")),
        SaddleEvidence::Directional if c.saddle.orthogonal_curvatures_ev_per_angstrom2.iter().any(|x|*x < -req.curvature_tolerance_ev_per_angstrom2)=>return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-005","negative sampled orthogonal curvature")),
        SaddleEvidence::Directional=>{},
    }
    let signs=[c.downhill_endpoints[0].displacement_sign,c.downhill_endpoints[1].displacement_sign];if !((signs[0]==-1&&signs[1]==1)||(signs[0]==1&&signs[1]==-1)){return Err(Status::simple(StatusCode::InvalidSaddle,"saddle","SADDLE-003","endpoint signs must be the set {-1,+1}"));}for e in &c.downhill_endpoints{e.geometry.validate_fixed_against(&req.origin.system)?;}
    Ok(())
}

pub fn canonical_search_id(req:&SaddleSearchRequest)->Result<String,Status>{let mut active=req.active_region_hint.clone();if let Some(indices)=&mut active{indices.sort_unstable();}let bytes=crate::checkpoint::canonical_json_bytes(&(req.origin.state_id(),req.search_family.as_str(),req.searcher_identity.as_str(),req.searcher_version.as_str(),req.search_index,active,req.rng_substream_digest.as_str(),req.calculator_identity.as_str(),req.minimizer_identity.as_str(),req.force_tolerance_ev_per_angstrom.to_bits(),req.curvature_tolerance_ev_per_angstrom2.to_bits(),req.evaluation_limit))?;Ok(format!("search:sha256:{}",crate::identity::hex_sha256(&bytes)))}
fn validate_search_request(req:&SaddleSearchRequest)->Result<(),Status>{req.origin.validate()?;validate_state(&req.rng_substream)?;let rng=Philox::from_state(req.rng_substream.clone())?;if req.search_family.is_empty()||req.searcher_identity.is_empty()||req.searcher_version.is_empty()||req.evaluation_limit==0||!req.force_tolerance_ev_per_angstrom.is_finite()||req.force_tolerance_ev_per_angstrom<0.0||!req.curvature_tolerance_ev_per_angstrom2.is_finite()||req.curvature_tolerance_ev_per_angstrom2<0.0||substream_digest(&rng)!=req.rng_substream_digest{return Err(Status::simple(StatusCode::InvalidInput,"saddle","SADDLE-001","invalid searcher/RNG/tolerance resource binding"));}if let Some(active)=&req.active_region_hint{let mut x=active.clone();x.sort_unstable();if x.iter().any(|&i|i>=req.origin.system.atom_ids.len())||x.windows(2).any(|w|w[0]==w[1]){return Err(Status::simple(StatusCode::InvalidInput,"saddle","SADDLE-001","active-region indices invalid"));}}if canonical_search_id(req)?!=req.search_id{return Err(Status::simple(StatusCode::InvalidInput,"saddle","SADDLE-001","search ID is not canonical for request"));}Ok(())}

fn normalize_calc_error(mut s:Status)->Status{if s.status!=StatusCode::Cancelled&&s.status!=StatusCode::ResourceLimit&&s.status!=StatusCode::NonfiniteResult{s.status=StatusCode::CalculatorFailure;s.severity=StatusCode::CalculatorFailure.severity();}s}
fn calc_fail(m:&str)->Status{Status::simple(StatusCode::CalculatorFailure,"calculator","CALC-006",m)}
