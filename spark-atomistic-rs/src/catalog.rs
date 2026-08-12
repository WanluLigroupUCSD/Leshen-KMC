// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::callbacks::{validate_saddle_candidate,ValidatedRelaxation,ValidatedSaddle};
use crate::checkpoint::canonical_json_bytes;
use crate::identity::{canonical_geometry_digest,hex_sha256,match_states,IdentityConfig,MatchTolerances};
use crate::model::{CommittedState,DirectedEvent,GeometryMatch,RelaxationProvenance,ValidationEvidence};
use crate::rate::common_prefactor_pair;
use crate::resource::ResourceLedger;
use crate::status::{Status,StatusCode};
use serde::{Deserialize,Serialize};
use std::collections::BTreeMap;

#[derive(Clone,Debug)]
pub struct EventTolerances{pub state:MatchTolerances,pub saddle_rms_angstrom:f64,pub saddle_max_angstrom:f64,pub saddle_energy_ev:f64,pub barrier_ev:f64,pub direction_cosine_abs_min:f64,pub epsilon_db:f64}
impl EventTolerances{pub fn from_identity(c:&IdentityConfig)->Result<Self,Status>{c.validate()?;Ok(Self{state:c.match_tolerances(),saddle_rms_angstrom:c.saddle_rms_angstrom,saddle_max_angstrom:c.saddle_max_angstrom,saddle_energy_ev:c.saddle_energy_ev,barrier_ev:c.barrier_tolerance_ev,direction_cosine_abs_min:c.direction_cosine_abs_min,epsilon_db:c.detailed_balance_epsilon})}fn validate_bound(&self,c:&IdentityConfig)->Result<(),Status>{let expected=Self::from_identity(c)?;if self.state!=expected.state||[self.saddle_rms_angstrom,self.saddle_max_angstrom,self.saddle_energy_ev,self.barrier_ev,self.direction_cosine_abs_min,self.epsilon_db].map(f64::to_bits)!=[expected.saddle_rms_angstrom,expected.saddle_max_angstrom,expected.saddle_energy_ev,expected.barrier_ev,expected.direction_cosine_abs_min,expected.epsilon_db].map(f64::to_bits){return Err(Status::simple(StatusCode::CatalogIncompatible,"catalog","CAT-006","event tolerances are not bound to canonical identity configuration"));}Ok(())}}

#[derive(Clone,Debug)]
struct ValidationReceipt{digest:String}

#[derive(Clone,Debug)]
struct EndpointReceipt{digest:String}
#[derive(Clone,Debug)]
pub struct ValidatedEndpoint{state:CommittedState,relaxation:ValidatedRelaxation,displacement_sign:i8,receipt:EndpointReceipt}

#[derive(Clone,Debug,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UnvalidatedDiagnostic{pub search_id:String,pub status:StatusCode,pub reason:String}
#[derive(Clone,Debug)]
pub struct CandidateStaging{allow_unvalidated:bool,diagnostics:Vec<UnvalidatedDiagnostic>}
impl CandidateStaging{pub fn new(allow_unvalidated:bool)->Self{Self{allow_unvalidated,diagnostics:Vec::new()}}pub fn retain(&mut self,diagnostic:UnvalidatedDiagnostic)->Result<(),Status>{if !self.allow_unvalidated{return Err(Status::simple(StatusCode::InvalidSaddle,"catalog","CAT-005","allow_unvalidated=false forbids diagnostic retention"));}if diagnostic.search_id.is_empty()||diagnostic.reason.is_empty()||matches!(diagnostic.status,StatusCode::Ok|StatusCode::DuplicateEvent){return Err(Status::simple(StatusCode::InvalidInput,"catalog","CAT-005","invalid unvalidated diagnostic"));}self.diagnostics.push(diagnostic);Ok(())}pub fn diagnostics(&self)->&[UnvalidatedDiagnostic]{&self.diagnostics}}

pub fn bind_validated_endpoint(saddle:&ValidatedSaddle,endpoint:&crate::callbacks::EndpointRequest,relaxation:ValidatedRelaxation)->Result<ValidatedEndpoint,Status>{saddle.validate_receipt()?;relaxation.validate_receipt()?;let req=saddle.request();let state=relaxation.state().clone();if relaxation.request.state!=endpoint.geometry||relaxation.request.calculator_identity!=req.calculator_identity||relaxation.minimizer_identity!=req.minimizer_identity||state.identity_config!=req.origin.identity_config||!(endpoint.displacement_sign==1||endpoint.displacement_sign == -1){return Err(Status::simple(StatusCode::InvalidState,"catalog","SADDLE-004","endpoint relaxation receipt is not bound to endpoint/search/callback request"));}let digest=endpoint_receipt_digest(saddle,endpoint,&relaxation)?;Ok(ValidatedEndpoint{state,relaxation,displacement_sign:endpoint.displacement_sign,receipt:EndpointReceipt{digest}})}
fn endpoint_receipt_digest(saddle:&ValidatedSaddle,endpoint:&crate::callbacks::EndpointRequest,relaxation:&ValidatedRelaxation)->Result<String,Status>{let bytes=canonical_json_bytes(&("spark-endpoint-receipt/1",saddle,endpoint,relaxation))?;Ok(format!("sha256:{}",hex_sha256(&bytes)))}

#[derive(Clone,Debug)]
pub struct ValidatedPair{forward:DirectedEvent,reverse:DirectedEvent,origin:CommittedState,destination:CommittedState,receipt:ValidationReceipt}
impl ValidatedPair{pub fn events(&self)->(&DirectedEvent,&DirectedEvent){(&self.forward,&self.reverse)}}

#[derive(Clone,Debug)]
pub struct CatalogCommitOutcome{pub(crate)status:StatusCode,pub(crate)event_id:String,pub(crate)log_rate:f64,pub(crate)validation_digest:String}
impl CatalogCommitOutcome{pub fn status(&self)->StatusCode{self.status}pub fn event_id(&self)->&str{&self.event_id}}

pub fn validate_reversible_candidate(saddle:ValidatedSaddle,endpoints:[ValidatedEndpoint;2],
    environment_key:String,environment_version:String,active_atom_mapping:Vec<[usize;2]>,
    discovery_provenance:BTreeMap<String,serde_json::Value>,temperature_k:f64,prefactor_per_s:f64,tol:&EventTolerances)->Result<ValidatedPair,Status>{
    saddle.validate_receipt()?;let req=saddle.request();let candidate=saddle.candidate();validate_saddle_candidate(req,candidate)?;for e in &endpoints{e.relaxation.validate_receipt()?;let request=candidate.downhill_endpoints.iter().find(|x|x.displacement_sign==e.displacement_sign).ok_or_else(||Status::simple(StatusCode::InvalidState,"catalog","SADDLE-004","endpoint sign absent from candidate"))?;if endpoint_receipt_digest(&saddle,request,&e.relaxation)?!=e.receipt.digest{return Err(Status::simple(StatusCode::InvalidState,"catalog","SADDLE-004","internal endpoint validation receipt mismatch"));}}
    tol.validate_bound(&req.origin.identity_config)?;
    for e in &endpoints{e.state.system.validate_fixed_against(&req.origin.system)?;}
    let m0=try_match(&req.origin,&endpoints[0].state,&tol.state)?;let m1=try_match(&req.origin,&endpoints[1].state,&tol.state)?;
    let (origin_match,destination,destination_match)=match(m0,m1){
        (None,None)=>return Err(Status::simple(StatusCode::SaddleWrongBasin,"catalog","SADDLE-004","neither endpoint matches origin")),
        (Some(_),Some(_))=>return Err(Status::simple(StatusCode::EndpointCollapsed,"catalog","SADDLE-004","both endpoints match origin")),
        (Some(m),None)=>(m,&endpoints[1].state,match_states(&endpoints[1].state,&endpoints[1].state,&tol.state)?),
        (None,Some(m))=>(m,&endpoints[0].state,match_states(&endpoints[0].state,&endpoints[0].state,&tol.state)?),
    };
    if try_match(&req.origin,destination,&tol.state)?.is_some(){return Err(Status::simple(StatusCode::EndpointCollapsed,"catalog","SADDLE-004","destination equals origin"));}
    let rate=common_prefactor_pair(req.origin.energy_ev,destination.energy_ev,candidate.saddle.energy_ev,temperature_k,prefactor_per_s,tol.barrier_ev,tol.epsilon_db)?;
    let bf=candidate.saddle.energy_ev-req.origin.energy_ev;let br=candidate.saddle.energy_ev-destination.energy_ev;
    let mut endpoint_ids=[req.origin.state_id.clone(),destination.state_id.clone()];endpoint_ids.sort();let saddle_geometry_digest=canonical_geometry_digest(&candidate.saddle.geometry,&req.origin.identity_config)?;let direction_bits:Vec<u64>=candidate.saddle.unstable_direction.iter().flatten().map(|x|if *x==0.0{0}else{x.to_bits()}).collect();let neg_direction_bits:Vec<u64>=candidate.saddle.unstable_direction.iter().flatten().map(|x|if *x==0.0{0}else{(-*x).to_bits()}).collect();let canonical_direction=if direction_bits<=neg_direction_bits{direction_bits}else{neg_direction_bits};let pair_payload=canonical_json_bytes(&("spark-event-pair/1",endpoint_ids,saddle_geometry_digest,candidate.saddle.energy_ev.to_bits(),canonical_direction))?;
    let pair_id=format!("pair:sha256:{}",hex_sha256(&pair_payload));
    let f_id=format!("event:sha256:{}",hex_sha256(&canonical_json_bytes(&(pair_id.as_str(),req.origin.state_id.as_str(),destination.state_id.as_str()))?));
    let r_id=format!("event:sha256:{}",hex_sha256(&canonical_json_bytes(&(pair_id.as_str(),destination.state_id.as_str(),req.origin.state_id.as_str()))?));
    let validation=ValidationEvidence{method:candidate.saddle.evidence_level,origin_match,destination_match,unstable_mode_count:1,
        full_endpoint_relaxations:true,calculator_model_digest:req.origin.system.calculator_model_digest.clone(),constraint_digest:req.origin.constraint_digest.clone(),search_id:req.search_id.clone(),rng_substream_digest:req.rng_substream_digest.clone(),calculator_callback_identity:req.calculator_identity.clone(),minimizer_callback_identity:req.minimizer_identity.clone(),endpoint_states:[endpoints[0].state.clone(),endpoints[1].state.clone()],endpoint_receipt_digests:[endpoints[0].receipt.digest.clone(),endpoints[1].receipt.digest.clone()],saddle_search:saddle.clone(),endpoint_relaxations:[endpoints[0].relaxation.clone(),endpoints[1].relaxation.clone()]};
    let mut reverse_validation=validation.clone();std::mem::swap(&mut reverse_validation.origin_match,&mut reverse_validation.destination_match);
    let identity_digest=req.origin.identity_config.digest()?;let tolerance_digest=identity_digest.clone();let schema_digest=format!("sha256:{}",hex_sha256(crate::IR_SCHEMA.as_bytes()));
    let common=(candidate.saddle.clone(),active_atom_mapping,environment_key,environment_version,discovery_provenance,validation,
        req.origin.system.calculator_model_digest.clone(),identity_digest,tolerance_digest,schema_digest);
    let forward=DirectedEvent{event_id:f_id.clone(),reverse_pair_id:r_id.clone(),origin_state_id:req.origin.state_id.clone(),destination_state_id:destination.state_id.clone(),
        saddle:common.0.clone(),barrier_ev:bf,reverse_barrier_ev:br,rate_model:rate.clone(),active_atom_mapping:common.1.clone(),environment_key:common.2.clone(),environment_version:common.3.clone(),discovery_provenance:common.4.clone(),validation:common.5.clone(),calculator_digest:common.6.clone(),identity_digest:common.7.clone(),tolerance_digest:common.8.clone(),schema_digest:common.9.clone()};
    let mut reverse_rate=rate;std::mem::swap(&mut reverse_rate.log_forward_rate_per_s,&mut reverse_rate.log_reverse_rate_per_s);reverse_rate.detailed_balance_residual=-reverse_rate.detailed_balance_residual;
    let reverse=DirectedEvent{event_id:r_id,reverse_pair_id:f_id,origin_state_id:destination.state_id.clone(),destination_state_id:req.origin.state_id.clone(),
        saddle:common.0,barrier_ev:br,reverse_barrier_ev:bf,rate_model:reverse_rate,active_atom_mapping:common.1.into_iter().map(|x|[x[1],x[0]]).collect(),environment_key:common.2,environment_version:common.3,discovery_provenance:common.4,validation:reverse_validation,calculator_digest:common.6,identity_digest:common.7,tolerance_digest:common.8,schema_digest:common.9};
    let digest=pair_receipt_digest(&forward,&reverse)?;Ok(ValidatedPair{forward,reverse,origin:req.origin.clone(),destination:destination.clone(),receipt:ValidationReceipt{digest}})
}

fn try_match(a:&CommittedState,b:&CommittedState,t:&MatchTolerances)->Result<Option<GeometryMatch>,Status>{match match_states(a,b,t){Ok(m)=>Ok(Some(m)),Err(s) if s.status==StatusCode::InvalidState=>Ok(None),Err(s)=>Err(s)}}
fn pair_receipt_digest(a:&DirectedEvent,b:&DirectedEvent)->Result<String,Status>{let bytes=canonical_json_bytes(&(a,b))?;Ok(hex_sha256(&bytes))}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogSnapshot{pub epoch:u64,pub digest:String,pub states:BTreeMap<String,CommittedState>,pub events:Vec<DirectedEvent>,pub duplicate_multiplicity:BTreeMap<String,u64>}

#[derive(Clone,Debug)]
pub struct Catalog{epoch:u64,digest:String,states:BTreeMap<String,CommittedState>,events:Vec<DirectedEvent>,duplicate_multiplicity:BTreeMap<String,u64>}
impl Catalog{
    pub fn empty()->Self{let mut c=Self{epoch:0,digest:String::new(),states:BTreeMap::new(),events:Vec::new(),duplicate_multiplicity:BTreeMap::new()};c.digest=c.compute_digest().expect("empty catalog serialization is finite");c}
    pub fn snapshot(&self)->CatalogSnapshot{CatalogSnapshot{epoch:self.epoch,digest:self.digest.clone(),states:self.states.clone(),events:self.events.clone(),duplicate_multiplicity:self.duplicate_multiplicity.clone()}}
    pub fn from_snapshot(s:CatalogSnapshot)->Result<Self,Status>{
        for(id,state)in &s.states{state.validate()?;if id!=&state.state_id{return Err(incompatible("catalog state key/ID mismatch"));}}
        if let Some(reference)=s.states.values().next(){for state in s.states.values(){state.system.validate_fixed_against(&reference.system)?;}}
        let mut ids=BTreeMap::new();for e in &s.events{validate_event_record(e,&s.states)?;if ids.insert(&e.event_id,true).is_some(){return Err(incompatible("duplicate event ID in catalog snapshot"));}}
        for e in &s.events{if !s.states.contains_key(&e.origin_state_id)||!s.states.contains_key(&e.destination_state_id){return Err(incompatible("event references absent committed state"));}let reciprocal=s.events.iter().find(|x|x.event_id==e.reverse_pair_id).ok_or_else(||incompatible("missing reciprocal event"))?;if reciprocal.reverse_pair_id!=e.event_id||reciprocal.origin_state_id!=e.destination_state_id||reciprocal.destination_state_id!=e.origin_state_id||reciprocal.saddle!=e.saddle||reciprocal.barrier_ev.to_bits()!=e.reverse_barrier_ev.to_bits()||reciprocal.reverse_barrier_ev.to_bits()!=e.barrier_ev.to_bits()||reciprocal.rate_model.log_forward_rate_per_s.to_bits()!=e.rate_model.log_reverse_rate_per_s.to_bits()||reciprocal.rate_model.log_reverse_rate_per_s.to_bits()!=e.rate_model.log_forward_rate_per_s.to_bits(){return Err(incompatible("invalid reciprocal scientific relationship"));}}
        let c=Self{epoch:s.epoch,digest:s.digest,states:s.states,events:s.events,duplicate_multiplicity:s.duplicate_multiplicity};if c.compute_digest()?!=c.digest{return Err(incompatible("catalog digest mismatch"));}Ok(c)
    }
    pub fn events(&self)->&[DirectedEvent]{&self.events}pub fn states(&self)->&BTreeMap<String,CommittedState>{&self.states}pub fn epoch(&self)->u64{self.epoch}pub fn digest(&self)->&str{&self.digest}
    pub fn rate_cache(&self,temperature_k:f64)->Result<RateCache,Status>{if !temperature_k.is_finite()||temperature_k<=0.0{return Err(Status::simple(StatusCode::RateInvalid,"catalog","RATE-008","invalid cache temperature"));}let mut rates=BTreeMap::new();for e in &self.events{let origin=self.states.get(&e.origin_state_id).ok_or_else(||incompatible("rate-cache origin missing"))?;let destination=self.states.get(&e.destination_state_id).ok_or_else(||incompatible("rate-cache destination missing"))?;let pref=e.rate_model.common_prefactor_per_s.ok_or_else(||incompatible("common prefactor missing"))?;let config=&origin.identity_config;let model=common_prefactor_pair(origin.energy_ev,destination.energy_ev,e.saddle.energy_ev,temperature_k,pref,config.barrier_tolerance_ev,config.detailed_balance_epsilon)?;rates.insert(e.event_id.clone(),crate::rate::exp_selectable(model.log_forward_rate_per_s)?);}Ok(RateCache{catalog_epoch:self.epoch,catalog_digest:self.digest.clone(),temperature_k,numeric_rates:rates})}
    pub fn commit_pair(&mut self,pair:ValidatedPair,tol:&EventTolerances,ledger:&mut ResourceLedger)->Result<CatalogCommitOutcome,Status>{
        if pair_receipt_digest(&pair.forward,&pair.reverse)?!=pair.receipt.digest{return Err(Status::simple(StatusCode::InternalError,"catalog","CAT-007","validation receipt mismatch"));}
        let pair_states=BTreeMap::from([(pair.origin.state_id.clone(),pair.origin.clone()),(pair.destination.state_id.clone(),pair.destination.clone())]);validate_event_record(&pair.forward,&pair_states)?;validate_event_record(&pair.reverse,&pair_states)?;
        let mut candidate_endpoints=[pair.forward.origin_state_id.as_str(),pair.forward.destination_state_id.as_str()];candidate_endpoints.sort();
        for existing in self.events.iter().filter(|x|x.event_id.as_str()<x.reverse_pair_id.as_str()){let mut existing_endpoints=[existing.origin_state_id.as_str(),existing.destination_state_id.as_str()];existing_endpoints.sort();if existing_endpoints!=candidate_endpoints{continue}let candidate_direction=if existing.origin_state_id==pair.forward.origin_state_id{&pair.forward}else{&pair.reverse};
            if let Some(geometry_match)=saddle_geometry_match(existing,candidate_direction,tol)?{
                if (existing.saddle.energy_ev-pair.forward.saddle.energy_ev).abs()>tol.saddle_energy_ev{return Err(Status::simple(StatusCode::CatalogConflict,"catalog","CAT-004","same saddle geometry has inconsistent energy"));}
                if direction_same(existing,candidate_direction,&geometry_match.mapping_a_to_b,tol.direction_cosine_abs_min)&&existing.active_atom_mapping==candidate_direction.active_atom_mapping{
                    let v=self.duplicate_multiplicity.get(&existing.event_id).copied().unwrap_or(1).checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"catalog","CAT-003","duplicate count overflow"))?;
                    let mut next=self.clone();next.duplicate_multiplicity.insert(existing.event_id.clone(),v);next.advance_epoch()?;let out=CatalogCommitOutcome{status:StatusCode::DuplicateEvent,event_id:existing.event_id.clone(),log_rate:existing.rate_model.log_forward_rate_per_s,validation_digest:pair.receipt.digest};*self=next;return Ok(out);
                }
            }
        }
        for state in [&pair.origin,&pair.destination]{if let Some(existing)=self.states.get(&state.state_id){if existing!=state{return Err(Status::simple(StatusCode::CatalogConflict,"catalog","CAT-004","state ID maps to inconsistent committed record"));}}}
        let event_id=pair.forward.event_id.clone();let log_rate=pair.forward.rate_model.log_forward_rate_per_s;let validation_digest=pair.receipt.digest.clone();let mut next=self.clone();next.states.insert(pair.origin.state_id.clone(),pair.origin);next.states.insert(pair.destination.state_id.clone(),pair.destination);next.events.push(pair.forward);next.events.push(pair.reverse);next.events.sort_by(|a,b|a.event_id.cmp(&b.event_id));next.advance_epoch()?;ledger.reserve_catalog_event(2)?;*self=next;Ok(CatalogCommitOutcome{status:StatusCode::Ok,event_id,log_rate,validation_digest})
    }
    fn advance_epoch(&mut self)->Result<(),Status>{self.epoch=self.epoch.checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"catalog","CAT-007","catalog epoch overflow"))?;self.digest=self.compute_digest()?;Ok(())}
    fn compute_digest(&self)->Result<String,Status>{let bytes=canonical_json_bytes(&(self.epoch,&self.states,&self.events,&self.duplicate_multiplicity))?;Ok(format!("sha256:{}",hex_sha256(&bytes)))}
}

fn validate_event_record(e:&DirectedEvent,states:&BTreeMap<String,CommittedState>)->Result<(),Status>{
    if e.event_id.is_empty()||e.reverse_pair_id.is_empty()||e.origin_state_id.is_empty()||e.destination_state_id.is_empty()||e.origin_state_id==e.destination_state_id{return Err(incompatible("invalid event identity/endpoints"));}
    let nums=[e.saddle.energy_ev,e.barrier_ev,e.reverse_barrier_ev,e.rate_model.temperature_k,e.rate_model.log_forward_rate_per_s,e.rate_model.log_reverse_rate_per_s,e.rate_model.detailed_balance_residual];if !nums.iter().all(|x|x.is_finite()){return Err(incompatible("invalid event numeric data"));}
    e.validation.saddle_search.validate_receipt()?;let persisted_request=e.validation.saddle_search.request();let persisted_candidate=e.validation.saddle_search.candidate();validate_saddle_candidate(persisted_request,persisted_candidate)?;if persisted_candidate.saddle!=e.saddle||e.validation.search_id!=persisted_request.search_id||e.validation.rng_substream_digest!=persisted_request.rng_substream_digest||e.validation.calculator_callback_identity!=persisted_request.calculator_identity||e.validation.minimizer_callback_identity!=persisted_request.minimizer_identity{return Err(incompatible("persisted full saddle request/search receipt does not bind event"));}
    if !e.validation.full_endpoint_relaxations||e.validation.unstable_mode_count!=1||e.validation.calculator_model_digest!=e.calculator_digest||e.validation.search_id!=e.saddle.search_id||e.validation.endpoint_receipt_digests.iter().any(|x|x.is_empty()){return Err(incompatible("event lacks full internal validation evidence"));}
    let origin=states.get(&e.origin_state_id).ok_or_else(||incompatible("event origin state missing"))?;let destination=states.get(&e.destination_state_id).ok_or_else(||incompatible("event destination state missing"))?;let config=&origin.identity_config;if destination.identity_config!=*config{return Err(incompatible("event endpoints use different identity configuration"));}let tol=EventTolerances::from_identity(config)?;
    e.saddle.geometry.validate_fixed_against(&origin.system)?;if e.saddle.energy_ev-origin.energy_ev!=e.barrier_ev||e.saddle.energy_ev-destination.energy_ev!=e.reverse_barrier_ev{return Err(incompatible("barriers are not exact differences from one saddle"));}let pref=e.rate_model.common_prefactor_per_s.ok_or_else(||incompatible("P0 event lacks common prefactor"))?;let recomputed=common_prefactor_pair(origin.energy_ev,destination.energy_ev,e.saddle.energy_ev,e.rate_model.temperature_k,pref,tol.barrier_ev,tol.epsilon_db)?;if recomputed!=e.rate_model{return Err(incompatible("rate/detailed-balance record does not recompute"));}
    if e.identity_digest!=config.digest()?||e.tolerance_digest!=e.identity_digest||e.schema_digest!=format!("sha256:{}",hex_sha256(crate::IR_SCHEMA.as_bytes()))||e.calculator_digest!=origin.system.calculator_model_digest{return Err(incompatible("event digest binding mismatch"));}
    if persisted_request.origin.state_id!=e.origin_state_id&&persisted_request.origin.state_id!=e.destination_state_id{return Err(incompatible("saddle request origin is outside event pair"));}let n=origin.system.atom_ids.len();let mut seen_left=BTreeMap::new();let mut seen_right=BTreeMap::new();for m in &e.active_atom_mapping{if m[0]>=n||m[1]>=n||seen_left.insert(m[0],m[1]).is_some()||seen_right.insert(m[1],m[0]).is_some()||origin.system.species[m[0]]!=destination.system.species[m[1]]{return Err(incompatible("active-atom mapping is not an in-range species-preserving partial bijection"));}}
    for i in 0..2{let relaxation=&e.validation.endpoint_relaxations[i];relaxation.validate_receipt()?;if relaxation.state()!=&e.validation.endpoint_states[i]{return Err(incompatible("endpoint state differs from persisted relaxation receipt"));}let mut receipt_matches=0;for endpoint in &persisted_candidate.downhill_endpoints{if endpoint_receipt_digest(&e.validation.saddle_search,endpoint,relaxation)?==e.validation.endpoint_receipt_digests[i]{receipt_matches+=1;}}if receipt_matches!=1{return Err(incompatible("endpoint search/relaxation receipt does not uniquely recompute"));}}
    for endpoint in &e.validation.endpoint_states{endpoint.validate()?;endpoint.system.validate_fixed_against(&origin.system)?;if endpoint.identity_config!=*config{return Err(incompatible("endpoint identity config mismatch"));}}
    let origin_matches:Vec<_>=e.validation.endpoint_states.iter().filter_map(|x|match_states(origin,x,&tol.state).ok()).collect();let destination_matches:Vec<_>=e.validation.endpoint_states.iter().filter_map(|x|match_states(destination,x,&tol.state).ok()).collect();if origin_matches.len()!=1||destination_matches.len()!=1||origin_matches[0]!=e.validation.origin_match||destination_matches[0]!=e.validation.destination_match{return Err(incompatible("endpoint geometry/mapping evidence does not recompute"));}Ok(())
}
fn saddle_geometry_match(a:&DirectedEvent,b:&DirectedEvent,t:&EventTolerances)->Result<Option<GeometryMatch>,Status>{
    let mk=|e:&DirectedEvent|CommittedState::try_new(e.saddle.geometry.clone(),0.0,e.saddle.forces_ev_per_angstrom.clone(),f64::MAX,e.validation.constraint_digest.clone(),RelaxationProvenance{minimizer_name:"validation".into(),minimizer_version:"1".into(),minimizer_callback_identity:e.validation.minimizer_callback_identity.clone(),calculator_callback_identity:e.validation.calculator_callback_identity.clone(),calculator_evaluations:0,steps:0,termination_reason:"catalog geometry comparison".into()},e.validation.endpoint_states[0].identity_config.clone());
    let mut mt=t.state.clone();mt.rms_angstrom=t.saddle_rms_angstrom;mt.max_angstrom=t.saddle_max_angstrom;
    match match_states(&mk(a)?,&mk(b)?,&mt){Ok(m)=>Ok(Some(m)),Err(s)if s.status==StatusCode::InvalidState=>Ok(None),Err(s)=>Err(s)}
}
fn direction_same(a:&DirectedEvent,b:&DirectedEvent,mapping:&[usize],min_abs_cos:f64)->bool{if mapping.len()!=a.saddle.unstable_direction.len(){return false}let mut dot=0.0;let mut aa=0.0;let mut bb=0.0;for(i,&j)in mapping.iter().enumerate(){if j>=b.saddle.unstable_direction.len(){return false}for k in 0..3{let x=a.saddle.unstable_direction[i][k];let y=b.saddle.unstable_direction[j][k];dot+=x*y;aa+=x*x;bb+=y*y;}}aa>0.0&&bb>0.0&&(dot/(aa.sqrt()*bb.sqrt())).abs()>=min_abs_cos}
fn incompatible(m:&str)->Status{Status::simple(StatusCode::CatalogIncompatible,"catalog","CAT-006",m)}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RateCache{pub(crate)catalog_epoch:u64,pub(crate)catalog_digest:String,pub(crate)temperature_k:f64,pub(crate)numeric_rates:BTreeMap<String,f64>}
impl RateCache{pub fn validate_for(&self,c:&Catalog)->Result<(),Status>{if self.catalog_epoch!=c.epoch||self.catalog_digest!=c.digest||!self.temperature_k.is_finite()||self.temperature_k<=0.0{return Err(incompatible("immutable rate-cache epoch/digest/temperature mismatch"));}let expected=c.rate_cache(self.temperature_k)?;if self.numeric_rates!=expected.numeric_rates{return Err(Status::simple(StatusCode::RateInvalid,"catalog","RATE-004","cached rates do not recompute from validated catalog"));}Ok(())}pub fn temperature_k(&self)->f64{self.temperature_k}}
