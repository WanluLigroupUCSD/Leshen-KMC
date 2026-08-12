// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::identity::match_states;
use crate::model::{CommittedState,DirectedEvent};
use crate::rate::neumaier_sum;
use crate::catalog::{Catalog,RateCache};
use crate::rng::{validate_state,Philox,PhiloxState};
use crate::status::{Status,StatusCode};
use crate::parity::{RateTablePayload,RateTableSnapshot};
use serde::{Deserialize,Serialize};

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KmcStepRecord{pub checkpoint_sequence:u64,pub log_sequence:u64,pub post_state_id:String,pub pre_state_id:String,pub rate_table_snapshot:RateTableSnapshot,pub selected_event_id:String,pub selected_rate_per_s:f64,pub selection_uniform:f64,pub step_index:u64,pub time_increment_s:f64,pub time_uniform:f64,pub total_rate_per_s:f64}

#[derive(Clone,Debug)]
pub struct KmcTrajectory{current:CommittedState,simulation_time_s:f64,step_index:u64,checkpoint_sequence:u64,log_sequence:u64,rng:Philox,records:Vec<KmcStepRecord>}

impl KmcTrajectory{
    pub fn new(initial:CommittedState,rng:Philox,catalog:&Catalog,cache:&RateCache)->Result<Self,Status>{initial.validate()?;cache.validate_for(catalog)?;validate_state(&rng.state())?;if catalog.states().get(initial.state_id())!=Some(&initial){return Err(Status::simple(StatusCode::CatalogIncompatible,"kmc","KMC-001","initial state is not from validated catalog epoch"));}Ok(Self{current:initial,simulation_time_s:0.0,step_index:0,checkpoint_sequence:0,log_sequence:0,rng,records:Vec::new()})}
    pub(crate) fn from_replayed(current:CommittedState,simulation_time_s:f64,step_index:u64,checkpoint_sequence:u64,log_sequence:u64,rng:Philox,records:Vec<KmcStepRecord>,catalog:&Catalog,cache:&RateCache)->Result<Self,Status>{current.validate()?;cache.validate_for(catalog)?;validate_state(&rng.state())?;if catalog.states().get(current.state_id())!=Some(&current)||!simulation_time_s.is_finite()||simulation_time_s<0.0||step_index!=records.len()as u64||log_sequence!=step_index||rng.state().consumed_uniforms!=step_index.checked_mul(2).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"kmc","KMC-005","resume uniform counter overflow"))?||records.last().map_or(step_index!=0,|r|r.post_state_id!=current.state_id){return Err(Status::simple(StatusCode::CheckpointCorrupt,"kmc","CKPT-004","replayed trajectory parts are inconsistent"));}Ok(Self{current,simulation_time_s,step_index,checkpoint_sequence,log_sequence,rng,records})}
    pub fn current(&self)->&CommittedState{&self.current}pub fn simulation_time_s(&self)->f64{self.simulation_time_s}pub fn step_index(&self)->u64{self.step_index}pub fn checkpoint_sequence(&self)->u64{self.checkpoint_sequence}pub fn log_sequence(&self)->u64{self.log_sequence}pub fn records(&self)->&[KmcStepRecord]{&self.records}
    pub fn step<F>(&mut self,catalog:&Catalog,cache:&RateCache,mut apply_and_reconverge:F)->Result<&KmcStepRecord,Status>
    where F:FnMut(&DirectedEvent,&CommittedState)->Result<CommittedState,Status>{
        self.current.validate()?;cache.validate_for(catalog)?;if catalog.states().get(self.current.state_id())!=Some(&self.current){return Err(Status::simple(StatusCode::CatalogIncompatible,"kmc","KMC-001","trajectory state is not from validated catalog epoch"));}let mut enabled:Vec<&DirectedEvent>=catalog.events().iter().filter(|e|e.origin_state_id==self.current.state_id).collect();enabled.sort_by(|a,b|a.event_id.cmp(&b.event_id));
        if enabled.is_empty(){return Err(Status::simple(StatusCode::NoEnabledEvent,"kmc","KMC-002","no enabled event"));}
        let mut rates=Vec::with_capacity(enabled.len());for e in &enabled{rates.push(*cache.numeric_rates.get(&e.event_id).ok_or_else(||Status::simple(StatusCode::CatalogIncompatible,"kmc","KMC-001","validated rate cache lacks event"))?);}let total=neumaier_sum(&rates)?;if total==0.0{return Err(Status::simple(StatusCode::NoEnabledEvent,"kmc","KMC-002","zero total rate"));}if !total.is_finite(){return Err(Status::simple(StatusCode::RateInvalid,"kmc","KMC-002","nonfinite total rate"));}let rate_table_snapshot=RateTableSnapshot::new(RateTablePayload{destination_state_ids:enabled.iter().map(|e|e.destination_state_id.clone()).collect(),event_ids:enabled.iter().map(|e|e.event_id.clone()).collect(),log_rates:enabled.iter().map(|e|e.rate_model.log_forward_rate_per_s).collect(),lost_rate_log_upper_bound:None,origin_state_id:self.current.state_id.clone(),rates:rates.clone(),schema:"spark-atomistic-rate-table-snapshot/1".into(),total_rate_per_s:total}).map_err(|_|Status::simple(StatusCode::RateInvalid,"kmc","E2-KMC-002","historical rate snapshot invalid"))?;
        let(u_select,u_time,next_rng)=self.rng.two_uniforms_atomic()?;let threshold=u_select*total;let mut cumulative=0.0;let mut chosen=enabled.len()-1;
        for(i,r)in rates.iter().enumerate(){cumulative+=*r;if cumulative>threshold{chosen=i;break;}}
        let event=enabled[chosen];let expected=catalog.states().get(&event.destination_state_id).ok_or_else(||Status::simple(StatusCode::CatalogIncompatible,"kmc","EVENT-004","validated destination record missing"))?;
        let recovered=apply_and_reconverge(event,&self.current).map_err(|s|application_failure(s,event))?;
        let tol=expected.identity_config().match_tolerances();if match_states(expected,&recovered,&tol).is_err(){return Err(Status::simple(StatusCode::EventApplicationFailed,"kmc","EVENT-004","reconverged endpoint does not match validated destination"));}
        let dt=-u_time.ln()/total;if !dt.is_finite()||dt<=0.0{return Err(Status::simple(StatusCode::RateInvalid,"kmc","KMC-002","invalid residence time"));}
        let new_time=self.simulation_time_s+dt;if !new_time.is_finite(){return Err(Status::simple(StatusCode::RateInvalid,"kmc","KMC-002","simulation time overflow"));}
        let new_step=self.step_index.checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"kmc","KMC-005","step index overflow"))?;
        let new_log=self.log_sequence.checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"kmc","KMC-004","log sequence overflow"))?;
        let record=KmcStepRecord{checkpoint_sequence:self.checkpoint_sequence,log_sequence:new_log,post_state_id:expected.state_id.clone(),pre_state_id:self.current.state_id.clone(),rate_table_snapshot,selected_event_id:event.event_id.clone(),selected_rate_per_s:rates[chosen],selection_uniform:u_select,step_index:new_step,time_increment_s:dt,time_uniform:u_time,total_rate_per_s:total};
        self.current=expected.clone();self.simulation_time_s=new_time;self.step_index=new_step;self.log_sequence=new_log;self.rng=next_rng;self.records.push(record);Ok(self.records.last().expect("just pushed"))
    }
    pub fn rng_state(&self)->PhiloxState{self.rng.state()}
}
fn application_failure(mut s:Status,e:&DirectedEvent)->Status{if s.status!=StatusCode::Cancelled&&s.status!=StatusCode::ResourceLimit{s.status=StatusCode::EventApplicationFailed;s.severity=StatusCode::EventApplicationFailed.severity();s.context.requirement_id="EVENT-004".into();s.context.search_or_event_id=Some(e.event_id.clone());}s}
