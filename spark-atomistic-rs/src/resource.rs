// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::status::{Status,StatusCode};
use serde::{Deserialize,Serialize};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool,Ordering};
use std::sync::Arc;
use std::time::{Duration,Instant};

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResourceLimits{
    pub wall_time_s:f64,pub total_calculator_evaluations:u64,pub evaluations_per_relaxation:u64,
    pub evaluations_per_saddle_attempt:u64,pub saddle_attempts_per_state:u64,pub catalog_events:u64,
    pub resident_memory_bytes:u64,pub output_bytes:u64,pub callback_timeout_s:f64,
    pub retry_count:u32,pub retry_backoff_s:f64,
}
impl ResourceLimits{pub fn validate(&self)->Result<(),Status>{validate_limits(self)}}

#[derive(Clone,Debug,Default,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResourceCounters{
    pub calculator_evaluations:u64,pub relaxation_evaluations:BTreeMap<String,u64>,
    pub saddle_evaluations:BTreeMap<String,u64>,pub saddle_attempts_by_state:BTreeMap<String,u64>,
    pub catalog_events:u64,pub resident_memory_bytes:u64,pub output_bytes:u64,
    pub retry_history:Vec<RetryRecord>,pub baseline_wall_elapsed_s:f64,
    pub budget_overruns:Vec<BudgetOverrun>,
}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RetryRecord{pub callback_id:String,pub attempt:u32,pub retryable:bool,pub identical_request_digest:String,pub substream_digest:Option<String>,pub backoff_s:f64,pub outcome:String}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BudgetOverrun{pub component:String,pub request_id:String,pub reservation:u64,pub actual:u64,pub total_counter_after:u64,pub scoped_counter_after:u64}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResourceSnapshot{pub limits:ResourceLimits,pub counters:ResourceCounters}
impl ResourceSnapshot{pub fn validate(&self)->Result<(),Status>{validate_limits(&self.limits)?;let c=&self.counters;let l=&self.limits;let evaluations:Vec<_>=c.budget_overruns.iter().filter(|x|x.component=="relaxation"||x.component=="saddle-search").collect();let evidence=evaluations.first().copied();let overrun_invalid=c.budget_overruns.iter().any(|x|x.component.is_empty()||x.request_id.is_empty()||x.actual<=x.reservation||match x.component.as_str(){"relaxation"=>x.total_counter_after!=c.calculator_evaluations||c.relaxation_evaluations.get(&x.request_id)!=Some(&x.scoped_counter_after)||x.scoped_counter_after<x.actual,"saddle-search"=>x.total_counter_after!=c.calculator_evaluations||c.saddle_evaluations.get(&x.request_id)!=Some(&x.scoped_counter_after)||x.scoped_counter_after<x.actual,"wall-time-ms"=>x.request_id!="resource-ledger"||x.total_counter_after!=x.actual||x.scoped_counter_after!=x.actual,_=>true});let total_unexplained=c.calculator_evaluations>l.total_calculator_evaluations&&evidence.map_or(true,|x|x.total_counter_after!=c.calculator_evaluations);let relax_unexplained=c.relaxation_evaluations.iter().any(|(id,v)|*v>l.evaluations_per_relaxation&&evidence.map_or(true,|x|x.component!="relaxation"||x.request_id.as_str()!=id.as_str()||x.scoped_counter_after!=*v));let saddle_unexplained=c.saddle_evaluations.iter().any(|(id,v)|*v>l.evaluations_per_saddle_attempt&&evidence.map_or(true,|x|x.component!="saddle-search"||x.request_id.as_str()!=id.as_str()||x.scoped_counter_after!=*v));if !c.baseline_wall_elapsed_s.is_finite()||c.baseline_wall_elapsed_s<0.0||c.baseline_wall_elapsed_s>l.wall_time_s||evaluations.len()>1||overrun_invalid||total_unexplained||relax_unexplained||saddle_unexplained||c.catalog_events>l.catalog_events||c.resident_memory_bytes>l.resident_memory_bytes||c.output_bytes>l.output_bytes||c.saddle_attempts_by_state.values().any(|x|*x>l.saddle_attempts_per_state)||!c.retry_history.is_empty(){return Err(Status::simple(StatusCode::CheckpointCorrupt,"resource","RES-004","persisted resource limits/counters lack exact scoped overrun evidence"));}Ok(())}}

#[derive(Clone)] pub struct CancelToken(Arc<AtomicBool>);
impl CancelToken{pub fn new()->Self{Self(Arc::new(AtomicBool::new(false)))}pub fn cancel(&self){self.0.store(true,Ordering::SeqCst)}pub fn is_cancelled(&self)->bool{self.0.load(Ordering::SeqCst)}}
impl Default for CancelToken{fn default()->Self{Self::new()}}

pub struct ResourceLedger{pub limits:ResourceLimits,pub counters:ResourceCounters,started:Instant,cancel:CancelToken}
impl ResourceLedger{
    pub fn new(limits:ResourceLimits,counters:ResourceCounters,cancel:CancelToken)->Result<Self,Status>{
        validate_limits(&limits)?;ResourceSnapshot{limits:limits.clone(),counters:counters.clone()}.validate()?;
        Ok(Self{limits,counters,started:Instant::now(),cancel})
    }
    pub fn elapsed_s(&self)->f64{self.counters.baseline_wall_elapsed_s+self.started.elapsed().as_secs_f64()}
    pub fn check(&self)->Result<(),Status>{
        if self.cancel.is_cancelled(){return Err(Status::simple(StatusCode::Cancelled,"resource","CKPT-005","cancellation requested"));}
        if self.elapsed_s()>self.limits.wall_time_s{return Err(limit("wall time exhausted"));}Ok(())
    }
    pub fn callback_deadline(&self)->Result<Instant,Status>{self.check()?;let remaining=(self.limits.wall_time_s-self.elapsed_s()).max(0.0);Ok(Instant::now()+Duration::from_secs_f64(remaining.min(self.limits.callback_timeout_s)))}
    pub fn reserve_calculator(&mut self,n:u64,relaxation:Option<&str>,saddle:Option<&str>)->Result<(),Status>{
        self.check()?;let total=self.counters.calculator_evaluations.checked_add(n).ok_or_else(||limit("calculator counter overflow"))?;
        if total>self.limits.total_calculator_evaluations{return Err(limit("total calculator evaluation limit"));}
        let relaxation_next=if let Some(id)=relaxation{let v=self.counters.relaxation_evaluations.get(id).copied().unwrap_or(0).checked_add(n).ok_or_else(||limit("relaxation counter overflow"))?;if v>self.limits.evaluations_per_relaxation{return Err(limit("per-relaxation evaluation limit"));}Some((id.to_owned(),v))}else{None};let saddle_next=if let Some(id)=saddle{let v=self.counters.saddle_evaluations.get(id).copied().unwrap_or(0).checked_add(n).ok_or_else(||limit("saddle counter overflow"))?;if v>self.limits.evaluations_per_saddle_attempt{return Err(limit("per-saddle evaluation limit"));}Some((id.to_owned(),v))}else{None};self.counters.calculator_evaluations=total;if let Some((id,v))=relaxation_next{self.counters.relaxation_evaluations.insert(id,v);}if let Some((id,v))=saddle_next{self.counters.saddle_evaluations.insert(id,v);}Ok(())
    }
    pub fn refund_calculator(&mut self,n:u64,relaxation:Option<&str>,saddle:Option<&str>)->Result<(),Status>{
        self.counters.calculator_evaluations=self.counters.calculator_evaluations.checked_sub(n).ok_or_else(||ledger_underflow("calculator ledger underflow"))?;
        if let Some(id)=relaxation{let v=self.counters.relaxation_evaluations.get(id).copied().unwrap_or(0).checked_sub(n).ok_or_else(||ledger_underflow("relaxation ledger underflow"))?;self.counters.relaxation_evaluations.insert(id.to_owned(),v);}
        if let Some(id)=saddle{let v=self.counters.saddle_evaluations.get(id).copied().unwrap_or(0).checked_sub(n).ok_or_else(||ledger_underflow("saddle ledger underflow"))?;self.counters.saddle_evaluations.insert(id.to_owned(),v);}Ok(())
    }
    pub fn finish_evaluation_reservation(&mut self,reserved:u64,actual:u64,component:&str,relaxation:Option<&str>,saddle:Option<&str>)->Result<(),Status>{let(request_id,is_relax)=match(relaxation,saddle,component){(Some(id),None,"relaxation")=>(id,true),(None,Some(id),"saddle-search")=>(id,false),_=>return Err(ledger_underflow("evaluation reservation scope/request binding invalid"))};if actual<=reserved{return self.refund_calculator(reserved-actual,relaxation,saddle)}let extra=actual-reserved;let total=self.counters.calculator_evaluations.checked_add(extra).ok_or_else(||limit("actual calculator counter overflow"))?;let scoped=if is_relax{self.counters.relaxation_evaluations.get(request_id).copied().unwrap_or(0).checked_add(extra).ok_or_else(||limit("actual relaxation counter overflow"))?}else{self.counters.saddle_evaluations.get(request_id).copied().unwrap_or(0).checked_add(extra).ok_or_else(||limit("actual saddle counter overflow"))?};self.counters.calculator_evaluations=total;if is_relax{self.counters.relaxation_evaluations.insert(request_id.to_owned(),scoped);}else{self.counters.saddle_evaluations.insert(request_id.to_owned(),scoped);}self.counters.budget_overruns.push(BudgetOverrun{component:component.into(),request_id:request_id.into(),reservation:reserved,actual,total_counter_after:total,scoped_counter_after:scoped});Err(limit("callback exceeded pre-reserved evaluation budget; exact actual counts committed"))}
    pub fn reserve_saddle_attempt(&mut self,state_id:&str)->Result<(),Status>{self.check()?;let v=self.counters.saddle_attempts_by_state.get(state_id).copied().unwrap_or(0).checked_add(1).ok_or_else(||limit("saddle-attempt counter overflow"))?;if v>self.limits.saddle_attempts_per_state{return Err(limit("per-state saddle-attempt limit"));}self.counters.saddle_attempts_by_state.insert(state_id.to_owned(),v);Ok(())}
    pub fn reserve_catalog_event(&mut self,directed:u64)->Result<(),Status>{self.check()?;let v=self.counters.catalog_events.checked_add(directed).ok_or_else(||limit("catalog counter overflow"))?;if v>self.limits.catalog_events{return Err(limit("catalog event limit"));}self.counters.catalog_events=v;Ok(())}
    pub fn reserve_memory(&mut self,bytes:u64)->Result<(),Status>{self.check()?;let v=self.counters.resident_memory_bytes.checked_add(bytes).ok_or_else(||limit("memory counter overflow"))?;if v>self.limits.resident_memory_bytes{return Err(limit("resident memory limit"));}self.counters.resident_memory_bytes=v;Ok(())}
    pub fn release_memory(&mut self,bytes:u64)->Result<(),Status>{self.counters.resident_memory_bytes=self.counters.resident_memory_bytes.checked_sub(bytes).ok_or_else(||ledger_underflow("memory ledger underflow"))?;Ok(())}
    pub fn reserve_output(&mut self,bytes:u64)->Result<(),Status>{self.check()?;let v=self.counters.output_bytes.checked_add(bytes).ok_or_else(||limit("output counter overflow"))?;if v>self.limits.output_bytes{return Err(limit("output byte limit"));}self.counters.output_bytes=v;Ok(())}
    pub fn refund_output(&mut self,bytes:u64)->Result<(),Status>{self.counters.output_bytes=self.counters.output_bytes.checked_sub(bytes).ok_or_else(||ledger_underflow("output ledger underflow"))?;Ok(())}
    pub fn checkpoint_counters(&self)->ResourceCounters{let mut c=self.counters.clone();let actual=self.elapsed_s();c.baseline_wall_elapsed_s=actual.min(self.limits.wall_time_s);if actual>self.limits.wall_time_s{let reservation=(self.limits.wall_time_s*1000.0).floor().max(0.0)as u64;let observed=(actual*1000.0).ceil().max(reservation as f64+1.0)as u64;if !c.budget_overruns.iter().any(|x|x.component=="wall-time-ms"){c.budget_overruns.push(BudgetOverrun{component:"wall-time-ms".into(),request_id:"resource-ledger".into(),reservation,actual:observed,total_counter_after:observed,scoped_counter_after:observed});}}c}
    pub fn snapshot(&self)->ResourceSnapshot{ResourceSnapshot{limits:self.limits.clone(),counters:self.checkpoint_counters()}}
}
fn validate_limits(l:&ResourceLimits)->Result<(),Status>{let times=[l.wall_time_s,l.callback_timeout_s,l.retry_backoff_s];if times.iter().any(|x|!x.is_finite())||l.wall_time_s<=0.0||l.callback_timeout_s<=0.0||l.retry_backoff_s!=0.0||l.retry_count!=0||l.total_calculator_evaluations==0||l.evaluations_per_relaxation==0||l.evaluations_per_saddle_attempt==0||l.saddle_attempts_per_state==0||l.catalog_events==0||l.resident_memory_bytes==0||l.output_bytes==0{return Err(invalid("resource limits must be finite/positive; retry_count and backoff must be 0"));}Ok(())}
fn invalid(m:&str)->Status{Status::simple(StatusCode::InvalidInput,"resource","RES-001",m)}
fn limit(m:&str)->Status{Status::simple(StatusCode::ResourceLimit,"resource","RES-002",m)}
fn ledger_underflow(m:&str)->Status{Status::simple(StatusCode::InternalError,"resource","RES-001",m)}
