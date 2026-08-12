// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::status::{Status,StatusCode};
use serde::{Deserialize,Serialize};
use std::collections::{BTreeMap,BTreeSet};
use crate::rng::derive_substream;
use crate::catalog::CatalogCommitOutcome;

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SearchClass{pub name:String,pub probability:f64,pub targeted:bool}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiscoveryConfig{
    pub search_classes:Vec<SearchClass>,pub minimum_successful_searches:u64,pub consecutive_redundant_successes:u64,
    pub maximum_searches:u64,pub maximum_evaluations:u64,pub relevance_rate_window_per_s:f64,
    pub strict:bool,pub alpha:Option<f64>,pub alpha_calibration:Option<String>,
}
impl DiscoveryConfig{pub fn validate(&self)->Result<(),Status>{
    if self.search_classes.is_empty()||!self.search_classes.iter().any(|x|!x.targeted&&x.probability>0.0){return Err(invalid("at least one non-targeted search class must have nonzero probability"));}
    let mut names=BTreeSet::new();if self.search_classes.iter().any(|x|!names.insert(x.name.as_str())){return Err(invalid("duplicate search class name"));}
    if self.search_classes.iter().any(|x|x.name.is_empty()||!x.probability.is_finite()||x.probability<0.0)
        ||!self.relevance_rate_window_per_s.is_finite()||self.relevance_rate_window_per_s<=0.0||self.minimum_successful_searches==0||self.consecutive_redundant_successes==0||self.maximum_searches==0||self.maximum_evaluations==0{return Err(invalid("invalid discovery budgets, probabilities, or relevance window"));}
    let sum:f64=self.search_classes.iter().map(|x|x.probability).sum();if (sum-1.0).abs()>1e-12{return Err(invalid("search-class probabilities must sum to 1"));}
    match(self.alpha,self.alpha_calibration.as_ref()){(Some(a),Some(c))if a.is_finite()&&a>0.0&&a<=1.0&&!c.is_empty()=>{},(None,None)=>{},_=>return Err(invalid("alpha requires a declared calibration and must be in (0,1]"))};Ok(())}}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiscoveryStats{
    pub(crate) state_id:String,pub(crate) attempts:u64,pub(crate) successes:u64,pub(crate) evaluations:u64,pub(crate) failures_by_status:BTreeMap<StatusCode,u64>,
    pub(crate) duplicates:u64,pub(crate) consecutive_redundant_successes:u64,pub(crate) event_log_rates:BTreeMap<String,f64>,
    pub(crate) config_digest:String,
    pub(crate) relevance_rate_window_per_s:f64,pub(crate) alpha:Option<f64>,pub(crate) alpha_calibration:Option<String>,
    pub(crate) heuristic_confidence:Option<f64>,pub(crate) stopping_state:DiscoveryStoppingState,pub(crate) permanently_incomplete_catalog:bool,
}

pub fn choose_search_class(run_seed:&[u8],state_id:&str,search_index:u64,c:&DiscoveryConfig)->Result<String,Status>{
    c.validate()?;let mut rng=derive_substream(run_seed,state_id,"class-selection",search_index);let u=rng.next_uniform()?;let mut cumulative=0.0;for class in &c.search_classes{cumulative+=class.probability;if cumulative>u{return Ok(class.name.clone())}}Ok(c.search_classes.last().expect("validated nonempty").name.clone())
}
#[derive(Clone,Copy,Debug,Eq,PartialEq,Serialize,Deserialize)]
#[serde(rename_all="SCREAMING_SNAKE_CASE")]
pub enum DiscoveryStoppingState{Running,ConvergedHeuristic,Incomplete}

impl DiscoveryStats{
    pub fn new(state_id:String,c:&DiscoveryConfig)->Result<Self,Status>{c.validate()?;if state_id.is_empty(){return Err(invalid("empty discovery state ID"));}Ok(Self{state_id,attempts:0,successes:0,evaluations:0,failures_by_status:BTreeMap::new(),duplicates:0,consecutive_redundant_successes:0,event_log_rates:BTreeMap::new(),config_digest:discovery_config_digest(c)?,relevance_rate_window_per_s:c.relevance_rate_window_per_s,alpha:c.alpha,alpha_calibration:c.alpha_calibration.clone(),heuristic_confidence:None,stopping_state:DiscoveryStoppingState::Running,permanently_incomplete_catalog:false})}
    pub fn record_catalog(&mut self,outcome:&CatalogCommitOutcome,evaluations_used:u64,c:&DiscoveryConfig)->Result<StatusCode,Status>{
        self.validate_binding(c)?;self.require_running()?;let mut next=self.clone();let result=next.record_catalog_inner(outcome,evaluations_used,c)?;*self=next;Ok(result)
    }
    fn record_catalog_inner(&mut self,outcome:&CatalogCommitOutcome,evaluations_used:u64,c:&DiscoveryConfig)->Result<StatusCode,Status>{
        self.attempts=self.attempts.checked_add(1).ok_or_else(||limit("discovery attempt counter overflow"))?;
        self.evaluations=self.evaluations.checked_add(evaluations_used).ok_or_else(||limit("discovery evaluation counter overflow"))?;
        if !outcome.log_rate.is_finite()||outcome.event_id.is_empty()||outcome.validation_digest.is_empty(){return Err(Status::simple(StatusCode::InternalError,"discovery","DISC-008","catalog outcome lacks internal event/rate/validation binding"));}
        match outcome.status{
            StatusCode::Ok=>{self.successes=self.successes.checked_add(1).ok_or_else(||limit("discovery success counter overflow"))?;self.consecutive_redundant_successes=0;self.event_log_rates.insert(outcome.event_id.clone(),outcome.log_rate);},
            StatusCode::DuplicateEvent=>{self.successes=self.successes.checked_add(1).ok_or_else(||limit("discovery success counter overflow"))?;self.duplicates=self.duplicates.checked_add(1).ok_or_else(||limit("duplicate counter overflow"))?;if outcome.log_rate>=c.relevance_rate_window_per_s.ln(){self.consecutive_redundant_successes=self.consecutive_redundant_successes.checked_add(1).ok_or_else(||limit("redundancy counter overflow"))?;}else{self.consecutive_redundant_successes=0;}},
            _=>return Err(Status::simple(StatusCode::InternalError,"discovery","DISC-008","catalog outcome has non-catalog status")),
        }
        self.after(c,outcome.status)
    }
    pub fn record_failure(&mut self,outcome:StatusCode,evaluations_used:u64,c:&DiscoveryConfig)->Result<StatusCode,Status>{
        self.validate_binding(c)?;self.require_running()?;let mut next=self.clone();let result=next.record_failure_inner(outcome,evaluations_used,c)?;*self=next;Ok(result)
    }
    fn record_failure_inner(&mut self,outcome:StatusCode,evaluations_used:u64,c:&DiscoveryConfig)->Result<StatusCode,Status>{
        if matches!(outcome,StatusCode::Ok|StatusCode::DuplicateEvent){return Err(Status::simple(StatusCode::InternalError,"discovery","DISC-008","success requires unforgeable catalog outcome"));}self.attempts=self.attempts.checked_add(1).ok_or_else(||limit("discovery attempt counter overflow"))?;self.evaluations=self.evaluations.checked_add(evaluations_used).ok_or_else(||limit("discovery evaluation counter overflow"))?;let v=self.failures_by_status.get(&outcome).copied().unwrap_or(0).checked_add(1).ok_or_else(||limit("failure counter overflow"))?;self.failures_by_status.insert(outcome,v);self.consecutive_redundant_successes=0;self.after(c,outcome)
    }
    fn require_running(&self)->Result<(),Status>{if self.stopping_state!=DiscoveryStoppingState::Running{return Err(Status::simple(StatusCode::InvalidInput,"discovery","DISC-008","terminal discovery statistics are immutable"));}Ok(())}
    fn validate_binding(&self,c:&DiscoveryConfig)->Result<(),Status>{c.validate()?;if self.config_digest!=discovery_config_digest(c)?||self.relevance_rate_window_per_s.to_bits()!=c.relevance_rate_window_per_s.to_bits()||self.alpha.map(f64::to_bits)!=c.alpha.map(f64::to_bits)||self.alpha_calibration!=c.alpha_calibration{return Err(Status::simple(StatusCode::CheckpointIncompatible,"discovery","DISC-008","discovery statistics/config digest binding mismatch"));}Ok(())}
    pub fn attempts(&self)->u64{self.attempts}pub fn config_digest(&self)->&str{&self.config_digest}
    fn after(&mut self,c:&DiscoveryConfig,fallback:StatusCode)->Result<StatusCode,Status>{
        if self.successes>=c.minimum_successful_searches&&self.consecutive_redundant_successes>=c.consecutive_redundant_successes{
            self.stopping_state=DiscoveryStoppingState::ConvergedHeuristic;self.heuristic_confidence=c.alpha.map(|a|1.0-1.0/(a*self.consecutive_redundant_successes as f64));return Ok(StatusCode::DiscoveryConvergedHeuristic);
        }
        if self.attempts>=c.maximum_searches||self.evaluations>=c.maximum_evaluations{
            self.stopping_state=DiscoveryStoppingState::Incomplete;if !c.strict{self.permanently_incomplete_catalog=true;}return Ok(StatusCode::DiscoveryIncomplete);
        }Ok(fallback)
    }
}
fn discovery_config_digest(c:&DiscoveryConfig)->Result<String,Status>{let bytes=crate::checkpoint::canonical_json_bytes(c)?;Ok(format!("sha256:{}",crate::identity::hex_sha256(&bytes)))}
fn invalid(m:&str)->Status{Status::simple(StatusCode::InvalidInput,"discovery","DISC-001",m)}
fn limit(m:&str)->Status{Status::simple(StatusCode::ResourceLimit,"discovery","DISC-008",m)}
