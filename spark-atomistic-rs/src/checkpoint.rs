// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::basin::BasinCheckpoint;
use crate::catalog::{Catalog,CatalogSnapshot,RateCache};
use crate::discovery::DiscoveryStats;
use crate::identity::hex_sha256;
use crate::kmc::{KmcStepRecord,KmcTrajectory};
use crate::model::CommittedState;
use crate::resource::{ResourceSnapshot,ResourceLedger};
use crate::rate::neumaier_sum;
use crate::rng::{validate_state,Philox,RngRecord};
use crate::status::{Status,StatusCode};
use serde::de::{self,MapAccess,SeqAccess,Visitor};
use serde::{Deserialize,Deserializer,Serialize};
use serde_json::{Map,Number,Value};
use std::collections::BTreeSet;
use std::fmt;
use std::fs::{self,File,OpenOptions};
use std::io::{Read,Write};
use std::path::{Path,PathBuf};
use std::sync::atomic::{AtomicU64,Ordering};

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CompletionFlags{pub complete:bool,pub incomplete_catalog:bool,pub cancelled:bool,pub resource_limited:bool}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CheckpointPayload{
    pub schema:String,pub config_digest:String,pub model_digest:String,pub discovery_config_digest:String,pub initial_state:CommittedState,pub current_state:CommittedState,
    pub simulation_time_s:f64,pub step_index:u64,pub catalog:CatalogSnapshot,pub rate_cache:RateCache,pub discovery_statistics:Vec<DiscoveryStats>,
    pub rng:RngRecord,pub basin:Option<BasinCheckpoint>,pub resources:ResourceSnapshot,pub trajectory:Vec<KmcStepRecord>,
    pub log_sequence:u64,pub checkpoint_sequence:u64,pub flags:CompletionFlags,pub cleanup_reasons:Vec<String>,
}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
struct CheckpointEnvelope{payload_sha256:String,payload:CheckpointPayload}

pub fn canonical_json_bytes<T:Serialize>(value:&T)->Result<Vec<u8>,Status>{
    let raw=serde_json::to_value(value).map_err(|e|corrupt(&format!("serialization failed: {e}")))?;let sorted=sort_and_validate(raw)?;let mut out=String::new();write_canonical(&sorted,&mut out)?;Ok(out.into_bytes())
}

fn write_canonical(v:&Value,out:&mut String)->Result<(),Status>{match v{Value::Null=>out.push_str("null"),Value::Bool(x)=>out.push_str(if *x{"true"}else{"false"}),Value::String(s)=>{out.push('"');for c in s.chars(){match c{'"'=>out.push_str("\\\""),'\\'=>out.push_str("\\\\"),c if (c as u32)<=0x1f=>out.push_str(&format!("\\u{:04x}",c as u32)),_=>out.push(c)}}out.push('"')},Value::Number(n)=>out.push_str(&canonical_number(n)?),Value::Array(a)=>{out.push('[');for(i,x)in a.iter().enumerate(){if i>0{out.push(',')}write_canonical(x,out)?}out.push(']')},Value::Object(m)=>{out.push('{');let mut keys:Vec<_>=m.keys().collect();keys.sort();for(i,k)in keys.iter().enumerate(){if i>0{out.push(',')}write_canonical(&Value::String((*k).clone()),out)?;out.push(':');write_canonical(&m[*k],out)?}out.push('}')}}Ok(())}
fn canonical_number(n:&Number)->Result<String,Status>{if n.is_i64(){let x=n.as_i64().ok_or_else(||corrupt("invalid integer"))?;if x.unsigned_abs()>9_007_199_254_740_991{return Err(corrupt("integer outside portable domain"))}return Ok(x.to_string())}if n.is_u64(){let x=n.as_u64().ok_or_else(||corrupt("invalid unsigned integer"))?;if x>9_007_199_254_740_991{return Err(corrupt("integer outside portable domain"))}return Ok(x.to_string())}let x=n.as_f64().ok_or_else(||corrupt("unrepresentable binary64"))?;if !x.is_finite(){return Err(corrupt("nonfinite binary64"))}if x==0.0{return Ok("0".into())}let negative=x.is_sign_negative();let magnitude=x.abs();let raw=magnitude.to_string();let(e_mant,e_exp)=match raw.find(|c|c=='e'||c=='E'){Some(i)=>(&raw[..i],raw[i+1..].parse::<i32>().map_err(|_|corrupt("invalid binary64 exponent"))?),None=>(raw.as_str(),0)};let decimal=e_mant.find('.').unwrap_or(e_mant.len())as i32;let mut digits:String=e_mant.chars().filter(|c|*c!='.').collect();let leading=digits.bytes().take_while(|b|*b==b'0').count();digits.drain(..leading);let scale=decimal+e_exp-leading as i32;while digits.ends_with('0'){digits.pop();}if digits.is_empty(){return Ok("0".into())}let exponent=scale-1;let body=if exponent>=-6&&exponent<21&&magnitude<9_007_199_254_740_992.0{if scale<=0{format!("0.{}{}","0".repeat((-scale)as usize),digits)}else if scale as usize>=digits.len(){format!("{}{}",digits,"0".repeat(scale as usize-digits.len()))}else{format!("{}.{}",&digits[..scale as usize],&digits[scale as usize..])}}else{let mantissa=if digits.len()==1{digits}else{format!("{}.{}",&digits[..1],&digits[1..])};format!("{}e{}{}",mantissa,if exponent>=0{"+"}else{""},exponent)};Ok(if negative{format!("-{body}")}else{body})}

fn sort_and_validate(v:Value)->Result<Value,Status>{match v{
    Value::Null|Value::Bool(_)|Value::String(_)=>Ok(v),
    Value::Number(ref n)if n.as_f64().map_or(false,|x|x.is_finite())=>Ok(v),
    Value::Number(_)=>Err(corrupt("nonfinite/unrepresentable JSON number")),
    Value::Array(a)=>Ok(Value::Array(a.into_iter().map(sort_and_validate).collect::<Result<_,_>>()?)),
    Value::Object(m)=>{let mut pairs:Vec<_>=m.into_iter().collect();pairs.sort_by(|a,b|a.0.cmp(&b.0));let mut out=Map::new();for(k,v)in pairs{out.insert(k,sort_and_validate(v)?);}Ok(Value::Object(out))}
}}

pub fn parse_strict_json(bytes:&[u8])->Result<Value,Status>{
    if let Ok(text)=std::str::from_utf8(bytes){if portable_integer_violation(text){return Err(corrupt("integer outside portable domain"))}}
    let mut de=serde_json::Deserializer::from_slice(bytes);let value=NoDupValue::deserialize(&mut de).map_err(|e|corrupt(&format!("strict JSON parse failed: {e}")))?.0;de.end().map_err(|e|corrupt(&format!("trailing JSON data: {e}")))?;Ok(value)
}

/// `E2-JSON-002`: "Every integer is in `[-9007199254740991,9007199254740991]`. A syntactically valid
/// integer outside this domain returns `INVALID_INPUT` before schema validation."
///
/// The visitor-level guards below reach only literals that fit `i64`/`u64`; `serde_json` routes any
/// wider literal to `visit_f64`, which cannot see that the source token was an integer. The domain is
/// therefore decided here, on the source text, exactly as E2-JSON-002 words it: a number token with
/// no fraction and no exponent is an integer, and any other number token is an `E2-JSON-003` binary64.
/// Tokens inside strings are not numbers and are skipped.
pub fn portable_integer_violation(text:&str)->bool{
    let b=text.as_bytes();let mut i=0;let mut in_string=false;
    while i<b.len(){let c=b[i];
        if in_string{if c==b'\\'{i+=2;continue}if c==b'"'{in_string=false}i+=1;continue}
        if c==b'"'{in_string=true;i+=1;continue}
        if c!=b'-'&&!c.is_ascii_digit(){i+=1;continue}
        if c==b'-'{i+=1}
        let start=i;while i<b.len()&&b[i].is_ascii_digit(){i+=1}
        let digits=&b[start..i];let mut integer_form=true;
        if i<b.len()&&b[i]==b'.'{integer_form=false;i+=1;while i<b.len()&&b[i].is_ascii_digit(){i+=1}}
        if i<b.len()&&(b[i]==b'e'||b[i]==b'E'){integer_form=false;i+=1;if i<b.len()&&(b[i]==b'+'||b[i]==b'-'){i+=1}while i<b.len()&&b[i].is_ascii_digit(){i+=1}}
        if integer_form&&!digits.is_empty()&&integer_digits_out_of_domain(digits){return true}
    }
    false
}
fn integer_digits_out_of_domain(digits:&[u8])->bool{
    let significant=match digits.iter().position(|x|*x!=b'0'){Some(k)=>&digits[k..],None=>b"0" as &[u8]};
    significant.len()>16||(significant.len()==16&&significant>b"9007199254740991" as &[u8])
}

struct NoDupValue(Value);
impl<'de>Deserialize<'de> for NoDupValue{fn deserialize<D:Deserializer<'de>>(d:D)->Result<Self,D::Error>{struct V;impl<'de>Visitor<'de>for V{type Value=NoDupValue;fn expecting(&self,f:&mut fmt::Formatter<'_>)->fmt::Result{f.write_str("strict JSON value")}
fn visit_bool<E:de::Error>(self,v:bool)->Result<Self::Value,E>{Ok(NoDupValue(Value::Bool(v)))}fn visit_i64<E:de::Error>(self,v:i64)->Result<Self::Value,E>{if v.unsigned_abs()>9_007_199_254_740_991{return Err(E::custom("integer outside portable domain"))}Ok(NoDupValue(Value::Number(v.into())))}fn visit_u64<E:de::Error>(self,v:u64)->Result<Self::Value,E>{if v>9_007_199_254_740_991{return Err(E::custom("integer outside portable domain"))}Ok(NoDupValue(Value::Number(v.into())))}fn visit_f64<E:de::Error>(self,v:f64)->Result<Self::Value,E>{if !v.is_finite(){return Err(E::custom("nonfinite number"))}Ok(NoDupValue(Value::Number(Number::from_f64(v).ok_or_else(||E::custom("unrepresentable number"))?)))}fn visit_str<E:de::Error>(self,v:&str)->Result<Self::Value,E>{Ok(NoDupValue(Value::String(v.to_owned())))}fn visit_string<E:de::Error>(self,v:String)->Result<Self::Value,E>{Ok(NoDupValue(Value::String(v)))}fn visit_none<E:de::Error>(self)->Result<Self::Value,E>{Ok(NoDupValue(Value::Null))}fn visit_unit<E:de::Error>(self)->Result<Self::Value,E>{Ok(NoDupValue(Value::Null))}fn visit_some<D2:Deserializer<'de>>(self,d:D2)->Result<Self::Value,D2::Error>{NoDupValue::deserialize(d)}
fn visit_seq<A:SeqAccess<'de>>(self,mut a:A)->Result<Self::Value,A::Error>{let mut out=Vec::new();while let Some(v)=a.next_element::<NoDupValue>()?{out.push(v.0)}Ok(NoDupValue(Value::Array(out)))}
fn visit_map<A:MapAccess<'de>>(self,mut a:A)->Result<Self::Value,A::Error>{let mut out=Map::new();let mut seen=BTreeSet::new();while let Some(k)=a.next_key::<String>()?{if !seen.insert(k.clone()){return Err(de::Error::custom(format!("duplicate key: {k}")))}let v=a.next_value::<NoDupValue>()?;out.insert(k,v.0);}Ok(NoDupValue(Value::Object(out)))}}
d.deserialize_any(V)}}

pub fn validate_payload(p:&CheckpointPayload)->Result<(),Status>{
    if p.schema!="spark-atomistic-checkpoint/1"||p.config_digest.is_empty()||p.model_digest.is_empty()||p.discovery_config_digest.is_empty(){return Err(incompatible("checkpoint schema/digests invalid"));}
    p.initial_state.validate()?;p.current_state.validate()?;p.resources.validate()?;if !p.simulation_time_s.is_finite()||p.simulation_time_s<0.0{return Err(corrupt("invalid simulation time"));}
    let catalog=Catalog::from_snapshot(p.catalog.clone()).map_err(|_|corrupt("catalog snapshot failed recursive validation"))?;p.rate_cache.validate_for(&catalog).map_err(|_|corrupt("rate cache failed epoch/digest/rate recomputation"))?;if catalog.states().get(p.initial_state.state_id())!=Some(&p.initial_state)||catalog.states().get(p.current_state.state_id())!=Some(&p.current_state){return Err(corrupt("initial/current committed state absent or inconsistent in catalog"));}validate_state(&p.rng.trajectory)?;if p.rng.substream_map.iter().any(|(k,v)|k.is_empty()||v.is_empty()){return Err(corrupt("invalid RNG substream map"));}
    if p.flags.complete&&(p.flags.incomplete_catalog||p.flags.cancelled||p.flags.resource_limited){return Err(corrupt("inconsistent completion flags"));}
    if !p.resources.counters.budget_overruns.is_empty()&&!p.flags.resource_limited{return Err(corrupt("resource overrun evidence requires resource_limited=true"));}
    if p.step_index!=p.trajectory.len() as u64||p.log_sequence!=p.trajectory.len() as u64{return Err(corrupt("step/log sequence does not match full trajectory"));}
    if p.rng.trajectory.consumed_uniforms!=p.step_index.checked_mul(2).ok_or_else(||corrupt("step-to-uniform counter overflow"))?{return Err(corrupt("trajectory RNG must consume exactly two uniforms per committed step"));}
    let mut replay_rng=Philox::new(p.rng.trajectory.key,p.rng.trajectory.initial_counter);let mut replay_state=p.initial_state.clone();let mut time=0.0;
    for(i,r)in p.trajectory.iter().enumerate(){let nums=[r.selection_uniform,r.time_uniform,r.total_rate_per_s,r.selected_rate_per_s,r.time_increment_s];if !nums.iter().all(|x|x.is_finite())||!(0.0<r.selection_uniform&&r.selection_uniform<1.0&&0.0<r.time_uniform&&r.time_uniform<1.0)||r.total_rate_per_s<=0.0||r.selected_rate_per_s<=0.0||r.time_increment_s<=0.0{return Err(corrupt("invalid trajectory record"));}if r.step_index!=i as u64+1||r.log_sequence!=i as u64+1||r.checkpoint_sequence>p.checkpoint_sequence||r.pre_state_id!=replay_state.state_id{return Err(corrupt("noncontiguous trajectory sequence/state chain"));}
        r.rate_table_snapshot.validate().map_err(|_|corrupt("historical rate snapshot hash/shape invalid"))?;let table=&r.rate_table_snapshot.payload;if table.origin_state_id!=replay_state.state_id||table.event_ids.is_empty(){return Err(corrupt("historical rate snapshot origin/events invalid"));}let rates=&table.rates;let total=neumaier_sum(rates).map_err(|_|corrupt("historical replay rate sum invalid"))?;if total.to_bits()!=table.total_rate_per_s.to_bits()||total.to_bits()!=r.total_rate_per_s.to_bits(){return Err(corrupt("historical replay total rate differs bit-exactly"));}
        let(us,ut,next_rng)=replay_rng.two_uniforms_atomic().map_err(|_|corrupt("replay RNG failed"))?;if us.to_bits()!=r.selection_uniform.to_bits()||ut.to_bits()!=r.time_uniform.to_bits(){return Err(corrupt("replay uniforms differ bit-exactly"));}let threshold=us*total;let mut cumulative=0.0;let mut chosen=rates.len()-1;for(j,rate)in rates.iter().enumerate(){cumulative+=*rate;if cumulative>threshold{chosen=j;break;}}let event=catalog.events().iter().find(|e|e.event_id==table.event_ids[chosen]).ok_or_else(||corrupt("historical selected event absent from immutable catalog"))?;let post=catalog.states().get(&table.destination_state_ids[chosen]).ok_or_else(||corrupt("historical destination state missing"))?;if event.destination_state_id!=post.state_id||event.event_id!=r.selected_event_id||rates[chosen].to_bits()!=r.selected_rate_per_s.to_bits()||post.state_id!=r.post_state_id{return Err(corrupt("historical replay selected event/rate/state differs"));}let dt=-ut.ln()/total;if !dt.is_finite()||dt<=0.0||dt.to_bits()!=r.time_increment_s.to_bits(){return Err(corrupt("replay residence time differs bit-exactly"));}time+=dt;if !time.is_finite(){return Err(corrupt("replay simulation time overflow"));}replay_state=post.clone();replay_rng=next_rng;
    }
    if time.to_bits()!=p.simulation_time_s.to_bits()||replay_state!=p.current_state||replay_rng.state()!=p.rng.trajectory{return Err(corrupt("replayed time/state/RNG does not equal checkpoint"));}
    let state_ids:BTreeSet<_>=p.catalog.states.keys().map(String::as_str).collect();for d in &p.discovery_statistics{if d.state_id.is_empty()||d.config_digest!=p.discovery_config_digest||!d.relevance_rate_window_per_s.is_finite()||d.event_log_rates.values().any(|x|!x.is_finite()){return Err(corrupt("recursive discovery statistics do not match exact checkpoint discovery-config digest"));}if !state_ids.is_empty()&&!state_ids.contains(d.state_id.as_str())&&d.state_id!=p.current_state.state_id{return Err(corrupt("discovery statistics reference unknown state"));}}
    for d in &p.discovery_statistics{let failures=d.failures_by_status.values().try_fold(0_u64,|a,x|a.checked_add(*x)).ok_or_else(||corrupt("discovery failure counter overflow"))?;if d.successes.checked_add(failures)!=Some(d.attempts)||d.duplicates>d.successes||d.consecutive_redundant_successes>d.duplicates||d.evaluations>p.resources.counters.calculator_evaluations||d.permanently_incomplete_catalog&&d.stopping_state!=crate::discovery::DiscoveryStoppingState::Incomplete{return Err(corrupt("discovery counters/stopping relations invalid"));}if d.event_log_rates.keys().any(|id|!p.catalog.events.iter().any(|e|&e.event_id==id)){return Err(corrupt("discovery rate references absent catalog event"));}}
    if p.flags.incomplete_catalog!=p.discovery_statistics.iter().any(|d|d.permanently_incomplete_catalog){return Err(corrupt("checkpoint incomplete-catalog flag disagrees with recursive discovery state"));}
    if let Some(b)=&p.basin{b.validate()?;}
    if p.resources.counters.retry_history.iter().any(|r|r.attempt!=0||r.retryable||r.identical_request_digest.is_empty()||!r.backoff_s.is_finite()||r.backoff_s!=0.0){return Err(corrupt("retry history invalid while retry implementation is disabled"));}Ok(())
}

pub fn resume_trajectory(p:&CheckpointPayload)->Result<KmcTrajectory,Status>{validate_payload(p)?;let catalog=Catalog::from_snapshot(p.catalog.clone())?;let rng=Philox::from_state(p.rng.trajectory.clone())?;KmcTrajectory::from_replayed(p.current_state.clone(),p.simulation_time_s,p.step_index,p.checkpoint_sequence,p.log_sequence,rng,p.trajectory.clone(),&catalog,&p.rate_cache)}

pub fn encode_checkpoint(p:&CheckpointPayload)->Result<Vec<u8>,Status>{validate_payload(p)?;let payload_bytes=canonical_json_bytes(p)?;let env=CheckpointEnvelope{payload_sha256:format!("sha256:{}",hex_sha256(&payload_bytes)),payload:p.clone()};canonical_json_bytes(&env)}
pub fn decode_checkpoint(bytes:&[u8],expected_config:&str,expected_model:&str,expected_discovery_config:&str)->Result<CheckpointPayload,Status>{let value=parse_strict_json(bytes)?;let env:CheckpointEnvelope=serde_json::from_value(value).map_err(|e|corrupt(&format!("checkpoint shape invalid: {e}")))?;let actual=format!("sha256:{}",hex_sha256(&canonical_json_bytes(&env.payload)?));if actual!=env.payload_sha256{return Err(corrupt("checkpoint payload hash mismatch"));}if env.payload.config_digest!=expected_config||env.payload.model_digest!=expected_model||env.payload.discovery_config_digest!=expected_discovery_config{return Err(incompatible("checkpoint config/model/discovery-config digest mismatch"));}validate_payload(&env.payload)?;Ok(env.payload)}

static TEMP_SEQUENCE:AtomicU64=AtomicU64::new(0);
pub fn write_crash_safe(path:&Path,p:&CheckpointPayload,has_valid_checkpoint:bool,ledger:&mut ResourceLedger)->Result<(),Status>{
    let bytes=encode_checkpoint(p)?;let reserved:u64=bytes.len().try_into().map_err(|_|Status::simple(StatusCode::ResourceLimit,"checkpoint","RES-001","checkpoint size exceeds u64"))?;let parent=path.parent().ok_or_else(||io_status("checkpoint has no parent directory",has_valid_checkpoint,None))?;let name=path.file_name().and_then(|x|x.to_str()).ok_or_else(||io_status("non-UTF8 checkpoint filename",has_valid_checkpoint,None))?;ledger.reserve_output(reserved)?;let mut temp=None;let mut file=None;for _ in 0..16{let seq=TEMP_SEQUENCE.fetch_add(1,Ordering::SeqCst);let candidate=parent.join(format!(".{name}.tmp.{}.{}",std::process::id(),seq));match OpenOptions::new().write(true).create_new(true).open(&candidate){Ok(f)=>{temp=Some(candidate);file=Some(f);break},Err(e)if e.kind()==std::io::ErrorKind::AlreadyExists=>continue,Err(e)=>{let _=ledger.refund_output(reserved);return Err(io_status(&format!("checkpoint temporary create failed: {e}"),has_valid_checkpoint,None))}}}let temp=match temp{Some(x)=>x,None=>{let _=ledger.refund_output(reserved);return Err(io_status("checkpoint unique temporary create retries exhausted",has_valid_checkpoint,None))}};let mut renamed=false;
    let result=(||->std::io::Result<()>{let mut f=file.take().expect("unique temporary file opened");f.write_all(&bytes)?;f.sync_all()?;fs::rename(&temp,path)?;renamed=true;File::open(parent)?.sync_all()?;Ok(())})();
    match result{Ok(())=>Ok(()),Err(e)=>{if !renamed{let _=ledger.refund_output(reserved);}let cleanup=if renamed{Some("atomic replace completed but parent-directory flush failed; destination retained".to_owned())}else{match fs::remove_file(&temp){Ok(())=>Some("temporary checkpoint removed after failure".to_owned()),Err(c)=>Some(format!("temporary checkpoint retained; cleanup failed: {c}"))}};Err(io_status(&format!("checkpoint write transaction failed: {e}"),has_valid_checkpoint,cleanup))}}
}
pub fn read_checkpoint(path:&Path,expected_config:&str,expected_model:&str,expected_discovery_config:&str)->Result<CheckpointPayload,Status>{let mut f=File::open(path).map_err(|e|corrupt(&format!("checkpoint open failed: {e}")))?;let mut b=Vec::new();f.read_to_end(&mut b).map_err(|e|corrupt(&format!("checkpoint read failed: {e}")))?;decode_checkpoint(&b,expected_config,expected_model,expected_discovery_config)}
pub fn resolve_path(input:&Path,p:&Path)->Result<PathBuf,Status>{let base=input.parent().ok_or_else(||Status::simple(StatusCode::InvalidInput,"io","IO-003","input path has no parent"))?;Ok(if p.is_absolute(){p.to_owned()}else{base.join(p)})}
fn corrupt(m:&str)->Status{Status::simple(StatusCode::CheckpointCorrupt,"checkpoint","CKPT-004",m)}
fn incompatible(m:&str)->Status{Status::simple(StatusCode::CheckpointIncompatible,"checkpoint","CKPT-004",m)}
fn io_status(m:&str,has_valid:bool,cleanup:Option<String>)->Status{let code=if has_valid{StatusCode::ResourceLimit}else{StatusCode::CheckpointCorrupt};let mut s=Status::simple(code,"checkpoint","ERR-004",m);if let Some(c)=cleanup{let cause=Status::simple(StatusCode::InternalError,"checkpoint-cleanup","ERR-004",&c);s.context.details.insert("cleanup_reason".into(),Value::String(c));s=s.caused_by(cause);}s}
