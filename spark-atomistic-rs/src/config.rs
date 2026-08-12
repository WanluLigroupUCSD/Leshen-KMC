// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::checkpoint::{canonical_json_bytes,parse_strict_json,resolve_path};
use crate::discovery::DiscoveryConfig;
use crate::identity::hex_sha256;
use crate::identity::IdentityConfig;
use crate::model::AtomicSystem;
use crate::resource::ResourceLimits;
use crate::status::{Status,StatusCode};
use serde::{Deserialize,Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::{Component,Path,PathBuf};

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config{
    pub schema:SchemaConfig,pub system:SystemConfig,pub calculator:CalculatorConfig,pub relaxation:RelaxationConfig,
    pub saddle_search:SaddleConfig,pub discovery:DiscoveryConfig,pub kinetics:KineticsConfig,pub resources:ResourceLimits,
    pub output:OutputConfig,#[serde(default)]pub metadata:BTreeMap<String,Value>,
}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SchemaConfig{pub ir:String,pub backend:String,pub validated:bool,pub production:bool,pub release:bool}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SystemConfig{pub atoms:AtomicSystem,pub identity:IdentityConfig}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CalculatorConfig{pub adapter:String,pub callback_identity:String,pub model_name:String,pub model_version:String,pub model_digest:String,pub deterministic_expected:bool}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RelaxationConfig{pub minimizer_identity:String,pub force_tolerance_ev_per_angstrom:f64,pub step_limit:u64,pub evaluation_limit:u64}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SaddleConfig{pub searcher_identity:String,pub force_tolerance_ev_per_angstrom:f64,pub curvature_tolerance_ev_per_angstrom2:f64,pub allow_unvalidated:bool}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KineticsConfig{pub mode:String,pub temperature_k:f64,pub common_prefactor_per_s:f64,pub strict:bool,pub requested_absorbing_success:bool,pub detailed_balance_epsilon:f64,pub basin_acceleration_enabled:bool}
#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputConfig{pub checkpoint_path:String,pub trajectory_path:String,pub summary_path:String,pub overwrite:bool,pub resume:bool,pub checkpoint_every_steps:u64,pub checkpoint_every_wall_s:f64}

#[derive(Clone,Debug)]pub struct ValidatedPaths{pub checkpoint:PathBuf,pub trajectory:PathBuf,pub summary:PathBuf}

pub fn parse_config(bytes:&[u8])->Result<Config,Status>{let value=parse_strict_json(bytes).map_err(|mut s|{s.status=StatusCode::InvalidInput;s.severity=StatusCode::InvalidInput.severity();s.context.requirement_id="IO-002".into();s})?;serde_json::from_value(value).map_err(|e|Status::simple(StatusCode::InvalidInput,"config","IO-002",&format!("configuration shape invalid: {e}")))}
pub fn validate_config(config:&Config,input_path:&Path)->Result<(ValidatedPaths,String),Status>{
    if config.schema.ir!=crate::IR_SCHEMA{return Err(Status::simple(StatusCode::SchemaUnsupported,"config","IO-001","IR schema must be spark-atomistic-model/1"));}
    if config.schema.backend!=crate::BACKEND{return Err(invalid("backend must be exact string rust"));}
    if config.schema.validated||config.schema.production||config.schema.release{return Err(invalid("unvalidated implementation must declare validated/production/release=false"));}
    config.system.atoms.validate()?;config.system.identity.validate()?;if config.calculator.model_digest!=config.system.atoms.calculator_model_digest||config.calculator.adapter.is_empty()||config.calculator.callback_identity.is_empty(){return Err(invalid("calculator identity/model digest mismatch"));}
    let nums=[config.relaxation.force_tolerance_ev_per_angstrom,config.saddle_search.force_tolerance_ev_per_angstrom,config.saddle_search.curvature_tolerance_ev_per_angstrom2,config.kinetics.temperature_k,config.kinetics.common_prefactor_per_s,config.kinetics.detailed_balance_epsilon,config.output.checkpoint_every_wall_s];
    if !nums.iter().all(|x|x.is_finite())||config.relaxation.force_tolerance_ev_per_angstrom<0.0||config.saddle_search.force_tolerance_ev_per_angstrom<0.0||config.saddle_search.curvature_tolerance_ev_per_angstrom2<0.0||config.kinetics.temperature_k<=0.0||config.kinetics.common_prefactor_per_s<=0.0||config.kinetics.detailed_balance_epsilon<0.0||config.output.checkpoint_every_wall_s<=0.0||config.output.checkpoint_every_steps==0{return Err(invalid("invalid finite positive tolerance/rate/checkpoint setting"));}
    if config.kinetics.mode!="COMMON_PREFACTOR"{return Err(invalid("P0 implements only COMMON_PREFACTOR"));}
    if config.kinetics.basin_acceleration_enabled{return Err(invalid("basin_acceleration_enabled must be false until an internal completeness validator exists"));}
    config.discovery.validate()?;config.resources.validate()?;
    if [&config.output.checkpoint_path,&config.output.trajectory_path,&config.output.summary_path].iter().any(|p|p.is_empty()||Path::new(p).file_name().and_then(|x|x.to_str()).map_or(true,str::is_empty)){return Err(invalid("output paths require nonempty UTF-8 basenames"));}
    let paths=ValidatedPaths{checkpoint:normalize_output(&resolve_path(input_path,Path::new(&config.output.checkpoint_path))?)?,trajectory:normalize_output(&resolve_path(input_path,Path::new(&config.output.trajectory_path))?)?,summary:normalize_output(&resolve_path(input_path,Path::new(&config.output.summary_path))?)?};
    if same_file(&paths.checkpoint,&paths.trajectory)?||same_file(&paths.checkpoint,&paths.summary)?||same_file(&paths.trajectory,&paths.summary)?{return Err(invalid("checkpoint, trajectory, and summary paths must be distinct files"));}
    for p in [&paths.checkpoint,&paths.trajectory,&paths.summary]{if p.file_name().and_then(|x|x.to_str()).map_or(true,str::is_empty){return Err(invalid("resolved output paths require nonempty UTF-8 basenames"));}if p.exists()&&!config.output.overwrite&&!(p==&paths.checkpoint&&config.output.resume){return Err(Status::simple(StatusCode::OutputExists,"io","IO-004","output exists and overwrite is false"));}}
    if config.output.resume&&!paths.checkpoint.exists(){return Err(invalid("resume=true but checkpoint does not exist"));}
    let digest=format!("sha256:{}",hex_sha256(&canonical_json_bytes(config)?));Ok((paths,digest))
}
fn normalize_output(path:&Path)->Result<PathBuf,Status>{let mut lexical=PathBuf::new();for part in path.components(){match part{Component::CurDir=>{},Component::ParentDir=>match lexical.components().next_back(){Some(Component::Normal(_))=>{lexical.pop();},Some(Component::RootDir)|Some(Component::Prefix(_))=>{},_ if !path.is_absolute()=>lexical.push(".."),_=>{}},_=>lexical.push(part.as_os_str())}}if lexical.exists(){return std::fs::canonicalize(&lexical).map_err(|e|invalid(&format!("cannot canonicalize existing output target: {e}")));}if let(Some(parent),Some(name))=(lexical.parent(),lexical.file_name()){if parent.exists(){let canonical_parent=std::fs::canonicalize(parent).map_err(|e|invalid(&format!("cannot canonicalize existing output parent: {e}")))?;return Ok(canonical_parent.join(name));}}Ok(lexical)}
fn same_file(a:&Path,b:&Path)->Result<bool,Status>{if a==b{return Ok(true)}if !(a.exists()&&b.exists()){return Ok(false)}let ma=std::fs::metadata(a).map_err(|e|invalid(&format!("cannot inspect existing output target: {e}")))?;let mb=std::fs::metadata(b).map_err(|e|invalid(&format!("cannot inspect existing output target: {e}")))?;#[cfg(unix)]{use std::os::unix::fs::MetadataExt;Ok(ma.dev()==mb.dev()&&ma.ino()==mb.ino())}#[cfg(not(unix))]{let _=(ma,mb);Ok(false)}}
fn invalid(m:&str)->Status{Status::simple(StatusCode::InvalidInput,"config","IO-001",m)}
