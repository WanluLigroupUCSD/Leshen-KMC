// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::identity::hex_sha256;
use crate::status::{Status, StatusCode};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;

pub const RNG_ALGORITHM: &str = "Philox4x32-10:errata-1-midpoint52";
const M0:u32=0xD251_1F53; const M1:u32=0xCD9E_8D57;
const W0:u32=0x9E37_79B9; const W1:u32=0xBB67_AE85;

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhiloxState {
    pub algorithm:String,
    pub key:[u32;2],
    pub initial_counter:[u32;4],
    pub next_counter:[u32;4],
    pub buffered_block:Option<[u32;4]>,
    pub next_pair:u8,
    pub consumed_blocks:u64,
    pub consumed_uniforms:u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Philox { state:PhiloxState }

impl Philox {
    pub fn new(key:[u32;2],counter:[u32;4])->Self{Self{state:PhiloxState{algorithm:RNG_ALGORITHM.to_owned(),key,initial_counter:counter,next_counter:counter,buffered_block:None,next_pair:0,consumed_blocks:0,consumed_uniforms:0}}}
    pub fn from_state(state:PhiloxState)->Result<Self,Status>{validate_state(&state)?;Ok(Self{state})}
    pub fn state(&self)->PhiloxState{self.state.clone()}

    pub fn next_uniform(&mut self)->Result<f64,Status>{
        if self.state.buffered_block.is_none(){
            let block=philox4x32_10(self.state.next_counter,self.state.key);
            increment_counter(&mut self.state.next_counter)?;
            self.state.consumed_blocks=self.state.consumed_blocks.checked_add(1).ok_or_else(||rng_limit("consumed-block counter overflow"))?;
            self.state.buffered_block=Some(block);self.state.next_pair=0;
        }
        let block=self.state.buffered_block.ok_or_else(||rng_corrupt("missing buffered block"))?;
        let (a,b)=match self.state.next_pair{0=>(block[0],block[1]),1=>(block[2],block[3]),_=>return Err(rng_corrupt("invalid pair index"))};
        self.state.consumed_uniforms=self.state.consumed_uniforms.checked_add(1).ok_or_else(||rng_limit("uniform counter overflow"))?;
        self.state.next_pair+=1;
        if self.state.next_pair==2{self.state.buffered_block=None;self.state.next_pair=0;}
        uniform_from_words(a,b)
    }

    pub fn two_uniforms_atomic(&self)->Result<(f64,f64,Philox),Status>{
        let mut next=self.clone();let a=next.next_uniform()?;let b=next.next_uniform()?;Ok((a,b,next))
    }
}

pub fn uniform_from_words(a:u32,b:u32)->Result<f64,Status>{let q=((a as u64)<<20)|((b as u64)>>12);let u=((2*q+1) as f64)*2f64.powi(-53);if !(u>0.0&&u<1.0){return Err(rng_corrupt("Errata-1 midpoint escaped open interval"));}Ok(u)}

pub fn validate_state(s:&PhiloxState)->Result<(),Status>{
    if s.algorithm!=RNG_ALGORITHM{return Err(rng_corrupt("RNG algorithm mismatch"));}
    let expected_blocks=s.consumed_uniforms.checked_add(1).ok_or_else(||rng_corrupt("uniform count overflow"))?/2;
    if s.consumed_blocks!=expected_blocks{return Err(rng_corrupt("consumed block/uniform relationship invalid"));}
    if add_blocks(s.initial_counter,s.consumed_blocks)?!=s.next_counter{return Err(rng_corrupt("persistent Philox counter does not equal initial counter plus consumed blocks"));}
    let odd=s.consumed_uniforms%2==1;
    if odd!=(s.buffered_block.is_some()) || (odd && s.next_pair!=1) || (!odd && s.next_pair!=0){
        return Err(rng_corrupt("counter/buffer/consumed relationship invalid"));
    }
    if odd{let previous=subtract_one(s.next_counter)?;if s.buffered_block!=Some(philox4x32_10(previous,s.key)){return Err(rng_corrupt("buffered words are not the actual previous Philox block"));}}
    Ok(())
}

pub fn derive_substream(run_seed:&[u8],state_id:&str,search_class:&str,search_index:u64)->Philox{
    let mut h=Sha256::new();h.update(b"spark-saddle-substream/1\0");h.update(run_seed);h.update([0]);h.update(state_id.as_bytes());h.update([0]);h.update(search_class.as_bytes());h.update(search_index.to_be_bytes());
    let d=h.finalize();
    let word=|i:usize|u32::from_be_bytes([d[i],d[i+1],d[i+2],d[i+3]]);
    Philox::new([word(0),word(4)],[word(8),word(12),word(16),word(20)])
}

pub fn derive_trajectory_stream(run_seed:u64)->Philox{derive_e2(b"spark-trajectory-stream/2\0",run_seed,&[])}
pub fn derive_saddle_substream(run_seed:u64,state_id:&str,search_class:&str,search_index:u64)->Result<Philox,Status>{let state_len:u32=state_id.len().try_into().map_err(|_|rng_limit("state ID exceeds len32"))?;let class_len:u32=search_class.len().try_into().map_err(|_|rng_limit("search class exceeds len32"))?;let mut rest=Vec::new();rest.extend(state_len.to_be_bytes());rest.extend(state_id.as_bytes());rest.extend(class_len.to_be_bytes());rest.extend(search_class.as_bytes());rest.extend(search_index.to_be_bytes());Ok(derive_e2(b"spark-saddle-substream/2\0",run_seed,&rest))}
fn derive_e2(label:&[u8],run_seed:u64,rest:&[u8])->Philox{let mut h=Sha256::new();h.update(label);h.update(run_seed.to_be_bytes());h.update(rest);let d=h.finalize();let word=|i:usize|u32::from_be_bytes([d[i],d[i+1],d[i+2],d[i+3]]);Philox::new([word(0),word(4)],[word(8),word(12),word(16),word(20)])}

pub fn substream_digest(p:&Philox)->String{
    let s=p.state();let bytes=crate::checkpoint::canonical_json_bytes(&s).expect("validated finite integer RNG state serializes");format!("sha256:{}",hex_sha256(&bytes))
}

pub fn philox4x32_10(mut c:[u32;4],mut k:[u32;2])->[u32;4]{
    for round in 0..10{
        let p0=(M0 as u64)*(c[0] as u64);let p1=(M1 as u64)*(c[2] as u64);
        c=[(p1>>32) as u32^c[1]^k[0],p1 as u32,(p0>>32) as u32^c[3]^k[1],p0 as u32];
        if round!=9{k[0]=k[0].wrapping_add(W0);k[1]=k[1].wrapping_add(W1);}
    }c
}

fn increment_counter(c:&mut[u32;4])->Result<(),Status>{for x in c.iter_mut(){let(v,carry)=x.overflowing_add(1);*x=v;if !carry{return Ok(());}}Err(rng_limit("Philox counter exhausted"))}
fn subtract_one(mut c:[u32;4])->Result<[u32;4],Status>{for x in &mut c{let(v,borrow)=x.overflowing_sub(1);*x=v;if !borrow{return Ok(c)}}Err(rng_corrupt("buffered block precedes zero counter"))}
fn add_blocks(c:[u32;4],blocks:u64)->Result<[u32;4],Status>{let value=(c[0]as u128)|((c[1]as u128)<<32)|((c[2]as u128)<<64)|((c[3]as u128)<<96);let out=value.checked_add(blocks as u128).ok_or_else(||rng_corrupt("persistent counter addition overflow"))?;Ok([out as u32,(out>>32)as u32,(out>>64)as u32,(out>>96)as u32])}
fn rng_corrupt(m:&str)->Status{Status::simple(StatusCode::CheckpointCorrupt,"rng","DET-001",m)}
fn rng_limit(m:&str)->Status{Status::simple(StatusCode::ResourceLimit,"rng","DET-001",m)}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RngRecord { pub trajectory:PhiloxState, pub substream_map:BTreeMap<String,String> }
