// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::model::{determinant, AtomicSystem, Cell, CommittedState, GeometryMatch, Vec3};
use crate::status::{Status, StatusCode};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use serde::{Deserialize,Serialize};

pub const IDENTITY_VERSION: &str = "spark-state-identity/1:pbc-translation-same-species-permutation";

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IdentityConfig{
    pub version:String,pub rms_angstrom:f64,pub max_angstrom:f64,pub energy_ev_per_atom:f64,
    pub saddle_rms_angstrom:f64,pub saddle_max_angstrom:f64,pub saddle_energy_ev:f64,
    pub barrier_tolerance_ev:f64,pub direction_cosine_abs_min:f64,pub detailed_balance_epsilon:f64,
    pub closest_vector_budget:u64,pub permutation_budget:u64,pub singularity_threshold:f64,
    pub rotation_invariant:bool,pub reflection_invariant:bool,
}
impl Default for IdentityConfig{fn default()->Self{Self{version:IDENTITY_VERSION.into(),rms_angstrom:1e-3,max_angstrom:5e-3,energy_ev_per_atom:1e-6,saddle_rms_angstrom:1e-3,saddle_max_angstrom:5e-3,saddle_energy_ev:1e-5,barrier_tolerance_ev:1e-10,direction_cosine_abs_min:1.0-1e-8,detailed_balance_epsilon:1e-8,closest_vector_budget:1_000_000,permutation_budget:1_000_000,singularity_threshold:1e-12,rotation_invariant:false,reflection_invariant:false}}}
impl IdentityConfig{
    pub fn validate(&self)->Result<(),Status>{let xs=[self.rms_angstrom,self.max_angstrom,self.energy_ev_per_atom,self.saddle_rms_angstrom,self.saddle_max_angstrom,self.saddle_energy_ev,self.barrier_tolerance_ev,self.direction_cosine_abs_min,self.detailed_balance_epsilon,self.singularity_threshold];if self.version!=IDENTITY_VERSION||xs.iter().any(|x|!x.is_finite())||self.rms_angstrom<0.0||self.max_angstrom<self.rms_angstrom||self.energy_ev_per_atom<0.0||self.saddle_rms_angstrom<0.0||self.saddle_max_angstrom<self.saddle_rms_angstrom||self.saddle_energy_ev<0.0||self.barrier_tolerance_ev<0.0||!(0.0..=1.0).contains(&self.direction_cosine_abs_min)||self.detailed_balance_epsilon<0.0||self.singularity_threshold<=0.0||self.closest_vector_budget==0||self.permutation_budget==0||self.rotation_invariant||self.reflection_invariant{return Err(Status::simple(StatusCode::InvalidInput,"identity","STATE-005","invalid or unsupported identity configuration"));}Ok(())}
    pub fn digest(&self)->Result<String,Status>{self.validate()?;let bytes=crate::checkpoint::canonical_json_bytes(self)?;Ok(format!("sha256:{}",hex_sha256(&bytes)))}
    pub fn match_tolerances(&self)->MatchTolerances{MatchTolerances{rms_angstrom:self.rms_angstrom,max_angstrom:self.max_angstrom,energy_ev_per_atom:self.energy_ev_per_atom,closest_vector_budget:self.closest_vector_budget,permutation_budget:self.permutation_budget,singularity_threshold:self.singularity_threshold}}
    pub fn saddle_match_tolerances(&self)->MatchTolerances{MatchTolerances{rms_angstrom:self.saddle_rms_angstrom,max_angstrom:self.saddle_max_angstrom,energy_ev_per_atom:self.energy_ev_per_atom,closest_vector_budget:self.closest_vector_budget,permutation_budget:self.permutation_budget,singularity_threshold:self.singularity_threshold}}
}

#[derive(Clone, Debug,PartialEq)]
pub struct MatchTolerances {
    pub rms_angstrom: f64,
    pub max_angstrom: f64,
    pub energy_ev_per_atom: f64,
    pub closest_vector_budget: u64,
    pub permutation_budget: u64,
    pub singularity_threshold: f64,
}

impl Default for MatchTolerances {
    fn default() -> Self {
        Self { rms_angstrom: 1e-3, max_angstrom: 5e-3, energy_ev_per_atom: 1e-6,
            closest_vector_budget: 1_000_000, permutation_budget: 1_000_000,
            singularity_threshold: 1e-12 }
    }
}

#[derive(Clone, Debug)]
pub struct ClosestVector {
    pub displacement: Vec3,
    pub lattice_shift: [i64; 3],
    pub norm2: f64,
    pub candidates_examined: u64,
}

pub fn closest_periodic_vector(
    delta: Vec3,
    cell: Cell,
    periodic: [bool; 3],
    budget: u64,
    singularity_threshold: f64,
) -> Result<ClosestVector, Status> {
    if !delta.iter().chain(cell.iter().flatten()).all(|x| x.is_finite()) {
        return Err(Status::simple(StatusCode::NonfiniteResult, "identity", "STATE-008", "nonfinite closest-vector input"));
    }
    let inv = inverse(cell, singularity_threshold)?;
    let fractional = row_mul(delta, inv);
    let mut initial = [0_i64; 3];
    for k in 0..3 { if periodic[k] { initial[k] = (-fractional[k]).round() as i64; } }
    let mut best_shift = initial;
    let mut best_vec = add_lattice(delta, cell, initial);
    let mut best2 = dot(best_vec, best_vec);
    let inv_frobenius = inv.iter().flatten().map(|x| x*x).sum::<f64>().sqrt();
    let radius = inv_frobenius * best2.sqrt() + 1e-12;
    if !radius.is_finite() || radius > i64::MAX as f64 / 4.0 {
        return Err(Status::simple(StatusCode::InvalidState, "identity", "STATE-005", "pathological cell prevents bounded closest-vector search"));
    }
    let mut lo = [0_i64; 3]; let mut hi = [0_i64; 3];
    for k in 0..3 {
        if periodic[k] {
            lo[k] = (-fractional[k] - radius).floor() as i64;
            hi[k] = (-fractional[k] + radius).ceil() as i64;
        }
    }
    let mut examined = 0_u64;
    for i in lo[0]..=hi[0] { for j in lo[1]..=hi[1] { for k in lo[2]..=hi[2] {
        examined = examined.checked_add(1).ok_or_else(|| cvp_budget_status("closest-vector counter overflow"))?;
        if examined > budget { return Err(cvp_budget_status("exact closest-vector budget exhausted")); }
        let shift = [i,j,k];
        let v = add_lattice(delta, cell, shift);
        let n2 = dot(v,v);
        if n2 < best2 || (n2 == best2 && shift < best_shift) { best2=n2; best_vec=v; best_shift=shift; }
    }}}
    Ok(ClosestVector { displacement: best_vec, lattice_shift: best_shift, norm2: best2, candidates_examined: examined })
}

fn cvp_budget_status(message: &str) -> Status {
    Status::simple(StatusCode::ResourceLimit, "identity", "STATE-005", message)
}

pub fn discrete_identity(state:&AtomicSystem,config:&IdentityConfig)->Result<String,Status>{
    state.validate()?;
    config.validate()?;
    let mut labels: Vec<String> = state.species.iter().zip(&state.movable)
        .map(|(s,m)| format!("{s}:{}{}{}", m[0] as u8,m[1] as u8,m[2] as u8)).collect();
    labels.sort();
    let cell_bits:[[u64;3];3]=state.cell_angstrom.map(|r|r.map(normalized_bits));
    let payload=crate::checkpoint::canonical_json_bytes(&(config.version.as_str(),labels,cell_bits,state.periodic,&state.electronic,&state.calculator_model_digest))?;
    Ok(format!("sha256:{}",hex_sha256(&payload)))
}

pub fn canonical_geometry_digest(state:&AtomicSystem,config:&IdentityConfig)->Result<String,Status>{state.validate()?;config.validate()?;let discrete=discrete_identity(state,config)?;let mut candidates=Vec::with_capacity(state.atom_ids.len());for anchor in 0..state.atom_ids.len(){let mut atoms=Vec::with_capacity(state.atom_ids.len());for j in 0..state.atom_ids.len(){let cv=closest_periodic_vector(sub(state.positions_angstrom[j],state.positions_angstrom[anchor]),state.cell_angstrom,state.periodic,config.closest_vector_budget,config.singularity_threshold)?;atoms.push((state.species[j].clone(),state.movable[j],cv.displacement.map(normalized_bits)));}atoms.sort();candidates.push(crate::checkpoint::canonical_json_bytes(&atoms)?);}candidates.sort();let payload=crate::checkpoint::canonical_json_bytes(&(config.digest()?,discrete,&candidates[0]))?;Ok(format!("sha256:{}",hex_sha256(&payload)))}
pub fn canonical_state_id(state:&AtomicSystem,config:&IdentityConfig,constraint_digest:&str)->Result<String,Status>{if constraint_digest.is_empty(){return Err(Status::simple(StatusCode::InvalidState,"identity","STATE-003","constraint digest is empty"));}let bytes=crate::checkpoint::canonical_json_bytes(&(canonical_geometry_digest(state,config)?,constraint_digest))?;Ok(format!("state:sha256:{}",hex_sha256(&bytes)))}
fn normalized_bits(x:f64)->u64{if x==0.0{0}else{x.to_bits()}}

pub fn match_states(a: &CommittedState, b: &CommittedState, tol: &MatchTolerances) -> Result<GeometryMatch, Status> {
    a.validate()?; b.validate()?; b.system.validate_fixed_against(&a.system)?;
    if a.constraint_digest!=b.constraint_digest{return Err(Status::simple(StatusCode::InvalidState,"identity","STATE-007","constraint digest differs"));}if a.identity_config!=b.identity_config||(*tol!=a.identity_config.match_tolerances()&&*tol!=a.identity_config.saddle_match_tolerances()){return Err(Status::simple(StatusCode::CatalogIncompatible,"identity","STATE-007","caller tolerance is not the committed canonical identity configuration"));}
    if discrete_identity(&a.system,&a.identity_config)? != discrete_identity(&b.system,&b.identity_config)? {
        return Err(Status::simple(StatusCode::InvalidState, "identity", "STATE-006", "discrete identities differ"));
    }
    let n = a.system.species.len();
    let allowed_energy = tol.energy_ev_per_atom * n as f64;
    let energy_difference = (a.energy_ev-b.energy_ev).abs();
    if energy_difference > allowed_energy {
        return Err(Status::simple(StatusCode::InvalidState, "identity", "STATE-007", "energy tolerance failed"));
    }
    let mut best: Option<(f64,f64,Vec<usize>,Vec3)> = None;
    let anchor_a = 0;
    let candidates_b: Vec<usize> = (0..n).filter(|&j| same_label(&a.system, anchor_a, &b.system, j)).collect();
    let mut work = 0_u64;
    for anchor_b in candidates_b {
        let translation = sub(a.system.positions_angstrom[anchor_a], b.system.positions_angstrom[anchor_b]);
        let mut costs = vec![vec![f64::INFINITY; n]; n];
        for i in 0..n { for j in 0..n {
            if same_label(&a.system,i,&b.system,j) {
                let shifted = add(b.system.positions_angstrom[j], translation);
                costs[i][j] = closest_periodic_vector(sub(shifted,a.system.positions_angstrom[i]),
                    a.system.cell_angstrom,a.system.periodic,tol.closest_vector_budget,tol.singularity_threshold)?.norm2;
            }
        }}
        let (_, _, map) = exact_assignment(&costs, tol.permutation_budget, &mut work)?;
        let (translation,sum,max2)=refine_translation(&a.system,&b.system,&map,translation,tol)?;
        match &best {
            None => best=Some((sum,max2,map,translation)),
            Some((bs,bm,bmap,_)) if sum < *bs || (sum == *bs && (max2 < *bm || (max2 == *bm && map < *bmap))) =>
                best=Some((sum,max2,map,translation)),
            _ => {}
        }
    }
    let (sum,max2,mapping,translation) = best.ok_or_else(|| Status::simple(StatusCode::InvalidState,"identity","STATE-007","no same-species mapping"))?;
    let rms = (sum/n as f64).sqrt(); let max = max2.sqrt();
    if rms > tol.rms_angstrom || max > tol.max_angstrom {
        return Err(Status::simple(StatusCode::InvalidState,"identity","STATE-007","geometry tolerance failed"));
    }
    Ok(GeometryMatch { mapping_a_to_b: mapping, translation_angstrom: translation,
        rms_displacement_angstrom:rms, max_displacement_angstrom:max,
        energy_difference_ev:energy_difference, identity_version:IDENTITY_VERSION.to_owned() })
}

fn refine_translation(a:&AtomicSystem,b:&AtomicSystem,map:&[usize],mut translation:Vec3,tol:&MatchTolerances)->Result<(Vec3,f64,f64),Status>{
    for _ in 0..8{let mut mean=[0.0;3];for(i,&j)in map.iter().enumerate(){let shifted=add(b.positions_angstrom[j],translation);let cv=closest_periodic_vector(sub(shifted,a.positions_angstrom[i]),a.cell_angstrom,a.periodic,tol.closest_vector_budget,tol.singularity_threshold)?;for k in 0..3{mean[k]+=cv.displacement[k];}}for x in &mut mean{*x/=map.len() as f64;}translation=[translation[0]-mean[0],translation[1]-mean[1],translation[2]-mean[2]];if dot(mean,mean)<=1e-30{break}}
    let mut sum=0.0;let mut max2=0.0_f64;for(i,&j)in map.iter().enumerate(){let cv=closest_periodic_vector(sub(add(b.positions_angstrom[j],translation),a.positions_angstrom[i]),a.cell_angstrom,a.periodic,tol.closest_vector_budget,tol.singularity_threshold)?;sum+=cv.norm2;max2=max2.max(cv.norm2);}Ok((translation,sum,max2))
}

fn exact_assignment(costs: &[Vec<f64>], budget: u64, work: &mut u64) -> Result<(f64,f64,Vec<usize>),Status> {
    fn rec(row:usize,costs:&[Vec<f64>],used:&mut[bool],map:&mut Vec<usize>,sum:f64,max2:f64,
        best:&mut Option<(f64,f64,Vec<usize>)>,budget:u64,work:&mut u64)->Result<(),Status>{
        *work=work.checked_add(1).ok_or_else(|| cvp_budget_status("permutation counter overflow"))?;
        if *work>budget{return Err(cvp_budget_status("joint translation/permutation budget exhausted"));}
        if best.as_ref().map_or(false,|x|sum>x.0){return Ok(());}
        if row==costs.len(){
            let candidate=(sum,max2,map.clone());
            if best.as_ref().map_or(true,|x|candidate.0<x.0 || (candidate.0==x.0 && (candidate.1<x.1 || (candidate.1==x.1 && candidate.2<x.2)))){*best=Some(candidate);}
            return Ok(());
        }
        for col in 0..costs.len(){if !used[col] && costs[row][col].is_finite(){used[col]=true;map.push(col);
            rec(row+1,costs,used,map,sum+costs[row][col],max2.max(costs[row][col]),best,budget,work)?;
            map.pop();used[col]=false;}}
        Ok(())
    }
    let mut best=None; let mut used=vec![false;costs.len()]; let mut map=Vec::new();
    rec(0,costs,&mut used,&mut map,0.0,0.0,&mut best,budget,work)?;
    best.ok_or_else(|| Status::simple(StatusCode::InvalidState,"identity","STATE-007","no valid assignment"))
}

fn same_label(a:&AtomicSystem,i:usize,b:&AtomicSystem,j:usize)->bool { a.species[i]==b.species[j] && a.movable[i]==b.movable[j] }
fn dot(a:Vec3,b:Vec3)->f64{a[0]*b[0]+a[1]*b[1]+a[2]*b[2]}
fn add(a:Vec3,b:Vec3)->Vec3{[a[0]+b[0],a[1]+b[1],a[2]+b[2]]}
fn sub(a:Vec3,b:Vec3)->Vec3{[a[0]-b[0],a[1]-b[1],a[2]-b[2]]}
fn add_lattice(d:Vec3,c:Cell,n:[i64;3])->Vec3{[
    d[0]+n[0] as f64*c[0][0]+n[1] as f64*c[1][0]+n[2] as f64*c[2][0],
    d[1]+n[0] as f64*c[0][1]+n[1] as f64*c[1][1]+n[2] as f64*c[2][1],
    d[2]+n[0] as f64*c[0][2]+n[1] as f64*c[1][2]+n[2] as f64*c[2][2]]}
fn row_mul(v:Vec3,m:Cell)->Vec3{[
    v[0]*m[0][0]+v[1]*m[1][0]+v[2]*m[2][0],
    v[0]*m[0][1]+v[1]*m[1][1]+v[2]*m[2][1],
    v[0]*m[0][2]+v[1]*m[1][2]+v[2]*m[2][2]]}

fn inverse(m:Cell,threshold:f64)->Result<Cell,Status>{
    let d=determinant(m); let scale=m.iter().flatten().fold(0.0_f64,|a,x|a.max(x.abs())).max(1.0);
    if !d.is_finite() || d.abs()<=threshold*scale.powi(3){return Err(Status::simple(StatusCode::InvalidState,"identity","STATE-005","singular or pathological triclinic cell"));}
    let inv=[
        [(m[1][1]*m[2][2]-m[1][2]*m[2][1])/d,(m[0][2]*m[2][1]-m[0][1]*m[2][2])/d,(m[0][1]*m[1][2]-m[0][2]*m[1][1])/d],
        [(m[1][2]*m[2][0]-m[1][0]*m[2][2])/d,(m[0][0]*m[2][2]-m[0][2]*m[2][0])/d,(m[0][2]*m[1][0]-m[0][0]*m[1][2])/d],
        [(m[1][0]*m[2][1]-m[1][1]*m[2][0])/d,(m[0][1]*m[2][0]-m[0][0]*m[2][1])/d,(m[0][0]*m[1][1]-m[0][1]*m[1][0])/d]];
    Ok(inv)
}

pub fn hex_sha256(bytes:&[u8])->String{let digest=Sha256::digest(bytes);digest.iter().map(|b|format!("{b:02x}")).collect()}

#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnvironmentDefinition {
    pub center_selection:String,pub radial_extent_angstrom:f64,pub neighbor_rule:String,
    pub ambiguity_band_angstrom:f64,pub element_labels:bool,pub periodic_image_rule:String,
    pub rotation_invariant:bool,pub reflection_invariant:bool,pub identity_version:String,
    pub closest_vector_budget:u64,pub singularity_threshold:f64,
}

pub fn environment_key(state:&AtomicSystem,def:&EnvironmentDefinition,center:usize)->Result<String,Status>{
    if center>=state.positions_angstrom.len() || !def.radial_extent_angstrom.is_finite()
        || !def.ambiguity_band_angstrom.is_finite() || !def.singularity_threshold.is_finite() || def.radial_extent_angstrom<=0.0 || def.ambiguity_band_angstrom<0.0 || def.closest_vector_budget==0 || def.singularity_threshold<=0.0 || def.center_selection.is_empty() || def.neighbor_rule.is_empty() || def.periodic_image_rule.is_empty() || def.identity_version.is_empty() {
        return Err(Status::simple(StatusCode::InvalidInput,"environment","ENV-001","invalid environment definition"));
    }
    let definition_digest=format!("sha256:{}",hex_sha256(&crate::checkpoint::canonical_json_bytes(def)?));let mut neighbors:BTreeMap<String,Vec<u64>>=BTreeMap::new();
    for j in 0..state.positions_angstrom.len(){if j!=center{
        let cv=closest_periodic_vector(sub(state.positions_angstrom[j],state.positions_angstrom[center]),state.cell_angstrom,state.periodic,def.closest_vector_budget,def.singularity_threshold)?;
        let r=cv.norm2.sqrt();
        if (r-def.radial_extent_angstrom).abs()<=def.ambiguity_band_angstrom{return Err(Status::simple(StatusCode::EnvironmentAmbiguous,"environment","ENV-004","neighbor lies in ambiguity band"));}
        if r<def.radial_extent_angstrom { neighbors.entry(if def.element_labels{state.species[j].clone()}else{"*".into()}).or_default().push(normalized_bits(r)); }
    }}
    for v in neighbors.values_mut(){v.sort();}
    let payload=crate::checkpoint::canonical_json_bytes(&(definition_digest,if def.element_labels{state.species[center].as_str()}else{"*"},neighbors))?;Ok(format!("sha256:{}",hex_sha256(&payload)))
}

#[derive(Clone,Debug,Default)]
pub struct EnvironmentReuseRegistry{pub statistics:EnvironmentStatistics,refinement_by_key:BTreeMap<String,u64>,collision_provenance:BTreeMap<String,Vec<String>>}
#[derive(Clone,Debug,Default,serde::Serialize,serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnvironmentStatistics{pub key_count:u64,pub reuse_attempts:u64,pub successful_reconvergences:u64,pub rejected_mappings:u64,pub ambiguous_identities:u64,pub fresh_searches:u64}
impl EnvironmentReuseRegistry{
    pub fn record_key(&mut self,key:&str)->Result<(),Status>{if !self.refinement_by_key.contains_key(key){self.statistics.key_count=self.statistics.key_count.checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"environment","ENV-008","key counter overflow"))?;self.refinement_by_key.insert(key.to_owned(),0);}Ok(())}
    pub fn record_ambiguous(&mut self)->Result<(),Status>{self.statistics.ambiguous_identities=self.statistics.ambiguous_identities.checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"environment","ENV-008","ambiguity counter overflow"))?;Ok(())}
    pub fn split_collision(&mut self,key:&str,reason:String)->Result<String,Status>{self.record_key(key)?;let generation=self.refinement_by_key.get(key).copied().unwrap_or(0).checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"environment","ENV-006","refinement version overflow"))?;self.refinement_by_key.insert(key.to_owned(),generation);self.collision_provenance.entry(key.to_owned()).or_default().push(reason);Ok(format!("{key}:refinement-{generation}"))}
    pub fn begin_reuse(&mut self)->Result<(),Status>{self.statistics.reuse_attempts=self.statistics.reuse_attempts.checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"environment","ENV-008","reuse counter overflow"))?;Ok(())}
    pub fn finish_reuse(&mut self,reconverged_and_validated:bool)->Result<(),Status>{let x=if reconverged_and_validated{&mut self.statistics.successful_reconvergences}else{&mut self.statistics.rejected_mappings};*x=x.checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"environment","ENV-008","reuse result counter overflow"))?;Ok(())}
    pub fn record_fresh(&mut self)->Result<(),Status>{self.statistics.fresh_searches=self.statistics.fresh_searches.checked_add(1).ok_or_else(||Status::simple(StatusCode::ResourceLimit,"environment","ENV-008","fresh-search counter overflow"))?;Ok(())}
}
