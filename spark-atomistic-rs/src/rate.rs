// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use crate::model::RateModelRecord;
use crate::status::{Status,StatusCode};

pub const KB_EV_PER_K:f64=8.617_333_262_145e-5;

pub fn common_prefactor_pair(origin_ev:f64,destination_ev:f64,saddle_ev:f64,temperature_k:f64,
    prefactor_per_s:f64,barrier_tolerance_ev:f64,epsilon_db:f64)->Result<RateModelRecord,Status>{
    let values=[origin_ev,destination_ev,saddle_ev,temperature_k,prefactor_per_s,barrier_tolerance_ev,epsilon_db];
    if !values.iter().all(|x|x.is_finite())||temperature_k<=0.0||prefactor_per_s<=0.0||barrier_tolerance_ev<0.0||epsilon_db<0.0{return Err(rate_invalid("invalid rate input"));}
    let bf=saddle_ev-origin_ev;let br=saddle_ev-destination_ev;
    if bf < -barrier_tolerance_ev || br < -barrier_tolerance_ev{return Err(rate_invalid("negative barrier below tolerance"));}
    let beta=1.0/(KB_EV_PER_K*temperature_k);let log_nu=prefactor_per_s.ln();
    let lf=log_nu-beta*bf;let lr=log_nu-beta*br;
    if !lf.is_finite()||!lr.is_finite(){return Err(rate_invalid("nonfinite log rate"));}
    // NO DETAILED-BALANCE VERIFICATION HAPPENS HERE, and none can under this rate model.
    //
    // `RATE-005` tests `|ln(k_ij/k_ji) + beta(F_j - F_i)| <= epsilon_DB` with `(F=E) in
    // common-prefactor energy mode`. `RATE-002` fixes BOTH directions to one saddle and one shared
    // `nu`: `k_ij = nu e^{-beta(E_s-E_i)}`, `k_ji = nu e^{-beta(E_s-E_j)}`. `E2-RATE-001` repeats it
    // for the record: "Raw barriers are exactly `b_f=E_s-E_i` and `b_r=E_s-E_j`", and
    // `E2-EVENT-001` calls them "Raw same-saddle differences". Substituting gives
    // `ln(k_ij/k_ji) = beta(E_i - E_j)`, so the residual below is `beta(E_i-E_j) + beta(E_j-E_i)`:
    // an algebraic identity of `RATE-002`, identically zero in exact arithmetic. What the
    // comparison can still detect is binary64 ROUNDING, not physics. Measured in this crate by
    // `detailed_balance_gate_measures_rounding_not_physics` (`tests/conformance_spec.rs`): over
    // 200,000 randomized `(E_i, E_j, E_s, T, nu)` inputs the worst `|residual|` is
    // `7.275957614183426e-12` against the `E2-SCHEMA-008` default `epsilon_DB = 1e-8`, so the gate
    // is unreachable there; 90059 of those inputs have a NONZERO residual and every one of them is
    // rejected once `epsilon_DB` is pushed to `0`, which is the whole content of the test. (The
    // figure reported for the twin backend, `1.7053025658242404e-13`, is the same quantity over a
    // narrower input range; the residual scale is set by `beta * |E_i - E_j|`, so a wider energy
    // and lower temperature range raises it. Neither reaches `1e-8`.)
    //
    // The two ways out are both refused by the frozen text, so neither was taken:
    //   * computing `log_reverse` from an INDEPENDENTLY obtained reverse barrier is forbidden --
    //     `E2-RATE-001`/`E2-EVENT-001` require the raw same-saddle difference, and no other
    //     quantity may be invented to keep the gate alive;
    //   * deleting the recorded residual is forbidden -- `E2-EVENT-003` lists
    //     `detailed_balance_residual` among the exactly six `rate_model` fields, `E2-CKPT-007`(5)
    //     recomputes it, and `E2-PAR-002`(7) makes a detailed-balance-failure fixture mandatory;
    //     `RATE-006` also requires `DETAILED_BALANCE_VIOLATION` to exist as an answer.
    // The comparison is therefore kept exactly as `RATE-005` words it, and is documented here as a
    // rounding guard rather than a verification. `RATE-001` allows `HARMONIC_TST` or an externally
    // supplied free-energy model, whose reverse log rate is NOT derived from the forward one; there
    // the same test has real force. `E2-SCOPE-002` excludes those from this frozen P0 subset.
    let residual=(lf-lr)+beta*(destination_ev-origin_ev);
    if !residual.is_finite()||residual.abs()>epsilon_db{return Err(Status::simple(StatusCode::DetailedBalanceViolation,"rate","RATE-005","detailed-balance residual exceeds tolerance"));}
    Ok(RateModelRecord{model:"COMMON_PREFACTOR".to_owned(),temperature_k,common_prefactor_per_s:Some(prefactor_per_s),log_forward_rate_per_s:lf,log_reverse_rate_per_s:lr,detailed_balance_residual:residual})
}

pub fn exp_selectable(log_rate:f64)->Result<f64,Status>{
    if !log_rate.is_finite(){return Err(rate_invalid("nonfinite log rate"));}
    let r=log_rate.exp();if !r.is_finite()||r<=0.0{return Err(rate_invalid("rate overflow or underflow; no undeclared clipping/cutoff"));}Ok(r)
}

pub fn log_sum_exp(log_rates:&[f64])->Result<f64,Status>{
    if log_rates.is_empty(){return Err(Status::simple(StatusCode::NoEnabledEvent,"rate","KMC-002","no enabled event"));}
    if !log_rates.iter().all(|x|x.is_finite()){return Err(rate_invalid("nonfinite log rate in sum"));}
    let m=log_rates.iter().copied().fold(f64::NEG_INFINITY,f64::max);
    let terms:Vec<f64>=log_rates.iter().map(|x|(x-m).exp()).collect();let s=neumaier_sum(&terms)?;
    let out=m+s.ln();if !out.is_finite(){return Err(rate_invalid("nonfinite total log rate"));}Ok(out)
}

pub fn neumaier_sum(values:&[f64])->Result<f64,Status>{
    let mut sum=0.0;let mut correction=0.0;
    for &x in values{if !x.is_finite()||x<0.0{return Err(rate_invalid("invalid linear rate"));}let t=sum+x;if sum.abs()>=x.abs(){correction+=(sum-t)+x}else{correction+=(x-t)+sum}sum=t;}
    let total=sum+correction;if !total.is_finite(){return Err(rate_invalid("nonfinite compensated sum"));}Ok(total)
}
fn rate_invalid(m:&str)->Status{Status::simple(StatusCode::RateInvalid,"rate","RATE-004",m)}

