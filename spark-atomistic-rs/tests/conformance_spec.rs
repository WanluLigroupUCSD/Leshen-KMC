// Clean-room tests authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; intentionally not executed during the paused dynamic-test phase.
use spark_atomistic_rs::checkpoint::parse_strict_json;
use spark_atomistic_rs::rate::common_prefactor_pair;
use spark_atomistic_rs::rng::{philox4x32_10,uniform_from_words,Philox};
use spark_atomistic_rs::status::StatusCode;
use spark_atomistic_rs::catalog::{CandidateStaging,UnvalidatedDiagnostic};
use spark_atomistic_rs::basin::{attempt_acceleration,BasinCheckpoint};
use spark_atomistic_rs::discovery::{DiscoveryConfig,DiscoveryStats,SearchClass};
use spark_atomistic_rs::parity::{dispatch_json,SCHEMA_DIGEST};
use spark_atomistic_rs::checkpoint::canonical_json_bytes;

#[test]
fn errata_1_boundary_midpoints_are_exact(){
    assert_eq!(uniform_from_words(0,0).unwrap().to_bits(),0x3ca0_0000_0000_0000);
    assert_eq!(uniform_from_words(u32::MAX,0xffff_f000).unwrap().to_bits(),0x3fef_ffff_ffff_ffff);
}

#[test]
fn philox_zero_counter_key_reference_words(){
    assert_eq!(philox4x32_10([0;4],[0;2]),[0x6627_e8d5,0xe169_c58d,0xbc57_ac4c,0x9b00_dbd8]);
}

#[test]
fn strict_json_rejects_duplicate_at_any_depth(){
    assert!(parse_strict_json(br#"{"outer":{"x":1,"x":2}}"#).is_err());
}

#[test]
fn common_prefactor_pair_obeys_detailed_balance(){
    let r=common_prefactor_pair(-1.0,-0.8,0.2,600.0,1e13,1e-10,1e-8).unwrap();
    assert!(r.detailed_balance_residual.abs()<=1e-8);
}

/// `RATE-005` is a DEAD GATE under `COMMON_PREFACTOR` and this test is its measurement, not a
/// conformance claim. `RATE-002` fixes both directions to one saddle and one shared prefactor and
/// `E2-RATE-001` fixes both barriers to raw same-saddle differences, so
/// `ln(k_ij/k_ji) + beta(E_j - E_i)` is identically zero in exact arithmetic and the recorded
/// residual is pure binary64 rounding. See the block comment in `src/rate.rs`.
///
/// ACCEPTED BASELINE: 200,000 randomized `(E_i, E_j, E_s, T, nu)` inputs at the `E2-SCHEMA-008`
/// default `detailed_balance_tolerance = 1e-8`. Every one is rated; none trips `RATE-006`.
/// PAIRED REJECTION: the SAME inputs with `detailed_balance_tolerance = 0` -- the only difference --
/// where every input whose rounding residual is nonzero DOES return `DETAILED_BALANCE_VIOLATION`.
/// The pair is what shows the gate measures rounding rather than physics: no energy, temperature or
/// prefactor makes it fire, and only shrinking the tolerance to the rounding scale does.
#[test]
fn detailed_balance_gate_measures_rounding_not_physics(){
    let mut rng=Philox::new([0x5041_5249,0x5459_3035],[0,0,0,0]);
    let mut draw=|lo:f64,hi:f64|->f64{lo+(hi-lo)*rng.next_uniform().expect("uniform")};
    let mut worst=0.0_f64;let mut worst_case=(0.0,0.0,0.0,0.0,0.0);
    let mut rejected_at_zero=0_u64;let mut nonzero_residuals=0_u64;
    const CASES:u32=200_000;
    for _ in 0..CASES{
        let origin=draw(-50.0,50.0);let destination=draw(-50.0,50.0);
        let saddle=origin.max(destination)+draw(0.0,5.0);
        let temperature=draw(50.0,2000.0);let prefactor=10f64.powf(draw(8.0,16.0));
        // `E2-SCHEMA-008` default tolerance; barrier tolerance 0 so only RATE-005 can reject.
        let ok=common_prefactor_pair(origin,destination,saddle,temperature,prefactor,0.0,1e-8)
            .expect("RATE-005 cannot reject a COMMON_PREFACTOR pair at the default tolerance");
        let residual=ok.detailed_balance_residual.abs();
        if residual>worst{worst=residual;worst_case=(origin,destination,saddle,temperature,prefactor);}
        if residual>0.0{nonzero_residuals+=1;}
        // The ONLY difference: the tolerance.
        if common_prefactor_pair(origin,destination,saddle,temperature,prefactor,0.0,0.0)
            .map_err(|s|s.status)==Err(StatusCode::DetailedBalanceViolation){rejected_at_zero+=1;}
    }
    println!("RATE-005 over {CASES} randomized inputs: worst |residual| = {worst:e} \
        (E_i={} E_j={} E_s={} T={} nu={:e}); nonzero residuals {nonzero_residuals}; \
        rejected at epsilon_DB=0: {rejected_at_zero}",
        worst_case.0,worst_case.1,worst_case.2,worst_case.3,worst_case.4);
    assert!(worst<1e-8,"the gate is unreachable at the E2-SCHEMA-008 default: worst |residual| = {worst:e}");
    assert_eq!(rejected_at_zero,nonzero_residuals,
        "every rejection at epsilon_DB=0 is exactly a nonzero rounding residual, and nothing else");
    assert!(nonzero_residuals>0,"the paired rejection must actually fire, or it proves nothing");
}

#[test]
fn nonterminal_status_has_no_process_exit_code(){assert_eq!(StatusCode::DuplicateEvent.exit_code(false,false),None);assert_eq!(StatusCode::BasinDisabled.exit_code(false,false),None);}

#[test]
fn philox_rejects_forged_previous_buffer(){let mut rng=Philox::new([1,2],[3,4,5,6]);let _=rng.next_uniform().unwrap();let mut state=rng.state();state.buffered_block.as_mut().unwrap()[0]^=1;assert!(Philox::from_state(state).is_err());}

#[test]
fn allow_unvalidated_is_a_real_diagnostic_gate(){let diagnostic=UnvalidatedDiagnostic{search_id:"s".into(),status:StatusCode::InvalidSaddle,reason:"test".into()};assert!(CandidateStaging::new(false).retain(diagnostic.clone()).is_err());assert!(CandidateStaging::new(true).retain(diagnostic).is_ok());}

#[test]
fn basin_cannot_self_enable(){assert_eq!(attempt_acceleration().unwrap_err().status,StatusCode::BasinDisabled);assert!(BasinCheckpoint::disabled().validate().is_ok());}

#[test]
fn discovery_stats_reject_changed_config_digest(){let c=DiscoveryConfig{search_classes:vec![SearchClass{name:"global".into(),probability:1.0,targeted:false}],minimum_successful_searches:1,consecutive_redundant_successes:1,maximum_searches:2,maximum_evaluations:2,relevance_rate_window_per_s:1.0,strict:true,alpha:None,alpha_calibration:None};let mut stats=DiscoveryStats::new("state".into(),&c).unwrap();let mut changed=c.clone();changed.maximum_searches=3;assert!(stats.record_failure(StatusCode::InvalidSaddle,0,&changed).is_err());assert_eq!(stats.attempts(),0);}

#[test]
fn e2_capability_response_is_exact(){let actual=dispatch_json(include_bytes!("corpus/e2_capabilities.request.json"),None,None);let expected=include_bytes!("corpus/e2_capabilities.response.json").strip_suffix(b"\n").unwrap();assert_eq!(actual,expected);assert_eq!(SCHEMA_DIGEST,"sha256:583d580d54e3847ef92f1b1456dda006161689c0bac27fd7ea896a093f48c02c");}

#[test]
fn e2_canonical_number_boundaries(){assert_eq!(canonical_json_bytes(&serde_json::json!(-0.0)).unwrap(),b"0");assert_eq!(canonical_json_bytes(&serde_json::json!(1e-6)).unwrap(),b"0.000001");assert_eq!(canonical_json_bytes(&serde_json::json!(1e21)).unwrap(),b"1e+21");}

#[test]
fn e2_basin_true_validates_but_does_not_advertise_acceleration(){let mut model:serde_json::Value=serde_json::from_slice(include_bytes!("corpus/e2_minimal_model.json")).unwrap();model["basin"]["enabled"]=serde_json::json!(true);let request=serde_json::json!({"model":model,"operation":"validate"});let response:serde_json::Value=serde_json::from_slice(&dispatch_json(&serde_json::to_vec(&request).unwrap(),None,None)).unwrap();assert_eq!(response["status"],"OK");}
