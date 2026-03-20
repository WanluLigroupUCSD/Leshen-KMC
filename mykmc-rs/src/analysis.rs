/// Analysis tools: trajectory recording, steady-state detection.

use std::collections::HashMap;
use crate::engine::KMCEngine;

/// Records simulation trajectory at regular intervals.
pub struct TrajectoryRecorder {
    pub times: Vec<f64>,
    pub steps: Vec<u64>,
    pub coverages: Vec<HashMap<String, f64>>,
    pub tofs: Vec<HashMap<String, f64>>,
}

impl TrajectoryRecorder {
    pub fn new() -> Self {
        Self {
            times: Vec::new(),
            steps: Vec::new(),
            coverages: Vec::new(),
            tofs: Vec::new(),
        }
    }

    pub fn record(&mut self, engine: &mut KMCEngine) {
        self.times.push(engine.kmc_time);
        self.steps.push(engine.kmc_step);
        self.coverages.push(engine.get_coverage());
        let tof = engine.get_tof();
        if !tof.is_empty() {
            self.tofs.push(tof);
        }
    }

    /// Get coverage time series for a specific species.
    pub fn get_coverage_series(&self, species: &str) -> Vec<f64> {
        self.coverages.iter()
            .map(|c| *c.get(species).unwrap_or(&0.0))
            .collect()
    }

    /// Get TOF time series for a specific observable.
    pub fn get_tof_series(&self, observable: &str) -> Vec<f64> {
        self.tofs.iter()
            .map(|t| *t.get(observable).unwrap_or(&0.0))
            .collect()
    }
}

/// Check if a time series has reached steady state.
pub fn check_steady_state(values: &[f64], window: usize, rtol: f64) -> bool {
    if values.len() < 2 * window {
        return false;
    }
    let n = values.len();
    let recent: f64 = values[n - window..].iter().sum::<f64>() / window as f64;
    let earlier: f64 = values[n - 2 * window..n - window].iter().sum::<f64>() / window as f64;

    if earlier.abs() < 1e-30 {
        return recent.abs() < 1e-30;
    }
    ((recent - earlier) / earlier).abs() < rtol
}

/// Run simulation until steady state.
pub fn run_to_steady_state(
    engine: &mut KMCEngine,
    max_steps: u64,
    sample_interval: u64,
    rtol: f64,
) -> TrajectoryRecorder {
    let mut recorder = TrajectoryRecorder::new();
    let mut total = 0u64;

    while total < max_steps {
        engine.do_steps(sample_interval, false);
        total += sample_interval;
        recorder.record(engine);

        if recorder.coverages.len() >= 100 {
            let cov = engine.get_coverage();
            let all_steady = cov.iter()
                .filter(|(_, &v)| v > 0.01)
                .all(|(sp, _)| {
                    let series = recorder.get_coverage_series(sp);
                    check_steady_state(&series, 50, rtol)
                });
            if all_steady {
                eprintln!("Steady state at step={}, time={:.6e} s",
                          engine.kmc_step, engine.kmc_time);
                return recorder;
            }
        }
    }
    eprintln!("Warning: steady state not reached after {} steps", max_steps);
    recorder
}
