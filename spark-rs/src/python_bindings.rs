//! Phase C: PyO3 bindings exposing the Rust dimer + minimize hot loop to Python.
//!
//! Design:
//!   * Rust holds the loop state and arithmetic.
//!   * Force evaluation crosses the Rust→Python boundary exactly ONCE per
//!     evaluation, via a user-supplied Python callable
//!     `force_callback(positions_ndarray) -> (energy: float, forces_ndarray)`.
//!   * The user passes ASE conventions (FORCES, not gradient); Rust negates
//!     internally to match the Calculator trait's gradient semantics.
//!
//! All PyO3 surface stays in this file so `lib.rs` doesn't need to know about
//! ndarray internals.

#![cfg(feature = "python")]

use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

use crate::offlattice::calc::Calculator;
use crate::offlattice::minimize::{minimize, FireParams};
use crate::offlattice::saddle::{find_saddle, DimerParams, DimerStatus};

// ---------------------------------------------------------------------------
//                         PYTHON-CALLBACK CALCULATOR
// ---------------------------------------------------------------------------
//
// Wraps a Python callable into the Rust Calculator trait. Each `eval()` call:
//   1. Acquires the GIL (we're on the same thread as the caller, so this is
//      effectively a no-op acquire-then-release).
//   2. Converts &[[f64; 3]] -> numpy (N, 3) array.
//   3. Calls the Python callable.
//   4. Unpacks (float, ndarray) tuple.
//   5. Negates forces to gradient.

/// Calculator that delegates to a Python callable.
struct PyCalculator {
    callback: Py<PyAny>,
    frozen: Option<Vec<bool>>,
}

impl Calculator for PyCalculator {
    fn eval(&mut self, positions: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
        Python::with_gil(|py| {
            // (N, 3) numpy array — copy so the Python side can't keep a
            // reference into our Vec.
            let n = positions.len();
            let flat: Vec<f64> = positions.iter()
                .flat_map(|p| p.iter().copied())
                .collect();
            let pos_arr = numpy::ndarray::Array2::from_shape_vec((n, 3), flat)
                .expect("positions reshape");
            let pos_py = pos_arr.into_pyarray(py);

            // Call the Python force callback. Expected return: (energy, forces).
            let result = self.callback
                .bind(py)
                .call1((pos_py.unbind(),))
                .expect("force callback raised");

            let tuple: Bound<PyTuple> = result.extract().expect("not a tuple");
            let energy: f64 = tuple.get_item(0).unwrap()
                .extract().expect("energy not a float");
            let forces_obj = tuple.get_item(1).unwrap();
            let forces_arr: PyReadonlyArray2<f64> = forces_obj
                .extract().expect("forces not a (N, 3) f64 ndarray");
            let forces_view = forces_arr.as_array();

            // Convert forces -> gradient (negate). Apply frozen mask.
            let mut grad = vec![[0.0; 3]; n];
            for i in 0..n {
                let frozen_i = self.frozen.as_ref()
                    .map(|m| m[i])
                    .unwrap_or(false);
                if frozen_i {
                    grad[i] = [0.0; 3];
                } else {
                    grad[i] = [
                        -forces_view[[i, 0]],
                        -forces_view[[i, 1]],
                        -forces_view[[i, 2]],
                    ];
                }
            }
            (energy, grad)
        })
    }

    fn frozen_mask(&self) -> Option<&[bool]> {
        self.frozen.as_deref()
    }
}

// ---------------------------------------------------------------------------
//                              PYTHON HELPERS
// ---------------------------------------------------------------------------

/// Convert a (N, 3) numpy array to Vec<[f64; 3]>.
fn pyarray_to_vec3(arr: &PyReadonlyArray2<f64>) -> Vec<[f64; 3]> {
    let view = arr.as_array();
    let n = view.shape()[0];
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        out.push([view[[i, 0]], view[[i, 1]], view[[i, 2]]]);
    }
    out
}

/// Convert Vec<[f64; 3]> to a Python (N, 3) numpy array.
fn vec3_to_pyarray<'py>(py: Python<'py>, v: &[[f64; 3]]) -> Bound<'py, PyArray2<f64>> {
    let n = v.len();
    let flat: Vec<f64> = v.iter().flat_map(|p| p.iter().copied()).collect();
    numpy::ndarray::Array2::from_shape_vec((n, 3), flat)
        .expect("vec3 reshape")
        .into_pyarray(py)
}

/// Status enum to a Python string matching the Python saddle.py convention.
fn status_to_str(s: DimerStatus) -> &'static str {
    match s {
        DimerStatus::Success => "success",
        DimerStatus::Convex => "convex",
        DimerStatus::Collision => "collision",
        DimerStatus::MaxIter => "maxiter",
    }
}

// ---------------------------------------------------------------------------
//                        PUBLIC PyO3 ENTRY POINTS
// ---------------------------------------------------------------------------

/// Dimer saddle search (Henkelman-Jónsson 1999).
///
/// # Python signature
/// ```python
/// def dimer_find_saddle(
///     positions, axis, force_callback,
///     dimer_sep=0.005, f_tol=0.05, max_steps=300,
///     max_rotor_steps=10, rotor_tol=0.01,
///     trust=0.1, trust_grow=1.1, trust_shrink=0.5, trust_max=0.5,
///     convex_max=5, frozen_mask=None,
///     history_sp=None, theta_tol=None,
///     debug=False,
/// ) -> dict:
///     # ... returns {'positions', 'axis', 'energy', 'curvature', 'status'}
/// ```
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (
    positions, axis, force_callback,
    dimer_sep=0.005, f_tol=0.05, max_steps=300,
    max_rotor_steps=10, rotor_tol=0.01,
    trust=0.1, trust_grow=1.1, trust_shrink=0.5, trust_max=0.5,
    convex_max=5, frozen_mask=None,
    history_sp=None, theta_tol=None,
    debug=false,
))]
pub fn dimer_find_saddle<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<f64>,
    axis: PyReadonlyArray2<f64>,
    force_callback: Py<PyAny>,
    dimer_sep: f64,
    f_tol: f64,
    max_steps: usize,
    max_rotor_steps: usize,
    rotor_tol: f64,
    trust: f64,
    trust_grow: f64,
    trust_shrink: f64,
    trust_max: f64,
    convex_max: usize,
    frozen_mask: Option<Vec<bool>>,
    history_sp: Option<Vec<PyReadonlyArray2<f64>>>,
    theta_tol: Option<f64>,
    debug: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let pos_vec = pyarray_to_vec3(&positions);
    let axis_vec = pyarray_to_vec3(&axis);

    let mut calc = PyCalculator {
        callback: force_callback,
        frozen: frozen_mask,
    };

    let params = DimerParams {
        dimer_sep,
        f_tol,
        max_steps,
        max_rotor_steps,
        rotor_tol,
        trust,
        trust_grow,
        trust_shrink,
        trust_max,
        convex_max,
        debug,
    };

    let history_vec: Option<Vec<Vec<[f64; 3]>>> = history_sp.map(|hs| {
        hs.iter().map(|sp| pyarray_to_vec3(sp)).collect()
    });

    let result = py.allow_threads(|| {
        find_saddle(
            &mut calc,
            &pos_vec,
            &axis_vec,
            &params,
            history_vec.as_deref(),
            theta_tol,
        )
    });

    let dict = PyDict::new(py);
    dict.set_item("positions", vec3_to_pyarray(py, &result.positions))?;
    dict.set_item("axis",      vec3_to_pyarray(py, &result.axis))?;
    dict.set_item("energy",    result.energy)?;
    dict.set_item("curvature", result.curvature)?;
    dict.set_item("status",    status_to_str(result.status))?;
    Ok(dict)
}

/// FIRE minimization (Bitzek 2006).
///
/// # Python signature
/// ```python
/// def fire_minimize(
///     positions, force_callback,
///     f_tol=0.01, max_steps=500, max_step=0.2,
///     dt_start=0.05, dt_max=0.5, dt_min=0.001,
///     alpha_start=0.1, f_inc=1.1, f_dec=0.5, f_alpha=0.99,
///     n_min=5, frozen_mask=None, debug=False,
/// ) -> dict:
///     # ... returns {'positions', 'energy', 'converged', 'n_steps', 'final_max_force'}
/// ```
#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (
    positions, force_callback,
    f_tol=0.01, max_steps=500, max_step=0.2,
    dt_start=0.05, dt_max=0.5, dt_min=0.001,
    alpha_start=0.1, f_inc=1.1, f_dec=0.5, f_alpha=0.99,
    n_min=5, frozen_mask=None, debug=false,
))]
pub fn fire_minimize<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<f64>,
    force_callback: Py<PyAny>,
    f_tol: f64,
    max_steps: usize,
    max_step: f64,
    dt_start: f64,
    dt_max: f64,
    dt_min: f64,
    alpha_start: f64,
    f_inc: f64,
    f_dec: f64,
    f_alpha: f64,
    n_min: usize,
    frozen_mask: Option<Vec<bool>>,
    debug: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let pos_vec = pyarray_to_vec3(&positions);

    let mut calc = PyCalculator {
        callback: force_callback,
        frozen: frozen_mask,
    };

    let params = FireParams {
        f_tol,
        max_steps,
        max_step,
        dt_start,
        dt_max,
        dt_min,
        alpha_start,
        f_inc,
        f_dec,
        f_alpha,
        n_min,
        debug,
    };

    let result = py.allow_threads(|| minimize(&mut calc, &pos_vec, &params));

    let dict = PyDict::new(py);
    dict.set_item("positions",       vec3_to_pyarray(py, &result.positions))?;
    dict.set_item("energy",          result.energy)?;
    dict.set_item("converged",       result.converged)?;
    dict.set_item("n_steps",         result.n_steps)?;
    dict.set_item("final_max_force", result.final_max_force)?;
    Ok(dict)
}
