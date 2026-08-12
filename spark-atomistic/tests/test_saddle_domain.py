# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Characterised domain of `DirectionalDimerSearcher` on the analytic double well.

`V = sum_i (x_i^2-1)^2 + y_i^2 + z_i^2` (eV, Angstrom): minima at x = +-1 with E = 0,
first-order saddle at x = 0 with E = 1.0 eV, forward and reverse barriers both exactly
1.000 eV.

These tests record a PARAMETER-DOMAIN LIMITATION, not a defect. `SADDLE-002` lets the
interface carry any first-order-saddle solver, `SADDLE-005` defines the acceptance
evidence, and `SADDLE-006`/`SADDLE-007` plus the `E2-STATUS-002` severity of
`SADDLE_NOT_FOUND` ("candidate reject") make a non-converging attempt a legitimate,
non-fatal outcome. No requirement obliges a minimum-mode follower to converge from a
launch point where the Hessian has no negative eigenvalue, so the algorithm is left
unchanged and the boundary is measured instead.
"""

from __future__ import annotations

import math
import unittest

import fixtures
from fixtures import (StubCalculator, analytic_hessian_eigenvalues, double_well_state,
                      unlimited_ledger)

from spark_atomistic.errors import DomainFailure
from spark_atomistic.rng import derive_saddle_stream
from spark_atomistic.solvers import DirectionalDimerSearcher, SteepestDescentMinimizer


BASE_SADDLE_CONFIG = {"curvature_tolerance": 1e-6, "endpoint_displacement": 0.01,
                      "force_tolerance": 0.05, "max_iterations": 2000,
                      "method": "directional-dimer", "orthogonal_directions": 1}


def search(x0: float, *, force_tolerance: float = 0.05, endpoint_displacement: float = 0.01,
           max_iterations: int = 2000, search_index: int = 0):
    """Drive the solver in process. The subprocess transport is out of scope here."""
    ledger = unlimited_ledger()
    calculator = StubCalculator(ledger)
    relaxation = {"force_tolerance": force_tolerance, "max_steps": 200_000,
                  "max_evaluations": 200_000}
    saddle = dict(BASE_SADDLE_CONFIG, force_tolerance=force_tolerance,
                  endpoint_displacement=endpoint_displacement, max_iterations=max_iterations)
    minimizer = SteepestDescentMinimizer(calculator, relaxation)
    searcher = DirectionalDimerSearcher(calculator, minimizer, saddle)
    origin = double_well_state(x0)
    stream = derive_saddle_stream(0, origin.state_id, "global", search_index)
    try:
        candidate = searcher.search(origin, stream, search_id="search:domain")
    except DomainFailure as exc:
        return {"status": exc.outcome.status, "evaluations": ledger.calculator_reserved,
                "candidate": None, "calculator": calculator,
                "details": exc.outcome.context["details"]}
    return {"status": "OK", "evaluations": ledger.calculator_reserved,
            "candidate": candidate, "calculator": calculator, "details": {}}


class AnalyticGeometryTests(unittest.TestCase):
    def test_launch_point_curvature_classification(self) -> None:
        # d2V/dx2 = 12x^2 - 4 ; d2V/dy2 = d2V/dz2 = 2.
        self.assertEqual(analytic_hessian_eigenvalues(0.0), (-4.0, 2.0, 2.0))
        self.assertEqual(analytic_hessian_eigenvalues(1.0), (8.0, 2.0, 2.0))
        # Negative curvature along the reaction coordinate exists only for |x| < 1/sqrt(3).
        self.assertAlmostEqual(1 / math.sqrt(3), 0.5773502691896258, places=15)
        self.assertLess(analytic_hessian_eigenvalues(0.577)[0], 0.0)
        self.assertGreater(analytic_hessian_eigenvalues(0.578)[0], 0.0)
        # The reaction coordinate stops being the LOWEST mode above |x| = 1/sqrt(2).
        self.assertAlmostEqual(1 / math.sqrt(2), 0.7071067811865475, places=15)
        self.assertLess(analytic_hessian_eigenvalues(0.70)[0], 2.0)
        self.assertGreater(analytic_hessian_eigenvalues(0.72)[0], 2.0)
        # The reported failing launch point is convex in every direction.
        eigenvalues = analytic_hessian_eigenvalues(0.9)
        self.assertAlmostEqual(eigenvalues[0], 5.72, places=12)
        self.assertTrue(all(value > 0.0 for value in eigenvalues))
        self.assertEqual(min(eigenvalues), 2.0)  # the lowest mode is transverse


class ConvergentDomainTests(unittest.TestCase):
    """ACCEPTED baseline: inside its domain the solver reproduces the analytic answer."""

    def test_dimer_finds_the_analytic_saddle_and_unit_barriers(self) -> None:
        result = search(0.30, force_tolerance=1e-4, endpoint_displacement=0.05)
        self.assertEqual(result["status"], "OK")
        candidate = result["candidate"]
        self.assertEqual(candidate.evidence_level, "DIRECTIONAL")  # SADDLE-005
        self.assertLess(candidate.curvature, -1e-6)
        self.assertTrue(all(value >= -1e-6 for value in candidate.orthogonal_curvatures))
        # Saddle at x = 0 with E = 1.0 eV.
        self.assertAlmostEqual(candidate.positions[0][0], 0.0, places=4)
        self.assertAlmostEqual(candidate.energy, 1.0, places=8)
        self.assertAlmostEqual(candidate.curvature, -4.0, places=3)
        self.assertAlmostEqual(candidate.orthogonal_curvatures[0], 2.0, places=9)
        # Endpoints relax to the two minima; both barriers are 1.000 eV.
        endpoints = sorted((candidate.plus_endpoint.positions[0][0],
                            candidate.minus_endpoint.positions[0][0]))
        self.assertAlmostEqual(abs(endpoints[0]), 1.0, places=4)
        self.assertAlmostEqual(abs(endpoints[1]), 1.0, places=4)
        forward = candidate.energy - candidate.plus_endpoint.energy
        reverse = candidate.energy - candidate.minus_endpoint.energy
        self.assertAlmostEqual(forward, 1.0, places=8)
        self.assertAlmostEqual(reverse, 1.0, places=8)

    def test_accuracy_is_bounded_by_the_configured_force_tolerance(self) -> None:
        # Near the saddle V ~= 1 - 2x^2 and |F| ~= 4|x|, so the force-tolerance ball has
        # radius f/4 and admits an energy error of at most 2*(f/4)^2 = f^2/8. The
        # measured error tracks that bound over five decades, which is why the loose
        # tolerance -- not the algorithm -- is what keeps the barrier away from 1.000 eV.
        measured = {}
        for tolerance in (5e-2, 1e-2, 1e-3, 1e-4, 1e-5):
            result = search(0.30, force_tolerance=tolerance, endpoint_displacement=0.05)
            self.assertEqual(result["status"], "OK")
            error = abs(result["candidate"].energy - 1.0)
            measured[tolerance] = error
            self.assertLess(error, 4.0 * tolerance * tolerance / 8.0)
        self.assertGreater(measured[5e-2], 1e-4)     # 2.9e-4 eV at the corpus tolerance
        self.assertLess(measured[1e-5], 1e-9)        # 2.4e-11 eV once tightened
        self.assertGreater(measured[5e-2] / measured[1e-5], 1e6)

    def test_loose_tolerance_can_produce_a_negative_forward_barrier(self) -> None:
        # `E2-RATE-001` returns `RATE_INVALID` when a raw barrier is below
        # -barrier_tolerance. With relaxation.force_tolerance = 0.05 eV/A the entire
        # region |x| <= 0.0125 A is force-converged, so an endpoint can be committed on
        # the barrier top and the raw forward barrier goes negative. The shared
        # e2_minimal_model.json pairs that force tolerance with barrier_tolerance =
        # 1e-10 eV, which is 3.1e6 times smaller than the admitted error f^2/8.
        result = search(0.50, force_tolerance=0.05, endpoint_displacement=0.01)
        self.assertEqual(result["status"], "OK")
        candidate = result["candidate"]
        forward = candidate.energy - candidate.plus_endpoint.energy
        self.assertLess(forward, 0.0)
        self.assertGreater(forward, -0.05 * 0.05 / 8.0)
        self.assertLess(forward, -1e-10)  # below the corpus barrier_tolerance


class DivergentDomainTests(unittest.TestCase):
    """REJECTED counterpart: outside the domain the solver reports a candidate reject."""

    def test_launch_without_a_negative_eigenvalue_exhausts_the_iteration_budget(self) -> None:
        result = search(0.9, max_iterations=2000)
        self.assertEqual(result["status"], "SADDLE_NOT_FOUND")
        self.assertEqual(result["evaluations"], 6000)  # 3 evaluations per dimer iteration
        self.assertEqual(result["details"]["iterations"], 2000)
        self.assertEqual(result["details"]["termination_reason"],
                         "maximum dimer iterations exhausted")  # SADDLE-007

    def test_the_failure_is_a_transverse_runaway_not_slow_convergence(self) -> None:
        # The lowest mode at x = 0.9 is transverse (eigenvalue 2, versus 5.72 along the
        # reaction coordinate), so minimum-mode following locks onto y/z. The reflected
        # force then climbs that transverse mode while the centre slides into the
        # minimum at x = +1. Growth per iteration is exactly 1 + dt*lambda_perp with the
        # hard-coded translation step dt = 0.01 and lambda_perp = 2 eV/A^2.
        result = search(0.9, max_iterations=2000)
        centres = result["calculator"].centres[::3]
        self.assertGreater(len(centres), 1000)
        self.assertAlmostEqual(centres[200][0], 1.0, places=6)      # slid into the minimum
        transverse = [math.hypot(centre[1], centre[2]) for centre in centres]
        self.assertLess(transverse[0], 1e-2)
        self.assertGreater(transverse[1800], 1e13)                  # 15 orders of growth
        ratios = [transverse[k + 1] / transverse[k] for k in range(600, 1000)]
        self.assertAlmostEqual(sum(ratios) / len(ratios), 1.02, places=9)

    def test_more_iterations_convert_the_reject_into_a_calculator_failure(self) -> None:
        # The runaway is not benign. Past ~2010 dimer iterations the transverse
        # coordinate overflows binary64, the analytic calculator exits nonzero, and
        # `CALC-006` maps that to `CALCULATOR_FAILURE`. Unlike `SADDLE_NOT_FOUND`, that
        # status is not in the engine's candidate-reject continue set, so a longer
        # budget turns one rejected candidate into a terminated run (exit 69).
        result = search(0.9, max_iterations=20000)
        self.assertEqual(result["status"], "CALCULATOR_FAILURE")
        self.assertEqual(result["evaluations"], 6031)
        self.assertEqual(result["calculator"].overflow_at, 6031)

    def test_the_failure_is_reproducible_across_independent_substreams(self) -> None:
        # `DISC-002`/`E2-RNG-004`: each search index is an independent substream. The
        # outcome is a property of the launch point, not of one unlucky random mode.
        for index in range(6):
            with self.subTest(search_index=index):
                result = search(0.9, max_iterations=200, search_index=index)
                self.assertEqual(result["status"], "SADDLE_NOT_FOUND")

    def test_domain_boundary_is_probabilistic_not_a_clean_threshold(self) -> None:
        # Honest characterisation, measured over eight independent search substreams per
        # launch point rather than one. A single scan of x0 produces an interleaved
        # pattern that looks like a sharp threshold only because each x0 draws a
        # different initial mode; the reproducible quantity is the success RATE.
        #
        #   x0        0.0 0.1 0.2 0.3 0.4 | 0.5 0.55 0.6 0.65 0.70 0.75 | 0.8 0.9 1.0
        #   OK / 8     8   8   8   8   8  |  7   8    6   6    4    1   |  0   0   0
        #
        # 1/sqrt(3) = 0.5774 (curvature along x turns positive) and 1/sqrt(2) = 0.7071
        # (x stops being the LOWEST mode) explain the decline; neither is a threshold.
        def rate(x0: float) -> int:
            return sum(1 for index in range(8)
                       if search(x0, max_iterations=400, search_index=index)["status"] == "OK")

        for x0 in (0.0, 0.1, 0.2, 0.3, 0.4):
            with self.subTest(region="every substream converges", x0=x0):
                self.assertEqual(rate(x0), 8)
        for x0 in (0.8, 0.9, 1.0):
            with self.subTest(region="no substream converges", x0=x0):
                self.assertEqual(rate(x0), 0)
        for x0 in (0.6, 0.65, 0.7, 0.75):
            with self.subTest(region="mixed", x0=x0):
                self.assertTrue(0 < rate(x0) < 8)
        # The decline is monotone across the mode-crossing at 1/sqrt(2).
        self.assertGreater(rate(0.6), rate(0.75))

    def test_outcome_flips_on_a_one_ulp_change_of_the_launch_coordinate(self) -> None:
        # The initial dimer mode is drawn from a substream keyed by `state_id`, and
        # `state_id` hashes the binary64 energy. Two adjacent doubles for the same
        # launch coordinate therefore draw unrelated modes and can land on opposite
        # sides of the outcome. This is why no threshold in x0 can be quoted.
        low = search(0.59, max_iterations=400)["status"]
        high = search(0.5900000000000001, max_iterations=400)["status"]
        self.assertNotEqual(low, high)
        self.assertIn("OK", (low, high))
        self.assertIn("SADDLE_NOT_FOUND", (low, high))


if __name__ == "__main__":
    unittest.main()
