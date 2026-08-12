# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""`SADDLE-005` saddle-order evidence: the orthogonal curvature must be a bound.

`SADDLE-005`: "Without a Hessian, negative curvature along the reported mode plus
nonnegative sampled orthogonal curvatures is required and the evidence level MUST be
`DIRECTIONAL`". A Rayleigh quotient `d^T H d` along ONE random `d` perpendicular to the
mode is a positively weighted average of the whole orthogonal spectrum, so it is
dominated by the stiff modes and cannot distinguish an index-1 saddle from an index-7
one. `E2-EVENT-002` then records `evidence_level = "DIRECTIONAL"` on that evidence, which
is a claim the measurement does not support.

`DirectionalDimerSearcher._orthogonal_curvatures` therefore reports, per configured
restart, the MINIMUM of the Rayleigh quotient found by a second dimer-style rotation
confined to `d` perpendicular to `mode`. Every reported number is still a directional
curvature that was measured by central differences, but it is an upper bound on the
lowest curvature of the complement instead of one draw from it.

The exact quadratic fixture is used because the answer is then known in closed form: the
gate must reject exactly the spectra that carry a second negative eigenvalue.
"""

from __future__ import annotations

import math
import unittest

import fixtures
from fixtures import quadratic_saddle_probe, unlimited_ledger

from spark_atomistic.rng import derive_saddle_stream
from spark_atomistic.solvers import (DirectionalDimerSearcher, SteepestDescentMinimizer,
                                     _flat_normalize, _flat_project_out)


ATOMS = 5           # 15 degrees of freedom, so the complement of the mode has 14
CURVATURE_TOLERANCE = 1e-6
RELAXATION = {"force_tolerance": 1e-4, "max_steps": 1000, "max_evaluations": 10000}
TRIALS = 50

# One soft mode in the complement (index 2 overall): must be rejected.
INDEX_TWO = [-1.0, -0.5] + [5.0] * 13
# Six soft modes in the complement (index 7 overall): must be rejected.
INDEX_SEVEN = [-1.0] + [-0.5] * 6 + [5.0] * 8
# No soft mode in the complement (index 1): must be accepted.
INDEX_ONE = [-1.0] + [5.0] * 14
# Every direction downhill (a local maximum): must be rejected.
LOCAL_MAXIMUM = [-1.0] + [-0.5] * 14


def searcher_for(spectrum, directions: int):
    request, positions, mode, factory = quadratic_saddle_probe(spectrum, atoms=ATOMS)
    ledger = unlimited_ledger()
    calculator = factory(ledger)
    config = {"curvature_tolerance": CURVATURE_TOLERANCE, "endpoint_displacement": 0.01,
              "force_tolerance": 0.05, "max_iterations": 2000,
              "method": "directional-dimer", "orthogonal_directions": directions}
    searcher = DirectionalDimerSearcher(
        calculator, SteepestDescentMinimizer(calculator, RELAXATION), config)
    return searcher, request, positions, mode, ledger


def bounds(spectrum, *, directions: int = 1, trial: int = 0):
    """`(reported curvatures, calculator evaluations)` for one independent restart set."""
    searcher, request, positions, mode, ledger = searcher_for(spectrum, directions)
    stream = derive_saddle_stream(0, "state:sha256:%d" % trial, "global", trial)
    values = searcher._orthogonal_curvatures(request, positions, mode, [True] * ATOMS,
                                             stream, "search:order")
    return values, ledger.calculator_reserved


def start_sample(spectrum, *, trial: int = 0):
    """The number the deleted implementation reported: one random-direction sample."""
    searcher, request, positions, mode, _ledger = searcher_for(spectrum, 1)
    stream = derive_saddle_stream(0, "state:sha256:%d" % trial, "global", trial)
    raw = [2.0 * stream.uniform() - 1.0 for _atom in range(ATOMS) for _axis in range(3)]
    direction = _flat_normalize(_flat_project_out(raw, mode))
    curvature, _product = searcher._curvature_along(request, positions, direction,
                                                    [True] * ATOMS, "search:order", 1e-3)
    return curvature


def accepted(values) -> bool:
    """The acceptance predicate of `DirectionalDimerSearcher._search`."""
    return all(item >= -CURVATURE_TOLERANCE for item in values)


def acceptance_rate(spectrum, *, directions: int = 1, trials: int = TRIALS):
    hits = 0
    lowest = math.inf
    for trial in range(trials):
        values, _evaluations = bounds(spectrum, directions=directions, trial=trial)
        lowest = min(lowest, min(values))
        hits += int(accepted(values))
    return hits, lowest


class SaddleOrderGateTests(unittest.TestCase):
    """One accepted baseline and three rejections, all on the same 15-dof fixture."""

    def test_index_one_saddle_is_accepted(self) -> None:
        # ACCEPTED BASELINE. It differs from the index-2 case below in exactly one
        # eigenvalue, so the rejections cannot be blamed on the fixture or the harness.
        hits, lowest = acceptance_rate(INDEX_ONE)
        self.assertEqual(hits, TRIALS)
        self.assertAlmostEqual(lowest, 5.0, places=9)
        # The complement is isotropic, so the rotation stops at once: 2 evaluations.
        _values, evaluations = bounds(INDEX_ONE)
        self.assertEqual(evaluations, 2)

    def test_index_two_saddle_is_rejected(self) -> None:
        # Before the repair this construction was accepted 500/500 at
        # `orthogonal_directions` 1, 3 and 10, and the most negative curvature ever
        # sampled was +2.843395: the true lambda_2 = -0.5 was never approached.
        for directions in (1, 3):
            with self.subTest(orthogonal_directions=directions):
                hits, lowest = acceptance_rate(INDEX_TWO, directions=directions)
                self.assertEqual(hits, 0)
                self.assertAlmostEqual(lowest, -0.5, places=6)

    def test_index_seven_saddle_is_rejected(self) -> None:
        for directions in (1, 3):
            with self.subTest(orthogonal_directions=directions):
                hits, lowest = acceptance_rate(INDEX_SEVEN, directions=directions)
                self.assertEqual(hits, 0)
                self.assertAlmostEqual(lowest, -0.5, places=6)

    def test_a_local_maximum_is_still_rejected(self) -> None:
        # POSITIVE CONTROL carried over from the pre-repair measurement: the old gate
        # could detect this one case, and the repair must not lose it.
        hits, lowest = acceptance_rate(LOCAL_MAXIMUM)
        self.assertEqual(hits, 0)
        self.assertAlmostEqual(lowest, -0.5, places=9)

    def test_a_single_random_direction_cannot_see_the_soft_mode(self) -> None:
        # This is the measurement that makes the repair necessary, and it is what fails
        # if `_orthogonal_curvatures` ever regresses to reporting its random start.
        # For the index-2 spectrum a negative sample needs the random direction to place
        # c1^2 > 5/5.5 = 0.909 of its norm on the single soft mode; 200000 draws produced
        # none.
        samples = [start_sample(INDEX_TWO, trial=trial) for trial in range(TRIALS)]
        self.assertTrue(all(item > 1.5 for item in samples))
        self.assertTrue(all(item >= -CURVATURE_TOLERANCE for item in samples))
        for trial in range(TRIALS):
            values, _evaluations = bounds(INDEX_TWO, trial=trial)
            with self.subTest(trial=trial):
                # Same start, same fixture: the minimisation only ever lowers the number.
                self.assertLessEqual(min(values), samples[trial])
                self.assertLess(min(values), -CURVATURE_TOLERANCE)

    def test_the_reported_value_never_claims_more_than_it_measured(self) -> None:
        # A bound must not invent curvature it did not find: for every fixture the
        # reported number stays at or above the true lowest eigenvalue of the complement,
        # and at or below the sample it started from.
        for label, spectrum in (("index-1", INDEX_ONE), ("index-2", INDEX_TWO),
                                ("index-7", INDEX_SEVEN), ("local-max", LOCAL_MAXIMUM)):
            true_minimum = min(spectrum[1:])
            for trial in range(5):
                values, _evaluations = bounds(spectrum, directions=2, trial=trial)
                with self.subTest(spectrum=label, trial=trial):
                    self.assertEqual(len(values), 2)  # one per configured restart
                    self.assertGreaterEqual(min(values), true_minimum - 1e-6)
                    self.assertLessEqual(min(values), start_sample(spectrum, trial=trial) + 1e-9)

    def test_the_rotation_budget_is_bounded_and_only_spent_while_it_pays(self) -> None:
        # `RES-001`: the new evidence must stay inside a declarable budget. A cycle costs
        # four evaluations and a restart starts with two, so a restart can never exceed
        # 2 + 4 * 40 = 162 evaluations; an isotropic complement stops at 2.
        for label, spectrum, expected in (("isotropic", INDEX_ONE, 2),
                                          ("local-max", LOCAL_MAXIMUM, 2)):
            with self.subTest(spectrum=label):
                _values, evaluations = bounds(spectrum)
                self.assertEqual(evaluations, expected)
        _values, evaluations = bounds(INDEX_TWO)
        self.assertGreater(evaluations, 2)
        self.assertLessEqual(evaluations, 2 + 4 * 40)


if __name__ == "__main__":
    unittest.main()
