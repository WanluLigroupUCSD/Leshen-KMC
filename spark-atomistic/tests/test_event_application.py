# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""`EVENT-004` event application: the verification relaxation must recover, not restate.

`EVENT-004`: "Event application MUST use the validated destination minimum, followed by
one verification relaxation. Failure to recover the destination within state tolerances
returns `EVENT_APPLICATION_FAILED`; no kinetic time is committed."

A relaxation launched at the destination's OWN coordinates cannot recover anything: the
destination was committed at `max_movable_force <= relaxation.force_tolerance`
(`RELAX-003`), so the minimiser returns at step zero and the match compares the
destination with itself. Measured before the repair: all nine (current, destination)
combinations of three mutually non-matching committed states were accepted after exactly
one calculator evaluation, the verdict was a pure function of `destination_state_id`, and
the event id `event:fabricated` -- which is in no catalog -- was accepted too.

The application geometry is now rebuilt the way `catalog.validate_candidate` built the
endpoints: the committed saddle geometry displaced along the committed unstable direction
by `saddle_search.endpoint_displacement`. Both signs are admissible because no
requirement fixes the sign of `saddle.unstable_direction` relative to the destination
(`E2-EVENT-006` canonicalises it only inside the `pair_id` hash).
"""

from __future__ import annotations

import dataclasses
import tempfile
import unittest
from unittest import mock

import fixtures
from fixtures import (StubCalculator, analytic_energy_forces, checkpoint_model,
                      double_well_chain_state, unlimited_ledger)

from spark_atomistic.canonical import deep_thaw
from spark_atomistic.catalog import match_states
from spark_atomistic.engine import ReferenceEngine
from spark_atomistic.errors import DomainFailure
from spark_atomistic.model import validate_model
from spark_atomistic.solvers import SaddleCandidate, SteepestDescentMinimizer


EXTENSION = {"calculator_command": ["/bin/true"]}
# Three minima of the analytic double well that are not related by whole-cell translation
# or same-species permutation, so they are three distinct committed states.
GEOMETRIES = {"A": (-1.0, -1.0, -1.0), "B": (-1.0, -1.0, 1.0), "C": (-1.0, 1.0, 1.0)}
# The saddle that connects A and B: atom 2 sits on the barrier top, and the unstable
# direction is that atom's x axis.
SADDLE_POSITIONS = ((-1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 0.0))
SADDLE_MODE = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0))


def build_engine():
    """An engine wired to the in-process analytic double well, with real budget scopes."""
    directory = tempfile.mkdtemp(prefix="spark-e2-apply-")
    model = deep_thaw(checkpoint_model(directory))
    model["relaxation"] = {"force_tolerance": 1e-4, "max_steps": 2000,
                           "max_evaluations": 2000}
    model["resources"]["evaluations_per_relaxation"] = 2000
    model["resources"]["total_calculator_evaluations"] = 100000
    model["resources"]["catalog_events"] = 64
    config = validate_model(model)
    engine = ReferenceEngine(config, extension=EXTENSION)
    ledger = unlimited_ledger()
    calculator = StubCalculator(ledger)
    engine.calculator = calculator
    engine.minimizer = SteepestDescentMinimizer(calculator, config["relaxation"])
    states = {name: double_well_chain_state(coordinates)
              for name, coordinates in GEOMETRIES.items()}
    for state in states.values():
        engine.catalog.states[state.state_id] = state
    engine.current_state_id = engine.initial_state_id = states["A"].state_id
    return engine, ledger, states, config


def commit_forward_event(engine, states, config):
    """The genuine A -> B event, produced by the real `Catalog.validate_candidate`."""
    candidate = SaddleCandidate(
        positions=SADDLE_POSITIONS, energy=analytic_energy_forces(SADDLE_POSITIONS)[0],
        forces=((0.0, 0.0, 0.0),) * 3, mode=SADDLE_MODE, curvature=-4.0,
        orthogonal_curvatures=(2.0,), evidence_level="DIRECTIONAL", evaluations=7,
        iterations=3, plus_endpoint=states["B"], minus_endpoint=states["A"],
        provenance={"method": "directional-dimer/1", "search_id": "search:apply",
                    "termination_reason": "full-force-curvature-orthogonal-gates"})
    outcome, event, _destination = engine.catalog.validate_candidate(
        states["A"], candidate, config["kinetics"],
        {"rng_substream_digest": "sha256:sub", "search_class": "global",
         "search_id": "search:apply", "search_index": 0})
    assert outcome.status == "OK" and event is not None
    return event


class EventApplicationTests(unittest.TestCase):

    def setUp(self) -> None:
        self.engine, self.ledger, self.states, self.config = build_engine()
        self.event = commit_forward_event(self.engine, self.states, self.config)

    def verify(self, event, destination_state_id):
        """`(verdict, calculator evaluations)` for one application attempt."""
        self.engine.catalog.events[event.event_id] = event
        before = self.ledger.calculator_reserved
        try:
            state = self.engine._verify_application(destination_state_id, event.event_id)
            return state.state_id, self.ledger.calculator_reserved - before
        except DomainFailure as exc:
            return exc.outcome.status, self.ledger.calculator_reserved - before

    def test_the_fixture_states_are_mutually_non_matching(self) -> None:
        # FIXTURE DISCRIMINATION CONTROL. If any two of these matched, every rejection
        # below could be an accident of the geometry rather than of the rule.
        kinetics = self.config["kinetics"]
        names = sorted(self.states)
        for left in names:
            for right in names:
                if left < right:
                    with self.subTest(pair=(left, right)):
                        self.assertFalse(match_states(self.states[left], self.states[right],
                                                      kinetics).equal)
        self.assertEqual(self.event.destination_state_id, self.states["B"].state_id)
        self.assertEqual(self.event.origin_state_id, self.states["A"].state_id)

    def test_a_real_event_is_accepted(self) -> None:
        # ACCEPTED BASELINE for every rejection below.
        verdict, evaluations = self.verify(self.event, self.states["B"].state_id)
        self.assertEqual(verdict, self.states["B"].state_id)
        # Before the repair this cost exactly 1 evaluation, because the minimiser was
        # launched on an already-converged geometry and returned at step 0. Real
        # reconvergence from the saddle displacement cannot be that cheap.
        self.assertGreater(evaluations, 1)

    def test_the_relaxation_starts_at_the_saddle_displacement_not_at_the_destination(self) -> None:
        # The sharpest pin on the defect: capture the geometry actually handed to the
        # minimiser. It must be `saddle +- endpoint_displacement * unstable_direction`
        # and must NOT be the destination's own coordinates.
        distance = self.config["saddle_search"]["endpoint_displacement"]
        expected = [[coordinate + distance * component
                     for coordinate, component in zip(position, direction)]
                    for position, direction in zip(SADDLE_POSITIONS, SADDLE_MODE)]
        seen: list[list[list[float]]] = []
        original = SteepestDescentMinimizer.minimize

        def recording(minimizer, request, *, object_id):
            seen.append([list(item) for item in request["positions"]])
            return original(minimizer, request, object_id=object_id)

        with mock.patch.object(SteepestDescentMinimizer, "minimize", recording):
            verdict, _evaluations = self.verify(self.event, self.states["B"].state_id)
        self.assertEqual(verdict, self.states["B"].state_id)
        self.assertEqual(len(seen), 1)  # `EVENT-004`: one verification relaxation
        for row, expected_row in zip(seen[0], expected):
            for value, expected_value in zip(row, expected_row):
                self.assertAlmostEqual(value, expected_value, places=15)
        destination_positions = [list(item) for item in self.states["B"].positions]
        self.assertNotEqual(seen[0], destination_positions)

    def test_a_fabricated_destination_is_rejected(self) -> None:
        # REJECTED counterpart, one field away from the accepted baseline: the event
        # record now names C, which the saddle displacement does not lead to.
        fabricated = dataclasses.replace(self.event, event_id="event:fabricated-destination",
                                         destination_state_id=self.states["C"].state_id)
        verdict, evaluations = self.verify(fabricated, self.states["C"].state_id)
        self.assertEqual(verdict, "EVENT_APPLICATION_FAILED")
        # Both signs were tried before refusing, so the refusal is not a sign accident.
        self.assertGreater(evaluations, 1)

    def test_a_fabricated_unstable_direction_is_rejected(self) -> None:
        # REJECTED counterpart, one field away: the mode no longer connects to B.
        fabricated = dataclasses.replace(
            self.event, event_id="event:fabricated-direction",
            unstable_direction=((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)))
        verdict, _evaluations = self.verify(fabricated, self.states["B"].state_id)
        self.assertEqual(verdict, "EVENT_APPLICATION_FAILED")

    def test_a_fabricated_saddle_geometry_is_rejected(self) -> None:
        # REJECTED counterpart, one field away: the saddle is moved into the A basin, so
        # both displaced endpoints relax back to A.
        fabricated = dataclasses.replace(
            self.event, event_id="event:fabricated-saddle",
            saddle_positions=((-1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (-1.0, 0.0, 0.0)))
        verdict, _evaluations = self.verify(fabricated, self.states["B"].state_id)
        self.assertEqual(verdict, "EVENT_APPLICATION_FAILED")

    def test_an_event_id_that_is_in_no_catalog_is_rejected(self) -> None:
        # Before the repair `event:fabricated` was accepted for all nine
        # (current, destination) combinations after one evaluation. There is no geometry
        # to rebuild from an unknown event, so it is refused before any calculator work.
        for current in sorted(self.states):
            for destination in sorted(self.states):
                with self.subTest(current=current, destination=destination):
                    self.engine.current_state_id = self.states[current].state_id
                    before = self.ledger.calculator_reserved
                    with self.assertRaises(DomainFailure) as failure:
                        self.engine._verify_application(self.states[destination].state_id,
                                                        "event:fabricated")
                    self.assertEqual(failure.exception.outcome.status,
                                     "EVENT_APPLICATION_FAILED")
                    self.assertEqual(self.ledger.calculator_reserved - before, 0)

    def test_an_event_that_does_not_leave_the_current_state_is_rejected(self) -> None:
        # `E2-KMC-001` only ever offers events whose origin is the current state; an
        # event that does not is an inconsistency, not an application.
        self.engine.current_state_id = self.states["C"].state_id
        verdict, evaluations = self.verify(self.event, self.states["B"].state_id)
        self.assertEqual(verdict, "EVENT_APPLICATION_FAILED")
        self.assertEqual(evaluations, 0)


if __name__ == "__main__":
    unittest.main()
