# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""What the `E2-EVENT-003` detailed-balance residual can and cannot demonstrate.

`E2-RATE-001` fixes both legs of a directed pair to ONE saddle energy ("Raw barriers are
exactly `b_f=E_s-E_i` and `b_r=E_s-E_j`"; `E2-EVENT-001` calls them "Raw same-saddle
differences") and `E2-RATE-002` fixes the rates to those barriers under the
`COMMON_PREFACTOR` equations. Therefore
`log_forward - log_reverse = beta*(E_i - E_j)` identically, and the residual
`(log_forward - log_reverse) + beta*(E_j - E_i)` is algebraically zero for every input.
No independently obtained reverse barrier exists to compare against: the freshly relaxed
reverse endpoint energy is not admissible because `E2-EVENT-001` requires
`destination_state_id` to be a COMMITTED state and `E2-RATE-001` pins `b_r` to that
committed state's energy, which is the quantity already used. Option (a) of the repair
is therefore unreachable without inventing a second physical quantity, so option (b) was
taken: the commit-time gate is deleted. The `detailed_balance_residual` FIELD stays,
because `E2-EVENT-003` fixes the `rate_model` key set exactly and deleting it would break
both the schema and the checkpoint contract.

The verification that does exist is at restore, where the two sides are independent:
`E2-CKPT-007`(5) makes `checkpoint._restore_event` recompute the residual from the stored
barriers and energies and refuse a stored value that disagrees. That path is exercised
here with an accepted baseline and a single-field tamper.
"""

from __future__ import annotations

import math
import random
import tempfile
import unittest

import fixtures
from fixtures import (CELL, CHECKPOINT_MODEL_DIGEST, build_grown_catalog_checkpoint,
                      clone_payload)

from spark_atomistic.canonical import digest
from spark_atomistic.catalog import Catalog
from spark_atomistic.checkpoint import read_checkpoint, validate_checkpoint_payload
from spark_atomistic.errors import DomainFailure
from spark_atomistic.model import state_from_relaxation
from spark_atomistic.solvers import SaddleCandidate


BASE_KINETICS = {
    "temperature": 300.0, "prefactor": 1e13, "barrier_tolerance": 1e-10,
    "detailed_balance_tolerance": 1e-8, "log_rate_cutoff": -700.0,
    "state_rms_tolerance": 1e-3, "state_max_tolerance": 5e-3,
    "state_energy_tolerance_per_atom": 1e-6, "saddle_rms_tolerance": 1e-3,
    "saddle_max_tolerance": 5e-3, "saddle_energy_tolerance": 1e-5,
}
BOLTZMANN_EV_PER_K = 8.617333262145e-5


def _state(positions, energy: float):
    request = {"schema": "spark-atomistic-model/1", "atom_ids": ["a0", "a1"],
               "species": ["H", "H"], "positions": positions, "cell": CELL,
               "pbc": [False, False, False], "movable": [True, True],
               "constraints": {"kind": "fixed-mask"}, "charge": 0.0, "spin": 0.0,
               "calculator_model_digest": CHECKPOINT_MODEL_DIGEST}
    return state_from_relaxation(request, energy, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                                 0.05, {"calculator_evaluations": 1, "steps": 0,
                                        "calculator_identity": "fixture@1",
                                        "minimizer_identity": "fixture-minimizer/1",
                                        "termination_reason": "force_tolerance"})


def commit(origin_energy: float, destination_energy: float, saddle_energy: float,
           *, keep: bool = False, **overrides):
    """Drive one candidate through the real `Catalog.validate_candidate` rate block."""
    kinetics = dict(BASE_KINETICS, **overrides)
    origin = _state([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], origin_energy)
    destination = _state([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], destination_energy)
    catalog = Catalog("sha256:model", "sha256:config", "sha256:tolerance", 64,
                      "sha256:identity")
    catalog.add_state(origin, kinetics)
    candidate = SaddleCandidate(
        positions=((0.0, 0.0, 0.0), (1.5, 0.0, 0.0)), energy=saddle_energy,
        forces=((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
        mode=((1.0, 0.0, 0.0), (0.0, 0.0, 0.0)), curvature=-4.0,
        orthogonal_curvatures=(2.0,), evidence_level="DIRECTIONAL", evaluations=7,
        iterations=3, plus_endpoint=origin, minus_endpoint=destination,
        provenance={"method": "directional-dimer/1", "search_id": "search:x",
                    "termination_reason": "full-force-curvature-orthogonal-gates"})
    result = catalog.validate_candidate(origin, candidate, kinetics,
                                        {"rng_substream_digest": "sha256:sub",
                                         "search_class": "global", "search_id": "search:x",
                                         "search_index": 0})
    return (*result, catalog) if keep else result


def sweep(samples: int, **overrides):
    """Physically valid candidates over the domain the defect report measured."""
    rng = random.Random(20260811)
    statuses: dict[str, int] = {}
    residuals: list[float] = []
    for _ in range(samples):
        temperature = rng.choice((100.0, 300.0, 700.0, 1200.0))
        prefactor = rng.choice((1e12, 1e13, 6.2e12))
        origin_energy = rng.uniform(-500.0, 500.0)
        destination_energy = origin_energy + rng.uniform(-3.0, 3.0)
        saddle_energy = max(origin_energy, destination_energy) + rng.uniform(0.0, 3.0)
        try:
            outcome, event, _destination = commit(
                origin_energy, destination_energy, saddle_energy,
                temperature=temperature, prefactor=prefactor, **overrides)
        except DomainFailure as exc:
            statuses[exc.outcome.status] = statuses.get(exc.outcome.status, 0) + 1
            continue
        statuses[outcome.status] = statuses.get(outcome.status, 0) + 1
        if event is not None:
            residuals.append(event.detailed_balance_residual)
    return statuses, residuals


class CommitTimeResidualTests(unittest.TestCase):
    """The residual written into every event record carries no information."""

    SAMPLES = 750

    def test_the_sweep_actually_reaches_the_rate_model_block(self) -> None:
        # ANTI-VACUITY CONTROL. Every assertion below is about what the rate block does
        # NOT reject; that is worthless unless this harness provably executes that block.
        # Two controls, both differing from the sweep only in the energies:
        # (1) a saddle below the origin by more than `barrier_tolerance` must return
        #     `RATE_INVALID` (`E2-RATE-001`), which is raised inside the block;
        # (2) a committed event must carry the exact `E2-RATE-002` COMMON_PREFACTOR
        #     numbers, which are computed AFTER the deleted gate's position, so the
        #     sweep runs through the whole block and not merely into it.
        with self.assertRaises(DomainFailure) as failure:
            commit(0.0, -1.0, -0.5)
        self.assertEqual(failure.exception.outcome.status, "RATE_INVALID")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-RATE-001")

        outcome, event, _destination = commit(-10.0, -10.5, -9.0)
        self.assertEqual(outcome.status, "OK")
        beta = 1.0 / (BOLTZMANN_EV_PER_K * BASE_KINETICS["temperature"])
        self.assertEqual(event.barrier, 1.0)
        self.assertEqual(event.reverse_barrier, 1.5)
        self.assertEqual(event.log_rate,
                         math.log(BASE_KINETICS["prefactor"]) - beta * 1.0)
        self.assertEqual(event.reverse_log_rate,
                         math.log(BASE_KINETICS["prefactor"]) - beta * 1.5)
        self.assertEqual(event.record()["rate_model"]["model"], "COMMON_PREFACTOR")

    def test_residual_is_an_algebraic_identity_over_the_physical_domain(self) -> None:
        # Measured: 750 randomised physically valid candidates here, 20000 in the report; T in
        # {100,300,700,1200} K, nu in {1e12,1e13,6.2e12} s^-1, |dE| <= 3 eV, |E| <= 500 eV.
        statuses, residuals = sweep(self.SAMPLES)
        self.assertEqual(statuses, {"OK": self.SAMPLES})
        self.assertEqual(len(residuals), self.SAMPLES)
        worst = max(abs(item) for item in residuals)
        # 1.1368683772161603e-13 measured on CPython 3.12.3; that is binary64 rounding of
        # ~5.8e4-magnitude intermediates, five orders below the configured 1e-8 tolerance.
        self.assertLess(worst, 1e-12)
        self.assertLess(max(residuals) - min(residuals), 1e-12)
        # The column is constant: knowing the residual tells a reader nothing about the
        # energies, the temperature, or the prefactor that produced it.
        self.assertGreater(statuses["OK"], 0)

    def test_a_residual_tolerance_gate_could_only_reject_valid_physics(self) -> None:
        # REGRESSION CATCHER for "the gate returns as an identity". The deleted gate is
        # `abs(residual) > kinetics["detailed_balance_tolerance"]`. Because the residual
        # is rounding noise around zero, that comparison cannot reject a detailed-balance
        # violation at any tolerance; the only thing it can do is reject valid physics
        # once the tolerance drops below the rounding floor. Measured with the gate still
        # present: 11829 of 20000 identical candidates were rejected with
        # `DETAILED_BALANCE_VIOLATION` at `detailed_balance_tolerance = 1e-18`.
        # ACCEPTED baseline first: the corpus tolerance.
        loose, _residuals = sweep(self.SAMPLES, detailed_balance_tolerance=1e-8)
        self.assertEqual(loose, {"OK": self.SAMPLES})
        # Same candidates, tolerance below the rounding floor, still every one accepted.
        tight, _tight_residuals = sweep(self.SAMPLES, detailed_balance_tolerance=1e-18)
        self.assertEqual(tight, {"OK": self.SAMPLES})
        self.assertNotIn("DETAILED_BALANCE_VIOLATION", tight)

    def test_the_reciprocal_record_negates_the_same_identity(self) -> None:
        # `E2-EVENT-003`: "For the reciprocal record, forward/reverse logs swap and the
        # residual changes sign." The field is mandatory and is still emitted; only the
        # gate is gone.
        outcome, event, _destination, catalog = commit(-10.0, -10.5, -9.0, keep=True)
        self.assertEqual(outcome.status, "OK")
        self.assertEqual(set(event.record()["rate_model"]),
                         {"common_prefactor_per_s", "detailed_balance_residual",
                          "log_forward_rate_per_s", "log_reverse_rate_per_s",
                          "model", "temperature_k"})
        reverse = catalog.events[event.reverse_event_id]
        self.assertEqual(reverse.log_rate, event.reverse_log_rate)
        self.assertEqual(reverse.reverse_log_rate, event.log_rate)
        self.assertEqual(reverse.detailed_balance_residual, -event.detailed_balance_residual)
        self.assertLess(abs(event.detailed_balance_residual), 1e-12)


class RestoreTimeDetailedBalanceTests(unittest.TestCase):
    """`E2-CKPT-007`(5): at restore the stored residual is an independent input."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.directory = tempfile.mkdtemp(prefix="spark-e2-balance-")
        cls.built = build_grown_catalog_checkpoint(cls.directory)
        cls.payload = read_checkpoint(cls.built.checkpoint_path)
        cls.arguments = cls.built.restore_arguments()

    def restore(self, payload):
        return validate_checkpoint_payload(payload, **self.arguments)

    @staticmethod
    def reseal(payload):
        """Re-derive the catalog digest so `E2-CKPT-006` cannot answer for `E2-EVENT-003`.

        Without this the catalog digest fires first and every tamper below would pass for
        the wrong reason, never reaching the detailed-balance recomputation it claims to
        exercise.
        """
        catalog = payload["catalog"]
        catalog["digest"] = digest({key: catalog[key]
                                    for key in ("events", "multiplicity", "schema", "states")})
        return payload

    def test_untampered_checkpoint_restores(self) -> None:
        # ACCEPTED baseline for the rejections below.
        restored = self.restore(clone_payload(self.payload))
        self.assertEqual(len(restored["catalog"].events), 6)
        for event in restored["catalog"].events.values():
            self.assertLess(abs(event.detailed_balance_residual), 1e-12)

    def test_reseal_alone_is_not_what_makes_a_tamper_pass_or_fail(self) -> None:
        # ANTI-VACUITY CONTROL for `reseal`: resealing an untouched payload must still
        # restore, and an unresealed tamper must fail on the DIGEST, not on the residual.
        self.restore(self.reseal(clone_payload(self.payload)))
        payload = clone_payload(self.payload)
        event_id = sorted(payload["catalog"]["events"])[0]
        payload["catalog"]["events"][event_id]["rate_model"]["detailed_balance_residual"] = 1e-12
        with self.assertRaises(DomainFailure) as failure:
            self.restore(payload)
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-CKPT-006")

    def test_a_tampered_stored_residual_is_refused(self) -> None:
        # REJECTED counterpart, differing from the accepted baseline in ONE field and
        # resealed so the catalog digest cannot answer for it. The tampered value is
        # 1e-12, INSIDE the `detailed_balance_tolerance` of 1e-8, so the vacuous
        # `abs(expected_residual) > tolerance` term cannot be what catches it: the
        # catching term is the recomputation-equality, the only part of the restore check
        # with any power.
        payload = clone_payload(self.payload)
        event_id = sorted(payload["catalog"]["events"])[0]
        rate = payload["catalog"]["events"][event_id]["rate_model"]
        self.assertLess(abs(rate["detailed_balance_residual"]), 1e-12)
        rate["detailed_balance_residual"] = 1e-12
        self.assertLess(abs(rate["detailed_balance_residual"]),
                        self.built.config["kinetics"]["detailed_balance_tolerance"])
        with self.assertRaises(DomainFailure) as failure:
            self.restore(self.reseal(payload))
        self.assertEqual(failure.exception.outcome.status, "CHECKPOINT_CORRUPT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-EVENT-003")

    def test_a_reciprocal_residual_that_is_not_the_negation_is_refused(self) -> None:
        # `E2-EVENT-003` reciprocity. Both records are moved to the SAME value, so each
        # one still fails its own recomputation; the point of the pair is that a
        # sign-consistent tamper is not available without breaking the recomputation too.
        payload = clone_payload(self.payload)
        events = payload["catalog"]["events"]
        event_id = sorted(events)[0]
        reverse_id = events[event_id]["reverse_event_id"]
        self.assertNotEqual(reverse_id, event_id)
        events[event_id]["rate_model"]["detailed_balance_residual"] = 1e-12
        events[reverse_id]["rate_model"]["detailed_balance_residual"] = 1e-12
        with self.assertRaises(DomainFailure) as failure:
            self.restore(self.reseal(payload))
        self.assertEqual(failure.exception.outcome.status, "CHECKPOINT_CORRUPT")


if __name__ == "__main__":
    unittest.main()
