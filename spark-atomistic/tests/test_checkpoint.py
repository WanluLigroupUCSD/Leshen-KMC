# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Errata-2 section 10/11 checkpoint closeout: E2-CKPT-* and the E2-KMC-005 replay.

Every rejection fixture below is paired with the accepted baseline it was derived
from, because a fixture that is rejected for the wrong reason proves nothing.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
import tempfile
import unittest

import fixtures
from fixtures import build_grown_catalog_checkpoint, clone_payload

from spark_atomistic.canonical import digest
from spark_atomistic.checkpoint import read_checkpoint, validate_checkpoint_payload
from spark_atomistic.errors import DomainFailure
from spark_atomistic.kinetics import build_rate_table


class HistoricalRateSnapshotReplayTests(unittest.TestCase):
    """`E2-KMC-005`: restore/replay MUST use the per-step historical snapshot,
    never a later expanded catalog or current rate cache."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.directory = tempfile.mkdtemp(prefix="spark-e2-ckpt-")
        cls.built = build_grown_catalog_checkpoint(cls.directory)
        cls.payload = read_checkpoint(cls.built.checkpoint_path)
        cls.arguments = cls.built.restore_arguments()

    def restore(self, payload):
        return validate_checkpoint_payload(payload, **self.arguments)

    def test_fixture_actually_discriminates_historical_from_current(self) -> None:
        # A pass is only evidence if the two tables genuinely disagree. Quantify the
        # disagreement before asserting anything about the restore path.
        built = self.built
        snapshot = self.payload["trajectory"][0]["rate_table_snapshot"]["payload"]
        self.assertEqual(len(snapshot["event_ids"]), 2)
        selectable_now = [event for event in built.engine.catalog.events.values()
                          if event.origin_state_id == built.origin.state_id and event.selectable]
        self.assertEqual(len(selectable_now), 3)

        self.assertNotEqual(built.historical.total_rate, built.current.total_rate)
        self.assertAlmostEqual(built.historical.total_rate, 17.762303659592, places=9)
        self.assertAlmostEqual(built.current.total_rate, 138.121865606364, places=9)

        # Discriminator 1: the selected event differs.
        self.assertNotEqual(built.selection.event_id, built.current_selected_event_id)
        # Discriminator 2: the time increment differs by 87.1 %, i.e. ~1.7e14 times
        # the `PAR-005` replay tolerance of relative 5e-15.
        relative = (abs(built.current_time_increment - built.selection.delta_time)
                    / built.selection.delta_time)
        self.assertGreater(relative, 0.87)
        self.assertGreater(relative / 5e-15, 1e14)

    def test_replay_accepts_the_historical_snapshot_after_catalog_growth(self) -> None:
        # ACCEPTED baseline. An implementation that rebuilt the table from the
        # checkpointed catalog would compute total rate 138.12 s^-1 instead of
        # 17.76 s^-1, select a different event, and reject this checkpoint.
        restored = self.restore(clone_payload(self.payload))
        built = self.built
        self.assertEqual(restored["step_index"], 1)
        self.assertEqual(restored["current_state_id"], built.selection.destination_state_id)
        self.assertEqual(restored["simulation_time"], built.selection.delta_time)
        self.assertEqual(len(restored["catalog"].events), 6)
        rebuilt = build_rate_table(tuple(restored["catalog"].events.values()),
                                   built.origin.state_id)
        self.assertEqual(len(rebuilt.event_ids), 3)
        self.assertNotEqual(rebuilt.total_rate,
                            self.payload["trajectory"][0]["total_rate_per_s"])

    def test_replay_rejects_the_time_increment_a_current_catalog_would_produce(self) -> None:
        # REJECTED counterpart 1: only `time_increment_s` is moved to the value the
        # current catalog implies. If replay silently used the current catalog this
        # corruption would be accepted.
        payload = clone_payload(self.payload)
        payload["trajectory"][0]["time_increment_s"] = self.built.current_time_increment
        payload["simulation_time_s"] = self.built.current_time_increment
        with self.assertRaises(DomainFailure) as failure:
            self.restore(payload)
        self.assertEqual(failure.exception.outcome.status, "CHECKPOINT_CORRUPT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-KMC-005")

    def test_replay_rejects_the_selection_a_current_catalog_would_produce(self) -> None:
        # REJECTED counterpart 2: the selection is moved to the current-catalog answer
        # while the historical snapshot stays in the record.
        payload = clone_payload(self.payload)
        step = payload["trajectory"][0]
        step["selected_rate_per_s"] = self.built.current_selected_rate
        step["total_rate_per_s"] = self.built.current.total_rate
        with self.assertRaises(DomainFailure) as failure:
            self.restore(payload)
        self.assertEqual(failure.exception.outcome.status, "CHECKPOINT_CORRUPT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-KMC-005")

    def test_replay_rejects_a_snapshot_replaced_by_the_grown_catalog(self) -> None:
        # REJECTED counterpart 3: the snapshot itself is replaced by the three-event
        # table. This pins that the snapshot, not the catalog, drives the replay: the
        # recorded uniforms no longer reproduce the recorded selection or time.
        payload = clone_payload(self.payload)
        payload["trajectory"][0]["rate_table_snapshot"] = self.built.current.snapshot()
        with self.assertRaises(DomainFailure) as failure:
            self.restore(payload)
        self.assertEqual(failure.exception.outcome.status, "CHECKPOINT_CORRUPT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-KMC-005")

    def test_snapshot_envelope_shape_and_hash_are_enforced(self) -> None:
        # `E2-KMC-002`: the envelope has exactly `payload` and `payload_sha256`.
        payload = clone_payload(self.payload)
        snapshot = payload["trajectory"][0]["rate_table_snapshot"]
        snapshot["payload"]["total_rate_per_s"] = snapshot["payload"]["total_rate_per_s"] * 2.0
        with self.assertRaises(DomainFailure) as failure:
            self.restore(payload)
        self.assertEqual(failure.exception.outcome.status, "CHECKPOINT_CORRUPT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-KMC-002")

    def test_snapshot_rates_must_agree_with_the_catalog_log_rates(self) -> None:
        # `E2-KMC-002` parallel arrays: a snapshot rate that is not exp(log_rate) of the
        # named event is corrupt even when its own hash is recomputed.
        payload = clone_payload(self.payload)
        snapshot = payload["trajectory"][0]["rate_table_snapshot"]
        snapshot["payload"]["rates"][0] = snapshot["payload"]["rates"][0] * 1.5
        snapshot["payload_sha256"] = digest(snapshot["payload"])
        with self.assertRaises(DomainFailure) as failure:
            self.restore(payload)
        self.assertEqual(failure.exception.outcome.status, "CHECKPOINT_CORRUPT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-KMC-002")


class DeepCheckpointInvariantTests(unittest.TestCase):
    """`E2-CKPT-001..009` recursive validate-before-mutation restore."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.directory = tempfile.mkdtemp(prefix="spark-e2-ckpt-deep-")
        cls.built = build_grown_catalog_checkpoint(cls.directory)
        cls.payload = read_checkpoint(cls.built.checkpoint_path)
        cls.arguments = cls.built.restore_arguments()

    def restore(self, payload):
        return validate_checkpoint_payload(payload, **self.arguments)

    def assert_baseline_accepted(self) -> None:
        self.assertEqual(self.restore(clone_payload(self.payload))["step_index"], 1)

    def assert_rejected(self, mutate, status: str, requirement: str) -> None:
        self.assert_baseline_accepted()
        payload = clone_payload(self.payload)
        mutate(payload)
        with self.assertRaises(DomainFailure) as failure:
            self.restore(payload)
        self.assertEqual(failure.exception.outcome.status, status)
        self.assertEqual(failure.exception.outcome.context["requirement_id"], requirement)

    def test_envelope_hash_and_canonical_form(self) -> None:
        # `E2-CKPT-001` plus `E2-CAN-001`: the on-disk bytes are canonical and the
        # envelope hash covers the canonical payload bytes.
        raw = Path(self.built.checkpoint_path).read_bytes()
        self.assertFalse(raw.endswith(b"\n"))
        self.assertEqual(read_checkpoint(self.built.checkpoint_path), self.payload)
        broken = self.built.checkpoint_path + ".broken"
        Path(broken).write_bytes(raw.replace(b'"payload_sha256":"sha256:', b'"payload_sha256":"sha256:0', 1))
        with self.assertRaises(DomainFailure) as failure:
            read_checkpoint(broken)
        self.assertEqual(failure.exception.outcome.status, "CHECKPOINT_CORRUPT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-CKPT-001")

    def test_payload_key_set_is_exact(self) -> None:
        # `E2-CKPT-002`: exactly fifteen keys and the `spark-atomistic-checkpoint/2` schema.
        self.assertEqual(set(self.payload), {
            "basin", "catalog", "checkpoint_sequence", "current_state", "digests",
            "discovery_statistics", "flags", "initial_state", "log_sequence", "resources",
            "rng", "schema", "simulation_time_s", "step_index", "trajectory"})
        self.assertEqual(self.payload["schema"], "spark-atomistic-checkpoint/2")
        self.assert_rejected(lambda payload: payload.pop("basin"),
                             "CHECKPOINT_CORRUPT", "E2-CKPT-002")

    def test_digest_mismatch_is_incompatible_not_corrupt(self) -> None:
        # A valid but mismatched run contract is `CHECKPOINT_INCOMPATIBLE`, a different
        # terminal status from damage. The governing requirement is `E2-CKPT-007` step 3
        # ("schema/config/model/tolerance digests"), NOT `E2-CKPT-003`: the latter fixes only
        # the key set of `digests`, which the mutation below leaves intact. Asserting the key
        # set first makes that explicit -- it is the precondition under which E2-CKPT-003 is
        # satisfied and E2-CKPT-007 is the only rule left to fire. Corrected 2026-08-11; this
        # was the last cross-language divergence among the 46 mandatory fixtures.
        self.assertEqual(set(self.payload["digests"]),
                         {"config", "model", "schema", "tolerances"})

        def mutate(payload):
            payload["digests"]["config"] = "sha256:" + "0" * 64

        self.assert_rejected(mutate, "CHECKPOINT_INCOMPATIBLE", "E2-CKPT-007")

    def test_basin_record_is_frozen_disabled(self) -> None:
        # `E2-CKPT-003` / `E2-BASIN-001`.
        self.assertEqual(self.payload["basin"], {"enabled": False, "reason": "v1-disabled"})
        self.assert_rejected(lambda payload: payload["basin"].update({"enabled": True}),
                             "CHECKPOINT_CORRUPT", "E2-BASIN-001")

    def test_flags_key_set_and_derived_consistency(self) -> None:
        # `E2-CKPT-004`: exact booleans and an exact status token, and the flags must
        # agree with what the restored content implies.
        self.assertEqual(set(self.payload["flags"]),
                         {"cancelled", "complete", "incomplete_catalog", "last_status",
                          "resource_limited"})
        self.assertEqual(self.payload["flags"]["last_status"], "OK")

        def mutate(payload):
            payload["flags"]["incomplete_catalog"] = True

        self.assert_rejected(mutate, "CHECKPOINT_CORRUPT", "E2-CKPT-004")

    def test_resource_counters_and_zero_retry_history(self) -> None:
        # `E2-CKPT-005`: retry history is empty in v1 because configured retry count is
        # zero, and `saddle_attempts_by_state` must equal the discovery attempt counts.
        resources = self.payload["resources"]
        self.assertEqual(set(resources), {
            "calculator_evaluations", "catalog_events", "output_bytes",
            "resident_memory_bytes", "retry_history", "saddle_attempts_by_state",
            "wall_elapsed_s"})
        self.assertEqual(resources["retry_history"], [])
        self.assertEqual(resources["catalog_events"], 6)
        self.assertEqual(resources["saddle_attempts_by_state"],
                         {self.built.origin.state_id: 3})
        self.assert_rejected(
            lambda payload: payload["resources"].__setitem__("retry_history", [{"attempt": 1}]),
            "CHECKPOINT_CORRUPT", "E2-CKPT-005")
        self.assert_rejected(
            lambda payload: payload["resources"].__setitem__("catalog_events", 4),
            "CHECKPOINT_CORRUPT", "E2-CKPT-005")

    def test_catalog_digest_covers_the_object_without_its_digest(self) -> None:
        # `E2-CKPT-006`.
        catalog = self.payload["catalog"]
        self.assertEqual(set(catalog), {"digest", "events", "multiplicity", "schema", "states"})
        self.assertEqual(catalog["schema"], "spark-atomistic-catalog/2")
        self.assertEqual(catalog["digest"], digest(
            {key: catalog[key] for key in ("events", "multiplicity", "schema", "states")}))
        self.assert_rejected(
            lambda payload: payload["catalog"]["multiplicity"].update(
                {next(iter(payload["catalog"]["events"])): 7}),
            "CHECKPOINT_CORRUPT", "E2-CKPT-006")

    def test_reciprocal_event_pair_is_recomputed_not_trusted(self) -> None:
        # `E2-EVENT-001` / `E2-EVENT-006`: identifiers are recomputed from content, so a
        # tampered stored ID cannot survive restore.
        def mutate(payload):
            events = payload["catalog"]["events"]
            victim = sorted(events)[0]
            events[victim]["reverse_event_id"] = "event:sha256:" + "0" * 64
            payload["catalog"]["digest"] = digest(
                {key: payload["catalog"][key]
                 for key in ("events", "multiplicity", "schema", "states")})

        self.assert_rejected(mutate, "CHECKPOINT_CORRUPT", "E2-EVENT-006")

    def test_barrier_and_rate_recomputation(self) -> None:
        # `E2-RATE-001`: raw barriers are exactly E_s - E_i, recomputed at restore.
        def mutate(payload):
            events = payload["catalog"]["events"]
            victim = sorted(events)[0]
            events[victim]["barrier_ev"] = events[victim]["barrier_ev"] + 0.05
            payload["catalog"]["digest"] = digest(
                {key: payload["catalog"][key]
                 for key in ("events", "multiplicity", "schema", "states")})

        self.assert_rejected(mutate, "CHECKPOINT_CORRUPT", "E2-RATE-001")

    def test_step_and_log_sequences_equal_trajectory_length(self) -> None:
        # `E2-CKPT-008`.
        self.assertEqual(self.payload["step_index"], len(self.payload["trajectory"]))
        self.assertEqual(self.payload["log_sequence"], len(self.payload["trajectory"]))
        self.assert_rejected(lambda payload: payload.__setitem__("log_sequence", 2),
                             "CHECKPOINT_CORRUPT", "E2-CKPT-008")

    def test_trajectory_rng_consumed_uniforms_equal_two_per_step(self) -> None:
        # `E2-CKPT-008`: consumed uniforms equal 2 * step_index, and the replayed
        # stream must equal the stored stream byte for byte.
        self.assertEqual(self.payload["rng"]["trajectory"]["consumed_uniforms"],
                         2 * self.payload["step_index"])
        self.assertEqual(self.payload["rng"]["trajectory"]["algorithm"],
                         "Philox4x32-10:errata-1-midpoint52")

        def mutate(payload):
            payload["rng"]["trajectory"]["consumed_uniforms"] = 4

        # Philox internal count relations are checked first (`E2-RNG-002`).
        self.assert_rejected(mutate, "CHECKPOINT_CORRUPT", "E2-RNG-002")

    def test_substream_map_must_be_complete_and_exact(self) -> None:
        # `E2-CKPT-007` item 8 plus `E2-RNG-005`: one class-selection stream and one
        # search stream per recorded attempt, and no extra keys.
        substreams = self.payload["rng"]["substream_map"]
        self.assertEqual(len(substreams), 6)
        self.assertEqual(sum(1 for key in substreams if key.startswith("class-selection:")), 3)

        def drop(payload):
            key = next(key for key in payload["rng"]["substream_map"]
                       if key.startswith("class-selection:"))
            payload["rng"]["substream_map"].pop(key)

        self.assert_rejected(drop, "CHECKPOINT_CORRUPT", "E2-RNG-005")

        def add(payload):
            payload["rng"]["substream_map"]["search:sha256:" + "0" * 64] = (
                payload["rng"]["trajectory"])

        self.assert_rejected(add, "CHECKPOINT_CORRUPT", "E2-CKPT-007")

    def test_discovery_counter_identities(self) -> None:
        # `E2-DISC-004` / `E2-DISC-005`: attempts = successes + failures, and duplicates
        # bound the consecutive redundant count.
        statistics = self.payload["discovery_statistics"][self.built.origin.state_id]
        self.assertEqual(statistics["attempts"], 3)
        self.assertEqual(statistics["successes"], 3)
        self.assertEqual(statistics["stopping_state"], "RUNNING")
        self.assertEqual(statistics["heuristic_confidence"], "UNAVAILABLE")
        self.assert_rejected(
            lambda payload: payload["discovery_statistics"][
                self.built.origin.state_id].__setitem__("successes", 2),
            "CHECKPOINT_CORRUPT", "E2-DISC-004")

    def test_crash_safe_write_is_atomic_replacement(self) -> None:
        # `E2-CKPT-009`: canonical bytes with no trailing newline, written through a
        # sibling temporary and an atomic replacement.
        raw = Path(self.built.checkpoint_path).read_bytes()
        self.assertFalse(raw.endswith(b"\n"))
        self.assertTrue(raw.startswith(b'{"payload":{'))
        leftovers = [name for name in os.listdir(self.directory) if ".tmp-" in name]
        self.assertEqual(leftovers, [])


class CheckpointRestoreOrderTests(unittest.TestCase):
    """Records a KNOWN DEVIATION from the `E2-CKPT-007` verification order.

    `E2-CKPT-007` lists step 8 as "Philox key/counter/buffer/count relations and
    complete substream map" and step 9 as "trajectory sequence/state chain using each
    historical rate snapshot". `validate_checkpoint_payload` performs the trajectory
    replay before the substream-map check, so a checkpoint damaged in BOTH places
    reports the step-9 requirement instead of the step-8 requirement.

    The terminal `status` and `exit_code` are unaffected (both paths are
    `CHECKPOINT_CORRUPT` / 74); only `context.requirement_id` differs. The ordering is
    deliberately NOT changed here: `context.requirement_id` is part of the public
    response that `E2-PAR-003` requires to be byte-identical across backends, so a
    unilateral Python change could create the very mismatch it aims to remove. The
    resolution belongs to a cross-language decision.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls.directory = tempfile.mkdtemp(prefix="spark-e2-ckpt-order-")
        cls.built = build_grown_catalog_checkpoint(cls.directory)
        cls.payload = read_checkpoint(cls.built.checkpoint_path)
        cls.arguments = cls.built.restore_arguments()

    def test_dual_defect_reports_the_trajectory_requirement(self) -> None:
        payload = clone_payload(self.payload)
        payload["trajectory"][0]["time_increment_s"] = self.built.current_time_increment
        payload["simulation_time_s"] = self.built.current_time_increment
        key = next(key for key in payload["rng"]["substream_map"]
                   if key.startswith("class-selection:"))
        payload["rng"]["substream_map"].pop(key)
        with self.assertRaises(DomainFailure) as failure:
            validate_checkpoint_payload(payload, **self.arguments)
        self.assertEqual(failure.exception.outcome.status, "CHECKPOINT_CORRUPT")
        # Spec order would report the substream-map failure (E2-RNG-005) first.
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-KMC-005")


if __name__ == "__main__":
    unittest.main()
