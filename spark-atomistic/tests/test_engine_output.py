# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Output byte reservation and fail-closed output transactions.

Governing requirements: `RES-001`/`RES-002`/`RES-004` (limits, limit-hit behaviour,
counters), `IO-003`/`IO-004` (path and overwrite policy validated before calculator
work), `E2-SCHEMA-009`/`E2-SCHEMA-010` (resource and output objects), and
`E2-CKPT-009` (sibling temporary, flush, atomic replacement, parent-directory flush).
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import fixtures
from fixtures import build_grown_catalog_checkpoint, checkpoint_model

from spark_atomistic import checkpoint as checkpoint_module
from spark_atomistic.canonical import deep_thaw
from spark_atomistic.engine import ReferenceEngine
from spark_atomistic.errors import DomainFailure
from spark_atomistic.model import validate_model
from spark_atomistic.resources import ResourceLedger


EXTENSION = {"calculator_command": ["/bin/true"]}


class OutputReservationTests(unittest.TestCase):
    """`RES-001`/`RES-004`: output bytes are a limited resource, reserved before use."""

    def setUp(self) -> None:
        self.directory = tempfile.mkdtemp(prefix="spark-e2-output-")
        self.built = build_grown_catalog_checkpoint(self.directory)

    def test_reservation_precedes_the_write_and_matches_it_exactly(self) -> None:
        engine = self.built.engine
        order: list[tuple[str, int]] = []
        real_reserve = ResourceLedger.reserve_output
        real_write = checkpoint_module._atomic_write

        def traced_reserve(ledger, count):
            order.append(("reserve_output", count))
            return real_reserve(ledger, count)

        def traced_write(path, encoded):
            order.append(("_atomic_write", len(encoded)))
            return real_write(path, encoded)

        reserved_before = engine.ledger.output_bytes_reserved
        written_before = engine.ledger.output_bytes_written
        with mock.patch.object(ResourceLedger, "reserve_output", traced_reserve), \
             mock.patch.object(checkpoint_module, "_atomic_write", traced_write):
            engine.write_checkpoint()

        self.assertEqual([name for name, _count in order],
                         ["reserve_output", "_atomic_write"])
        # The reservation is not an estimate: the self-referential `output_bytes` field
        # inside the payload is driven to a fixed point before the reservation is taken,
        # so reserved bytes equal written bytes exactly.
        self.assertEqual(order[0][1], order[1][1])
        self.assertEqual(engine.ledger.output_bytes_reserved - reserved_before, order[0][1])
        self.assertEqual(engine.ledger.output_bytes_written - written_before, order[1][1])
        self.assertEqual(engine.ledger.output_bytes_reserved,
                         engine.ledger.output_bytes_written)

    def test_recorded_output_bytes_account_for_the_checkpoint_being_written(self) -> None:
        # `E2-CKPT-005`: the counter stored inside the checkpoint includes that
        # checkpoint, which is what makes the fixed point necessary in the first place.
        payload = checkpoint_module.read_checkpoint(self.built.checkpoint_path)
        self.assertEqual(payload["resources"]["output_bytes"],
                         Path(self.built.checkpoint_path).stat().st_size)

    def test_output_limit_refuses_the_write_and_leaves_no_artifact(self) -> None:
        # `RES-002`: a limit hit returns `RESOURCE_LIMIT` and aborts the uncommitted
        # transaction. Paired with the accepted baseline above, which does write.
        limited_directory = tempfile.mkdtemp(prefix="spark-e2-output-limited-")
        model = deep_thaw(self.built.config)
        model["resources"]["output_bytes"] = 1024
        for key, name in (("checkpoint_path", "checkpoint.json"),
                          ("summary_path", "summary.json"),
                          ("trajectory_path", "trajectory.json")):
            model["output"][key] = str(Path(limited_directory) / name)
        engine = ReferenceEngine(validate_model(model), extension=EXTENSION)
        source = self.built.engine
        engine.catalog = source.catalog
        engine.current_state_id = source.current_state_id
        engine.initial_state_id = source.initial_state_id
        engine.trajectory_log = source.trajectory_log
        engine.substreams = source.substreams
        engine.trajectory_rng = source.trajectory_rng
        engine.step_index = 1
        engine.log_sequence = 1
        engine.simulation_time = source.simulation_time
        engine.ledger.per_state_saddle_attempts[source.initial_state_id] = 3

        sequence_before = engine.checkpoint_sequence
        with self.assertRaises(DomainFailure) as failure:
            engine.write_checkpoint()
        self.assertEqual(failure.exception.outcome.status, "RESOURCE_LIMIT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "RES-002")
        # Fail-closed: sequence rolled back, nothing reserved, nothing on disk at all.
        self.assertEqual(engine.checkpoint_sequence, sequence_before)
        self.assertEqual(engine.ledger.output_bytes_reserved, 0)
        self.assertEqual(engine.ledger.output_bytes_written, 0)
        self.assertEqual(sorted(os.listdir(limited_directory)), [])

    def test_failed_atomic_replace_preserves_the_previous_artifact(self) -> None:
        # `E2-CKPT-009`: crash-safe write is sibling temporary, flush, atomic
        # replacement, parent flush. A failure at the replacement boundary must leave
        # the previous checkpoint byte-identical and no temporary behind.
        engine = self.built.engine
        target = Path(self.built.checkpoint_path)
        original = target.read_bytes()
        sequence_before = engine.checkpoint_sequence
        with mock.patch.object(checkpoint_module.os, "replace",
                               side_effect=OSError(28, "No space left on device")):
            with self.assertRaises(DomainFailure) as failure:
                engine.write_checkpoint()
        self.assertEqual(failure.exception.outcome.status, "INTERNAL_ERROR")
        self.assertEqual(failure.exception.outcome.context["requirement_id"], "E2-CKPT-009")
        self.assertEqual(target.read_bytes(), original)
        self.assertEqual(engine.checkpoint_sequence, sequence_before)
        self.assertEqual([name for name in os.listdir(self.directory) if ".tmp-" in name], [])

    def test_failed_content_write_preserves_the_previous_artifact(self) -> None:
        # Same invariant one stage earlier: the failure happens while the sibling
        # temporary is being filled, so the destination is never touched.
        engine = self.built.engine
        target = Path(self.built.checkpoint_path)
        original = target.read_bytes()
        with mock.patch.object(checkpoint_module.os, "fsync",
                               side_effect=OSError(5, "I/O error")):
            with self.assertRaises(DomainFailure) as failure:
                engine.write_checkpoint()
        self.assertEqual(failure.exception.outcome.status, "INTERNAL_ERROR")
        self.assertEqual(target.read_bytes(), original)
        self.assertEqual([name for name in os.listdir(self.directory) if ".tmp-" in name], [])

    def test_auxiliary_outputs_reserve_before_writing(self) -> None:
        # `RES-004`: summary and trajectory are reserved and completed like the
        # checkpoint, and each byte count is exact rather than estimated.
        engine = self.built.engine
        reserved_before = engine.ledger.output_bytes_reserved
        engine._write_auxiliary_outputs()
        summary = Path(engine.config["output"]["summary_path"])
        trajectory = Path(engine.config["output"]["trajectory_path"])
        total = summary.stat().st_size + trajectory.stat().st_size
        self.assertEqual(engine.ledger.output_bytes_reserved - reserved_before, total)
        self.assertEqual(engine.ledger.output_bytes_reserved,
                         engine.ledger.output_bytes_written)
        self.assertEqual(json.loads(summary.read_bytes())["step_index"], 1)


class OutputCollisionPolicyTests(unittest.TestCase):
    """`IO-004`: only an existing checkpoint under `resume=true` may be tolerated."""

    def setUp(self) -> None:
        self.directory = tempfile.mkdtemp(prefix="spark-e2-collision-")

    def preflight(self, present, *, resume: bool, overwrite: bool):
        for name in os.listdir(self.directory):
            os.unlink(Path(self.directory) / name)
        for name in present:
            (Path(self.directory) / name).write_bytes(b"{}")
        config = checkpoint_model(self.directory, resume=resume, overwrite=overwrite)
        engine = ReferenceEngine(config, extension=EXTENSION)
        return engine._validate_output_preflight()

    def test_resume_tolerates_only_the_checkpoint_collision(self) -> None:
        # ACCEPTED baseline: an existing checkpoint with `resume=true` resumes.
        self.assertTrue(self.preflight(["checkpoint.json"], resume=True, overwrite=False))
        # ACCEPTED baseline: a clean directory is not a collision at all.
        self.assertFalse(self.preflight([], resume=True, overwrite=False))

        # REJECTED: "every other collision returns OUTPUT_EXISTS" (IO-004). Before the
        # 2026-08-11 fix the resume exemption was applied to all three paths at once,
        # so a pre-existing summary or trajectory was silently overwritten at the end
        # of a resumed run without `overwrite=true`.
        for present in (["checkpoint.json", "summary.json"],
                        ["checkpoint.json", "trajectory.json"],
                        ["checkpoint.json", "summary.json", "trajectory.json"],
                        ["summary.json"],
                        ["trajectory.json"]):
            with self.subTest(present=present):
                with self.assertRaises(DomainFailure) as failure:
                    self.preflight(present, resume=True, overwrite=False)
                self.assertEqual(failure.exception.outcome.status, "OUTPUT_EXISTS")
                self.assertEqual(failure.exception.outcome.context["requirement_id"],
                                 "E2-SCHEMA-010")

    def test_checkpoint_collision_without_resume_is_refused(self) -> None:
        with self.assertRaises(DomainFailure) as failure:
            self.preflight(["checkpoint.json"], resume=False, overwrite=False)
        self.assertEqual(failure.exception.outcome.status, "OUTPUT_EXISTS")

    def test_overwrite_permits_every_collision(self) -> None:
        # `IO-004` first sentence: `overwrite=true` is the blanket permission.
        self.assertTrue(self.preflight(
            ["checkpoint.json", "summary.json", "trajectory.json"],
            resume=False, overwrite=True))

    def test_missing_output_directory_is_rejected_before_calculator_work(self) -> None:
        # `IO-003`: output paths are validated before calculator work starts.
        config = checkpoint_model(Path(self.directory) / "absent")
        engine = ReferenceEngine(config, extension=EXTENSION)
        with self.assertRaises(DomainFailure) as failure:
            engine._validate_output_preflight()
        self.assertEqual(failure.exception.outcome.status, "INVALID_INPUT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"],
                         "E2-SCHEMA-010")


class AdapterBoundaryTests(unittest.TestCase):
    """`E2-SCOPE-004`: adapter transport settings live only under run-request
    `extension` and never enter the model or any digest."""

    def test_extension_is_excluded_from_the_config_digest(self) -> None:
        directory = tempfile.mkdtemp(prefix="spark-e2-extension-")
        config = checkpoint_model(directory)
        plain = ReferenceEngine(config, extension={"calculator_command": ["/bin/true"]})
        other = ReferenceEngine(config, extension={"calculator_command": ["/bin/false"],
                                                   "callback_stdout_bytes": 4096})
        self.assertEqual(plain.config_digest, other.config_digest)

    def test_relative_calculator_command_is_refused(self) -> None:
        directory = tempfile.mkdtemp(prefix="spark-e2-extension-bad-")
        config = checkpoint_model(directory)
        with self.assertRaises(DomainFailure) as failure:
            ReferenceEngine(config, extension={"calculator_command": ["true"]})
        self.assertEqual(failure.exception.outcome.status, "INVALID_INPUT")
        self.assertEqual(failure.exception.outcome.context["requirement_id"],
                         "E2-SCOPE-004")


if __name__ == "__main__":
    unittest.main()
