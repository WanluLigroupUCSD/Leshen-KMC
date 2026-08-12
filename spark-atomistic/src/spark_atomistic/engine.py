# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Errata-2 deterministic serial fixed-composition transaction engine."""

from __future__ import annotations

import math
import os
from pathlib import Path
import time
from typing import Any, Mapping, Sequence

from .basin import DISABLED_CHECKPOINT_RECORD, basin_disabled_outcome
from .calculator import ProcessCalculator
from .canonical import deep_freeze, deep_thaw, digest
from .catalog import Catalog, match_states
from .checkpoint import (checkpoint_encoded_size, read_checkpoint,
                         canonical_output_size, validate_checkpoint_payload, write_canonical_output,
                         write_checkpoint)
from .errors import DomainFailure, Outcome
from .kinetics import Selection, build_rate_table, propose_serial_step
from .model import (IR, SCHEMA_SHA256, AtomicState, identity_digest,
                    state_request_from_system)
from .resources import ResourceLedger
from .rng import PhiloxStream, derive_saddle_stream, derive_trajectory_stream
from .solvers import DirectionalDimerSearcher, SteepestDescentMinimizer


class ReferenceEngine:
    def __init__(self, config: Mapping[str, Any], *, extension: Mapping[str, Any]) -> None:
        self.config = config
        self.extension = extension
        behavior = deep_thaw(config)
        behavior.pop("metadata", None)
        self.config_digest = digest(behavior)
        self.model_digest = config["calculator"]["model_digest"]
        tolerance_payload = {key: config["kinetics"][key] for key in (
            "barrier_tolerance", "detailed_balance_tolerance", "log_rate_cutoff",
            "saddle_energy_tolerance", "saddle_max_tolerance", "saddle_rms_tolerance",
            "state_energy_tolerance_per_atom", "state_max_tolerance", "state_rms_tolerance")}
        self.tolerance_digest = digest(tolerance_payload)
        self.identity_digest = identity_digest(config["kinetics"])
        self.ledger = ResourceLedger.start(config["resources"])
        command = extension.get("calculator_command")
        if (not isinstance(command, (list, tuple)) or not command
                or any(not isinstance(item, str) or not item for item in command)
                or not os.path.isabs(command[0])):
            raise DomainFailure("INVALID_INPUT", "input invalid", component="adapter",
                                requirement="E2-SCOPE-004",
                                details={"extension_field": "calculator_command"})
        stdout_limit = extension.get("callback_stdout_bytes", config["resources"]["output_bytes"])
        stderr_limit = extension.get("callback_stderr_bytes", config["resources"]["output_bytes"])
        if (type(stdout_limit) is not int or stdout_limit < 1
                or type(stderr_limit) is not int or stderr_limit < 1):
            raise DomainFailure("INVALID_INPUT", "input invalid", component="adapter",
                                requirement="E2-SCOPE-004")
        adapter_config = {"command": list(command), "model_name": config["calculator"]["model_name"],
                          "model_version": config["calculator"]["model_version"],
                          "model_digest": self.model_digest,
                          "deterministic": config["calculator"]["deterministic"]}
        self.calculator = ProcessCalculator(
            adapter_config, config["resources"]["callback_timeout_s"],
            stdout_limit, stderr_limit, self.ledger)
        self.minimizer = SteepestDescentMinimizer(self.calculator, config["relaxation"])
        self.searcher = DirectionalDimerSearcher(self.calculator, self.minimizer,
                                                 config["saddle_search"])
        self.catalog = Catalog(self.model_digest, self.config_digest, self.tolerance_digest,
                               config["resources"]["catalog_events"], self.identity_digest)
        self.current_state_id = ""
        self.initial_state_id = ""
        self.simulation_time = 0.0
        self.step_index = 0
        self.log_sequence = 0
        self.checkpoint_sequence = 0
        self.trajectory_log: list[Mapping[str, Any]] = []
        self.substreams: dict[str, PhiloxStream] = {}
        self.trajectory_rng = derive_trajectory_stream(config["kinetics"]["run_seed"])
        self.last_status = "OK"
        self.completed = False
        self.cancelled = False
        self.resource_limited = False
        self.last_checkpoint_monotonic = time.monotonic()
        self.commit_generation = 0
        self.checkpointed_generation = -1

    def _validate_output_preflight(self) -> bool:
        paths = [Path(self.config["output"][key]) for key in
                 ("checkpoint_path", "summary_path", "trajectory_path")]
        for path in paths:
            if not path.parent.exists() or not path.parent.is_dir():
                raise DomainFailure("INVALID_INPUT", "input invalid", component="output",
                                    requirement="E2-SCHEMA-010",
                                    details={"path": str(path)})
        exists = [path.exists() for path in paths]
        # IO-004: "Existing outputs MUST NOT be overwritten unless `overwrite=true`.
        # Existing compatible checkpoint plus `resume=true` resumes; every other
        # collision returns `OUTPUT_EXISTS`." The resume exemption therefore covers the
        # checkpoint path only. A pre-existing summary or trajectory is "every other
        # collision" even when the checkpoint is resumable, because those two artifacts
        # are rewritten at completion and resuming would silently overwrite them.
        tolerated = (exists[0] and self.config["output"]["resume"], False, False)
        if (any(present and not allowed for present, allowed in zip(exists, tolerated))
                and not self.config["output"]["overwrite"]):
            raise DomainFailure("OUTPUT_EXISTS", "output exists", component="output",
                                requirement="E2-SCHEMA-010")
        return exists[0]

    def initialize_or_resume(self) -> None:
        checkpoint_exists = self._validate_output_preflight()
        if checkpoint_exists and self.config["output"]["resume"]:
            raw = read_checkpoint(self.config["output"]["checkpoint_path"])
            restored = validate_checkpoint_payload(
                raw, expected_config_digest=self.config_digest,
                expected_model_digest=self.model_digest,
                expected_tolerance_digest=self.tolerance_digest,
                expected_identity_digest=self.identity_digest,
                kinetics=self.config["kinetics"], saddle_config=self.config["saddle_search"],
                relaxation_config=self.config["relaxation"],
                discovery_config=self.config["discovery"],
                resource_config=self.config["resources"],
                maximum_events=self.config["resources"]["catalog_events"])
            ledger = ResourceLedger.start(self.config["resources"])
            ledger.restore_counts(restored["resources"])
            expected_attempts = {state_id: record.attempts
                                 for state_id, record in restored["catalog"].discovery.items()
                                 if record.attempts}
            if ledger.per_state_saddle_attempts != expected_attempts:
                raise DomainFailure("CHECKPOINT_CORRUPT", "checkpoint corrupt", component="checkpoint",
                                    requirement="E2-CKPT-007")
            self.catalog = restored["catalog"]
            self.current_state_id = restored["current_state_id"]
            self.initial_state_id = restored["initial_state_id"]
            self.simulation_time = restored["simulation_time"]
            self.step_index = restored["step_index"]
            self.log_sequence = restored["log_sequence"]
            self.checkpoint_sequence = restored["checkpoint_sequence"]
            self.trajectory_log = restored["trajectory"]
            self.trajectory_rng = restored["trajectory_rng"]
            self.substreams = restored["substreams"]
            self.last_status = restored["flags"]["last_status"]
            self.completed = restored["flags"]["complete"]
            self.cancelled = restored["flags"]["cancelled"]
            self.resource_limited = restored["flags"]["resource_limited"]
            self.commit_generation = self.step_index
            self.checkpointed_generation = self.step_index
            self.ledger = ledger
            self.calculator.ledger = ledger
            return
        request = state_request_from_system(self.config["system"])
        with self.calculator.evaluation_scope(
                "initial-relaxation", self.config["resources"]["evaluations_per_relaxation"]):
            result = self.minimizer.minimize(request, object_id="initial-relaxation")
        initial = self.catalog.add_state(result.state, self.config["kinetics"])
        self.current_state_id = initial.state_id
        self.initial_state_id = initial.state_id

    def _class_for_attempt(self, state_id: str, index: int) -> tuple[Mapping[str, Any], PhiloxStream]:
        stream = derive_saddle_stream(self.config["kinetics"]["run_seed"], state_id,
                                      "class-selection", index)
        uniform = stream.uniform()
        cumulative = 0.0
        classes = self.config["discovery"]["classes"]
        for entry in classes:
            cumulative += entry["probability"]
            if cumulative > uniform:
                return entry, stream
        return classes[-1], stream

    def _search_id(self, state_id: str, search_class: str, index: int) -> str:
        return "search:" + digest({"run_seed": self.config["kinetics"]["run_seed"],
                                    "search_class": search_class, "search_index": index,
                                    "state_id": state_id})

    def discover_current_state(self) -> Outcome:
        state = self.catalog.states[self.current_state_id]
        discovery = self.config["discovery"]
        record = self.catalog.discovery_record(state.state_id, discovery["mode"], discovery)
        if record.stopping_state == "CONVERGED_HEURISTIC":
            return Outcome("DISCOVERY_CONVERGED_HEURISTIC", "", {
                "component": "discovery", "details": {}, "requirement_id": "E2-DISC-004",
                "retryable": False, "search_or_event_id": None, "state_id": state.state_id})
        while (record.attempts < discovery["maximum_attempts"]
               and record.evaluations < discovery["maximum_evaluations"]):
            index = record.attempts
            self.ledger.reserve_saddle_attempt(state.state_id)
            search_class, class_stream = self._class_for_attempt(state.state_id, index)
            search_id = self._search_id(state.state_id, search_class["name"], index)
            stream = derive_saddle_stream(self.config["kinetics"]["run_seed"], state.state_id,
                                          search_class["name"], index)
            before = self.ledger.calculator_reserved
            candidate = None
            failure: DomainFailure | None = None
            try:
                active_hint: Sequence[int] | None = None
                if search_class["kind"] in {"local", "targeted"}:
                    movable = [atom for atom, active in enumerate(state.movable) if active]
                    active_hint = [movable[index % len(movable)]]
                remaining = discovery["maximum_evaluations"] - record.evaluations
                attempt_limit = min(self.config["resources"]["evaluations_per_saddle_attempt"], remaining)
                with self.calculator.evaluation_scope(search_id, attempt_limit):
                    candidate = self.searcher.search(state, stream, search_id=search_id,
                                                     active_hint=active_hint)
            except DomainFailure as exc:
                failure = exc
            finally:
                record.attempts += 1
                record.evaluations += self.ledger.calculator_reserved - before
                self.substreams["class-selection:" + search_id] = class_stream.clone()
                self.substreams[search_id] = stream.clone()
            if failure is not None:
                status = failure.outcome.status
                record.failures_by_status[status] = record.failures_by_status.get(status, 0) + 1
                record.consecutive_redundant_successes = 0
                if status in {"SADDLE_NOT_FOUND", "INVALID_SADDLE", "SADDLE_WRONG_BASIN",
                              "ENDPOINT_COLLAPSED", "RELAX_NOT_CONVERGED"}:
                    continue
                raise failure
            if candidate is None:
                raise DomainFailure("INTERNAL_ERROR", "internal error", component="discovery",
                                    requirement="E2-DISC-005", state_id=state.state_id)
            try:
                outcome, event, _destination = self.catalog.validate_candidate(
                    state, candidate, self.config["kinetics"],
                    {"rng_substream_digest": digest(stream.checkpoint()),
                     "search_class": search_class["name"], "search_id": search_id,
                     "search_index": index})
            except DomainFailure as exc:
                record.failures_by_status[exc.outcome.status] = record.failures_by_status.get(exc.outcome.status, 0) + 1
                record.consecutive_redundant_successes = 0
                raise
            if outcome.status in {"SADDLE_WRONG_BASIN", "ENDPOINT_COLLAPSED"}:
                record.failures_by_status[outcome.status] = record.failures_by_status.get(outcome.status, 0) + 1
                record.consecutive_redundant_successes = 0
                continue
            record.successes += 1
            if outcome.status == "DUPLICATE_EVENT":
                record.duplicates += 1
                duplicate = self.catalog.events[outcome.context["search_or_event_id"]]
                relevant = (record.relevance_rate_min == 0.0
                            or duplicate.log_rate >= math.log(record.relevance_rate_min))
                record.consecutive_redundant_successes = (
                    record.consecutive_redundant_successes + 1 if relevant else 0)
            elif event is not None:
                if (record.relevance_rate_min == 0.0
                        or event.log_rate >= math.log(record.relevance_rate_min)):
                    record.event_log_rates[event.event_id] = event.log_rate
                record.consecutive_redundant_successes = 0
            if (record.successes >= discovery["minimum_successful"]
                    and record.consecutive_redundant_successes >= discovery["consecutive_redundant"]):
                record.stopping_state = "CONVERGED_HEURISTIC"
                if record.alpha is not None:
                    record.heuristic_confidence = 1.0 - 1.0 / (
                        record.alpha * record.consecutive_redundant_successes)
                return Outcome("DISCOVERY_CONVERGED_HEURISTIC", "", {
                    "component": "discovery", "details": {}, "requirement_id": "E2-DISC-004",
                    "retryable": False, "search_or_event_id": None, "state_id": state.state_id})
        record.stopping_state = "INCOMPLETE"
        record.permanently_incomplete_catalog = discovery["mode"] == "exploratory"
        return Outcome("DISCOVERY_INCOMPLETE", "", {
            "component": "discovery", "details": {}, "requirement_id": "E2-DISC-004",
            "retryable": False, "search_or_event_id": None, "state_id": state.state_id})

    def _verify_application(self, destination_state_id: str, event_id: str) -> AtomicState:
        """Reconverge the event's own saddle displacement and match the destination.

        `EVENT-004`: "Event application MUST use the validated destination minimum,
        followed by one verification relaxation. Failure to recover the destination
        within state tolerances returns `EVENT_APPLICATION_FAILED`." The validated
        destination minimum is what the relaxation must RECOVER, not what it starts
        from: a relaxation launched at the destination's own coordinates returns at step
        zero, because `RELAX-003` already committed that state at
        `max_movable_force <= relaxation.force_tolerance`, and the match then compares
        the destination with itself. The application geometry is therefore rebuilt the
        way `E2-EVENT-002`/`catalog.validate_candidate` built the endpoints in the first
        place, from the committed saddle geometry displaced along the committed unstable
        direction by `saddle_search.endpoint_displacement`.

        Both signs are tried because no requirement fixes the sign of
        `saddle.unstable_direction` relative to the destination: `E2-EVENT-006`
        canonicalises the sign only inside the `pair_id` hash, and
        `catalog.validate_candidate` picks whichever endpoint did NOT match the origin
        without recording which one that was. The first sign that recovers the
        destination is the verification relaxation `EVENT-004` names; the second is
        attempted only to resolve that unrecorded sign.
        """
        destination = self.catalog.states[destination_state_id]
        event = self.catalog.events.get(event_id)
        if (event is None or event.destination_state_id != destination_state_id
                or event.origin_state_id != self.current_state_id):
            raise DomainFailure("EVENT_APPLICATION_FAILED", "event application failed",
                                component="kinetics", requirement="E2-KMC-003",
                                state_id=self.current_state_id, object_id=event_id)
        origin = self.catalog.states[event.origin_state_id]
        distance = self.config["saddle_search"]["endpoint_displacement"]
        system = _state_as_system(origin)
        causal: str | None = None
        for sign, suffix in ((1.0, ":plus"), (-1.0, ":minus")):
            object_id = "apply-" + event_id + suffix
            system["positions"] = [
                [coordinate + sign * distance * component
                 for coordinate, component in zip(position, direction)]
                for position, direction in zip(event.saddle_positions, event.unstable_direction)]
            request = state_request_from_system(system)
            try:
                with self.calculator.evaluation_scope(
                        object_id, self.config["resources"]["evaluations_per_relaxation"]):
                    verified = self.minimizer.minimize(request, object_id=object_id).state
            except DomainFailure as exc:
                if exc.outcome.status in {"CALCULATOR_FAILURE", "RESOURCE_LIMIT", "CANCELLED"}:
                    raise
                causal = causal or exc.outcome.status
                continue
            if match_states(destination, verified, self.config["kinetics"]).equal:
                return destination
        raise DomainFailure("EVENT_APPLICATION_FAILED", "event application failed",
                            component="kinetics", requirement="E2-KMC-003",
                            state_id=self.current_state_id, object_id=event_id,
                            causal_status=causal)

    def _serial_proposal(self) -> Selection:
        table = build_rate_table(tuple(self.catalog.events.values()), self.current_state_id)
        return propose_serial_step(table, self.trajectory_rng)

    def commit_one_step(self) -> None:
        pre_state = self.current_state_id
        selection = self._serial_proposal()
        destination = self._verify_application(selection.destination_state_id, selection.event_id)
        new_time = self.simulation_time + selection.delta_time
        if not math.isfinite(new_time) or new_time <= self.simulation_time:
            raise DomainFailure("RATE_INVALID", "rate invalid", component="kinetics",
                                requirement="E2-KMC-003", state_id=pre_state,
                                object_id=selection.event_id)
        entry = deep_freeze({"checkpoint_sequence": self.checkpoint_sequence,
                 "log_sequence": self.log_sequence + 1,
                 "post_state_id": destination.state_id, "pre_state_id": pre_state,
                 "rate_table_snapshot": selection.rate_table_snapshot,
                 "selected_event_id": selection.event_id,
                 "selected_rate_per_s": selection.selected_rate,
                 "selection_uniform": selection.selection_uniform,
                 "step_index": self.step_index + 1,
                 "time_increment_s": selection.delta_time,
                 "time_uniform": selection.time_uniform,
                 "total_rate_per_s": selection.total_rate})
        self.trajectory_log.append(entry)
        self.current_state_id = destination.state_id
        self.simulation_time = new_time
        self.trajectory_rng.commit_from(selection.rng_after)
        self.step_index += 1
        self.log_sequence += 1
        self.commit_generation += 1

    def checkpoint_payload(self, resource_record: Mapping[str, Any] | None = None) -> dict[str, Any]:
        incomplete = any(record.stopping_state == "INCOMPLETE"
                         for record in self.catalog.discovery.values())
        return {"basin": dict(DISABLED_CHECKPOINT_RECORD),
                "catalog": self.catalog.checkpoint(),
                "checkpoint_sequence": self.checkpoint_sequence,
                "current_state": self.catalog.states[self.current_state_id].record(),
                "digests": {"config": self.config_digest, "model": self.model_digest,
                            "schema": SCHEMA_SHA256, "tolerances": self.tolerance_digest},
                "discovery_statistics": {key: self.catalog.discovery[key].record()
                                         for key in sorted(self.catalog.discovery)},
                "flags": {"cancelled": self.cancelled, "complete": self.completed,
                          "incomplete_catalog": incomplete, "last_status": self.last_status,
                          "resource_limited": self.resource_limited},
                "initial_state": self.catalog.states[self.initial_state_id].record(),
                "log_sequence": self.log_sequence,
                "resources": (self.ledger.checkpoint(catalog_events=len(self.catalog.events))
                              if resource_record is None else dict(resource_record)),
                "rng": {"run_seed": self.config["kinetics"]["run_seed"],
                        "substream_map": {key: self.substreams[key].checkpoint()
                                          for key in sorted(self.substreams)},
                        "trajectory": self.trajectory_rng.checkpoint()},
                "schema": "spark-atomistic-checkpoint/2",
                "simulation_time_s": self.simulation_time,
                "step_index": self.step_index,
                "trajectory": list(self.trajectory_log)}

    def write_checkpoint(self) -> None:
        self.checkpoint_sequence += 1
        try:
            wall = self.ledger.wall_elapsed()
            size = 0
            for _ in range(12):
                resources = self.ledger.checkpoint(
                    wall_elapsed=wall, projected_written=self.ledger.output_bytes_written + size,
                    catalog_events=len(self.catalog.events))
                candidate_size = checkpoint_encoded_size(self.checkpoint_payload(resources))
                if candidate_size == size:
                    break
                size = candidate_size
            else:
                raise DomainFailure("INTERNAL_ERROR", "internal error", component="checkpoint",
                                    requirement="E2-CKPT-009")
            self.ledger.reserve_output(size)
            resources = self.ledger.checkpoint(
                wall_elapsed=wall, projected_written=self.ledger.output_bytes_written + size,
                catalog_events=len(self.catalog.events))
            payload = self.checkpoint_payload(resources)
            written = write_checkpoint(self.config["output"]["checkpoint_path"], payload,
                                       byte_limit=self.config["resources"]["output_bytes"])
            self.ledger.complete_output(written)
            self.last_checkpoint_monotonic = time.monotonic()
            self.checkpointed_generation = self.commit_generation
        except BaseException:
            self.checkpoint_sequence -= 1
            raise

    def _write_auxiliary_outputs(self) -> None:
        summary = self.public_summary()
        trajectory = self.trajectory_log
        for path, value in ((self.config["output"]["summary_path"], summary),
                            (self.config["output"]["trajectory_path"], trajectory)):
            size = canonical_output_size(value)
            self.ledger.reserve_output(size)
            written = write_canonical_output(path, value,
                                             byte_limit=self.config["resources"]["output_bytes"])
            if written != size:
                raise DomainFailure("INTERNAL_ERROR", "internal error", component="output",
                                    requirement="E2-CKPT-009")
            self.ledger.complete_output(written)

    def has_valid_last_checkpoint(self) -> bool:
        return self.checkpointed_generation == self.commit_generation

    def run(self) -> Outcome:
        self.initialize_or_resume()
        if self.completed:
            return Outcome("OK", "", {"component": "engine", "details": {},
                           "requirement_id": "E2-CKPT-008", "retryable": False,
                           "search_or_event_id": None, "state_id": self.current_state_id})
        if self.config["basin"]["enabled"]:
            _disabled = basin_disabled_outcome(self.current_state_id)
        while self.step_index < self.config["kinetics"]["maximum_steps"]:
            try:
                self.ledger.check_wall_time()
            except DomainFailure:
                self.resource_limited = True
                self.last_status = "RESOURCE_LIMIT"
                raise
            discovery = self.discover_current_state()
            self.last_status = discovery.status
            if discovery.status == "DISCOVERY_INCOMPLETE" and self.config["discovery"]["mode"] == "strict":
                self.write_checkpoint()
                return discovery
            try:
                self.commit_one_step()
            except DomainFailure as exc:
                if exc.outcome.status == "NO_ENABLED_EVENT" and self.config["kinetics"]["absorbing_ok"]:
                    self.last_status = "NO_ENABLED_EVENT"
                    self.completed = True
                    self.write_checkpoint()
                    self._write_auxiliary_outputs()
                    return exc.outcome
                raise
            self.last_status = "OK"
            if (self.step_index % self.config["output"]["checkpoint_every_steps"] == 0
                    or time.monotonic() - self.last_checkpoint_monotonic
                    >= self.config["output"]["checkpoint_wall_time_s"]):
                self.write_checkpoint()
        self.completed = True
        self.last_status = "OK"
        self.write_checkpoint()
        self._write_auxiliary_outputs()
        return Outcome("OK", "", {"component": "engine", "details": {},
                       "requirement_id": "E2-KMC-004", "retryable": False,
                       "search_or_event_id": None, "state_id": self.current_state_id})

    def public_summary(self) -> dict[str, Any]:
        incomplete = any(record.stopping_state == "INCOMPLETE"
                         for record in self.catalog.discovery.values())
        return {"checkpoint_sequence": self.checkpoint_sequence,
                "current_state_id": self.current_state_id,
                "incomplete_catalog": incomplete,
                "simulation_time_s": self.simulation_time,
                "step_index": self.step_index}


def _state_as_system(state: AtomicState) -> dict[str, Any]:
    return {"atom_ids": list(state.atom_ids), "species": list(state.species),
            "positions": [list(item) for item in state.positions],
            "cell": [list(item) for item in state.cell], "pbc": list(state.pbc),
            "movable": list(state.movable), "constraints": {"kind": "fixed-mask"},
            "charge": state.charge, "spin": state.spin,
            "calculator_model_digest": state.calculator_model_digest}
