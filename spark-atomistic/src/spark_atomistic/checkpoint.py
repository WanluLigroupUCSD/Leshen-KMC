# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Errata-2 canonical checkpoint wire and recursive validate-before-mutation restore."""

from __future__ import annotations

import math
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping

from .canonical import canonical_bytes, deep_freeze, digest, parse_json
from .catalog import Catalog, DiscoveryRecord, Event
from .errors import DomainFailure, STATUSES
from .geometry import max_movable_force
from .model import (AtomicState, SCHEMA_SHA256, geometry_certificate,
                    validate_state_record)
from .rng import PhiloxStream, derive_saddle_stream, derive_trajectory_stream


CHECKPOINT_SCHEMA = "spark-atomistic-checkpoint/2"


def _failure(status: str, requirement: str, *, causal: str | None = None) -> DomainFailure:
    return DomainFailure(status, "", component="checkpoint", requirement=requirement,
                         causal_status=causal)


def checkpoint_encoded_size(payload: Mapping[str, Any]) -> int:
    return len(canonical_bytes({"payload": dict(payload),
                                "payload_sha256": digest(payload)}))


def canonical_output_size(value: Any) -> int:
    return len(canonical_bytes(value))


def _atomic_write(path: str, encoded: bytes) -> int:
    destination = Path(path)
    descriptor = -1
    temporary = ""
    try:
        descriptor, temporary = tempfile.mkstemp(prefix=destination.name + ".tmp-",
                                                 dir=str(destination.parent))
        with os.fdopen(descriptor, "wb", closefd=True) as handle:
            descriptor = -1
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        temporary = ""
        directory_fd = os.open(destination.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        return len(encoded)
    except OSError as exc:
        raise DomainFailure("INTERNAL_ERROR", "internal error", component="checkpoint",
                            requirement="E2-CKPT-009",
                            details={"io_error": type(exc).__name__}) from exc
    finally:
        if descriptor >= 0:
            try:
                os.close(descriptor)
            except OSError:
                pass
        if temporary:
            try:
                os.unlink(temporary)
            except OSError:
                pass


def write_checkpoint(path: str, payload: Mapping[str, Any], *, byte_limit: int) -> int:
    encoded = canonical_bytes({"payload": dict(payload),
                               "payload_sha256": digest(payload)})
    if len(encoded) > byte_limit:
        raise DomainFailure("RESOURCE_LIMIT", "resource limit reached", component="checkpoint",
                            requirement="E2-SCHEMA-009")
    return _atomic_write(path, encoded)


def write_canonical_output(path: str, value: Any, *, byte_limit: int) -> int:
    encoded = canonical_bytes(value)
    if len(encoded) > byte_limit:
        raise DomainFailure("RESOURCE_LIMIT", "resource limit reached", component="output",
                            requirement="E2-SCHEMA-009")
    return _atomic_write(path, encoded)


def read_checkpoint(path: str) -> Any:
    try:
        raw = Path(path).read_bytes()
    except OSError as exc:
        raise DomainFailure("INTERNAL_ERROR", "internal error", component="checkpoint",
                            requirement="E2-CKPT-009",
                            details={"io_error": type(exc).__name__}) from exc
    try:
        envelope = parse_json(raw)
        if canonical_bytes(envelope) != raw:
            raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-001")
        if (not isinstance(envelope, dict) or set(envelope) != {"payload", "payload_sha256"}
                or not isinstance(envelope["payload_sha256"], str)
                or envelope["payload_sha256"] != digest(envelope["payload"])):
            raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-001")
        return envelope["payload"]
    except DomainFailure as exc:
        if exc.outcome.status in {"INVALID_INPUT", "NONFINITE_RESULT"}:
            raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-001",
                           causal=exc.outcome.status) from exc
        raise


def _finite(value: Any, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-007")
    result = float(value)
    if nonnegative and result < 0.0:
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-007")
    return result


def _vectors(value: Any, count: int) -> tuple[tuple[float, float, float], ...]:
    if not isinstance(value, list) or len(value) != count:
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-007")
    output = []
    for row in value:
        if not isinstance(row, list) or len(row) != 3:
            raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-007")
        output.append(tuple(_finite(item) for item in row))
    return tuple(output)  # type: ignore[return-value]


def _restore_event(value: Any, states: Mapping[str, AtomicState],
                   model_digest: str, tolerance_digest: str, identity_digest: str,
                   kinetics: Mapping[str, Any], saddle_config: Mapping[str, Any]) -> Event:
    required = {"active_atom_mapping", "barrier_ev", "calculator_digest",
                "destination_state_id", "discovery_provenance", "environment_key",
                "environment_version", "event_id", "identity_digest", "origin_state_id",
                "pair_id", "rate_model", "reverse_barrier_ev", "reverse_event_id",
                "saddle", "schema", "schema_digest", "selectable", "tolerance_digest",
                "validation"}
    if not isinstance(value, dict) or set(value) != required:
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-001")
    origin_id, destination_id = value["origin_state_id"], value["destination_state_id"]
    if (not isinstance(origin_id, str) or not isinstance(destination_id, str)
            or origin_id not in states or destination_id not in states or origin_id == destination_id
            or value["schema"] != "spark-atomistic-directed-event/2"
            or value["schema_digest"] != SCHEMA_SHA256
            or value["calculator_digest"] != model_digest
            or value["tolerance_digest"] != tolerance_digest
            or value["identity_digest"] != identity_digest
            or value["environment_key"] != "disabled"
            or value["environment_version"] != "none/1"):
        raise _failure("CHECKPOINT_INCOMPATIBLE", "E2-EVENT-001")
    origin, destination = states[origin_id], states[destination_id]
    saddle = value["saddle"]
    saddle_keys = {"curvature_ev_per_angstrom2", "energy_ev", "evaluation_count",
                   "evidence_level", "forces_ev_per_angstrom",
                   "orthogonal_curvatures_ev_per_angstrom2", "positions", "search_id",
                   "termination_reason", "unstable_direction"}
    count = len(origin.atom_ids)
    if (not isinstance(saddle, dict) or set(saddle) != saddle_keys
            or saddle["evidence_level"] not in {"HESSIAN", "DIRECTIONAL"}
            or type(saddle["evaluation_count"]) is not int or saddle["evaluation_count"] < 0
            or not isinstance(saddle["search_id"], str) or not saddle["search_id"]
            or not isinstance(saddle["termination_reason"], str) or not saddle["termination_reason"]):
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-002")
    positions = _vectors(saddle["positions"], count)
    forces = _vectors(saddle["forces_ev_per_angstrom"], count)
    mode = _vectors(saddle["unstable_direction"], count)
    if abs(math.sqrt(sum(item * item for row in mode for item in row)) - 1.0) > 1e-10:
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-002")
    curvature = _finite(saddle["curvature_ev_per_angstrom2"])
    orthogonal_raw = saddle["orthogonal_curvatures_ev_per_angstrom2"]
    if (curvature >= -saddle_config["curvature_tolerance"]
            or not isinstance(orthogonal_raw, list)
            or len(orthogonal_raw) != saddle_config["orthogonal_directions"]
            or any(_finite(item) < -saddle_config["curvature_tolerance"] for item in orthogonal_raw)
            or max_movable_force(forces, origin.movable) > saddle_config["force_tolerance"]):
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-002")
    energy = _finite(saddle["energy_ev"])
    barrier = _finite(value["barrier_ev"])
    reverse_barrier = _finite(value["reverse_barrier_ev"])
    if (barrier != energy - origin.energy or reverse_barrier != energy - destination.energy
            or barrier < -kinetics["barrier_tolerance"]
            or reverse_barrier < -kinetics["barrier_tolerance"]):
        raise _failure("CHECKPOINT_CORRUPT", "E2-RATE-001")
    rate = value["rate_model"]
    rate_keys = {"common_prefactor_per_s", "detailed_balance_residual",
                 "log_forward_rate_per_s", "log_reverse_rate_per_s", "model", "temperature_k"}
    if not isinstance(rate, dict) or set(rate) != rate_keys or rate["model"] != "COMMON_PREFACTOR":
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-003")
    beta = 1.0 / (8.617333262145e-5 * kinetics["temperature"])
    expected_forward = math.log(kinetics["prefactor"]) - beta * barrier
    expected_reverse = math.log(kinetics["prefactor"]) - beta * reverse_barrier
    # `E2-CKPT-007`(5) requires restore to recompute detailed balance. Here that
    # recomputation HAS power, unlike the commit-time one deleted from
    # `catalog.validate_candidate`: `rate["detailed_balance_residual"]` is an independent
    # value read from the file and is compared against the value implied by the stored
    # barriers and energies. The `abs(expected_residual) > tolerance` term immediately
    # below is NOT independent — the two `!=` checks four lines up already pinned
    # `barrier`/`reverse_barrier` to the same-saddle differences of `E2-RATE-001`, which
    # makes `expected_residual` the same algebraic identity, bounded by binary64 rounding
    # (measured 1.1368683772161603e-13 worst case). It is retained because the tolerance
    # is the specification's own and costs nothing, but it must not be read as evidence
    # that detailed balance was checked; the `!=` on the stored residual is that evidence.
    expected_residual = expected_forward - expected_reverse + beta * (destination.energy - origin.energy)
    if (rate["common_prefactor_per_s"] != kinetics["prefactor"]
            or rate["temperature_k"] != kinetics["temperature"]
            or rate["log_forward_rate_per_s"] != expected_forward
            or rate["log_reverse_rate_per_s"] != expected_reverse
            or rate["detailed_balance_residual"] != expected_residual
            or abs(expected_residual) > kinetics["detailed_balance_tolerance"]
            or type(value["selectable"]) is not bool
            or value["selectable"] != (expected_forward >= kinetics["log_rate_cutoff"])):
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-003")
    active_raw = value["active_atom_mapping"]
    if (not isinstance(active_raw, list)
            or any(not isinstance(item, list) or len(item) != 2
                   or any(type(index) is not int for index in item) for item in active_raw)):
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-001")
    active = tuple((item[0], item[1]) for item in active_raw)
    if (list(active) != sorted(set(active))
            or any(not 0 <= left < count or not 0 <= right < count
                   or origin.species[left] != destination.species[right] for left, right in active)
            or len({left for left, _right in active}) != len(active)
            or len({right for _left, right in active}) != len(active)):
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-001")
    provenance = value["discovery_provenance"]
    if (not isinstance(provenance, dict)
            or set(provenance) != {"rng_substream_digest", "search_class", "search_id", "search_index"}
            or provenance["search_id"] != saddle["search_id"]
            or not isinstance(provenance["rng_substream_digest"], str)
            or not isinstance(provenance["search_class"], str)
            or type(provenance["search_index"]) is not int or provenance["search_index"] < 0):
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-004")
    validation = value["validation"]
    match_keys = {"atom_mapping", "energy_difference_ev", "max_displacement_angstrom",
                  "rms_displacement_angstrom"}
    validation_keys = {"calculator_model_digest", "constraint_digest", "destination_match",
                       "full_endpoint_relaxations", "method", "origin_match", "unstable_mode_count"}
    if (not isinstance(validation, dict) or set(validation) != validation_keys
            or validation["calculator_model_digest"] != model_digest
            or validation["constraint_digest"] != origin.constraint_digest
            or validation["full_endpoint_relaxations"] is not True
            or validation["unstable_mode_count"] != 1
            or not isinstance(validation["method"], str) or not validation["method"]
            or any(not isinstance(validation[key], dict) or set(validation[key]) != match_keys
                   for key in ("destination_match", "origin_match"))):
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-005")
    for match in (validation["destination_match"], validation["origin_match"]):
        atom_mapping = match["atom_mapping"]
        if (not isinstance(atom_mapping, list) or len(atom_mapping) != count
                or any(type(index) is not int or not 0 <= index < count for index in atom_mapping)
                or sorted(atom_mapping) != list(range(count))
                or _finite(match["energy_difference_ev"], nonnegative=True) < 0.0
                or _finite(match["max_displacement_angstrom"], nonnegative=True) < 0.0
                or _finite(match["rms_displacement_angstrom"], nonnegative=True) < 0.0):
            raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-005")
    certificate = geometry_certificate({"positions": origin.positions, "species": origin.species,
                                        "movable": origin.movable, "cell": origin.cell,
                                        "pbc": origin.pbc}, positions)
    saddle_geometry_digest = digest(certificate)
    flat = [component for vector in mode for component in vector]
    canonical_mode = min((flat, [-component for component in flat]), key=canonical_bytes)
    endpoints = sorted((origin_id, destination_id))
    oriented_active = active if origin_id == endpoints[0] else tuple(sorted((right, left) for left, right in active))
    pair_id = "pair:" + digest({"active_atom_mapping": [list(item) for item in oriented_active],
                                 "endpoint_state_ids": endpoints, "saddle_energy_ev": energy,
                                 "saddle_geometry_digest": saddle_geometry_digest,
                                 "schema": "spark-atomistic-event-pair/2",
                                 "unstable_direction": canonical_mode})
    event_id = "event:" + digest({"destination_state_id": destination_id,
                                   "origin_state_id": origin_id, "pair_id": pair_id,
                                   "schema": "spark-atomistic-directed-event/2"})
    reverse_id = "event:" + digest({"destination_state_id": origin_id,
                                     "origin_state_id": destination_id, "pair_id": pair_id,
                                     "schema": "spark-atomistic-directed-event/2"})
    if value["pair_id"] != pair_id or value["event_id"] != event_id or value["reverse_event_id"] != reverse_id:
        raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-006")
    return Event(event_id, reverse_id, pair_id, origin_id, destination_id, positions, energy,
                 forces, mode, curvature, tuple(float(item) for item in orthogonal_raw),
                 saddle["evidence_level"], saddle["evaluation_count"], saddle["search_id"],
                 saddle["termination_reason"], barrier, reverse_barrier, expected_forward,
                 expected_reverse, expected_residual, kinetics["prefactor"],
                 kinetics["temperature"], active, deep_freeze(provenance),
                 deep_freeze(validation), model_digest, identity_digest, tolerance_digest,
                 value["selectable"])


def _restore_discovery(value: Any, state_id: str, config_digest: str,
                       config: Mapping[str, Any], catalog: Catalog) -> DiscoveryRecord:
    required = {"alpha", "alpha_calibration", "attempts", "config_digest",
                "consecutive_redundant_successes", "duplicates", "evaluations",
                "event_log_rates", "failures_by_status", "heuristic_confidence",
                "permanently_incomplete_catalog", "relevance_rate_min", "state_id",
                "stopping_state", "successes"}
    if not isinstance(value, dict) or set(value) != required:
        raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-004")
    if (not isinstance(value["failures_by_status"], dict)
            or not isinstance(value["event_log_rates"], dict)):
        raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-004")
    integer_keys = ("attempts", "consecutive_redundant_successes", "duplicates",
                    "evaluations", "successes")
    if (value["state_id"] != state_id or value["config_digest"] != config_digest
            or any(type(value[key]) is not int or value[key] < 0 for key in integer_keys)
            or type(value["permanently_incomplete_catalog"]) is not bool
            or value["stopping_state"] not in {"RUNNING", "CONVERGED_HEURISTIC", "INCOMPLETE"}
            or value["attempts"] != value["successes"] + sum(value["failures_by_status"].values())
            or value["duplicates"] > value["successes"]
            or value["consecutive_redundant_successes"] > value["duplicates"]
            or value["attempts"] > config["maximum_attempts"]
            or value["evaluations"] > config["maximum_evaluations"]
            or value["relevance_rate_min"] != config["relevance_rate_min"]
            or value["alpha"] != config["alpha"]
            or value["alpha_calibration"] != (None if config["alpha_calibration"] is None
                                               else dict(config["alpha_calibration"]))
            or any(status not in STATUSES or type(count) is not int or count < 0
                   for status, count in value["failures_by_status"].items())
            or any(not isinstance(event_id, str) or not event_id
                   for event_id in value["event_log_rates"])):
        raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-004")
    confidence = value["heuristic_confidence"]
    if confidence != "UNAVAILABLE" and not math.isfinite(_finite(confidence)):
        raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-004")
    if value["stopping_state"] == "CONVERGED_HEURISTIC":
        if (value["successes"] < config["minimum_successful"]
                or value["consecutive_redundant_successes"] < config["consecutive_redundant"]
                or value["permanently_incomplete_catalog"]):
            raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-004")
    elif value["stopping_state"] == "INCOMPLETE":
        if (confidence != "UNAVAILABLE"
                or (value["attempts"] < config["maximum_attempts"]
                    and value["evaluations"] < config["maximum_evaluations"])
                or value["permanently_incomplete_catalog"] != (config["mode"] == "exploratory")):
            raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-004")
    elif confidence != "UNAVAILABLE" or value["permanently_incomplete_catalog"]:
        raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-004")
    for event_id, log_rate in value["event_log_rates"].items():
        event = catalog.events.get(event_id)
        if event is None or event.origin_state_id != state_id or event.log_rate != log_rate:
            raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-005")
    return DiscoveryRecord(state_id, config_digest, value["relevance_rate_min"], value["alpha"],
                           None if value["alpha_calibration"] is None else deep_freeze(value["alpha_calibration"]),
                           value["attempts"], value["successes"], dict(value["failures_by_status"]),
                           value["duplicates"], value["consecutive_redundant_successes"],
                           dict(value["event_log_rates"]), confidence, value["evaluations"],
                           value["stopping_state"], value["permanently_incomplete_catalog"])


def _restore_rate_snapshot(value: Any, catalog: Catalog, origin_id: str) -> dict[str, Any]:
    if (not isinstance(value, dict) or set(value) != {"payload", "payload_sha256"}
            or value["payload_sha256"] != digest(value["payload"])):
        raise _failure("CHECKPOINT_CORRUPT", "E2-KMC-002")
    payload = value["payload"]
    required = {"destination_state_ids", "event_ids", "log_rates", "lost_rate_log_upper_bound",
                "origin_state_id", "rates", "schema", "total_rate_per_s"}
    if (not isinstance(payload, dict) or set(payload) != required
            or payload["schema"] != "spark-atomistic-rate-table-snapshot/1"
            or payload["origin_state_id"] != origin_id):
        raise _failure("CHECKPOINT_CORRUPT", "E2-KMC-002")
    ids, destinations = payload["event_ids"], payload["destination_state_ids"]
    logs, rates = payload["log_rates"], payload["rates"]
    if (not isinstance(ids, list) or not ids or ids != sorted(set(ids))
            or not isinstance(destinations, list) or not isinstance(logs, list)
            or not isinstance(rates, list)
            or not len(ids) == len(destinations) == len(logs) == len(rates)):
        raise _failure("CHECKPOINT_CORRUPT", "E2-KMC-002")
    total = 0.0
    correction = 0.0
    for index, event_id in enumerate(ids):
        event = catalog.events.get(event_id)
        rate = _finite(rates[index])
        if (event is None or event.origin_state_id != origin_id or not event.selectable
                or destinations[index] != event.destination_state_id
                or logs[index] != event.log_rate or rate != math.exp(event.log_rate)):
            raise _failure("CHECKPOINT_CORRUPT", "E2-KMC-002")
        updated = total + rate
        correction += ((total - updated) + rate) if abs(total) >= abs(rate) else ((rate - updated) + total)
        total = updated
    if payload["total_rate_per_s"] != total + correction:
        raise _failure("CHECKPOINT_CORRUPT", "E2-KMC-002")
    if payload["lost_rate_log_upper_bound"] is not None:
        _finite(payload["lost_rate_log_upper_bound"])
    return payload


def validate_checkpoint_payload(value: Any, *, expected_config_digest: str,
                                expected_model_digest: str,
                                expected_tolerance_digest: str,
                                expected_identity_digest: str,
                                kinetics: Mapping[str, Any],
                                saddle_config: Mapping[str, Any],
                                relaxation_config: Mapping[str, Any],
                                discovery_config: Mapping[str, Any],
                                resource_config: Mapping[str, Any],
                                maximum_events: int) -> dict[str, Any]:
    required = {"basin", "catalog", "checkpoint_sequence", "current_state", "digests",
                "discovery_statistics", "flags", "initial_state", "log_sequence",
                "resources", "rng", "schema", "simulation_time_s", "step_index", "trajectory"}
    if not isinstance(value, dict) or set(value) != required or value["schema"] != CHECKPOINT_SCHEMA:
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-002")
    expected_digests = {"config": expected_config_digest, "model": expected_model_digest,
                        "schema": SCHEMA_SHA256, "tolerances": expected_tolerance_digest}
    if value["digests"] != expected_digests:
        # E2-CKPT-003 fixes only the KEY SET of `digests` ("has exactly `config`, `model`,
        # `schema`, and `tolerances`"), which a well-formed but wrong digest still satisfies.
        # A mismatched value is a run-contract mismatch, caught by E2-CKPT-007 step 3
        # ("schema/config/model/tolerance digests"). Cited as E2-CKPT-003 until 2026-08-11;
        # this was the last cross-language divergence in the mandatory fixture corpus.
        raise _failure("CHECKPOINT_INCOMPATIBLE", "E2-CKPT-007")
    if value["basin"] != {"enabled": False, "reason": "v1-disabled"}:
        raise _failure("CHECKPOINT_CORRUPT", "E2-BASIN-001")
    catalog_raw = value["catalog"]
    if (not isinstance(catalog_raw, dict)
            or set(catalog_raw) != {"digest", "events", "multiplicity", "schema", "states"}
            or catalog_raw["schema"] != "spark-atomistic-catalog/2"):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-006")
    catalog_payload = {key: catalog_raw[key] for key in ("events", "multiplicity", "schema", "states")}
    if catalog_raw["digest"] != digest(catalog_payload):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-006")
    if (not isinstance(catalog_raw["states"], dict) or not isinstance(catalog_raw["events"], dict)
            or len(catalog_raw["events"]) > maximum_events):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-006")
    catalog = Catalog(expected_model_digest, expected_config_digest,
                      expected_tolerance_digest, maximum_events, expected_identity_digest)
    for state_id in sorted(catalog_raw["states"]):
        state = validate_state_record(catalog_raw["states"][state_id], checkpoint=True)
        if (state.state_id != state_id or state.calculator_model_digest != expected_model_digest
                or state.force_tolerance != relaxation_config["force_tolerance"]):
            raise _failure("CHECKPOINT_INCOMPATIBLE", "E2-CKPT-007")
        if catalog.resolve_state(state, kinetics) is not state:
            raise _failure("CHECKPOINT_CORRUPT", "E2-ID-005")
        catalog.states[state_id] = state
    for event_id in sorted(catalog_raw["events"]):
        event = _restore_event(catalog_raw["events"][event_id], catalog.states,
                               expected_model_digest, expected_tolerance_digest,
                               expected_identity_digest, kinetics, saddle_config)
        if event.event_id != event_id:
            raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-001")
        catalog.events[event_id] = event
    for event in catalog.events.values():
        reverse = catalog.events.get(event.reverse_event_id)
        if (reverse is None or reverse.reverse_event_id != event.event_id
                or reverse.pair_id != event.pair_id
                or reverse.origin_state_id != event.destination_state_id
                or reverse.destination_state_id != event.origin_state_id
                or reverse.barrier != event.reverse_barrier
                or reverse.reverse_barrier != event.barrier
                or reverse.log_rate != event.reverse_log_rate
                or reverse.reverse_log_rate != event.log_rate
                or reverse.detailed_balance_residual != -event.detailed_balance_residual
                or reverse.saddle_energy != event.saddle_energy
                or reverse.curvature != event.curvature
                or reverse.orthogonal_curvatures != event.orthogonal_curvatures
                or reverse.evidence_level != event.evidence_level
                or reverse.evaluation_count != event.evaluation_count
                or reverse.search_id != event.search_id
                or reverse.termination_reason != event.termination_reason
                or reverse.discovery_provenance != event.discovery_provenance):
            raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-001")
        mapping = event.validation["destination_match"]["atom_mapping"]
        expected_positions: list[Any] = [None] * len(mapping)
        expected_forces: list[Any] = [None] * len(mapping)
        expected_mode: list[Any] = [None] * len(mapping)
        for source, destination in enumerate(mapping):
            expected_positions[destination] = event.saddle_positions[source]
            expected_forces[destination] = event.saddle_forces[source]
            expected_mode[destination] = tuple(-component for component in event.unstable_direction[source])
        if (tuple(expected_positions) != reverse.saddle_positions
                or tuple(expected_forces) != reverse.saddle_forces
                or tuple(expected_mode) != reverse.unstable_direction):
            raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-001")
    multiplicity = catalog_raw["multiplicity"]
    if (not isinstance(multiplicity, dict) or set(multiplicity) != set(catalog.events)
            or any(type(count) is not int or count < 1 for count in multiplicity.values())
            or any(multiplicity[event.event_id] != multiplicity[event.reverse_event_id]
                   for event in catalog.events.values())):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-006")
    catalog.multiplicity = dict(multiplicity)
    statistics = value["discovery_statistics"]
    if not isinstance(statistics, dict):
        raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-004")
    for state_id in sorted(statistics):
        if state_id not in catalog.states:
            raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-004")
        catalog.discovery[state_id] = _restore_discovery(
            statistics[state_id], state_id, expected_config_digest,
            discovery_config, catalog)
    for event in catalog.events.values():
        provenance = event.discovery_provenance
        matching_sources = []
        for source_id in (event.origin_state_id, event.destination_state_id):
            expected_search_id = "search:" + digest({
                "run_seed": kinetics["run_seed"],
                "search_class": provenance["search_class"],
                "search_index": provenance["search_index"],
                "state_id": source_id})
            if expected_search_id == provenance["search_id"]:
                matching_sources.append(source_id)
        if len(matching_sources) != 1:
            raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-002")
        source_record = catalog.discovery.get(matching_sources[0])
        if source_record is None or provenance["search_index"] >= source_record.attempts:
            raise _failure("CHECKPOINT_CORRUPT", "E2-DISC-005")
    resources = value["resources"]
    resource_keys = {"calculator_evaluations", "catalog_events", "output_bytes",
                     "resident_memory_bytes", "retry_history", "saddle_attempts_by_state",
                     "wall_elapsed_s"}
    if (not isinstance(resources, dict) or set(resources) != resource_keys
            or resources["catalog_events"] != len(catalog.events)
            or resources["retry_history"] != []
            or any(type(resources[key]) is not int or resources[key] < 0 for key in
                   ("calculator_evaluations", "catalog_events", "output_bytes",
                    "resident_memory_bytes"))
            or resources["calculator_evaluations"] > resource_config["total_calculator_evaluations"]
            or resources["output_bytes"] > resource_config["output_bytes"]
            or resources["resident_memory_bytes"] > resource_config["resident_memory_bytes"]
            or not isinstance(resources["saddle_attempts_by_state"], dict)
            or any(not isinstance(key, str) or not key or type(count) is not int
                   or not 0 <= count <= resource_config["saddle_attempts_per_state"]
                   for key, count in resources["saddle_attempts_by_state"].items())):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-005")
    _finite(resources["wall_elapsed_s"], nonnegative=True)
    expected_attempts = {state_id: record.attempts for state_id, record in catalog.discovery.items()
                         if record.attempts}
    if resources["saddle_attempts_by_state"] != expected_attempts:
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-005")
    flags = value["flags"]
    flag_keys = {"cancelled", "complete", "incomplete_catalog", "last_status", "resource_limited"}
    if (not isinstance(flags, dict) or set(flags) != flag_keys
            or any(type(flags[key]) is not bool for key in
                   ("cancelled", "complete", "incomplete_catalog", "resource_limited"))
            or flags["last_status"] not in STATUSES):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-004")
    actual_incomplete = any(record.stopping_state == "INCOMPLETE"
                            for record in catalog.discovery.values())
    if flags["incomplete_catalog"] != actual_incomplete:
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-004")
    initial = validate_state_record(value["initial_state"], checkpoint=True)
    current = validate_state_record(value["current_state"], checkpoint=True)
    if initial.state_id not in catalog.states or current.state_id not in catalog.states:
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-007")
    if (initial.record() != catalog.states[initial.state_id].record()
            or current.record() != catalog.states[current.state_id].record()):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-007")
    for key in ("checkpoint_sequence", "log_sequence", "step_index"):
        if type(value[key]) is not int or value[key] < 0:
            raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-008")
    if value["checkpoint_sequence"] < 1:
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-008")
    trajectory = value["trajectory"]
    if (not isinstance(trajectory, list) or len(trajectory) != value["step_index"]
            or value["log_sequence"] != len(trajectory)):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-008")
    replay = derive_trajectory_stream(kinetics["run_seed"])
    cursor = initial.state_id
    cumulative_time = 0.0
    prior_checkpoint_sequence = -1
    step_keys = {"checkpoint_sequence", "log_sequence", "post_state_id", "pre_state_id",
                 "rate_table_snapshot", "selected_event_id", "selected_rate_per_s",
                 "selection_uniform", "step_index", "time_increment_s", "time_uniform",
                 "total_rate_per_s"}
    for index, step in enumerate(trajectory, 1):
        if (not isinstance(step, dict) or set(step) != step_keys
                or step["step_index"] != index or step["log_sequence"] != index
                or step["pre_state_id"] != cursor
                or type(step["checkpoint_sequence"]) is not int
                or not prior_checkpoint_sequence <= step["checkpoint_sequence"]
                       < value["checkpoint_sequence"]):
            raise _failure("CHECKPOINT_CORRUPT", "E2-KMC-004")
        prior_checkpoint_sequence = step["checkpoint_sequence"]
        event = catalog.events.get(step["selected_event_id"])
        if (event is None or event.origin_state_id != cursor
                or event.destination_state_id != step["post_state_id"]):
            raise _failure("CHECKPOINT_CORRUPT", "E2-KMC-004")
        table = _restore_rate_snapshot(step["rate_table_snapshot"], catalog, cursor)
        selection_uniform = replay.uniform()
        time_uniform = replay.uniform()
        if (step["selection_uniform"] != selection_uniform or step["time_uniform"] != time_uniform):
            raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-008")
        threshold = selection_uniform * table["total_rate_per_s"]
        cumulative = 0.0
        selected = len(table["rates"]) - 1
        for rate_index, rate in enumerate(table["rates"]):
            cumulative += rate
            if cumulative > threshold:
                selected = rate_index
                break
        expected_time = -math.log(time_uniform) / table["total_rate_per_s"]
        if (table["event_ids"][selected] != event.event_id
                or step["selected_rate_per_s"] != table["rates"][selected]
                or step["total_rate_per_s"] != table["total_rate_per_s"]
                or not math.isclose(step["time_increment_s"], expected_time,
                                    rel_tol=5e-15, abs_tol=1e-18)):
            raise _failure("CHECKPOINT_CORRUPT", "E2-KMC-005")
        cumulative_time += step["time_increment_s"]
        cursor = event.destination_state_id
    simulation_time = _finite(value["simulation_time_s"], nonnegative=True)
    if (cursor != current.state_id
            or not math.isclose(simulation_time, cumulative_time, rel_tol=5e-15, abs_tol=1e-18)):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-008")
    rng = value["rng"]
    if (not isinstance(rng, dict) or set(rng) != {"run_seed", "substream_map", "trajectory"}
            or rng["run_seed"] != kinetics["run_seed"]
            or not isinstance(rng["substream_map"], dict)):
        raise _failure("CHECKPOINT_INCOMPATIBLE", "E2-CKPT-003")
    trajectory_rng = PhiloxStream.restore(rng["trajectory"])
    if (trajectory_rng.checkpoint() != replay.checkpoint()
            or trajectory_rng.consumed_uniforms != 2 * value["step_index"]):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-008")
    substreams = {key: PhiloxStream.restore(item)
                  for key, item in rng["substream_map"].items()
                  if isinstance(key, str) and key}
    expected_keys: set[str] = set()
    for state_id, record in catalog.discovery.items():
        for search_index in range(record.attempts):
            class_stream = derive_saddle_stream(kinetics["run_seed"], state_id,
                                                 "class-selection", search_index)
            uniform = class_stream.uniform()
            cumulative = 0.0
            selected_class = discovery_config["classes"][-1]
            for entry in discovery_config["classes"]:
                cumulative += entry["probability"]
                if cumulative > uniform:
                    selected_class = entry
                    break
            search_id = "search:" + digest({"run_seed": kinetics["run_seed"],
                                            "search_class": selected_class["name"],
                                            "search_index": search_index,
                                            "state_id": state_id})
            class_key = "class-selection:" + search_id
            expected_keys.update({class_key, search_id})
            if (class_key not in substreams
                    or substreams[class_key].checkpoint() != class_stream.checkpoint()):
                raise _failure("CHECKPOINT_CORRUPT", "E2-RNG-005")
            expected_search = derive_saddle_stream(kinetics["run_seed"], state_id,
                                                   selected_class["name"], search_index)
            actual_search = substreams.get(search_id)
            if (actual_search is None or actual_search.key != expected_search.key
                    or actual_search.initial_counter != expected_search.initial_counter):
                raise _failure("CHECKPOINT_CORRUPT", "E2-RNG-004")
    if set(substreams) != expected_keys:
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-007")
    for event in catalog.events.values():
        search_id = event.discovery_provenance["search_id"]
        if (search_id not in substreams
                or event.discovery_provenance["rng_substream_digest"]
                   != digest(substreams[search_id].checkpoint())):
            raise _failure("CHECKPOINT_CORRUPT", "E2-EVENT-004")
    if value["step_index"] > kinetics["maximum_steps"]:
        raise _failure("CHECKPOINT_INCOMPATIBLE", "E2-CKPT-008")
    expected_complete = (value["step_index"] == kinetics["maximum_steps"]
                         or (flags["last_status"] == "NO_ENABLED_EVENT"
                             and kinetics["absorbing_ok"]))
    if (flags["complete"] != expected_complete
            or flags["cancelled"] != (flags["last_status"] == "CANCELLED")
            or flags["resource_limited"] != (flags["last_status"] == "RESOURCE_LIMIT")):
        raise _failure("CHECKPOINT_CORRUPT", "E2-CKPT-004")
    return {"catalog": catalog, "current_state_id": current.state_id,
            "initial_state_id": initial.state_id, "simulation_time": simulation_time,
            "step_index": value["step_index"], "log_sequence": value["log_sequence"],
            "checkpoint_sequence": value["checkpoint_sequence"],
            "trajectory": [deep_freeze(item) for item in trajectory],
            "trajectory_rng": trajectory_rng, "substreams": substreams,
            "resources": resources, "flags": dict(flags)}
