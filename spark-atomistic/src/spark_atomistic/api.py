# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Three-operation Errata-2 canonical JSON API."""

from __future__ import annotations

from typing import Any, Mapping

from .canonical import canonical_text, deep_freeze, deep_thaw, digest, parse_json
from .engine import ReferenceEngine
from .errors import DomainFailure, Outcome, exit_code
from .model import EXPECTED_SCHEMA_SHA256, IR, SCHEMA_SHA256, validate_model


_CAPABILITY_VALUE = deep_freeze({
    "api": "spark-atomistic-json/1",
    "basin_acceleration": "disabled",
    "conformance": "unvalidated",
    "features": {"common_prefactor": True, "deterministic_checkpoint": True,
                 "fixed_composition_off_lattice": True, "harmonic_tst": False,
                 "local_environment_generic_reuse": False, "serial_kmc": True,
                 "variable_composition": False},
    "ir": IR,
    "operations": ["capabilities", "validate", "run"],
    "production": False, "release": False, "validated": False,
})


def capabilities() -> Mapping[str, Any]:
    return _CAPABILITY_VALUE


def _outcome(status: str, requirement: str, *, component: str = "api",
             state_id: str | None = None, object_id: str | None = None,
             details: Mapping[str, Any] | None = None,
             causal_status: str | None = None, value: Any = None) -> Outcome:
    return Outcome(status, "", {"component": component, "details": dict(details or {}),
                                "requirement_id": requirement, "retryable": False,
                                "search_or_event_id": object_id, "state_id": state_id},
                   value=value, causal_status=causal_status)


def _response(operation: str, outcome: Outcome, *, exploratory: bool = False,
              absorbing_ok: bool = False, checkpoint_lost: bool = False,
              value: Any = None) -> str:
    final_exit = exit_code(
        outcome.status, exploratory=exploratory, absorbing_ok=absorbing_ok,
        checkpoint_lost=checkpoint_lost)
    record = outcome.record(terminal_exit_code=final_exit)
    record["operation"] = operation
    record["value"] = value if final_exit == 0 else None
    return canonical_text(record)


def _strict_request(value: Any, operation: str) -> dict[str, Any]:
    required = ({"operation"} if operation == "capabilities" else
                {"operation", "model"} if operation == "validate" else
                {"allow_unvalidated", "extension", "model", "operation"})
    if not isinstance(value, dict) or set(value) != required:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="api",
                            requirement="E2-API-001",
                            details={"required_keys": sorted(required)})
    return value


def process_atomistic_json(request_json: str | bytes | bytearray, *,
                           source_path: str | None = None) -> str:
    operation = "unknown"
    engine: ReferenceEngine | None = None
    exploratory = False
    absorbing_ok = False
    try:
        request = parse_json(request_json)
        if isinstance(request, dict) and isinstance(request.get("operation"), str):
            operation = request["operation"]
        if operation not in {"capabilities", "validate", "run"}:
            raise DomainFailure("INVALID_INPUT", "input invalid", component="api",
                                requirement="E2-API-001")
        request = _strict_request(request, operation)
        if operation == "capabilities":
            return _response(operation, _outcome("OK", "E2-API-001"),
                             value=deep_thaw(_CAPABILITY_VALUE))
        if operation == "run" and request["allow_unvalidated"] is not True:
            raise DomainFailure("INVALID_INPUT", "input invalid", component="api",
                                requirement="E2-API-004")
        if operation == "run" and not isinstance(request["extension"], dict):
            raise DomainFailure("INVALID_INPUT", "input invalid", component="api",
                                requirement="E2-API-003")
        model = validate_model(request["model"], source_path=source_path)
        behavior = deep_thaw(model)
        behavior.pop("metadata", None)
        config_digest = digest(behavior)
        if SCHEMA_SHA256 != EXPECTED_SCHEMA_SHA256:
            raise DomainFailure("INTERNAL_ERROR", "internal error", component="schema",
                                requirement="E2-CAN-006")
        if operation == "validate":
            value = {"config_digest": config_digest, "ir": IR,
                     "schema_digest": SCHEMA_SHA256}
            return _response(operation, _outcome("OK", "E2-API-002"), value=value)
        exploratory = model["discovery"]["mode"] == "exploratory"
        absorbing_ok = model["kinetics"]["absorbing_ok"]
        engine = ReferenceEngine(model, extension=request["extension"])
        result = engine.run()
        summary = engine.public_summary()
        return _response(operation, result, exploratory=exploratory,
                         absorbing_ok=absorbing_ok, value=summary)
    except KeyboardInterrupt:
        outcome = _outcome("CANCELLED", "E2-STATUS-002",
                           state_id=engine.current_state_id if engine else None)
    except DomainFailure as exc:
        outcome = exc.outcome
    except Exception:
        outcome = _outcome("INTERNAL_ERROR", "E2-STATUS-002")
    checkpoint_failure: Outcome | None = None
    if engine is not None and engine.current_state_id:
        engine.last_status = outcome.status
        engine.cancelled = outcome.status == "CANCELLED"
        engine.resource_limited = outcome.status == "RESOURCE_LIMIT"
        try:
            engine.write_checkpoint()
        except DomainFailure as exc:
            checkpoint_failure = exc.outcome
        except Exception:
            checkpoint_failure = _outcome("INTERNAL_ERROR", "E2-CKPT-009",
                                          component="checkpoint", causal_status=outcome.status)
    if checkpoint_failure is not None:
        details = dict(outcome.context["details"])
        details["checkpoint_failure"] = checkpoint_failure.record()
        outcome = _outcome(outcome.status, str(outcome.context["requirement_id"]),
                           component=str(outcome.context["component"]),
                           state_id=outcome.context["state_id"],
                           object_id=outcome.context["search_or_event_id"],
                           details=details, causal_status=outcome.causal_status)
    checkpoint_io = (outcome.status == "INTERNAL_ERROR"
                     and outcome.context["component"] == "checkpoint")
    valid_last = engine is not None and engine.has_valid_last_checkpoint()
    return _response(operation, outcome, exploratory=exploratory,
                     absorbing_ok=absorbing_ok,
                     checkpoint_lost=(checkpoint_failure is not None or checkpoint_io)
                     and not valid_last)


def validate_atomistic_model_json(request_json: str | bytes | bytearray) -> str:
    return process_atomistic_json(request_json)


def run_atomistic_reference_unvalidated_json(
        request_json: str | bytes | bytearray, *, allow_unvalidated: bool = False) -> str:
    # The normative trust gate is inside the JSON request. The legacy keyword
    # cannot authorize a request and is intentionally ignored.
    _ = allow_unvalidated
    return process_atomistic_json(request_json)
