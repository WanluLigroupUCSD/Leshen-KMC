# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Errata-2 exact status, context, message, severity, and exit contract."""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Mapping


_STATUS_ROWS = {
    "OK": ("success", "transaction committed"),
    "DISCOVERY_CONVERGED_HEURISTIC": ("success-with-qualification", "heuristic discovery criterion passed"),
    "DUPLICATE_EVENT": ("candidate reject", "duplicate event rejected"),
    "SADDLE_NOT_FOUND": ("candidate reject", "saddle not found"),
    "INVALID_SADDLE": ("candidate reject", "saddle validation failed"),
    "SADDLE_WRONG_BASIN": ("candidate reject", "neither endpoint matches origin"),
    "ENDPOINT_COLLAPSED": ("candidate reject", "both endpoints match origin"),
    "ENVIRONMENT_AMBIGUOUS": ("recoverable", "environment identity ambiguous"),
    "BASIN_DISABLED": ("recoverable", "basin acceleration disabled"),
    "DISCOVERY_INCOMPLETE": ("pause/qualified", "discovery budget exhausted"),
    "RELAX_NOT_CONVERGED": ("transaction fail", "relaxation not converged"),
    "EVENT_APPLICATION_FAILED": ("transaction fail", "event application failed"),
    "CALCULATOR_FAILURE": ("transaction fail", "calculator callback failed"),
    "NONFINITE_RESULT": ("fatal", "nonfinite value rejected"),
    "INVALID_INPUT": ("fatal", "input invalid"),
    "SCHEMA_UNSUPPORTED": ("fatal", "schema unsupported"),
    "INVALID_STATE": ("fatal", "state invalid"),
    "RATE_INVALID": ("fatal in strict mode", "rate invalid"),
    "DETAILED_BALANCE_VIOLATION": ("fatal in strict mode", "detailed balance violated"),
    "CATALOG_CONFLICT": ("fatal", "catalog conflict"),
    "CATALOG_INCOMPATIBLE": ("fatal", "catalog incompatible"),
    "ATOM_COUNT_CHANGE_UNSUPPORTED": ("fatal", "atom count change unsupported"),
    "NO_ENABLED_EVENT": ("terminal-success if requested, else fatal", "no enabled event"),
    "RESOURCE_LIMIT": ("pause", "resource limit reached"),
    "OUTPUT_EXISTS": ("fatal", "output exists"),
    "CHECKPOINT_CORRUPT": ("fatal", "checkpoint corrupt"),
    "CHECKPOINT_INCOMPATIBLE": ("fatal", "checkpoint incompatible"),
    "CANCELLED": ("pause", "run cancelled"),
    "INTERNAL_ERROR": ("fatal", "internal error"),
}
STATUSES = frozenset(_STATUS_ROWS)
SEVERITY = MappingProxyType({key: row[0] for key, row in _STATUS_ROWS.items()})
MESSAGE = MappingProxyType({key: row[1] for key, row in _STATUS_ROWS.items()})
NONTERMINAL = frozenset({
    "DUPLICATE_EVENT", "SADDLE_NOT_FOUND", "INVALID_SADDLE", "SADDLE_WRONG_BASIN",
    "ENDPOINT_COLLAPSED", "ENVIRONMENT_AMBIGUOUS", "BASIN_DISABLED",
})


def exit_code(status: str, *, exploratory: bool = False, absorbing_ok: bool = False,
              checkpoint_lost: bool = False) -> int:
    if checkpoint_lost:
        return 74
    if status in {"OK", "DISCOVERY_CONVERGED_HEURISTIC"}:
        return 0
    if status == "DISCOVERY_INCOMPLETE":
        return 0 if exploratory else 75
    if status == "NO_ENABLED_EVENT":
        return 0 if absorbing_ok else 65
    if status in {"INVALID_INPUT", "SCHEMA_UNSUPPORTED", "OUTPUT_EXISTS"}:
        return 64
    if status in {"CHECKPOINT_CORRUPT", "CHECKPOINT_INCOMPATIBLE"}:
        return 74
    if status == "CALCULATOR_FAILURE":
        return 69
    if status == "INTERNAL_ERROR":
        return 70
    if status in {"RESOURCE_LIMIT", "CANCELLED", "DISCOVERY_INCOMPLETE"}:
        return 75
    return 65


def _normalize_context(value: Mapping[str, Any]) -> Mapping[str, Any]:
    fixed = {"component", "details", "requirement_id", "retryable",
             "search_or_event_id", "state_id"}
    details = dict(value.get("details", {})) if isinstance(value.get("details", {}), Mapping) else {}
    for key, item in value.items():
        if key not in fixed:
            details[key] = item
    context = {
        "component": str(value.get("component", "internal")),
        "details": details,
        "requirement_id": str(value.get("requirement_id", "E2-STATUS-001")),
        "retryable": value.get("retryable", False) is True,
        "search_or_event_id": value.get("search_or_event_id"),
        "state_id": value.get("state_id"),
    }
    return MappingProxyType(context)


@dataclass(frozen=True, slots=True)
class Outcome:
    status: str
    message: str
    context: Mapping[str, Any]
    value: Any = None
    causal_status: str | None = None

    def __post_init__(self) -> None:
        if self.status not in STATUSES:
            raise ValueError("unknown internal status")
        if self.causal_status is not None and self.causal_status not in STATUSES:
            raise ValueError("unknown causal status")
        object.__setattr__(self, "message", MESSAGE[self.status])
        object.__setattr__(self, "context", _normalize_context(self.context))

    @property
    def ok(self) -> bool:
        return self.status in {"OK", "DISCOVERY_CONVERGED_HEURISTIC"}

    def record(self, *, terminal_exit_code: int | None = None) -> dict[str, Any]:
        return {
            "causal_status": self.causal_status,
            "context": {"component": self.context["component"],
                        "details": dict(self.context["details"]),
                        "requirement_id": self.context["requirement_id"],
                        "retryable": self.context["retryable"],
                        "search_or_event_id": self.context["search_or_event_id"],
                        "state_id": self.context["state_id"]},
            "exit_code": terminal_exit_code,
            "message": self.message,
            "severity": SEVERITY[self.status],
            "status": self.status,
        }


class DomainFailure(Exception):
    def __init__(self, status: str, message: str, *, component: str,
                 requirement: str, retryable: bool = False,
                 state_id: str | None = None, object_id: str | None = None,
                 causal_status: str | None = None, calculator_evaluations: int = 0,
                 iterations: int = 0, termination_reason: str | None = None,
                 details: Mapping[str, Any] | None = None) -> None:
        super().__init__(MESSAGE.get(status, message))
        detail_map = dict(details or {})
        if calculator_evaluations or iterations or termination_reason is not None:
            detail_map.update({"calculator_evaluations": calculator_evaluations,
                               "iterations": iterations,
                               "termination_reason": termination_reason or status})
        self.outcome = Outcome(
            status, message,
            {"component": component, "details": detail_map,
             "requirement_id": requirement, "retryable": retryable,
             "search_or_event_id": object_id, "state_id": state_id},
            causal_status=causal_status)
