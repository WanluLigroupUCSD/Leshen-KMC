# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Evaluation reservation and persistent resource counters."""

from __future__ import annotations

from dataclasses import dataclass, field
import resource
import sys
import time
from typing import Any, Mapping

from .errors import DomainFailure


@dataclass(slots=True)
class ResourceLedger:
    total_limit: int
    wall_time_limit: float
    resident_memory_limit: int
    output_limit: int
    saddle_attempt_limit: int
    started_monotonic: float
    calculator_reserved: int = 0
    calculator_completed: int = 0
    calculator_failed: int = 0
    retry_count: int = 0
    output_bytes_reserved: int = 0
    output_bytes_written: int = 0
    per_state_saddle_attempts: dict[str, int] = field(default_factory=dict)
    retry_history: list[dict[str, Any]] = field(default_factory=list)
    attempt_history: list[dict[str, Any]] = field(default_factory=list)
    accumulated_wall_s: float = 0.0

    @classmethod
    def start(cls, config: Mapping[str, Any]) -> "ResourceLedger":
        return cls(config["total_calculator_evaluations"], config["wall_time_s"],
                   config["resident_memory_bytes"], config["output_bytes"],
                   config["saddle_attempts_per_state"], time.monotonic())

    def check_wall_time(self) -> None:
        if self.accumulated_wall_s + time.monotonic() - self.started_monotonic >= self.wall_time_limit:
            raise DomainFailure("RESOURCE_LIMIT", "wall-time limit reached",
                                component="resources", requirement="RES-002")
        maximum_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        resident_bytes = maximum_rss if sys.platform == "darwin" else maximum_rss * 1024
        if resident_bytes >= self.resident_memory_limit:
            raise DomainFailure("RESOURCE_LIMIT", "resident-memory limit reached",
                                component="resources", requirement="RES-002")

    def reserve_evaluation(self) -> int:
        self.check_wall_time()
        if self.calculator_reserved >= self.total_limit:
            raise DomainFailure("RESOURCE_LIMIT", "calculator evaluation limit reached",
                                component="resources", requirement="RES-002")
        self.calculator_reserved += 1
        return self.calculator_reserved

    def complete(self, success: bool) -> None:
        self.calculator_completed += 1
        if not success:
            self.calculator_failed += 1

    def reserve_saddle_attempt(self, state_id: str) -> None:
        count = self.per_state_saddle_attempts.get(state_id, 0)
        if count >= self.saddle_attempt_limit:
            raise DomainFailure("RESOURCE_LIMIT", "per-state saddle-attempt limit reached",
                                component="resources", requirement="RES-002", state_id=state_id)
        self.per_state_saddle_attempts[state_id] = count + 1

    def reserve_output(self, byte_count: int) -> None:
        if byte_count <= 0 or self.output_bytes_reserved + byte_count > self.output_limit:
            raise DomainFailure("RESOURCE_LIMIT", "cumulative output byte limit reached",
                                component="resources", requirement="RES-002")
        self.output_bytes_reserved += byte_count

    def complete_output(self, byte_count: int) -> None:
        if byte_count <= 0 or self.output_bytes_written + byte_count > self.output_bytes_reserved:
            raise DomainFailure("INTERNAL_ERROR", "output-byte accounting invariant failed",
                                component="resources", requirement="RES-004")
        self.output_bytes_written += byte_count

    def record_attempt(self, kind: str, attempt_id: str, reserved_before: int,
                       final_status: str, termination_reason: str) -> None:
        reserved_after = self.calculator_reserved
        if kind not in {"callback", "solver", "endpoint"} or not attempt_id:
            raise DomainFailure("INTERNAL_ERROR", "invalid attempt audit record",
                                component="resources", requirement="RES-004")
        self.attempt_history.append({
            "sequence": len(self.attempt_history) + 1, "kind": kind,
            "attempt_id": attempt_id, "calculator_reserved_before": reserved_before,
            "calculator_reserved_after": reserved_after,
            "calculator_evaluations": reserved_after - reserved_before,
            "final_status": final_status, "termination_reason": termination_reason,
        })

    def wall_elapsed(self) -> float:
        return self.accumulated_wall_s + time.monotonic() - self.started_monotonic

    @staticmethod
    def resident_memory() -> int:
        maximum_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return maximum_rss if sys.platform == "darwin" else maximum_rss * 1024

    def checkpoint(self, *, wall_elapsed: float | None = None,
                   projected_reserved: int | None = None,
                   projected_written: int | None = None,
                   catalog_events: int = 0,
                   resident_memory_bytes: int | None = None) -> dict[str, Any]:
        output_bytes = (self.output_bytes_written if projected_written is None
                        else projected_written)
        return {
            "calculator_evaluations": self.calculator_reserved,
            "catalog_events": catalog_events,
            "output_bytes": output_bytes,
            "resident_memory_bytes": (self.resident_memory() if resident_memory_bytes is None
                                      else resident_memory_bytes),
            "retry_history": list(self.retry_history),
            "saddle_attempts_by_state": dict(sorted(self.per_state_saddle_attempts.items())),
            "wall_elapsed_s": self.wall_elapsed() if wall_elapsed is None else wall_elapsed,
        }

    def restore_counts(self, record: Any) -> None:
        required = {"calculator_evaluations", "catalog_events", "output_bytes",
                    "resident_memory_bytes", "retry_history", "saddle_attempts_by_state",
                    "wall_elapsed_s"}
        integer_keys = {"calculator_evaluations", "catalog_events", "output_bytes",
                        "resident_memory_bytes"}
        if (not isinstance(record, dict) or set(record) != required
                or any(type(record[key]) is not int or record[key] < 0 for key in integer_keys)
                or isinstance(record["wall_elapsed_s"], bool)
                or not isinstance(record["wall_elapsed_s"], (int, float))
                or not 0.0 <= float(record["wall_elapsed_s"]) < float("inf")):
            raise DomainFailure("CHECKPOINT_CORRUPT", "invalid resource counters",
                                component="checkpoint", requirement="CKPT-004")
        if (record["calculator_evaluations"] > self.total_limit
                or record["catalog_events"] < 0
                or record["output_bytes"] > self.output_limit
                or record["resident_memory_bytes"] > self.resident_memory_limit):
            raise DomainFailure("CHECKPOINT_INCOMPATIBLE", "checkpoint exceeds resource configuration",
                                component="checkpoint", requirement="E2-CKPT-005")
        if (record["retry_history"] != []
                or not isinstance(record["saddle_attempts_by_state"], dict)
                or any(not isinstance(key, str) or not key or type(item) is not int
                       or not 0 <= item <= self.saddle_attempt_limit
                       for key, item in record["saddle_attempts_by_state"].items())):
            raise DomainFailure("CHECKPOINT_CORRUPT", "resource history/accounting mismatch",
                                component="checkpoint", requirement="E2-CKPT-005")
        self.calculator_reserved = record["calculator_evaluations"]
        self.calculator_completed = record["calculator_evaluations"]
        self.calculator_failed = 0
        self.retry_count = 0
        self.retry_history = []
        self.attempt_history = []
        self.per_state_saddle_attempts = dict(record["saddle_attempts_by_state"])
        self.output_bytes_reserved = record["output_bytes"]
        self.output_bytes_written = record["output_bytes"]
        self.accumulated_wall_s = float(record["wall_elapsed_s"])
        self.started_monotonic = time.monotonic()
