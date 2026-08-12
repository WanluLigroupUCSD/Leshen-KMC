# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Process-isolated, timeout-enforced JSON-lines calculator transport."""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
import selectors
import signal
import subprocess
import time
from contextlib import contextmanager
from typing import Any, Mapping

from .canonical import canonical_bytes, digest, parse_json
from .errors import DomainFailure
from .geometry import Vector
from .resources import ResourceLedger


@dataclass(frozen=True, slots=True)
class Evaluation:
    energy: float
    forces: tuple[Vector, ...]
    evaluation_id: str
    deterministic: bool
    request_digest: str


class ProcessCalculator:
    """One fresh child process per request; timeout kills the entire child process."""

    def __init__(self, config: Mapping[str, Any], timeout: float, stdout_limit: int,
                 stderr_limit: int, ledger: ResourceLedger) -> None:
        self.command = tuple(config["command"])
        self.model_name = config["model_name"]
        self.model_version = config["model_version"]
        self.model_digest = config["model_digest"]
        self.declared_deterministic = config["deterministic"]
        self.timeout = timeout
        self.stdout_limit = stdout_limit
        self.stderr_limit = stderr_limit
        self.ledger = ledger
        self._scope: tuple[str, int, int] | None = None

    @contextmanager
    def evaluation_scope(self, scope_id: str, limit: int):
        if self._scope is not None:
            raise DomainFailure("INTERNAL_ERROR", "nested calculator budget scope",
                                component="resources", requirement="RES-001")
        self._scope = (scope_id, limit, 0)
        try:
            yield
        finally:
            self._scope = None

    def evaluate(self, state_request: Mapping[str, Any], *, component: str,
                 object_id: str | None = None) -> Evaluation:
        before = self.ledger.calculator_reserved
        completed_before = self.ledger.calculator_completed
        final_status = "INTERNAL_ERROR"
        termination_reason = "callback aborted unexpectedly"
        try:
            result = self._evaluate(state_request, component=component, object_id=object_id)
            final_status = "OK"
            termination_reason = "validated calculator response"
            return result
        except DomainFailure as exc:
            final_status = exc.outcome.status
            termination_reason = str(exc.outcome.context["details"].get(
                "termination_reason", exc.outcome.status))
            raise
        except BaseException:
            final_status = "CANCELLED"
            termination_reason = "callback interrupted"
            raise
        finally:
            if (self.ledger.calculator_reserved - before
                    > self.ledger.calculator_completed - completed_before):
                self.ledger.complete(False)
            self.ledger.record_attempt("callback", object_id or component, before,
                                       final_status, termination_reason)

    def _evaluate(self, state_request: Mapping[str, Any], *, component: str,
                  object_id: str | None = None) -> Evaluation:
        if self._scope is not None:
            scope_id, limit, count = self._scope
            if count >= limit:
                raise DomainFailure("RESOURCE_LIMIT", "scoped calculator evaluation limit reached",
                                    component="resources", requirement="RES-002", object_id=scope_id)
            self._scope = (scope_id, limit, count + 1)
        self.ledger.reserve_evaluation()
        payload = {
            "schema": "spark-atomistic-calculator-request/1",
            "state": dict(state_request),
            "properties": ["energy", "forces"],
            "units": {"energy": "eV", "forces": "eV/angstrom", "length": "angstrom"},
            "model_digest": self.model_digest,
        }
        before = digest(payload)
        payload["request_digest"] = before
        request_bytes = canonical_bytes(payload) + b"\n"
        try:
            returncode, stdout, _stderr = self._invoke(request_bytes, component, object_id)
            if returncode != 0:
                raise DomainFailure("CALCULATOR_FAILURE", "calculator process returned nonzero",
                                    component=component, requirement="CALC-006", object_id=object_id)
            if b"\n" in stdout.rstrip(b"\n"):
                raise DomainFailure("CALCULATOR_FAILURE", "calculator returned more than one JSON record",
                                    component=component, requirement="CALC-006", object_id=object_id)
            try:
                response = parse_json(stdout)
            except DomainFailure as exc:
                if exc.outcome.status == "NONFINITE_RESULT":
                    raise
                raise DomainFailure("CALCULATOR_FAILURE", "calculator returned malformed JSON",
                                    component=component, requirement="CALC-006",
                                    object_id=object_id, causal_status=exc.outcome.status) from exc
            result = self._validate_response(response, state_request, before, component, object_id)
            after_payload = {
                "schema": "spark-atomistic-calculator-request/1",
                "state": dict(state_request), "properties": ["energy", "forces"],
                "units": {"energy": "eV", "forces": "eV/angstrom", "length": "angstrom"},
                "model_digest": self.model_digest,
            }
            if digest(after_payload) != before:
                raise DomainFailure("CALCULATOR_FAILURE", "calculator request mutated during evaluation",
                                    component=component, requirement="CALC-003", object_id=object_id)
            self.ledger.complete(True)
            return result
        except DomainFailure:
            self.ledger.complete(False)
            raise
        except (OSError, ValueError, TypeError, MemoryError) as exc:
            self.ledger.complete(False)
            raise DomainFailure("CALCULATOR_FAILURE", "calculator process/response failure",
                                component=component, requirement="CALC-006", object_id=object_id) from exc

    @staticmethod
    def _kill_group_and_wait(process: subprocess.Popen[bytes]) -> None:
        deadline = time.monotonic() + 2.0
        while True:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                if process.poll() is None:
                    try:
                        process.kill()
                    except ProcessLookupError:
                        pass
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                raise DomainFailure("CALCULATOR_FAILURE", "calculator process group could not be reaped within bound",
                                    component="calculator", requirement="CALC-006")
            try:
                process.wait(timeout=remaining)
                return
            except subprocess.TimeoutExpired:
                try:
                    process.kill()
                except ProcessLookupError:
                    pass

    def _invoke(self, request: bytes, component: str,
                object_id: str | None) -> tuple[int, bytes, bytes]:
        process = subprocess.Popen(
            self.command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, bufsize=0, start_new_session=True,
        )
        if process.stdin is None or process.stdout is None or process.stderr is None:
            self._kill_group_and_wait(process)
            raise DomainFailure("CALCULATOR_FAILURE", "calculator pipes unavailable",
                                component=component, requirement="CALC-006", object_id=object_id)
        streams = selectors.DefaultSelector()
        stdout = bytearray()
        stderr = bytearray()
        offset = 0
        deadline = time.monotonic() + self.timeout
        try:
            for pipe in (process.stdin, process.stdout, process.stderr):
                os.set_blocking(pipe.fileno(), False)
            streams.register(process.stdin, selectors.EVENT_WRITE, "stdin")
            streams.register(process.stdout, selectors.EVENT_READ, "stdout")
            streams.register(process.stderr, selectors.EVENT_READ, "stderr")
            while streams.get_map():
                remaining = deadline - time.monotonic()
                if remaining <= 0.0:
                    self._kill_group_and_wait(process)
                    raise DomainFailure("CALCULATOR_FAILURE", "calculator callback timed out",
                                        component=component, requirement="CALC-006", object_id=object_id)
                ready = streams.select(remaining)
                if not ready:
                    continue
                for key, _mask in ready:
                    pipe = key.fileobj
                    if key.data == "stdin":
                        try:
                            offset += os.write(pipe.fileno(), request[offset:])
                        except BrokenPipeError:
                            offset = len(request)
                        if offset >= len(request):
                            streams.unregister(pipe)
                            pipe.close()
                    else:
                        target = stdout if key.data == "stdout" else stderr
                        limit = self.stdout_limit if key.data == "stdout" else self.stderr_limit
                        chunk = os.read(pipe.fileno(), min(65536, limit - len(target) + 1))
                        if not chunk:
                            streams.unregister(pipe)
                            pipe.close()
                            continue
                        target.extend(chunk)
                        if len(target) > limit:
                            self._kill_group_and_wait(process)
                            raise DomainFailure("CALCULATOR_FAILURE", f"calculator {key.data} byte limit exceeded",
                                                component=component, requirement="CALC-006",
                                                object_id=object_id)
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                self._kill_group_and_wait(process)
                raise DomainFailure("CALCULATOR_FAILURE", "calculator callback timed out",
                                    component=component, requirement="CALC-006", object_id=object_id)
            try:
                returncode = process.wait(timeout=remaining)
            except subprocess.TimeoutExpired as exc:
                self._kill_group_and_wait(process)
                raise DomainFailure("CALCULATOR_FAILURE", "calculator callback timed out",
                                    component=component, requirement="CALC-006", object_id=object_id) from exc
            return returncode, bytes(stdout), bytes(stderr)
        finally:
            streams.close()
            if process.poll() is None:
                self._kill_group_and_wait(process)

    def _validate_response(self, value: Any, state: Mapping[str, Any], request_digest: str,
                           component: str, object_id: str | None) -> Evaluation:
        required = {"status", "energy", "forces", "units", "model_name", "model_version",
                    "model_digest", "evaluation_id", "deterministic", "request_digest"}
        if not isinstance(value, dict) or set(value) != required or value.get("status") != "OK":
            raise DomainFailure("CALCULATOR_FAILURE", "malformed calculator response record",
                                component=component, requirement="CALC-002", object_id=object_id)
        units = value["units"]
        if units != {"energy": "eV", "forces": "eV/angstrom"}:
            raise DomainFailure("CALCULATOR_FAILURE", "calculator units mismatch",
                                component=component, requirement="CALC-006", object_id=object_id)
        if (value["model_name"] != self.model_name or value["model_version"] != self.model_version
                or value["model_digest"] != self.model_digest):
            raise DomainFailure("CALCULATOR_FAILURE", "calculator model identity mismatch",
                                component=component, requirement="CALC-006", object_id=object_id)
        if value["request_digest"] != request_digest:
            raise DomainFailure("CALCULATOR_FAILURE", "calculator request digest mismatch",
                                component=component, requirement="CALC-006", object_id=object_id)
        if not isinstance(value["evaluation_id"], str) or not value["evaluation_id"]:
            raise DomainFailure("CALCULATOR_FAILURE", "missing calculator evaluation ID",
                                component=component, requirement="CALC-002", object_id=object_id)
        if type(value["deterministic"]) is not bool:
            raise DomainFailure("CALCULATOR_FAILURE", "invalid deterministic flag",
                                component=component, requirement="CALC-002", object_id=object_id)
        energy = value["energy"]
        if isinstance(energy, bool) or not isinstance(energy, (int, float)):
            raise DomainFailure("CALCULATOR_FAILURE", "calculator energy type is malformed",
                                component=component, requirement="CALC-002", object_id=object_id)
        if not math.isfinite(float(energy)):
            raise DomainFailure("NONFINITE_RESULT", "calculator energy is nonfinite",
                                component=component, requirement="STATE-008", object_id=object_id)
        forces = value["forces"]
        count = len(state["atom_ids"])
        if not isinstance(forces, list) or len(forces) != count:
            raise DomainFailure("CALCULATOR_FAILURE", "calculator force shape mismatch",
                                component=component, requirement="CALC-006", object_id=object_id)
        parsed: list[Vector] = []
        for row in forces:
            if (not isinstance(row, list) or len(row) != 3
                    or any(isinstance(item, bool) or not isinstance(item, (int, float)) for item in row)):
                raise DomainFailure("CALCULATOR_FAILURE", "calculator force shape/type mismatch",
                                    component=component, requirement="CALC-006", object_id=object_id)
            vector = tuple(float(item) for item in row)
            if any(not math.isfinite(item) for item in vector):
                raise DomainFailure("NONFINITE_RESULT", "calculator force is nonfinite",
                                    component=component, requirement="STATE-008", object_id=object_id)
            parsed.append(vector)  # type: ignore[arg-type]
        deterministic = bool(value["deterministic"] and self.declared_deterministic)
        return Evaluation(float(energy), tuple(parsed), value["evaluation_id"], deterministic, request_digest)
