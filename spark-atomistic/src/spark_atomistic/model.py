# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Errata-2 exact model schema and canonical state identities."""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

from .canonical import JSON_INTEGER_LIMIT, canonical_bytes, canonical_text, deep_freeze, digest
from .errors import DomainFailure
from .geometry import Cell, Vector, closest_periodic_displacement, max_movable_force


IR = "spark-atomistic-model/1"
SCHEMA_REVISION = 2
IDENTITY_VERSION = "spark-state-identity/2"
STATE_SCHEMA = "spark-atomistic-state/2"
SCHEMA_DESCRIPTOR = {
    "base_spec_sha256": "8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84",
    "errata_1_sha256": "52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40",
    "ir": IR,
    "revision": 2,
}
SCHEMA_SHA256 = digest(SCHEMA_DESCRIPTOR)
EXPECTED_SCHEMA_SHA256 = "sha256:583d580d54e3847ef92f1b1456dda006161689c0bac27fd7ea896a093f48c02c"
# Errata 2 froze the wire contract but deliberately did NOT enter SCHEMA_DESCRIPTOR: the
# descriptor is hashed into SCHEMA_SHA256, which E2-PAR-003 requires to stay byte-identical
# across backends. The provenance SHA is therefore carried as a separate machine-readable
# constant, mirroring `ERRATA_2_SHA256` in the Rust backend.
BASE_SPEC_SHA256 = SCHEMA_DESCRIPTOR["base_spec_sha256"]
ERRATA_1_SHA256 = SCHEMA_DESCRIPTOR["errata_1_sha256"]
ERRATA_2_SHA256 = "eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995"


def _object(value: Any, path: str, required: set[str], optional: set[str] = frozenset()) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema",
                            requirement="E2-SCOPE-003", details={"path": path})
    missing = required - set(value)
    unknown = set(value) - required - optional
    if missing or unknown:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema",
                            requirement="E2-SCOPE-003",
                            details={"path": path, "missing": sorted(missing),
                                     "unknown": sorted(unknown)})
    return value


def _string(value: Any, path: str) -> str:
    if not isinstance(value, str) or not value:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema",
                            requirement="E2-SCHEMA-001", details={"path": path})
    return value


def _bool(value: Any, path: str) -> bool:
    if type(value) is not bool:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema",
                            requirement="E2-SCHEMA-001", details={"path": path})
    return value


def _integer(value: Any, path: str, minimum: int = 0) -> int:
    if type(value) is not int or not minimum <= value <= JSON_INTEGER_LIMIT:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema",
                            requirement="E2-JSON-002", details={"path": path})
    return value


def _number(value: Any, path: str, *, positive: bool = False,
            nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema",
                            requirement="E2-SCHEMA-001", details={"path": path})
    result = float(value)
    if not math.isfinite(result):
        raise DomainFailure("NONFINITE_RESULT", "nonfinite value rejected", component="schema",
                            requirement="E2-JSON-003", details={"path": path})
    if positive and result <= 0.0 or nonnegative and result < 0.0:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema",
                            requirement="E2-SCHEMA-001", details={"path": path})
    return result


def _vector(value: Any, path: str) -> Vector:
    if not isinstance(value, list) or len(value) != 3:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema",
                            requirement="E2-SCHEMA-003", details={"path": path})
    return tuple(_number(item, f"{path}[{index}]") for index, item in enumerate(value))  # type: ignore[return-value]


def _determinant(cell: Sequence[Sequence[float]]) -> float:
    return (cell[0][0] * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1])
            - cell[0][1] * (cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0])
            + cell[0][2] * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0]))


def geometry_certificate(request: Mapping[str, Any],
                         positions: Sequence[Sequence[float]] | None = None) -> dict[str, Any]:
    coordinates = request["positions"] if positions is None else positions
    certificates: list[list[list[Any]]] = []
    for anchor_position in coordinates:
        rows: list[list[Any]] = []
        for species, movable, position in zip(request["species"], request["movable"], coordinates):
            displacement, _shift = closest_periodic_displacement(
                (position[0] - anchor_position[0], position[1] - anchor_position[1],
                 position[2] - anchor_position[2]), request["cell"], request["pbc"])
            rows.append([species, movable, displacement[0], displacement[1], displacement[2]])
        rows.sort(key=lambda row: tuple(canonical_text(item) for item in row))
        certificates.append(rows)
    rows = min(certificates, key=canonical_bytes)
    return {"cell": [list(item) for item in request["cell"]], "pbc": list(request["pbc"]),
            "rows": rows, "version": "anchor-minimum-closest-image/1"}


def fixed_contract_digest(request: Mapping[str, Any]) -> str:
    atom_contracts = sorted(([species, movable]
                             for species, movable in zip(request["species"], request["movable"])),
                            key=canonical_bytes)
    return digest({"atom_contracts": atom_contracts,
                   "calculator_model_digest": request["calculator_model_digest"],
                   "cell": [list(item) for item in request["cell"]],
                   "charge": request["charge"], "constraints": {"kind": "fixed-mask"},
                   "pbc": list(request["pbc"]), "spin": request["spin"]})


def identity_digest(kinetics: Mapping[str, Any]) -> str:
    return digest({"reflection_invariant": False, "rotation_invariant": False,
                   "state_energy_tolerance_per_atom": kinetics["state_energy_tolerance_per_atom"],
                   "state_max_tolerance": kinetics["state_max_tolerance"],
                   "state_rms_tolerance": kinetics["state_rms_tolerance"],
                   "version": IDENTITY_VERSION})


@dataclass(frozen=True, slots=True)
class AtomicState:
    state_id: str
    candidate_identity: str
    candidate_identity_version: str
    fixed_contract_digest: str
    atom_ids: tuple[str, ...]
    species: tuple[str, ...]
    positions: tuple[Vector, ...]
    cell: Cell
    pbc: tuple[bool, bool, bool]
    movable: tuple[bool, ...]
    constraints: Mapping[str, Any]
    charge: float
    spin: float
    calculator_model_digest: str
    energy: float
    forces: tuple[Vector, ...]
    maximum_movable_force: float
    force_tolerance: float
    constraint_digest: str
    relaxation_provenance: Mapping[str, Any]

    def record(self) -> dict[str, Any]:
        return {
            "atom_ids": list(self.atom_ids),
            "calculator_model_digest": self.calculator_model_digest,
            "candidate_identity": self.candidate_identity,
            "cell": [list(item) for item in self.cell], "charge": self.charge,
            "constraint_digest": self.constraint_digest,
            "constraints": {"kind": "fixed-mask"}, "energy_ev": self.energy,
            "fixed_contract_digest": self.fixed_contract_digest,
            "force_tolerance_ev_per_angstrom": self.force_tolerance,
            "forces_ev_per_angstrom": [list(item) for item in self.forces],
            "identity_version": self.candidate_identity_version,
            "max_movable_force_ev_per_angstrom": self.maximum_movable_force,
            "movable": list(self.movable), "pbc": list(self.pbc),
            "positions": [list(item) for item in self.positions],
            "relaxation_provenance": dict(self.relaxation_provenance),
            "schema": STATE_SCHEMA, "species": list(self.species), "spin": self.spin,
            "state_id": self.state_id,
        }


def state_request_from_system(system: Mapping[str, Any]) -> dict[str, Any]:
    return {"schema": IR, "atom_ids": list(system["atom_ids"]),
            "species": list(system["species"]), "positions": [list(item) for item in system["positions"]],
            "cell": [list(item) for item in system["cell"]], "pbc": list(system["pbc"]),
            "movable": list(system["movable"]), "constraints": {"kind": "fixed-mask"},
            "charge": system["charge"], "spin": system["spin"],
            "calculator_model_digest": system["calculator_model_digest"]}


def state_from_relaxation(request: Mapping[str, Any], energy: float, forces: Sequence[Vector],
                          tolerance: float, provenance: Mapping[str, Any]) -> AtomicState:
    fixed_digest = fixed_contract_digest(request)
    certificate = geometry_certificate(request)
    candidate = "candidate:" + digest({"fixed_contract_digest": fixed_digest,
                                        "geometry": certificate, "version": IDENTITY_VERSION})
    state_id = "state:" + digest({"candidate_identity": candidate, "energy_ev": float(energy),
                                  "fixed_contract_digest": fixed_digest,
                                  "version": IDENTITY_VERSION})
    constraint = digest({"constraints": {"kind": "fixed-mask"},
                         "movable": list(request["movable"])})
    exact_provenance = {
        "calculator_evaluations": int(provenance.get("calculator_evaluations", 0)),
        "calculator_identity": str(provenance.get("calculator_identity",
                                                   request["calculator_model_digest"])),
        "minimizer_identity": str(provenance.get("minimizer_identity",
                                                  provenance.get("method", "reference-minimizer/1"))),
        "steps": int(provenance.get("steps", 0)),
        "termination_reason": str(provenance.get("termination_reason", "force_tolerance")),
    }
    return AtomicState(
        state_id, candidate, IDENTITY_VERSION, fixed_digest, tuple(request["atom_ids"]),
        tuple(request["species"]), tuple(tuple(item) for item in request["positions"]),
        tuple(tuple(item) for item in request["cell"]), tuple(request["pbc"]),
        tuple(request["movable"]), deep_freeze({"kind": "fixed-mask"}),
        float(request["charge"]), float(request["spin"]), request["calculator_model_digest"],
        float(energy), tuple(tuple(item) for item in forces),
        max_movable_force(forces, request["movable"]), float(tolerance), constraint,
        deep_freeze(exact_provenance))


def validate_state_record(value: Any, *, checkpoint: bool = False) -> AtomicState:
    required = {"atom_ids", "calculator_model_digest", "candidate_identity", "cell", "charge",
                "constraint_digest", "constraints", "energy_ev", "fixed_contract_digest",
                "force_tolerance_ev_per_angstrom", "forces_ev_per_angstrom", "identity_version",
                "max_movable_force_ev_per_angstrom", "movable", "pbc", "positions",
                "relaxation_provenance", "schema", "species", "spin", "state_id"}
    try:
        obj = _object(value, "state", required)
        if obj["schema"] != STATE_SCHEMA or obj["identity_version"] != IDENTITY_VERSION:
            raise DomainFailure("INVALID_STATE", "state invalid", component="state",
                                requirement="E2-ID-006")
        count = len(obj["atom_ids"]) if isinstance(obj["atom_ids"], list) else -1
        if (count < 1 or len(set(obj["atom_ids"])) != count
                or any(not isinstance(item, str) or not item for item in obj["atom_ids"])
                or not isinstance(obj["species"], list) or len(obj["species"]) != count
                or any(not isinstance(item, str) or not item for item in obj["species"])):
            raise DomainFailure("INVALID_STATE", "state invalid", component="state",
                                requirement="E2-ID-006")
        positions = tuple(_vector(item, "state.positions") for item in obj["positions"])
        cell = tuple(_vector(item, "state.cell") for item in obj["cell"])
        if len(positions) != count or len(cell) != 3 or _determinant(cell) == 0.0:
            raise DomainFailure("INVALID_STATE", "state invalid", component="state",
                                requirement="E2-ID-006")
        if (not isinstance(obj["pbc"], list) or len(obj["pbc"]) != 3
                or any(type(item) is not bool for item in obj["pbc"])
                or not isinstance(obj["movable"], list) or len(obj["movable"]) != count
                or any(type(item) is not bool for item in obj["movable"])):
            raise DomainFailure("INVALID_STATE", "state invalid", component="state",
                                requirement="E2-ID-006")
        if obj["constraints"] != {"kind": "fixed-mask"}:
            raise DomainFailure("INVALID_STATE", "state invalid", component="state",
                                requirement="E2-ID-006")
        forces = tuple(_vector(item, "state.forces") for item in obj["forces_ev_per_angstrom"])
        if len(forces) != count:
            raise DomainFailure("INVALID_STATE", "state invalid", component="state",
                                requirement="E2-ID-006")
        request = {"atom_ids": obj["atom_ids"], "species": obj["species"],
                   "positions": positions, "cell": cell, "pbc": obj["pbc"],
                   "movable": obj["movable"], "constraints": obj["constraints"],
                   "charge": _number(obj["charge"], "state.charge"),
                   "spin": _number(obj["spin"], "state.spin"),
                   "calculator_model_digest": _string(obj["calculator_model_digest"], "state.calculator_model_digest")}
        energy = _number(obj["energy_ev"], "state.energy_ev")
        tolerance = _number(obj["force_tolerance_ev_per_angstrom"], "state.force_tolerance", positive=True)
        actual_maximum = max_movable_force(forces, obj["movable"])
        stored_maximum = _number(obj["max_movable_force_ev_per_angstrom"], "state.max_force", nonnegative=True)
        if stored_maximum != actual_maximum or stored_maximum > tolerance:
            raise DomainFailure("INVALID_STATE", "state invalid", component="state",
                                requirement="E2-CKPT-007")
        provenance = _object(obj["relaxation_provenance"], "state.relaxation_provenance",
                             {"calculator_evaluations", "calculator_identity", "minimizer_identity",
                              "steps", "termination_reason"})
        _integer(provenance["calculator_evaluations"], "state.provenance.evaluations")
        _integer(provenance["steps"], "state.provenance.steps")
        for key in ("calculator_identity", "minimizer_identity", "termination_reason"):
            _string(provenance[key], f"state.provenance.{key}")
        candidate = state_from_relaxation(request, energy, forces, tolerance, provenance)
        if (obj["fixed_contract_digest"] != candidate.fixed_contract_digest
                or obj["constraint_digest"] != candidate.constraint_digest
                or obj["candidate_identity"] != candidate.candidate_identity
                or obj["state_id"] != candidate.state_id):
            raise DomainFailure("INVALID_STATE", "state invalid", component="state",
                                requirement="E2-ID-004")
        return candidate
    except DomainFailure as exc:
        if checkpoint and exc.outcome.status not in {"CHECKPOINT_CORRUPT", "CHECKPOINT_INCOMPATIBLE"}:
            raise DomainFailure("CHECKPOINT_CORRUPT", "checkpoint corrupt", component="checkpoint",
                                requirement="E2-CKPT-007", causal_status=exc.outcome.status) from exc
        raise


def validate_model(value: Any, *, source_path: str | None = None) -> Mapping[str, Any]:
    root = _object(value, "$", {"schema", "system", "calculator", "relaxation", "saddle_search",
                                  "discovery", "kinetics", "resources", "output", "basin"},
                   {"metadata"})
    schema = _object(root["schema"], "$.schema", {"id"})
    if schema["id"] != IR:
        raise DomainFailure("SCHEMA_UNSUPPORTED", "schema unsupported", component="schema",
                            requirement="E2-SCHEMA-002")
    system = _object(root["system"], "$.system", {"atom_ids", "species", "positions", "cell", "pbc",
                                                                    "movable", "constraints", "charge", "spin",
                                                                    "calculator_model_digest"})
    count = len(system["atom_ids"]) if isinstance(system["atom_ids"], list) else -1
    if (count < 1 or len(set(system["atom_ids"])) != count
            or any(not isinstance(item, str) or not item for item in system["atom_ids"])):
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-003")
    if (not isinstance(system["species"], list) or len(system["species"]) != count
            or any(not isinstance(item, str) or not item for item in system["species"])
            or not isinstance(system["positions"], list) or len(system["positions"]) != count):
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-003")
    system["positions"] = [_vector(item, "$.system.positions") for item in system["positions"]]
    if not isinstance(system["cell"], list) or len(system["cell"]) != 3:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-003")
    system["cell"] = [_vector(item, "$.system.cell") for item in system["cell"]]
    if _determinant(system["cell"]) == 0.0:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-003")
    if (not isinstance(system["pbc"], list) or len(system["pbc"]) != 3
            or any(type(item) is not bool for item in system["pbc"])
            or not isinstance(system["movable"], list) or len(system["movable"]) != count
            or any(type(item) is not bool for item in system["movable"]) or not any(system["movable"])
            or system["constraints"] != {"kind": "fixed-mask"}):
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-003")
    system["charge"] = _number(system["charge"], "$.system.charge")
    system["spin"] = _number(system["spin"], "$.system.spin")
    _string(system["calculator_model_digest"], "$.system.calculator_model_digest")

    calculator = _object(root["calculator"], "$.calculator",
                         {"deterministic", "model_digest", "model_name", "model_version"})
    _bool(calculator["deterministic"], "$.calculator.deterministic")
    for key in ("model_digest", "model_name", "model_version"):
        _string(calculator[key], f"$.calculator.{key}")
    if calculator["model_digest"] != system["calculator_model_digest"]:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-004")

    relaxation = _object(root["relaxation"], "$.relaxation",
                         {"force_tolerance", "max_evaluations", "max_steps"})
    relaxation["force_tolerance"] = _number(relaxation["force_tolerance"], "$.relaxation.force_tolerance", positive=True)
    relaxation["max_evaluations"] = _integer(relaxation["max_evaluations"], "$.relaxation.max_evaluations", 1)
    relaxation["max_steps"] = _integer(relaxation["max_steps"], "$.relaxation.max_steps", 1)

    saddle = _object(root["saddle_search"], "$.saddle_search",
                     {"curvature_tolerance", "endpoint_displacement", "force_tolerance",
                      "max_iterations", "method", "orthogonal_directions"})
    if saddle["method"] != "directional-dimer":
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-006")
    for key in ("curvature_tolerance", "endpoint_displacement", "force_tolerance"):
        saddle[key] = _number(saddle[key], f"$.saddle_search.{key}", positive=True)
    saddle["max_iterations"] = _integer(saddle["max_iterations"], "$.saddle_search.max_iterations", 1)
    saddle["orthogonal_directions"] = _integer(saddle["orthogonal_directions"], "$.saddle_search.orthogonal_directions", 1)

    discovery = _object(root["discovery"], "$.discovery",
                        {"mode", "classes", "minimum_successful", "consecutive_redundant",
                         "maximum_attempts", "maximum_evaluations", "relevance_rate_min",
                         "alpha", "alpha_calibration"})
    if discovery["mode"] not in {"strict", "exploratory"} or not isinstance(discovery["classes"], list) or not discovery["classes"]:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-007")
    names: set[str] = set()
    probability_sum = 0.0
    non_targeted = False
    for index, item in enumerate(discovery["classes"]):
        entry = _object(item, f"$.discovery.classes[{index}]", {"kind", "name", "probability"})
        name = _string(entry["name"], f"$.discovery.classes[{index}].name")
        if name in names or entry["kind"] not in {"global", "local", "targeted"}:
            raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-007")
        names.add(name)
        entry["probability"] = _number(entry["probability"], f"$.discovery.classes[{index}].probability", positive=True)
        probability_sum += entry["probability"]
        non_targeted = non_targeted or entry["kind"] != "targeted"
    if abs(probability_sum - 1.0) > 1e-12 or not non_targeted:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-007")
    for key in ("minimum_successful", "consecutive_redundant", "maximum_attempts", "maximum_evaluations"):
        discovery[key] = _integer(discovery[key], f"$.discovery.{key}", 1)
    discovery["relevance_rate_min"] = _number(discovery["relevance_rate_min"], "$.discovery.relevance_rate_min", nonnegative=True)
    alpha = discovery["alpha"]
    calibration = discovery["alpha_calibration"]
    if alpha is None:
        if calibration is not None:
            raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-007")
    else:
        discovery["alpha"] = _number(alpha, "$.discovery.alpha", positive=True)
        if discovery["alpha"] > 1.0:
            raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-007")
        calibration = _object(calibration, "$.discovery.alpha_calibration", {"id", "source", "version"})
        for key in ("id", "source", "version"):
            _string(calibration[key], f"$.discovery.alpha_calibration.{key}")

    kinetics = _object(root["kinetics"], "$.kinetics",
                       {"temperature", "rate_model", "prefactor", "barrier_tolerance",
                        "detailed_balance_tolerance", "log_rate_cutoff", "absorbing_ok",
                        "maximum_steps", "run_seed", "state_rms_tolerance", "state_max_tolerance",
                        "state_energy_tolerance_per_atom", "saddle_rms_tolerance",
                        "saddle_max_tolerance", "saddle_energy_tolerance"})
    if kinetics["rate_model"] != "COMMON_PREFACTOR":
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-008")
    for key in ("temperature", "prefactor", "detailed_balance_tolerance", "state_rms_tolerance",
                "state_max_tolerance", "state_energy_tolerance_per_atom", "saddle_rms_tolerance",
                "saddle_max_tolerance", "saddle_energy_tolerance"):
        kinetics[key] = _number(kinetics[key], f"$.kinetics.{key}", positive=True)
    kinetics["barrier_tolerance"] = _number(kinetics["barrier_tolerance"], "$.kinetics.barrier_tolerance", nonnegative=True)
    kinetics["log_rate_cutoff"] = _number(kinetics["log_rate_cutoff"], "$.kinetics.log_rate_cutoff")
    kinetics["absorbing_ok"] = _bool(kinetics["absorbing_ok"], "$.kinetics.absorbing_ok")
    kinetics["maximum_steps"] = _integer(kinetics["maximum_steps"], "$.kinetics.maximum_steps", 1)
    kinetics["run_seed"] = _integer(kinetics["run_seed"], "$.kinetics.run_seed", 0)
    beta = 1.0 / (8.617333262145e-5 * kinetics["temperature"])
    log_prefactor = math.log(kinetics["prefactor"])
    if (not math.isfinite(beta) or not math.isfinite(log_prefactor)
            or not math.log(math.nextafter(0.0, 1.0)) <= kinetics["log_rate_cutoff"] <= math.log(sys.float_info.max)):
        raise DomainFailure("RATE_INVALID", "rate invalid", component="schema", requirement="E2-SCHEMA-008")

    resources = _object(root["resources"], "$.resources",
                        {"callback_timeout_s", "catalog_events", "evaluations_per_relaxation",
                         "evaluations_per_saddle_attempt", "output_bytes", "resident_memory_bytes",
                         "retry_backoff_s", "retry_count", "saddle_attempts_per_state",
                         "total_calculator_evaluations", "wall_time_s"})
    resources["callback_timeout_s"] = _number(resources["callback_timeout_s"], "$.resources.callback_timeout_s", positive=True)
    resources["wall_time_s"] = _number(resources["wall_time_s"], "$.resources.wall_time_s", positive=True)
    for key in ("catalog_events", "evaluations_per_relaxation", "evaluations_per_saddle_attempt",
                "output_bytes", "resident_memory_bytes", "saddle_attempts_per_state",
                "total_calculator_evaluations"):
        resources[key] = _integer(resources[key], f"$.resources.{key}", 1)
    if resources["retry_count"] != 0 or resources["retry_backoff_s"] != 0:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-009")
    resources["retry_count"] = _integer(resources["retry_count"], "$.resources.retry_count", 0)
    resources["retry_backoff_s"] = _number(resources["retry_backoff_s"], "$.resources.retry_backoff_s", nonnegative=True)

    output = _object(root["output"], "$.output",
                     {"checkpoint_every_steps", "checkpoint_path", "checkpoint_wall_time_s",
                      "overwrite", "resume", "summary_path", "trajectory_path"})
    output["checkpoint_every_steps"] = _integer(output["checkpoint_every_steps"], "$.output.checkpoint_every_steps", 1)
    output["checkpoint_wall_time_s"] = _number(output["checkpoint_wall_time_s"], "$.output.checkpoint_wall_time_s", positive=True)
    output["overwrite"] = _bool(output["overwrite"], "$.output.overwrite")
    output["resume"] = _bool(output["resume"], "$.output.resume")
    resolved: list[str] = []
    for key in ("checkpoint_path", "summary_path", "trajectory_path"):
        path = _string(output[key], f"$.output.{key}")
        if not os.path.isabs(path):
            if source_path is None:
                raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-010")
            path = str((Path(source_path).resolve().parent / path).resolve())
        output[key] = path
        resolved.append(path)
    if len(set(resolved)) != 3:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-010")

    basin = _object(root["basin"], "$.basin", {"enabled"})
    basin["enabled"] = _bool(basin["enabled"], "$.basin.enabled")
    if relaxation["max_evaluations"] > resources["evaluations_per_relaxation"]:
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-012")
    if (discovery["maximum_attempts"] > resources["saddle_attempts_per_state"]
            or discovery["maximum_evaluations"] > resources["total_calculator_evaluations"]):
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCHEMA-012")
    if "metadata" in root and not isinstance(root["metadata"], dict):
        raise DomainFailure("INVALID_INPUT", "input invalid", component="schema", requirement="E2-SCOPE-003")
    return deep_freeze(root)
