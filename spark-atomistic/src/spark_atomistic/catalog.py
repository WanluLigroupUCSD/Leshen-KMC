# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Errata-2 directed event transactions and discovery statistics."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping

from .canonical import canonical_bytes, deep_freeze, digest
from .errors import DomainFailure, Outcome
from .geometry import MatchReport, Vector, closest_periodic_displacement, match_positions, norm
from .model import (AtomicState, SCHEMA_SHA256, geometry_certificate,
                    identity_digest as make_identity_digest)
from .solvers import SaddleCandidate


def dot3(left: Vector, right: Vector) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def same_fixed_contract(left: AtomicState, right: AtomicState) -> bool:
    return left.fixed_contract_digest == right.fixed_contract_digest


def match_states(left: AtomicState, right: AtomicState, kinetics: Mapping[str, Any]) -> MatchReport:
    if not same_fixed_contract(left, right):
        return MatchReport(False, (), math.inf, math.inf, math.inf, (0.0, 0.0, 0.0))
    labels_left = tuple(f"{species}\0{int(movable)}"
                        for species, movable in zip(left.species, left.movable))
    labels_right = tuple(f"{species}\0{int(movable)}"
                         for species, movable in zip(right.species, right.movable))
    return match_positions(
        labels_left, left.positions, left.energy, labels_right, right.positions, right.energy,
        left.cell, left.pbc, rms_tolerance=kinetics["state_rms_tolerance"],
        max_tolerance=kinetics["state_max_tolerance"],
        energy_tolerance_per_atom=kinetics["state_energy_tolerance_per_atom"])


def _state_request(state: AtomicState) -> dict[str, Any]:
    return {"atom_ids": list(state.atom_ids), "species": list(state.species),
            "positions": [list(item) for item in state.positions],
            "cell": [list(item) for item in state.cell], "pbc": list(state.pbc),
            "movable": list(state.movable), "constraints": {"kind": "fixed-mask"},
            "charge": state.charge, "spin": state.spin,
            "calculator_model_digest": state.calculator_model_digest}


def _match_record(report: MatchReport) -> dict[str, Any]:
    return {"atom_mapping": list(report.atom_mapping),
            "energy_difference_ev": report.energy_difference,
            "max_displacement_angstrom": report.maximum_displacement,
            "rms_displacement_angstrom": report.rms_displacement}


def _inverse_match_record(report: MatchReport) -> dict[str, Any]:
    inverse = [0] * len(report.atom_mapping)
    for source, destination in enumerate(report.atom_mapping):
        inverse[destination] = source
    return {"atom_mapping": inverse,
            "energy_difference_ev": report.energy_difference,
            "max_displacement_angstrom": report.maximum_displacement,
            "rms_displacement_angstrom": report.rms_displacement}


@dataclass(frozen=True, slots=True)
class Event:
    event_id: str
    reverse_event_id: str
    pair_id: str
    origin_state_id: str
    destination_state_id: str
    saddle_positions: tuple[Vector, ...]
    saddle_energy: float
    saddle_forces: tuple[Vector, ...]
    unstable_direction: tuple[Vector, ...]
    curvature: float
    orthogonal_curvatures: tuple[float, ...]
    evidence_level: str
    evaluation_count: int
    search_id: str
    termination_reason: str
    barrier: float
    reverse_barrier: float
    log_rate: float
    reverse_log_rate: float
    detailed_balance_residual: float
    prefactor: float
    temperature: float
    active_atom_mapping: tuple[tuple[int, int], ...]
    discovery_provenance: Mapping[str, Any]
    validation: Mapping[str, Any]
    calculator_digest: str
    identity_digest: str
    tolerance_digest: str
    selectable: bool

    @property
    def enabled(self) -> bool:
        return self.selectable

    def record(self) -> dict[str, Any]:
        return {
            "active_atom_mapping": [list(item) for item in self.active_atom_mapping],
            "barrier_ev": self.barrier,
            "calculator_digest": self.calculator_digest,
            "destination_state_id": self.destination_state_id,
            "discovery_provenance": dict(self.discovery_provenance),
            "environment_key": "disabled", "environment_version": "none/1",
            "event_id": self.event_id, "identity_digest": self.identity_digest,
            "origin_state_id": self.origin_state_id, "pair_id": self.pair_id,
            "rate_model": {"common_prefactor_per_s": self.prefactor,
                           "detailed_balance_residual": self.detailed_balance_residual,
                           "log_forward_rate_per_s": self.log_rate,
                           "log_reverse_rate_per_s": self.reverse_log_rate,
                           "model": "COMMON_PREFACTOR", "temperature_k": self.temperature},
            "reverse_barrier_ev": self.reverse_barrier,
            "reverse_event_id": self.reverse_event_id,
            "saddle": {"curvature_ev_per_angstrom2": self.curvature,
                       "energy_ev": self.saddle_energy,
                       "evaluation_count": self.evaluation_count,
                       "evidence_level": self.evidence_level,
                       "forces_ev_per_angstrom": [list(item) for item in self.saddle_forces],
                       "orthogonal_curvatures_ev_per_angstrom2": list(self.orthogonal_curvatures),
                       "positions": [list(item) for item in self.saddle_positions],
                       "search_id": self.search_id,
                       "termination_reason": self.termination_reason,
                       "unstable_direction": [list(item) for item in self.unstable_direction]},
            "schema": "spark-atomistic-directed-event/2",
            "schema_digest": SCHEMA_SHA256, "selectable": self.selectable,
            "tolerance_digest": self.tolerance_digest,
            "validation": dict(self.validation),
        }


@dataclass(slots=True)
class DiscoveryRecord:
    state_id: str
    config_digest: str
    relevance_rate_min: float
    alpha: float | None
    alpha_calibration: Mapping[str, Any] | None
    attempts: int = 0
    successes: int = 0
    failures_by_status: dict[str, int] = field(default_factory=dict)
    duplicates: int = 0
    consecutive_redundant_successes: int = 0
    event_log_rates: dict[str, float] = field(default_factory=dict)
    heuristic_confidence: float | str = "UNAVAILABLE"
    evaluations: int = 0
    stopping_state: str = "RUNNING"
    permanently_incomplete_catalog: bool = False

    def record(self) -> dict[str, Any]:
        return {"alpha": self.alpha,
                "alpha_calibration": (None if self.alpha_calibration is None
                                      else dict(self.alpha_calibration)),
                "attempts": self.attempts, "config_digest": self.config_digest,
                "consecutive_redundant_successes": self.consecutive_redundant_successes,
                "duplicates": self.duplicates, "evaluations": self.evaluations,
                "event_log_rates": dict(sorted(self.event_log_rates.items())),
                "failures_by_status": dict(sorted(self.failures_by_status.items())),
                "heuristic_confidence": self.heuristic_confidence,
                "permanently_incomplete_catalog": self.permanently_incomplete_catalog,
                "relevance_rate_min": self.relevance_rate_min, "state_id": self.state_id,
                "stopping_state": self.stopping_state, "successes": self.successes}


class Catalog:
    def __init__(self, model_digest: str, config_digest: str, tolerance_digest: str,
                 maximum_events: int, identity_digest: str | None = None) -> None:
        self.model_digest = model_digest
        self.config_digest = config_digest
        self.tolerance_digest = tolerance_digest
        self.identity_digest = identity_digest or ""
        self.maximum_events = maximum_events
        self.states: dict[str, AtomicState] = {}
        self.events: dict[str, Event] = {}
        self.discovery: dict[str, DiscoveryRecord] = {}
        self.multiplicity: dict[str, int] = {}

    def add_state(self, candidate: AtomicState, kinetics: Mapping[str, Any]) -> AtomicState:
        resolved = self.resolve_state(candidate, kinetics)
        if resolved is candidate:
            self.states[candidate.state_id] = candidate
        return resolved

    def resolve_state(self, candidate: AtomicState, kinetics: Mapping[str, Any]) -> AtomicState:
        for state_id in sorted(self.states):
            if match_states(self.states[state_id], candidate, kinetics).equal:
                return self.states[state_id]
        return candidate

    def discovery_record(self, state_id: str, _mode: str,
                         discovery: Mapping[str, Any]) -> DiscoveryRecord:
        if state_id not in self.discovery:
            self.discovery[state_id] = DiscoveryRecord(
                state_id, self.config_digest, discovery["relevance_rate_min"],
                discovery["alpha"], discovery["alpha_calibration"])
        return self.discovery[state_id]

    def validate_candidate(self, origin: AtomicState, candidate: SaddleCandidate,
                           kinetics: Mapping[str, Any], provenance: Mapping[str, Any]) -> tuple[Outcome, Event | None, AtomicState | None]:
        if (not same_fixed_contract(origin, candidate.plus_endpoint)
                or not same_fixed_contract(origin, candidate.minus_endpoint)):
            raise DomainFailure("ATOM_COUNT_CHANGE_UNSUPPORTED", "atom count change unsupported",
                                component="catalog", requirement="E2-ID-003",
                                state_id=origin.state_id)
        plus_match = match_states(origin, candidate.plus_endpoint, kinetics)
        minus_match = match_states(origin, candidate.minus_endpoint, kinetics)
        if plus_match.equal and minus_match.equal:
            return Outcome("ENDPOINT_COLLAPSED", "", {"component": "catalog", "details": {},
                           "requirement_id": "E2-EVENT-005", "retryable": False,
                           "search_or_event_id": provenance.get("search_id"),
                           "state_id": origin.state_id}), None, None
        if not plus_match.equal and not minus_match.equal:
            return Outcome("SADDLE_WRONG_BASIN", "", {"component": "catalog", "details": {},
                           "requirement_id": "E2-EVENT-005", "retryable": False,
                           "search_or_event_id": provenance.get("search_id"),
                           "state_id": origin.state_id}), None, None
        endpoint = candidate.minus_endpoint if plus_match.equal else candidate.plus_endpoint
        origin_report = plus_match if plus_match.equal else minus_match
        destination = self.resolve_state(endpoint, kinetics)
        if destination.state_id == origin.state_id:
            return Outcome("ENDPOINT_COLLAPSED", "", {"component": "catalog", "details": {},
                           "requirement_id": "E2-EVENT-005", "retryable": False,
                           "search_or_event_id": provenance.get("search_id"),
                           "state_id": origin.state_id}), None, None
        destination_report = match_states(endpoint, destination, kinetics)
        if not destination_report.equal:
            raise DomainFailure("INVALID_STATE", "state invalid", component="catalog",
                                requirement="E2-ID-005", state_id=origin.state_id)
        raw_forward = candidate.energy - origin.energy
        raw_reverse = candidate.energy - destination.energy
        tolerance = kinetics["barrier_tolerance"]
        if raw_forward < -tolerance or raw_reverse < -tolerance:
            raise DomainFailure("RATE_INVALID", "rate invalid", component="rates",
                                requirement="E2-RATE-001", state_id=origin.state_id)
        beta = 1.0 / (8.617333262145e-5 * kinetics["temperature"])
        log_prefactor = math.log(kinetics["prefactor"])
        log_forward = log_prefactor - beta * raw_forward
        log_reverse = log_prefactor - beta * raw_reverse
        # `detailed_balance_residual` is mandatory (`E2-EVENT-003` fixes the `rate_model`
        # key set exactly), but on this path it carries NO information and NO gate is
        # applied to it. `E2-RATE-001` fixes both legs to one saddle energy ("Raw
        # barriers are exactly b_f=E_s-E_i and b_r=E_s-E_j"; `E2-EVENT-001` calls them
        # "Raw same-saddle differences"), and `E2-RATE-002` fixes the rates to those
        # barriers under the COMMON_PREFACTOR equations. Hence
        # log_forward - log_reverse = beta*(E_i - E_j) identically and the residual below
        # is algebraically zero for every input; measured worst |residual| over 20000
        # randomised physical inputs (T in {100,300,700,1200} K, nu in {1e12,1e13,6.2e12},
        # |dE| <= 3 eV, |E| <= 500 eV) is 1.1368683772161603e-13, pure binary64 rounding,
        # against a configured `detailed_balance_tolerance` of 1e-8.
        # No independent reverse barrier exists to compare against: the reverse endpoint
        # relaxation energy is not admissible here because `E2-EVENT-001` requires
        # `destination_state_id` to be a COMMITTED state and `E2-RATE-001` pins b_r to
        # that committed state's energy, which is exactly the quantity already used.
        # Manufacturing a second, independently obtained reverse barrier would require a
        # separate reverse saddle search or a Hessian-derived asymmetric prefactor, i.e.
        # new physics that no requirement asks for. A tolerance gate here therefore has
        # zero power to reject a detailed-balance violation and can only reject valid
        # physics on rounding noise: with `detailed_balance_tolerance` set below the
        # rounding floor the deleted gate rejected 11829 of the same 20000 valid inputs.
        # Detailed balance IS verified where the two sides are independent inputs: at
        # checkpoint restore (`E2-CKPT-007`(5)), `checkpoint._restore_event` recomputes
        # the residual from the stored barriers and energies and refuses any stored value
        # that disagrees.
        residual = (log_forward - log_reverse) + beta * (destination.energy - origin.energy)
        if not all(math.isfinite(item) for item in (log_forward, log_reverse, residual)):
            raise DomainFailure("RATE_INVALID", "rate invalid", component="rates",
                                requirement="E2-RATE-001", state_id=origin.state_id)
        mapping = tuple(destination_report.atom_mapping)
        active_pairs = tuple(sorted(
            (source, mapping[source]) for source, (left, right) in enumerate(
                zip(origin.positions, endpoint.positions))
            if norm(closest_periodic_displacement(
                (left[0] - right[0], left[1] - right[1], left[2] - right[2]),
                origin.cell, origin.pbc)[0]) > kinetics["state_rms_tolerance"]))
        reverse_active = tuple(sorted((destination_index, origin_index)
                                      for origin_index, destination_index in active_pairs))
        saddle_geometry_digest = digest(geometry_certificate(
            _state_request(origin), candidate.positions))
        flat_mode = [component for vector in candidate.mode for component in vector]
        negated_mode = [-component for component in flat_mode]
        canonical_mode = min((flat_mode, negated_mode), key=canonical_bytes)
        endpoint_ids = sorted((origin.state_id, destination.state_id))
        oriented_active = active_pairs if origin.state_id == endpoint_ids[0] else reverse_active
        pair_id = "pair:" + digest({"active_atom_mapping": [list(item) for item in oriented_active],
                                     "endpoint_state_ids": endpoint_ids,
                                     "saddle_energy_ev": candidate.energy,
                                     "saddle_geometry_digest": saddle_geometry_digest,
                                     "schema": "spark-atomistic-event-pair/2",
                                     "unstable_direction": canonical_mode})
        forward_id = "event:" + digest({"destination_state_id": destination.state_id,
                                        "origin_state_id": origin.state_id,
                                        "pair_id": pair_id,
                                        "schema": "spark-atomistic-directed-event/2"})
        reverse_id = "event:" + digest({"destination_state_id": origin.state_id,
                                        "origin_state_id": destination.state_id,
                                        "pair_id": pair_id,
                                        "schema": "spark-atomistic-directed-event/2"})
        duplicate = self._find_duplicate(origin.state_id, destination.state_id,
                                         candidate.positions, candidate.energy,
                                         candidate.mode, active_pairs, kinetics)
        if duplicate is not None:
            reverse_duplicate = self.events[duplicate].reverse_event_id
            self.multiplicity[duplicate] += 1
            self.multiplicity[reverse_duplicate] += 1
            return Outcome("DUPLICATE_EVENT", "", {"component": "catalog", "details": {},
                           "requirement_id": "E2-DISC-005", "retryable": False,
                           "search_or_event_id": duplicate, "state_id": origin.state_id}), None, destination
        if len(self.events) + 2 > self.maximum_events:
            raise DomainFailure("RESOURCE_LIMIT", "resource limit reached", component="catalog",
                                requirement="E2-SCHEMA-009", state_id=origin.state_id)
        search_id = str(provenance["search_id"])
        discovery_provenance = deep_freeze({
            "rng_substream_digest": provenance["rng_substream_digest"],
            "search_class": provenance["search_class"], "search_id": search_id,
            "search_index": provenance["search_index"]})
        validation = deep_freeze({
            "calculator_model_digest": self.model_digest,
            "constraint_digest": origin.constraint_digest,
            "destination_match": _match_record(destination_report),
            "full_endpoint_relaxations": True, "method": "full-endpoint-relaxation/1",
            "origin_match": _match_record(origin_report), "unstable_mode_count": 1})
        reverse_validation = deep_freeze({
            "calculator_model_digest": self.model_digest,
            "constraint_digest": destination.constraint_digest,
            "destination_match": _inverse_match_record(origin_report),
            "full_endpoint_relaxations": True, "method": "full-endpoint-relaxation/1",
            "origin_match": _inverse_match_record(destination_report),
            "unstable_mode_count": 1})
        reverse_positions: list[Vector | None] = [None] * len(candidate.positions)
        reverse_forces: list[Vector | None] = [None] * len(candidate.forces)
        reverse_mode: list[Vector | None] = [None] * len(candidate.mode)
        for source, target in enumerate(mapping):
            reverse_positions[target] = candidate.positions[source]
            reverse_forces[target] = candidate.forces[source]
            reverse_mode[target] = tuple(-component for component in candidate.mode[source])
        if any(item is None for item in reverse_positions + reverse_forces + reverse_mode):
            raise DomainFailure("INVALID_STATE", "state invalid", component="catalog",
                                requirement="E2-EVENT-001", state_id=origin.state_id)
        identity_value = self.identity_digest or make_identity_digest(kinetics)
        forward = Event(
            forward_id, reverse_id, pair_id, origin.state_id, destination.state_id,
            candidate.positions, candidate.energy, candidate.forces, candidate.mode,
            candidate.curvature, candidate.orthogonal_curvatures, candidate.evidence_level,
            candidate.evaluations, search_id,
            str(candidate.provenance.get("termination_reason", "full-force-curvature-orthogonal-gates")),
            raw_forward, raw_reverse, log_forward, log_reverse, residual,
            kinetics["prefactor"], kinetics["temperature"], active_pairs,
            discovery_provenance, validation, self.model_digest, identity_value,
            self.tolerance_digest, log_forward >= kinetics["log_rate_cutoff"])
        reverse = Event(
            reverse_id, forward_id, pair_id, destination.state_id, origin.state_id,
            tuple(reverse_positions), candidate.energy, tuple(reverse_forces),
            tuple(reverse_mode),
            candidate.curvature, candidate.orthogonal_curvatures, candidate.evidence_level,
            candidate.evaluations, search_id,
            str(candidate.provenance.get("termination_reason", "full-force-curvature-orthogonal-gates")),
            raw_reverse, raw_forward, log_reverse, log_forward, -residual,
            kinetics["prefactor"], kinetics["temperature"], reverse_active,
            discovery_provenance, reverse_validation, self.model_digest, identity_value,
            self.tolerance_digest, log_reverse >= kinetics["log_rate_cutoff"])
        if destination.state_id not in self.states:
            self.states[destination.state_id] = destination
        self.events.update({forward_id: forward, reverse_id: reverse})
        self.multiplicity.update({forward_id: 1, reverse_id: 1})
        return Outcome("OK", "", {"component": "catalog", "details": {},
                       "requirement_id": "E2-EVENT-001", "retryable": False,
                       "search_or_event_id": forward_id, "state_id": origin.state_id}), forward, destination

    def _find_duplicate(self, origin_state_id: str, destination_state_id: str,
                        saddle_positions: tuple[Vector, ...], saddle_energy: float,
                        unstable_direction: tuple[Vector, ...],
                        active_atom_mapping: tuple[tuple[int, int], ...],
                        kinetics: Mapping[str, Any]) -> str | None:
        origin = self.states[origin_state_id]
        labels = tuple(f"{species}\0{int(movable)}"
                       for species, movable in zip(origin.species, origin.movable))
        for event_id in sorted(self.events):
            event = self.events[event_id]
            if (event.origin_state_id != origin_state_id
                    or event.destination_state_id != destination_state_id):
                continue
            report = match_positions(labels, event.saddle_positions, 0.0,
                                     labels, saddle_positions, 0.0, origin.cell, origin.pbc,
                                     rms_tolerance=kinetics["saddle_rms_tolerance"],
                                     max_tolerance=kinetics["saddle_max_tolerance"],
                                     energy_tolerance_per_atom=0.0)
            if not report.equal:
                continue
            if abs(event.saddle_energy - saddle_energy) > kinetics["saddle_energy_tolerance"]:
                raise DomainFailure("CATALOG_CONFLICT", "catalog conflict", component="catalog",
                                    requirement="E2-EVENT-006", object_id=event_id)
            mapped_dot = sum(dot3(event.unstable_direction[source], unstable_direction[target])
                             for source, target in enumerate(report.atom_mapping))
            if abs(mapped_dot) >= 1.0 - 1e-6 and event.active_atom_mapping == active_atom_mapping:
                return event_id
        return None

    def checkpoint(self) -> dict[str, Any]:
        payload = {"events": {key: self.events[key].record() for key in sorted(self.events)},
                   "multiplicity": dict(sorted(self.multiplicity.items())),
                   "schema": "spark-atomistic-catalog/2",
                   "states": {key: self.states[key].record() for key in sorted(self.states)}}
        return {"digest": digest(payload), **payload}
