# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Reference steepest-descent minimizer and directional dimer search."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Any, Mapping, Sequence

from .calculator import Evaluation, ProcessCalculator
from .errors import DomainFailure
from .geometry import Vector, dot, max_movable_force, norm, vscale
from .model import AtomicState, state_from_relaxation
from .rng import PhiloxStream


def _copy_request(state: Mapping[str, Any], positions: Sequence[Vector]) -> dict[str, Any]:
    return {
        "schema": state["schema"], "atom_ids": list(state["atom_ids"]),
        "species": list(state["species"]), "positions": [list(item) for item in positions],
        "cell": [list(item) for item in state["cell"]], "pbc": list(state["pbc"]),
        "movable": list(state["movable"]), "constraints": dict(state["constraints"]),
        "charge": state["charge"], "spin": state["spin"],
        "calculator_model_digest": state["calculator_model_digest"],
    }


def _state_request(state: AtomicState) -> dict[str, Any]:
    record = state.record()
    return {key: record[key] for key in (
        "atom_ids", "species", "positions", "cell", "pbc", "movable", "constraints",
        "charge", "spin", "calculator_model_digest"
    )} | {"schema": "spark-atomistic-model/1"}


def _flatten(vectors: Sequence[Vector], movable: Sequence[bool]) -> list[float]:
    result: list[float] = []
    for vector, active in zip(vectors, movable):
        result.extend(vector if active else (0.0, 0.0, 0.0))
    return result


def _unflatten(values: Sequence[float]) -> tuple[Vector, ...]:
    return tuple((values[index], values[index + 1], values[index + 2])
                 for index in range(0, len(values), 3))


def _flat_dot(a: Sequence[float], b: Sequence[float]) -> float:
    return sum(left * right for left, right in zip(a, b))


def _flat_normalize(values: Sequence[float]) -> list[float]:
    length = math.sqrt(_flat_dot(values, values))
    if not math.isfinite(length) or length == 0.0:
        raise DomainFailure("INVALID_SADDLE", "dimer direction is zero/nonfinite",
                            component="saddle", requirement="SADDLE-005")
    return [item / length for item in values]


def _flat_project_out(values: Sequence[float], axis: Sequence[float]) -> list[float]:
    coefficient = _flat_dot(values, axis)
    return [item - coefficient * component for item, component in zip(values, axis)]


# Rotation budget for the orthogonal Rayleigh-quotient minimisation. Each cycle costs
# four calculator evaluations (two curvature probes), so a restart costs at most
# 2 + 4 * _ORTHOGONAL_ROTATIONS = 162 evaluations and the search stays bounded. The
# budget is a cap, not a cost: a cycle that buys less than `curvature_tolerance` ends
# the minimisation, so an isotropic complement stops after the first two evaluations.
# Measured on a 15-dof quadratic with a deliberately clustered complement
# (lambda_min = -0.05 against a +0.3 second mode and a +9.2 stiffest mode) the cap sets
# how tight the bound gets: cap 6 -> 4/200 rejected, 12 -> 52/200, 24 -> 155/200,
# 40 -> 186/200, 60 -> 196/200. The reported value is an upper bound at every cap, so a
# too-small cap under-reports the negative curvature and never invents one.
_ORTHOGONAL_ROTATIONS = 40
# Trial rotation angle for the Fourier curvature model. It must not be a multiple of
# pi/2, or cos(2*theta) - 1 vanishes and the model cannot be solved.
_ORTHOGONAL_TRIAL_ANGLE = 0.25


@dataclass(frozen=True, slots=True)
class RelaxationResult:
    state: AtomicState
    steps: int
    evaluations: int
    termination_reason: str


class SteepestDescentMinimizer:
    def __init__(self, calculator: ProcessCalculator, config: Mapping[str, Any]) -> None:
        self.calculator = calculator
        self.config = config

    def minimize(self, request: Mapping[str, Any], *, object_id: str) -> RelaxationResult:
        before = self.calculator.ledger.calculator_reserved
        final_status = "INTERNAL_ERROR"
        termination_reason = "minimizer aborted unexpectedly"
        kind = ("endpoint" if object_id.endswith(":plus") or object_id.endswith(":minus")
                or object_id.startswith("apply-") else "solver")
        try:
            result = self._minimize(request, object_id=object_id)
            final_status = "OK"
            termination_reason = result.termination_reason
            return replace(result, evaluations=self.calculator.ledger.calculator_reserved - before)
        except DomainFailure as exc:
            context = exc.outcome.context
            final_status = exc.outcome.status
            detail = context["details"]
            termination_reason = str(detail.get("termination_reason", exc.outcome.status))
            raise DomainFailure(
                exc.outcome.status, str(exc), component=str(context["component"]),
                requirement=str(context["requirement_id"]), retryable=bool(context["retryable"]),
                state_id=context["state_id"], object_id=context["search_or_event_id"],
                causal_status=exc.outcome.causal_status,
                calculator_evaluations=self.calculator.ledger.calculator_reserved - before,
                iterations=int(detail.get("iterations", 0)), termination_reason=termination_reason) from exc
        except BaseException:
            final_status = "CANCELLED"
            termination_reason = "minimizer interrupted"
            raise
        finally:
            self.calculator.ledger.record_attempt(kind, object_id, before,
                                                  final_status, termination_reason)

    def _minimize(self, request: Mapping[str, Any], *, object_id: str) -> RelaxationResult:
        positions = tuple(tuple(item) for item in request["positions"])
        evaluations = 0
        evaluation = self.calculator.evaluate(_copy_request(request, positions),
                                              component="relaxation", object_id=object_id)
        evaluations += 1
        tolerance = self.config["force_tolerance"]
        for step in range(self.config["max_steps"] + 1):
            maximum = max_movable_force(evaluation.forces, request["movable"])
            if maximum <= tolerance:
                relaxed_request = _copy_request(request, positions)
                provenance = {
                    "method": "steepest-descent-backtracking/1", "steps": step,
                    "calculator_evaluations": evaluations,
                    "calculator_identity": (self.calculator.model_name + "@"
                                            + self.calculator.model_version + ":"
                                            + self.calculator.model_digest),
                    "minimizer_identity": "steepest-descent-backtracking/1",
                    "termination_reason": "force_tolerance",
                    "evaluation_id": evaluation.evaluation_id,
                    "deterministic": evaluation.deterministic,
                }
                state = state_from_relaxation(relaxed_request, evaluation.energy,
                                              evaluation.forces, tolerance, provenance)
                return RelaxationResult(state, step, evaluations, "force_tolerance")
            if step == self.config["max_steps"] or evaluations >= self.config["max_evaluations"]:
                break
            direction = tuple(vscale(force, 1.0 if active else 0.0)
                              for force, active in zip(evaluation.forces, request["movable"]))
            trial_step = 0.05
            accepted: tuple[tuple[Vector, ...], Evaluation] | None = None
            while trial_step >= 1e-8 and evaluations < self.config["max_evaluations"]:
                trial = tuple((position[0] + trial_step * force[0],
                               position[1] + trial_step * force[1],
                               position[2] + trial_step * force[2])
                              for position, force in zip(positions, direction))
                trial_evaluation = self.calculator.evaluate(
                    _copy_request(request, trial), component="relaxation", object_id=object_id)
                evaluations += 1
                if trial_evaluation.energy < evaluation.energy:
                    accepted = trial, trial_evaluation
                    break
                trial_step *= 0.5
            if accepted is None:
                break
            positions, evaluation = accepted
        raise DomainFailure("RELAX_NOT_CONVERGED", "minimizer exhausted its budget",
                            component="relaxation", requirement="RELAX-003", object_id=object_id,
                            calculator_evaluations=evaluations, iterations=step,
                            termination_reason="step/evaluation/backtracking budget exhausted")


@dataclass(frozen=True, slots=True)
class SaddleCandidate:
    positions: tuple[Vector, ...]
    energy: float
    forces: tuple[Vector, ...]
    mode: tuple[Vector, ...]
    curvature: float
    orthogonal_curvatures: tuple[float, ...]
    evidence_level: str
    evaluations: int
    iterations: int
    plus_endpoint: AtomicState
    minus_endpoint: AtomicState
    provenance: Mapping[str, Any]


class DirectionalDimerSearcher:
    def __init__(self, calculator: ProcessCalculator, minimizer: SteepestDescentMinimizer,
                 config: Mapping[str, Any]) -> None:
        self.calculator = calculator
        self.minimizer = minimizer
        self.config = config

    def _initial_mode(self, atom_count: int, movable: Sequence[bool], rng: PhiloxStream,
                      active_hint: Sequence[int] | None) -> list[float]:
        hinted = set(active_hint or ())
        values: list[float] = []
        for atom in range(atom_count):
            active = movable[atom]
            proposal_scale = 2.0 if atom in hinted else 1.0
            for _ in range(3):
                values.append(proposal_scale * (2.0 * rng.uniform() - 1.0) if active else 0.0)
        return _flat_normalize(values)

    def _offset(self, positions: Sequence[Vector], mode: Sequence[float], distance: float) -> tuple[Vector, ...]:
        vectors = _unflatten(mode)
        return tuple((position[0] + distance * vector[0],
                      position[1] + distance * vector[1],
                      position[2] + distance * vector[2])
                     for position, vector in zip(positions, vectors))

    def search(self, origin: AtomicState, rng: PhiloxStream, *, search_id: str,
               active_hint: Sequence[int] | None = None) -> SaddleCandidate:
        before = self.calculator.ledger.calculator_reserved
        final_status = "INTERNAL_ERROR"
        termination_reason = "dimer search aborted unexpectedly"
        try:
            candidate = self._search(origin, rng, search_id=search_id, active_hint=active_hint)
            evaluations = self.calculator.ledger.calculator_reserved - before
            provenance = dict(candidate.provenance)
            provenance["calculator_evaluations"] = evaluations
            final_status = "OK"
            termination_reason = str(provenance["termination_reason"])
            return replace(candidate, evaluations=evaluations, provenance=provenance)
        except DomainFailure as exc:
            context = exc.outcome.context
            final_status = exc.outcome.status
            detail = context["details"]
            termination_reason = str(detail.get("termination_reason", exc.outcome.status))
            raise DomainFailure(
                exc.outcome.status, str(exc), component=str(context["component"]),
                requirement=str(context["requirement_id"]), retryable=bool(context["retryable"]),
                state_id=context["state_id"], object_id=context["search_or_event_id"],
                causal_status=exc.outcome.causal_status,
                calculator_evaluations=self.calculator.ledger.calculator_reserved - before,
                iterations=int(detail.get("iterations", 0)), termination_reason=termination_reason) from exc
        except BaseException:
            final_status = "CANCELLED"
            termination_reason = "dimer search interrupted"
            raise
        finally:
            self.calculator.ledger.record_attempt("solver", search_id, before,
                                                  final_status, termination_reason)

    def _search(self, origin: AtomicState, rng: PhiloxStream, *, search_id: str,
                active_hint: Sequence[int] | None = None) -> SaddleCandidate:
        request = _state_request(origin)
        mode = self._initial_mode(len(origin.atom_ids), origin.movable, rng, active_hint)
        positions = self._offset(origin.positions, mode, self.config["endpoint_displacement"])
        separation = 1e-3
        evaluations = 0
        center_eval: Evaluation | None = None
        curvature = math.inf
        for iteration in range(1, self.config["max_iterations"] + 1):
            center_eval = self.calculator.evaluate(_copy_request(request, positions),
                                                  component="saddle", object_id=search_id)
            plus_eval = self.calculator.evaluate(
                _copy_request(request, self._offset(positions, mode, separation)),
                component="saddle", object_id=search_id)
            minus_eval = self.calculator.evaluate(
                _copy_request(request, self._offset(positions, mode, -separation)),
                component="saddle", object_id=search_id)
            evaluations += 3
            plus_force = _flatten(plus_eval.forces, origin.movable)
            minus_force = _flatten(minus_eval.forces, origin.movable)
            hessian_mode = [(minus - plus) / (2.0 * separation)
                            for minus, plus in zip(minus_force, plus_force)]
            curvature = _flat_dot(mode, hessian_mode)
            physical_force = _flatten(center_eval.forces, origin.movable)
            if (max_movable_force(center_eval.forces, origin.movable) <= self.config["force_tolerance"]
                    and curvature < -self.config["curvature_tolerance"]):
                before_orthogonal = self.calculator.ledger.calculator_reserved
                orthogonal = self._orthogonal_curvatures(request, positions, mode, origin.movable,
                                                         rng, search_id)
                # The orthogonal minimisation runs a variable number of rotation cycles,
                # so the audited count is the ledger delta, never a per-direction constant.
                evaluations += self.calculator.ledger.calculator_reserved - before_orthogonal
                if all(item >= -self.config["curvature_tolerance"] for item in orthogonal):
                    break
            rotational = [component - curvature * axis
                          for component, axis in zip(hessian_mode, mode)]
            rotated = [axis - 0.1 * component
                       for axis, component in zip(mode, rotational)]
            mode = _flat_normalize(rotated)
            along = _flat_dot(physical_force, mode)
            climbing_force = [component - 2.0 * along * axis
                              for component, axis in zip(physical_force, mode)]
            step_vectors = _unflatten([0.01 * item for item in climbing_force])
            positions = tuple((position[0] + step[0], position[1] + step[1], position[2] + step[2])
                              for position, step in zip(positions, step_vectors))
        else:
            raise DomainFailure("SADDLE_NOT_FOUND", "dimer search exhausted iterations",
                                component="saddle", requirement="SADDLE-007", object_id=search_id,
                                calculator_evaluations=evaluations,
                                iterations=self.config["max_iterations"],
                                termination_reason="maximum dimer iterations exhausted")
        if center_eval is None:
            raise DomainFailure("INTERNAL_ERROR", "dimer ended without a center evaluation",
                                component="saddle", requirement="SADDLE-007", object_id=search_id)
        mode_vectors = _unflatten(mode)
        endpoint_distance = self.config["endpoint_displacement"]
        plus_request = _copy_request(request, self._offset(positions, mode, endpoint_distance))
        minus_request = _copy_request(request, self._offset(positions, mode, -endpoint_distance))
        plus_endpoint = self.minimizer.minimize(plus_request, object_id=search_id + ":plus").state
        minus_endpoint = self.minimizer.minimize(minus_request, object_id=search_id + ":minus").state
        provenance = {
            "method": "directional-dimer/1", "search_id": search_id,
            "active_hint_role": "initial-proposal-only", "iterations": iteration,
            "calculator_evaluations": evaluations,
            "termination_reason": "full-force-curvature-orthogonal-gates",
        }
        return SaddleCandidate(
            positions, center_eval.energy, center_eval.forces, mode_vectors, curvature,
            tuple(orthogonal), "DIRECTIONAL", evaluations, iteration,
            plus_endpoint, minus_endpoint, provenance,
        )

    def _curvature_along(self, request: Mapping[str, Any], positions: Sequence[Vector],
                         direction: Sequence[float], movable: Sequence[bool], search_id: str,
                         separation: float) -> tuple[float, list[float]]:
        """Rayleigh quotient `d^T H d` and the product `H d`, by central differences."""
        plus = self.calculator.evaluate(
            _copy_request(request, self._offset(positions, direction, separation)),
            component="saddle", object_id=search_id)
        minus = self.calculator.evaluate(
            _copy_request(request, self._offset(positions, direction, -separation)),
            component="saddle", object_id=search_id)
        plus_force = _flatten(plus.forces, movable)
        minus_force = _flatten(minus.forces, movable)
        product = [(left - right) / (2.0 * separation)
                   for left, right in zip(minus_force, plus_force)]
        return _flat_dot(direction, product), product

    def _minimum_orthogonal_curvature(self, request: Mapping[str, Any],
                                      positions: Sequence[Vector], start: Sequence[float],
                                      mode: Sequence[float], movable: Sequence[bool],
                                      search_id: str) -> float:
        """Minimise the Rayleigh quotient over the complement of `mode`.

        A second dimer-style rotation, confined to `d` perpendicular to `mode`: at each
        cycle the projected residual `P(Hd) - lambda d` gives the steepest ascent
        direction inside the complement, one trial rotation fixes the Fourier model
        `C(t) = a0/2 + a1 cos 2t + b1 sin 2t` (with `b1 = |residual|` analytically), and
        the model's minimising angle is taken. The returned number is therefore an upper
        BOUND on the lowest curvature of the complement, not one sample of it: a single
        random direction is a positively weighted average of the whole orthogonal
        spectrum and is dominated by the stiff modes, so it cannot distinguish index 1
        from index 7.
        """
        separation = 1e-3
        direction = _flat_normalize(_flat_project_out(start, mode))
        curvature, product = self._curvature_along(request, positions, direction, movable,
                                                   search_id, separation)
        lowest = curvature
        tolerance = self.config["curvature_tolerance"]
        for _cycle in range(_ORTHOGONAL_ROTATIONS):
            # `residual` is the component of `H d` perpendicular to both `d` and `mode`;
            # it is the steepest-ASCENT direction of the Rayleigh quotient on the unit
            # sphere of the complement, and it vanishes exactly at an eigendirection.
            residual = _flat_project_out(
                [item - curvature * axis for item, axis in zip(product, direction)], mode)
            gradient = math.sqrt(_flat_dot(residual, residual))
            if not math.isfinite(gradient) or gradient <= tolerance:
                break
            ascent = [item / gradient for item in residual]
            angle = _ORTHOGONAL_TRIAL_ANGLE
            trial = _flat_normalize(_flat_project_out(
                [axis * math.cos(angle) + step * math.sin(angle)
                 for axis, step in zip(direction, ascent)], mode))
            trial_curvature, _trial_product = self._curvature_along(
                request, positions, trial, movable, search_id, separation)
            lowest = min(lowest, trial_curvature)
            # In the plane spanned by (direction, ascent) the curvature is exactly
            # C(t) = a0/2 + a1 cos 2t + b1 sin 2t with C(0) = a0/2 + a1 = curvature and
            # dC/dt(0) = 2 * |residual|, so b1 is known analytically and the single trial
            # angle fixes a1; the model's minimising angle is then taken directly.
            denominator = math.cos(2.0 * angle) - 1.0
            a1 = (trial_curvature - curvature - gradient * math.sin(2.0 * angle)) / denominator
            if not math.isfinite(a1):
                break
            optimum = 0.5 * math.atan2(gradient, a1) + 0.5 * math.pi
            rotated = _flat_normalize(_flat_project_out(
                [axis * math.cos(optimum) + step * math.sin(optimum)
                 for axis, step in zip(direction, ascent)], mode))
            rotated_curvature, rotated_product = self._curvature_along(
                request, positions, rotated, movable, search_id, separation)
            lowest = min(lowest, rotated_curvature)
            # Stop as soon as a whole cycle buys less than the curvature tolerance. An
            # isotropic complement (every orthogonal curvature equal) leaves at the
            # `gradient <= tolerance` test above after two evaluations, so the rotation
            # budget is spent only while the bound is still moving.
            if rotated_curvature >= curvature - tolerance:
                break
            direction, curvature, product = rotated, rotated_curvature, rotated_product
        return lowest

    def _orthogonal_curvatures(self, request: Mapping[str, Any], positions: Sequence[Vector],
                               mode: Sequence[float], movable: Sequence[bool], rng: PhiloxStream,
                               search_id: str) -> list[float]:
        """One bound per configured restart, each from an independent random start.

        `SADDLE-005` requires "negative curvature along the reported mode plus
        nonnegative sampled orthogonal curvatures" without a Hessian. Each entry here is
        the minimum of the Rayleigh quotient found from its restart, which is a strictly
        stronger statement than the value at the random start itself, so the
        `DIRECTIONAL` evidence level of `E2-EVENT-002` is not over-claimed.
        """
        results: list[float] = []
        for _ in range(self.config["orthogonal_directions"]):
            trial = [2.0 * rng.uniform() - 1.0 if active else 0.0
                     for active in movable for _coordinate in range(3)]
            results.append(self._minimum_orthogonal_curvature(request, positions, trial,
                                                              mode, movable, search_id))
        return results
