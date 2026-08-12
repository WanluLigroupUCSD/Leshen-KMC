# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Executable fixture builders shared by the test modules.

Importing this module also puts `src/` on `sys.path`, so
`python3 -m unittest discover -s tests -p 'test_*.py'` works from the project root
without an install step. Every test module imports this module first.

Nothing here is a golden artifact. The shared Errata-2 golden corpus lives at
`../spark-atomistic-rs/tests/corpus/` and is referenced, never copied
(`E2-PAR-001` requires both backends to consume the same canonical corpus).
"""

from __future__ import annotations

import contextlib
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

CORPUS = ROOT.parent / "spark-atomistic-rs" / "tests" / "corpus"

from spark_atomistic.calculator import Evaluation  # noqa: E402
from spark_atomistic.canonical import deep_freeze, deep_thaw, digest  # noqa: E402
from spark_atomistic.catalog import Event  # noqa: E402
from spark_atomistic.engine import ReferenceEngine  # noqa: E402
from spark_atomistic.errors import DomainFailure  # noqa: E402
from spark_atomistic.kinetics import RateTable, build_rate_table, propose_serial_step  # noqa: E402
from spark_atomistic.model import AtomicState, state_from_relaxation, validate_model  # noqa: E402
from spark_atomistic.resources import ResourceLedger  # noqa: E402
from spark_atomistic.rng import PhiloxStream, derive_saddle_stream  # noqa: E402
from spark_atomistic.solvers import SaddleCandidate  # noqa: E402


CELL = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
CHECKPOINT_MODEL_DIGEST = "sha256:test-model"
# `corpus/mock_calculator.py` advertises this digest for the analytic double well.
DOUBLE_WELL_MODEL_DIGEST = "7c799e3c0c25eb952d433430027d3d73de8d9f8f3d06064b6374f4b6eab4dd47"


# --------------------------------------------------------------------------------------
# Analytic double well, identical to corpus/mock_calculator.py
# --------------------------------------------------------------------------------------

def analytic_energy_forces(positions: Sequence[Sequence[float]]):
    """V = sum_i (x_i^2-1)^2 + y_i^2 + z_i^2 in eV with positions in Angstrom.

    Minima at x = +-1 (E = 0); first-order saddle at x = 0 (E = 1.0 eV);
    forward and reverse barriers are both exactly 1.000 eV.
    """
    energy = 0.0
    forces = []
    for x, y, z in positions:
        energy += (x * x - 1.0) ** 2 + y * y + z * z
        forces.append((-4.0 * x * (x * x - 1.0), -2.0 * y, -2.0 * z))
    return energy, tuple(forces)


def analytic_hessian_eigenvalues(x: float) -> tuple[float, float, float]:
    """Exact eigenvalues of the double-well Hessian for one atom at (x, 0, 0)."""
    return (12.0 * x * x - 4.0, 2.0, 2.0)


class StubCalculator:
    """In-process test double with the `ProcessCalculator` surface used by the solvers.

    The subprocess transport is deliberately replaced so that a solver experiment
    measures the SOLVER. Overflow is mapped the way the real transport maps it:
    `corpus/mock_calculator.py` raises `OverflowError`, exits nonzero, and
    `ProcessCalculator` turns a nonzero exit into `CALCULATOR_FAILURE` (`CALC-006`).
    """

    model_name = "analytic-double-well"
    model_version = "1"
    model_digest = DOUBLE_WELL_MODEL_DIGEST

    def __init__(self, ledger: ResourceLedger) -> None:
        self.ledger = ledger
        self.centres: list[tuple[float, float, float]] = []
        self.overflow_at: int | None = None
        self._scope: tuple[str, int, int] | None = None

    @contextlib.contextmanager
    def evaluation_scope(self, scope_id: str, limit: int):
        """Same scoping contract as `ProcessCalculator`, so a scoped caller is tested."""
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
        if self._scope is not None:
            scope_id, limit, count = self._scope
            if count >= limit:
                raise DomainFailure("RESOURCE_LIMIT", "scoped calculator evaluation limit reached",
                                    component="resources", requirement="RES-002",
                                    object_id=scope_id)
            self._scope = (scope_id, limit, count + 1)
        self.ledger.reserve_evaluation()
        self.ledger.complete(True)
        positions = [tuple(item) for item in state_request["positions"]]
        try:
            energy, forces = analytic_energy_forces(positions)
        except OverflowError as exc:
            self.overflow_at = self.ledger.calculator_reserved
            raise DomainFailure("CALCULATOR_FAILURE", "calculator process returned nonzero",
                                component=component, requirement="CALC-006",
                                object_id=object_id) from exc
        if not math.isfinite(energy) or any(not math.isfinite(item)
                                            for row in forces for item in row):
            self.overflow_at = self.ledger.calculator_reserved
            raise DomainFailure("NONFINITE_RESULT", "calculator energy is nonfinite",
                                component=component, requirement="STATE-008",
                                object_id=object_id)
        if component == "saddle":
            self.centres.append(positions[0])
        return Evaluation(energy, forces, "eval-%d" % self.ledger.calculator_reserved,
                          True, "fixture-request-digest")


def double_well_state(x: float, *, force_tolerance: float = 10.0) -> AtomicState:
    """One movable H atom at (x, 0, 0) in a 10 A non-periodic cubic cell."""
    request = {
        "schema": "spark-atomistic-model/1", "atom_ids": ["a0"], "species": ["H"],
        "positions": [[x, 0.0, 0.0]], "cell": CELL, "pbc": [False, False, False],
        "movable": [True], "constraints": {"kind": "fixed-mask"},
        "charge": 0.0, "spin": 0.0, "calculator_model_digest": DOUBLE_WELL_MODEL_DIGEST,
    }
    energy, forces = analytic_energy_forces(request["positions"])
    return state_from_relaxation(request, energy, forces, force_tolerance,
                                 {"calculator_evaluations": 1, "steps": 0,
                                  "calculator_identity": "fixture@1",
                                  "minimizer_identity": "fixture-minimizer/1",
                                  "termination_reason": "force_tolerance"})


def double_well_chain_state(coordinates: Sequence[float], *,
                            force_tolerance: float = 1e-4) -> AtomicState:
    """N movable H atoms at (x_i, 0, 0) in the same analytic double well.

    Every `x_i = +-1` configuration is a true minimum with zero force and energy 0, and
    configurations that are not related by whole-cell translation plus same-species
    permutation are distinct committed states.
    """
    request = {
        "schema": "spark-atomistic-model/1",
        "atom_ids": ["a%d" % index for index in range(len(coordinates))],
        "species": ["H"] * len(coordinates),
        "positions": [[float(value), 0.0, 0.0] for value in coordinates], "cell": CELL,
        "pbc": [False, False, False], "movable": [True] * len(coordinates),
        "constraints": {"kind": "fixed-mask"}, "charge": 0.0, "spin": 0.0,
        "calculator_model_digest": DOUBLE_WELL_MODEL_DIGEST,
    }
    energy, forces = analytic_energy_forces(request["positions"])
    return state_from_relaxation(request, energy, forces, force_tolerance,
                                 {"calculator_evaluations": 1, "steps": 0,
                                  "calculator_identity": "fixture@1",
                                  "minimizer_identity": "fixture-minimizer/1",
                                  "termination_reason": "force_tolerance"})


# --------------------------------------------------------------------------------------
# Exact quadratic model with a prescribed Hessian spectrum
# --------------------------------------------------------------------------------------

def hessian_from_spectrum(spectrum: Sequence[float]):
    """`(H, first eigenvector)` for `H = Q diag(spectrum) Q^T`, `Q` a Householder frame.

    The eigenbasis is deliberately NOT the coordinate basis: an axis-aligned Hessian
    would let a direction-handling bug pass unnoticed.
    """
    dimension = len(spectrum)
    raw = [math.sin(1.0 + 0.7 * index) for index in range(dimension)]
    length = math.sqrt(sum(item * item for item in raw))
    unit = [item / length for item in raw]
    frame = [[(1.0 if row == column else 0.0) - 2.0 * unit[row] * unit[column]
              for column in range(dimension)] for row in range(dimension)]
    matrix = [[sum(frame[row][index] * spectrum[index] * frame[column][index]
                   for index in range(dimension))
               for column in range(dimension)] for row in range(dimension)]
    return matrix, [frame[row][0] for row in range(dimension)]


class QuadraticCalculator:
    """`V(x) = 0.5 (x-c)^T H (x-c)`; exact forces, so central differences recover `H`."""

    model_name = "quadratic"
    model_version = "1"
    model_digest = "sha256:quadratic"

    def __init__(self, ledger: ResourceLedger, hessian, centre: Sequence[float]) -> None:
        self.ledger = ledger
        self.hessian = hessian
        self.centre = list(centre)

    def evaluate(self, state_request: Mapping[str, Any], *, component: str,
                 object_id: str | None = None) -> Evaluation:
        self.ledger.reserve_evaluation()
        self.ledger.complete(True)
        flat = [item for row in state_request["positions"] for item in row]
        delta = [left - right for left, right in zip(flat, self.centre)]
        product = [sum(self.hessian[row][column] * delta[column]
                       for column in range(len(delta))) for row in range(len(delta))]
        energy = 0.5 * sum(left * right for left, right in zip(delta, product))
        forces = tuple(tuple(-product[3 * atom + axis] for axis in range(3))
                       for atom in range(len(delta) // 3))
        return Evaluation(energy, forces, "eval-%d" % self.ledger.calculator_reserved,
                          True, "fixture-request-digest")


def quadratic_saddle_probe(spectrum: Sequence[float], *, atoms: int):
    """`(request, positions, mode, calculator factory)` at the stationary point of `H`."""
    hessian, mode = hessian_from_spectrum(spectrum)
    centre = [0.0] * (3 * atoms)
    positions = tuple((0.0, 0.0, 0.0) for _ in range(atoms))
    request = {"schema": "spark-atomistic-model/1",
               "atom_ids": ["a%d" % index for index in range(atoms)],
               "species": ["H"] * atoms, "positions": [list(item) for item in positions],
               "cell": CELL, "pbc": [False, False, False], "movable": [True] * atoms,
               "constraints": {"kind": "fixed-mask"}, "charge": 0.0, "spin": 0.0,
               "calculator_model_digest": "sha256:quadratic"}
    return request, positions, mode, (lambda ledger: QuadraticCalculator(ledger, hessian, centre))


def unlimited_ledger(total_evaluations: int = 1_000_000) -> ResourceLedger:
    return ResourceLedger(total_evaluations, 1e9, 1 << 62, 1 << 40, 1000, 0.0)


# --------------------------------------------------------------------------------------
# Checkpoint fixture: two committed events, one KMC step, then catalog growth
# --------------------------------------------------------------------------------------

_ZERO_FORCES_2 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

# (destination energy, saddle x of atom 1, saddle energy, unstable direction)
# Energies are chosen so that the event committed AFTER the step dominates the rate
# sum: its inclusion changes both the selected event and the time increment, which is
# what gives the E2-KMC-005 fixture its discriminating power.
_SADDLE_PLAN = (
    (-0.8, 1.5, 1.25, -0.20, ((1.0, 0.0, 0.0), (0.0, 0.0, 0.0))),
    (-0.9, 2.0, 1.60, -0.30, ((0.0, 1.0, 0.0), (0.0, 0.0, 0.0))),
    (-0.7, 2.5, 1.90, -0.35, ((0.0, 0.0, 1.0), (0.0, 0.0, 0.0))),
)


@dataclass
class GrownCatalogCheckpoint:
    engine: ReferenceEngine
    config: Mapping[str, Any]
    origin: AtomicState
    historical: RateTable
    current: RateTable
    selection: Any
    checkpoint_path: str
    current_selected_event_id: str
    current_selected_rate: float
    current_time_increment: float

    def restore_arguments(self) -> dict[str, Any]:
        engine, config = self.engine, self.config
        return {"expected_config_digest": engine.config_digest,
                "expected_model_digest": engine.model_digest,
                "expected_tolerance_digest": engine.tolerance_digest,
                "expected_identity_digest": engine.identity_digest,
                "kinetics": config["kinetics"], "saddle_config": config["saddle_search"],
                "relaxation_config": config["relaxation"],
                "discovery_config": config["discovery"],
                "resource_config": config["resources"],
                "maximum_events": config["resources"]["catalog_events"]}


def checkpoint_model(directory: str | Path, **output_overrides: Any) -> Mapping[str, Any]:
    """A validated model derived from the shared minimal model.

    Three values are raised relative to `e2_minimal_model.json` because that fixture
    cannot host three searches from one state, and because its
    `resources.resident_memory_bytes` (1 MB) is below the resident set of any Python
    process, so `ResourceLedger.check_wall_time` would raise `RESOURCE_LIMIT`
    (`RES-002`) before the first calculator evaluation. See TEST_SPEC.md.
    """
    model = json.loads((CORPUS / "e2_minimal_model.json").read_bytes())
    model["system"] = {
        "atom_ids": ["a0", "a1"], "species": ["H", "H"],
        "positions": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], "cell": CELL,
        "pbc": [False, False, False], "movable": [True, True],
        "constraints": {"kind": "fixed-mask"}, "charge": 0.0, "spin": 0.0,
        "calculator_model_digest": CHECKPOINT_MODEL_DIGEST,
    }
    model["discovery"]["maximum_attempts"] = 3
    model["resources"]["saddle_attempts_per_state"] = 3
    model["resources"]["resident_memory_bytes"] = 1 << 40
    model["resources"]["output_bytes"] = 1 << 22
    directory = Path(directory)
    model["output"]["checkpoint_path"] = str(directory / "checkpoint.json")
    model["output"]["summary_path"] = str(directory / "summary.json")
    model["output"]["trajectory_path"] = str(directory / "trajectory.json")
    model["output"].update(output_overrides)
    return validate_model(model)


def _committed_state(positions, energy: float, tolerance: float) -> AtomicState:
    request = {"schema": "spark-atomistic-model/1", "atom_ids": ["a0", "a1"],
               "species": ["H", "H"], "positions": positions, "cell": CELL,
               "pbc": [False, False, False], "movable": [True, True],
               "constraints": {"kind": "fixed-mask"}, "charge": 0.0, "spin": 0.0,
               "calculator_model_digest": CHECKPOINT_MODEL_DIGEST}
    return state_from_relaxation(request, energy, _ZERO_FORCES_2, tolerance,
                                 {"calculator_evaluations": 1, "steps": 0,
                                  "calculator_identity": "fixture@1:" + CHECKPOINT_MODEL_DIGEST,
                                  "minimizer_identity": "fixture-minimizer/1",
                                  "termination_reason": "force_tolerance"})


def build_grown_catalog_checkpoint(directory: str | Path) -> GrownCatalogCheckpoint:
    """Commit two events, take one KMC step, then grow the catalog, then checkpoint.

    This is the `kmc-historical-snapshot-growth` fixture required by `E2-PAR-002`
    item 9. The step-1 rate-table snapshot is the two-event table that existed when
    the step was taken; the checkpointed catalog holds three selectable events from
    the same origin. `E2-KMC-005` requires replay to use the snapshot.
    """
    config = checkpoint_model(directory)
    engine = ReferenceEngine(config, extension={"calculator_command": ["/bin/true"]})
    kinetics = config["kinetics"]
    tolerance = config["relaxation"]["force_tolerance"]

    origin = _committed_state([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], -1.0, tolerance)
    engine.catalog.add_state(origin, kinetics)
    engine.current_state_id = engine.initial_state_id = origin.state_id
    record = engine.catalog.discovery_record(origin.state_id, "strict", config["discovery"])

    def commit(index: int) -> Event:
        energy, destination_x, saddle_x, saddle_energy, mode = _SADDLE_PLAN[index]
        destination = _committed_state([[0.0, 0.0, 0.0], [destination_x, 0.0, 0.0]],
                                       energy, tolerance)
        class_stream = derive_saddle_stream(kinetics["run_seed"], origin.state_id,
                                            "class-selection", index)
        class_stream.uniform()
        search_id = "search:" + digest({"run_seed": kinetics["run_seed"],
                                        "search_class": "global", "search_index": index,
                                        "state_id": origin.state_id})
        stream = derive_saddle_stream(kinetics["run_seed"], origin.state_id, "global", index)
        stream.uniform()
        candidate = SaddleCandidate(
            positions=((0.0, 0.0, 0.0), (saddle_x, 0.0, 0.0)), energy=saddle_energy,
            forces=((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), mode=mode, curvature=-4.0,
            orthogonal_curvatures=(2.0,), evidence_level="DIRECTIONAL", evaluations=7,
            iterations=3, plus_endpoint=origin, minus_endpoint=destination,
            provenance={"method": "directional-dimer/1", "search_id": search_id,
                        "termination_reason": "full-force-curvature-orthogonal-gates"})
        _outcome, event, _destination = engine.catalog.validate_candidate(
            origin, candidate, kinetics,
            {"rng_substream_digest": digest(stream.checkpoint()), "search_class": "global",
             "search_id": search_id, "search_index": index})
        engine.substreams["class-selection:" + search_id] = class_stream.clone()
        engine.substreams[search_id] = stream.clone()
        record.attempts += 1
        record.successes += 1
        record.evaluations += 2
        if event is not None:
            record.event_log_rates[event.event_id] = event.log_rate
        assert event is not None
        return event

    commit(0)
    commit(1)

    historical = build_rate_table(tuple(engine.catalog.events.values()), origin.state_id)
    selection = propose_serial_step(historical, engine.trajectory_rng)
    engine.trajectory_log.append(deep_freeze({
        "checkpoint_sequence": engine.checkpoint_sequence, "log_sequence": 1,
        "post_state_id": selection.destination_state_id, "pre_state_id": origin.state_id,
        "rate_table_snapshot": selection.rate_table_snapshot,
        "selected_event_id": selection.event_id,
        "selected_rate_per_s": selection.selected_rate,
        "selection_uniform": selection.selection_uniform, "step_index": 1,
        "time_increment_s": selection.delta_time, "time_uniform": selection.time_uniform,
        "total_rate_per_s": selection.total_rate}))
    engine.current_state_id = selection.destination_state_id
    engine.simulation_time = selection.delta_time
    engine.trajectory_rng.commit_from(selection.rng_after)
    engine.step_index = 1
    engine.log_sequence = 1
    engine.commit_generation = 1

    commit(2)  # catalog growth, strictly after the step was committed

    current = build_rate_table(tuple(engine.catalog.events.values()), origin.state_id)
    threshold = selection.selection_uniform * current.total_rate
    cumulative = 0.0
    index = len(current.rates) - 1
    for position, rate in enumerate(current.rates):
        cumulative += rate
        if cumulative > threshold:
            index = position
            break

    engine.ledger.per_state_saddle_attempts[origin.state_id] = record.attempts
    engine.last_status = "OK"
    engine.write_checkpoint()

    return GrownCatalogCheckpoint(
        engine=engine, config=config, origin=origin, historical=historical,
        current=current, selection=selection,
        checkpoint_path=config["output"]["checkpoint_path"],
        current_selected_event_id=current.event_ids[index],
        current_selected_rate=current.rates[index],
        current_time_increment=-math.log(selection.time_uniform) / current.total_rate)


def clone_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    """A mutable deep copy, so a corruption fixture cannot leak into its baseline."""
    return deep_thaw(payload)
