# Independently authored from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1 + ERRATA_2.
# Base: 8e21e556838d997aa257db1a194c3ed88018ad3d77dcf1d811a82f7194f49c84
# Erratum 1: 52c9c8bc9e6839ec04709da99aa7cca39963f4efd6fc420a865d6bb0b911ec40
# Erratum 2: eecc5f13c72d90ae2b5c6ff99385cb51f69c0d0500b563eb3a8df2e7d9863995
"""Triclinic periodic geometry and permutation-aware state matching."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence

from .errors import DomainFailure


Vector = tuple[float, float, float]
Cell = tuple[Vector, Vector, Vector]


def vadd(a: Vector, b: Vector) -> Vector:
    return a[0] + b[0], a[1] + b[1], a[2] + b[2]


def vsub(a: Vector, b: Vector) -> Vector:
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]


def vscale(a: Vector, scalar: float) -> Vector:
    return a[0] * scalar, a[1] * scalar, a[2] * scalar


def dot(a: Vector, b: Vector) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm2(a: Vector) -> float:
    return dot(a, a)


def norm(a: Vector) -> float:
    return math.sqrt(norm2(a))


def _basis_qr(cell: Cell, pbc: tuple[bool, bool, bool]) -> tuple[list[Vector], list[list[float]], list[Vector]]:
    basis = [cell[index] for index, periodic in enumerate(pbc) if periodic]
    count = len(basis)
    q: list[Vector] = []
    r = [[0.0] * count for _ in range(count)]
    for column, vector in enumerate(basis):
        residual = vector
        for row, unit in enumerate(q):
            coefficient = dot(unit, vector)
            r[row][column] = coefficient
            residual = vsub(residual, vscale(unit, coefficient))
        diagonal = norm(residual)
        scale = max(1.0, *(norm(item) for item in basis))
        if not math.isfinite(diagonal) or diagonal <= 1e-14 * scale:
            raise DomainFailure("INVALID_STATE", "periodic cell basis is singular",
                                component="geometry", requirement="STATE-001")
        r[column][column] = diagonal
        q.append(vscale(residual, 1.0 / diagonal))
    return q, r, basis


def closest_periodic_displacement(delta: Vector, cell: Cell,
                                  pbc: tuple[bool, bool, bool]) -> tuple[Vector, tuple[int, ...]]:
    """Exact p<=3 closest-lattice-vector search by QR sphere decoding."""
    q_basis, upper, lattice_basis = _basis_qr(cell, pbc)
    dimension = len(lattice_basis)
    if dimension == 0:
        return delta, ()
    projected = [dot(unit, delta) for unit in q_basis]

    babai = [0] * dimension
    for row in range(dimension - 1, -1, -1):
        tail = sum(upper[row][column] * babai[column]
                   for column in range(row + 1, dimension))
        babai[row] = math.floor((projected[row] - tail) / upper[row][row] + 0.5)

    def projected_error(indices: Sequence[int]) -> float:
        return sum((projected[row] - sum(upper[row][column] * indices[column]
                                         for column in range(row, dimension))) ** 2
                   for row in range(dimension))

    best_indices = tuple(babai)
    best = projected_error(best_indices)
    working = [0] * dimension

    def visit(row: int, partial: float) -> None:
        nonlocal best, best_indices
        if row < 0:
            candidate = tuple(working)
            tolerance = 32.0 * math.ulp(max(1.0, best, partial))
            if partial < best - tolerance or (abs(partial - best) <= tolerance and candidate < best_indices):
                best = partial
                best_indices = candidate
            return
        tail = sum(upper[row][column] * working[column]
                   for column in range(row + 1, dimension))
        center = (projected[row] - tail) / upper[row][row]
        remaining = max(0.0, best - partial)
        radius = math.sqrt(remaining) / abs(upper[row][row])
        low = math.ceil(center - radius - 8.0 * math.ulp(max(1.0, abs(center), radius)))
        high = math.floor(center + radius + 8.0 * math.ulp(max(1.0, abs(center), radius)))
        candidates = list(range(low, high + 1))
        candidates.sort(key=lambda integer: (abs(integer - center), integer))
        for integer in candidates:
            residual = projected[row] - tail - upper[row][row] * integer
            next_partial = partial + residual * residual
            if next_partial <= best + 32.0 * math.ulp(max(1.0, best)):
                working[row] = integer
                visit(row - 1, next_partial)

    visit(dimension - 1, 0.0)
    lattice = (0.0, 0.0, 0.0)
    for coefficient, vector in zip(best_indices, lattice_basis):
        lattice = vadd(lattice, vscale(vector, coefficient))
    return vsub(delta, lattice), best_indices


def _hungarian(cost: list[list[float]]) -> tuple[int, ...]:
    """Deterministic minimum-cost square assignment."""
    size = len(cost)
    if any(len(row) != size for row in cost):
        raise DomainFailure("INVALID_STATE", "assignment matrix is not square",
                            component="geometry", requirement="STATE-005")
    u = [0.0] * (size + 1)
    v = [0.0] * (size + 1)
    p = [0] * (size + 1)
    way = [0] * (size + 1)
    for i in range(1, size + 1):
        p[0] = i
        minimum = [math.inf] * (size + 1)
        used = [False] * (size + 1)
        column0 = 0
        while True:
            used[column0] = True
            row0 = p[column0]
            delta = math.inf
            column1 = 0
            for column in range(1, size + 1):
                if not used[column]:
                    current = cost[row0 - 1][column - 1] - u[row0] - v[column]
                    if current < minimum[column]:
                        minimum[column] = current
                        way[column] = column0
                    if minimum[column] < delta or (minimum[column] == delta and column < column1):
                        delta = minimum[column]
                        column1 = column
            for column in range(size + 1):
                if used[column]:
                    u[p[column]] += delta
                    v[column] -= delta
                else:
                    minimum[column] -= delta
            column0 = column1
            if p[column0] == 0:
                break
        while True:
            column1 = way[column0]
            p[column0] = p[column1]
            column0 = column1
            if column0 == 0:
                break
    assignment = [0] * size
    for column in range(1, size + 1):
        assignment[p[column] - 1] = column - 1
    return tuple(assignment)


@dataclass(frozen=True, slots=True)
class MatchReport:
    equal: bool
    atom_mapping: tuple[int, ...]
    rms_displacement: float
    maximum_displacement: float
    energy_difference: float
    translation: Vector


def match_positions(species_a: Sequence[str], positions_a: Sequence[Vector], energy_a: float,
                    species_b: Sequence[str], positions_b: Sequence[Vector], energy_b: float,
                    cell: Cell, pbc: tuple[bool, bool, bool], *, rms_tolerance: float,
                    max_tolerance: float, energy_tolerance_per_atom: float) -> MatchReport:
    """Jointly refine translation, periodic images, and same-species assignment."""
    count = len(positions_a)
    if count == 0 or len(positions_b) != count or len(species_a) != count or len(species_b) != count:
        return MatchReport(False, (), math.inf, math.inf, math.inf, (0.0, 0.0, 0.0))
    groups_a = {element: [i for i, item in enumerate(species_a) if item == element]
                for element in sorted(set(species_a))}
    groups_b = {element: [i for i, item in enumerate(species_b) if item == element]
                for element in sorted(set(species_b))}
    if {key: len(value) for key, value in groups_a.items()} != {key: len(value) for key, value in groups_b.items()}:
        return MatchReport(False, (), math.inf, math.inf, math.inf, (0.0, 0.0, 0.0))

    seeds: list[Vector] = [(0.0, 0.0, 0.0)]
    rarest = min(groups_a, key=lambda element: (len(groups_a[element]), element))
    for source in groups_a[rarest]:
        for target in groups_b[rarest]:
            seed, _ = closest_periodic_displacement(vsub(positions_a[source], positions_b[target]), cell, pbc)
            seeds.append(seed)

    best: tuple[float, float, tuple[int, ...], Vector, list[Vector]] | None = None
    for initial in seeds:
        translation = initial
        previous_mapping: tuple[int, ...] | None = None
        mapping = [-1] * count
        displacements: list[Vector] = []
        for _ in range(32):
            for element in sorted(groups_a):
                left = groups_a[element]
                right = groups_b[element]
                costs: list[list[float]] = []
                for source in left:
                    row: list[float] = []
                    for target in right:
                        displacement, _ = closest_periodic_displacement(
                            vsub(positions_a[source], vadd(positions_b[target], translation)), cell, pbc)
                        row.append(norm2(displacement))
                    costs.append(row)
                assignment = _hungarian(costs)
                for local_source, local_target in enumerate(assignment):
                    mapping[left[local_source]] = right[local_target]
            current_mapping = tuple(mapping)
            unshifted = []
            for source, target in enumerate(current_mapping):
                raw = vsub(positions_a[source], positions_b[target])
                residual, _ = closest_periodic_displacement(vsub(raw, translation), cell, pbc)
                unshifted.append(vadd(translation, residual))
            refined = (
                sum(vector[0] for vector in unshifted) / count,
                sum(vector[1] for vector in unshifted) / count,
                sum(vector[2] for vector in unshifted) / count,
            )
            refined, _ = closest_periodic_displacement(refined, cell, pbc)
            shift_change = norm(vsub(refined, translation))
            translation = refined
            if current_mapping == previous_mapping and shift_change <= 1e-14:
                break
            previous_mapping = current_mapping
        displacements = [closest_periodic_displacement(
            vsub(positions_a[source], vadd(positions_b[target], translation)), cell, pbc)[0]
            for source, target in enumerate(mapping)]
        squared = sum(norm2(vector) for vector in displacements)
        maximum = max(norm(vector) for vector in displacements)
        candidate = (squared, maximum, tuple(mapping), translation, displacements)
        if best is None or candidate[:4] < best[:4]:
            best = candidate
    if best is None:
        raise DomainFailure("INTERNAL_ERROR", "state matcher produced no candidate",
                            component="geometry", requirement="STATE-006")
    rms = math.sqrt(best[0] / count)
    difference = abs(energy_a - energy_b)
    equal = (rms <= rms_tolerance and best[1] <= max_tolerance
             and difference <= energy_tolerance_per_atom * count)
    return MatchReport(equal, best[2], rms, best[1], difference, best[3])


def max_movable_force(forces: Sequence[Vector], movable: Sequence[bool]) -> float:
    selected = [norm(force) for force, active in zip(forces, movable) if active]
    return max(selected, default=0.0)


def normalized(vector: Vector) -> Vector:
    length = norm(vector)
    if not math.isfinite(length) or length == 0.0:
        raise DomainFailure("INVALID_SADDLE", "zero or nonfinite mode vector",
                            component="saddle", requirement="SADDLE-005")
    return vscale(vector, 1.0 / length)
