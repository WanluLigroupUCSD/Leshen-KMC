"""
Local environment detection, fingerprinting, and matching.

Ported from openFLY env/geometry.hpp + env/heuristics.hpp.

Three-step environment matching pipeline:
  1. Graph canonical hash (fast topology check)
  2. Fingerprint equivalence (sorted interatomic distances)
  3. Full geometry matching (Kabsch + greedy permutation)
"""

import numpy as np
from hashlib import sha256
from itertools import combinations


# ---------------------------------------------------------------------------
# Utility: Kabsch alignment
# ---------------------------------------------------------------------------

def kabsch(X, Y):
    """
    Kabsch algorithm: find optimal rotation O minimizing ||O @ X - Y||.

    Parameters
    ----------
    X, Y : ndarray, shape (n, 3)
        Centered point sets (centroid at origin).

    Returns
    -------
    O : ndarray, shape (3, 3)
        Optimal orthogonal matrix.
    rmsd : float
        RMSD after alignment.
    """
    H = X.T @ Y  # covariance matrix
    U, S, Vt = np.linalg.svd(H)
    # Correct for reflection
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1.0, 1.0, np.sign(d)])
    O = Vt.T @ sign_matrix @ U.T
    # RMSD
    diff = O @ X.T - Y.T  # (3, n)
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=0)))
    return O, rmsd


def grmsd(O, X, Y):
    """Generalized RMSD: sqrt(mean(||O @ x_i - y_i||^2))."""
    diff = (O @ X.T - Y.T)  # (3, n)
    return np.sqrt(np.mean(np.sum(diff**2, axis=0)))


# ---------------------------------------------------------------------------
# Fingerprint
# ---------------------------------------------------------------------------

class Fingerprint:
    """
    Fast fingerprint for local environment equivalence testing.

    Stores sorted center-to-neighbor distances (r_0j) and sorted
    inter-neighbor distances (r_ij). Two environments with different
    fingerprints cannot be equivalent.

    Parameters
    ----------
    positions : ndarray, shape (n, 3)
        Positions relative to the central atom (center at origin).
    """

    __slots__ = ('r_0j', 'r_ij', '_r_min')

    def __init__(self, positions):
        n = len(positions)
        # Distances from center (index 0) to neighbors
        if n > 1:
            dr = positions[1:]  # center is at origin
            self.r_0j = np.sort(np.linalg.norm(dr, axis=1))
        else:
            self.r_0j = np.array([], dtype=np.float64)

        # Inter-neighbor distances (sorted)
        if n > 2:
            dists = []
            for i, j in combinations(range(1, n), 2):
                dists.append(np.linalg.norm(positions[i] - positions[j]))
            self.r_ij = np.sort(dists)
        else:
            self.r_ij = np.array([], dtype=np.float64)

        self._r_min = self.r_0j[0] if len(self.r_0j) > 0 else 1.0

    def r_min(self):
        """Minimum center-neighbor distance."""
        return self._r_min

    def equiv(self, other, delta):
        """
        Fast equivalence check against another fingerprint.

        Returns True if all sorted distances match within delta * sqrt(2).
        This is a necessary but not sufficient condition for geometric
        equivalence — false positives are resolved by full permutation.

        Parameters
        ----------
        other : Fingerprint
        delta : float
            Tolerance (L2 norm).

        Returns
        -------
        bool
        """
        tol = delta * np.sqrt(2.0)

        if len(self.r_0j) != len(other.r_0j):
            return False
        if np.any(np.abs(self.r_0j - other.r_0j) > tol):
            return False
        if len(self.r_ij) != len(other.r_ij):
            return False
        if len(self.r_ij) > 0 and np.any(np.abs(self.r_ij - other.r_ij) > tol):
            return False
        return True


# ---------------------------------------------------------------------------
# Graph hash (replaces nauty canonicalization)
# ---------------------------------------------------------------------------

def _graph_hash(positions, colours, r_edge):
    """
    Compute a canonical hash of the local environment graph.

    Builds a connectivity graph (edge if dist < r_edge), encodes with
    sorted adjacency + colour information, and hashes with SHA-256.

    This replaces the nauty-based canonicalization in openFLY with a
    simpler but effective approach using sorted invariants.

    Parameters
    ----------
    positions : ndarray, shape (n, 3)
        Atom positions (center at origin).
    colours : ndarray of int, shape (n,)
        Atom type encoding: 2 * type_id + is_frozen.
    r_edge : float
        Bonding distance for graph edges.

    Returns
    -------
    int
        Hash value.
    """
    n = len(positions)
    if n == 0:
        return 0

    # Build adjacency per atom: sorted list of (colour, dist) to neighbors
    adj_invariants = []
    for i in range(n):
        neighbors = []
        for j in range(n):
            if i == j:
                continue
            d = np.linalg.norm(positions[i] - positions[j])
            if d < r_edge:
                neighbors.append(colours[j])
        neighbors.sort()
        # Per-atom invariant: (own colour, degree, sorted neighbor colours)
        adj_invariants.append((int(colours[i]), len(neighbors), tuple(neighbors)))

    # Sort to get canonical ordering
    adj_invariants.sort()

    # Centre atom colour + colour counts
    colour_counts = np.bincount(colours.astype(int), minlength=1)

    # Combine into hashable representation
    hash_input = (
        int(colours[0]),  # center colour
        tuple(colour_counts.tolist()),
        tuple(adj_invariants),
    )

    h = sha256(repr(hash_input).encode()).hexdigest()
    return int(h[:16], 16)  # 64-bit hash


# ---------------------------------------------------------------------------
# Geometry: local environment representation
# ---------------------------------------------------------------------------

class Geometry:
    """
    Local environment around a central atom.

    The center atom is always at index 0 with position at the origin.
    Neighbour atoms are stored with their displacement vectors from center.

    Parameters
    ----------
    positions : ndarray, shape (n, 3)
        Positions of atoms in local frame (center at origin).
    colours : ndarray of int, shape (n,)
        Colour encoding: 2 * type_id + is_frozen.
    indices : ndarray of int, shape (n,)
        Global atom indices in the parent system.
    """

    def __init__(self, positions, colours, indices):
        self.positions = np.asarray(positions, dtype=np.float64)
        self.colours = np.asarray(colours, dtype=int)
        self.indices = np.asarray(indices, dtype=int)
        self._fingerprint = None
        self._hash = None

    @property
    def n_atoms(self):
        return len(self.positions)

    @property
    def center_index(self):
        """Global index of the center atom."""
        return int(self.indices[0])

    def center(self):
        """Shift so centroid is at origin."""
        c = np.mean(self.positions, axis=0)
        self.positions -= c
        self._fingerprint = None
        self._hash = None

    def fingerprint(self):
        """Get or compute the Fingerprint for this geometry."""
        if self._fingerprint is None:
            self._fingerprint = Fingerprint(self.positions)
        return self._fingerprint

    def graph_hash(self, r_edge):
        """Compute the graph canonical hash."""
        if self._hash is None:
            self._hash = _graph_hash(self.positions, self.colours, r_edge)
        return self._hash

    def permute_onto(self, reference, delta):
        """
        Find a permutation + rotation mapping this geometry onto reference.

        Uses greedy recursive matching with Kabsch alignment.
        Returns (rotation_matrix, rmsd, permutation) or None if no match.

        Parameters
        ----------
        reference : Geometry
            Reference environment (from catalogue).
        delta : float
            Maximum allowed RMSD for a match.

        Returns
        -------
        result : dict or None
            If matched: {'O': rotation_matrix, 'rmsd': rmsd, 'perm': permutation}
            None if no matching permutation found within tolerance.
        """
        n = self.n_atoms
        if n != reference.n_atoms:
            return None

        # Build distance-based compatibility
        # For greedy matching: for each atom in self, which atoms in reference
        # could it map to?
        best_result = [None]
        best_rmsd = [delta]  # only accept matches below delta

        def _greedy_match(assigned_self, assigned_ref, perm, depth):
            if depth == n:
                # All atoms assigned — compute final alignment
                X = self.positions[perm]
                Y = reference.positions
                O, rmsd_val = kabsch(X, Y)
                if rmsd_val < best_rmsd[0]:
                    best_rmsd[0] = rmsd_val
                    best_result[0] = {
                        'O': O.copy(),
                        'rmsd': rmsd_val,
                        'perm': list(perm),
                    }
                return

            # Index in reference we're trying to match
            ref_idx = depth

            # Try each unassigned atom in self
            ref_col = reference.colours[ref_idx]
            for self_idx in range(n):
                if self_idx in assigned_self:
                    continue
                # Must match colour
                if self.colours[self_idx] != ref_col:
                    continue

                # Distance pruning: check already-assigned pairs
                ok = True
                for k in range(depth):
                    d_self = np.linalg.norm(
                        self.positions[perm[k]] - self.positions[self_idx])
                    d_ref = np.linalg.norm(
                        reference.positions[k] - reference.positions[ref_idx])
                    if abs(d_self - d_ref) > delta * np.sqrt(2.0):
                        ok = False
                        break
                if not ok:
                    continue

                perm.append(self_idx)
                assigned_self.add(self_idx)
                _greedy_match(assigned_self, assigned_ref, perm, depth + 1)
                assigned_self.discard(self_idx)
                perm.pop()

                # Early exit if we found a match below threshold
                if best_result[0] is not None:
                    return

        _greedy_match(set(), set(), [], 0)
        return best_result[0]

    def self_symmetries(self, delta):
        """
        Find all symmetries (rotation + permutation) of this geometry.

        Returns list of (O, perm) tuples where O is a rotation matrix
        and perm is the atom permutation that maps the geometry onto itself.

        Parameters
        ----------
        delta : float
            Tolerance for self-equivalence.

        Returns
        -------
        syms : list of (ndarray, list)
            Each element is (rotation_matrix_3x3, permutation_list).
        """
        n = self.n_atoms
        syms = []

        def _find_syms(assigned_self, perm, depth):
            if depth == n:
                X = self.positions[perm]
                Y = self.positions
                O, rmsd_val = kabsch(X, Y)
                if rmsd_val < delta:
                    syms.append((O.copy(), list(perm)))
                return

            ref_idx = depth
            ref_col = self.colours[ref_idx]

            for self_idx in range(n):
                if self_idx in assigned_self:
                    continue
                if self.colours[self_idx] != ref_col:
                    continue

                ok = True
                for k in range(depth):
                    d_self = np.linalg.norm(
                        self.positions[perm[k]] - self.positions[self_idx])
                    d_ref = np.linalg.norm(
                        self.positions[k] - self.positions[ref_idx])
                    if abs(d_self - d_ref) > delta * np.sqrt(2.0):
                        ok = False
                        break
                if not ok:
                    continue

                perm.append(self_idx)
                assigned_self.add(self_idx)
                _find_syms(assigned_self, perm, depth + 1)
                assigned_self.discard(self_idx)
                perm.pop()

        _find_syms(set(), [], 0)
        return syms


# ---------------------------------------------------------------------------
# Build local environment from atoms + neighbor list
# ---------------------------------------------------------------------------

def build_geometry(center_idx, positions, types, frozen, neighbor_indices,
                   r_env):
    """
    Build a Geometry for the environment around a central atom.

    Parameters
    ----------
    center_idx : int
        Index of the central atom.
    positions : ndarray, shape (n_atoms, 3)
        All atom positions.
    types : ndarray of int, shape (n_atoms,)
        Atom type IDs (species index).
    frozen : ndarray of bool, shape (n_atoms,)
        Frozen atom flags.
    neighbor_indices : list of int
        Indices of atoms neighboring the center within r_env.
    r_env : float
        Environment radius cutoff.

    Returns
    -------
    Geometry
    """
    center_pos = positions[center_idx]

    local_pos = [np.zeros(3)]  # center at origin
    local_col = [2 * int(types[center_idx]) + int(frozen[center_idx])]
    local_idx = [center_idx]

    for j in neighbor_indices:
        if j == center_idx:
            continue
        dr = positions[j] - center_pos
        dist = np.linalg.norm(dr)
        if dist < r_env:
            local_pos.append(dr)
            local_col.append(2 * int(types[j]) + int(frozen[j]))
            local_idx.append(j)

    geo = Geometry(
        np.array(local_pos),
        np.array(local_col),
        np.array(local_idx),
    )
    geo.center()
    return geo


def build_neighbor_list(positions, cell, pbc, r_cut):
    """
    Build neighbor list using ASE or grid-based approach.

    Parameters
    ----------
    positions : ndarray, shape (n_atoms, 3)
    cell : ndarray, shape (3, 3)
        Unit cell matrix.
    pbc : array-like of bool, shape (3,)
        Periodic boundary conditions.
    r_cut : float
        Cutoff radius.

    Returns
    -------
    neighbors : list of list of int
        neighbors[i] = list of indices of atoms within r_cut of atom i.
    """
    try:
        from ase.neighborlist import neighbor_list
        from ase import Atoms

        atoms = Atoms(positions=positions, cell=cell, pbc=pbc)
        i_list, j_list = neighbor_list('ij', atoms, r_cut)

        n = len(positions)
        neighbors = [[] for _ in range(n)]
        for i, j in zip(i_list, j_list):
            neighbors[i].append(j)
        return neighbors

    except ImportError:
        # Fallback: brute force with minimum image convention
        n = len(positions)
        neighbors = [[] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                dr = positions[j] - positions[i]
                # Minimum image (orthorhombic only)
                if any(pbc):
                    cell_diag = np.diag(cell) if cell.ndim == 2 else cell
                    for dim in range(3):
                        if pbc[dim] and cell_diag[dim] > 0:
                            dr[dim] -= round(dr[dim] / cell_diag[dim]) * cell_diag[dim]
                dist = np.linalg.norm(dr)
                if dist < r_cut:
                    neighbors[i].append(j)
        return neighbors
