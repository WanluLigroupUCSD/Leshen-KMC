"""
Geometry minimizer wrapping scipy L-BFGS-B.

Ported from openFLY minimise/LBFGS.
"""

import numpy as np
from scipy.optimize import minimize as scipy_minimize


class Minimizer:
    """
    Local geometry minimizer using L-BFGS-B via scipy.

    Parameters
    ----------
    calc_adapter : CalculatorAdapter
        Provides energy and gradient for atomic configurations.
    f_tol : float
        Force convergence criterion: max|f| < f_tol in eV/Angstrom.
    max_steps : int
        Maximum number of L-BFGS-B iterations.
    max_step : float
        Maximum step size in Angstrom.
    debug : bool
        Print convergence info.
    """

    def __init__(self, calc_adapter, f_tol=0.01, max_steps=500,
                 max_step=0.2, debug=False):
        self.calc = calc_adapter
        self.f_tol = f_tol
        self.max_steps = max_steps
        self.max_step = max_step
        self.debug = debug
        self._n_calls = 0

    def minimize(self, atoms):
        """
        Relax atomic positions to a local energy minimum.

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic configuration (modified in-place).

        Returns
        -------
        atoms : ase.Atoms
            Relaxed configuration (same object, modified).
        energy : float
            Final potential energy in eV.
        converged : bool
            Whether force criterion was met.
        n_steps : int
            Number of iterations taken.
        """
        frozen = getattr(self.calc, 'frozen_mask', None)
        pos0 = atoms.get_positions().copy()
        n = len(atoms)

        # Determine free DOFs
        if frozen is not None:
            free_mask = ~frozen
        else:
            free_mask = np.ones(n, dtype=bool)

        free_indices = np.where(free_mask)[0]
        n_free = len(free_indices)

        if n_free == 0:
            e, _ = self.calc.gradient(atoms)
            return atoms, e, True, 0

        def _objective(x_flat):
            """Objective for scipy: x_flat contains free atom positions."""
            pos = pos0.copy()
            pos[free_indices] = x_flat.reshape(-1, 3)
            atoms.set_positions(pos)
            e, grad = self.calc.gradient(atoms)
            grad_free = grad[free_indices].ravel()
            self._n_calls += 1
            return e, grad_free

        x0 = pos0[free_indices].ravel()

        # Run L-BFGS-B
        result = scipy_minimize(
            _objective,
            x0,
            method='L-BFGS-B',
            jac=True,
            options={
                'maxiter': self.max_steps,
                'ftol': 1e-15,  # use force criterion, not energy
                'gtol': self.f_tol,
                'maxcor': 10,
                'maxls': 20,
            },
        )

        # Apply final positions
        pos_final = pos0.copy()
        pos_final[free_indices] = result.x.reshape(-1, 3)
        atoms.set_positions(pos_final)

        # Check convergence via force norm
        e_final, grad_final = self.calc.gradient(atoms)
        max_force = np.max(np.linalg.norm(
            -grad_final[free_indices], axis=1))
        converged = max_force < self.f_tol

        if self.debug:
            status = "converged" if converged else "NOT converged"
            print(f"MIN: {status} in {result.nit} steps, "
                  f"E={e_final:.6f} eV, max|f|={max_force:.4f} eV/A")

        return atoms, e_final, converged, result.nit

    def minimize_from_positions(self, atoms, positions):
        """
        Set positions and minimize.

        Parameters
        ----------
        atoms : ase.Atoms
            Template atoms (with calculator, cell, etc.).
        positions : ndarray, shape (n_atoms, 3)
            Starting positions.

        Returns
        -------
        Same as minimize().
        """
        atoms_copy = atoms.copy()
        atoms_copy.calc = self.calc.calculator
        atoms_copy.set_positions(positions)
        return self.minimize(atoms_copy)
