"""
Calculator adapter — wraps ASE calculators for the off-lattice KMC engine.

Provides a unified interface for energy, gradient, and Hessian calculations.
Supports any ASE-compatible calculator (VASP, ORCA, EMT, MACE, etc.).
"""

import numpy as np
from copy import deepcopy


class CalculatorAdapter:
    """
    Wraps an ASE calculator to provide energy, gradient, and Hessian.

    Parameters
    ----------
    calculator : ase.calculators.calculator.Calculator
        Any ASE-compatible calculator instance.
    frozen_mask : ndarray of bool, shape (n_atoms,), optional
        If provided, atoms where frozen_mask[i] is True have zero gradient.
    """

    def __init__(self, calculator, frozen_mask=None):
        self.calculator = calculator
        self.frozen_mask = frozen_mask

    def energy(self, atoms):
        """
        Compute potential energy.

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic configuration.

        Returns
        -------
        float
            Potential energy in eV.
        """
        atoms = atoms.copy()
        atoms.calc = self.calculator
        return atoms.get_potential_energy()

    def gradient(self, atoms):
        """
        Compute potential energy gradient (negative forces).

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic configuration.

        Returns
        -------
        energy : float
            Potential energy in eV.
        grad : ndarray, shape (n_atoms, 3)
            Energy gradient in eV/Angstrom. Frozen atoms have zero gradient.
        """
        atoms = atoms.copy()
        atoms.calc = self.calculator
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        grad = -forces

        if self.frozen_mask is not None:
            grad[self.frozen_mask] = 0.0

        return energy, grad

    def hessian(self, atoms, dx=0.01):
        """
        Compute Hessian matrix via finite differences of forces.

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic configuration at a stationary point.
        dx : float
            Finite difference step size in Angstrom.

        Returns
        -------
        H : ndarray, shape (3*n_atoms, 3*n_atoms)
            Hessian matrix in eV/Angstrom^2.
        """
        n = len(atoms)
        H = np.zeros((3 * n, 3 * n))
        pos0 = atoms.get_positions().copy()

        for i in range(3 * n):
            atom_idx = i // 3
            cart_idx = i % 3

            if self.frozen_mask is not None and self.frozen_mask[atom_idx]:
                continue

            # Forward displacement
            atoms_plus = atoms.copy()
            atoms_plus.calc = self.calculator
            pos_plus = pos0.copy()
            pos_plus[atom_idx, cart_idx] += dx
            atoms_plus.set_positions(pos_plus)
            f_plus = atoms_plus.get_forces().ravel()

            # Backward displacement
            atoms_minus = atoms.copy()
            atoms_minus.calc = self.calculator
            pos_minus = pos0.copy()
            pos_minus[atom_idx, cart_idx] -= dx
            atoms_minus.set_positions(pos_minus)
            f_minus = atoms_minus.get_forces().ravel()

            # Central difference: H_ij = -df_j / dx_i
            H[i, :] = -(f_plus - f_minus) / (2.0 * dx)

        # Symmetrize
        H = 0.5 * (H + H.T)

        # Zero out frozen DOFs
        if self.frozen_mask is not None:
            for i in range(n):
                if self.frozen_mask[i]:
                    for c in range(3):
                        idx = 3 * i + c
                        H[idx, :] = 0.0
                        H[:, idx] = 0.0

        return H

    def mass_weighted_hessian(self, atoms, dx=0.01):
        """
        Compute mass-weighted Hessian: H_mw[i,j] = H[i,j] / sqrt(m_i * m_j).

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic configuration.
        dx : float
            Finite difference step size.

        Returns
        -------
        H_mw : ndarray, shape (3*n_atoms, 3*n_atoms)
            Mass-weighted Hessian.
        """
        H = self.hessian(atoms, dx)
        masses = atoms.get_masses()
        n = len(atoms)
        inv_sqrt_m = np.zeros(3 * n)
        for i in range(n):
            m = masses[i]
            if m > 0:
                inv_sqrt_m[3*i:3*i+3] = 1.0 / np.sqrt(m)
        # H_mw[i,j] = H[i,j] * inv_sqrt_m[i] * inv_sqrt_m[j]
        H_mw = H * np.outer(inv_sqrt_m, inv_sqrt_m)
        return H_mw

    @property
    def r_cut(self):
        """Cutoff radius of the potential, if available."""
        if hasattr(self.calculator, 'r_cut'):
            return self.calculator.r_cut
        if hasattr(self.calculator, 'rc'):
            return self.calculator.rc
        # Fallback for many ML potentials
        if hasattr(self.calculator, 'cutoff'):
            return self.calculator.cutoff
        return 5.0  # Default fallback
