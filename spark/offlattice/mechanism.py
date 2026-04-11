"""
Mechanism data class — represents a single min→saddle→min transition.

Ported from openFLY env/mechanisms.hpp.
"""

import numpy as np
import json


class Mechanism:
    """
    A transition mechanism between two local minima via a saddle point.

    Stores the energy barrier, energy change, kinetic pre-factor,
    and the local displacements (initial→saddle and initial→final).

    Parameters
    ----------
    barrier : float
        Activation energy E_sp - E_0 in eV. Must be positive.
    delta : float
        Energy change E_final - E_0 in eV. Can be negative (exothermic).
    kinetic_pre : float
        Arrhenius pre-factor (attempt frequency) in Hz.
    delta_sp : ndarray, shape (n_local, 3)
        Displacement vectors from initial to saddle point for each atom
        in the local environment.
    delta_fwd : ndarray, shape (n_local, 3)
        Displacement vectors from initial to final state for each atom
        in the local environment.
    axis : ndarray, shape (n_local, 3), optional
        Dimer axis at the saddle point (normalized). Used for reconstruction
        hints in subsequent searches.
    err_fwd : float
        Reconstruction error for the forward (final) state in Angstrom.
    err_sp : float
        Reconstruction error for the saddle point in Angstrom.
    poison_fwd : bool
        If True, the forward state cannot be reliably reconstructed.
    poison_sp : bool
        If True, the saddle point cannot be reliably reconstructed.
    """

    __slots__ = (
        'barrier', 'delta', 'kinetic_pre',
        'delta_sp', 'delta_fwd', 'axis',
        'err_fwd', 'err_sp',
        'poison_fwd', 'poison_sp',
    )

    def __init__(self, barrier, delta, kinetic_pre,
                 delta_sp, delta_fwd, axis=None,
                 err_fwd=0.0, err_sp=0.0,
                 poison_fwd=False, poison_sp=False):
        self.barrier = float(barrier)
        self.delta = float(delta)
        self.kinetic_pre = float(kinetic_pre)
        self.delta_sp = np.asarray(delta_sp, dtype=np.float64)
        self.delta_fwd = np.asarray(delta_fwd, dtype=np.float64)
        self.axis = np.asarray(axis, dtype=np.float64) if axis is not None else None
        self.err_fwd = float(err_fwd)
        self.err_sp = float(err_sp)
        self.poison_fwd = bool(poison_fwd)
        self.poison_sp = bool(poison_sp)

    @property
    def n_atoms(self):
        """Number of atoms in the local environment."""
        return self.delta_fwd.shape[0]

    @property
    def reverse_barrier(self):
        """Barrier for the reverse transition: E_sp - E_final."""
        return self.barrier - self.delta

    def rate(self, temperature):
        """
        Compute Arrhenius rate at given temperature.

        k = kinetic_pre * exp(-barrier / (k_B * T))

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin.

        Returns
        -------
        float
            Rate in s^-1.
        """
        INV_BOLTZ = 16021766340.0 / 1380649.0  # 1/eV, matches openFLY
        return self.kinetic_pre * np.exp(-self.barrier * INV_BOLTZ / temperature)

    def to_dict(self):
        """Serialize to dictionary."""
        d = {
            'barrier': self.barrier,
            'delta': self.delta,
            'kinetic_pre': self.kinetic_pre,
            'delta_sp': self.delta_sp.tolist(),
            'delta_fwd': self.delta_fwd.tolist(),
            'err_fwd': self.err_fwd,
            'err_sp': self.err_sp,
            'poison_fwd': self.poison_fwd,
            'poison_sp': self.poison_sp,
        }
        if self.axis is not None:
            d['axis'] = self.axis.tolist()
        return d

    @classmethod
    def from_dict(cls, d):
        """Deserialize from dictionary."""
        return cls(
            barrier=d['barrier'],
            delta=d['delta'],
            kinetic_pre=d['kinetic_pre'],
            delta_sp=np.array(d['delta_sp']),
            delta_fwd=np.array(d['delta_fwd']),
            axis=np.array(d['axis']) if 'axis' in d else None,
            err_fwd=d.get('err_fwd', 0.0),
            err_sp=d.get('err_sp', 0.0),
            poison_fwd=d.get('poison_fwd', False),
            poison_sp=d.get('poison_sp', False),
        )

    def __repr__(self):
        return (f"Mechanism(barrier={self.barrier:.4f} eV, "
                f"delta={self.delta:.4f} eV, "
                f"pre={self.kinetic_pre:.2e} Hz, "
                f"n_atoms={self.n_atoms})")
