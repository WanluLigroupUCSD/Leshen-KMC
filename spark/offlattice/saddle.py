"""
Saddle point finding — dimer method + master coordinator.

Ported from openFLY saddle/dimer.cpp + saddle/find.cpp.

The dimer method locates first-order saddle points by:
  1. Random perturbation from a minimum
  2. Rotation to align with lowest curvature mode
  3. Translation along effective gradient toward saddle point
  4. Collision detection against known saddle points
"""

import numpy as np
from .mechanism import Mechanism


class DimerSearch:
    """
    Dimer method for saddle point finding.

    Parameters
    ----------
    calc_adapter : CalculatorAdapter
        Provides energy and gradient.
    dimer_sep : float
        Dimer half-separation in Angstrom.
    f_tol : float
        Force convergence criterion for saddle point.
    max_steps : int
        Maximum number of dimer translation steps.
    max_rotor_steps : int
        Maximum rotor (rotation) iterations per translation step.
    rotor_tol : float
        Convergence tolerance for dimer rotation (torque norm).
    trust : float
        Initial trust radius for translation step.
    trust_grow : float
        Factor to grow trust radius on good steps.
    trust_shrink : float
        Factor to shrink trust radius on bad steps.
    convex_max : int
        Maximum consecutive convex curvature steps before giving up.
    debug : bool
        Print convergence info.
    """

    def __init__(self, calc_adapter, dimer_sep=0.005, f_tol=0.05,
                 max_steps=300, max_rotor_steps=10, rotor_tol=0.01,
                 trust=0.1, trust_grow=1.1, trust_shrink=0.5,
                 convex_max=5, debug=False):
        self.calc = calc_adapter
        self.dimer_sep = dimer_sep
        self.f_tol = f_tol
        self.max_steps = max_steps
        self.max_rotor_steps = max_rotor_steps
        self.rotor_tol = rotor_tol
        self.trust = trust
        self.trust_grow = trust_grow
        self.trust_shrink = trust_shrink
        self.convex_max = convex_max
        self.debug = debug

    def find_saddle(self, atoms, axis=None, history_sp=None, theta_tol=None):
        """
        Search for a saddle point starting from atoms configuration.

        Parameters
        ----------
        atoms : ase.Atoms
            Starting configuration (perturbed from minimum).
        axis : ndarray, shape (n_atoms, 3), optional
            Initial dimer axis. If None, uses random direction.
        history_sp : list of ndarray, optional
            Previous saddle point positions for collision detection.
        theta_tol : float, optional
            Angular tolerance for collision detection (radians).

        Returns
        -------
        result : dict or None
            If successful:
                'positions': saddle point positions
                'axis': dimer axis at saddle point
                'energy': energy at saddle point
                'status': 'success', 'convex', 'collision', 'maxiter'
            None if search failed entirely.
        """
        frozen = getattr(self.calc, 'frozen_mask', None)
        n = len(atoms)
        pos = atoms.get_positions().copy()

        # Initialize dimer axis
        if axis is None:
            axis = np.random.randn(n, 3)
            if frozen is not None:
                axis[frozen] = 0.0
            norm = np.linalg.norm(axis)
            if norm > 0:
                axis /= norm

        trust = self.trust
        convex_count = 0
        status = 'maxiter'

        for step in range(self.max_steps):
            # --- Rotation phase: align axis with lowest curvature mode ---
            axis, curvature = self._rotate_dimer(atoms, pos, axis, frozen)

            # --- Compute effective gradient ---
            e_center, g_center = self.calc.gradient(atoms)
            g_eff = self._effective_gradient(g_center, axis, curvature)

            if frozen is not None:
                g_eff[frozen] = 0.0

            # --- Check convergence ---
            g_norm = np.linalg.norm(g_eff)
            if g_norm < self.f_tol:
                status = 'success'
                if self.debug:
                    print(f"DIMER: Converged at step {step}, "
                          f"E={e_center:.6f}, |g|={g_norm:.6f}, "
                          f"curv={curvature:.4f}")
                break

            # --- Convex check ---
            if curvature > 0:
                convex_count += 1
                if convex_count >= self.convex_max:
                    status = 'convex'
                    if self.debug:
                        print(f"DIMER: Convex exit at step {step}")
                    break
            else:
                convex_count = 0

            # --- Collision detection ---
            if history_sp and theta_tol is not None:
                for sp_prev in history_sp:
                    disp = pos - sp_prev
                    d = np.linalg.norm(disp)
                    if d > 0:
                        cos_theta = np.sum(disp * axis.ravel()) / (
                            d * np.linalg.norm(axis))
                        if abs(cos_theta) > np.cos(theta_tol):
                            status = 'collision'
                            if self.debug:
                                print(f"DIMER: Collision at step {step}")
                            return {'positions': pos, 'axis': axis,
                                    'energy': e_center, 'status': status}

            # --- Translation step ---
            step_vec = -g_eff  # Descend effective gradient
            step_norm = np.linalg.norm(step_vec)
            if step_norm > 0:
                step_vec *= min(1.0, trust / step_norm)

            pos += step_vec
            atoms.set_positions(pos)

            # --- Trust radius adaptation ---
            if step > 0:
                proj = np.sum(step_vec.ravel() * g_eff.ravel())
                if proj < -0.5 * trust:
                    trust *= self.trust_shrink
                elif proj > 0.5 * trust:
                    trust = min(trust * self.trust_grow, 0.5)

        e_final, _ = self.calc.gradient(atoms)
        return {
            'positions': pos.copy(),
            'axis': axis.copy(),
            'energy': e_final,
            'status': status,
        }

    def _rotate_dimer(self, atoms, pos, axis, frozen):
        """
        Rotate dimer to align with lowest curvature mode.

        Returns updated axis and curvature.
        """
        n = len(atoms)
        h = self.dimer_sep

        for rot_step in range(self.max_rotor_steps):
            # Evaluate at pos + h*axis and pos - h*axis
            atoms_plus = atoms.copy()
            atoms_plus.calc = self.calc.calculator
            atoms_plus.set_positions(pos + h * axis)
            _, g_plus = self.calc.gradient(atoms_plus)

            atoms_minus = atoms.copy()
            atoms_minus.calc = self.calc.calculator
            atoms_minus.set_positions(pos - h * axis)
            _, g_minus = self.calc.gradient(atoms_minus)

            # Curvature estimate
            dg = g_plus - g_minus
            curvature = np.sum(dg * axis) / (2.0 * h)

            # Rotor "torque" — we want to MINIMIZE curvature
            #     κ(N̂) = dg · N̂ / (2h),
            # whose gradient in N̂-space is dg/(2h). To DECREASE κ we step in
            # direction -dg/(2h), projected perpendicular to N̂ to keep |N̂|=1.
            #
            # Sign fix 2026-05-01 (a3d6b14 in Rust, mirrored here): the
            # original openFLY-derived port had a missing minus sign. The
            # rotor was rotating axis toward HIGHER curvature, so most random
            # initial axes hit convex exits before finding a saddle. SPARK
            # production OTF runs on Cu slab benchmark show 1/20 success rate
            # without this fix; 20/20 with it. See spark-rs/src/offlattice/
            # saddle.rs for the Rust counterpart and Phase D bench results.
            dg_per_2h = dg / (2.0 * h)
            proj = np.sum(dg_per_2h * axis)
            torque = -(dg_per_2h - proj * axis)

            if frozen is not None:
                torque[frozen] = 0.0

            torque_norm = np.linalg.norm(torque)
            if torque_norm < self.rotor_tol:
                break

            # Rotate axis toward torque direction
            # Use finite-difference approach
            theta = torque / torque_norm
            # Small rotation angle from torque magnitude
            angle = min(0.1, torque_norm * 0.1)  # damped rotation

            # Rodrigues rotation: rotate axis toward theta
            new_axis = axis * np.cos(angle) + theta * np.sin(angle)

            if frozen is not None:
                new_axis[frozen] = 0.0

            norm = np.linalg.norm(new_axis)
            if norm > 0:
                axis = new_axis / norm

        return axis, curvature

    def _effective_gradient(self, gradient, axis, curvature):
        """
        Compute effective gradient for translation.

        If curvature < 0 (saddle region):
            g_eff = g - 2*(g.axis)*axis  (invert parallel component)
        If curvature > 0 (convex):
            g_eff = -(g.axis)*axis  (only climb along axis)
        """
        g_par = np.sum(gradient * axis) * axis

        if curvature < 0:
            return gradient - 2.0 * g_par
        else:
            return -g_par


class SaddleMaster:
    """
    Master coordinator for saddle point searches.

    Manages parallel searches, mechanism extraction, symmetry exploitation,
    and duplicate detection.

    Parameters
    ----------
    dimer : DimerSearch
        The dimer method instance.
    minimizer : Minimizer
        For relaxing trial structures after dimer finds saddle.
    calc_adapter : CalculatorAdapter
        For energy/gradient evaluations.
    r_pert : float
        Perturbation envelope radius in Angstrom.
    stddev : float
        Perturbation standard deviation.
    max_searches : int
        Maximum number of saddle point searches per environment.
    max_failed : int
        Maximum consecutive failed searches before stopping.
    batch_size : int
        Number of parallel searches per batch.
    mech_tol : float
        RMSD tolerance for distinguishing mechanisms.
    basin_tol : float
        Distance tolerance for same-basin detection.
    nudge_frac : float
        Fraction of displacement for nudging from SP to minima.
    capture_r_tol : float
        Position tolerance for mechanism reconstruction validation.
    capture_e_tol : float
        Energy tolerance for mechanism reconstruction validation.
    debug : bool
    """

    def __init__(self, dimer, minimizer, calc_adapter,
                 r_pert=3.0, stddev=0.3, max_searches=50,
                 max_failed=20, batch_size=4, mech_tol=0.2,
                 basin_tol=0.3, nudge_frac=0.1,
                 capture_r_tol=0.5, capture_e_tol=0.1,
                 debug=False):
        self.dimer = dimer
        self.minimizer = minimizer
        self.calc = calc_adapter
        self.r_pert = r_pert
        self.stddev = stddev
        self.max_searches = max_searches
        self.max_failed = max_failed
        self.batch_size = batch_size
        self.mech_tol = mech_tol
        self.basin_tol = basin_tol
        self.nudge_frac = nudge_frac
        self.capture_r_tol = capture_r_tol
        self.capture_e_tol = capture_e_tol
        self.debug = debug

    def find_mechanisms(self, atoms, center_idx, geometry, catalogue=None):
        """
        Find all transition mechanisms from the environment around center_idx.

        Parameters
        ----------
        atoms : ase.Atoms
            Current minimized configuration.
        center_idx : int
            Index of the central atom.
        geometry : Geometry
            Local environment geometry.
        catalogue : Catalogue, optional
            If provided, used for reconstruction validation.

        Returns
        -------
        mechanisms : list of Mechanism
            Discovered mechanisms, sorted by barrier height.
        """
        mechanisms = []
        history_sp = []
        n_failed = 0
        n_total = 0

        pos_min = atoms.get_positions().copy()
        e_min = self.calc.energy(atoms)

        while n_total < self.max_searches and n_failed < self.max_failed:
            # Adaptive collision tolerance
            theta_tol = ((30 - 2.5) * np.exp(-0.02 * n_total) + 2.5) \
                        / 360.0 * 2.0 * np.pi

            # Generate perturbation
            perturbed_atoms = atoms.copy()
            perturbed_atoms.calc = self.calc.calculator
            perturbed_pos = self._perturb(pos_min, center_idx, geometry)
            perturbed_atoms.set_positions(perturbed_pos)

            # Random initial axis
            axis = self._random_axis(len(atoms), center_idx, geometry)

            # Run dimer search
            result = self.dimer.find_saddle(
                perturbed_atoms,
                axis=axis,
                history_sp=history_sp,
                theta_tol=theta_tol,
            )

            n_total += 1

            if result is None or result['status'] not in ('success',):
                n_failed += 1
                continue

            # Found a saddle point — extract mechanism
            sp_pos = result['positions']
            sp_axis = result['axis']
            e_sp = result['energy']

            # Check it's actually a saddle (barrier > 0)
            barrier = e_sp - e_min
            if barrier <= 0.005:
                n_failed += 1
                continue

            # Nudge from SP to find forward and reverse minima
            fwd_mech = self._extract_mechanism(
                atoms, pos_min, sp_pos, sp_axis, e_min, e_sp, geometry)

            if fwd_mech is None:
                n_failed += 1
                continue

            # Check if this is a new mechanism
            if self._is_new(fwd_mech, mechanisms):
                mechanisms.append(fwd_mech)
                history_sp.append(sp_pos.copy())
                n_failed = 0  # Reset failure counter

                if self.debug:
                    print(f"MASTER: Found mechanism #{len(mechanisms)}: "
                          f"barrier={fwd_mech.barrier:.4f} eV, "
                          f"delta={fwd_mech.delta:.4f} eV")

                # Generate symmetric copies
                syms = geometry.self_symmetries(
                    geometry.fingerprint().r_min() * 0.4)
                for O, perm in syms[1:]:  # Skip identity
                    sym_mech = self._apply_symmetry(fwd_mech, O, perm)
                    if self._is_new(sym_mech, mechanisms):
                        mechanisms.append(sym_mech)
            else:
                n_failed += 1

        # Sort by barrier
        mechanisms.sort(key=lambda m: m.barrier)

        if self.debug:
            print(f"MASTER: Found {len(mechanisms)} mechanisms in "
                  f"{n_total} searches ({n_failed} failed)")

        return mechanisms

    def _perturb(self, positions, center_idx, geometry):
        """
        Generate a perturbation of the minimum configuration.

        Gaussian displacement weighted by distance from center atom,
        with exponential envelope decay.
        """
        pos = positions.copy()
        n = len(pos)
        center_pos = pos[center_idx]

        # Perturbation for center atom
        stddev = np.random.normal(self.stddev, self.stddev / 3.0)
        stddev = max(0.01, abs(stddev))

        for j in range(geometry.n_atoms):
            global_idx = int(geometry.indices[j])
            r = np.linalg.norm(pos[global_idx] - center_pos)
            weight = max(0.0, 1.0 - r / self.r_pert)
            displacement = np.random.randn(3) * stddev * weight
            pos[global_idx] += displacement

        # Respect frozen atoms
        frozen = getattr(self.calc, 'frozen_mask', None)
        if frozen is not None:
            pos[frozen] = positions[frozen]

        return pos

    def _random_axis(self, n_atoms, center_idx, geometry):
        """Generate random initial dimer axis."""
        axis = np.zeros((n_atoms, 3))
        for j in range(geometry.n_atoms):
            global_idx = int(geometry.indices[j])
            axis[global_idx] = np.random.randn(3)

        frozen = getattr(self.calc, 'frozen_mask', None)
        if frozen is not None:
            axis[frozen] = 0.0

        norm = np.linalg.norm(axis)
        if norm > 0:
            axis /= norm
        return axis

    def _extract_mechanism(self, atoms, pos_min, pos_sp, axis_sp,
                           e_min, e_sp, geometry):
        """
        Extract a Mechanism from minimum, saddle point, and axis.

        Nudges from SP along ±axis, minimizes both sides, validates.
        """
        n = len(atoms)

        # Compute nudge displacement
        disp = pos_sp - pos_min
        disp_norm = np.linalg.norm(disp)
        if disp_norm < 1e-10:
            return None
        nudge = self.nudge_frac * disp_norm

        # Forward: SP + nudge * axis
        pos_fwd = pos_sp + nudge * axis_sp
        atoms_fwd = atoms.copy()
        atoms_fwd.calc = self.calc.calculator
        atoms_fwd.set_positions(pos_fwd)
        atoms_fwd, e_fwd, conv_fwd, _ = self.minimizer.minimize(atoms_fwd)
        pos_fwd = atoms_fwd.get_positions()

        # Reverse: SP - nudge * axis
        pos_rev = pos_sp - nudge * axis_sp
        atoms_rev = atoms.copy()
        atoms_rev.calc = self.calc.calculator
        atoms_rev.set_positions(pos_rev)
        atoms_rev, e_rev, conv_rev, _ = self.minimizer.minimize(atoms_rev)
        pos_rev = atoms_rev.get_positions()

        if not (conv_fwd and conv_rev):
            return None

        # Determine which is forward and which is reverse
        # Reverse should be close to initial minimum
        d_fwd_min = np.sqrt(np.mean(np.sum((pos_fwd - pos_min)**2, axis=1)))
        d_rev_min = np.sqrt(np.mean(np.sum((pos_rev - pos_min)**2, axis=1)))

        if d_fwd_min < d_rev_min:
            # Swap: what we called fwd is actually closer to initial
            pos_fwd, pos_rev = pos_rev, pos_fwd
            e_fwd, e_rev = e_rev, e_fwd

        # Validate: reverse should be close to initial minimum
        d_rev_min = np.sqrt(np.mean(np.sum((pos_rev - pos_min)**2, axis=1)))
        if d_rev_min > self.basin_tol:
            # Reverse is not the same basin — skip
            if self.debug:
                print(f"MASTER: Reverse min too far from initial: "
                      f"d={d_rev_min:.4f}")
            return None

        # Forward should be different from initial
        d_fwd_min = np.sqrt(np.mean(np.sum((pos_fwd - pos_min)**2, axis=1)))
        if d_fwd_min < self.basin_tol * 0.1:
            # Same basin — not a real transition
            return None

        # Build local displacements
        delta_sp = np.zeros((geometry.n_atoms, 3))
        delta_fwd = np.zeros((geometry.n_atoms, 3))
        mech_axis = np.zeros((geometry.n_atoms, 3))

        for j in range(geometry.n_atoms):
            gi = int(geometry.indices[j])
            delta_sp[j] = pos_sp[gi] - pos_min[gi]
            delta_fwd[j] = pos_fwd[gi] - pos_min[gi]
            mech_axis[j] = axis_sp[gi]

        barrier = e_sp - e_min
        delta_e = e_fwd - e_min

        # Estimate kinetic pre-factor (simplified TST)
        # Full version would use Hessian analysis at min and SP
        kBT_per_h = 6.25e12  # kB*300K/h in Hz (order of magnitude)
        kinetic_pre = kBT_per_h

        return Mechanism(
            barrier=barrier,
            delta=delta_e,
            kinetic_pre=kinetic_pre,
            delta_sp=delta_sp,
            delta_fwd=delta_fwd,
            axis=mech_axis,
        )

    def _is_new(self, mech, existing):
        """Check if mechanism is distinct from all existing ones."""
        for other in existing:
            if mech.n_atoms != other.n_atoms:
                continue

            rmsd_sp = np.sqrt(np.mean(np.sum(
                (mech.delta_sp - other.delta_sp)**2, axis=1)))
            rmsd_fwd = np.sqrt(np.mean(np.sum(
                (mech.delta_fwd - other.delta_fwd)**2, axis=1)))

            if rmsd_sp < self.mech_tol and rmsd_fwd < self.mech_tol:
                return False

            # Also check barrier similarity
            if abs(mech.barrier - other.barrier) < 0.01 and rmsd_fwd < self.mech_tol * 2:
                return False

        return True

    def _apply_symmetry(self, mech, O, perm):
        """
        Generate a mechanism copy by applying a symmetry operation.

        Parameters
        ----------
        mech : Mechanism
        O : ndarray, shape (3, 3)
            Rotation matrix.
        perm : list of int
            Atom permutation.

        Returns
        -------
        Mechanism
        """
        n = mech.n_atoms
        delta_sp = np.zeros_like(mech.delta_sp)
        delta_fwd = np.zeros_like(mech.delta_fwd)
        axis = np.zeros_like(mech.axis) if mech.axis is not None else None

        for j in range(n):
            delta_sp[j] = O @ mech.delta_sp[perm[j]]
            delta_fwd[j] = O @ mech.delta_fwd[perm[j]]
            if axis is not None:
                axis[j] = O @ mech.axis[perm[j]]

        return Mechanism(
            barrier=mech.barrier,
            delta=mech.delta,
            kinetic_pre=mech.kinetic_pre,
            delta_sp=delta_sp,
            delta_fwd=delta_fwd,
            axis=axis,
            err_fwd=mech.err_fwd,
            err_sp=mech.err_sp,
            poison_fwd=mech.poison_fwd,
            poison_sp=mech.poison_sp,
        )
