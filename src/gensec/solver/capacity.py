# ---------------------------------------------------------------------------
# GenSec — Copyright (c) 2026 Andrea Albero
#
# This file is part of GenSec.
#
# GenSec is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# GenSec is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public
# License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with GenSec.  If not, see <https://www.gnu.org/licenses/>.
# ---------------------------------------------------------------------------
r"""
Resistance surface generator — N-Mx-My interaction domain.

Phase 2: generates a 3D point cloud in (N, Mx, My) space by
scanning ultimate strain configurations across all curvature
directions.

For uniaxial analysis (chi_y=0), produces the classic N-M diagram.
"""

import numpy as np
from .integrator import FiberSolver


class NMDiagram:
    r"""
    Material-agnostic resistance domain generator.

    For **uniaxial** bending, generates the N-Mx interaction diagram
    (chi_y = 0). For **biaxial** bending, generates a 3D point cloud
    in (N, Mx, My) space by scanning curvature directions in the
    :math:`(\chi_x, \chi_y)` plane.

    The strain plane at any point is:

    .. math::

        \varepsilon(x, y) = \varepsilon_0
            + \chi_x\,(y - y_{\text{ref}})
            - \chi_y\,(x - x_{\text{ref}})

    For the 3D scan, the curvature direction angle :math:`\theta` is:

    .. math::

        \chi_x = \chi \cos\theta, \quad \chi_y = \chi \sin\theta

    and the curvature magnitude :math:`\chi` is determined by the
    strain limits at the extreme fibers for that direction.

    Parameters
    ----------
    solver : FiberSolver
    """

    def __init__(self, solver):
        self.solver = solver

    def _collect_strain_limits(self):
        """
        Global strain envelope from all materials.

        Returns
        -------
        eps_min_global, eps_max_global, eps_min_bulk, eps_max_bulk
        """
        sec = self.solver.sec
        b = sec.bulk_material
        emin_g, emax_g = b.eps_min, b.eps_max
        for r in sec.rebars:
            emin_g = min(emin_g, r.material.eps_min)
            emax_g = max(emax_g, r.material.eps_max)
        return emin_g, emax_g, b.eps_min, b.eps_max

    def _ultimate_strain_configs_1d(self, direction='x', n_points=200):
        r"""
        Generate ``(eps_inf, eps_sup)`` pairs for single-axis bending.

        For ``direction='x'`` (bending about x-axis, curvature
        :math:`\chi_x`), the section depth is ``H`` and the reference
        lever arm is ``y_ref``.

        For ``direction='y'`` (bending about y-axis, curvature
        :math:`\chi_y`), the section depth is ``B`` and the reference
        lever arm is ``x_ref``.

        Parameters
        ----------
        direction : ``'x'`` or ``'y'``
        n_points : int

        Returns
        -------
        list of (eps0, chi)
            Strain configurations as ``(eps0, curvature)`` pairs.
        """
        emg, exg, emb, _ = self._collect_strain_limits()
        sec = self.solver.sec

        if direction == 'x':
            depth = sec.H
            ref = self.solver.y_ref
        else:
            depth = sec.B
            ref = self.solver.x_ref

        configs = []

        def _to_eps0_chi(ei, es):
            """Convert edge strains to (eps0, chi)."""
            chi = (es - ei) / depth if depth > 0 else 0.0
            eps0 = ei + chi * ref
            return eps0, chi

        # Positive curvature: "top" (max coord) compressed
        for ei in np.linspace(exg, emg, n_points):
            configs.append(_to_eps0_chi(ei, emb))
        for es in np.linspace(emb, emb * 0.5, n_points // 2):
            configs.append(_to_eps0_chi(emb, es))

        # Negative curvature: "bottom" (min coord) compressed
        for es in np.linspace(exg, emg, n_points):
            configs.append(_to_eps0_chi(emb, es))
        for ei in np.linspace(emb, emb * 0.5, n_points // 2):
            configs.append(_to_eps0_chi(ei, emb))

        # Pure tension
        for ev in np.linspace(0, exg, n_points // 4):
            configs.append(_to_eps0_chi(ev, 0.0))
            configs.append(_to_eps0_chi(0.0, ev))
            configs.append(_to_eps0_chi(exg, ev))
            configs.append(_to_eps0_chi(ev, exg))

        # Near-uniform compression
        for ev in np.linspace(emb * 0.5, emb, n_points // 4):
            configs.append(_to_eps0_chi(ev, ev))
            span = abs(emb) * 0.15
            for d in np.linspace(-span, span, 5):
                configs.append(_to_eps0_chi(ev + d, ev - d))

        return configs

    def generate(self, n_points=300, direction='x'):
        r"""
        Generate the N-M interaction diagram for a single bending
        direction.

        For ``direction='x'``: N-Mx diagram (:math:`\chi_y = 0`).
        For ``direction='y'``: N-My diagram (:math:`\chi_x = 0`).

        Parameters
        ----------
        n_points : int, optional
            Base resolution. Default 300.
        direction : ``'x'`` or ``'y'``, optional
            Bending direction. Default ``'x'``.

        Returns
        -------
        dict
            ``N``, ``M`` [N, N*mm], ``N_kN``, ``M_kNm``.
            Also ``Mx`` or ``My`` alias depending on direction.
        """
        configs = self._ultimate_strain_configs_1d(direction, n_points)
        Nl, Ml = [], []
        for eps0, chi in configs:
            if direction == 'x':
                N, Mx, My = self.solver.integrate(eps0, chi, 0.0)
                Nl.append(N)
                Ml.append(Mx)
            else:
                N, Mx, My = self.solver.integrate(eps0, 0.0, chi)
                Nl.append(N)
                Ml.append(My)
        Na, Ma = np.array(Nl), np.array(Ml)
        M_key = "Mx" if direction == 'x' else "My"
        return {
            "N": Na, "M": Ma, M_key: Ma,
            "N_kN": Na / 1e3, "M_kNm": Ma / 1e6,
            f"{M_key}_kNm": Ma / 1e6,
            "direction": direction,
        }

    def generate_biaxial(self, n_angles=72, n_points_per_angle=200):
        r"""
        Generate the 3D resistance surface (N, Mx, My).

        Scans curvature directions :math:`\theta` in
        :math:`[0, 2\pi)` and, for each direction, scans curvature
        magnitudes through ultimate strain configurations.

        Parameters
        ----------
        n_angles : int, optional
            Number of curvature direction angles. Default 36 (every 10°).
        n_points_per_angle : int, optional
            Strain configurations per angle. Default 200.

        Returns
        -------
        dict
            ``N`` [N], ``Mx``, ``My`` [N*mm], and ``_kN``/``_kNm``
            variants.
        """
        sec = self.solver.sec
        emg, exg, emb, _ = self._collect_strain_limits()

        # All bulk fiber positions relative to reference
        lx = sec.x_fibers - self.solver.x_ref
        ly = sec.y_fibers - self.solver.y_ref

        # Also rebar positions
        lx_r = sec.x_rebars - self.solver.x_ref
        ly_r = sec.y_rebars - self.solver.y_ref

        # All lever arms combined
        all_lx = np.concatenate([lx, lx_r])
        all_ly = np.concatenate([ly, ly_r])

        thetas = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)

        Nl, Mxl, Myl = [], [], []

        for theta in thetas:
            cos_t, sin_t = np.cos(theta), np.sin(theta)

            # Project all fiber positions onto curvature direction.
            # Strain field: eps = eps0 + chi_x*ly - chi_y*lx
            # With chi_x = chi*cos(t), chi_y = chi*sin(t):
            # eps = eps0 + chi*(ly*cos_t - lx*sin_t)
            proj = all_ly * cos_t - all_lx * sin_t

            # Extreme projections
            p_max = proj.max()
            p_min = proj.min()

            if abs(p_max - p_min) < 1e-10:
                continue

            # Scan strain configurations: similar to uniaxial but
            # along the projected direction.
            # The "top" extreme (max projection) gets eps_top,
            # the "bottom" extreme (min projection) gets eps_bot.
            # eps0 + chi * p_max = eps_top
            # eps0 + chi * p_min = eps_bot
            # => chi = (eps_top - eps_bot) / (p_max - p_min)
            # => eps0 = eps_bot - chi * p_min

            span = p_max - p_min

            def _scan(eps_bot, eps_top):
                chi = (eps_top - eps_bot) / span
                eps0 = eps_bot - chi * p_min
                chi_x = chi * cos_t
                chi_y = chi * sin_t
                return self.solver.integrate(eps0, chi_x, chi_y)

            n = n_points_per_angle

            # Branch 1: "top" at bulk crush limit
            for eb in np.linspace(exg, emg, n):
                N, Mx, My = _scan(eb, emb)
                Nl.append(N); Mxl.append(Mx); Myl.append(My)

            # Branch 2: "bottom" at bulk crush limit
            for et in np.linspace(exg, emg, n):
                N, Mx, My = _scan(emb, et)
                Nl.append(N); Mxl.append(Mx); Myl.append(My)

            # Branch 3: near-uniform compression
            for ev in np.linspace(emb * 0.5, emb, n // 4):
                N, Mx, My = _scan(ev, ev)
                Nl.append(N); Mxl.append(Mx); Myl.append(My)
                sp = abs(emb) * 0.15
                for d in np.linspace(-sp, sp, 3):
                    N, Mx, My = _scan(ev + d, ev - d)
                    Nl.append(N); Mxl.append(Mx); Myl.append(My)

            # Branch 4: pure tension
            for ev in np.linspace(0, exg, n // 4):
                N, Mx, My = _scan(ev, 0.0)
                Nl.append(N); Mxl.append(Mx); Myl.append(My)
                N, Mx, My = _scan(0.0, ev)
                Nl.append(N); Mxl.append(Mx); Myl.append(My)

        Na = np.array(Nl)
        Mxa = np.array(Mxl)
        Mya = np.array(Myl)

        return {
            "N": Na, "Mx": Mxa, "My": Mya,
            "N_kN": Na / 1e3,
            "Mx_kNm": Mxa / 1e6,
            "My_kNm": Mya / 1e6,
        }

    def generate_mx_my(self, N_fixed, n_angles=72,
                       n_points_per_angle=200):
        r"""
        Generate the Mx-My interaction contour at a fixed axial force.

        For each curvature direction :math:`\theta \in [0, 2\pi)`, the
        algorithm scans curvature magnitudes :math:`\chi` from 0 to the
        ultimate value, at each step solving for :math:`\varepsilon_0`
        such that :math:`N = N_{\text{fixed}}`.  The envelope point
        is the :math:`(M_x, M_y)` at the **maximum** curvature before
        any fiber reaches a material strain limit.

        This approach is **direct** (no interpolation of strain
        configurations) and relies on the robust
        :meth:`_solve_eps0_for_N` for the 1-D equilibrium at each
        step.

        Parameters
        ----------
        N_fixed : float
            Fixed axial force [N].
        n_angles : int, optional
            Number of curvature directions. Default 72 (every 5°).
        n_points_per_angle : int, optional
            Curvature steps per direction. Default 200.

        Returns
        -------
        dict
            ``Mx`` [N*mm], ``My`` [N*mm], ``Mx_kNm``, ``My_kNm``,
            ``N_fixed_kN``.
        """
        sec = self.solver.sec
        emg, exg, emb, _ = self._collect_strain_limits()

        thetas = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
        Mxl, Myl = [], []

        # Warm-start eps0 across angles
        eps0_warm = 0.0

        for theta in thetas:
            cos_t, sin_t = np.cos(theta), np.sin(theta)

            # Estimate chi_max from section depth in this direction
            chi_max_candidates = []
            if abs(cos_t) > 0.01 and sec.H > 0:
                chi_max_candidates.append(
                    abs(emb) / (sec.H * 0.25) * 1.5)
            if abs(sin_t) > 0.01 and sec.B > 0:
                chi_max_candidates.append(
                    abs(emb) / (sec.B * 0.25) * 1.5)
            chi_max = (max(chi_max_candidates)
                       if chi_max_candidates else 1e-4)

            # Scan chi from 0 to chi_max, solving for eps0 at each step
            best_Mx = 0.0
            best_My = 0.0
            best_M_mag = 0.0
            eps0_prev = eps0_warm

            for chi_mag in np.linspace(0, chi_max,
                                       n_points_per_angle):
                chi_x = chi_mag * cos_t
                chi_y = chi_mag * sin_t

                eps0 = self._solve_eps0_for_N(
                    self.solver, N_fixed, chi_x, chi_y,
                    eps0_prev, emb)
                eps0_prev = eps0

                # Check material limits
                eb, er = self.solver.strain_field(
                    eps0, chi_x, chi_y)
                all_eps = (np.concatenate([eb, er])
                           if len(er) > 0 else eb)

                if (all_eps.min() <= emb * 0.99
                        or all_eps.max() >= exg * 0.99):
                    # At or past ultimate — record this last point
                    # and stop scanning further
                    N, Mx, My = self.solver.integrate(
                        eps0, chi_x, chi_y)
                    M_mag = np.sqrt(Mx**2 + My**2)
                    if M_mag > best_M_mag:
                        best_M_mag = M_mag
                        best_Mx = Mx
                        best_My = My
                    break

                # Within limits — check if this is the best so far
                N, Mx, My = self.solver.integrate(
                    eps0, chi_x, chi_y)
                M_mag = np.sqrt(Mx**2 + My**2)
                if M_mag > best_M_mag:
                    best_M_mag = M_mag
                    best_Mx = Mx
                    best_My = My

            Mxl.append(best_Mx)
            Myl.append(best_My)
            eps0_warm = eps0_prev

        Mxa = np.array(Mxl)
        Mya = np.array(Myl)

        return {
            "Mx": Mxa, "My": Mya,
            "Mx_kNm": Mxa / 1e6, "My_kNm": Mya / 1e6,
            "N_fixed_kN": N_fixed / 1e3,
        }

    def generate_polar_ductility(self, N_fixed, n_angles=72,
                                 n_points=400):
        r"""
        Compute ultimate curvature as a function of bending direction.

        For each curvature direction :math:`\theta` in the
        :math:`(\chi_x, \chi_y)` plane, the section is loaded
        incrementally (increasing curvature magnitude) at constant
        axial force :math:`N = N_{\text{fixed}}`.  The ultimate
        curvature :math:`\chi_u(\theta)` is the **largest** curvature
        at which all fibers remain within their material strain
        limits.

        Once the coarse scan identifies the failure step, a
        **bisection refinement** between the last safe step and the
        first failed step pins down the ultimate curvature to high
        accuracy.

        Parameters
        ----------
        N_fixed : float
            Axial force [N].
        n_angles : int, optional
            Number of curvature directions. Default 72.
        n_points : int, optional
            Curvature steps per direction. Default 400.

        Returns
        -------
        dict
            ``thetas`` [rad], ``chi_u`` [1/mm],
            ``chi_u_km`` [1/km], ``N_fixed_kN``.
        """
        sec = self.solver.sec
        emg, exg, emb, _ = self._collect_strain_limits()

        thetas = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
        chi_ultimates = np.zeros(n_angles)

        eps0_warm = 0.0

        for i, theta in enumerate(thetas):
            cos_t, sin_t = np.cos(theta), np.sin(theta)

            # Estimate chi_max from section depth in this direction
            chi_max_candidates = []
            if abs(cos_t) > 0.01 and sec.H > 0:
                chi_max_candidates.append(
                    abs(emb) / (sec.H * 0.25) * 1.5)
            if abs(sin_t) > 0.01 and sec.B > 0:
                chi_max_candidates.append(
                    abs(emb) / (sec.B * 0.25) * 1.5)
            chi_max = (max(chi_max_candidates)
                       if chi_max_candidates else 1e-4)

            # --- Coarse scan: find the first chi that violates limits ---
            chi_safe = 0.0
            chi_fail = None
            eps0_prev = eps0_warm

            for chi_mag in np.linspace(0, chi_max, n_points):
                chi_x = chi_mag * cos_t
                chi_y = chi_mag * sin_t

                eps0 = self._solve_eps0_for_N(
                    self.solver, N_fixed, chi_x, chi_y,
                    eps0_prev, emb)
                eps0_prev = eps0

                eb, er = self.solver.strain_field(
                    eps0, chi_x, chi_y)
                all_eps = (np.concatenate([eb, er])
                           if len(er) > 0 else eb)

                if (all_eps.min() <= emb * 0.99
                        or all_eps.max() >= exg * 0.99):
                    chi_fail = chi_mag
                    break
                chi_safe = chi_mag

            # --- Bisection refinement between chi_safe and chi_fail ---
            if chi_fail is not None and chi_fail > chi_safe:
                a, b = chi_safe, chi_fail
                for _ in range(20):
                    mid = 0.5 * (a + b)
                    chi_x = mid * cos_t
                    chi_y = mid * sin_t
                    eps0 = self._solve_eps0_for_N(
                        self.solver, N_fixed, chi_x, chi_y,
                        eps0_prev, emb)
                    eb, er = self.solver.strain_field(
                        eps0, chi_x, chi_y)
                    all_eps = (np.concatenate([eb, er])
                               if len(er) > 0 else eb)
                    if (all_eps.min() <= emb * 0.99
                            or all_eps.max() >= exg * 0.99):
                        b = mid
                    else:
                        a = mid
                        eps0_prev = eps0
                chi_ultimates[i] = a
            else:
                # Never failed → chi_max is conservative estimate
                chi_ultimates[i] = chi_safe

            eps0_warm = eps0_prev

        return {
            "thetas": thetas,
            "chi_u": chi_ultimates,
            "chi_u_km": chi_ultimates * 1e6,
            "N_fixed_kN": N_fixed / 1e3,
        }
    def generate_moment_curvature(self, N_fixed, chi_max=None,
                                  n_points=200, direction='x'):
        r"""
        Generate the moment-curvature diagram at fixed axial force.

        Scans curvature :math:`\chi` from 0 to the ultimate value
        (where a material strain limit is reached), tracking the
        moment response. Key points (first yield, ultimate) are
        identified.

        Parameters
        ----------
        N_fixed : float
            Fixed axial force [N]. Negative = compression.
        chi_max : float, optional
            Maximum curvature to scan [1/mm]. If ``None``, computed
            automatically from strain limits.
        n_points : int, optional
            Number of curvature steps. Default 200.
        direction : str, optional
            ``'x'`` for Mx-chi_x (default) or ``'y'`` for My-chi_y.

        Returns
        -------
        dict
            ``chi`` [1/mm], ``M`` [N*mm], ``M_kNm``, ``chi_km``
            [1/km], ``eps_min``, ``eps_max`` (extreme strains at each
            step), ``yield_index`` (first yield step or None),
            ``ultimate_index`` (ultimate step or None),
            ``N_fixed_kN``, ``direction``.
        """
        sec = self.solver.sec
        emg, _, emb, _ = self._collect_strain_limits()

        # Determine chi_max from geometry and strain limits.
        # Use the full section height (not half) as lever arm to
        # ensure the ultimate strain is reached even when the
        # neutral axis is far from the centroid.
        if chi_max is None:
            if direction == 'x':
                d_max = sec.H
            else:
                d_max = sec.B
            if d_max > 0:
                chi_max = abs(emb) / (d_max * 0.3) * 1.5
            else:
                chi_max = 1e-4

        # Collect rebar yield strains for yield detection
        rebar_eps_yd = []
        for rb in sec.rebars:
            if hasattr(rb.material, 'eps_yd'):
                rebar_eps_yd.append(rb.material.eps_yd)

        # Scan both positive and negative curvature
        results_pos = self._scan_chi(
            N_fixed, 0, chi_max, n_points, direction, rebar_eps_yd)
        results_neg = self._scan_chi(
            N_fixed, 0, -chi_max, n_points, direction, rebar_eps_yd)

        # Merge: negative reversed + positive
        chi_all = np.concatenate([
            results_neg["chi"][::-1], results_pos["chi"][1:]])
        M_all = np.concatenate([
            results_neg["M"][::-1], results_pos["M"][1:]])
        eps_min_all = np.concatenate([
            results_neg["eps_min"][::-1], results_pos["eps_min"][1:]])
        eps_max_all = np.concatenate([
            results_neg["eps_max"][::-1], results_pos["eps_max"][1:]])

        return {
            "chi": chi_all,
            "M": M_all,
            "M_kNm": M_all / 1e6,
            "chi_km": chi_all * 1e6,  # 1/mm -> 1/km
            "eps_min": eps_min_all,
            "eps_max": eps_max_all,
            "N_fixed_kN": N_fixed / 1e3,
            "direction": direction,
            "yield_chi_pos": results_pos.get("yield_chi"),
            "yield_M_pos": results_pos.get("yield_M"),
            "ultimate_chi_pos": results_pos.get("ultimate_chi"),
            "ultimate_M_pos": results_pos.get("ultimate_M"),
            "yield_chi_neg": results_neg.get("yield_chi"),
            "yield_M_neg": results_neg.get("yield_M"),
            "ultimate_chi_neg": results_neg.get("ultimate_chi"),
            "ultimate_M_neg": results_neg.get("ultimate_M"),
        }

    def _scan_chi(self, N_fixed, chi_start, chi_end, n_points,
                  direction, rebar_eps_yd):
        """
        Internal: scan curvature from chi_start to chi_end,
        solving for eps0 at each step to maintain N = N_fixed.

        Uses Newton iteration with warm-start from previous step,
        falling back to bisection if Newton fails.
        """
        sec = self.solver.sec
        emg, exg, emb, _ = self._collect_strain_limits()
        sv = self.solver

        chis = np.linspace(chi_start, chi_end, n_points)
        Ms = np.zeros(n_points)
        eps_mins = np.zeros(n_points)
        eps_maxs = np.zeros(n_points)

        yield_chi = None
        yield_M = None
        ultimate_chi = None
        ultimate_M = None

        # Initial eps0 estimate from the equilibrium solver at chi=0
        if abs(chi_start) < 1e-15:
            # Use a rough estimate: N_fixed ~ sigma_mean * A
            A_gross = getattr(sec, 'gross_area', sec.B * sec.H)
            eps0_guess = N_fixed / (A_gross * 15000)  # rough Ec/2
            eps0_guess = np.clip(eps0_guess, emb, -emb)
        else:
            eps0_guess = 0.0

        for k, chi in enumerate(chis):
            if direction == 'x':
                chi_x, chi_y = chi, 0.0
            else:
                chi_x, chi_y = 0.0, chi

            # Solve for eps0 such that N(eps0, chi_x, chi_y) = N_fixed
            eps0 = self._solve_eps0_for_N(
                sv, N_fixed, chi_x, chi_y, eps0_guess, emb)

            N, Mx, My = sv.integrate(eps0, chi_x, chi_y)
            M = Mx if direction == 'x' else My
            Ms[k] = M
            eps0_guess = eps0  # warm-start

            # Track extreme strains
            eb, er = sv.strain_field(eps0, chi_x, chi_y)
            all_eps = np.concatenate([eb, er]) if len(er) > 0 else eb
            eps_mins[k] = all_eps.min()
            eps_maxs[k] = all_eps.max()

            # Detect first yield
            if yield_chi is None and len(rebar_eps_yd) > 0 and abs(chi) > 0:
                for j, rb in enumerate(sec.rebars):
                    if hasattr(rb.material, 'eps_yd'):
                        if abs(er[j]) >= rb.material.eps_yd * 0.99:
                            yield_chi = chi
                            yield_M = M
                            break

            # Detect ultimate (any fiber at material limit).
            # emb is negative: emb * 0.99 is less negative, so the
            # check triggers when strain reaches (or just exceeds)
            # the limit, not 1% past it.
            if ultimate_chi is None and abs(chi) > 0:
                if (all_eps.min() <= emb * 0.99
                        or all_eps.max() >= exg * 0.99):
                    ultimate_chi = chi
                    ultimate_M = M

        return {
            "chi": chis, "M": Ms,
            "eps_min": eps_mins, "eps_max": eps_maxs,
            "yield_chi": yield_chi, "yield_M": yield_M,
            "ultimate_chi": ultimate_chi, "ultimate_M": ultimate_M,
        }

    @staticmethod
    def _solve_eps0_for_N(sv, N_target, chi_x, chi_y, eps0_init, emb,
                          tol=1.0):
        r"""
        Solve for :math:`\varepsilon_0` at fixed curvature such that
        :math:`N = N_{\text{target}}`.

        The function :math:`N(\varepsilon_0)` at fixed
        :math:`(\chi_x, \chi_y)` is generally non-decreasing in the
        working strain range, but can be non-monotone near the
        material strain limits.

        Strategy:

        1.  **Newton** with warm-start from ``eps0_init``.  Fast when
            the initial guess is close (typical in sequential scans).
        2.  **Scan + bisect fallback**: if Newton fails, scan
            :math:`N(\varepsilon_0)` on a coarse grid, find a local
            crossing bracket, and bisect within it.

        Parameters
        ----------
        sv : FiberSolver
        N_target : float
            Target axial force [N].
        chi_x, chi_y : float
            Fixed curvatures [1/mm].
        eps0_init : float
            Warm-start for Newton.
        emb : float
            Concrete ultimate strain (``eps_cu2``, negative).
        tol : float, optional
            Force tolerance [N].  Default 1.0.

        Returns
        -------
        float
            The :math:`\varepsilon_0` value that satisfies equilibrium.
        """
        # ---- Phase 1: Newton with warm-start ----
        eps0 = eps0_init
        for _ in range(30):
            N, _, _ = sv.integrate(eps0, chi_x, chi_y)
            r = N - N_target
            if abs(r) < tol:
                return eps0
            de = max(abs(eps0) * 1e-7, 1e-9)
            N1, _, _ = sv.integrate(eps0 + de, chi_x, chi_y)
            dNde = (N1 - N) / de
            if abs(dNde) > 1.0:
                step = -r / dNde
                # Adaptive clamp: allow larger steps for larger strains
                max_step = max(abs(eps0) * 0.5, 5e-4)
                step = np.clip(step, -max_step, max_step)
                eps0 += step
            else:
                # Near-zero stiffness — nudge toward compression
                eps0 -= 1e-5
                continue

        # ---- Phase 2: scan + bisect fallback ----
        # Build strain range from material limits
        sec = sv.sec
        eps_lo = sec.bulk_material.eps_min
        eps_hi = 0.0
        for rb in sec.rebars:
            eps_lo = min(eps_lo, rb.material.eps_min)
            eps_hi = max(eps_hi, rb.material.eps_max)
        eps_lo *= 1.01
        eps_hi *= 1.01

        n_scan = 80
        eps_scan = np.linspace(eps_lo, eps_hi, n_scan)
        N_scan = np.empty(n_scan)
        for k in range(n_scan):
            N_scan[k], _, _ = sv.integrate(eps_scan[k], chi_x, chi_y)

        # Find all crossings
        crossings = []
        for k in range(n_scan - 1):
            if ((N_scan[k] - N_target) * (N_scan[k + 1] - N_target) <= 0
                    and abs(N_scan[k + 1] - N_scan[k]) > 0.01):
                crossings.append(k)

        if not crossings:
            # Return best Newton result as last resort
            return eps0

        # Pick crossing closest to eps0=0 (working region)
        best_k = min(crossings,
                     key=lambda k: abs(eps_scan[k] + eps_scan[k + 1]))

        a = eps_scan[best_k]
        b = eps_scan[best_k + 1]
        Na = N_scan[best_k]

        for _ in range(60):
            mid = 0.5 * (a + b)
            Nm, _, _ = sv.integrate(mid, chi_x, chi_y)
            if abs(Nm - N_target) < tol:
                return mid
            if (Nm - N_target) * (Na - N_target) <= 0:
                b = mid
            else:
                a = mid
                Na = Nm

        return 0.5 * (a + b)
