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
Demand verification against the resistance domain.

Supports both 2D (N-M) and 3D (N-Mx-My) domains using
:class:`scipy.spatial.ConvexHull`.
"""

import numpy as np
from scipy.spatial import ConvexHull


class DemandChecker:
    r"""
    Verify load demands against the resistance domain.

    Works in both 2D ``(N, M)`` and 3D ``(N, Mx, My)`` mode,
    depending on the data passed at construction.

    The **utilization ratio** :math:`\eta` is computed by ray-casting
    from a reference point through the demand point to the hull
    boundary:

    .. math::

        \eta = \frac{|\mathbf{d} - \mathbf{d}_{\text{ref}}|}
                    {|\mathbf{r} - \mathbf{d}_{\text{ref}}|}

    Parameters
    ----------
    nm_data : dict
        Output of :meth:`NMDiagram.generate` (2D, keys ``N``, ``M``)
        or :meth:`NMDiagram.generate_biaxial` (3D, keys ``N``, ``Mx``,
        ``My``).

    Attributes
    ----------
    ndim : int
        Dimensionality (2 or 3).
    hull : scipy.spatial.ConvexHull
    """

    def __init__(self, nm_data):
        if "My" in nm_data and nm_data.get("My") is not None:
            Mx = nm_data.get("Mx", nm_data.get("M"))
            pts = np.column_stack([nm_data["N"], Mx, nm_data["My"]])
            self.ndim = 3
        else:
            M = nm_data.get("M", nm_data.get("Mx"))
            pts = np.column_stack([M, nm_data["N"]])
            self.ndim = 2

        self.hull = ConvexHull(pts)
        self._equations = self.hull.equations

    def is_inside(self, N, Mx_or_M, My=None):
        """
        Check if a point is inside the resistance domain.

        Parameters
        ----------
        N : float [N]
        Mx_or_M : float [N*mm]
        My : float or None [N*mm]
            Required for 3D mode.

        Returns
        -------
        bool
        """
        p = self._make_point(N, Mx_or_M, My)
        return bool(np.all(
            self._equations[:, :-1] @ p + self._equations[:, -1] <= 1e-6
        ))

    def utilization_ratio(self, N, Mx_or_M, My=None,
                          N_ref=0.0, M_ref=0.0, My_ref=0.0):
        """
        Utilization ratio by ray-casting from reference to demand.

        Parameters
        ----------
        N, Mx_or_M : float
            Demand forces [N, N*mm].
        My : float or None
        N_ref, M_ref, My_ref : float, optional
            Reference point. Default: origin.

        Returns
        -------
        float
            Utilization ratio. < 1 = verified.
        """
        p = self._make_point(N, Mx_or_M, My)
        o = self._make_point(N_ref, M_ref, My_ref)
        direction = p - o
        length = np.linalg.norm(direction)
        if length < 1e-12:
            return 0.0 if self.is_inside(N, Mx_or_M, My) else float('inf')

        A = self._equations[:, :-1]
        c = self._equations[:, -1]
        denom = A @ direction
        numer = -(A @ o + c)
        valid = np.abs(denom) > 1e-15
        t = np.full_like(denom, np.inf)
        t[valid] = numer[valid] / denom[valid]
        t[t <= 1e-12] = np.inf
        t_bnd = float(np.min(t))

        if t_bnd == np.inf:
            return float('inf') if not self.is_inside(N, Mx_or_M, My) else 0.0
        return 1.0 / t_bnd

    def check_demands(self, demands, N_ref=0.0, M_ref=0.0, My_ref=0.0):
        """
        Batch verification of demands.

        Parameters
        ----------
        demands : list of dict
            Each dict: ``name``, ``N`` [N], ``M`` or ``Mx`` [N*mm],
            optionally ``My`` [N*mm].
        N_ref, M_ref, My_ref : float

        Returns
        -------
        list of dict
            Each: ``name``, ``N_kN``, ``M_kNm`` (and ``My_kNm``),
            ``inside``, ``utilization``, ``verified``.
        """
        results = []
        for d in demands:
            N = d["N"]
            Mx = d.get("Mx", d.get("M", 0.0))
            My = d.get("My", 0.0)
            inside = self.is_inside(N, Mx, My if self.ndim == 3 else None)
            eta = self.utilization_ratio(
                N, Mx, My if self.ndim == 3 else None,
                N_ref, M_ref, My_ref)
            r = {
                "name": d["name"],
                "N_kN": N / 1e3,
                "M_kNm": Mx / 1e6,
                "inside": inside,
                "utilization": eta,
                "verified": eta <= 1.0,
            }
            if self.ndim == 3:
                r["My_kNm"] = My / 1e6
            results.append(r)
        return results

    def _make_point(self, N, Mx_or_M, My=None):
        """Build point array consistent with hull dimensionality."""
        if self.ndim == 3:
            return np.array([N, Mx_or_M, My if My is not None else 0.0])
        else:
            return np.array([Mx_or_M, N])
