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
"""
Tabulated material — arbitrary stress-strain curve from data points.
"""

import numpy as np
from .base import Material


class TabulatedMaterial(Material):
    r"""
    Material defined by a tabulated :math:`(\varepsilon, \sigma)` curve.

    Piecewise-linear interpolation between data points. Strains
    outside the table range return zero stress (material failure).

    Parameters
    ----------
    strains : array_like
        Strictly increasing strain values.
    stresses : array_like
        Corresponding stress values [MPa].
    name : str, optional
        Label. Default ``"Tabulated"``.

    Raises
    ------
    ValueError
        Bad dimensions or non-monotonic strains.

    Examples
    --------
    >>> mat = TabulatedMaterial([-0.01, -0.002, 0, 0.002, 0.01],
    ...                        [-400, -400, 0, 400, 400], "EPP")
    >>> mat.stress(0.001)
    200.0
    """

    def __init__(self, strains, stresses, name="Tabulated"):
        strains = np.asarray(strains, dtype=float)
        stresses = np.asarray(stresses, dtype=float)
        if len(strains) != len(stresses):
            raise ValueError("strains/stresses length mismatch")
        if len(strains) < 2:
            raise ValueError("Need >= 2 points")
        if not np.all(np.diff(strains) > 0):
            raise ValueError("strains must be strictly increasing")
        self._strains = strains
        self._stresses = stresses
        self._name = name

    @property
    def eps_min(self):
        return float(self._strains[0])

    @property
    def eps_max(self):
        return float(self._strains[-1])

    @property
    def name(self):
        """Material label."""
        return self._name

    def stress(self, eps):
        if eps < self._strains[0] or eps > self._strains[-1]:
            return 0.0
        return float(np.interp(eps, self._strains, self._stresses))

    def stress_array(self, eps):
        sigma = np.interp(eps, self._strains, self._stresses)
        sigma[eps < self._strains[0]] = 0.0
        sigma[eps > self._strains[-1]] = 0.0
        return sigma
