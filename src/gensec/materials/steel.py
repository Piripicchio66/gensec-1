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
Reinforcing steel constitutive law — elastic-plastic with optional hardening.
"""

import numpy as np
from dataclasses import dataclass
from .base import Material


@dataclass
class Steel(Material):
    r"""
    Elastic-plastic steel with optional linear hardening.

    .. math::

        \sigma_s(\varepsilon) =
        \begin{cases}
            E_s\,\varepsilon
                & |\varepsilon| \le \varepsilon_{yd} \\
            \mathrm{sign}(\varepsilon)\!\left[f_{yd}
                + (f_{td}-f_{yd})
                \dfrac{|\varepsilon|-\varepsilon_{yd}}
                      {\varepsilon_{su}-\varepsilon_{yd}}\right]
                & \varepsilon_{yd} < |\varepsilon| \le \varepsilon_{su} \\
            0 & |\varepsilon| > \varepsilon_{su}
        \end{cases}

    Parameters
    ----------
    fyk : float
        Characteristic yield strength [MPa].
    gamma_s : float, optional
        Partial safety factor. Default 1.15.
    Es : float, optional
        Young's modulus [MPa]. Default 200000.
    k_hardening : float, optional
        :math:`f_t/f_y` ratio. Default 1.0 (perfectly plastic).
    eps_su : float, optional
        Ultimate strain. Default 0.01.
    works_in_compression : bool, optional
        If False, compressive stress is zero. Default True.

    Attributes
    ----------
    fyd : float
        Design yield strength [MPa].
    ftd : float
        Design ultimate strength [MPa].
    eps_yd : float
        Yield strain.
    """

    fyk: float = 450.0
    gamma_s: float = 1.15
    Es: float = 200000.0
    k_hardening: float = 1.0
    eps_su: float = 0.01
    works_in_compression: bool = True

    def __post_init__(self):
        self.fyd = self.fyk / self.gamma_s
        self.ftd = self.fyd * self.k_hardening
        self.eps_yd = self.fyd / self.Es

    @property
    def eps_min(self):
        return -self.eps_su if self.works_in_compression else 0.0

    @property
    def eps_max(self):
        return self.eps_su

    def stress(self, eps):
        if eps >= 0:
            if eps <= self.eps_yd:
                return self.Es * eps
            elif eps <= self.eps_su:
                return self.fyd + (self.ftd - self.fyd) * (
                    (eps - self.eps_yd) / (self.eps_su - self.eps_yd))
            return 0.0
        if not self.works_in_compression:
            return 0.0
        ea = abs(eps)
        if ea <= self.eps_yd:
            return self.Es * eps
        elif ea <= self.eps_su:
            return -(self.fyd + (self.ftd - self.fyd) * (
                (ea - self.eps_yd) / (self.eps_su - self.eps_yd)))
        return 0.0

    def stress_array(self, eps):
        sigma = np.zeros_like(eps)
        m_te = (eps >= 0) & (eps <= self.eps_yd)
        sigma[m_te] = self.Es * eps[m_te]
        m_tp = (eps > self.eps_yd) & (eps <= self.eps_su)
        sigma[m_tp] = self.fyd + (self.ftd - self.fyd) * (
            (eps[m_tp] - self.eps_yd) / (self.eps_su - self.eps_yd))
        if self.works_in_compression:
            ea = np.abs(eps)
            m_ce = (eps < 0) & (ea <= self.eps_yd)
            sigma[m_ce] = self.Es * eps[m_ce]
            m_cp = (eps < 0) & (ea > self.eps_yd) & (ea <= self.eps_su)
            sigma[m_cp] = -(self.fyd + (self.ftd - self.fyd) * (
                (ea[m_cp] - self.eps_yd) / (self.eps_su - self.eps_yd)))
        return sigma
