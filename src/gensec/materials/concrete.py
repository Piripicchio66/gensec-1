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
Concrete constitutive law — parabola-rectangle (EC2 3.1.7).
"""

import numpy as np
from dataclasses import dataclass
from .base import Material


@dataclass
class Concrete(Material):
    r"""
    Parabola-rectangle concrete (EC2 3.1.7 / NTC 2018).

    .. math::

        \sigma_c(\varepsilon) =
        \begin{cases}
            0 & \varepsilon > 0 \\
            -f_{cd}\!\left[1 - \left(1 - \dfrac{\varepsilon}
                {\varepsilon_{c2}}\right)^{\!n}\right]
                & \varepsilon_{c2} \le \varepsilon \le 0 \\
            -f_{cd}
                & \varepsilon_{cu2} \le \varepsilon < \varepsilon_{c2} \\
            0 & \varepsilon < \varepsilon_{cu2}
        \end{cases}

    where :math:`f_{cd} = \alpha_{cc}\,f_{ck}/\gamma_c`.

    Parameters
    ----------
    fck : float
        Characteristic cylinder strength [MPa].
    gamma_c : float, optional
        Partial safety factor. Default 1.5.
    alpha_cc : float, optional
        Long-term coefficient. Default 0.85.
    n_parabola : float, optional
        Parabolic exponent. Default 2.0.
    eps_c2 : float, optional
        Peak-stress strain. Default -0.002.
    eps_cu2 : float, optional
        Ultimate compressive strain. Default -0.0035.

    Attributes
    ----------
    fcd : float
        Design compressive strength [MPa].
    """

    fck: float = 25.0
    gamma_c: float = 1.5
    alpha_cc: float = 0.85
    n_parabola: float = 2.0
    eps_c2: float = -0.002
    eps_cu2: float = -0.0035

    def __post_init__(self):
        self.fcd = self.alpha_cc * self.fck / self.gamma_c

    @property
    def eps_min(self):
        return self.eps_cu2

    @property
    def eps_max(self):
        return 0.0

    def stress(self, eps):
        if eps > 0 or eps < self.eps_cu2:
            return 0.0
        if eps >= self.eps_c2:
            eta = eps / self.eps_c2
            return -self.fcd * (1.0 - (1.0 - eta) ** self.n_parabola)
        return -self.fcd

    def stress_array(self, eps):
        sigma = np.zeros_like(eps)
        m1 = (eps <= 0) & (eps >= self.eps_c2)
        eta = eps[m1] / self.eps_c2
        sigma[m1] = -self.fcd * (1.0 - (1.0 - eta) ** self.n_parabola)
        m2 = (eps < self.eps_c2) & (eps >= self.eps_cu2)
        sigma[m2] = -self.fcd
        return sigma
