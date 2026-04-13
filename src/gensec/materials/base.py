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
Abstract base class for constitutive material laws.
"""

import numpy as np
from abc import ABC, abstractmethod


class Material(ABC):
    r"""
    Abstract base class for all constitutive material laws.

    Every material must expose:

    - ``stress(eps)`` — scalar stress evaluation.
    - ``stress_array(eps)`` — vectorized evaluation over an array.
    - ``eps_min`` / ``eps_max`` — admissible strain range.

    The strain limits are consumed by the N-M diagram generator to
    automatically determine the scan range for ultimate configurations.

    Attributes
    ----------
    eps_min : float
        Most compressive admissible strain (typically negative).
    eps_max : float
        Most tensile admissible strain (typically positive).
    """

    @property
    @abstractmethod
    def eps_min(self):
        """Most compressive admissible strain."""
        ...

    @property
    @abstractmethod
    def eps_max(self):
        """Most tensile admissible strain."""
        ...

    @abstractmethod
    def stress(self, eps):
        """
        Compute stress for a single strain value.

        Parameters
        ----------
        eps : float
            Strain.

        Returns
        -------
        float
            Stress [MPa].
        """
        ...

    def stress_array(self, eps):
        """
        Vectorized stress computation.

        Default implementation loops over :meth:`stress`. Subclasses
        should override for performance.

        Parameters
        ----------
        eps : numpy.ndarray

        Returns
        -------
        numpy.ndarray
        """
        return np.array([self.stress(e) for e in eps])
