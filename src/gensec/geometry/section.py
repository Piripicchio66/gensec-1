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
Section assembly — backward-compatible wrappers.

The :class:`RectSection` class is retained for backward compatibility
and YAML files that specify simple rectangular sections. Internally,
it creates a :class:`GenericSection` with a rectangular polygon.

For new code, use :class:`GenericSection` directly with any
Shapely polygon, or use the factory functions in
:mod:`gensec.geometry.primitives`.
"""

import numpy as np
from dataclasses import dataclass
from typing import List

from .fiber import RebarLayer
from .geometry import GenericSection
from .primitives import rect_poly
from ..materials.base import Material


@dataclass
class RectSection:
    r"""
    Rectangular cross-section — backward-compatible wrapper.

    Constructs a :class:`GenericSection` with a rectangular polygon
    of size :math:`B \times H`. All attributes expected by
    :class:`~gensec.solver.FiberSolver` are delegated to the inner
    ``GenericSection``.

    For new code, prefer using :class:`GenericSection` directly.

    Parameters
    ----------
    B : float
        Width (x-direction) [mm].
    H : float
        Height (y-direction) [mm].
    bulk_material : Material
        Bulk material (concrete, timber, ...).
    rebars : list of RebarLayer
        Point fibers.
    n_fibers_y : int, optional
        Fiber rows. Default 100.
    n_fibers_x : int, optional
        Fiber columns. Default 1 (uniaxial mode).

    Attributes
    ----------
    _generic : GenericSection
        Inner section object that does all the work.
    """

    B: float
    H: float
    bulk_material: Material
    rebars: List[RebarLayer]
    n_fibers_y: int = 100
    n_fibers_x: int = 1

    def __post_init__(self):
        # Compute mesh_size from the requested grid resolution.
        # The grid mesher in GenericSection uses ceil(B/s) and
        # ceil(H/s), so we pick s to reproduce the exact nx, ny.
        dx = self.B / max(self.n_fibers_x, 1)
        dy = self.H / max(self.n_fibers_y, 1)
        # Use the smaller of the two so neither axis is coarser
        # than requested.
        mesh_size = min(dx, dy)

        poly = rect_poly(self.B, self.H)
        self._generic = GenericSection(
            polygon=poly,
            bulk_material=self.bulk_material,
            rebars=self.rebars,
            mesh_size=mesh_size,
            mesh_method="grid",
        )

        # Expose grid dimensions (may differ slightly due to ceil)
        self.n_fibers_x = self._generic.n_fibers_x
        self.n_fibers_y = self._generic.n_fibers_y
        self.n_fibers = self._generic.n_fibers
        self.dx = self.B / self.n_fibers_x
        self.dy = self.H / self.n_fibers_y

    # ------------------------------------------------------------------
    #  Delegated array attributes
    # ------------------------------------------------------------------

    @property
    def x_fibers(self):
        """x-coordinates of bulk fiber centroids [mm]."""
        return self._generic.x_fibers

    @property
    def y_fibers(self):
        """y-coordinates of bulk fiber centroids [mm]."""
        return self._generic.y_fibers

    @property
    def A_fibers(self):
        """Bulk fiber areas [mm²]."""
        return self._generic.A_fibers

    @property
    def x_rebars(self):
        """x-coordinates of rebar layers [mm]."""
        return self._generic.x_rebars

    @property
    def y_rebars(self):
        """y-coordinates of rebar layers [mm]."""
        return self._generic.y_rebars

    @property
    def A_rebars(self):
        """Rebar areas [mm²]."""
        return self._generic.A_rebars

    @property
    def embedded_rebars(self):
        """Embedded flags for rebars."""
        return self._generic.embedded_rebars

    @property
    def x_centroid(self):
        r"""Gross centroid x-coordinate: :math:`B/2` [mm]."""
        return self._generic.x_centroid

    @property
    def y_centroid(self):
        r"""Gross centroid y-coordinate: :math:`H/2` [mm]."""
        return self._generic.y_centroid

    @property
    def polygon(self):
        """Underlying Shapely polygon."""
        return self._generic.polygon
