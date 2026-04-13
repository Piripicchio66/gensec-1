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
Point fiber definitions (rebars, FRP strips, tendons).

Phase 2: each fiber has both x and y coordinates for biaxial bending.
"""

from dataclasses import dataclass
from typing import Optional
from ..materials.base import Material


@dataclass
class RebarLayer:
    r"""
    A point fiber (rebar, FRP strip, tendon, etc.).

    Each fiber carries its own :class:`Material` reference, allowing
    mixed-material sections. For biaxial bending, both ``x`` and ``y``
    coordinates are needed.

    Parameters
    ----------
    y : float
        Vertical coordinate from bottom edge [mm].
    As : float
        Cross-sectional area [mm^2].
    material : Material
        Constitutive law.
    x : float, optional
        Horizontal coordinate from left edge [mm]. If ``None``,
        defaults to the section centroid x-coordinate (set during
        section assembly). For uniaxial bending this is irrelevant.
    embedded : bool, optional
        If ``True`` (default), the fiber is embedded within the bulk
        material. The integrator will subtract the bulk material
        contribution at this location to avoid double-counting the
        area. Set to ``False`` for external elements (e.g. external
        FRP strips, steel truss chords outside the concrete).
    n_bars : int, optional
        Number of bars (informational). Default 1.
    diameter : float, optional
        Bar diameter [mm] (informational). Default 0.
    """

    y: float
    As: float
    material: Material
    x: Optional[float] = None
    embedded: bool = True
    n_bars: int = 1
    diameter: float = 0.0
