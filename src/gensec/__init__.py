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
GenSec - Generic Section Calculator
=====================================

Modular Python package for fiber-based cross-section analysis of
composite structural members under combined axial force and bending.

Subpackages
-----------
materials
    Constitutive laws (abstract base, concrete, steel, tabulated).
geometry
    Section definition (patches, point fibers, assemblies).
solver
    Strain integration, equilibrium solver, N-M diagram generation.
output
    Reporting and plotting utilities.
"""

#__version__ = "0.6.0"

# subfolders_tool/__init__.py
from importlib.metadata import version as _pkg_version, PackageNotFoundError

try:
    __version__ = _pkg_version("gensec")
except PackageNotFoundError:
    # Fallback per PyInstaller o esecuzione da sorgente
    try:
        from ._version import __version__
    except ImportError:
        __version__ = "X.X.X"