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
YAML input loader for GenSec.

Reads a YAML file describing materials, section geometry, and
(optionally) load demands, and returns fully constructed GenSec
objects ready for analysis.

Section geometry format
-----------------------
The ``section`` block supports two modes:

**Legacy rectangular** (backward-compatible):

.. code-block:: yaml

    section:
      B: 300
      H: 600
      bulk_material: concrete_1
      n_fibers_y: 100
      n_fibers_x: 1
      rebars:
        - y: 40
          As: 942.5
          material: steel_1

**Generic section** (new):

.. code-block:: yaml

    section:
      shape: tee            # or: rect, circle, annulus, h, box,
                             #     single_tee, double_tee, custom
      params:
        bf: 800
        hf: 150
        bw: 300
        hw: 450
      bulk_material: concrete_1
      mesh_size: 15
      mesh_method: grid      # or: triangle
      rebars:
        - y: 40
          x: 150
          As: 942.5
          material: steel_1

**Custom polygon** (arbitrary vertex list):

.. code-block:: yaml

    section:
      shape: custom
      params:
        exterior: [[0,0], [300,0], [300,600], [0,600]]
        holes:
          - [[50,50], [250,50], [250,150], [50,150]]
      bulk_material: concrete_1
      mesh_size: 10
      mesh_method: triangle
      rebars: []

The YAML parser detects which mode to use:

- If ``shape`` is present → generic section.
- If ``B`` and ``H`` are present without ``shape`` → legacy
  rectangular (wrapped via :class:`RectSection`).
"""

import yaml
import numpy as np

from .materials import Concrete, Steel, TabulatedMaterial
from .geometry.fiber import RebarLayer
from .geometry.section import RectSection
from .geometry.geometry import GenericSection
from .geometry import primitives as prim
from .materials.ec2_bridge import concrete_from_class, concrete_from_ec2


# ---- Material builders (unchanged) ----

_MATERIAL_BUILDERS = {
    "concrete": {
        "cls": Concrete,
        "params": ["fck", "gamma_c", "alpha_cc", "n_parabola",
                    "eps_c2", "eps_cu2"],
    },
    "steel": {
        "cls": Steel,
        "params": ["fyk", "gamma_s", "Es", "k_hardening", "eps_su",
                    "works_in_compression"],
    },
    "tabulated": {
        "cls": TabulatedMaterial,
        "params": ["strains", "stresses", "name"],
    },
}


def _build_material(name, spec):
    """
    Build a Material instance from a YAML specification dict.

    Supported types: ``concrete``, ``concrete_ec2``, ``steel``,
    ``tabulated``.

    Parameters
    ----------
    name : str
        Key used in the YAML ``materials`` block.
    spec : dict
        Must contain a ``'type'`` key.

    Returns
    -------
    Material

    Raises
    ------
    ValueError
        Unknown material type.
    """
    mat_type = spec.get("type", "").lower()

    if mat_type == "concrete_ec2":
        conc_class = spec.get("class")
        if conc_class:
            return concrete_from_class(
                conc_class,
                ls=spec.get("ls", "F"),
                loadtype=spec.get("loadtype", "slow"),
                TypeConc=spec.get("TypeConc", "R"),
                NA=spec.get("NA", "French"),
                time=spec.get("time", 28),
            )
        fck = spec.get("fck")
        if fck is None:
            raise ValueError(
                f"Material '{name}': concrete_ec2 requires 'class' "
                f"(e.g. 'C30/37') or 'fck'."
            )
        return concrete_from_ec2(
            fck=float(fck),
            ls=spec.get("ls", "F"),
            loadtype=spec.get("loadtype", "slow"),
            TypeConc=spec.get("TypeConc", "R"),
            NA=spec.get("NA", "French"),
            time=spec.get("time", 28),
        )

    if mat_type not in _MATERIAL_BUILDERS:
        raise ValueError(
            f"Unknown material type '{mat_type}' for '{name}'. "
            f"Valid: {list(_MATERIAL_BUILDERS.keys())} + 'concrete_ec2'"
        )

    builder = _MATERIAL_BUILDERS[mat_type]
    cls = builder["cls"]
    kwargs = {}
    for p in builder["params"]:
        if p in spec:
            val = spec[p]
            if isinstance(val, list):
                val = np.array(val, dtype=float)
            kwargs[p] = val
    return cls(**kwargs)


# ---- Shape factory dispatch ----

_SHAPE_FACTORIES = {
    "rect": lambda p: prim.rect_poly(p["B"], p["H"]),
    "circle": lambda p: prim.circle_poly(
        p["D"], resolution=p.get("resolution", 64)),
    "annulus": lambda p: prim.annulus_poly(
        p["D_ext"], p["D_int"],
        resolution=p.get("resolution", 64)),
    "tee": lambda p: prim.tee_poly(
        p["bf"], p["hf"], p["bw"], p["hw"]),
    "inv_tee": lambda p: prim.inv_tee_poly(
        p["bf"], p["hf"], p["bw"], p["hw"]),
    "h": lambda p: prim.h_poly(
        p["bf"], p["hf_top"], p["hf_bot"], p["bw"], p["hw"]),
    "box": lambda p: prim.box_poly(
        p["B"], p["H"], p["tw"], p["tf_top"],
        tf_bot=p.get("tf_bot")),
    "single_tee": lambda p: prim.single_tee_slab_poly(
        p["b_top"], p["h_top"], p["bw"], p["hw"]),
    "double_tee": lambda p: prim.double_tee_slab_poly(
        p["b_top"], p["h_top"], p["bw"], p["hw"],
        p["stem_spacing"]),
    "custom": lambda p: prim.custom_poly(
        p["exterior"], holes=p.get("holes")),
}


def _build_polygon(sec_spec):
    r"""
    Build a Shapely polygon from the ``section`` YAML block.

    Parameters
    ----------
    sec_spec : dict
        The ``section`` block from YAML.

    Returns
    -------
    shapely.geometry.Polygon

    Raises
    ------
    ValueError
        If the shape type is not recognized.
    """
    shape = sec_spec["shape"].lower()
    params = sec_spec.get("params", {})

    if shape not in _SHAPE_FACTORIES:
        raise ValueError(
            f"Unknown section shape '{shape}'. "
            f"Valid: {list(_SHAPE_FACTORIES.keys())}"
        )

    return _SHAPE_FACTORIES[shape](params)


# ---- Main loader ----

def load_yaml(filepath):
    r"""
    Load a GenSec input file and return constructed objects.

    Detects whether the section block uses the legacy rectangular
    format (``B`` + ``H``) or the new generic format (``shape``).

    Parameters
    ----------
    filepath : str or pathlib.Path

    Returns
    -------
    dict
        Keys: ``'materials'``, ``'section'`` (GenericSection or
        RectSection), ``'demands'``, ``'combinations'``,
        ``'output_options'``.
    """
    with open(filepath, 'r') as f:
        data = yaml.safe_load(f)

    # ---- Materials ----
    materials = {}
    for mat_name, mat_spec in data.get("materials", {}).items():
        materials[mat_name] = _build_material(mat_name, mat_spec)

    # ---- Section ----
    sec_spec = data["section"]
    bulk_name = sec_spec["bulk_material"]
    if bulk_name not in materials:
        raise ValueError(
            f"Bulk material '{bulk_name}' not found in materials."
        )

    # Parse rebars (common to both modes)
    rebars = _parse_rebars(sec_spec, materials)

    if "shape" in sec_spec:
        # ---- New generic mode ----
        polygon = _build_polygon(sec_spec)

        # Optional multi-material zones
        bulk_materials = []
        for zone_spec in sec_spec.get("material_zones", []):
            zone_poly = _SHAPE_FACTORIES[
                zone_spec["shape"].lower()](zone_spec.get("params", {}))
            zone_mat_name = zone_spec["material"]
            if zone_mat_name not in materials:
                raise ValueError(
                    f"Zone material '{zone_mat_name}' not found."
                )
            bulk_materials.append((zone_poly, materials[zone_mat_name]))

        section = GenericSection(
            polygon=polygon,
            bulk_material=materials[bulk_name],
            rebars=rebars,
            mesh_size=float(sec_spec.get("mesh_size", 10)),
            mesh_method=sec_spec.get("mesh_method", "grid"),
            bulk_materials=bulk_materials,
        )
    else:
        # ---- Legacy rectangular mode ----
        section = RectSection(
            B=float(sec_spec["B"]),
            H=float(sec_spec["H"]),
            bulk_material=materials[bulk_name],
            rebars=rebars,
            n_fibers_y=int(sec_spec.get("n_fibers_y",
                            sec_spec.get("n_fibers", 100))),
            n_fibers_x=int(sec_spec.get("n_fibers_x", 1)),
        )

    # ---- Demands ----
    demands = [_parse_demand(d) for d in data.get("demands", [])]

    # ---- Combinations ----
    combinations = []
    for c_spec in data.get("combinations", []):
        combinations.append({
            "name": c_spec.get("name", "unnamed"),
            "demands": [_parse_demand(d)
                        for d in c_spec.get("demands", [])],
        })

    # ---- Output options ----
    output_opts = data.get("output", {})

    return {
        "materials": materials,
        "section": section,
        "demands": demands,
        "combinations": combinations,
        "output_options": output_opts,
    }


def _parse_rebars(sec_spec, materials):
    """
    Parse the ``rebars`` list from a section YAML block.

    Parameters
    ----------
    sec_spec : dict
        Section specification dict.
    materials : dict
        Material name → Material mapping.

    Returns
    -------
    list of RebarLayer
    """
    rebars = []
    for rb_spec in sec_spec.get("rebars", []):
        mat_name = rb_spec["material"]
        if mat_name not in materials:
            raise ValueError(
                f"Rebar material '{mat_name}' not found in materials."
            )
        rebars.append(RebarLayer(
            y=float(rb_spec["y"]),
            As=float(rb_spec["As"]),
            material=materials[mat_name],
            x=float(rb_spec["x"]) if "x" in rb_spec else None,
            embedded=bool(rb_spec.get("embedded", True)),
            n_bars=int(rb_spec.get("n_bars", 1)),
            diameter=float(rb_spec.get("diameter", 0)),
        ))
    return rebars


def _parse_demand(d_spec):
    """
    Parse a single demand triple from YAML.

    Accepts ``Mx_kNm`` / ``My_kNm`` (canonical) or legacy
    ``M_kNm`` (Mx only, My=0).

    Parameters
    ----------
    d_spec : dict

    Returns
    -------
    dict
        Keys: ``name``, ``N`` [N], ``Mx`` [N*mm], ``My`` [N*mm].
    """
    N = float(d_spec.get("N_kN", 0)) * 1e3

    if "Mx_kNm" in d_spec:
        Mx = float(d_spec["Mx_kNm"]) * 1e6
        My = float(d_spec.get("My_kNm", 0)) * 1e6
    elif "M_kNm" in d_spec:
        Mx = float(d_spec["M_kNm"]) * 1e6
        My = 0.0
    else:
        Mx = 0.0
        My = 0.0

    return {
        "name": d_spec.get("name", "unnamed"),
        "N": N,
        "Mx": Mx,
        "My": My,
    }
