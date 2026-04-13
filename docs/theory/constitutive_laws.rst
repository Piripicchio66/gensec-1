.. _constitutive_laws:

==================
Constitutive laws
==================

GenSec is material-agnostic: any stress-strain relationship that
implements the :class:`~gensec.materials.Material` interface can be used.
Three built-in laws and a generic tabulated law are provided.


Abstract interface
-------------------

Every material must implement:

- ``stress(eps)`` — evaluate :math:`\sigma(\varepsilon)` for a scalar
  strain.
- ``stress_array(eps)`` — vectorized evaluation over a NumPy array.
- ``eps_min`` / ``eps_max`` — admissible strain range, consumed by the
  N-M diagram generator to determine scan bounds.


Concrete — parabola-rectangle (EC2 3.1.7)
-------------------------------------------

The parabola-rectangle law per EN 1992-1-1, §3.1.7:

.. math::

   \sigma_c(\varepsilon) =
   \begin{cases}
       0
           & \varepsilon > 0
           \\[4pt]
       -f_{cd}\!\left[1 - \left(1 - \dfrac{\varepsilon}
           {\varepsilon_{c2}}\right)^{\!n}\right]
           & \varepsilon_{c2} \le \varepsilon \le 0
           \\[4pt]
       -f_{cd}
           & \varepsilon_{cu2} \le \varepsilon < \varepsilon_{c2}
           \\[4pt]
       0
           & \varepsilon < \varepsilon_{cu2}
   \end{cases}

where:

- :math:`f_{cd} = \alpha_{cc}\,f_{ck}/\gamma_c` is the design
  compressive strength,
- :math:`\varepsilon_{c2}` is the strain at peak stress
  (default :math:`-0.002`),
- :math:`\varepsilon_{cu2}` is the ultimate compressive strain
  (default :math:`-0.0035`),
- :math:`n` is the parabolic exponent (default 2.0).

Concrete carries **no tensile stress** (:math:`\sigma = 0` for
:math:`\varepsilon > 0`).

The stress is negative (compressive) by sign convention.  The strain
range is :math:`[\varepsilon_{cu2},\, 0]`.

For high-strength concrete (:math:`f_{ck} > 50\;\text{MPa}`), the
parameters :math:`\varepsilon_{c2}`, :math:`\varepsilon_{cu2}`, and
:math:`n` deviate from the standard values and must be taken from
EC2 Table 3.1.  The factory function
:func:`~gensec.materials.concrete_from_ec2` computes them automatically.


Reinforcing steel — elastic-plastic with hardening
----------------------------------------------------

.. math::

   \sigma_s(\varepsilon) =
   \begin{cases}
       E_s\,\varepsilon
           & |\varepsilon| \le \varepsilon_{yd}
           \\[4pt]
       \mathrm{sign}(\varepsilon)\!\left[f_{yd}
           + (f_{td} - f_{yd})\,
           \dfrac{|\varepsilon| - \varepsilon_{yd}}
                 {\varepsilon_{su} - \varepsilon_{yd}}\right]
           & \varepsilon_{yd} < |\varepsilon| \le \varepsilon_{su}
           \\[4pt]
       0
           & |\varepsilon| > \varepsilon_{su}
   \end{cases}

where:

- :math:`f_{yd} = f_{yk}/\gamma_s` is the design yield strength,
- :math:`f_{td} = k \cdot f_{yd}` with :math:`k =` ``k_hardening``
  (ratio :math:`f_t/f_y`; 1.0 = perfectly plastic),
- :math:`\varepsilon_{yd} = f_{yd}/E_s` is the yield strain,
- :math:`\varepsilon_{su}` is the ultimate strain.

The law is **symmetric in tension and compression** by default.  Set
``works_in_compression = False`` to suppress compressive stress
(useful for modelling FRP or tendons that buckle in compression).


Structural steel (EN 10025-2)
------------------------------

Uses the same elastic-plastic model as reinforcing steel, but yield
and ultimate strengths are computed from plate thickness according to
NF EN 10025-2 Table 7.  The factory function
:func:`~gensec.materials.steel_from_en10025` handles the lookup.

Supported grades: **S235**, **S275**, **S355**.

The hardening ratio :math:`k = f_u/f_y` is computed automatically
from the tabulated :math:`f_u` and thickness-dependent :math:`f_y`.


Tabulated material
-------------------

For materials not covered by the built-in laws (CFRP, timber, FRP bars,
shotcrete, etc.), :class:`~gensec.materials.TabulatedMaterial` accepts
an arbitrary set of :math:`(\varepsilon_i,\,\sigma_i)` data points
and interpolates linearly between them.

Strains outside the table range produce :math:`\sigma = 0` (material
failure / rupture).

.. code-block:: python

   from gensec.materials import TabulatedMaterial

   cfrp = TabulatedMaterial(
       strains=[0.0, 0.017],
       stresses=[0.0, 2800.0],
       name="CFRP",
   )

The tabulated law is fully compatible with all solvers and diagram
generators.  The strain limits ``eps_min`` / ``eps_max`` are inferred
from the first and last data points.
