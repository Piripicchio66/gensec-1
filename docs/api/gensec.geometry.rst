.. _api_geometry:

===========================
:mod:`gensec.geometry`
===========================

Section definition, fiber meshing, and parametric shape factories.

A cross-section in GenSec is represented as a collection of **bulk
fibers** (from the meshed polygon) and **point fibers** (rebars,
tendons, FRP strips).  The :class:`~gensec.geometry.geometry.GenericSection`
class is the primary section object; :class:`~gensec.geometry.section.RectSection`
is a backward-compatible wrapper.


Point fibers (rebars)
----------------------

.. automodule:: gensec.geometry.fiber
   :members:
   :show-inheritance:
   :no-index: 


Generic section with polygon meshing
--------------------------------------

.. automodule:: gensec.geometry.geometry
   :members:
   :show-inheritance:
   :no-index: 


Parametric section primitives
------------------------------

Factory functions that return :class:`shapely.geometry.Polygon` objects
for common section shapes.  Each polygon can be passed directly to
:class:`~gensec.geometry.geometry.GenericSection`.

.. automodule:: gensec.geometry.primitives
   :members:
   :show-inheritance:
   :no-index: 


Rectangular section (legacy wrapper)
--------------------------------------

.. automodule:: gensec.geometry.section
   :members:
   :show-inheritance:
   :no-index: 

..
   Package-level exports
   ----------------------

   .. automodule:: gensec.geometry
      :members:
      :show-inheritance:
      :no-index: 
