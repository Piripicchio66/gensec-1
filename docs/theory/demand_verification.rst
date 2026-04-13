.. _demand_verification:

=====================
Demand verification
=====================

Once the resistance domain is computed, GenSec verifies load demands
by measuring how close each demand point is to the domain boundary.
The :class:`~gensec.solver.DemandChecker` handles both 2-D
(:math:`N\text{-}M`) and 3-D (:math:`N\text{-}M_x\text{-}M_y`)
domains.


Utilization ratio
-----------------

The **utilization ratio** :math:`\eta` is a scalar measure of how
"used up" the section capacity is by a given demand.  GenSec computes
it via ray-casting from a reference point inside the domain through the
demand point to the convex hull boundary:

.. math::

   \eta =
       \frac{|\mathbf{d} - \mathbf{d}_{\text{ref}}|}
            {|\mathbf{r} - \mathbf{d}_{\text{ref}}|}

where:

- :math:`\mathbf{d}` is the demand point
  :math:`(N,\, M_x)` or :math:`(N,\, M_x,\, M_y)`,
- :math:`\mathbf{d}_{\text{ref}}` is a reference point inside the
  domain (typically the centroid of the convex hull),
- :math:`\mathbf{r}` is the intersection of the ray
  :math:`\mathbf{d}_{\text{ref}} \to \mathbf{d}` with the domain
  boundary.

Interpretation:

+--------------------+---------------------------------------------+
| :math:`\eta < 1`  | Demand inside the domain — **verified**.     |
+--------------------+---------------------------------------------+
| :math:`\eta = 1`  | Demand exactly on the boundary.              |
+--------------------+---------------------------------------------+
| :math:`\eta > 1`  | Demand outside the domain — **not verified**.|
+--------------------+---------------------------------------------+


Convex hull
-----------

The domain boundary is represented as a
:class:`scipy.spatial.ConvexHull`.  This is exact for convex domains
(which is the typical case for RC sections at ULS under the
parabola-rectangle law) and conservative for non-convex domains
(the convex hull encloses the true domain, so :math:`\eta` may be
slightly underestimated).

In 2-D mode, the hull is built on the :math:`(N, M_x)` point cloud.
In 3-D mode, on :math:`(N, M_x, M_y)`.  The ray-casting algorithm
solves for the intersection parameter along the half-line from the
reference point through the demand, testing against all hull facets.


Combinations
------------

When demands are grouped into **combinations** (see
:doc:`/user_guide/yaml_reference`), GenSec computes :math:`\eta` for
every triple in the combination and reports:

- :math:`\eta_i` for each individual demand,
- :math:`\eta_{\max} = \max_i \eta_i` for the combination as a whole.

This is useful for envelope verification of multiple columns sharing
the same section.


Per-fiber post-processing
--------------------------

For any demand point, GenSec can solve the inverse problem (find the
strain plane that equilibrates the demand) and then extract the
stress/strain state at every fiber and rebar.  This is available
through:

- :meth:`FiberSolver.solve_equilibrium` — inverse solver,
- :meth:`FiberSolver.get_fiber_results` — stress/strain extraction.

The results include bulk fiber strains, bulk fiber stresses, rebar
strains, rebar stresses, and lever arms.  These are exported to CSV
and plotted on the section geometry by the CLI.
