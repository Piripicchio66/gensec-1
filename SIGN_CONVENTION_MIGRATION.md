# Sign Convention Migration Guide
# ================================
#
# Summary of changes to the integrator (integrator.py):
#
# OLD convention:
#   strain:  eps = eps0 + chi_x * (y - yref) + chi_y * (x - xref)
#   forces:  My  = +sum(sigma * A * (x - xref))
#   meaning: My > 0 compressed LEFT edge (x=0)
#
# NEW convention (right-hand rule):
#   strain:  eps = eps0 + chi_x * (y - yref) - chi_y * (x - xref)
#   forces:  My  = -sum(sigma * A * (x - xref))
#   meaning: My > 0 compresses RIGHT edge (x=xmax)
#
# Mx is UNCHANGED:
#   strain:  eps = eps0 + chi_x * (y - yref) + ...
#   forces:  Mx  = +sum(sigma * A * (y - yref))
#   meaning: Mx > 0 compresses BOTTOM edge (y=0)
#
#
# How to update tests:
# ====================
#
# 1. Any test that calls solver.integrate(eps0, chi_x, chi_y)
#    and checks the returned My value:
#    -> The NEW My has OPPOSITE SIGN to the OLD My.
#    -> If old test expected My = +X, new test should expect My = -X.
#    -> N and Mx are unchanged.
#
# 2. Any test that calls solver.strain_field(eps0, chi_x, chi_y)
#    and checks strains at specific (x, y) locations:
#    -> The strain contribution from chi_y is now NEGATED.
#    -> OLD: eps_at(x,y) = eps0 + chi_x*(y-yref) + chi_y*(x-xref)
#    -> NEW: eps_at(x,y) = eps0 + chi_x*(y-yref) - chi_y*(x-xref)
#    -> For a fiber at x > xref with chi_y > 0:
#       OLD strain was MORE POSITIVE (tensile), NEW is MORE NEGATIVE (compressive).
#
# 3. Any test that calls solver.solve_equilibrium(N, Mx, My)
#    -> The My INPUT meaning has changed.
#    -> If old test imposed My = +X (meaning "compress left"),
#       the same physical configuration is now My = -X.
#    -> OR: keep the same My input value, but expect different
#       eps0, chi_x, chi_y output (chi_y will have opposite sign).
#
# 4. Neutral axis position tests:
#    -> For biaxial bending with chi_y != 0, the neutral axis
#       orientation is mirrored in x. If old test expected NA
#       at x=352.5, new test should check the recalculated position.
#
# 5. Tests that only use uniaxial bending (chi_y = 0):
#    -> UNCHANGED. N, Mx, eps are all the same.
#
# 6. The Jacobian:
#    -> J[:, 2] (partial derivatives w.r.t. chi_y) changes sign.
#    -> J[2, :] (partial derivatives of My) changes sign.
#    -> J[2, 2] is UNCHANGED (double negation).
#
#
# Quick fix pattern:
# ==================
#
# For tests that verify My output:
#   OLD: self.assertAlmostEqual(My, +6082123, delta=100)
#   NEW: self.assertAlmostEqual(My, -6082123, delta=100)
#
# For tests that verify strain at (x, y) with chi_y != 0:
#   OLD: expected_eps = eps0 + chi_x*(y-yref) + chi_y*(x-xref)
#   NEW: expected_eps = eps0 + chi_x*(y-yref) - chi_y*(x-xref)
#
# For tests that impose a physical configuration via My:
#   OLD: sol = solver.solve_equilibrium(N, Mx, My=+100e6)  # compress left
#   NEW: sol = solver.solve_equilibrium(N, Mx, My=-100e6)  # same physical config
