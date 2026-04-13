"""
Unit tests for material constitutive laws.

Covers: Concrete, Steel, TabulatedMaterial, EC2 bridge, EN 10025 bridge.
"""

import unittest
import numpy as np
from gensec.materials import (
    Concrete, Steel, TabulatedMaterial,
    concrete_from_ec2, concrete_from_class, steel_from_en10025,
    fben2,
)


class TestConcreteStressStrain(unittest.TestCase):
    """Pointwise verification of the parabola-rectangle law."""

    def setUp(self):
        self.c = Concrete(fck=25.0, gamma_c=1.5, alpha_cc=0.85)
        self.fcd = 0.85 * 25.0 / 1.5  # 14.1667

    def test_zero_strain(self):
        self.assertAlmostEqual(self.c.stress(0.0), 0.0)

    def test_tension_zero(self):
        self.assertAlmostEqual(self.c.stress(0.001), 0.0)
        self.assertAlmostEqual(self.c.stress(0.01), 0.0)

    def test_peak_stress(self):
        """At eps_c2, sigma must equal -fcd."""
        self.assertAlmostEqual(self.c.stress(-0.002), -self.fcd, places=3)

    def test_plateau(self):
        """Between eps_c2 and eps_cu2, sigma = -fcd (constant)."""
        self.assertAlmostEqual(self.c.stress(-0.003), -self.fcd, places=3)
        self.assertAlmostEqual(self.c.stress(-0.0035), -self.fcd, places=3)

    def test_crushing(self):
        """Beyond eps_cu2, sigma = 0."""
        self.assertAlmostEqual(self.c.stress(-0.004), 0.0)
        self.assertAlmostEqual(self.c.stress(-0.01), 0.0)

    def test_parabolic_midpoint(self):
        r"""At eps = eps_c2/2 = -0.001, eta=0.5,
        sigma = -fcd * (1 - (1-0.5)^2) = -fcd * 0.75."""
        expected = -self.fcd * 0.75
        self.assertAlmostEqual(self.c.stress(-0.001), expected, places=3)

    def test_parabolic_quarter(self):
        r"""At eps = eps_c2/4 = -0.0005, eta=0.25,
        sigma = -fcd * (1 - 0.75^2) = -fcd * 0.4375."""
        expected = -self.fcd * (1 - 0.75**2)
        self.assertAlmostEqual(self.c.stress(-0.0005), expected, places=3)

    def test_stress_array_matches_scalar(self):
        eps = np.array([0.001, 0.0, -0.0005, -0.001, -0.002,
                        -0.003, -0.0035, -0.004])
        arr = self.c.stress_array(eps)
        for i, e in enumerate(eps):
            self.assertAlmostEqual(arr[i], self.c.stress(e), places=8)

    def test_eps_limits(self):
        self.assertEqual(self.c.eps_min, -0.0035)
        self.assertEqual(self.c.eps_max, 0.0)


class TestSteelStressStrain(unittest.TestCase):
    """Pointwise verification of elastic-plastic steel law."""

    def setUp(self):
        self.s = Steel(fyk=450.0, gamma_s=1.15, Es=200000.0,
                       k_hardening=1.0, works_in_compression=True)
        self.fyd = 450.0 / 1.15
        self.eps_yd = self.fyd / 200000.0

    def test_zero(self):
        self.assertAlmostEqual(self.s.stress(0.0), 0.0)

    def test_elastic_tension(self):
        eps = self.eps_yd / 2
        self.assertAlmostEqual(self.s.stress(eps), 200000 * eps, places=2)

    def test_at_yield(self):
        self.assertAlmostEqual(self.s.stress(self.eps_yd), self.fyd,
                               places=2)

    def test_plastic_tension(self):
        """Perfectly plastic: sig = fyd for eps_yd < eps <= eps_su."""
        self.assertAlmostEqual(self.s.stress(0.005), self.fyd, places=2)
        self.assertAlmostEqual(self.s.stress(0.01), self.fyd, places=2)

    def test_rupture(self):
        self.assertAlmostEqual(self.s.stress(0.02), 0.0)

    def test_elastic_compression(self):
        eps = -self.eps_yd / 2
        self.assertAlmostEqual(self.s.stress(eps), 200000 * eps, places=2)

    def test_plastic_compression(self):
        self.assertAlmostEqual(self.s.stress(-0.005), -self.fyd, places=2)

    def test_compression_rupture(self):
        self.assertAlmostEqual(self.s.stress(-0.02), 0.0)

    def test_tension_only(self):
        s_t = Steel(fyk=450.0, works_in_compression=False)
        self.assertAlmostEqual(s_t.stress(-0.005), 0.0)
        self.assertAlmostEqual(s_t.stress(-self.eps_yd), 0.0)
        self.assertGreater(s_t.stress(self.eps_yd), 0.0)

    def test_hardening(self):
        s_h = Steel(fyk=450.0, gamma_s=1.15, k_hardening=1.15,
                    eps_su=0.05)
        fyd = 450.0 / 1.15
        ftd = fyd * 1.15
        # At eps_su, should be at ftd
        self.assertAlmostEqual(s_h.stress(0.05), ftd, places=1)
        # Midway between yield and ultimate
        eps_mid = (s_h.eps_yd + 0.05) / 2
        sig_mid = fyd + (ftd - fyd) * 0.5
        self.assertAlmostEqual(s_h.stress(eps_mid), sig_mid, places=1)

    def test_eps_limits(self):
        self.assertEqual(self.s.eps_min, -0.01)
        self.assertEqual(self.s.eps_max, 0.01)
        s_t = Steel(fyk=450.0, works_in_compression=False)
        self.assertEqual(s_t.eps_min, 0.0)

    def test_stress_array_matches_scalar(self):
        eps = np.linspace(-0.015, 0.015, 100)
        arr = self.s.stress_array(eps)
        for i, e in enumerate(eps):
            self.assertAlmostEqual(arr[i], self.s.stress(e), places=8)


class TestTabulatedMaterial(unittest.TestCase):
    """Tests for user-defined tabulated material."""

    def setUp(self):
        self.mat = TabulatedMaterial(
            [-0.01, -0.002, 0.0, 0.002, 0.01],
            [-400.0, -400.0, 0.0, 400.0, 400.0],
            name="EPP_test",
        )

    def test_interpolation(self):
        self.assertAlmostEqual(self.mat.stress(0.001), 200.0)
        self.assertAlmostEqual(self.mat.stress(-0.001), -200.0)

    def test_plateau(self):
        self.assertAlmostEqual(self.mat.stress(0.005), 400.0)
        self.assertAlmostEqual(self.mat.stress(-0.005), -400.0)

    def test_outside_range(self):
        self.assertAlmostEqual(self.mat.stress(0.02), 0.0)
        self.assertAlmostEqual(self.mat.stress(-0.02), 0.0)

    def test_limits(self):
        self.assertEqual(self.mat.eps_min, -0.01)
        self.assertEqual(self.mat.eps_max, 0.01)

    def test_cfrp(self):
        frp = TabulatedMaterial([0.0, 0.017], [0.0, 2800.0], "CFRP")
        self.assertAlmostEqual(frp.stress(0.0085), 1400.0, places=0)
        self.assertEqual(frp.stress(-0.001), 0.0)
        self.assertEqual(frp.stress(0.02), 0.0)

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            TabulatedMaterial([1, 0], [0, 1])  # not increasing
        with self.assertRaises(ValueError):
            TabulatedMaterial([0], [0])  # too few points

    def test_vectorized(self):
        eps = np.array([-0.02, -0.005, 0.0, 0.001, 0.005, 0.02])
        arr = self.mat.stress_array(eps)
        for i, e in enumerate(eps):
            self.assertAlmostEqual(arr[i], self.mat.stress(e), places=8)


class TestEC2Bridge(unittest.TestCase):
    """Tests for concrete_from_ec2 and concrete_from_class."""

    def test_c25_standard(self):
        c = concrete_from_ec2(fck=25, ls='F')
        self.assertAlmostEqual(c.fcd, 0.85 * 25 / 1.5, places=2)
        self.assertEqual(c.eps_cu2, -0.0035)
        self.assertEqual(c.n_parabola, 2.0)

    def test_c30_french_annex(self):
        c = concrete_from_ec2(fck=30, ls='F', loadtype='slow',
                              TypeConc='R', NA='French')
        self.assertAlmostEqual(c.fcd, 0.85 * 30 / 1.5, places=2)
        self.assertTrue(hasattr(c, 'ec2'))
        self.assertAlmostEqual(c.ec2.fcm, 38.0)

    def test_c70_high_strength(self):
        """fck > 50 must give different n, eps_c2, eps_cu2."""
        c = concrete_from_ec2(fck=70, ls='F')
        self.assertNotEqual(c.n_parabola, 2.0)
        self.assertNotEqual(c.eps_c2, -0.002)
        self.assertNotEqual(c.eps_cu2, -0.0035)
        # n should be ~1.44, eps_c2 ~-0.00242, eps_cu2 ~-0.00266
        self.assertAlmostEqual(c.n_parabola, 1.4 + 23.4*((90-70)/100)**4,
                               places=2)

    def test_from_class(self):
        c = concrete_from_class('C50/60', ls='F')
        self.assertEqual(c.fck, 50.0)

    def test_accidental_ls(self):
        c = concrete_from_ec2(fck=25, ls='A')
        self.assertAlmostEqual(c.ec2.gamma_c, 1.2)

    def test_unknown_class(self):
        with self.assertRaises(ValueError):
            concrete_from_class('C999')


class TestEN10025Bridge(unittest.TestCase):
    """Tests for steel_from_en10025."""

    def test_s355_thin(self):
        s = steel_from_en10025('S355', t=10)
        self.assertEqual(s.fyk, 355)

    def test_s355_thick(self):
        s = steel_from_en10025('S355', t=20)
        self.assertEqual(s.fyk, 345)  # thickness-dependent

    def test_s235(self):
        s = steel_from_en10025('S235', t=50)
        self.assertEqual(s.fyk, 215)

    def test_hardening_ratio(self):
        s = steel_from_en10025('S355', t=10)
        self.assertAlmostEqual(s.k_hardening, 470 / 355, places=3)


if __name__ == '__main__':
    unittest.main()
