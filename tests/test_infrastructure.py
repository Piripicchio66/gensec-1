"""
Unit tests for YAML loader, export, and DemandChecker.
"""

import unittest
import json
import csv
import os
import tempfile
import numpy as np
from gensec.materials import Concrete, Steel
from gensec.geometry import RebarLayer, RectSection
from gensec.solver import FiberSolver, NMDiagram, DemandChecker
from gensec.io_yaml import load_yaml
from pathlib import Path

_EXAMPLE_YAML = str(Path(__file__).parent.parent / "examples" / "example_input.yaml")
from gensec.output import (
    export_nm_domain_csv, export_nm_domain_json,
    export_demand_results_csv, export_demand_results_json,
    export_fiber_results_csv,
)


class TestYAMLLoader(unittest.TestCase):

    def test_loads_section(self):
        data = load_yaml(_EXAMPLE_YAML)
        self.assertEqual(data["section"].B, 300)
        self.assertEqual(data["section"].H, 600)

    def test_loads_demands(self):
        data = load_yaml(_EXAMPLE_YAML)
        self.assertEqual(len(data["demands"]), 2)
        self.assertEqual(data["demands"][0]["name"], "Gravity")

    def test_loads_materials(self):
        data = load_yaml(_EXAMPLE_YAML)
        self.assertIn("concrete_1", data["materials"])
        self.assertIsInstance(data["materials"]["concrete_1"], Concrete)

    def test_rebars_count(self):
        data = load_yaml(_EXAMPLE_YAML)
        self.assertEqual(len(data["section"].rebars), 3)


class TestExport(unittest.TestCase):

    def setUp(self):
        concrete = Concrete(fck=25.0)
        steel = Steel(fyk=450.0)
        rb = [RebarLayer(y=40, As=942.5, material=steel),
              RebarLayer(y=560, As=942.5, material=steel)]
        sec = RectSection(B=300, H=600, bulk_material=concrete,
                          rebars=rb, n_fibers_y=100, n_fibers_x=1)
        self.sv = FiberSolver(sec)
        self.nm = NMDiagram(self.sv).generate(200)
        self.td = tempfile.mkdtemp()

    def test_nm_csv(self):
        p = os.path.join(self.td, "nm.csv")
        export_nm_domain_csv(self.nm, p)
        with open(p) as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)
        self.assertEqual(header, ["N_kN", "Mx_kNm"])
        self.assertEqual(len(rows), len(self.nm["N"]))

    def test_nm_json(self):
        p = os.path.join(self.td, "nm.json")
        export_nm_domain_json(self.nm, p)
        with open(p) as f:
            data = json.load(f)
        self.assertEqual(data["n_points"], len(self.nm["N"]))

    def test_fiber_csv(self):
        fr = self.sv.get_fiber_results(-0.001, 5e-6)
        p = os.path.join(self.td, "fibers.csv")
        export_fiber_results_csv(fr, p)
        self.assertGreater(os.path.getsize(p), 100)

    def test_demand_csv(self):
        results = [{"name": "A", "N_kN": -1500, "M_kNm": 200,
                     "inside": True, "utilization": 0.7, "verified": True}]
        p = os.path.join(self.td, "dem.csv")
        export_demand_results_csv(results, p)
        with open(p) as f:
            reader = csv.reader(f)
            next(reader)  # header
            rows = list(reader)
        self.assertEqual(len(rows), 1)

    def test_demand_json(self):
        results = [{"name": "A", "N_kN": -1500, "M_kNm": 200,
                     "inside": True, "utilization": 0.7, "verified": True}]
        p = os.path.join(self.td, "dem.json")
        export_demand_results_json(results, p)
        with open(p) as f:
            data = json.load(f)
        self.assertEqual(len(data["demands"]), 1)


class TestDemandChecker2D(unittest.TestCase):

    def setUp(self):
        concrete = Concrete(fck=25.0)
        steel = Steel(fyk=450.0)
        rb = [RebarLayer(y=40, As=942.5, material=steel),
              RebarLayer(y=560, As=942.5, material=steel)]
        sec = RectSection(B=300, H=600, bulk_material=concrete,
                          rebars=rb, n_fibers_y=200, n_fibers_x=1)
        nm = NMDiagram(FiberSolver(sec)).generate(500)
        self.checker = DemandChecker(nm)

    def test_inside(self):
        self.assertTrue(self.checker.is_inside(-500e3, 50e6))

    def test_outside(self):
        self.assertFalse(self.checker.is_inside(-5000e3, 0))

    def test_eta_inside(self):
        eta = self.checker.utilization_ratio(-1500e3, 200e6)
        self.assertGreater(eta, 0)
        self.assertLess(eta, 1.0)

    def test_eta_outside(self):
        eta = self.checker.utilization_ratio(-1500e3, 500e6)
        self.assertGreater(eta, 1.0)

    def test_batch(self):
        results = self.checker.check_demands([
            {"name": "A", "N": -1500e3, "M": 200e6},
            {"name": "B", "N": -1500e3, "M": 500e6},
        ])
        self.assertTrue(results[0]["verified"])
        self.assertFalse(results[1]["verified"])


if __name__ == '__main__':
    unittest.main()


class TestYAMLCombinations(unittest.TestCase):
    """Test YAML loading with combinations and explicit Mx/My."""

    def test_loads_combinations(self):
        data = load_yaml(_EXAMPLE_YAML)
        self.assertIn("combinations", data)
        combos = data["combinations"]
        self.assertEqual(len(combos), 2)
        self.assertEqual(combos[0]["name"], "Seismic_X_envelope")
        self.assertEqual(len(combos[0]["demands"]), 5)

    def test_demand_has_Mx_My(self):
        data = load_yaml(_EXAMPLE_YAML)
        for d in data["demands"]:
            self.assertIn("Mx", d)
            self.assertIn("My", d)

    def test_combination_demand_has_Mx_My(self):
        data = load_yaml(_EXAMPLE_YAML)
        for c in data["combinations"]:
            for d in c["demands"]:
                self.assertIn("Mx", d)
                self.assertIn("My", d)

    def test_combination_demand_My_nonzero(self):
        data = load_yaml(_EXAMPLE_YAML)
        combo = data["combinations"][0]
        has_my = any(abs(d["My"]) > 0 for d in combo["demands"])
        self.assertTrue(has_my)


class TestYAMLConcreteEC2(unittest.TestCase):
    """Test concrete_ec2 material type in YAML."""

    def test_loads_ec2_class(self):
        data = load_yaml(
            str(Path(__file__).parent.parent / "examples" / "biaxial_column.yaml"))
        c = data["materials"]["C30"]
        self.assertAlmostEqual(c.fcd, 0.85 * 30 / 1.5, places=2)
        self.assertEqual(c.eps_cu2, -0.0035)
        self.assertEqual(c.n_parabola, 2.0)
        self.assertTrue(hasattr(c, 'ec2'))
        self.assertAlmostEqual(c.ec2.fcm, 38.0)


class TestSectionPlot(unittest.TestCase):
    """Test section drawing functionality."""

    def test_geometry_only(self):
        from gensec.output import plot_section
        import matplotlib; matplotlib.use('Agg')
        concrete = Concrete(fck=25.0)
        steel = Steel(fyk=450.0)
        sec = RectSection(B=300, H=600, bulk_material=concrete,
                          rebars=[
                              RebarLayer(y=40, As=942, material=steel,
                                         diameter=20),
                              RebarLayer(y=560, As=942, material=steel,
                                         diameter=20),
                          ], n_fibers_y=50, n_fibers_x=1)
        fig = plot_section(sec)
        self.assertIsNotNone(fig)

    def test_with_stress_field(self):
        from gensec.output import plot_section
        import matplotlib; matplotlib.use('Agg')
        concrete = Concrete(fck=25.0)
        steel = Steel(fyk=450.0)
        A20 = np.pi / 4 * 20**2
        sec = RectSection(B=300, H=600, bulk_material=concrete,
                          rebars=[
                              RebarLayer(y=40, x=40, As=A20, material=steel,
                                         diameter=20),
                              RebarLayer(y=560, x=260, As=A20, material=steel,
                                         diameter=20),
                          ], n_fibers_y=30, n_fibers_x=15)
        sv = FiberSolver(sec)
        fr = sv.get_fiber_results(-0.001, 5e-6, 2e-6)
        fig = plot_section(sec, fr)
        self.assertIsNotNone(fig)
