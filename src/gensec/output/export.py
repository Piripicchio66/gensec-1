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
Numerical export utilities — CSV and JSON.

Provides functions to export N-M domain points, demand verification
results, and per-fiber stress/strain states.
"""

import json
import csv
import numpy as np
import os


def export_nm_domain_csv(nm_data, filepath):
    """
    Export N-M domain point cloud to CSV.

    Columns: ``N_kN``, ``Mx_kNm``.

    Parameters
    ----------
    nm_data : dict
        Output of :meth:`NMDiagram.generate`.
    filepath : str
        Output CSV file path.
    """
    with open(filepath, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["N_kN", "Mx_kNm"])
        for n, m in zip(nm_data["N_kN"], nm_data["M_kNm"]):
            w.writerow([f"{n:.4f}", f"{m:.6f}"])


def export_nm_domain_json(nm_data, filepath):
    """
    Export N-M domain point cloud to JSON.

    Parameters
    ----------
    nm_data : dict
        Output of :meth:`NMDiagram.generate`.
    filepath : str
        Output JSON file path.
    """
    data = {
        "description": "N-Mx interaction domain point cloud",
        "units": {"N": "kN", "Mx": "kN*m"},
        "n_points": len(nm_data["N_kN"]),
        "N_kN": [round(float(x), 4) for x in nm_data["N_kN"]],
        "Mx_kNm": [round(float(x), 6) for x in nm_data["M_kNm"]],
    }
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)


def export_demand_results_csv(results, filepath):
    """
    Export demand verification summary to CSV.

    Parameters
    ----------
    results : list of dict
        Output of :meth:`DemandChecker.check_demands`.
    filepath : str
    """
    with open(filepath, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["name", "N_kN", "Mx_kNm", "My_kNm", "inside",
                     "utilization", "verified"])
        for r in results:
            w.writerow([
                r["name"],
                f"{r['N_kN']:.2f}",
                f"{r['M_kNm']:.4f}",
                f"{r.get('My_kNm', 0):.4f}",
                r["inside"],
                f"{r['utilization']:.4f}",
                r["verified"],
            ])


def export_demand_results_json(results, filepath):
    """
    Export demand verification summary to JSON.

    Parameters
    ----------
    results : list of dict
        Output of :meth:`DemandChecker.check_demands`.
    filepath : str
    """
    # Ensure serializable
    clean = []
    for r in results:
        clean.append({
            "name": r["name"],
            "N_kN": round(r["N_kN"], 2),
            "M_kNm": round(r["M_kNm"], 4),
            "inside": bool(r["inside"]),
            "utilization": round(r["utilization"], 4),
            "verified": bool(r["verified"]),
        })
    with open(filepath, 'w') as f:
        json.dump({"demands": clean}, f, indent=2)


def export_fiber_results_csv(fiber_results, filepath):
    """
    Export per-fiber strain/stress state to CSV.

    Two sections in the file: BULK and REBARS, separated by a blank line.

    Parameters
    ----------
    fiber_results : dict
        Output of :meth:`FiberSolver.get_fiber_results`.
    filepath : str
    """
    with open(filepath, 'w', newline='') as f:
        w = csv.writer(f)

        bulk = fiber_results["bulk"]
        has_x = "x" in bulk

        # Bulk fibers
        w.writerow(["# BULK FIBERS"])
        header_b = (["x_mm", "y_mm", "eps", "sigma_MPa", "dA_mm2"]
                    if has_x else
                    ["y_mm", "eps", "sigma_MPa", "dA_mm2"])
        w.writerow(header_b)
        for i in range(len(bulk["y"])):
            row = []
            if has_x:
                row.append(f"{bulk['x'][i]:.2f}")
            row += [
                f"{bulk['y'][i]:.2f}",
                f"{bulk['eps'][i]:.8f}",
                f"{bulk['sigma'][i]:.4f}",
                f"{bulk['dA'][i]:.2f}",
            ]
            w.writerow(row)

        w.writerow([])

        # Rebars
        rb = fiber_results["rebars"]
        has_x_r = "x" in rb
        has_net = "sigma_net" in rb
        w.writerow(["# REBAR LAYERS"])
        header_r = ["bar_#"]
        if has_x_r:
            header_r.append("x_mm")
        header_r += ["y_mm", "eps", "sigma_MPa"]
        if has_net:
            header_r.append("sigma_net_MPa")
        header_r += ["A_mm2", "F_net_kN"]
        w.writerow(header_r)
        for i in range(len(rb["y"])):
            sig_net = rb["sigma_net"][i] if has_net else rb["sigma"][i]
            F_net = sig_net * rb["A"][i] / 1e3
            row = [f"{i+1}"]
            if has_x_r:
                row.append(f"{rb['x'][i]:.2f}")
            row += [
                f"{rb['y'][i]:.2f}",
                f"{rb['eps'][i]:.8f}",
                f"{rb['sigma'][i]:.4f}",
            ]
            if has_net:
                row.append(f"{sig_net:.4f}")
            row += [
                f"{rb['A'][i]:.2f}",
                f"{F_net:.4f}",
            ]
            w.writerow(row)


def export_3d_surface_csv(nm_3d, filepath):
    """
    Export 3D resistance surface (N, Mx, My) to CSV.

    Parameters
    ----------
    nm_3d : dict
        Output of :meth:`NMDiagram.generate_biaxial`.
    filepath : str
    """
    with open(filepath, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["N_kN", "Mx_kNm", "My_kNm"])
        for n, mx, my in zip(nm_3d["N_kN"], nm_3d["Mx_kNm"],
                              nm_3d["My_kNm"]):
            w.writerow([f"{n:.4f}", f"{mx:.6f}", f"{my:.6f}"])


def export_3d_surface_json(nm_3d, filepath):
    """
    Export 3D resistance surface (N, Mx, My) to JSON.

    Parameters
    ----------
    nm_3d : dict
        Output of :meth:`NMDiagram.generate_biaxial`.
    filepath : str
    """
    data = {
        "description": "N-Mx-My resistance surface point cloud",
        "units": {"N": "kN", "Mx": "kN*m", "My": "kN*m"},
        "n_points": len(nm_3d["N_kN"]),
        "N_kN": [round(float(x), 4) for x in nm_3d["N_kN"]],
        "Mx_kNm": [round(float(x), 6) for x in nm_3d["Mx_kNm"]],
        "My_kNm": [round(float(x), 6) for x in nm_3d["My_kNm"]],
    }
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
