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
Command-line interface for GenSec.

Usage::

    uv run gensec input.yaml [--n-points 400] [--output-dir ./results]
"""

import argparse
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .io_yaml import load_yaml
from .solver import FiberSolver, NMDiagram, DemandChecker
from .output import (
    print_section_info, print_fiber_results,
    plot_nm_diagram, plot_stress_profile, plot_mx_my_diagram,
    plot_moment_curvature, plot_section, plot_section_state,
    plot_demand_heatmap, plot_3d_surface,
    plot_moment_curvature_bundle, plot_polar_ductility,
    plot_moment_curvature_surface,
    export_nm_domain_csv, export_nm_domain_json,
    export_demand_results_csv, export_demand_results_json,
    export_fiber_results_csv,
    export_3d_surface_csv, export_3d_surface_json,
)


def _print_verification_table(title, check_results):
    """Print a formatted verification summary."""
    print(f"\n{'=' * 75}")
    print(f"  {title}")
    print(f"{'=' * 75}")
    has_my = any("My_kNm" in r for r in check_results)
    if has_my:
        header = (f"  {'Name':>20} {'N[kN]':>8} {'Mx[kNm]':>9}"
                  f" {'My[kNm]':>9} {'eta':>7} {'Status':>8}")
    else:
        header = (f"  {'Name':>20} {'N[kN]':>8} {'Mx[kNm]':>9}"
                  f" {'eta':>7} {'Status':>8}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    all_ok = True
    for r in check_results:
        status = "OK" if r["verified"] else "FAIL"
        if not r["verified"]:
            all_ok = False
        if has_my:
            print(f"  {r['name']:>20} {r['N_kN']:>8.1f}"
                  f" {r['M_kNm']:>9.2f} {r.get('My_kNm', 0):>9.2f}"
                  f" {r['utilization']:>7.3f} {status:>8}")
        else:
            print(f"  {r['name']:>20} {r['N_kN']:>8.1f}"
                  f" {r['M_kNm']:>9.2f}"
                  f" {r['utilization']:>7.3f} {status:>8}")

    print("  " + "-" * (len(header) - 2))
    if all_ok:
        print("  All demands VERIFIED.")
    else:
        n_fail = sum(1 for r in check_results if not r["verified"])
        print(f"  *** {n_fail} demand(s) NOT VERIFIED ***")
    print("=" * 75)


def main(argv=None):
    """GenSec CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="gensec",
        description="GenSec — fiber-based cross-section analysis.",
    )
    parser.add_argument("input_file", help="YAML input file.")
    parser.add_argument("--n-points", type=int, default=400,
                        help="N-M diagram resolution (default: 400).")
    parser.add_argument("--output-dir", type=str, default=".",
                        help="Output directory (default: current dir).")
    args = parser.parse_args(argv)
    outdir = args.output_dir
    os.makedirs(outdir, exist_ok=True)

    # --- Load ---
    print(f"\nLoading: {args.input_file}")
    data = load_yaml(args.input_file)
    section = data["section"]
    demands = data["demands"]
    combinations = data.get("combinations", [])

    # Output flags from YAML
    output_opts = data.get("output_options", {})
    do_mx_my = output_opts.get("generate_mx_my", False)
    do_3d_surface = output_opts.get("generate_3d_surface", False)

    # N-level strategy for M-chi and Mx-My diagrams.
    # Options in YAML (under output):
    #   n_levels_mode: "demands"    -> one per demand N (default)
    #                  "auto"       -> auto-discretize from Nmin to Nmax
    #                  "explicit"   -> user-specified list
    #   n_levels_count: 10          -> number of levels for "auto"
    #   n_levels_values: [-3000, -1500, 0]  -> explicit N values [kN]
    n_levels_mode = output_opts.get("n_levels_mode", "demands")
    n_levels_count = int(output_opts.get("n_levels_count", 10))
    n_levels_explicit = output_opts.get("n_levels_values", [])
    # Angular resolution for Mx-My contour (higher = smoother)
    n_angles_mx_my = int(output_opts.get("n_angles_mx_my", 180))

    print_section_info(section)

    # --- Solver ---
    solver = FiberSolver(section)
    is_biaxial = section.n_fibers_x > 1

    # --- N-M diagram (always uniaxial N-Mx) ---
    print("\nGenerating N-Mx diagram...")
    nm_gen = NMDiagram(solver)
    nm_data = nm_gen.generate(n_points=args.n_points)
    print(f"  {len(nm_data['N'])} points")
    print(f"  N:  [{nm_data['N_kN'].min():.1f},"
          f" {nm_data['N_kN'].max():.1f}] kN")
    print(f"  Mx: [{nm_data['M_kNm'].min():.1f},"
          f" {nm_data['M_kNm'].max():.1f}] kN*m")

    export_nm_domain_csv(nm_data, os.path.join(outdir, "n_mx_domain.csv"))
    export_nm_domain_json(nm_data, os.path.join(outdir, "n_mx_domain.json"))

    # --- N-My diagram (only for biaxial sections) ---
    nm_data_y = None
    if is_biaxial:
        print("\nGenerating N-My diagram...")
        nm_data_y = nm_gen.generate(n_points=args.n_points, direction='y')
        print(f"  {len(nm_data_y['N'])} points")
        print(f"  N:  [{nm_data_y['N_kN'].min():.1f},"
              f" {nm_data_y['N_kN'].max():.1f}] kN")
        print(f"  My: [{nm_data_y['M_kNm'].min():.1f},"
              f" {nm_data_y['M_kNm'].max():.1f}] kN*m")
        export_nm_domain_csv(nm_data_y,
                             os.path.join(outdir, "n_my_domain.csv"))

    # --- Biaxial surface if section supports it ---
    nm_3d = None
    has_biaxial_demands = (
        any(abs(d.get("My", 0)) > 0 for d in demands)
        or any(
            any(abs(dd.get("My", 0)) > 0 for dd in c["demands"])
            for c in combinations
        )
    )

    if is_biaxial:
        print("\nGenerating 3D resistance surface (N-Mx-My)...")
        nm_3d = nm_gen.generate_biaxial(
            n_angles=72, n_points_per_angle=args.n_points // 2)
        print(f"  {len(nm_3d['N'])} points")

        # Export 3D surface if requested
        if do_3d_surface:
            export_3d_surface_csv(
                nm_3d, os.path.join(outdir, "surface_3d.csv"))
            export_3d_surface_json(
                nm_3d, os.path.join(outdir, "surface_3d.json"))
            print("  Exported: surface_3d.csv / .json")
    elif has_biaxial_demands:
        print("\n  WARNING: Demands with My != 0 found, but section has"
              " n_fibers_x=1 (uniaxial mode).")
        print("  Set n_fibers_x > 1 in the YAML for biaxial analysis.")
        print("  My will be IGNORED for verification; only N-Mx checked.")

    # --- Collect ALL demand triples for checking ---
    all_demands = list(demands)
    for combo in combinations:
        for i, d in enumerate(combo["demands"]):
            # Tag with combination name for reporting
            d_tagged = dict(d)
            d_tagged["name"] = f"{combo['name']}[{i+1}]"
            all_demands.append(d_tagged)

    # --- Demand verification ---
    if all_demands:
        # Choose 2D or 3D checker
        if nm_3d is not None:
            checker = DemandChecker(nm_3d)
        else:
            checker = DemandChecker(nm_data)

        check_results = checker.check_demands(all_demands)
        _print_verification_table("DEMAND VERIFICATION", check_results)

        # Export summary
        export_demand_results_csv(
            check_results, os.path.join(outdir, "demand_summary.csv"))
        export_demand_results_json(
            check_results, os.path.join(outdir, "demand_summary.json"))

        # --- Per-combination summary ---
        if combinations:
            print(f"\n{'=' * 75}")
            print("  COMBINATION SUMMARY")
            print(f"{'=' * 75}")
            print(f"  {'Combination':>25} {'n_demands':>10}"
                  f" {'eta_max':>8} {'Status':>8}")
            print("  " + "-" * 55)
            for combo in combinations:
                combo_results = [
                    r for r in check_results
                    if r["name"].startswith(combo["name"] + "[")
                ]
                if combo_results:
                    eta_max = max(r["utilization"] for r in combo_results)
                    n_d = len(combo_results)
                    status = "OK" if eta_max <= 1.0 else "FAIL"
                    print(f"  {combo['name']:>25} {n_d:>10}"
                          f" {eta_max:>8.3f} {status:>8}")
            print(f"{'=' * 75}")

    # --- Per-demand fiber details (individual demands only) ---
    # --- Demand utilization heatmap ---
    if all_demands and check_results:
        fig_hm = plot_demand_heatmap(check_results)
        p = os.path.join(outdir, "demand_utilization.png")
        fig_hm.savefig(p, dpi=150)
        plt.close(fig_hm)
        print(f"\n  Exported: {p}")

    demand_plot = []
    if demands:
        print("\n  Per-demand details:")
        for d in demands:
            N_d, Mx_d, My_d = d["N"], d["Mx"], d["My"]
            sol = solver.solve_equilibrium(N_d, Mx_d, My_d)
            label = d["name"]
            demand_plot.append((N_d / 1e3, Mx_d / 1e6, label))

            if sol["converged"]:
                fr = solver.get_fiber_results(
                    sol["eps0"], sol["chi_x"], sol["chi_y"])
                print(f"\n  --- {label} ---")
                print_fiber_results(fr, section)

                p = os.path.join(outdir, f"fibers_{label}.csv")
                export_fiber_results_csv(fr, p)
                print(f"    Exported: {p}")

                fig = plot_stress_profile(
                    fr, section,
                    title=f"{label}: N={N_d/1e3:.0f},"
                          f" Mx={Mx_d/1e6:.0f}, My={My_d/1e6:.0f} kN·m")
                p = os.path.join(outdir, f"state_{label}.png")
                fig.savefig(p, dpi=150)
                plt.close(fig)
                print(f"    Exported: {p}")
            else:
                print(f"\n  --- {label} --- SOLVER DID NOT CONVERGE")

    # --- N-M diagram plot ---
    fig_nm = plot_nm_diagram(nm_data, demand_plot)
    p = os.path.join(outdir, "n_mx_diagram.png")
    fig_nm.savefig(p, dpi=150)
    plt.close(fig_nm)
    print(f"\n  Exported: {p}")

    # --- N-My diagram plot (biaxial only) ---
    if nm_data_y is not None:
        demand_plot_y = [
            (d["N"] / 1e3, d["My"] / 1e6, d["name"])
            for d in demands if abs(d["My"]) > 0
        ]
        fig_nmy = plot_nm_diagram(nm_data_y, demand_plot_y or None,
                                  title="N-My Interaction Diagram")
        p = os.path.join(outdir, "n_my_diagram.png")
        fig_nmy.savefig(p, dpi=150)
        plt.close(fig_nmy)
        print(f"  Exported: {p}")

    # --- 3D resistance surface plot (biaxial only) ---
    if nm_3d is not None and do_3d_surface:
        print("\n  Generating 3D surface plot...")
        fig_3d = plot_3d_surface(nm_3d, demands=all_demands)
        p = os.path.join(outdir, "surface_3d.png")
        fig_3d.savefig(p, dpi=150)
        plt.close(fig_3d)
        print(f"  Exported: {p}")

    # ==============================================================
    #  Determine N levels for M-chi and Mx-My diagrams
    # ==============================================================
    if n_levels_mode == "explicit" and n_levels_explicit:
        # User-specified N values [kN] → convert to [N]
        N_levels_analysis = sorted(
            float(v) * 1e3 for v in n_levels_explicit)
        print(f"\n  N levels (explicit): "
              f"{[v/1e3 for v in N_levels_analysis]} kN")

    elif n_levels_mode == "auto":
        # Auto-discretize from N_min to N_max of the domain
        N_min_domain = float(nm_data["N"].min())
        N_max_domain = float(nm_data["N"].max())
        N_levels_analysis = sorted(
            np.linspace(N_min_domain, N_max_domain,
                        n_levels_count).tolist())
        print(f"\n  N levels (auto, {n_levels_count} values): "
              f"[{N_min_domain/1e3:.0f}, {N_max_domain/1e3:.0f}] kN")

    else:
        # Default: one per unique demand N
        N_levels_analysis = sorted(set(d["N"] for d in demands))
        print(f"\n  N levels (from demands): "
              f"{[v/1e3 for v in N_levels_analysis]} kN")

    # Always include N=0 (pure bending) if not already present
    if not any(abs(n) < 1.0 for n in N_levels_analysis):
        N_levels_analysis.append(0.0)
        N_levels_analysis.sort()
        print(f"  Added N=0 (pure bending)")

    # --- Mx-My diagrams at fixed N (biaxial sections) ---
    if is_biaxial and do_mx_my:
        for N_fixed in N_levels_analysis:
            print(f"  Generating Mx-My at N={N_fixed/1e3:.0f} kN...")
            mx_my = nm_gen.generate_mx_my(
                N_fixed, n_angles=n_angles_mx_my,
                n_points_per_angle=args.n_points // 2)

            mx_my_demands = [
                (d["Mx"] / 1e6, d["My"] / 1e6, d["name"])
                for d in all_demands
                if abs(d["N"] - N_fixed) < abs(N_fixed) * 0.05 + 1.0
            ]

            fig = plot_mx_my_diagram(mx_my, mx_my_demands or None)
            N_label = f"{N_fixed/1e3:.0f}".replace("-", "m")
            p = os.path.join(outdir, f"mx_my_N{N_label}.png")
            fig.savefig(p, dpi=150)
            plt.close(fig)
            print(f"    Exported: {p}")

    # --- Moment-curvature diagrams (Mx-chi and optionally My-chi) ---
    mc_collection_x = []
    mc_collection_y = []
    for N_fixed in N_levels_analysis:
        N_label = f"{N_fixed/1e3:.0f}".replace("-", "m")

        # Mx-chi (always)
        print(f"  Generating Mx-χ at N={N_fixed/1e3:.0f} kN...")
        mc_x = nm_gen.generate_moment_curvature(
            N_fixed, n_points=args.n_points, direction='x')
        mc_collection_x.append(mc_x)
        fig = plot_moment_curvature(mc_x)
        p = os.path.join(outdir, f"mx_chi_N{N_label}.png")
        fig.savefig(p, dpi=150)
        plt.close(fig)
        print(f"    Exported: {p}")

        # My-chi (biaxial only)
        if is_biaxial:
            print(f"  Generating My-χ at N={N_fixed/1e3:.0f} kN...")
            mc_y = nm_gen.generate_moment_curvature(
                N_fixed, n_points=args.n_points, direction='y')
            mc_collection_y.append(mc_y)
            fig = plot_moment_curvature(mc_y)
            p = os.path.join(outdir, f"my_chi_N{N_label}.png")
            fig.savefig(p, dpi=150)
            plt.close(fig)
            print(f"    Exported: {p}")

    # --- A) Bundle of Mx-chi curves at varying N ---
    if len(mc_collection_x) > 1:
        print("\n  Generating Mx-χ bundle...")
        fig = plot_moment_curvature_bundle(mc_collection_x, direction='x')
        p = os.path.join(outdir, "mx_chi_bundle.png")
        fig.savefig(p, dpi=150)
        plt.close(fig)
        print(f"  Exported: {p}")

    if len(mc_collection_y) > 1:
        print("  Generating My-χ bundle...")
        fig = plot_moment_curvature_bundle(mc_collection_y, direction='y')
        p = os.path.join(outdir, "my_chi_bundle.png")
        fig.savefig(p, dpi=150)
        plt.close(fig)
        print(f"  Exported: {p}")

    # --- B) Polar ductility diagrams ---
    if is_biaxial:
        for N_fixed in N_levels_analysis:
            N_label = f"{N_fixed/1e3:.0f}".replace("-", "m")
            print(f"  Generating polar ductility at N={N_fixed/1e3:.0f} kN...")
            fig = plot_polar_ductility(
                nm_gen, N_fixed, n_angles=n_angles_mx_my,
                n_points=args.n_points)
            p = os.path.join(outdir, f"polar_ductility_N{N_label}.png")
            fig.savefig(p, dpi=150)
            plt.close(fig)
            print(f"    Exported: {p}")

    # --- C) 3D surface Mx-chi-N ---
    if len(mc_collection_x) > 1:
        print("\n  Generating 3D Mx-χ-N surface...")
        fig = plot_moment_curvature_surface(mc_collection_x, direction='x')
        p = os.path.join(outdir, "mx_chi_N_surface.png")
        fig.savefig(p, dpi=150)
        plt.close(fig)
        print(f"  Exported: {p}")

    if len(mc_collection_y) > 1:
        print("  Generating 3D My-χ-N surface...")
        fig = plot_moment_curvature_surface(mc_collection_y, direction='y')
        p = os.path.join(outdir, "my_chi_N_surface.png")
        fig.savefig(p, dpi=150)
        plt.close(fig)
        print(f"  Exported: {p}")

    # --- Section drawing for each individual demand ---
    if demands:
        # Geometry-only section drawing
        fig_geom = plot_section(section, title="Section geometry")
        p = os.path.join(outdir, "section_geometry.png")
        fig_geom.savefig(p, dpi=150)
        plt.close(fig_geom)
        print(f"\n  Exported: {p}")

        for d in demands:
            N_d, Mx_d, My_d = d["N"], d["Mx"], d["My"]
            sol = solver.solve_equilibrium(N_d, Mx_d, My_d)
            if sol["converged"]:
                fr = solver.get_fiber_results(
                    sol["eps0"], sol["chi_x"], sol["chi_y"])
                label = d["name"]
                fig = plot_section(
                    section, fr,
                    title=f"{label}: N={N_d/1e3:.0f},"
                          f" Mx={Mx_d/1e6:.0f}, My={My_d/1e6:.0f} kNm")
                p = os.path.join(outdir, f"section_{label}.png")
                fig.savefig(p, dpi=150)
                plt.close(fig)
                print(f"  Exported: {p}")

    # --- Final summary ---
    print(f"\n{'=' * 75}")
    print(f"  All outputs in: {os.path.abspath(outdir)}/")
    files = sorted(os.listdir(outdir))
    for fn in files:
        size = os.path.getsize(os.path.join(outdir, fn))
        print(f"    {fn:45s}  {size:>8,} bytes")
    print("=" * 75)
    print("\nDone.")


if __name__ == "__main__":
    main()
