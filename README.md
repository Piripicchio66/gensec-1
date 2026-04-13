# GenSec — Generic Section Calculator

Fiber-based cross-section analysis for composite structural members under combined axial force and biaxial bending.

GenSec computes the **resistance domain** (N-Mx-My interaction surface) of any cross-section composed of any combination of materials, and verifies arbitrary load demands against it.

## Features

- **Material-agnostic**: concrete (parabola-rectangle EC2), reinforcing steel (elastic-plastic with hardening), structural steel (EN 10025-2), CFRP, timber, or any material via tabulated σ-ε curves.
- **Biaxial bending**: full (N, Mx, My) interaction, with 3D resistance surface generation and convex-hull-based demand verification.
- **Mx-My diagrams**: interaction contour at any fixed N, with demand points plotted.
- **Uniaxial bending**: classic N-M interaction diagram as a special case.
- **Utilization ratios**: for each load demand, computes the ratio η = distance to domain boundary.
- **Combinations**: named groups of (N, Mx, My) triples for envelope verification (e.g. many columns with the same section).
- **Per-fiber post-processing**: strain and stress at every fiber and rebar, for any given load combination.
- **YAML-driven input**: section geometry, materials, and load demands defined in a single YAML file.
- **CLI and API**: usable from the command line or as a Python library.
- **EC2/NTC 2018 integration**: automatic computation of all Table 3.1 parameters from fck, including high-strength concrete (fck > 50 MPa).

## Installation

Requires Python ≥ 3.10.

```bash
# With uv (recommended)
uv sync

# Or with pip
pip install -e .
```

## Quick start

### From the command line

```bash
# Run analysis on a YAML input file
uv run gensec examples/example_input.yaml --output-dir results

# Or equivalently
uv run python -m gensec examples/example_input.yaml --output-dir results
```

This produces:
- `nm_diagram.png` — N-M interaction diagram with demand points
- `nm_domain.csv` / `.json` — point cloud of the resistance domain
- `demand_summary.csv` / `.json` — verification table with utilization ratios
- `fibers_<name>.csv` — per-fiber strain/stress for each demand
- `stress_<name>.png` — strain and stress profile plots

### From Python

```python
from gensec.materials import Concrete, Steel
from gensec.geometry import RebarLayer, RectSection
from gensec.solver import FiberSolver, NMDiagram, DemandChecker
import numpy as np

# Define materials
concrete = Concrete(fck=25.0, gamma_c=1.5, alpha_cc=0.85)
steel = Steel(fyk=450.0, gamma_s=1.15)

# Define section: 300×600 mm column, 3Φ20 top + 3Φ20 bottom
A20 = np.pi / 4 * 20**2
rebars = [
    RebarLayer(y=40,  As=3*A20, material=steel, n_bars=3, diameter=20),
    RebarLayer(y=560, As=3*A20, material=steel, n_bars=3, diameter=20),
]
section = RectSection(B=300, H=600, bulk_material=concrete, rebars=rebars)

# Create solver and generate N-M diagram
solver = FiberSolver(section)
nm = NMDiagram(solver).generate(n_points=400)

# Verify a load demand
checker = DemandChecker(nm)
eta = checker.utilization_ratio(N=-1500e3, Mx_or_M=200e6)
print(f"Utilization ratio: {eta:.3f}")  # < 1.0 = verified

# Get detailed stress state for a specific demand
sol = solver.solve_equilibrium(N_target=-1500e3, Mx_target=200e6)
results = solver.get_fiber_results(sol["eps0"], sol["chi_x"])
```

### Using EC2 properties from fck

```python
from gensec.materials import concrete_from_ec2, concrete_from_class

# From fck directly (French National Annex)
c30 = concrete_from_ec2(fck=30, ls='F', loadtype='slow', TypeConc='R')
# c30.fcd, c30.eps_c2, c30.eps_cu2, c30.n_parabola are all computed

# From class name
c50 = concrete_from_class('C50/60', ls='F')

# Access full EC2 properties
print(c30.ec2.fcm)   # 38.0 MPa
print(c30.ec2.fctm)  # 2.9 MPa
print(c30.ec2.ecm)   # 32837 MPa
```

### Using structural steel (EN 10025-2)

```python
from gensec.materials import steel_from_en10025

s355 = steel_from_en10025('S355', t=20)  # 20mm plate
# s355.fyk = 345 (thickness-dependent)
# s355.k_hardening = 470/345 (fu/fy ratio)
```

### Biaxial bending

```python
# For biaxial, place rebars with explicit x coordinates
# and use a 2D fiber grid
rebars = [
    RebarLayer(y=40, x=40,  As=A20, material=steel),
    RebarLayer(y=40, x=260, As=A20, material=steel),
    RebarLayer(y=560, x=40,  As=A20, material=steel),
    RebarLayer(y=560, x=260, As=A20, material=steel),
]
section = RectSection(
    B=300, H=600, bulk_material=concrete,
    rebars=rebars, n_fibers_y=60, n_fibers_x=30,
)

solver = FiberSolver(section)
nm_3d = NMDiagram(solver).generate_biaxial(n_angles=36, n_points_per_angle=200)

# 3D demand checker
checker = DemandChecker(nm_3d)
eta = checker.utilization_ratio(N=-1000e3, Mx_or_M=100e6, My=50e6)
```

### Mixed materials (e.g., RC + CFRP strip)

```python
from gensec.materials import TabulatedMaterial

cfrp = TabulatedMaterial(
    strains=[0.0, 0.017],
    stresses=[0.0, 2800.0],
    name="CFRP",
)

rebars = [
    RebarLayer(y=40, As=942, material=steel, embedded=True),
    RebarLayer(y=560, As=942, material=steel, embedded=True),
    RebarLayer(y=5, As=150, material=cfrp, embedded=False),  # external strip
]
```

## YAML input format

```yaml
materials:
  concrete_1:
    type: concrete
    fck: 25.0
    gamma_c: 1.5
    alpha_cc: 0.85

  steel_1:
    type: steel
    fyk: 450.0
    gamma_s: 1.15
    works_in_compression: true

  cfrp:
    type: tabulated
    name: CFRP_sheet
    strains: [0.0, 0.017]
    stresses: [0.0, 2800.0]

section:
  B: 300
  H: 600
  bulk_material: concrete_1
  n_fibers_y: 100
  n_fibers_x: 1        # 1 for uniaxial, >1 for biaxial
  rebars:
    - y: 40
      As: 942.5
      material: steel_1
      embedded: true     # default: true
      n_bars: 3
      diameter: 20
    - y: 560
      As: 942.5
      material: steel_1
      x: 150             # optional: x-coordinate for biaxial
    - y: 5
      As: 150
      material: cfrp
      embedded: false     # external element

demands:
  - name: Gravity
    N_kN: -1500
    Mx_kNm: 200
    My_kNm: 0
  - name: Seismic_biax
    N_kN: -1200
    Mx_kNm: 200
    My_kNm: 80

# Combinations: groups of (N, Mx, My) triples.
# Use case: same section for many columns in a building.
combinations:
  - name: Seismic_X_envelope
    demands:
      - {N_kN: -1800, Mx_kNm: 180, My_kNm: 30}
      - {N_kN: -1500, Mx_kNm: 250, My_kNm: 45}
      - {N_kN: -1200, Mx_kNm: 300, My_kNm: 20}
```

**Demands** always use `Mx_kNm` and `My_kNm` (both required). The legacy field `M_kNm` is accepted as `Mx_kNm` with `My_kNm=0`.

**Combinations** group many triples under a name. GenSec reports the utilization ratio η for every triple, plus η_max per combination.

## Sign conventions

| Quantity | Positive | Negative |
|----------|----------|----------|
| ε (strain) | Tension | Compression |
| N (axial force) | Tension | Compression |
| Mx (moment about x) | Bottom edge compressed | Top edge compressed |
| My (moment about y) | Left edge compressed | Right edge compressed |

Coordinate system: origin at the bottom-left corner of the section, x-axis horizontal (width), y-axis vertical (height).

Moment convention follows directly from the integration:

    Mx = Σ σᵢ · Aᵢ · (yᵢ − y_ref)
    My = Σ σᵢ · Aᵢ · (xᵢ − x_ref)

## Embedded vs external rebars

For rebars **embedded** within the bulk material (the default), GenSec subtracts the displaced bulk material to avoid double-counting:

    F_net = [σ_rebar(ε) − σ_bulk(ε)] · As

For **external** elements (e.g., CFRP strips bonded to the surface, steel truss chords outside the concrete), set `embedded: false` — the full rebar stress is used:

    F_net = σ_rebar(ε) · As

## Running tests

```bash
uv run python -m pytest tests/ -v

# Or with unittest directly
uv run python -m unittest discover -s tests -v
```

106 tests covering:
- Material constitutive laws (pointwise verification)
- Analytical comparisons (pure compression, tension, bending)
- Solver convergence (systematic grid test)
- Sign conventions
- Embedded bar subtraction
- Biaxial bending with analytical rebar verification
- Asymmetric sections
- Mixed materials
- Mx-My diagram generation and symmetry
- Combinations loading from YAML
- Solver dispatch (My=0 on biaxial section → uniaxial fallback)
- YAML loading
- Export round-trips

## Architecture

```
src/gensec/
├── materials/
│   ├── base.py            # Abstract Material class
│   ├── concrete.py        # Parabola-rectangle (EC2)
│   ├── steel.py           # Elastic-plastic with hardening
│   ├── tabulated.py       # Arbitrary σ-ε from data points
│   ├── ec2_properties.py  # Full EC2 Table 3.1
│   ├── en10025_properties.py  # EN 10025-2 structural steel
│   └── ec2_bridge.py      # Factory functions: EC2 → GenSec
├── geometry/
│   ├── fiber.py           # RebarLayer (point fibers)
│   └── section.py         # RectSection (2D fiber grid)
├── solver/
│   ├── integrator.py      # FiberSolver: direct + inverse problem
│   ├── capacity.py        # NMDiagram: N-M and N-Mx-My surfaces
│   └── check.py           # DemandChecker: utilization ratios
├── output/
│   ├── report.py          # Terminal reporting
│   ├── plots.py           # matplotlib diagrams
│   └── export.py          # CSV / JSON export
├── io_yaml.py             # YAML input loader
├── cli.py                 # Command-line interface
└── __main__.py            # python -m gensec entry point
```

## Roadmap

- **Phase 3**: Multiple load combinations, life stages, seismic utilization vectors
- **Phase 4**: Mechanical and inertial properties (area, centroid, moments of inertia)
- **Phase 5**: Prestress (cumulative strain states)
- **Phase 6**: Shear and interface slip

## License

GenSec is licensed under the **GNU Affero General Public License v3.0**
(AGPL-3.0).  See [LICENSE](LICENSE) for the full text.

### What This Means for You

| Use case | License needed | Cost |
|---|---|---|
| Internal engineering calculations | AGPL-3.0 | Free |
| Academic research and teaching | AGPL-3.0 | Free |
| Open-source projects (AGPL-compatible) | AGPL-3.0 | Free |
| Embedding in proprietary software | Commercial | Negotiated |
| Offering as SaaS / cloud service | Commercial | Negotiated |
| Redistributing in commercial products | Commercial | Negotiated |

A **commercial license** is available for organizations that cannot comply
with the AGPL-3.0 terms.  See [COMMERCIAL_LICENSE.md](COMMERCIAL_LICENSE.md)
for details.

### Attribution

Any use of GenSec — including in reports, publications, and presentations —
must include the following attribution:

> Cross-section analysis performed with GenSec
> (https://alberoandrea.com) by Andrea Albero.

### Contributing

Contributions are welcome!  Please read [CONTRIBUTING.md](CONTRIBUTING.md)
and sign the [CLA](CLA.md) before submitting a pull request.

