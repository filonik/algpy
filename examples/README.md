# Example notebooks

These notebooks use algpy's Python/NumPy API to connect algebra operations to
geometry and differentiation. Run each notebook from top to bottom. Plotly
figures include interactive sliders, hover information, and 3D orbiting.

| Notebook | Demonstration |
| --- | --- |
| [01 — Algebras and functions](01_algebras_and_functions.ipynb) | Scalars and arrays, elementary functions, powers, multiplication matrices, all three five-generator families, and a fifth mixed derivative. |
| [02 — Projective geometry](02_projective_geometry.ipynb) | PGA2 line/point and PGA3 plane/line constructions, meets, joins, and motor interpolation. |
| [03 — Conformal geometry](03_conformal_geometry.ipynb) | CGA2 circle/line intersections and CGA3 sphere/plane intersections, null points, translation, rotation, dilation, and interpolation. |
| [04 — Curves and dual numbers](04_curves_and_dual_numbers.ipynb) | A planar curve's 2×1 Jacobian, tangent, normal, second derivatives, curvature, and osculating circle. |
| [05 — Surfaces and curvature](05_surfaces_and_curvature.ipynb) | A surface's 3×2 Jacobian, coordinate Hessians, fundamental forms, and Gaussian/mean curvature, with a vectorized torus example. |
| [06 — Exterior operations](06_exterior_opns.ipynb) | Lines through points and their intersections with `exterior3`, circles through three points and line-circle intersections with `exterior4`, drawn shader-style by evaluating incidence on whole pixel grids. |

## Setup

From the repository root, install the library and notebook dependencies into
your Python environment:

```sh
python -m pip install -e ".[examples]"
cd examples
python -m jupyterlab
```

Select that environment's Python kernel. Starting in `examples/` makes the
local helper modules available. The notebooks are saved without outputs;
**Run All** computes their numerical checks and interactive figures. Figures
use `plotly.offline.init_notebook_mode(connected=False)` and `fig.show(auto_play=False)`.
The shared `_plotting.slider_controls` helper selects precomputed frames with
zero transition duration. There are no playback buttons; the slider track has
no ticks or readout. Hover over the geometry to inspect numerical values.
All axes use fixed ranges, and the 3D figures use `aspectmode="cube"`.

The figures embed Plotly.js, so the sliders also work in exported HTML without
a running Python kernel or CDN connection. The 3D figures use WebGL.

After updating these examples, **restart the kernel and clear all outputs**, then
**Run All**. This reloads the shared display helper and removes outputs created
by the previous renderer.

## Mathematical conventions

- PGA uses the dual model: grade-one elements are lines in 2D and planes in
  3D. Geometric meets are exterior products, `a | b`, and joins are
  regressive products, `a & b`, using the PGA types' lattice API.
  The null generator comes first.
- CGA uses outer-product point constructions, `a | b`, and intersects them
  with the regressive product, `a & b`. The two special generators have
  squares `-1, +1` and come before the Euclidean generators. Round objects are
  decoded with the nondegenerate metric dual. Example intersections are finite
  and real.
- The geometry notebooks use small helpers in [_geometry.py](_geometry.py)
  that construct points, lines, planes, circles, and spheres and decode their
  coordinates, centers, and radii. All other operations come from the algebra
  types.
- Independent commuting dual generators carry repeated derivative directions.
  The notebooks explain the coefficient indices and the absence of factorial
  rescaling when extracting the displayed derivatives.
- Surface mean curvature uses the outward normal and shape operator `I⁻¹ II`.
  Its sign changes if the normal orientation is reversed; Gaussian curvature
  does not.

## Execute and check all examples

From the repository root, with the examples environment active:

```sh
python examples/execute.py
```

Each notebook runs in a fresh kernel, including its assertions against analytic
formulas and incidence identities. To retain executed copies without changing
the source notebooks, choose an output directory:

```sh
python examples/execute.py --output-dir /tmp/algpy-executed
```

Use `--kernel NAME` to select another installed Jupyter kernel. Individual
notebook paths can also be supplied. To export an executed notebook to HTML:

```sh
python -m jupyter nbconvert --to html /tmp/algpy-executed/03_conformal_geometry.ipynb
```
