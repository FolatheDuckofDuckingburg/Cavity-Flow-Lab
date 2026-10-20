# Validation protocol

## Current status

Cavity Flow Lab is an experimental two-dimensional incompressible lid-driven cavity solver. It is not yet a validated CFD package. Software tests, discrete mass conservation, reference-data agreement, grid convergence, and runtime are distinct questions.

## Repair validation on `fix/projection-operator-consistency`

The merged `main` at `c05c6b9a4938e48b5238b60be7a354d6c9a94272` failed
`tests/test_solver.py::test_projection_reduces_divergence_to_poisson_tolerance`.
The exact CI traceback was:

```text
assert 3.689083450379851 <= 1e-10
```

The failure was reproducible in a clean Python 3.12.3 virtual environment with
NumPy 2.5.4: `1 failed, 18 passed`. The root causes were (1) Jacobi stagnation
at that residual after 10,000 iterations and (2) a mismatch between the
pressure Laplacian and the pressure-gradient faces that were later overwritten
by the wall boundary assignment. The repair uses projected conjugate gradients
for the same discrete operator and constructs that operator as the divergence
of the free-face pressure gradient. The repaired suite reports `20 passed`.

Reproduce the software and projection checks with:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e "[dev]"
pytest -q
python scripts/benchmark.py --grids 8 --steps 2 --repeats 1 \
  --output results/ci-benchmark.json
```

For the exact random 8x8 projection case, using cell-centred pressure,
staggered face velocities, and the infinity norm over all cells, the repaired
run measured pre-projection `29.888522548610382`, post-projection
`4.760636329592671e-11`, and pressure residual
`4.760764040963687e-11` in 33 iterations. The discrete global mass sum was
within roundoff in the time-dependent runs. All pressure and velocity fields
were finite, and the prescribed wall values were retained.

The old `main` parent reproduces a post-projection value of `8.175570219322616`
for that same random case, not the previously reported `0.0944`. On the
repaired commit, the two-step 8x8 benchmark gives post-projection divergence
`8.649134410927814e-12`. Therefore `0.0944` was not independently reproduced
from the recorded configuration; it should not be treated as a validated
result. The discrepancy is explained qualitatively by the old boundary
overwrite and inconsistent/stagnated solve, but its original provenance remains
unknown.

The operator audit verifies, to about `1e-12` on a random 8x6 field, that the
implemented pressure Laplacian equals the divergence of the exact gradient
used for correction. This is a cell-centred finite-volume divergence with
`u` shaped `(ny, nx+1)`, `v` shaped `(ny+1, nx)`, and pressure shaped
`(ny, nx)`. Pressure correction is applied only to non-wall velocity faces;
the top-row tangential lid values and side-wall values are not subsequently
changed by projection.

## Re=100 evidence and limits

The Ghia, Ghia & Shin (1982) tabulations were compared at Re=100. The
published coordinates were linearly interpolated to the solver's cell-centre
coordinates; no extrapolation or point matching was used. On a 32x32 run with
adaptive time steps, pressure tolerance `1e-8`, pressure iteration limit
`3000`, 2,000 steps, and simulated time `10.960780857927581`, the final
centreline errors were:

```text
vertical u profile:   RMSE 0.03296268185927465, max |error| 0.106016292735054
horizontal v profile: RMSE 0.01621946281505418, max |error| 0.029577323717949
```

The final pressure residual was `8.682251155641196e-09` after 150 iterations;
post-projection divergence was `4.695799304954562e-11`; all fields were finite.
The run had not met the requested velocity-update steady criterion by the
2,000-step cap (`max(du_linf,dv_linf)=1.4223919926198425e-06`), so these are
transient/finite-grid comparison results, not a claim of benchmark validation.
The source tabulations are the u and v tables attributed to Ghia et al. in
the [u-data gist](https://gist.github.com/ivan-pi/3e9326d18a366ffe6a8e5bfda6353219)
and [v-data gist](https://gist.github.com/ivan-pi/caa6c6737d36a9140fbcf2ea59c78b3c),
with the publication identified there as *Journal of Computational Physics*
48(3), 387–411 (1982).

At matched physical time near `T=0.25`, adaptive runs used 18, 39, and 86
steps on 16x16, 32x32, and 64x64 grids. Kinetic energies were respectively
`0.0339712073`, `0.0223552403`, and `0.0176576565`; post-projection divergence
was `1.22e-10`, `4.46e-11`, and `2.81e-11`. On 32x32, fixed `dt=0.01` and
`dt=0.005` gave kinetic energies `0.0223552403` and `0.0222431116` at matched
time near 0.25. These sensitivity results are finite-time diagnostics only.

The evidence supports a software-level projection repair and discrete mass
conservation for the stated convention. It does not support the stronger claim
that this implementation is an accurate or validated Navier–Stokes solver.

## Reproduce

```bash
python -m pip install -e ".[dev]"
pytest -q
python scripts/benchmark.py --grids 32 64 --steps 100 --repeats 3 --output results/benchmark.json
python scripts/grid_study.py --grids 32 64 128 --steps 1000 --output results/grid-study.json
```

Record the Git commit, Python/NumPy versions, OS/CPU, command line, Reynolds number, grid, time-step policy, steps, and output artifacts. Compare timing runs on the same machine and report repeated-run variability. Fixed step counts do not guarantee equal simulated time.

## Gate A — Software invariants

Test array dimensions, finite values, configuration validation, time advancement, boundary assignment, and output reproducibility. These tests detect software regressions but do not establish discretization accuracy.

## Gate B — Pressure projection and mass conservation

For each run, record pre- and post-projection discrete divergence, pressure residual, iteration count, and whether the pressure solve met its tolerance. The benchmark script reports these values. Independently verify the discrete Laplacian, gradient, and divergence operators before interpreting a small residual as a valid projection.

A provisional screening threshold is `||div(u)||_infinity <= 1e-6` in the solver's nondimensional units. This is a project-level diagnostic threshold, not a universal CFD standard.

A previous diagnostic reported post-projection divergence near 0.0944 on an 8x8 grid, while other runs were reported near machine precision. The discrepancy must be reproduced and explained on a single commit and configuration before claiming the projection is validated.

## Gate C — Reference benchmark

At Reynolds number 100, compare centreline velocity profiles against the published lid-driven cavity data of Ghia, Ghia & Shin (1982), *Journal of Computational Physics*, 48(3), 387–411, DOI: 10.1016/0021-9991(82)90058-4.

Reference values must be transcribed and checked against the publication before being committed. Report the source, coordinates, interpolation, RMSE, maximum absolute error, and number of compared points. Do not invent or silently adjust reference values.

## Gate D — Grid/time-step sensitivity

Run at least 32x32, 64x64, and 128x128 grids where feasible. Compare on common physical coordinates and comparable physical times or justified steady-state conditions. Repeat selected cases with a smaller time step. Do not infer formal order of accuracy from a single pair of grids.

The supplied grid-study script records fixed-step results and centreline profiles. It intentionally warns that the simulated times may differ; it does not by itself prove grid convergence.

## Gate E — Performance

Measure multiple runs on the same machine and compare only equivalent configurations that pass the chosen numerical checks. A faster run that fails a numerical diagnostic is not an accepted optimization. The current Python runtime benchmark is not a hardware-independent score.

## Limitations

The solver is restricted to a simplified 2D lid-driven cavity implementation. Wall treatment, pressure projection, spatial discretization, and benchmark agreement need further independent verification. It does not validate three-dimensional flow, turbulence models, arbitrary geometries, experimental predictions, or any general mathematical claim about Navier–Stokes regularity.
