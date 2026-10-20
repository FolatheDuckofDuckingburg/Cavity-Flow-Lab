# Validation and benchmark protocol

## Status

This is an experimental 2D incompressible lid-driven cavity solver. Software tests, mass-conservation diagnostics, benchmark agreement, grid convergence, and speed are separate questions. Passing one does not imply the others.

## Reproduce a run

From a clean checkout, record the commit SHA, Python/NumPy versions, OS/CPU, command, Re, grid, time-step policy, step count, and output files.

```bash
python -m pip install -e ".[dev]"
pytest -q
python scripts/benchmark.py --grids 32 64 --steps 100 --repeats 3 --output results/benchmark.json
```

Archive the JSON with the commit SHA. Timing is machine-specific; use repeated runs and report median/spread. The benchmark currently reports individual timings and diagnostics, not a portable performance score. Python tracemalloc may omit native NumPy allocations.

## Gate A — Software invariants

Tests cover staggered array shapes, initial state, boundary values, finite fields, and a single step. These checks are necessary but do not prove that the discretization or pressure projection is correct.

## Gate B — Projection and mass conservation

Report pre- and post-projection divergence using the same discrete operator, and report the reduction ratio. Independently recompute the pressure-Poisson residual. Record whether the pressure solver reached its iteration limit. The current benchmark only records final divergence and the residual returned by the solver; it is a screening diagnostic, not a complete Gate B test.

A provisional engineering screen is post-projection `||div(u)||_infinity <= 1e-6` in the solver's nondimensional units. This is not a universal CFD standard. A small pressure residual alone is insufficient.

**Known unresolved concern:** a prior diagnostic reported post-projection divergence around 0.0944 on an 8x8 grid, while other runs reportedly gave values near machine precision. These have not been reconciled. Do not claim the projection is validated until the same commit, operator, and configuration reproduce and explain the result.

## Gate C — Published benchmark

At Re=100, compare the horizontal velocity profile on the vertical centreline and the vertical velocity profile on the horizontal centreline with Ghia, Ghia & Shin (1982), *Journal of Computational Physics* 48(3), 387–411, DOI: 10.1016/0021-9991(82)90058-4.

Store reference values in a separately documented data file only after checking the transcription against the source. Report the reference source, sample coordinates, interpolation procedure, RMSE, maximum absolute error, and number of points. Do not invent or silently alter reference values.

## Gate D — Grid and time-step studies

Compare at least 32x32, 64x64, and 128x128 grids where feasible, plus a smaller-time-step run. Compare profiles on common coordinates. Fixed step counts do not mean equal physical times when dt depends on the grid; compare at a common documented simulated time or after justified steady-state convergence. Do not claim a formal order of accuracy without a suitable refinement analysis.

## Gate E — Performance

Compare runtime only for matching configurations and environments, with multiple repeats. Record grid, Re, steps, simulated time, pressure iteration counts, divergence, pressure residual, and finite-field status. A faster run that fails a numerical gate is not an accepted optimization.

## Required results table

For every configuration record commit SHA, environment, Re, grid, dt policy, steps, simulated time, pressure iterations, pressure residual, pre/post-projection divergence, centreline RMSE/max error (when reference data is available), runtime median/spread, and gate statuses.

## Scope

This suite concerns a 2D lid-driven cavity at selected conditions. It does not validate 3D flow, arbitrary geometries, turbulence models, physical experiments, or the general Navier–Stokes regularity problem. Agreement with one benchmark is evidence only for the tested configuration.
