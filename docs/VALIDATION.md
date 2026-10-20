# Validation protocol

## Current status

Cavity Flow Lab is an experimental two-dimensional incompressible lid-driven cavity solver. It is not yet a validated CFD package. Software tests, discrete mass conservation, reference-data agreement, grid convergence, and runtime are distinct questions.

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
