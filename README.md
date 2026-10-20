# Cavity Flow Lab

**An open-source laboratory for reproducible numerical experiments with two-dimensional incompressible Navier–Stokes flow.**

Cavity Flow Lab studies the lid-driven square cavity at configurable Reynolds number. It is an educational research codebase for examining staggered-grid velocity fields, a pressure-Poisson solve, projection diagnostics, benchmark comparisons, and grid sensitivity.

> **Research status:** Experimental. The numerical discretization, pressure projection, and wall-boundary treatment require independent verification. Passing software tests does not establish a validated CFD solution.

## Current capabilities

- MAC-style staggered velocity arrays and cell-centred pressure.
- Explicit advection/diffusion predictor and pressure correction.
- Adaptive explicit time-step estimate, with an optional upper-bound time step.
- Pressure iteration diagnostics and pre/post-projection divergence.
- CSV centreline profiles, run metadata, speed contour, and streamlines.
- Fixed-step performance/diagnostic benchmark and grid-study script.
- Pytest suite and GitHub Actions workflow.

## Repository layout

```text
.github/workflows/tests.yml
docs/VALIDATION.md
scripts/benchmark.py
scripts/grid_study.py
src/cavity_flow/__init__.py
src/cavity_flow/solver.py
src/cavity_flow/cli.py
tests/test_solver.py
tests/test_solver_invariants.py
pyproject.toml
README.md
LICENSE
```

## Installation

Requires Python 3.10 or later.

```bash
git clone https://github.com/FolatheDuckofDuckingburg/Cavity-Flow-Lab.git
cd Cavity-Flow-Lab
python -m venv .venv
```

Activate the environment:

**Windows PowerShell**
```powershell
.venv\Scripts\Activate.ps1
```

**Linux/macOS**
```bash
source .venv/bin/activate
```

Install the package and test dependencies:

```bash
python -m pip install -e ".[dev]"
pytest -q
```

## Run a simulation

```bash
cavity-flow --re 100 --nx 64 --ny 64 --steps 3000 --output results/re100_64
```

Outputs include a history CSV, centreline velocity CSVs, plots, and `run_metadata.json`. For a velocity-update stopping heuristic, use `--until-steady`; this criterion is a numerical heuristic, not proof of physical steady state.

## Run the diagnostic benchmark

```bash
python scripts/benchmark.py --grids 32 64 --steps 100 --repeats 3 --output results/benchmark.json
```

The benchmark records wall time, simulated time, pressure residual/convergence, pre/post-projection divergence, kinetic energy, and finite-field checks. Runtime is machine-specific. A failed numerical gate is a result to investigate, not a reason to suppress the diagnostic.

## Run a grid study

```bash
python scripts/grid_study.py --grids 32 64 128 --steps 1000 --output results/grid-study.json
```

The script records centreline profiles and simulated time. Since time steps may depend on grid spacing, equal step counts do not necessarily represent equal physical times. Compare the recorded times before drawing convergence conclusions.

## Validation plan

See [docs/VALIDATION.md](docs/VALIDATION.md). In summary:

1. Verify software invariants and boundary assignments.
2. Independently verify pressure-Poisson, gradient, and divergence consistency.
3. Compare centreline profiles against Ghia, Ghia & Shin (1982), DOI: [10.1016/0021-9991(82)90058-4](https://doi.org/10.1016/0021-9991(82)90058-4).
4. Study grid and time-step sensitivity.
5. Benchmark performance only on comparable configurations that pass numerical checks.

Reference data must be transcribed and verified against the publication before use. The repository does not claim that benchmark agreement has already been achieved.

## Known limitations

- The current wall-value convention is simplified, including a lid/side-wall corner conflict.
- The pressure solve uses projected conjugate gradients and may be slow at finer grids.
- The discrete operators and boundary treatment are documented and regression-tested, but physical accuracy still requires independent verification.
- Grid-study output alone does not prove convergence.
- A small pressure residual does not guarantee a correct velocity field.
- This 2D cavity code does not validate 3D flows, arbitrary geometries, turbulence models, experimental predictions, or the general Navier–Stokes regularity problem.

## Contributing

Numerical changes should include tests and reproducible evidence. Preserve failed runs and document limitations rather than hiding discrepancies.

The current repair evidence, including the exact CI failure, operator audit,
Ghia-profile comparison, and matched-time grid/time-step studies, is recorded
in [docs/VALIDATION.md](docs/VALIDATION.md).

## License

MIT. See [LICENSE](LICENSE).

**Research principle:** the objective is not just to produce a plausible flow visualization; it is to establish which numerical claims the evidence actually supports.
