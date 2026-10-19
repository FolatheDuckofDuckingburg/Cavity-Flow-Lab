# Cavity Flow Lab

**A Python laboratory for numerical experiments with the incompressible Navier–Stokes equations.**

Cavity Flow Lab is an open-source computational fluid dynamics (CFD) project for simulating two-dimensional lid-driven cavity flow. It provides a foundation for investigating numerical discretization, pressure–velocity coupling, convergence, and grid-dependent error.

The goal is to make numerical methods inspectable, experiments reproducible, and scientific limitations explicit.

> **Research status:** This is an experimental solver. The pressure-projection implementation and boundary treatment require further verification before its results can be considered validated CFD solutions.

## Governing Equations

The incompressible Navier–Stokes equations describe how viscous fluids move under the influence of inertia, pressure, and viscosity.

The solver addresses two fundamental requirements:

- **Momentum conservation:** Fluid acceleration is governed by convection, pressure forces, and viscous diffusion.
- **Mass conservation:** For an incompressible fluid, the velocity field must satisfy the zero-divergence condition.

These equations form the mathematical foundation of the simulation.

## Lid-Driven Cavity Flow

The benchmark consists of a square cavity containing a viscous fluid. The top wall moves horizontally at a prescribed velocity, while the other walls remain stationary.

The moving lid drives the fluid inside the cavity, producing circulation and characteristic velocity profiles.

The Reynolds number describes the relative importance of inertial forces and viscous forces. The initial target case is Reynolds number 100.

## Features

- Two-dimensional incompressible flow simulation.
- Staggered-grid velocity representation.
- Explicit advection and diffusion predictor.
- Pressure-Poisson equation solved using Jacobi iterations.
- Velocity correction through a pressure-projection step.
- Configurable Reynolds number and grid dimensions.
- Numerical diagnostics and flow-field output.
- Grid-study script for comparing resolutions.
- Automated tests and GitHub Actions workflow.

## Project Structure

```text
Cavity-Flow-Lab/
├── .github/
│   └── workflows/
│       └── tests.yml
├── docs/
│   └── VALIDATION.md
├── scripts/
│   └── grid_study.py
├── src/
│   └── cavity_flow/
│       ├── __init__.py
│       ├── solver.py
│       └── cli.py
├── tests/
│   └── test_solver.py
├── .gitignore
├── LICENSE
├── pyproject.toml
└── README.md
```

## Requirements

- Python 3.10 or later.
- NumPy.
- Matplotlib.
- Development dependencies defined in `pyproject.toml`.

## Installation

Clone the repository:

```bash
git clone https://github.com/FolatheDuckofDuckingburg/Cavity-Flow-Lab.git
cd Cavity-Flow-Lab
```

Create a virtual environment:

```bash
python -m venv .venv
```

Activate the environment.

**Windows PowerShell:**

```powershell
.venv\Scripts\Activate.ps1
```

**Linux or macOS:**

```bash
source .venv/bin/activate
```

Install the project and development dependencies:

```bash
python -m pip install -e ".[dev]"
```

## Running a Simulation

Run the solver at Reynolds number 100 on a 64-by-64 grid:

```bash
cavity-flow --re 100 --nx 64 --ny 64 --steps 3000 --output results/re100_64
```

Display the available command-line options:

```bash
cavity-flow --help
```

Store results from different configurations separately so that experiments can be compared reproducibly.

## Running Tests

Execute the automated test suite:

```bash
pytest
```

The tests check basic solver behavior, array consistency, and selected numerical properties.

Passing the tests does not establish that the numerical solution is physically accurate or fully validated.

## Numerical Validation

A credible CFD result requires more than a plausible flow visualization or a small pressure-solver residual.

### 1. Discrete Consistency

Verify that the discrete divergence operator, pressure-gradient operator, pressure-Poisson equation, and boundary conditions are mathematically consistent.

### 2. Mass Conservation

Measure the remaining divergence of the velocity field after the pressure-correction step. Document the diagnostic, numerical tolerance, and results for each experiment.

### 3. Grid Convergence

Compare solutions on successively refined meshes. Quantify differences between grids rather than assuming that increasing resolution automatically improves accuracy.

### 4. Benchmark Comparison

Compare the computed centerline velocity profiles against established lid-driven cavity reference data.

A key reference is:

Ghia, U. K. N., Ghia, K. N., and Shin, C. T. (1982). *High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method.* Journal of Computational Physics, 48(3), 387–411.

DOI: https://doi.org/10.1016/0021-9991(82)90058-4

### 5. Reproducibility

Document the parameters, grid dimensions, time step, stopping criteria, residuals, and error metrics for each reported experiment.

See [`docs/VALIDATION.md`](docs/VALIDATION.md) for the validation checklist.

## Scientific Limitations

This implementation is an evolving numerical experiment, not a validated production CFD package.

The pressure-projection step and wall-boundary treatment require further numerical verification. A small pressure-solver residual does not, by itself, guarantee that the corrected velocity field satisfies mass conservation or that the computed flow is accurate.

Benchmark agreement, mass conservation, and grid convergence must be demonstrated independently before quantitative scientific conclusions are drawn.

## Project Goals

- Develop a transparent and reproducible CFD codebase.
- Investigate pressure–velocity coupling in incompressible flow.
- Evaluate numerical stability and convergence.
- Quantify errors against trusted reference solutions.
- Document numerical limitations, failures, and improvements.

## Contributing

Contributions involving numerical consistency, boundary-condition handling, convergence analysis, testing, documentation, and benchmark comparisons are welcome.

Numerical changes should include an explanation of the method and, where possible, tests or quantitative evidence supporting the change.

## License

Distributed under the MIT License. See [`LICENSE`](LICENSE) for details.

## Research Principle

**The objective is not merely to produce a flow visualization. It is to establish whether the computed solution is numerically credible.**
