"""Command-line runner and CSV/plot output for Cavity Flow Lab."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import platform
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .solver import CavityConfig, CavitySolver


def write_csv(path: Path, headers: list[str], rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the experimental lid-driven cavity solver")
    parser.add_argument("--re", type=float, default=100.0, dest="reynolds")
    parser.add_argument("--nx", type=int, default=64)
    parser.add_argument("--ny", type=int, default=64)
    parser.add_argument("--steps", type=int, default=3000)
    parser.add_argument("--dt", type=float, default=None)
    parser.add_argument("--cfl", type=float, default=0.25)
    parser.add_argument("--pressure-tolerance", type=float, default=1e-8)
    parser.add_argument("--pressure-max-iterations", type=int, default=2000)
    parser.add_argument("--log-every", type=int, default=100)
    parser.add_argument("--until-steady", action="store_true",
                        help="stop when velocity updates remain below tolerance")
    parser.add_argument("--steady-tolerance", type=float, default=1e-7)
    parser.add_argument("--min-steps", type=int, default=100)
    parser.add_argument("--output", type=Path, default=Path("results/run"))
    args = parser.parse_args()
    if args.steps < 1 or args.log_every < 1:
        parser.error("--steps and --log-every must be >= 1")

    args.output.mkdir(parents=True, exist_ok=True)
    cfg = CavityConfig(
        nx=args.nx, ny=args.ny, reynolds=args.reynolds, dt=args.dt, cfl=args.cfl,
        pressure_tolerance=args.pressure_tolerance,
        pressure_max_iterations=args.pressure_max_iterations,
    )
    solver = CavitySolver(cfg)
    if args.until_steady:
        history = solver.run_until_steady(args.steps, args.steady_tolerance, args.min_steps)
    else:
        history = solver.run(args.steps, args.log_every)

    headers = list(history[0].keys()) if history else []
    if history:
        write_csv(args.output / "history.csv", headers, ([row[k] for k in headers] for row in history))
    uc, vc, pressure = solver.cell_centered_fields()
    x = (np.arange(args.nx) + 0.5) * solver.dx
    y = (np.arange(args.ny) + 0.5) * solver.dy
    write_csv(args.output / "centerline_u.csv", ["y", "u"], zip(y, uc[:, args.nx // 2]))
    write_csv(args.output / "centerline_v.csv", ["x", "v"], zip(x, vc[args.ny // 2, :]))

    X, Y = np.meshgrid(x, y)
    speed = np.sqrt(uc**2 + vc**2)
    fig, ax = plt.subplots(figsize=(7, 6))
    contour = ax.contourf(X, Y, speed, levels=24)
    fig.colorbar(contour, ax=ax, label="speed")
    ax.set(xlabel="x", ylabel="y", title=f"Speed, Re={args.reynolds:g}")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(args.output / "velocity.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.streamplot(X, Y, uc, vc, density=1.4)
    ax.set(xlabel="x", ylabel="y", title=f"Streamlines, Re={args.reynolds:g}")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(args.output / "streamlines.png", dpi=160)
    plt.close(fig)

    final = history[-1] if history else {}
    metadata = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "parameters": {
            "nx": args.nx, "ny": args.ny, "reynolds": args.reynolds,
            "dt_requested": args.dt, "cfl": args.cfl, "steps_requested": args.steps,
            "until_steady": args.until_steady, "steady_tolerance": args.steady_tolerance,
            "pressure_tolerance": args.pressure_tolerance,
        },
        "results": {
            "iterations": solver.iteration, "simulated_time": solver.time,
            "divergence_linf": solver.divergence_linf(),
            "kinetic_energy": solver.kinetic_energy(),
            "last_pressure_residual_linf": final.get("pressure_residual_linf"),
            "last_pressure_converged": final.get("pressure_converged"),
        },
        "limitations": "Experimental 2D solver; not independently validated CFD.",
    }
    (args.output / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"iterations={solver.iteration} time={solver.time:.8g}")
    print(f"divergence_linf={solver.divergence_linf():.6e} kinetic_energy={solver.kinetic_energy():.8g}")
    print(f"results={args.output.resolve()}")


if __name__ == "__main__":
    main()
