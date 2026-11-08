"""Matched-physical-time grid sensitivity study for the lid-driven cavity solver.

This is a finite-time sensitivity diagnostic, not a formal grid-convergence
study and not proof of CFD validation. Every grid is advanced to the same
requested nondimensional time; the final step is shortened when necessary.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import replace
from pathlib import Path
import platform
import sys

import numpy as np

from cavity_flow import CavityConfig, CavitySolver


def run_to_time(n: int, reynolds: float, target_time: float, dt_cap: float,
                pressure_tolerance: float, pressure_max_iterations: int,
                max_steps: int) -> dict:
    config = CavityConfig(
        nx=n,
        ny=n,
        reynolds=reynolds,
        dt=dt_cap,
        pressure_tolerance=pressure_tolerance,
        pressure_max_iterations=pressure_max_iterations,
    )
    solver = CavitySolver(config)
    last: dict = {}
    while solver.time < target_time:
        if solver.iteration >= max_steps:
            raise RuntimeError(
                f"{n}x{n} did not reach T={target_time:g} within {max_steps} steps"
            )
        remaining = target_time - solver.time
        # Keep the requested cap while shortening only the final step if needed.
        solver.cfg = replace(solver.cfg, dt=min(dt_cap, remaining))
        last = solver.step()

    uc, vc, _ = solver.cell_centered_fields()
    y = (np.arange(n) + 0.5) / n
    x = (np.arange(n) + 0.5) / n
    return {
        "grid": n,
        "reynolds": reynolds,
        "target_time": target_time,
        "simulated_time": solver.time,
        "steps": solver.iteration,
        "dt_min": min(float(row["dt"]) for row in [last]),
        "dt_last": float(last["dt"]),
        "divergence_linf": solver.divergence_linf(),
        "pressure_residual_linf": float(last["pressure_residual_linf"]),
        "pressure_converged": bool(last["pressure_converged"]),
        "kinetic_energy": solver.kinetic_energy(),
        "u_y": y.tolist(),
        "u_centerline": uc[:, n // 2].tolist(),
        "v_x": x.tolist(),
        "v_centerline": vc[n // 2, :].tolist(),
        "all_fields_finite": bool(np.isfinite(uc).all() and np.isfinite(vc).all()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grids", nargs="+", type=int, default=[32, 64, 128])
    parser.add_argument("--re", type=float, default=100.0)
    parser.add_argument("--time", type=float, default=0.25,
                        help="common nondimensional physical time")
    parser.add_argument("--dt", type=float, default=5e-4,
                        help="maximum time step; stability estimate may reduce it")
    parser.add_argument("--pressure-tolerance", type=float, default=1e-8)
    parser.add_argument("--pressure-max-iterations", type=int, default=10000)
    parser.add_argument("--max-steps", type=int, default=20000)
    parser.add_argument("--output", type=Path,
                        default=Path("results/matched-time-grid-study"))
    args = parser.parse_args()

    if (not args.grids or any(n < 4 for n in args.grids)
            or args.re <= 0 or args.time <= 0 or args.dt <= 0
            or args.max_steps < 1):
        parser.error("grids must be >= 4; Re, time, dt, and max-steps must be positive")

    runs = []
    for n in args.grids:
        print(f"Running Re={args.re:g}, grid={n}x{n}, target T={args.time:g} ...")
        run = run_to_time(
            n, args.re, args.time, args.dt,
            args.pressure_tolerance, args.pressure_max_iterations,
            args.max_steps,
        )
        runs.append(run)
        print(
            f"  steps={run['steps']} T={run['simulated_time']:.12g} "
            f"div_inf={run['divergence_linf']:.3e} "
            f"p_res={run['pressure_residual_linf']:.3e} "
            f"p_ok={run['pressure_converged']}"
        )

    # Compare all profiles on a common interior coordinate set. The finest
    # requested grid is only a comparison reference, not assumed exact truth.
    reference = max(runs, key=lambda item: item["grid"])
    common = np.linspace(0.02, 0.98, 241)
    ref_u = np.interp(common, reference["u_y"], reference["u_centerline"])
    ref_v = np.interp(common, reference["v_x"], reference["v_centerline"])
    comparisons = []
    for run in runs:
        u = np.interp(common, run["u_y"], run["u_centerline"])
        v = np.interp(common, run["v_x"], run["v_centerline"])
        comparisons.append({
            "grid": run["grid"],
            "comparison_reference_grid": reference["grid"],
            "u_centerline_rmse_vs_finest_requested_grid": float(np.sqrt(np.mean((u-ref_u)**2))),
            "u_centerline_max_abs_difference_vs_finest_requested_grid": float(np.max(np.abs(u-ref_u))),
            "v_centerline_rmse_vs_finest_requested_grid": float(np.sqrt(np.mean((v-ref_v)**2))),
            "v_centerline_max_abs_difference_vs_finest_requested_grid": float(np.max(np.abs(v-ref_v))),
            "kinetic_energy_difference_vs_finest_requested_grid": float(
                abs(run["kinetic_energy"] - reference["kinetic_energy"])
            ),
        })

    args.output.mkdir(parents=True, exist_ok=True)
    # Keep the machine-readable file compact; profile arrays are saved separately.
    payload = {
        "schema_version": 1,
        "study": "matched physical time grid sensitivity",
        "interpretation_limit": (
            "Finite-time diagnostic only. The finest requested grid is a "
            "comparison reference, not exact truth; this is not formal order "
            "of accuracy or CFD validation."
        ),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "parameters": {
            "grids": args.grids,
            "reynolds": args.re,
            "target_time": args.time,
            "dt_cap": args.dt,
            "pressure_tolerance": args.pressure_tolerance,
            "pressure_max_iterations": args.pressure_max_iterations,
        },
        "runs": [
            {k: v for k, v in run.items()
             if k not in {"u_y", "u_centerline", "v_x", "v_centerline"}}
            for run in runs
        ],
        "comparisons": comparisons,
    }
    (args.output / "summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    with (args.output / "centerlines.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["grid", "coordinate", "u_vertical_centerline", "v_horizontal_centerline"])
        for run in runs:
            for coord, uval in zip(run["u_y"], run["u_centerline"]):
                writer.writerow([run["grid"], coord, uval, ""])
            for coord, vval in zip(run["v_x"], run["v_centerline"]):
                writer.writerow([run["grid"], coord, "", vval])

    print(f"Saved results to {args.output}")
    return 0 if all(
        run["all_fields_finite"] and run["pressure_converged"]
        for run in runs
    ) else 2


if __name__ == "__main__":
    raise SystemExit(main())
