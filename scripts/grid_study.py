"""Run a small grid-refinement comparison and write machine-readable results.

This script compares fixed-step snapshots, not necessarily steady-state solutions.
Use the recorded simulated time to avoid mistaking different physical times for
grid-convergence evidence.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from cavity_flow import CavityConfig, CavitySolver


def run_grid(n: int, re: float, steps: int) -> dict:
    solver = CavitySolver(CavityConfig(nx=n, ny=n, reynolds=re))
    history = solver.run(steps=steps, log_every=steps)
    uc, vc, _ = solver.cell_centered_fields()
    return {
        "grid": n,
        "reynolds": re,
        "steps": solver.iteration,
        "simulated_time": solver.time,
        "dt_last": history[-1]["dt"],
        "divergence_linf": solver.divergence_linf(),
        "pressure_residual_linf": history[-1]["pressure_residual_linf"],
        "pressure_converged": history[-1]["pressure_converged"],
        "kinetic_energy": solver.kinetic_energy(),
        "u_centerline": uc[:, n // 2].tolist(),
        "v_centerline": vc[n // 2, :].tolist(),
        "all_fields_finite": bool(np.isfinite(uc).all() and np.isfinite(vc).all()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grids", nargs="+", type=int, default=[32, 64, 128])
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--re", type=float, default=100.0)
    parser.add_argument("--output", type=Path, default=Path("results/grid-study.json"))
    args = parser.parse_args()
    if args.steps < 1 or any(n < 4 for n in args.grids) or args.re <= 0:
        parser.error("steps and grid sizes must be positive (grid sizes >= 4); Re must be positive")
    rows = [run_grid(n, args.re, args.steps) for n in args.grids]
    payload = {
        "schema_version": 1,
        "method": "same requested step count per grid; compare simulated times before interpreting differences",
        "parameters": {"grids": args.grids, "steps": args.steps, "re": args.re},
        "runs": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}; simulated times: {[round(x['simulated_time'], 6) for x in rows]}")
    return 0 if all(x["all_fields_finite"] and x["pressure_converged"] for x in rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
