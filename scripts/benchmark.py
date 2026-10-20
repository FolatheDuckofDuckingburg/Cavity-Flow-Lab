"""Repeatable fixed-step runtime and numerical-diagnostic benchmark."""
from __future__ import annotations

import argparse
import json
import platform
import statistics
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from cavity_flow import CavityConfig, CavitySolver


def run_once(n: int, re: float, steps: int, tolerance: float) -> dict:
    solver = CavitySolver(CavityConfig(nx=n, ny=n, reynolds=re))
    start = time.perf_counter()
    history = solver.run(steps=steps, log_every=steps)
    elapsed = time.perf_counter() - start
    last = history[-1]
    div = solver.divergence_linf()
    return {
        "nx": n, "ny": n, "reynolds": re, "requested_steps": steps,
        "completed_steps": solver.iteration, "simulated_time": solver.time,
        "wall_seconds": elapsed,
        "steps_per_second": solver.iteration / elapsed if elapsed > 0 else None,
        "final_divergence_linf": div, "divergence_tolerance": tolerance,
        "divergence_gate": "PASS" if div <= tolerance else "FAIL",
        "pre_projection_divergence_linf": last["divergence_pre_projection_linf"],
        "pressure_residual_linf": last["pressure_residual_linf"],
        "pressure_converged": last["pressure_converged"],
        "mean_pressure_iterations": last["pressure_iterations"],
        "final_kinetic_energy": solver.kinetic_energy(),
        "all_fields_finite": bool(np.isfinite(solver.u).all()
                                  and np.isfinite(solver.v).all()
                                  and np.isfinite(solver.p).all()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grids", nargs="+", type=int, default=[32, 64])
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--re", type=float, default=100.0)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--divergence-tolerance", type=float, default=1e-6)
    parser.add_argument("--output", type=Path, default=Path("results/benchmark.json"))
    args = parser.parse_args()
    if args.steps < 1 or args.repeats < 1 or any(n < 4 for n in args.grids):
        parser.error("steps/repeats must be positive and grid sizes must be >= 4")
    if args.re <= 0 or args.divergence_tolerance <= 0:
        parser.error("Reynolds number and divergence tolerance must be positive")

    runs = []
    for n in args.grids:
        repeats = []
        for repeat in range(1, args.repeats + 1):
            row = run_once(n, args.re, args.steps, args.divergence_tolerance)
            row["repeat"] = repeat
            runs.append(row)
            repeats.append(row["wall_seconds"])
            print(f"{n}x{n} run {repeat}: {row['wall_seconds']:.3f}s; "
                  f"div={row['final_divergence_linf']:.3e}; {row['divergence_gate']}")
        print(f"{n}x{n} median runtime: {statistics.median(repeats):.3f}s")
    payload = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "numpy_version": np.__version__,
        "method": "fixed step count; timings are environment-specific",
        "parameters": {"grids": args.grids, "steps": args.steps, "re": args.re,
                       "repeats": args.repeats, "divergence_tolerance": args.divergence_tolerance},
        "runs": runs,
        "overall_status": "PASS" if all(
            x["divergence_gate"] == "PASS" and x["all_fields_finite"] and x["pressure_converged"]
            for x in runs
        ) else "FAIL",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Diagnostic status: {payload['overall_status']}; report: {args.output}")
    return 0 if payload["overall_status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
