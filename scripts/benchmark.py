"""Fixed-step runtime and numerical diagnostic benchmark.

Example: python scripts/benchmark.py --grids 32 64 --steps 100 --repeats 3 --output results/benchmark.json
Timings are machine-specific. tracemalloc does not necessarily capture native NumPy allocations.
The divergence gate is a diagnostic, not complete CFD validation.
"""
import argparse
import json
import platform
import statistics
import sys
import time
import tracemalloc
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from cavity_flow import CavityConfig, CavitySolver


def run_once(n, re, steps, tolerance):
    solver = CavitySolver(CavityConfig(nx=n, ny=n, reynolds=re))
    tracemalloc.start()
    start = time.perf_counter()
    history = solver.run(steps=steps, log_every=1)
    elapsed = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    div = solver.divergence_linf()
    return {
        "nx": n, "ny": n, "reynolds": re, "requested_steps": steps,
        "completed_steps": solver.iteration, "simulated_time": solver.time,
        "wall_seconds": elapsed, "steps_per_second": solver.iteration / elapsed,
        "python_traced_peak_bytes": peak, "final_divergence_linf": div,
        "divergence_tolerance": tolerance,
        "divergence_gate": "PASS" if div <= tolerance else "FAIL",
        "final_pressure_residual_linf": history[-1]["pressure_residual_linf"],
        "mean_pressure_iterations": statistics.mean(x["pressure_iterations"] for x in history),
        "final_kinetic_energy": solver.kinetic_energy(),
        "all_fields_finite": bool(np.isfinite(solver.u).all() and np.isfinite(solver.v).all()
                                  and np.isfinite(solver.p).all()),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--grids", nargs="+", type=int, default=[32, 64])
    p.add_argument("--steps", type=int, default=100)
    p.add_argument("--re", type=float, default=100.0)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--divergence-tolerance", type=float, default=1e-6)
    p.add_argument("--output", type=Path, default=Path("results/benchmark.json"))
    a = p.parse_args()
    if a.steps < 1 or a.repeats < 1 or any(n < 4 for n in a.grids):
        p.error("steps/repeats must be positive and grid sizes must be >= 4")
    if a.re <= 0 or a.divergence_tolerance <= 0:
        p.error("Reynolds number and divergence tolerance must be positive")
    runs = []
    for n in a.grids:
        for repeat in range(1, a.repeats + 1):
            row = run_once(n, a.re, a.steps, a.divergence_tolerance)
            row["repeat"] = repeat
            runs.append(row)
            print(f"{n}x{n} repeat {repeat}: {row['wall_seconds']:.3f}s; "
                  f"div={row['final_divergence_linf']:.3e}; {row['divergence_gate']}")
    payload = {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version, "platform": platform.platform(),
        "numpy_version": np.__version__,
        "timing_note": "fixed number of steps; compare timings only on the same environment",
        "parameters": vars(a) | {"output": str(a.output)},
        "runs": runs,
        "overall_status": "PASS" if all(x["divergence_gate"] == "PASS" and x["all_fields_finite"]
                                         for x in runs) else "FAIL",
    }
    payload["parameters"]["output"] = str(a.output)
    payload["parameters"]["grids"] = a.grids
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"Overall diagnostic status: {payload['overall_status']}; report={a.output}")
    return 0 if payload["overall_status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
