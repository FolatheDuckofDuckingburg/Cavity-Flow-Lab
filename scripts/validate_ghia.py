#!/usr/bin/env python3
"""Validate the cavity solver against Ghia, Ghia & Shin (1982) at Re = 100.

Reference provenance
--------------------
Table I  (u along the vertical centreline x = 0.5)
Table II (v along the horizontal centreline y = 0.5)
Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for
incompressible flow using the Navier-Stokes equations and a multigrid method.
Journal of Computational Physics, 48(3), 387-411.
doi:10.1016/0021-9991(82)90058-4

Only the Re = 100 column is stored below, so that a different Reynolds number
cannot be compared by accident. The values match the transcriptions the
project's own docs/VALIDATION.md cites.

Sampling convention (explicit; no extrapolation, no point "matching")
---------------------------------------------------------------------
* solver cell-centre coordinates are x_i = (i+0.5)/nx, y_j = (j+0.5)/ny
* solver u at x = 0.5 is linearly interpolated from the two adjacent columns
* Ghia's tabulated profile is linearly interpolated onto the solver grid
* the wall points y = 0/1 and x = 0/1 are excluded: cell centres never reach
  the wall, and extrapolating to them would fabricate a comparison point
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cavity_flow.solver import CavityConfig, CavitySolver  # noqa: E402

GHIA_U_RE100 = [
    (1.0000, 1.00000), (0.9766, 0.84123), (0.9688, 0.78871), (0.9609, 0.73722),
    (0.9531, 0.68717), (0.8516, 0.23151), (0.7344, 0.00332), (0.6172, -0.13641),
    (0.5000, -0.20581), (0.4531, -0.21090), (0.2813, -0.15662),
    (0.1719, -0.10150), (0.1016, -0.06434), (0.0703, -0.04775),
    (0.0625, -0.04192), (0.0547, -0.03717), (0.0000, 0.00000),
]
GHIA_V_RE100 = [
    (1.00000, 0.00000), (0.9688, -0.05906), (0.9609, -0.07391),
    (0.9531, -0.08864), (0.9453, -0.10313), (0.9063, -0.16914),
    (0.8594, -0.22445), (0.8047, -0.24533), (0.5000, 0.05454),
    (0.2344, 0.17527), (0.2266, 0.17507), (0.1563, 0.16077),
    (0.0938, 0.12317), (0.0781, 0.10890), (0.0703, 0.10091),
    (0.0625, 0.09233), (0.00000, 0.00000),
]


def ghia_values_at(table, coords):
    """Linear interpolation of a Ghia table onto coords, no extrapolation."""
    src = np.array([p for p, _ in table], dtype=float)
    val = np.array([v for _, v in table], dtype=float)
    out = np.full(coords.shape, np.nan)
    inside = (coords >= src.min()) & (coords <= src.max())
    out[inside] = np.interp(coords[inside], src[::-1], val[::-1])
    return out


def profile_errors(solver_coord, solver_val, table, label):
    ref = ghia_values_at(table, solver_coord)
    ok = np.isfinite(ref)
    err = solver_val[ok] - ref[ok]
    n = int(ok.sum())
    return {
        "profile": label,
        "n_compared_points": n,
        "n_solver_points": int(solver_coord.size),
        "rmse": float(np.sqrt(np.mean(err**2))),
        "max_abs_error": float(np.max(np.abs(err))),
        "mean_abs_error": float(np.mean(np.abs(err))),
        "argmax_abs_error_coord": float(solver_coord[ok][np.argmax(np.abs(err))]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=64)
    ap.add_argument("--ny", type=int, default=64)
    ap.add_argument("--reynolds", type=float, default=100.0)
    ap.add_argument("--steps", type=int, default=4800)
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--max-pressure-iterations", type=int, default=3000)
    ap.add_argument("--output", default="ghia_re100.json")
    args = ap.parse_args()

    cfg = CavityConfig(
        nx=args.nx, ny=args.ny, reynolds=args.reynolds,
        pressure_tolerance=args.tol,
        pressure_max_iterations=args.max_pressure_iterations,
    )
    solver = CavitySolver(cfg)
    dt0 = solver.stable_dt()

    t0 = time.perf_counter()
    last = None
    update_trace = []
    for k in range(args.steps):
        last = solver.step()
        if (k + 1) % 500 == 0 or k == args.steps - 1:
            update_trace.append({
                "iteration": int(last["iteration"]), "time": last["time"],
                "max_du_dv": max(last["du_linf"], last["dv_linf"]),
            })
            print(f"  step {k+1}/{args.steps} t={last['time']:.4f} "
                  f"upd={max(last['du_linf'], last['dv_linf']):.3e} "
                  f"div={last['divergence_post_projection_linf']:.3e} "
                  f"piter={int(last['pressure_iterations'])}", flush=True)
    elapsed = time.perf_counter() - t0

    uc, vc, _ = solver.cell_centered_fields()
    xc = (np.arange(cfg.nx) + 0.5) * solver.dx
    yc = (np.arange(cfg.ny) + 0.5) * solver.dy

    # u along vertical centreline x = 0.5
    u_center = np.array([
        np.interp(0.5, xc, uc[j, :]) for j in range(cfg.ny)
    ])
    # v along horizontal centreline y = 0.5
    v_center = np.array([
        np.interp(0.5, yc, vc[:, i]) for i in range(cfg.nx)
    ])

    u_stats = profile_errors(yc, u_center, GHIA_U_RE100, "u along x=0.5")
    v_stats = profile_errors(xc, v_center, GHIA_V_RE100, "v along y=0.5")

    div = solver.divergence()
    report = {
        "config": {
            "nx": cfg.nx, "ny": cfg.ny, "reynolds": cfg.reynolds,
            "kinematic_viscosity": solver.nu, "dx": solver.dx, "dy": solver.dy,
            "pressure_tolerance": cfg.pressure_tolerance,
            "pressure_max_iterations": cfg.pressure_max_iterations,
            "requested_steps": args.steps,
        },
        "time_stepping": {
            "stable_dt_at_start": dt0,
            "final_dt": last["dt"],
            "simulated_time": last["time"],
            "note": ("dt is chosen adaptively by stable_dt(); a fixed dt is only "
                     "used if CavityConfig(dt=...) is set."),
        },
        "convergence": {
            "final_max_du_dv": max(last["du_linf"], last["dv_linf"]),
            "update_trace": update_trace,
            "diagnostic_only": "A small velocity update is a numerical heuristic, not proof of steady state.",
        },
        "projection": {
            "pre_projection_divergence_linf": last["divergence_pre_projection_linf"],
            "post_projection_divergence_linf": last["divergence_post_projection_linf"],
            "pressure_iterations_final": int(last["pressure_iterations"]),
            "pressure_residual_linf_final": last["pressure_residual_linf"],
            "pressure_converged": bool(last["pressure_converged"]),
        },
        "mass_conservation": {
            "global_divergence_sum": float(div.sum()),
            "global_divergence_sum_abs_max_cell_volume": float(div.sum() * solver.dx * solver.dy),
            "divergence_linf": float(np.max(np.abs(div))),
            "divergence_l2": float(np.sqrt(np.mean(div**2))),
        },
        "ghia_comparison": {
            "reference": "Ghia, Ghia & Shin (1982), Tables I and II, Re=100 column",
            "doi": "10.1016/0021-9991(82)90058-4",
            "interpolation": "linear, onto solver cell-centre coordinates; no extrapolation",
            "u_vertical_centerline": u_stats,
            "v_horizontal_centerline": v_stats,
        },
        "run": {
            "wall_clock_seconds": elapsed,
            "ms_per_step": 1000.0 * elapsed / max(1, args.steps),
            "numpy_version": np.__version__,
        },
        "limitations": [
            "Second-order centred advection on a collocated-in-face sense; no turbulence model.",
            "The lid corner convention overwrites u's top-row corner values.",
            "Agreement with a 129x129 reference solution from a different discretization is not proof of correctness.",
        ],
    }
    Path(args.output).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in
                      ("time_stepping", "projection", "mass_conservation")}, indent=2))
    print(json.dumps(report["ghia_comparison"], indent=2))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
