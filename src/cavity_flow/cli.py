"""CLI runner and plotting/output utilities."""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .solver import CavityConfig, CavitySolver


def write_csv(path: Path, headers: list[str], rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f); writer.writerow(headers); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the lid-driven cavity solver")
    parser.add_argument("--re", type=float, default=100.0, dest="reynolds")
    parser.add_argument("--nx", type=int, default=64); parser.add_argument("--ny", type=int, default=64)
    parser.add_argument("--steps", type=int, default=3000); parser.add_argument("--dt", type=float)
    parser.add_argument("--output", type=Path, default=Path("results/run"))
    parser.add_argument("--log-every", type=int, default=100)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    solver = CavitySolver(CavityConfig(nx=args.nx, ny=args.ny, reynolds=args.reynolds, dt=args.dt))
    hist = solver.run(args.steps, args.log_every)
    if hist:
        write_csv(args.output/"history.csv", list(hist[0]), [[r[k] for k in hist[0]] for r in hist])
    uc, vc, _ = solver.cell_centered_fields()
    x = (np.arange(args.nx)+0.5)*solver.dx; y = (np.arange(args.ny)+0.5)*solver.dy
    write_csv(args.output/"centerline_u.csv", ["y", "u"], zip(y, uc[:, args.nx//2]))
    write_csv(args.output/"centerline_v.csv", ["x", "v"], zip(x, vc[args.ny//2, :]))
    X, Y = np.meshgrid(x, y); speed = np.sqrt(uc**2 + vc**2)
    fig, ax = plt.subplots(figsize=(7, 6)); cf = ax.contourf(X, Y, speed, levels=24)
    fig.colorbar(cf, ax=ax, label="speed"); ax.set(xlabel="x", ylabel="y", title=f"Speed, Re={args.reynolds:g}")
    ax.set_aspect("equal"); fig.tight_layout(); fig.savefig(args.output/"velocity.png", dpi=160); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 6)); ax.streamplot(X, Y, uc, vc, density=1.4)
    ax.set(xlabel="x", ylabel="y", title=f"Streamlines, Re={args.reynolds:g}"); ax.set_aspect("equal")
    fig.tight_layout(); fig.savefig(args.output/"streamlines.png", dpi=160); plt.close(fig)
    print(f"iterations={solver.iteration} time={solver.time:.8g}")
    print(f"divergence_linf={solver.divergence_linf():.6e} kinetic_energy={solver.kinetic_energy():.8g}")
    print(f"results={args.output.resolve()}")

if __name__ == "__main__": main()
