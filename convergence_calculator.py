from pathlib import Path
import zipfile, textwrap

out = Path("/mnt/data/cavity_flow_convergence_calculator")
out.mkdir(parents=True, exist_ok=True)

script = r'''#!/usr/bin/env python3
"""
Cavity Flow Grid-Convergence Calculator
=======================================

Estimate observed convergence order from errors measured on successively
refined lid-driven-cavity grids.

CSV input format (one row per grid):
    grid_n,error_u,error_v
    32,0.01122,0.00862
    64,0.00545,0.00339
    128,0.00270,0.00170

The numbers above are EXAMPLES ONLY, not validated solver results.

By default, h is taken proportional to 1/(grid_n - 1), which is appropriate
when grid_n is the number of points spanning the same domain including both
boundaries. If grid_n instead means the number of cells, use --spacing cells,
which takes h proportional to 1/grid_n.

The input errors must all use the same definition, reference data, domain,
norm, and boundary treatment across grids. This script does not run the
Navier-Stokes solver and does not independently verify the reference solution.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class GridRecord:
    grid_n: int
    errors: dict[str, float]


def observed_order(error_coarse: float, error_fine: float,
                   h_coarse: float, h_fine: float) -> float:
    """Return p from E ~ C*h**p for two grids."""
    if error_coarse <= 0 or error_fine <= 0:
        raise ValueError("Errors must be strictly positive to estimate order.")
    if h_coarse <= h_fine or h_fine <= 0:
        raise ValueError("Grids must be ordered coarse-to-fine with h_coarse > h_fine.")
    return math.log(error_coarse / error_fine) / math.log(h_coarse / h_fine)


def grid_spacing_proxy(grid_n: int, spacing: str) -> float:
    if grid_n < 2:
        raise ValueError(f"grid_n must be at least 2; received {grid_n}.")
    if spacing == "points":
        return 1.0 / (grid_n - 1)
    if spacing == "cells":
        return 1.0 / grid_n
    raise ValueError(f"Unknown spacing convention: {spacing}")


def read_csv(path: Path, spacing: str) -> tuple[list[GridRecord], list[str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError("CSV has no header row.")
        fields = [name.strip() for name in reader.fieldnames]
        if "grid_n" not in fields:
            raise ValueError("CSV must contain a 'grid_n' column.")
        error_fields = [name for name in fields if name.startswith("error_")]
        if not error_fields:
            raise ValueError("CSV needs at least one error column, e.g. 'error_u'.")
        rows: list[GridRecord] = []
        for line_number, raw in enumerate(reader, start=2):
            row = {str(k).strip(): (v.strip() if isinstance(v, str) else v)
                   for k, v in raw.items() if k is not None}
            if not row or all(v in ("", None) for v in row.values()):
                continue
            try:
                n_float = float(row["grid_n"])
                if not n_float.is_integer():
                    raise ValueError
                n = int(n_float)
                errors = {key.removeprefix("error_"): float(row[key])
                          for key in error_fields}
            except (ValueError, TypeError, KeyError):
                raise ValueError(f"Invalid numeric value on CSV line {line_number}.")
            if any(not math.isfinite(v) or v <= 0 for v in errors.values()):
                raise ValueError(
                    f"Errors on CSV line {line_number} must be finite and > 0."
                )
            grid_spacing_proxy(n, spacing)
            rows.append(GridRecord(n, errors))
    if len(rows) < 2:
        raise ValueError("Provide at least two grid resolutions.")
    rows.sort(key=lambda item: grid_spacing_proxy(item.grid_n, spacing), reverse=True)
    ns = [r.grid_n for r in rows]
    if len(set(ns)) != len(ns):
        raise ValueError("Each grid_n must be unique.")
    variables = list(rows[0].errors)
    if any(list(r.errors) != variables for r in rows):
        raise ValueError("Every row must contain the same error columns.")
    return rows, variables


def format_report(rows: list[GridRecord], variables: list[str],
                  spacing: str, source: str) -> str:
    h_values = [grid_spacing_proxy(r.grid_n, spacing) for r in rows]
    lines = [
        "CAVITY FLOW GRID-CONVERGENCE REPORT",
        "=" * 38,
        f"Input: {source}",
        f"Grid convention: {spacing} (h proportional to "
        + ("1/(N-1)" if spacing == "points" else "1/N") + ")",
        "",
        "Input errors:",
    ]
    header = ["grid_n", "h_proxy"] + [f"error_{v}" for v in variables]
    lines.append("  ".join(f"{x:>14}" for x in header))
    for r, h in zip(rows, h_values):
        vals = [str(r.grid_n), f"{h:.8g}"] + [f"{r.errors[v]:.8g}" for v in variables]
        lines.append("  ".join(f"{x:>14}" for x in vals))

    lines.extend(["", "Pairwise observed orders (E proportional to h^p):"])
    if len(rows) == 2:
        lines.append("  Only one pair is available; use 3+ grids to assess consistency.")
    for i in range(len(rows) - 1):
        coarse, fine = rows[i], rows[i + 1]
        ratio = h_values[i] / h_values[i + 1]
        lines.append(
            f"  {coarse.grid_n} -> {fine.grid_n} "
            f"(h ratio {ratio:.4g})"
        )
        for variable in variables:
            p = observed_order(coarse.errors[variable], fine.errors[variable],
                               h_values[i], h_values[i + 1])
            lines.append(f"    {variable}: p = {p:.4f}")
            if p <= 0:
                lines.append("      WARNING: error did not decrease under refinement.")
    lines.extend([
        "",
        "INTERPRETATION / LIMITATIONS",
        "- Orders are meaningful only if each error uses the same norm and reference.",
        "- A benchmark/reference error is not the same as a difference between grids.",
        "- Pairwise order can be misleading outside the asymptotic convergence range.",
        "- Non-monotone errors, changing time-step error, or under-converged runs can",
        "  obscure spatial order. Keep solver stopping criteria consistent.",
        "- This report does not prove physical correctness or validate the reference data.",
        "- For a stronger study, use 3+ grids and separately test time-step sensitivity.",
    ])
    return "\n".join(lines)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Estimate observed grid-convergence order for cavity-flow errors."
    )
    parser.add_argument("csv_file", type=Path,
                        help="CSV with grid_n and one or more error_* columns")
    parser.add_argument("--spacing", choices=("points", "cells"), default="points",
                        help="Interpret grid_n as points (default) or cells")
    parser.add_argument("--output", type=Path,
                        help="Optional path to save the text report")
    args = parser.parse_args(argv)

    try:
        rows, variables = read_csv(args.csv_file, args.spacing)
        report = format_report(rows, variables, args.spacing, str(args.csv_file))
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print(report)
    if args.output:
        try:
            args.output.write_text(report + "\n", encoding="utf-8")
        except OSError as exc:
            print(f"Could not save report: {exc}", file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''
(out / "cavity_flow_convergence.py").write_text(script, encoding="utf-8")

sample = """# ILLUSTRATIVE DATA ONLY — replace with measured errors from your solver.
grid_n,error_u,error_v
32,0.01122,0.00862
64,0.00545,0.00339
128,0.00270,0.00170
"""
(out / "example_errors_ILLUSTRATIVE.csv").write_text(sample, encoding="utf-8")

readme = """# Cavity Flow Grid-Convergence Calculator

A small standalone Python CLI tool for estimating observed convergence order
from consistent error measurements across successively refined cavity-flow grids.

## Requirements
- Python 3.10+
- Standard library only; no packages to install.

## CSV format

Use one row per grid and one or more error columns prefixed with `error_`:

```csv
grid_n,error_u,error_v
32,REPLACE_WITH_REAL_U_ERROR,REPLACE_WITH_REAL_V_ERROR
64,REPLACE_WITH_REAL_U_ERROR,REPLACE_WITH_REAL_V_ERROR
128,REPLACE_WITH_REAL_U_ERROR,REPLACE_WITH_REAL_V_ERROR
```

The included example CSV has illustrative numbers only. Do not cite them as solver results.

## Run

```bash
python cavity_flow_convergence.py your_errors.csv
```

Save a report:

```bash
python cavity_flow_convergence.py your_errors.csv --output convergence_report.txt
```

By default, `grid_n` means grid points spanning the domain, so the spacing proxy is
`h = 1/(N-1)`. If your `grid_n` is a number of cells, use:

```bash
python cavity_flow_convergence.py your_errors.csv --spacing cells
```

## Interpretation

For each adjacent grid pair the script computes

`p = ln(E_coarse / E_fine) / ln(h_coarse / h_fine)`.

Errors must be positive and measured consistently using the same reference,
norm, domain, and boundary conditions. Prefer at least three grids. Keep the
solver's convergence tolerance and temporal-error contribution controlled.
The script estimates order; it does not run the flow solver or prove that the
reference data are correct.
"""
(out / "README.md").write_text(readme, encoding="utf-8")

# Run an actual smoke test using the illustrative CSV; record only that the script executes.
import subprocess, sys
run = subprocess.run(
    [sys.executable, str(out / "cavity_flow_convergence.py"),
     str(out / "example_errors_ILLUSTRATIVE.csv")],
    capture_output=True, text=True
)
(out / "smoke_test.txt").write_text(
    f"Exit code: {run.returncode}\n\nSTDOUT:\n{run.stdout}\nSTDERR:\n{run.stderr}",
    encoding="utf-8"
)

zip_path = Path("/mnt/data/cavity_flow_convergence_calculator.zip")
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for p in out.iterdir():
        zf.write(p, arcname=p.name)

print(f"Created: {out / 'cavity_flow_convergence.py'}")
print(f"Created: {zip_path}")
print(f"Smoke test exit code: {run.returncode}")
