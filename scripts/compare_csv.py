#!/usr/bin/env python3
"""Compare two regression CSVs cell-by-cell, report aggregate statistics.

Usage: compare_csv.py <baseline.csv> <test.csv>

Prints a one-line summary plus details about the largest disagreements.
Lines that aren't pure CSV numerics (e.g. "Simulating address N" headers
written by test_eq) are skipped and counted separately.
"""
import sys
from pathlib import Path


def parse_row(line):
    """Return list of floats, or None if the row isn't all numeric."""
    parts = line.strip().split(",")
    if not parts:
        return None
    out = []
    for p in parts:
        p = p.strip()
        if p == "":
            return None
        try:
            out.append(float(p))
        except ValueError:
            return None
    return out


def relerr(a, b):
    """Symmetric relative error; returns 0 when both are 0."""
    denom = max(abs(a), abs(b))
    if denom == 0:
        return 0.0
    return abs(a - b) / denom


def compare(baseline_path, test_path):
    base_rows = []
    test_rows = []
    skipped_base = 0
    skipped_test = 0

    with open(baseline_path) as f:
        for line in f:
            if not line.strip():
                continue
            row = parse_row(line)
            if row is None:
                skipped_base += 1
                continue
            base_rows.append(row)

    with open(test_path) as f:
        for line in f:
            if not line.strip():
                continue
            row = parse_row(line)
            if row is None:
                skipped_test += 1
                continue
            test_rows.append(row)

    if len(base_rows) != len(test_rows):
        print(f"  WARN: row counts differ "
              f"(baseline={len(base_rows)}, test={len(test_rows)})")

    n_rows = min(len(base_rows), len(test_rows))
    if n_rows == 0:
        print("  no comparable rows!")
        return

    total_cells = 0
    exact = 0
    close_1e9 = 0   # rel err < 1e-9 (matches 9-digit print precision)
    close_1e6 = 0   # rel err < 1e-6
    close_1e3 = 0   # rel err < 1e-3
    max_rel = 0.0
    max_loc = (0, 0)
    max_vals = (0.0, 0.0)

    for r in range(n_rows):
        rb = base_rows[r]
        rt = test_rows[r]
        if len(rb) != len(rt):
            print(f"  WARN: row {r} column count differs "
                  f"({len(rb)} vs {len(rt)})")
            continue
        for c, (a, b) in enumerate(zip(rb, rt)):
            total_cells += 1
            if a == b:
                exact += 1
                close_1e9 += 1
                close_1e6 += 1
                close_1e3 += 1
                continue
            re = relerr(a, b)
            if re < 1e-9:
                close_1e9 += 1
            if re < 1e-6:
                close_1e6 += 1
            if re < 1e-3:
                close_1e3 += 1
            if re > max_rel:
                max_rel = re
                max_loc = (r, c)
                max_vals = (a, b)

    print(f"  rows: {n_rows}, cells: {total_cells}")
    print(f"  exact match:    {exact}/{total_cells} "
          f"({100.0*exact/total_cells:.2f}%)")
    print(f"  rel < 1e-9:     {close_1e9}/{total_cells} "
          f"({100.0*close_1e9/total_cells:.2f}%)")
    print(f"  rel < 1e-6:     {close_1e6}/{total_cells} "
          f"({100.0*close_1e6/total_cells:.2f}%)")
    print(f"  rel < 1e-3:     {close_1e3}/{total_cells} "
          f"({100.0*close_1e3/total_cells:.2f}%)")
    print(f"  max rel error:  {max_rel:.3e}  "
          f"at row {max_loc[0]}, col {max_loc[1]}  "
          f"(baseline={max_vals[0]:.6e}, test={max_vals[1]:.6e})")
    if skipped_base or skipped_test:
        print(f"  skipped non-numeric lines: "
              f"baseline={skipped_base}, test={skipped_test}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    print(f"  {Path(sys.argv[1]).name} vs {Path(sys.argv[2]).name}")
    compare(sys.argv[1], sys.argv[2])
