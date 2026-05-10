#!/usr/bin/env python3
"""Boxplot: distribution of average amplification ratio across random
F/R pairing trials, per temperature, plus one extra box for "best
temperature per address" optimization.

Reads the rolling sweep.*.csv files produced by sweep_pairings.
Streaming group-by-trial_id parser keeps memory bounded.

Usage: scripts/plot_average_ratio.py [sweep_dir] [out_png]
Defaults: sweep_results/, average_ratio.png
"""
import csv
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sweep_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "sweep_results")
out_png   = Path(sys.argv[2] if len(sys.argv) > 2 else "average_ratio.png")

sweep_files = sorted(sweep_dir.glob("sweep.*.csv")) + \
              sorted(sweep_dir.glob("sweep_job_*/sweep.*.csv"))
if not sweep_files:
    sys.exit(f"no sweep.*.csv under {sweep_dir}")

with open(sweep_files[0]) as fp:
    header = next(csv.reader(fp))
temp_idxs = [i for i, c in enumerate(header) if re.fullmatch(r"\d+C", c)]
temps     = [int(header[i][:-1]) for i in temp_idxs]
n_temps   = len(temps)

per_temp_avg  = [[] for _ in range(n_temps)]
optimized_avg = []

buf = np.empty((4096, n_temps), dtype=np.float64)
n_rows = 0
current_trial = None

def flush():
    global n_rows
    if n_rows == 0:
        return
    arr = buf[:n_rows]
    col_mean = arr.mean(axis=0)
    for t, v in enumerate(col_mean):
        per_temp_avg[t].append(v)
    optimized_avg.append(arr.max(axis=1).mean())
    n_rows = 0

for path in sweep_files:
    with open(path, "rb") as fp:
        text = fp.read()
        os.posix_fadvise(fp.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
    lines = text.split(b"\n")
    for line in lines[1:]:
        if not line:
            continue
        parts = line.split(b",")
        tid = int(parts[0])
        if tid != current_trial:
            flush()
            current_trial = tid
        for j, idx in enumerate(temp_idxs):
            buf[n_rows, j] = float(parts[idx])
        n_rows += 1
flush()

n_trials = len(optimized_avg)

fig, ax = plt.subplots(figsize=(11, 6))
data   = per_temp_avg + [optimized_avg]
labels = [f"{t}°C" for t in temps] + ["best T\nper addr"]
ax.boxplot(data, labels=labels, whis=(0, 100), widths=0.6)
ax.set_yscale("log")
ax.set_ylabel("Average amplification ratio across addresses")
ax.set_title(f"Average amplification ratio across {n_trials} random F/R pairings")
ax.grid(True, axis="y", which="both", alpha=0.25)
ax.axvline(n_temps + 0.5, color="gray", linestyle=":", alpha=0.5)
fig.tight_layout()
fig.savefig(out_png, dpi=120)
plt.close(fig)
print(f"wrote {out_png}: {n_trials} trials, {n_temps + 1} boxes")
