#!/usr/bin/env python3
"""Boxplot: distribution of amplification-ratio standard deviation
across random F/R pairing trials, per temperature, plus one extra
box for "best temperature per address" optimization.

Memory-bounded — see plot_worst_address.py for the parser-buffer
rationale. Total RSS stays around 150 MB regardless of trial count.

Usage: scripts/plot_std_ratio.py [sweep_dir] [out_png]
Defaults: sweep_results/, std_ratio.png
"""
import csv
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sweep_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "sweep_results")
out_png   = Path(sys.argv[2] if len(sys.argv) > 2 else "std_ratio.png")

trial_re = re.compile(r"^trial_\d+\.csv$")
trial_files = sorted(p for p in sweep_dir.glob("sweep_job_[0-9][0-9]/trial_*.csv")
                     if trial_re.match(p.name))
if not trial_files:
    sys.exit(f"no trial CSVs under {sweep_dir}")

with open(trial_files[0]) as fp:
    header = next(csv.reader(fp))
temp_idxs = [i for i, c in enumerate(header) if re.fullmatch(r"\d+C", c)]
temps = [int(header[i][:-1]) for i in temp_idxs]
n_temps = len(temps)
n_trials = len(trial_files)

buf = np.empty((1500, n_temps), dtype=np.float64)

per_temp_std = np.zeros((n_trials, n_temps))
optimized_std = np.zeros(n_trials)

for ti, path in enumerate(trial_files):
    with open(path, "rb") as fp:
        text = fp.read()
        os.posix_fadvise(fp.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
    lines = text.split(b"\n")
    n_rows = 0
    for line in lines[1:]:
        if not line:
            continue
        parts = line.split(b",")
        for j, idx in enumerate(temp_idxs):
            buf[n_rows, j] = float(parts[idx])
        n_rows += 1
    arr = buf[:n_rows]
    per_temp_std[ti] = arr.std(axis=0, ddof=1)
    optimized_std[ti] = arr.max(axis=1).std(ddof=1)
    del text, lines

fig, ax = plt.subplots(figsize=(11, 6))
data   = [per_temp_std[:, c] for c in range(n_temps)] + [optimized_std]
labels = [f"{t}°C" for t in temps]                    + ["best T\nper addr"]
ax.boxplot(data, labels=labels, whis=(0, 100), widths=0.6)
ax.set_yscale("log")
ax.set_ylabel("Std. dev. of amplification ratio across addresses")
ax.set_title(f"Per-trial address-ratio spread across {n_trials} random F/R pairings")
ax.grid(True, axis="y", which="both", alpha=0.25)
ax.axvline(n_temps + 0.5, color="gray", linestyle=":", alpha=0.5)
fig.tight_layout()
fig.savefig(out_png, dpi=120)
plt.close(fig)
print(f"wrote {out_png}: {n_trials} trials, {n_temps + 1} boxes")
