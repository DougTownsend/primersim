#!/usr/bin/env python3
"""Boxplot: distribution of worst-address amplification ratio across
random F/R pairing trials, per temperature, plus one extra box for
"best temperature per address" optimization.

Reads the rolling sweep.*.csv files produced by sweep_pairings.
Each row is "trial_id,f_seq,r_seq,<T1>C,<T2>C,...". Rows of the
same trial are contiguous within a file (the writer mutex
guarantees it), so a streaming group-by-trial_id pass extracts
each trial's per-pair ratios with bounded memory.

Memory-bounded: posix_fadvise(POSIX_FADV_DONTNEED) drops each file
from page cache after read; per-trial work uses a reusable buffer.

Usage: scripts/plot_worst_address.py [sweep_dir] [out_png]
Defaults: sweep_results/, worst_address.png
"""
import csv
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sweep_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "sweep_results")
out_png   = Path(sys.argv[2] if len(sys.argv) > 2 else "worst_address.png")

# Pick up sweep.NNNN.csv either directly under sweep_dir or one level
# down (sweep_job_*/sweep.NNNN.csv) — works for the LSF array layout.
sweep_files = sorted(sweep_dir.glob("sweep.*.csv")) + \
              sorted(sweep_dir.glob("sweep_job_*/sweep.*.csv"))
if not sweep_files:
    sys.exit(f"no sweep.*.csv under {sweep_dir}")

with open(sweep_files[0]) as fp:
    header = next(csv.reader(fp))
temp_idxs = [i for i, c in enumerate(header) if re.fullmatch(r"\d+C", c)]
temps     = [int(header[i][:-1]) for i in temp_idxs]
n_temps   = len(temps)

per_temp_worst  = [[] for _ in range(n_temps)]   # [temp_i][trial]
optimized_worst = []                              # [trial]

# Reusable parse buffer — overwritten per trial.
buf = np.empty((4096, n_temps), dtype=np.float64)
n_rows = 0
current_trial = None

def flush():
    global n_rows
    if n_rows == 0:
        return
    arr = buf[:n_rows]
    col_min = arr.min(axis=0)
    for t, v in enumerate(col_min):
        per_temp_worst[t].append(v)
    optimized_worst.append(arr.max(axis=1).min())
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

n_trials = len(optimized_worst)

fig, ax = plt.subplots(figsize=(11, 6))
data   = per_temp_worst + [optimized_worst]
labels = [f"{t}°C" for t in temps] + ["best T\nper addr"]
ax.boxplot(data, labels=labels, whis=(0, 100), widths=0.6)
ax.set_yscale("log")
ax.set_ylabel("Worst-address amplification ratio")
ax.set_title(f"Worst-address ratio across {n_trials} random F/R pairings")
ax.grid(True, axis="y", which="both", alpha=0.25)
ax.axvline(n_temps + 0.5, color="gray", linestyle=":", alpha=0.5)
fig.tight_layout()
fig.savefig(out_png, dpi=120)
plt.close(fig)
print(f"wrote {out_png}: {n_trials} trials, {n_temps + 1} boxes")
