#!/usr/bin/env python3
"""Boxplot: distribution of average amplification ratio across random
F/R pairing trials, per temperature, plus one extra box for "best
temperature per address" optimization.

For each trial CSV (rows = addresses, columns = ratios at 50..60 °C):
  • Per-temperature box: the mean over addresses of the ratio at that
    temp — one value per trial.
  • Optimized box: for each address pick its maximum (best) ratio
    across temps, then take the mean over addresses — one value per
    trial. This is the average amplification if every address ran at
    its preferred annealing temperature.

Usage: scripts/plot_average_ratio.py [sweep_dir] [out_png]
Defaults: sweep_results/, average_ratio.png
"""
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sweep_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "sweep_results")
out_png   = Path(sys.argv[2] if len(sys.argv) > 2 else "average_ratio.png")

trial_files = sorted(sweep_dir.glob("sweep_job_*/trial_*.csv"))
if not trial_files:
    sys.exit(f"no trial CSVs under {sweep_dir}")

df0 = pd.read_csv(trial_files[0])
temp_cols = [c for c in df0.columns if re.fullmatch(r"\d+C", c)]
temps = [int(c[:-1]) for c in temp_cols]

per_temp_avg = {t: [] for t in temps}
optimized_avg = []

for path in trial_files:
    df = pd.read_csv(path)
    for col, t in zip(temp_cols, temps):
        per_temp_avg[t].append(df[col].mean())
    optimized_avg.append(df[temp_cols].max(axis=1).mean())

fig, ax = plt.subplots(figsize=(11, 6))
data   = [per_temp_avg[t] for t in temps] + [optimized_avg]
labels = [f"{t}°C" for t in temps]        + ["best T\nper addr"]
ax.boxplot(data, labels=labels, whis=(0, 100), widths=0.6)
ax.set_yscale("log")
ax.set_ylabel("Average amplification ratio across addresses")
ax.set_title(f"Average amplification ratio across {len(trial_files)} random F/R pairings")
ax.grid(True, axis="y", which="both", alpha=0.25)
ax.axvline(len(temps) + 0.5, color="gray", linestyle=":", alpha=0.5)
fig.tight_layout()
fig.savefig(out_png, dpi=120)
print(f"wrote {out_png}: {len(trial_files)} trials, {len(temps)+1} boxes")
