#!/usr/bin/env python3
"""Boxplot: distribution of worst-address amplification ratio across
random F/R pairing trials, per temperature, plus one extra box for
"best temperature per address" optimization.

For each trial CSV (rows = addresses, columns = ratios at 50..60 °C):
  • Per-temperature box: the minimum (worst) address ratio at that
    temp — one value per trial.
  • Optimized box: for each address pick its maximum (best) ratio
    across temperatures, then take the minimum over addresses —
    one value per trial. This is the best achievable worst-case if
    every address could be run at its preferred annealing temp.

Usage: scripts/plot_worst_address.py [sweep_dir] [out_png]
Defaults: sweep_results/, worst_address.png
"""
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sweep_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "sweep_results")
out_png   = Path(sys.argv[2] if len(sys.argv) > 2 else "worst_address.png")

trial_files = sorted(sweep_dir.glob("sweep_job_*/trial_*.csv"))
if not trial_files:
    sys.exit(f"no trial CSVs under {sweep_dir}")

df0 = pd.read_csv(trial_files[0])
temp_cols = [c for c in df0.columns if re.fullmatch(r"\d+C", c)]
temps = [int(c[:-1]) for c in temp_cols]

per_temp_worst = {t: [] for t in temps}
optimized_worst = []

for path in trial_files:
    df = pd.read_csv(path)
    for col, t in zip(temp_cols, temps):
        per_temp_worst[t].append(df[col].min())
    optimized_worst.append(df[temp_cols].max(axis=1).min())

fig, ax = plt.subplots(figsize=(11, 6))
data   = [per_temp_worst[t] for t in temps] + [optimized_worst]
labels = [f"{t}°C" for t in temps]          + ["best T\nper addr"]
ax.boxplot(data, labels=labels, whis=(0, 100), widths=0.6)
ax.set_yscale("log")
ax.set_ylabel("Worst-address amplification ratio")
ax.set_title(f"Worst-address ratio across {len(trial_files)} random F/R pairings")
ax.grid(True, axis="y", which="both", alpha=0.25)
ax.axvline(len(temps) + 0.5, color="gray", linestyle=":", alpha=0.5)
fig.tight_layout()
fig.savefig(out_png, dpi=120)
print(f"wrote {out_png}: {len(trial_files)} trials, {len(temps)+1} boxes")
