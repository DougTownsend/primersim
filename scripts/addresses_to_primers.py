#!/usr/bin/env python3
"""Convert addresses.csv (paired F,R format) into:

  primers.csv  — one primer per line (2N lines for N addresses)
  pairings.csv — one "f_idx,r_idx" line per pair (N lines)

Run once after addresses.csv is in place. The generated pairings.csv
is the trivial sequential pairing (0,1), (2,3), ..., reproducing the
original addresses.csv's pairing — useful for verifying that the
flat-primer refactor preserves the regression.

Usage: scripts/addresses_to_primers.py [addresses.csv] [primers.csv] [pairings.csv]
"""
import sys

addr_path     = sys.argv[1] if len(sys.argv) > 1 else "addresses.csv"
primers_path  = sys.argv[2] if len(sys.argv) > 2 else "primers.csv"
pairings_path = sys.argv[3] if len(sys.argv) > 3 else "pairings.csv"

VALID = set("ACGT")

primers = []
pairs = []
with open(addr_path) as f:
    for line in f:
        parts = line.strip().split(",")
        if len(parts) < 2:
            continue
        f_seq, r_seq = parts[0], parts[1]
        if not (f_seq and r_seq and set(f_seq) <= VALID and set(r_seq) <= VALID):
            continue
        f_idx = len(primers); primers.append(f_seq)
        r_idx = len(primers); primers.append(r_seq)
        pairs.append((f_idx, r_idx))

with open(primers_path, "w") as f:
    for p in primers:
        f.write(p + "\n")

with open(pairings_path, "w") as f:
    for f_idx, r_idx in pairs:
        f.write(f"{f_idx},{r_idx}\n")

print(f"wrote {len(primers)} primers to {primers_path}")
print(f"wrote {len(pairs)} pairings to {pairings_path}")
