#!/usr/bin/env python3
"""Convert a paired-primer CSV (column 1 = primer, column 2 =
reverse complement, optional metadata header line) into the input
format test_eq expects:

  primers.csv  — one primer sequence per line.
  pairings.csv — trivial sequential (0,1), (2,3), ... pairing.

The reverse-complement column is dropped because the simulator
recomputes RC internally from each primer. The sweep generates
random pairings at runtime, but test_eq still reads pairings.csv on
startup, so a trivial one is written for compatibility.

Usage: scripts/import_primers.py <input_csv> [primers.csv] [pairings.csv]
"""
import re
import sys
from pathlib import Path

if len(sys.argv) < 2:
    sys.exit(__doc__)

src       = Path(sys.argv[1])
primers_p = Path(sys.argv[2] if len(sys.argv) > 2 else "primers.csv")
pairs_p   = Path(sys.argv[3] if len(sys.argv) > 3 else "pairings.csv")

VALID = re.compile(r"^[ACGT]+$")

primers = []
with open(src) as f:
    for ln, line in enumerate(f, 1):
        first = line.strip().split(",")[0]
        if VALID.match(first):
            primers.append(first)
        # else: header / metadata / blank — silently skip

if len(primers) % 2:
    sys.exit(f"odd primer count ({len(primers)}); pairings need an even count")

with open(primers_p, "w") as f:
    for p in primers:
        f.write(p + "\n")

with open(pairs_p, "w") as f:
    for i in range(0, len(primers), 2):
        f.write(f"{i},{i+1}\n")

print(f"wrote {primers_p}: {len(primers)} primers")
print(f"wrote {pairs_p}: {len(primers)//2} sequential pairings")
