#!/usr/bin/env python3
"""Scan a regression CSV from test_eq, extract per-address cycle-30
amplification rates, and rank addresses by how poorly they amplify.

Output columns:
  addr, fwd_amp, rev_amp, fwd_total_30, rev_total_30, F_seq, R_seq
"""
import sys
import re
from pathlib import Path


def parse_regression(path):
    """Yield (addr, list_of_30_cycle_rows) tuples.

    Each cycle row is the 12-element list-of-floats for that cycle.
    """
    cur_addr = None
    cur_cycles = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = re.match(r"Simulating address (\d+)", line)
            if m:
                if cur_addr is not None:
                    yield cur_addr, cur_cycles
                cur_addr = int(m.group(1))
                cur_cycles = []
                continue
            parts = line.split(",")
            try:
                row = [float(p) for p in parts]
            except ValueError:
                continue
            cur_cycles.append(row)
        if cur_addr is not None:
            yield cur_addr, cur_cycles


def parse_addresses(path):
    """Return list of (F, R) tuples indexed by address number."""
    out = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < 2:
                continue
            out.append((parts[0], parts[1]))
    return out


def main():
    if len(sys.argv) < 2:
        print("usage: find_outliers.py <regression.csv> [addresses.csv]")
        sys.exit(1)

    reg_path = sys.argv[1]
    addr_path = sys.argv[2] if len(sys.argv) > 2 else "addresses.csv"

    addr_seqs = parse_addresses(addr_path)
    rows = []  # (addr, fwd_amp_30, rev_amp_30, fwd_total_30, rev_total_30)
    for addr, cycles in parse_regression(reg_path):
        if not cycles:
            continue
        # cycles[0] is cycle 00, cycles[30] is cycle 30 if all present
        # we want the last cycle
        last = cycles[-1]
        # columns: cycle, temp, fwd_amp, rev_amp, c0F, c0R, totalF, totalR,
        #          nonspec_f, nonspec_f_amp, nonspec_r, nonspec_r_amp
        if len(last) < 8:
            continue
        rows.append((addr, last[2], last[3], last[6], last[7]))

    # Combined "amplification quality": geometric mean of fwd & rev
    # (something both has to be high for PCR to work)
    def quality(r):
        f, rv = r[1], r[2]
        return min(f, rv)  # bottlenecked by whichever primer is slower

    rows.sort(key=quality)

    print(f"{'addr':>5} {'fwd_amp':>10} {'rev_amp':>10} {'fwd_tot30':>14} {'rev_tot30':>14}  {'F_seq':<22}{'R_seq':<22}")
    print("-" * 110)
    print("Worst 20:")
    for r in rows[:20]:
        addr = r[0]
        f_seq, r_seq = addr_seqs[addr] if addr < len(addr_seqs) else ("?","?")
        print(f"{addr:>5} {r[1]:>10.6f} {r[2]:>10.6f} {r[3]:>14.4e} {r[4]:>14.4e}  {f_seq:<22}{r_seq:<22}")

    print()
    print("Best 5:")
    for r in rows[-5:]:
        addr = r[0]
        f_seq, r_seq = addr_seqs[addr] if addr < len(addr_seqs) else ("?","?")
        print(f"{addr:>5} {r[1]:>10.6f} {r[2]:>10.6f} {r[3]:>14.4e} {r[4]:>14.4e}  {f_seq:<22}{r_seq:<22}")

    print()
    n_total = len(rows)
    n_below_50 = sum(1 for r in rows if quality(r) < 0.5)
    n_below_10 = sum(1 for r in rows if quality(r) < 0.1)
    print(f"total addresses: {n_total}")
    print(f"  cycle-30 min(fwd_amp, rev_amp) < 0.50: {n_below_50}  ({100.0*n_below_50/n_total:.1f}%)")
    print(f"  cycle-30 min(fwd_amp, rev_amp) < 0.10: {n_below_10}  ({100.0*n_below_10/n_total:.1f}%)")


if __name__ == "__main__":
    main()
