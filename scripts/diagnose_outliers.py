#!/usr/bin/env python3
"""For poorly-amplifying addresses, look at whether the failure is
F-side, R-side, or symmetric, and try to localize it to a primer
feature (hairpin / self-dimer / FR cross-binding).
"""
import sys
import re


def parse_regression(path):
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
    out = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < 2:
                continue
            out.append((parts[0], parts[1]))
    return out


def revcomp(seq):
    comp = {"A":"T","T":"A","G":"C","C":"G"}
    return "".join(comp[c] for c in reversed(seq))


def longest_self_dup(seq):
    """Return length of longest substring of seq that appears reversed-complemented elsewhere in seq.
    Crude proxy for hairpin / self-dimer potential."""
    rc = revcomp(seq)
    best = 0
    for L in range(len(seq), 1, -1):
        for i in range(len(seq) - L + 1):
            if seq[i:i+L] in rc:
                return L
    return 0


def main():
    reg_path = sys.argv[1]
    addr_path = sys.argv[2] if len(sys.argv) > 2 else "addresses.csv"
    addrs = parse_addresses(addr_path)

    # rows: addr, cycle1_F, cycle1_R, cycle30_F_amp, cycle30_R_amp,
    #       cycle30_F_total, cycle30_R_total
    rows = []
    for addr, cycles in parse_regression(reg_path):
        if len(cycles) < 31:
            continue
        c1 = cycles[1]   # cycle 01 (cycles[0] is cycle 00)
        c30 = cycles[-1]
        rows.append((addr, c1[2], c1[3], c30[2], c30[3], c30[6], c30[7]))

    # bucket by failure mode:
    #   "R-stuck" : cycle1 F good (>0.9), cycle1 R bad (<0.1)
    #   "F-stuck" : cycle1 F bad, cycle1 R good
    #   "both"    : both bad on cycle 1
    #   "ok"      : both good on cycle 1
    def bucket(r):
        f1, r1 = r[1], r[2]
        if f1 > 0.9 and r1 > 0.9: return "ok"
        if f1 > 0.9 and r1 < 0.1: return "R-stuck"
        if r1 > 0.9 and f1 < 0.1: return "F-stuck"
        return "both"

    counts = {"ok":0, "R-stuck":0, "F-stuck":0, "both":0}
    for r in rows:
        counts[bucket(r)] += 1
    print("cycle-1 failure modes:")
    for k, v in counts.items():
        print(f"  {k:>10}: {v:4d}  ({100.0*v/len(rows):.1f}%)")
    print()

    # show worst 15 by min(cycle30_F, cycle30_R), with sequence + dup-stretch
    rows.sort(key=lambda r: min(r[3], r[4]))
    print(f"{'addr':>5} {'mode':>8}  {'c1_F':>7} {'c1_R':>7}  {'c30_F':>8} {'c30_R':>8}  "
          f"{'F30/R30':>7}  {'F_dup':>5} {'R_dup':>5}  F_seq                R_seq")
    print("-" * 120)
    for r in rows[:15]:
        addr = r[0]
        f_seq, r_seq = addrs[addr] if addr < len(addrs) else ("?","?")
        f_dup = longest_self_dup(f_seq)
        r_dup = longest_self_dup(r_seq)
        ratio = (r[5]/r[6]) if r[6] != 0 else float('inf')
        print(f"{addr:>5} {bucket(r):>8}  {r[1]:>7.4f} {r[2]:>7.4f}  {r[3]:>8.4f} {r[4]:>8.4f}  "
              f"{ratio:>7.2f}  {f_dup:>5} {r_dup:>5}  {f_seq:<20} {r_seq:<20}")


if __name__ == "__main__":
    main()
