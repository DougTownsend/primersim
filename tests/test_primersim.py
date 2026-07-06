"""Phase-two acceptance tests for the primersim Python API.

Run from the repo root after `pip install .`:

    python3 tests/test_primersim.py

Covers, per the phase-two handoff:
  * Python <-> C++ parity against regression_new.csv (100 primers / 50 pairings)
  * Parallel determinism (threads=1 vs threads=8)
  * Item-9 concentration-matrix / ratio consistency
  * Address-token parsing + write_pairings-format round-trip
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

from primersim import Simulator, parse_token

# Exact params test_eq uses to generate regression_new.csv.
PRIMERS = os.path.join(REPO, "primers.csv")
PAIRINGS = os.path.join(REPO, "pairings.csv")
REGRESSION = os.path.join(REPO, "regression_new.csv")
DNA_TOTAL = 3e-15
PRIMER = 250e-9
MV, DV, DNTP = 100.0, 1.5, 0.2
CYCLES = 30
TEMP = 55.0


def read_pairings(path):
    """Parse pairings.csv into a list of (f_token, r_token) string pairs."""
    pairs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            a, b = line.split(",")
            pairs.append((a, b))
    return pairs


def regression_ratios(path):
    """Recover each address's amplification ratio from the cycle-30 log row.

    ratio = (target_F + target_R) / (nonspec_F + nonspec_R)
          = (col7 + col8) / (col9 + col11)   [1-indexed, final cycle]
    """
    ratios = {}
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("Simulating address"):
                cur = int(line.split()[-1])
            elif line and cur is not None and line[:2] == "%02d" % CYCLES:
                c = line.split(",")
                spec = float(c[6]) + float(c[7])
                nonspec = float(c[8]) + float(c[10])
                ratios[cur] = spec / nonspec
    return [ratios[i] for i in sorted(ratios)]


def build_sim(num_cpu):
    sim = Simulator(PRIMERS, mv=MV, dv=DV, dntp=DNTP, num_cpu=num_cpu)
    sim.set_addresses(read_pairings(PAIRINGS))
    # test_eq seeds each address with dna_total / N (uniform split).
    sim.set_address_conc(DNA_TOTAL / sim.num_addresses)
    sim.set_primer_conc(PRIMER)
    sim.set_temperature(TEMP)
    return sim


def relerr(a, b):
    denom = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / denom


def test_parity():
    expected = regression_ratios(REGRESSION)
    sim = build_sim(num_cpu=8)
    assert sim.num_primers == 100, sim.num_primers
    assert sim.num_addresses == 50, sim.num_addresses
    res = sim.simulate_all(cycles=CYCLES)
    assert len(res.ratios) == len(expected) == 50
    worst = max(relerr(g, e) for g, e in zip(res.ratios, expected))
    print(f"  parity: worst relative error vs regression = {worst:.3e}")
    # Extension is built without -mavx by default; the standalone
    # regression uses -mavx. Allow a small FP tolerance for codegen drift.
    assert worst < 1e-6, f"parity relerr too large: {worst:.3e}"


def test_determinism():
    r1 = build_sim(num_cpu=1).simulate_all(cycles=CYCLES, threads=1)
    r8 = build_sim(num_cpu=8).simulate_all(cycles=CYCLES, threads=8)
    assert r1.ratios == r8.ratios, "ratios differ across thread counts"
    assert r1.concentrations == r8.concentrations, "conc matrix differs across thread counts"
    print("  determinism: threads=1 and threads=8 bit-identical")


def test_item9_consistency():
    sim = build_sim(num_cpu=4)
    res = sim.simulate_all(cycles=CYCLES)
    worst = 0.0
    for a in range(sim.num_addresses):
        row = res.concentrations[a]
        target = row[a]
        offtarget = sum(row[i] for i in range(len(row)) if i != a)
        derived = target / offtarget
        worst = max(worst, relerr(derived, res.ratios[a]))
    print(f"  item9: worst relerr (target/off-target vs reported ratio) = {worst:.3e}")
    assert worst < 1e-9, f"concentration/ratio inconsistency: {worst:.3e}"


def test_token_parsing():
    assert parse_token("2r") == (2, True)
    assert parse_token("2R") == (2, True)
    assert parse_token(3) == (3, False)
    assert parse_token("7") == (7, False)

    sim = Simulator(PRIMERS, mv=MV, dv=DV, dntp=DNTP, num_cpu=1)
    sim.add_address(1, "2r")
    f_idx, r_idx, f_rc, r_rc = sim.addresses[0]
    assert (f_idx, r_idx, f_rc, r_rc) == (1, 2, False, True), sim.addresses[0]
    # Round-trip against the write_pairings "<idx>[r]" token format.
    assert sim.address_tokens()[0] == ("1", "2r"), sim.address_tokens()[0]
    print("  token parsing: add_address(1,'2r') -> (1,2,F,T), round-trips to ('1','2r')")


def main():
    if not os.path.exists(REGRESSION):
        sys.exit(f"missing {REGRESSION}; run `make && ./test_eq 24` first")
    failures = 0
    for name in ("test_parity", "test_determinism", "test_item9_consistency",
                 "test_token_parsing"):
        try:
            print(f"[{name}]")
            globals()[name]()
            print(f"  PASS\n")
        except AssertionError as e:
            failures += 1
            print(f"  FAIL: {e}\n")
    if failures:
        sys.exit(f"{failures} test(s) failed")
    print("all tests passed")


if __name__ == "__main__":
    main()
