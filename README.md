# primersim — Python interface

`primersim` simulates PCR amplification for a DNA-data-storage primer set. You
give it a pool of primers and a set of **addresses** (forward/reverse primer
pairs), and it reports, per address, how well the intended product amplifies
relative to off-target products.

The thermodynamics and the parallel numerics run in a compiled C++ core
(`primersim._core`); the Python layer (`primersim.Simulator`) is where you
describe the experiment and read the results.

---

## Installation

`pip install` builds the C++ extension automatically.

```bash
cd primersim
pip install .
```

For a faster (machine-specific, non-portable) build that enables native SIMD
codegen:

```bash
PRIMERSIM_NATIVE=1 pip install .
```

For development (builds the extension in place so it's importable from the repo
root too):

```bash
pip install -e .
```

**Import gotcha:** a plain `pip install .` puts the extension in site-packages,
but the repo root contains a `primersim/` source directory that will *shadow* it
if you `import primersim` while your working directory is the repo root. Run your
scripts from any other directory, or use `pip install -e .`.

Requirements: Python ≥ 3.10, a C++17 compiler, pybind11 ≥ 2.12 (pulled in as a
build dependency).

---

## Quick start

```python
from primersim import Simulator

# Build the simulator. The primers file is one ACGT sequence per line.
# mv/dv/dntp (mM) are baked into the thermodynamics cache built here — the
# one expensive step. num_cpu sets the default worker-thread count.
sim = Simulator("primers.csv", mv=100.0, dv=1.5, dntp=0.2, num_cpu=8)

# Define addresses = (forward primer, reverse primer). A slot is a primer
# index (literal sequence) or a "<idx>r" token (that primer's reverse
# complement), matching the pairings.csv convention ("0r,4").
sim.add_address(0, 1)          # primer 0 forward, primer 1 reverse
sim.add_address(1, "2r")       # primer 1 forward, RC of primer 2 reverse

# Reaction setup.
sim.set_address_conc(1.5e-15)  # initial template conc per address (mol/L)
sim.set_primer_conc(250e-9)    # primer conc (mol/L); F and R set equal
sim.set_temperature(55.0)      # constant anneal temp (deg C)

# Simulate every address in parallel.
res = sim.simulate_all(cycles=30)
print(res.ratios)              # one amplification ratio per address
print(res.concentrations)      # N x N matrix of final concentrations
```

---

## Concepts

### Primers and addresses

- **Primers** are loaded once from the file passed to the constructor (one
  sequence per line). Index them 0..N-1 in file order.
- An **address** is a forward/reverse primer pair — the thing PCR amplifies.
  Each slot is:
  - an `int` → use that primer's sequence as written, or
  - a `"<idx>r"` string (e.g. `"2r"`) → use that primer's **reverse complement**.
- Primers may be reused across many addresses; the thermodynamics cache is keyed
  by primer, so this costs no extra memory.

Ways to set addresses:

```python
sim.add_address(1, "2r")                       # append one
sim.set_addresses([(0, 1), (0, "2r"), (3, 4)]) # bulk replace
sim.enumerate_pairs()                          # all binom(N, 2) primer pairs
sim.enumerate_pairs([0, 1, 2, 5])              # binom(k, 2) over a subset
```

### Cations and the cache

`mv` (monovalent), `dv` (divalent), and `dntp` concentrations (all mM) enter the
nearest-neighbor entropy, so they are **baked into the dH/dS cache** built at
construction. Changing them requires rebuilding that cache:

```python
sim.set_cations(mv=50.0, dv=2.0, dntp=0.2)   # rebuilds the cache (expensive)
```

This is the one costly operation (a full pairwise thermodynamic sweep over
`num_cpu` threads) and it emits a warning. Do not call it in a loop. Temperature
is **not** part of the cache — it is applied per cycle and is free to change.

### Temperature (constant or touchdown)

```python
sim.set_temperature(58.0)                    # constant for all cycles
sim.set_temperature([60, 60, 59, 58, ...])   # per-cycle profile (touchdown PCR)
```

A per-cycle profile's length must equal the `cycles` argument at simulate time.

### Concentrations

- `set_address_conc(x)` — a scalar seeds every address with the same initial
  template concentration; a sequence gives a per-address value (its length must
  equal the number of addresses).
- `set_primer_conc(x)` sets forward and reverse primer to `x`;
  `set_primer_conc(f=..., r=...)` sets them independently (they are distinct,
  independently-depleted species in the model).

All concentrations are in mol/L.

---

## Running simulations

### One address

```python
r = sim.simulate_address(addr=0, cycles=30)
r.ratio            # float: target conc / sum of off-target concs
r.concentrations   # list of N floats: final total DNA (F + R strands) of
                   # every address present in this run
```

### All addresses (parallel)

```python
res = sim.simulate_all(cycles=30, threads=8)  # threads defaults to num_cpu
res.ratios          # list of N floats, one amplification ratio per address
res.concentrations  # N x N matrix; row a is address a's simulation, and
                    # entry [a][i] is address i's final F+R total in that run
```

Each worker thread owns its own mutable simulation state; the shared cache is
read-only during simulation, and the GIL is released for real parallelism.
Results are deterministic and independent of the thread count.

Before simulating you must have set: at least one address, the address
concentration, the primer concentration, and the temperature — otherwise a
`RuntimeError` explains what's missing.

---

## Introspection

```python
sim.num_primers        # count of loaded primers
sim.num_addresses      # count of configured addresses
sim.addresses          # [(f_idx, r_idx, f_rc, r_rc), ...]
sim.address_tokens()   # [("1", "2r"), ...]  round-trips the pairings.csv format

from primersim import parse_token
parse_token("2r")      # (2, True)
parse_token(3)         # (3, False)
```

---

## Worked example: shared reverse primer

A common design is many forward primers sharing one reverse primer. If the
reverse primer is the last line of the file (index N-1):

```python
from primersim import Simulator

sim = Simulator("primers.csv", mv=50.0, dv=2.0, dntp=0.2, num_cpu=8)
rev = sim.num_primers - 1
sim.set_addresses([(f, rev) for f in range(rev)])   # each fwd + shared rev

sim.set_address_conc(3e-15 / sim.num_addresses)
sim.set_primer_conc(250e-9)
sim.set_temperature(58.0)

res = sim.simulate_all(cycles=30)
for a, ratio in enumerate(res.ratios):
    print(a, ratio)
```

---

## Testing

Acceptance tests live in `tests/test_primersim.py`. They check parity against
the standalone C++ regression, thread determinism, concentration/ratio
consistency, and token parsing.

```bash
make && ./test_eq 24        # generate regression_new.csv (the parity baseline)
pip install .
python3 tests/test_primersim.py
```
