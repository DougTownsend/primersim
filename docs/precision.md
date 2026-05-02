# Precision sweep — full address-set workload

The codebase historically ran at `FLOAT_PREC = 1024` bits (~308 decimal
digits). [docs/suggestions.md](suggestions.md) flagged this as
overkill. This experiment sweeps mpfr at 1024 / 512 / 256 / 128 / 80
bits, a hardware-`long double` rewrite, and a `double` build,
simulating PCR for **every address** in `addresses.csv` (979
addresses × 30 cycles each) at 55 °C.

After the first sweep flagged a model bug ([k_RH·c[F]] term in the
c[R] equation, see [src/eq.cpp:90](src/eq.cpp#L90) and the related fix
note), the numbers below are from the post-fix run.

## Method

mpfr precision is overridden on the g++ command line:

```bash
g++ -Wall src/eq.cpp src/sim.cpp src/test_eq.cpp src/thal.cpp \
    -Iinclude -O3 -lm -lpthread -lmpfr -lgmp -mavx \
    -DFLOAT_PREC=$P -o bench_$P/test_eq
```

The hardware-FP versions live in [longdouble/](../longdouble) with the
underlying type controlled by a typedef in
[longdouble/include/eq.hpp](../longdouble/include/eq.hpp):

```cpp
#ifdef PSIM_USE_DOUBLE
typedef double Real;
#else
typedef long double Real;
#endif
```

Default builds use 80-bit `long double`; `-DPSIM_USE_DOUBLE` switches
to 64-bit `double`. mpfr is not linked in either.

All seven binaries were launched in parallel from per-precision dirs
so their output files (`regression_new.csv`) didn't collide. One run
per precision because the 1024-bit run alone takes ~50 minutes.

## Results

| version | time | vs 1024 mpfr | match against mpfr-1024 |
|---:|---:|---:|---|
| mpfr 1024 bit (~308 digits) | 50 m 33 s | — | (gold) |
| mpfr  512 bit (~154 digits) | 42 m 42 s | 1.18× | byte-identical (md5 match) |
| mpfr  256 bit ( ~77 digits) | 34 m 43 s | 1.46× | byte-identical |
| mpfr  128 bit ( ~38 digits) | 31 m 30 s | 1.61× | byte-identical |
| mpfr   80 bit ( ~24 digits) | 29 m 23 s | 1.72× | byte-identical |
| `long double` (80-bit hw, ~19 digits) | 16 m 56 s | 2.99× | byte-identical |
| **`double` (64-bit hw, ~16 digits)**  | **15 m 18 s** | **3.30×** | 99.98 % exact, max rel err 4.05e-10 |

For the mpfr precisions and long double, every output file has md5
`d7e02e82b50957301dc3c3810f675600`. The double run differs in 83 of
364,188 cells (0.02 %), all by less than 1e-9 relative — that's at
the print-precision floor, not a meaningful answer change.

Outlier diagnosis ([scripts/diagnose_outliers.py](../scripts/diagnose_outliers.py))
returns the same buckets across all seven precisions:

| mode | count | % |
|---|---:|---:|
| ok (both primers amplify) | 976 | 99.7 % |
| R-stuck | 0 | 0.0 % |
| F-stuck | 0 | 0.0 % |
| both fail | 3 | 0.3 % |

## Read

### mpfr precision sweep

The mpfr speedup curve is gentle and flattens fast. Each halving saves
progressively less wall time:

| step       | time saved | % of 1024 |
|---|---:|---:|
| 1024 → 512 |  7.9 min |  15.6 % |
| 512  → 256 |  8.0 min |  15.8 % |
| 256  → 128 |  3.2 min |   6.3 % |
| 128  →  80 |  2.1 min |   4.2 % |

Below ~256 the per-call mpfr_init/clear and dispatch overhead
dominates the per-bit limb arithmetic, so further bit reductions buy
less.

### long double vs mpfr

Native hardware `long double` is **2.99× faster than 1024-bit mpfr
and 1.74× faster than 80-bit mpfr at the same nominal precision** —
the same precision in software costs ~74 % more time than in
hardware. That overhead is mpfr's per-call function dispatch,
allocation, and rounding, none of which scale with limb count.

Going to long double also dropped two library dependencies (`-lmpfr
-lgmp`) and let `Psim_f` go away entirely:

```cpp
// before (mpfr)
eq.tmp[2].mul(k[K_FX], c[F]);
eq.tmp[2].add(eq.tmp[2], eq.tmp[3]);
eq.tmp[2].add_d(eq.tmp[2], 1.0);
c[X].div(c0[X], eq.tmp[2]);

// after (long double / double)
c[X] = c0[X] / (1.0 + k[K_FX]*c[F] + k[K_RX]*c[R]);
```

### double vs long double

Going from 80-bit hardware (long double) to 64-bit hardware (double)
saves another ~10 % wall time. On Linux x86-64 doubles run through
SSE2 vector units; long doubles drop back to the x87 FPU which is
slightly slower per op.

The 19→16-digit precision drop is invisible at the `%.9` output
format — 364,105 / 364,188 cells (99.98 %) are bit-identical to the
mpfr-1024 gold. The 83 cells that differ are all within ~1e-10
relative error, which is below 9-digit print precision. The outlier
analysis returns the same answer.

## Recommendation

**Switch the main tree to `double`.** It's a 3.30× speedup over
1024-bit mpfr, drops two library dependencies, simplifies the code,
and produces the same answer to 9+ digits across every cell. Both
upstream inputs (thal `dh`/`ds`) and the simulation magnitudes are
already double-friendly — there's nothing in the workload that needs
the extra precision.

`long double` would be a defensible alternative (still ~3× over
mpfr-1024, byte-identical to the gold), but the 10 % extra speed and
SSE2-friendly type alignment make `double` the better default. The
typedef in [eq.hpp](../longdouble/include/eq.hpp) makes both buildable
from the same source — pass `-DPSIM_USE_DOUBLE` for double, omit it
for long double.

If a future workload turns out to need more precision (very stiff
parameters, much longer cycle counts where rounding accumulates),
bringing mpfr back is mechanical — it's a `Psim_f`-shaped slot in the
codebase and the typedef points there. Today, no measurement justifies
that cost.

## Caveats

- All hardware-FP results are Linux x86-64 specific. On Windows MSVC
  long double = double (64-bit). On some ARM/POWER long double is
  128-bit IEEE quad, software-emulated.
- Byte-identical CSVs prove agreement to 9 print digits. Underlying
  math may differ in the 10th-19th digit; that's invisible in this
  output format. The double run's 83 different-last-digit cells
  confirm precision differences exist — they're just below the
  display threshold.
- The 80-bit mpfr vs 80-bit long double comparison is *not quite*
  apples-to-apples: mpfr's significand is 80 bits, long double's is
  64 bits (with a 16-bit exponent for 80 total). So mpfr-80 actually
  has slightly more precision than long double, yet long double still
  wins by 1.74× on speed.
