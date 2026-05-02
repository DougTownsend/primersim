# Code review

Suggestions for readability, performance, and accuracy on the current
tree (post-mpfr-removal). Each item cites file:line and notes whether
it's a bug, a cleanup, or a tradeoff.

> **Project state at the time of this revision**
>
> - `Psim_f` / mpfr removed entirely. Concentrations are `Real` (a
>   `typedef` in [include/eq.hpp](../include/eq.hpp)) — `double` by
>   default, `-DPSIM_USE_LONG_DOUBLE` swaps in 80-bit long double.
> - Code split across [src/eq.cpp](../src/eq.cpp) (EQ class /
>   solver), [src/sim.cpp](../src/sim.cpp) (sim_pcr and helpers),
>   [src/address_eval.cpp](../src/address_eval.cpp) (older address
>   eval / shuffle / primer assignment), [src/test_eq.cpp](../src/test_eq.cpp).
> - Several earlier suggestions are already done — see the tail of
>   this doc for the changelog.

---

## Correctness bugs (still open)

### 1. `read_addresses` allocates `r` / `r_rc` with `f`-length
[src/sim.cpp:81-84](../src/sim.cpp#L81-L84)

```cpp
a.f    = (char *)malloc(strlen(tmp_f)+1);
a.r    = (char *)malloc(strlen(tmp_f)+1);  // ← uses tmp_f length
a.f_rc = (char *)malloc(strlen(tmp_f)+1);
a.r_rc = (char *)malloc(strlen(tmp_f)+1);
```

If `tmp_r` is longer than `tmp_f` the next `strcpy(a.r, tmp_r)` writes
past the buffer. Today both inputs are 20 nt so it's latent. Fix:
size each buffer with its own `strlen`, or switch to `std::string`.

### 2. `eval_thread` reads uninitialized `eq.tmp[1]`
[src/address_eval.cpp:291-298](../src/address_eval.cpp#L291-L298)

```cpp
eq.tmp[0] = eq.c[FX] + eq.c[RX];
eq.tmp[2] = std::fmin(eq.tmp[0], eq.tmp[1]);   // ← eq.tmp[1] never set
...
eq.tmp[0] = std::fmax(eq.tmp[0], eq.tmp[1]);
```

When the Y-branch was dropped (A/B/Y removal) the line that wrote
`eq.tmp[1]` went with it but the readers stayed. The min/max read
whatever was last in `tmp[1]`. Either set `tmp[1] = 0` (degenerate
but defined — the lin/exp split collapses) or rewrite without
`min`/`max` now that there's only one nonspec pool.

This is the most likely reason the post-mpfr-removal address-eval
output is suspicious if you ever exercise it. `sim_pcr` (the active
path) doesn't go through this code — but `evaluate_addresses` /
`eval_thread` do.

### 3. `read_addresses` leaks on second call
[src/sim.cpp:67-91](../src/sim.cpp#L67-L91)

`addresses.clear()` drops the vector but doesn't `free` the four
`char *` inside each `address`. The TODO acknowledges this. Add a
per-entry free loop before the clear, or move to `std::string`.

---

## Cleanups

### 4. Right-size the c, c0, k arrays
[include/eq.hpp:106-117](../include/eq.hpp#L106-L117)

```cpp
Real c[19];
Real c0[6];
Real k[13];
```

The active enums use 10 / 3 / 7 of those slots respectively. Resize
to match. With `Real = double` each unused slot is only 8 bytes per
EQ instance (and only one EQ per thread), so the savings on EQ are
negligible — but the pattern propagates to `address_k_conc_vec`,
which does scale with addresses, and the oversized arrays make the
data layout misleading.

### 5. Replace bare-int constants with typed enums
[include/eq.hpp:44-61](../include/eq.hpp#L44-L61)

```cpp
const int F = 0;
const int R = 1;
...
const int K_FH = 0;
const int K_RH = 1;
...
```

Bare ints provide no type-checking. Two parallel sets of int
constants share index 0..whatever, so `c[K_FF]` compiles silently.
Two `enum class Species : int { F=0, ... };` and `enum class Rate :
int { K_FH=0, ... };` would catch the mixup at compile time, and
let `eq.hpp` be `#include`d without polluting the namespace.

### 6. Function-length budget
- `Primeanneal::sim_pcr` — ~220 lines [src/sim.cpp:203-419](../src/sim.cpp#L203-L419)
- `Primeanneal::eval_thread` — ~170 lines [src/address_eval.cpp:154-323](../src/address_eval.cpp#L154-L323)
- `Primeanneal::assign_addresses` — ~75 lines [src/address_eval.cpp:69-143](../src/address_eval.cpp#L69-L143)

Each has clear sub-phases (per-address dhds precompute, per-cycle
update, output). Extract those into named helpers — even just
splitting the inner `for(temp_c = 40 ... 80)` body of `eval_thread`
into a free function makes the outer skimmable.

### 7. Comment what `tmp[0..3]` mean per call site
[src/sim.cpp:158-160](../src/sim.cpp#L158-L160) has the only good
per-tmp documentation in the codebase. Other functions reuse
`eq.tmp[0..3]` with no commentary; one-line comments at the start of
each function (`tmp[0] = nonspec total, tmp[1] = ...`) make the body
legible.

### 8. Output CSV column documentation
The CSVs written by `sim_pcr` ([src/sim.cpp:243, 345, 379](../src/sim.cpp))
and `eval_thread` ([src/address_eval.cpp:317](../src/address_eval.cpp#L317))
have wide column lists with no header row. Either:
- emit a header on first write, or
- add a comment block above the `fprintf` enumerating columns.

### 9. Strip `read_primers_individual`'s `printf` debug spew
[src/address_eval.cpp:43, 49, 52-55](../src/address_eval.cpp#L43)

That function prints "primer is not 20 nt long" warnings to stdout.
Either remove or guard behind a `verbose` flag.

---

## Performance — likely wins

### 10. Cache `dhds_to_eq_const` results
[src/sim.cpp:166-198](../src/sim.cpp#L166-L198), 285-291

`calc_strand_bindings` is called 25× per address per cycle and
computes `dhds_to_eq_const` for the same `(address, primer-end,
primer-end)` pair many times. Hoist the per-cycle dhds → eq_const
conversions out of the end5/end3 loops into a small precomputed
table indexed by `[primer_f_or_r][addr_seq_kind]` once per cycle.

The same pattern applies in `eval_thread`'s temp-c loop.

### 11. Warm-start `solve_eq` across PCR cycles
[src/eq.cpp:57-58](../src/eq.cpp#L57-L58)

```cpp
c[F] = c0[F] / 2.0;   // discards the previous solution
c[R] = c0[R] / 2.0;
calc_cx();
```

Across PCR cycles `c[F]` and `c[R]` change slowly (primers stay
mostly free in the active address's case). Skipping this reset and
reusing the prior solution should drop the outer loop to 1-2
iterations in steady state. Caveat: `c0[F]` itself changes between
cycles, so guard with a "first call" flag and re-init when c0
shifts by more than a tolerance.

### 12. Vectorize `calc_strand_bindings`
[src/sim.cpp:134-201](../src/sim.cpp#L134-L201)

The function is called 5×5 = 25 times per address per cycle, each
time with similar but not identical math. With `Real = double` and
`-mavx2`, the compiler may auto-vectorize the per-(end5, end3)
accumulations, but the address-keyed indirection
(`eq.address_k_conc_vec[bind_addr5].rstrand[end5][end3]`) probably
inhibits it. Hoisting the per-address inner state into a flat array
and unrolling could yield 2-4× on this hotspot — but profile first.

---

## Performance — speculative

### 13. Multi-thread inside `sim_pcr`
The current threading parallelizes across addresses (`eval_addresses_thread`,
[src/sim.cpp:402](../src/sim.cpp#L402)), but a single `sim_pcr` call
is fully serial. The 5×5 loops over end5/end3 in `update_strand_concs`
and `calc_strand_bindings` are independent across address index `i`
— could be parallelized for large address sets. Only worth doing if
single-address simulation is the bottleneck.

### 14. Profile before tuning further
Most of the items above are at-most-2× wins. Run `perf record ./test_eq`
on a representative workload first; the hotspot may be in `thal`
calls during per-cycle dhds setup rather than in the equilibrium
math, since the FP arithmetic is now native and very fast.

---

## Investigations

### 15. Three genuinely failing addresses
After the [k_RH model fix](../src/eq.cpp#L66) (described in
[docs/precision.md](precision.md)), 3 of 979 addresses still fail
to amplify symmetrically — addresses 686, 803, 322. They show the
"both" failure mode (F1 < 0.9 *and* R1 < 0.9 on cycle 1), with
F30/R30 ratios of 0.87-1.05 — i.e. F and R are equally suppressed,
which is consistent with a real primer-design problem (strong
hairpins or cross-dimers) rather than a model artifact.

Worth running these primers through a standalone hairpin / dimer
predictor to confirm. If they're real-world bad primers, document
them. If not, there may still be a residual model issue.

---

## Cleanup of leftover files

Various `bench_$P/` directories may exist at the repo root from the
mpfr precision sweep (they hold per-precision regression CSVs). Same
for `regression_$P.csv` files. After the migration to native FP
these are reference data only — `rm -rf bench_*` plus
`rm regression_{1024,512,256,128,80}.csv` once you've decided you
don't need them.

---

## Done since previous revision (changelog)

The following items from the prior `suggestions.md` are now resolved
in the tree:

- Newton solver replaces binary-search `solve_cf_cr` (faster, simpler).
- `Psim_f` / mpfr removed entirely; native `double` / `long double`
  via `Real` typedef.
- File split: `eq.cpp` (solver) / `sim.cpp` (sim_pcr + helpers) /
  `address_eval.cpp` (older eval path).
- Dead code removed: `solve_cf_cr`, `calc_cf`, `calc_cr`,
  `bounds[]`, `solutions[]`, `saved_frc_conc`, `saved_rrc_conc`,
  `avg_nonspec_amp`, `address_k_conc::dhds[2][4]`, doubled
  `address_index = 0`, commented A/B/Y blocks, commented enum
  blocks, scratch `solve.cpp` at repo root, `CONTEXT.md`.
- Outer-loop iteration cap reduced from `FLOAT_PREC * 2 = 2048` to
  a flat `64`.
- Bool-cast bug in `evaluate_addresses(in_filename, dna_conc)` fixed
  (now passes `false`).
- Inconsistent `.method()`-form vs operator-form arithmetic resolved
  — there's no more wrapper-method form, all native operators.
- Model bug found and fixed: c[R] mass balance had `k_RH·c[F]`
  instead of `k_RH·c[R]`, which was sequestering free R primer in
  proportion to F primer concentration and producing 78 spurious
  "R-stuck" addresses out of 979. Fixed in
  [src/eq.cpp:66](../src/eq.cpp#L66) (k_RH folded into B1) and the
  Jacobian's J11. Post-fix the failure rate is 3/979 ≈ 0.3 % (real
  primer-design issues; see #15 above).
- Build dropped `-lmpfr -lgmp` and the `/share/tuck/dktownse/...`
  paths.
- Outer-loop convergence in `solve_eq` switched from bit-exact `==`
  to a relative-tolerance check
  ([src/eq.cpp:108-117](../src/eq.cpp#L108-L117)) at `16 * epsilon`.
  Output is bit-identical to the pre-change baseline across all 979
  addresses × 30 cycles (full diff via
  [scripts/compare_csv.py](../scripts/compare_csv.py): 364188/364188
  exact match), but the loop is no longer one-ULP-of-rounding away
  from oscillating up to the 64-iteration cap.
