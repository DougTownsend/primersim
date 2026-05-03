# Code review — sim_pcr / solve_eq

Suggestions for the active simulation path: `solve_eq`
([src/eq.cpp](../src/eq.cpp)), `sim_pcr` and its helpers
([src/sim.cpp](../src/sim.cpp)), and the shared classes in
[include/eq.hpp](../include/eq.hpp).

Items specific to the address-evaluation path
([src/address_eval.cpp](../src/address_eval.cpp)) live in
[docs/address_suggestions.md](address_suggestions.md).

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
[src/sim.cpp:88-91](../src/sim.cpp#L88-L91)

```cpp
a.f    = new char[strlen(tmp_f)+1];
a.r    = new char[strlen(tmp_f)+1];  // ← uses tmp_f length
a.f_rc = new char[strlen(tmp_f)+1];
a.r_rc = new char[strlen(tmp_f)+1];
```

If `tmp_r` is longer than `tmp_f` the next `strcpy(a.r, tmp_r)` writes
past the buffer. Today both inputs are 20 nt so it's latent. Fix:
size each buffer with its own `strlen`, or switch to `std::string`.

---

## Cleanups

### 2. Right-size the c, c0, k arrays
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

### 3. Replace bare-int constants with typed enums
[include/eq.hpp:45-62](../include/eq.hpp#L45-L62)

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
Two `enum class Species : int { F=0, ... };` and
`enum class Rate : int { K_FH=0, ... };` would catch the mixup at
compile time, and let `eq.hpp` be `#include`d without polluting the
namespace.

### 4. Function-length budget in `sim_pcr`
[src/sim.cpp:194-398](../src/sim.cpp#L194-L398)

`Primeanneal::sim_pcr` is ~205 lines. It has clear sub-phases
(per-address dhds precompute, K-cache populate, per-cycle update,
output). Extract those into named helpers — even just naming the
"populate per-cycle K cache" block as a function would make
`sim_pcr` skimmable.

### 5. Output CSV column documentation
The cycle rows written by `sim_pcr`
([src/sim.cpp:268, 362](../src/sim.cpp)) and the continuation row
([src/sim.cpp:380](../src/sim.cpp#L380)) have wide column lists
with no header row — the only existing hint is a comment above the
fprintf. Either:
- emit a header on first write, or
- expand the comment block to enumerate every column position.

---

## Performance — likely wins

### 6. Vectorize `calc_strand_bindings`
[src/sim.cpp:134-191](../src/sim.cpp#L134-L191)

The function is called 5×5 = 25 times per address per cycle, each
time with similar but not identical math. With `Real = double` and
`-mavx2`, the compiler may auto-vectorize the per-(end5, end3)
accumulations, but the address-keyed indirection
(`eq.address_k_conc_vec[bind_addr5].rstrand[end5][end3]`) probably
inhibits it. Hoisting the per-address inner state into a flat array
and unrolling could yield 2-4× on this hotspot — but profile first.

---

## Performance — speculative

### 7. Multi-thread inside `sim_pcr`
The current threading parallelizes across addresses
(`eval_addresses_thread`, [src/sim.cpp:400](../src/sim.cpp#L400)),
but a single `sim_pcr` call is fully serial. The outer i-loops in
the per-cycle code (calling `calc_strand_bindings` and
`update_strand_concs`) are independent across address index `i` —
could be parallelized for large address sets. Only worth doing if
single-address simulation is the bottleneck.

### 8. Profile before tuning further
Most of the items above are at-most-2× wins. Run `perf record
./test_eq` on a representative workload first; the hotspot may be
in `thal` calls during per-cycle dhds setup rather than in the
equilibrium math, since the FP arithmetic is now native and very
fast. (See for instance the `dhds_to_eq_const` cache in the
changelog: the work it eliminated was real but the wall-time impact
was only ~1.9 %.)

---

## Investigations

### 9. Three genuinely failing addresses
After the [k_RH model fix](../src/eq.cpp#L75) (described in
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

Various `bench_$P/` and `verify*/` directories may exist at the
repo root from the precision sweeps and verification runs (they
hold per-precision regression CSVs). After the migration to native
FP these are reference data only — `rm -rf bench_* verify*` plus
`rm regression_{1024,512,256,128,80}.csv` once you've decided you
don't need them.

---

## Done since previous revision (changelog)

Project-wide history of resolved items:

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
  [src/eq.cpp:75](../src/eq.cpp#L75) (k_RH folded into B1) and the
  Jacobian's J11. Post-fix the failure rate is 3/979 ≈ 0.3 % (real
  primer-design issues; see #9 above).
- Build dropped `-lmpfr -lgmp` and the `/share/tuck/dktownse/...`
  paths.
- Outer-loop convergence in `solve_eq` switched from bit-exact `==`
  to a relative-tolerance check
  ([src/eq.cpp:121-123](../src/eq.cpp#L121-L123)) at `16 * epsilon`.
  Output is bit-identical to the pre-change baseline across all 979
  addresses × 30 cycles (full diff via
  [scripts/compare_csv.py](../scripts/compare_csv.py): 364188/364188
  exact match), but the loop is no longer one-ULP-of-rounding away
  from oscillating up to the 64-iteration cap.
- `dhds_to_eq_const` results cached per (address, cycle): added
  `tmp_dhds_K[2][4]` and `dhds_K_primer_*_addr_{frc,rrc}` Real fields
  to `address_k_conc`, populated once per cycle in `sim_pcr`'s outer
  loop, read by `calc_strand_bindings` and `update_strand_concs`.
  Saves ~85 % of the `exp()` calls in the cycle loop. **Wall time
  impact: ~1.9 %** — the function was already cheap enough that
  caching barely registered. Output bit-identical to baseline.
- `solve_eq` warm-starts from the prior solution. Added
  `bool warm_start_valid = false` to EQ; first call resets
  `c[F]`/`c[R]` to `c0/2`, subsequent calls reuse the previous
  solution as the seed for `calc_cx`. Wall-time impact is in the
  noise (Newton converges in 1-2 outer iters either way); the win is
  more about the loop being more correct in spirit (we're solving
  perturbations of a known fixed point) than about speed. Output
  bit-identical to baseline.
- All `malloc` / `free` calls converted to `new[]` / `delete[]`
  (covers `read_addresses`, `read_primers_individual`, and the
  `~Primeanneal` cleanup). At the same time, fixed the
  `read_addresses` repeat-call leak (the `addresses.clear()` was
  dropping the four `char *` per address) and an analogous leak in
  `read_primers_individual`. `~Primeanneal` now also frees the
  `addresses` strings — previously it only freed `primers`, so
  destruction-time leaks were the norm whenever `read_addresses`
  was used.
- `EQ::tmp[4]` deleted. Every former use site now has a named local
  declared with a comment for what it holds. In the active `sim_pcr`
  path the locals are `nonspec_total` / `sum_f_weighted` /
  `sum_r_weighted` (also threaded into `calc_strand_bindings` as
  three new `Real &` parameters), `init_sum_fstrand` /
  `init_sum_rstrand`, `f_growth_ratio` / `r_growth_ratio`, and the
  `total_nonspec_*` / `last_nonspec_*` / `nonspec_*_growth` family
  for the per-cycle nonspec dump. Output bit-identical to baseline
  (md5 unchanged across 979 addresses × 30 cycles).
