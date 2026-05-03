# Performance — calc_dimer cache + parallel populate (symmetric)

This experiment caches every `calc_dimer` result that any `sim_pcr`
might need, computes them in parallel across `num_cpu` threads at
startup, and replaces the per-`sim_pcr` `thal` calls with cache
lookups. The cache exploits thal's near-symmetry: only one of
`thal(A, B)` and `thal(B, A)` is stored per pair, so populate work
and storage both halve.

## Setup

The cache is `Primeanneal::dimer_cache`, a packed triangular array of
`dh_ds_cache_entry` covering the canonical pairs `(a, b)` with
`a ≤ b` over `M = 4·N` "sequence slots" (each address has 4 kinds:
f, f_rc, r, r_rc). Index formula:

```
idx(a, b) = a*(2*M - a - 1)/2 + b   for a <= b
```

Total entries: `M*(M+1)/2 ≈ 8·N²` (down from `16·N²` in the
asymmetric version). `cached_dimer(addr1, ki, addr2, kj)` swaps to
canonical ordering before indexing, so both calling directions hit
the same entry.

`populate_dimer_cache` distributes the outer `a` axis cyclically
across `num_cpu` threads (`a % nthreads == t`) — a contiguous slice
would be heavily imbalanced because triangular row `a` has only
`M-a` entries. `sim_pcr` reads only.

`test_eq.cpp` calls `pa.populate_dimer_cache(...)` once after
`read_addresses` and before the simulation loop.

## Wall-time results (N = 979)

| build | wall | speedup vs baseline | regression md5 |
|---|---:|---:|---|
| baseline (no cache)                                          | 14 m 18 s | — | `0b8ca853…` |
| asymmetric cache + parallel populate (ThalReal=double)        | 1 m 27 s | 9.9×  | `0b8ca853…` (matches baseline) |
| asymmetric cache + parallel populate (ThalReal=float)         | 1 m 27 s | 9.9×  | `8a927c84…` |
| **symmetric cache + parallel populate (ThalReal=double)**     | **57.3 s** | **15.0×** | `55b677b1…` (drift, see accuracy) |

The symmetric cache is **34 % faster than the asymmetric version**
and **15× faster than baseline**. md5 differs from the asymmetric
run because the symmetric lookup picks one of `thal(A, B)` /
`thal(B, A)` and `thal` isn't bitwise symmetric — see the accuracy
section below.

CPU time vs wall time:

| build | wall | user | parallelism |
|---|---:|---:|---:|
| asymmetric cache (double)  | 87 s   | 1295 s | 14.9× (avg) |
| **symmetric cache (double)**  | **57.3 s** | **673 s** | **11.7× (avg)** |

The benchmark has two phases. Estimated breakdown from CPU/wall
arithmetic:

| build | populate (wall) | sim_pcr loop (wall) | total |
|---|---:|---:|---:|
| asymmetric cache  | ~39 s  | ~48 s  | 87 s |
| **symmetric cache**  | **~20 s** | **~37 s** | **57 s** |

The populate phase halves cleanly (half the unique pairs). The
sim_pcr loop also got faster — likely from better L1/L2 cache
behavior with a 122 MB working set vs the prior 245 MB, plus
slightly tighter inner-loop indexing.

## Profile (perf, `-O3 -g`, ThalReal=float)

Sampled at 999 Hz with `--call-graph dwarf`, 1.29 M samples over
the 87 s run. perf samples per-thread accurately, no instrumentation
overhead, and exposes hardware counters that gprof can't.

**Inclusive call graph** (time accumulating into children):

```
94.4 %  populate_dimer_cache lambdas (32 worker threads)
   ├── 91.0 %  thal()                            (89.7 % "self")
   │   ├── 86.1 %  fillMatrix_dimer              (inlined into thal)
   │   │   ├── 39.2 %  calc_bulge_internal_dimer (inlined)
   │   │   └── ~46 %   direct DP arithmetic / inlined helpers
   │   ├──  4.4 %  RSH (right-side helper)
   │   └──  1.7 %  LSH (left-side helper)
 1.5 %  sim_pcr (single-threaded cache-lookup + cycle math)
```

So the populate phase is **94 % of all CPU time** (across 32 cores),
and the single-threaded sim_pcr loop is only 1.5 %.

**Hardware counters** (`perf stat -d`):

| metric | value | reading |
|---|---:|---|
| wall time | 87.6 s | matches the optimized baseline |
| CPU utilization | 15.4× of 32 | parallel populate well-saturated; serial sim_pcr drags average |
| IPC | **1.53** | scalar with stalls; no SIMD |
| **Branch miss rate** | **7.61 %** | high — target is < 3 % for hot loops |
| L1 d-cache miss rate | 0.56 % | very low — hot data fits in L1 |
| LLC counters | not supported | this CPU doesn't expose them |

The branch miss rate is the headline. With 828 B retired branches ×
7.61 % = **63 B mispredicts × ~15 cycles each ≈ 17 % of total cycles
wasted on branch recovery**. That's the dominant inefficiency in
the thal DP. The L1 miss rate of 0.56 % means the parameter tables
and DP matrices fit comfortably in L1; we are not memory-bound.

This explains directly why **`-DTHAL_USE_FLOAT` produced no
speedup**: float would help if we were memory-bound (smaller working
set) or vectorizable (wider SIMD). We are neither — we're stalling
on branch recovery, which is independent of value width.

## Accuracy

Two comparisons matter here: ThalReal=float vs ThalReal=double under
the asymmetric cache (a value-width comparison), and the symmetric
cache vs the asymmetric cache (an orientation-of-thal comparison).
Both are vs the same set of 364,188 numeric cells in
`regression_new.csv`.

**ThalReal=float vs the asymmetric-double baseline:**

| tolerance | matching cells |
|---|---|
| exact     | 169,242 / 364,188 (46.5 %) |
| rel < 1e-9 | 194,013 / 364,188 (53.3 %) |
| rel < 1e-6 | 340,880 / 364,188 (93.6 %) |
| rel < 1e-3 | **364,188 / 364,188 (100 %)** |
| max rel error | **2.63e-4** |

**Symmetric cache vs the asymmetric-double baseline:**

| tolerance | matching cells |
|---|---|
| exact     | 359,535 / 364,188 (98.7 %) |
| rel < 1e-9 | 359,766 / 364,188 (98.8 %) |
| rel < 1e-6 | 361,781 / 364,188 (99.3 %) |
| rel < 1e-3 | 364,153 / 364,188 (99.99 %) |
| max rel error | **4.15e-3** |

The symmetric version's drift is concentrated in a handful of cells
— the worst is in column 11 (total\_nonspec\_r) at row 15627, where
the asymmetric and symmetric runs differ by 0.4 %. That comes from
empirically observing thal asymmetry: `thal(A, B)` differs from
`thal(B, A)` in the 7th significant digit on **41 % of pairs in a
50-address sample**. The values are equal at any reasonable
thermodynamic precision (~3 decimal digits matter for ΔG predictions),
but the per-pair pick of one orientation propagates through to
cumulative differences in some downstream cells. 99.99 % of cells
agree to better than 1e-3 relative.

The outlier diagnostic
([scripts/diagnose_outliers.py](../scripts/diagnose_outliers.py))
returns the same buckets across all three modes (asymmetric double,
asymmetric float, symmetric double): 976 ok, 0 R-stuck, 0 F-stuck,
3 both-fail (addresses 686, 803, 322 — the same three primer-design
issues each time).

---

## Where to look for further speedup

Three options remain after the symmetric cache. Listed in order of
expected wall-time gain at this point in the optimization curve.

### ✅ Done: symmetric cache

Implemented and measured: 87 s → 57 s (**1.52× speedup**), cache
memory halved (245 MB → 122 MB at N = 979 double). Output drift is
0.4 % max relative on 0.01 % of cells; outlier diagnosis unchanged.
See "Wall-time results" above.

### 1. Parallelize the `sim_pcr` loop in `test_eq.cpp`

**What**: post-populate, the `for (i = 0; i < N; i++) sim_pcr(i, ...)`
loop in `test_eq.cpp` runs single-threaded. With the symmetric cache
that loop accounts for ~37 s of the 57 s total wall. With 32
threads it would compress to ~1.2 s.

**How**: replace the loop with a thread-pool pattern like the
existing `eval_addresses_thread` in `address_eval.cpp`. Each thread
takes the next address index from an atomic counter, runs `sim_pcr`,
and appends its row to `regression_new.csv` under `outfile_mtx`. The
`dimer_cache` is read-only after populate so it's freely shared.

**Effect — wall time**:
- sim_pcr loop: ~37 s → ~1.2 s.
- Populate phase unchanged at ~20 s.
- Total: 57 s → **~22 s (≈61 % faster)**.

**Effect — memory**:
- Each worker thread has its own `EQ` (304 B + per-N `address_k_conc_vec`
  ≈ 1.18 MB at N = 979). 32 workers = ~38 MB — small relative to the
  cache.

**Caveats**:
- File output ordering changes (rows interleave by address index in
  whatever order threads finish). Already true in the existing
  `eval_addresses_thread` path. Fine for analysis as long as you
  don't rely on cycle-by-cycle ordering matching address index.
- Thread-safety in `sim_pcr` itself: each call constructs its own EQ
  and writes only to that EQ + the shared output file (under
  mutex). No issue.

**Risk level**: low. Standard pattern, already proven by
`eval_addresses_thread`.

### 2. Tighten the cache to entries `sim_pcr` actually reads

**What**: the symmetric cache holds entries for every kind pair
`(ki, kj) ∈ [0..4) × [0..4)`. But `sim_pcr` only reads first-kinds
`ki ∈ {0=f, 1=f_rc, 2=r}` — the `target.r_rc` first-arg is never
queried. That's 25 % of the populated entries that no caller reads.

**How**: parameterize the cache slot space by "useful kinds" instead
of all four. The canonical-pair index becomes more involved (the
useful-kind set isn't a clean prefix of the 4-kind set), but is
straightforward to implement.

**Effect — wall time**:
- Populate work drops ~25 %: 20 s → ~15 s.
- sim_pcr unchanged.
- Total: 57 s → **~52 s (≈9 % faster)** standalone, or in combination
  with #1, 22 s → ~17 s.

**Effect — memory**:
- Symmetric cache drops ~25 %: 122 MB → ~92 MB (double); 61 MB → 46 MB
  (float).

**Caveats**:
- None on accuracy — we're literally not computing entries no one
  reads. Output bit-identical to the current symmetric cache.
- Indexing is more fiddly than the symmetric formula already in
  place. Easy to get the canonical normalization wrong.

**Risk level**: zero accuracy risk; modest implementation risk
(indexing).

### 3. Reduce branch mispredictions inside thal's DP

**What**: 17 % of all cycles in this run are stalled on mispredicted
branch recovery (63 B mispredicts × ~15 cyc each, out of 5.6 T total
cycles). The mispredicts are concentrated in `calc_bulge_internal_dimer`
and `fillMatrix_dimer`'s inner loops, which evaluate per-cell DP
recurrences with data-dependent control flow.

**How** (in rough increasing order of effort):

(a) **Add `__builtin_expect` hints** on the most common branches.
Cheap, but only helps when the bias is strong and consistent.
Probably 2-4 % wall-time win; depends on whether the branches are
biased or genuinely 50/50.

(b) **Replace conditionals with `std::min`/`std::max`/`fmin`/`fmax`
arithmetic** where the DP is computing element-wise min/max over
candidates. CPUs implement these as branchless `cmov`/`vminps`
instructions. Eliminates the misprediction entirely for those sites.
Likely 5-10 % wall-time win.

(c) **Restructure the DP to be branchless / vectorizable**: pad the
recurrence to fixed length (no inner loop break), use SIMD over the
"j" inner axis. Significant rewrite; would also lift IPC from 1.53
toward ~3-4 (vectorized FP). Potential win: 30-50 % wall-time.

**Effect — wall time** (in the populate phase, where this matters):
- Option (a): 87 s → ~84 s.
- Option (b): 87 s → ~80 s.
- Option (c): 87 s → ~50-60 s.

**Caveats**:
- This is upstream primer3 code. Unless we fork primer3, every
  upstream rebase has to re-apply the changes.
- Branchless rewrites of DP code are notoriously easy to break.
  Need a regression test that's tighter than 9-digit print
  precision; the inner DP scores need bit-equal preservation or we
  risk picking different alignments.
- Only the populate phase benefits — `sim_pcr` is already cache
  lookups. Wall-time impact is bounded by the populate fraction
  (~45 % of total wall, post-caching).

**Risk level**: high for option (c); moderate for (a)/(b). Most
work, most fragile, vendored code, not all wins are guaranteed.

### Composability and recommendation

Combined ceiling, starting from the current symmetric cache (57 s):

| variant | wall (est.) | from baseline 14m18s | from current 57 s |
|---|---:|---:|---:|
| current (symmetric cache + parallel populate) | 57 s | 15.0× | — |
| + #1 (parallel sim_pcr loop)                  | ~22 s | 39× | 2.6× |
| + #1 + #2 (tighten kinds)                     | ~17 s | 50× | 3.4× |
| + #1 + #2 + #3(b) (branchless DP min/max)     | ~14 s | 61× | 4.1× |
| + #1 + #2 + #3(c) (vectorized DP)             | ~9 s | 95× | 6.3× |

Practical recommendation:

- **#1 first** — biggest single win at this point in the curve
  (~2.6×). Standard threading pattern, infrastructure already exists
  in `eval_addresses_thread`, no accuracy concerns.
- **#2 second** — small win on top of #1, no accuracy risk, gives
  another ~25 % memory headroom for larger N.
- **#3 only if you really need it** — primer3 is upstream code and
  rewriting its DP is a fork-maintenance burden. The 17 % cycles
  wasted on branch mispredicts are real but the engineering cost is
  significant. Worth it only if the workflow runs frequently enough
  that minutes saved over months exceed the maintenance cost.

## Caveats on the perf data

- Single-machine measurement on the 32-core dev box. Mileage will
  vary on different CPUs (branch predictor quality, cache sizes).
- `perf record --call-graph dwarf` produced ~10 GB of `perf.data` for
  this run. Use `--call-graph fp` if disk I/O matters; results are
  similar but slightly less reliable for inlined functions.
- `perf stat`'s `<not supported>` lines for LLC events on this CPU
  mean we couldn't directly measure cache thrashing past L1. Given
  the L1 miss rate of 0.56 % is already very low, L2/L3 wouldn't
  have changed the picture.
