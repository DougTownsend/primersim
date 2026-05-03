# Performance — calc_dimer cache + parallel populate

This experiment caches every `calc_dimer` result that any `sim_pcr`
might need, computes them in parallel across `num_cpu` threads at
startup, and replaces the per-`sim_pcr` `thal` calls with cache
lookups.

## Setup

The cache is `Primeanneal::dimer_cache`, a flat
`std::vector<dh_ds_cache_entry>` of size `16 * N * N` indexed as

```
dimer_cache[((i * 4 + ki) * N + j) * 4 + kj]
```

where `(i, ki)` and `(j, kj)` are `(address index, sequence kind)` and
`kind ∈ {0=f, 1=f_rc, 2=r, 3=r_rc}`. `populate_dimer_cache` divides
the outer `i` index across threads — each thread writes its own slice
of the cache, no mutex needed. `sim_pcr` reads only.

`test_eq.cpp` calls `pa.populate_dimer_cache(...)` once after
`read_addresses` and before the simulation loop.

## Wall-time results (N = 979)

| build | wall | speedup | regression md5 |
|---|---:|---:|---|
| baseline (no cache, no `-DTHAL_USE_FLOAT`)        | 14 m 18 s | — | `0b8ca853…` |
| **cache + parallel populate, ThalReal=double**    | **1 m 27 s** | **9.86×** | `0b8ca853…` (matches) |
| cache + parallel populate, ThalReal=float         | 1 m 27 s     | 9.86×    | `8a927c84…` (matches prior float baseline) |

Output md5 in `double` mode is bit-identical to the pre-cache
baseline. The float output matches the previous float experiment —
the cache rearranges the thal work but doesn't change values.

CPU time vs wall time tells the story of where time is spent:

|           | wall   | user  | parallelism |
|---|---:|---:|---:|
| double    | 87 s   | 1295 s | 14.9× (avg) |
| float     | 87 s   | 1334 s | 15.3× (avg) |

The benchmark has two phases:

| phase | wall (estimated) | CPU usage |
|---|---:|---|
| `populate_dimer_cache` | ~39 s | 32 threads, near-saturating |
| `sim_pcr` loop (979 calls) | ~48 s | 1 thread, cache lookups only |

So the previously-dominant 14 m of single-threaded thal compresses
to ~39 s spread across 32 cores, and the lookup-only `sim_pcr` loop
runs in another ~48 s.

## Profile (gprof, `-O2 -pg`)

`gprof`'s call counts and per-call times are unreliable in
multi-threaded code (the `mcount` counters race; thread time is
sampled inconsistently). Treat the *shape* of the profile as
informative; ignore absolute counts. With that caveat:

```
ThalReal = double                 ThalReal = float
─────────────────────────────     ─────────────────────────────
65.2% thal()                      64.2% thal()
29.2% calc_bulge_internal_dimer   30.0% calc_bulge_internal_dimer
 1.96% RSH                          2.09% RSH
 0.92% LSH                          0.97% LSH
 0.75% EQ::calc_cx                  0.67% EQ::calc_cx
 0.64% calc_strand_bindings         0.65% calc_strand_bindings
 0.44% update_strand_concs          0.44% update_strand_concs
 0.43% EQ::calc_bound_concs         0.48% EQ::calc_bound_concs
 0.29% sim_pcr (self)               0.25% sim_pcr (self)
 0.03% populate_dimer_cache lambda  0.04% populate_dimer_cache lambda
```

The shape is essentially unchanged from the pre-cache profile: thal
still dominates total CPU time. What's different is the **wall-time**
distribution — those CPU seconds are now spread across 32 cores
instead of 1, so wall time drops 10×.

The float profile is essentially identical to double. `thal()`
internal arithmetic does run on `float` (`-DTHAL_USE_FLOAT` switches
the `ThalReal` typedef), but on x86-64 the underlying scalar SSE2
ops are similar speed for `float` and `double` — the dimer DP is
branchy and table-lookup heavy, not a vectorizable arithmetic loop.

`-pg` profile binaries ran ~12 m in this experiment, vs 1m27s for
the optimized binary. Almost all of that overhead is in the
populate phase: gprof's per-function instrumentation fires on every
thal internal call (billions across 32 threads), and the contention
on the global `mcount` counter serializes much of it. **Use the `time`
output, not the gprof totals, to compare wall-clock performance.**

## Accuracy

ThalReal=float vs the double baseline (compared cell-by-cell across
all 364,188 numeric cells in `regression_new.csv`):

| tolerance | matching cells |
|---|---|
| exact | 169,242 / 364,188 (46.5 %) |
| rel < 1e-9 | 194,013 / 364,188 (53.3 %) |
| rel < 1e-6 | 340,880 / 364,188 (93.6 %) |
| rel < 1e-3 | **364,188 / 364,188 (100 %)** |
| max rel error | **2.63e-4** |

The outlier diagnostic
([scripts/diagnose_outliers.py](../scripts/diagnose_outliers.py))
returns the same buckets in both modes: 976 ok, 0 R-stuck, 0 F-stuck,
3 both-fail (the same three primer-design issues at addresses 686,
803, 322).

## Caveats

- The cache populate is now the dominant phase (~45 % of wall time).
  At larger N the populate scales as O(N²) thal calls — the wall
  cost grows quadratically.
- The cache memory grows O(N²) too. See
  [docs/memory_usage.md](memory_usage.md). Up to N ≈ 2 k it's a
  ~1 GB allocation; at N = 10 k it's 25 GB (double) and prohibitive.
- For N ≫ 1000, alternatives:
  1. Drop the cache entirely and parallelize the per-`sim_pcr`
     setup loop across threads. Same total thal work but better
     scaling.
  2. Exploit thal's near-symmetry — `thal(A, B)` differs from
     `thal(B, A)` only at the 7th significant digit (41 % of pairs
     in a 50-address sample). A symmetric cache halves the
     populate work and the cache memory at the cost of sub-9-digit
     precision drift.
  3. Move `populate_dimer_cache` onto a sparse / on-demand model
     — only fill entries that a sim_pcr actually requests. With
     `M < N` sim_pcr calls each touching 12·N pairs, the active
     set is `12·M·N`, smaller than the full `16·N²` when `M < N`.

## Recommendation

For workloads up to **N ≈ 1000-2000 addresses on a 32-core machine**,
the cache + parallel populate is a clear win: ~10× wall-time
speedup, no accuracy change, modest memory cost.

Above that, the O(N²) cache becomes the bottleneck and a different
strategy (parallel sim_pcrs, lazy/sparse cache, or symmetric cache)
is needed.

ThalReal=float vs ThalReal=double has no meaningful wall-time
impact under this caching strategy. If you want the smaller cache
footprint (122 MB vs 245 MB at N = 979), float is fine; otherwise
keep the default double.
