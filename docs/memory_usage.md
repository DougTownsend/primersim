# Memory usage as a function of address count

How RAM scales with `N = addresses.size()` in the active workload
(`test_eq.cpp` → `populate_dimer_cache` → `sim_pcr`). Two precision
knobs change the picture:

- `-DPSIM_USE_LONG_DOUBLE`: `Real` becomes 16-byte long double (default
  is `double`).
- `-DTHAL_USE_FLOAT`: `ThalReal` becomes 4-byte float (default is
  `double`). Affects every `dhds_*` field in `address_k_conc` and
  every entry in the `dimer_cache`.

All sizes are measured with `sizeof()` on x86-64 Linux, gcc -O3:

| flags | `Real` | `address_k_conc` | `EQ` minus vec | `dh_ds_cache_entry` |
|---|---:|---:|---:|---:|
| (default)                                  | 8  | **1184** | 304 | 16 |
| `-DTHAL_USE_FLOAT`                         | 8  | **1056** | 304 | 8  |
| `-DPSIM_USE_LONG_DOUBLE`                   | 16 | **2112** | 592 | 16 |
| `-DPSIM_USE_LONG_DOUBLE -DTHAL_USE_FLOAT`  | 16 | **1984** | 592 | 8  |

## Three pieces of N-dependent storage

```
  1.  addresses             40 B + 4 strings × ~24 B  ≈    150 B per N
  2.  EQ::address_k_conc_vec  1184 B (default), 1056 B (THAL_FLOAT) per N
  3.  dimer_cache           16 × N entries × {8|16} B per N²
```

The cache (3) is the new piece — it scales as **O(N²)**, the others
as O(N). For modest N it's a small constant; past N ≈ 1000 it
dominates.

### Per-address breakdown (the O(N) parts)

```
  addresses[i]         40 B   (struct: 4 char*, 1 double)
  4 string allocations ~96 B  (4 × ~24 B for 21-byte strings + alloc overhead)
  address_k_conc_vec[i] 1184 B / 1056 B (default / THAL_FLOAT)
                       -----
                       ~1320 / ~1190 B per address
```

`address_k_conc` is dominated by the four 5×5 `Real` matrices —
`fstrand`, `rstrand`, `fstrand_change`, `rstrand_change` — at 200 B
each, plus the dhds and per-cycle K-cache fields.

### `dimer_cache` (the O(N²) part)

`populate_dimer_cache` precomputes every `(i, ki, j, kj)` thal result
so `sim_pcr` never calls thal. The cache holds `16 N²` entries of
`{ThalReal dh, ThalReal ds}`, so:

| flag | per-entry | per-N² |
|---|---:|---:|
| (default)            | 16 B | **16 N²** = 16 N² B |
| `-DTHAL_USE_FLOAT`   | 8 B  | **8 N²** = 8 N² B |

For specific N:

| N | dimer_cache (double) | dimer_cache (float) |
|---:|---:|---:|
| 100      | 160 kB    | 80 kB     |
| 1 000    | 16 MB     | 8 MB      |
| **2 000**| **64 MB** | **32 MB** |
| 10 000   | 1.6 GB    | 800 MB    |
| 100 000  | 160 GB    | 80 GB     |

## Total memory (single-threaded `sim_pcr` path)

`test_eq.cpp` runs sequentially: one `EQ` is alive at a time,
`addresses` is built once and shared, `dimer_cache` is built once
and shared. Total:

```
  total ≈  ~150 × N    (addresses)
       +  ~1184 × N    (1 EQ × address_k_conc_vec)
       +  16 × N²      (dimer_cache, double mode)
```

| N | addresses | EQ + per-N | dimer_cache | **total (double)** | total (float) |
|---:|---:|---:|---:|---:|---:|
| 100      | 14 kB   | 119 kB   | 160 kB    | **293 kB**   | 213 kB |
| 1 000    | 136 kB  | 1.18 MB  | 16 MB     | **17.4 MB**  | 9.4 MB |
| 2 000    | 272 kB  | 2.37 MB  | 64 MB     | **66.6 MB**  | 34.6 MB |
| 5 000    | 680 kB  | 5.92 MB  | 400 MB    | **407 MB**   | 206 MB |
| 10 000   | 1.36 MB | 11.8 MB  | 1.6 GB    | **1.61 GB**  | 813 MB |
| 100 000  | 13.6 MB | 118 MB   | 160 GB    | **160 GB** ⚠ | 80 GB ⚠ |

Crossover where the cache exceeds everything else combined: around
**N = 100** (the cache is already 80 % of total at N = 1000).

Beyond N ≈ 5 k the double-mode cache pushes past 1 GB; at N = 10 k
it hits 1.6 GB and at N = 100 k it would need 160 GB — clearly
infeasible on a workstation. `-DTHAL_USE_FLOAT` halves the cache
size and buys roughly 2× headroom on N.

If you want to scale past N ≈ 5 k without dropping the cache, see
the alternatives in [docs/performance.md](performance.md):
on-demand sparse caching, parallel `sim_pcr`-without-cache, or
exploiting thal's near-symmetry to halve the cache.

## Multi-threaded path (`eval_addresses_thread`, num_cpu = 32)

This path isn't currently exercised by `test_eq.cpp`, but
`Primeanneal::eval_addresses` exists. Each worker thread holds its
own `EQ` (with its own `address_k_conc_vec`), but the `addresses`
and `dimer_cache` are shared. Total:

```
  shared addresses     ≈    150 × N
  shared dimer_cache   ≈ 16 × N²  (double) or 8 × N² (float)
  num_cpu × per-EQ     ≈ 32 × (304 + 1184 × N)
```

| N | shared (cache + addrs) | 32× EQ | total (double) |
|---:|---:|---:|---:|
| 100      | 174 kB    | 3.8 MB    | **4.0 MB** |
| 1 000    | 16.1 MB   | 38 MB     | **54 MB** |
| 10 000   | 1.6 GB    | 380 MB    | **2.0 GB** |
| 100 000  | 160 GB ⚠  | 3.8 GB    | **164 GB** ⚠ |

The cache dominates above N = 1 k. The 32× EQ multiplier matters
less than it used to.

## On-disk artifacts

`regression_new.csv` is appended row-by-row, only one line buffered
at a time:

| N | regression_new.csv |
|---:|---:|
| 100 | ~500 kB |
| 1 000 | ~5 MB |
| 10 000 | ~50 MB |
| 100 000 | ~500 MB |

`addresses.csv` is ~50 B per address.

## TL;DR

For small N (< 1 k), memory is dominated by the per-N pieces and is
trivial. **Above N ≈ 1 k, the new `dimer_cache` becomes the dominant
storage** and grows quadratically — N ≈ 5 k uses ~400 MB,
N = 10 k uses 1.6 GB, N = 100 k is 160 GB and not feasible.
`-DTHAL_USE_FLOAT` cuts the cache memory in half (and shrinks
`address_k_conc` by 11 %) for free, with no measurable accuracy
loss in the regression at 9 print digits.

If you ever need to scale past ~5 k addresses, the cache strategy
needs to change — see the alternatives discussed in
[docs/performance.md](performance.md).
