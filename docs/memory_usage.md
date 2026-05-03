# Memory usage as a function of address count

How RAM scales with `N = addresses.size()` in the active workload
(`test_eq.cpp` → `populate_dimer_cache` → `sim_pcr`). Two precision
knobs change the picture:

- `-DPSIM_USE_LONG_DOUBLE`: `Real` becomes 16-byte long double (default
  is `double`).
- `-DTHAL_USE_FLOAT`: `ThalReal` becomes 4-byte float (default is
  `double`). Affects every `dhds_*` field in `address_k_conc` and
  every entry in `dimer_cache`.

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
  2.  EQ::address_k_conc_vec  1184 B (default) per N (1056 with THAL_FLOAT)
  3.  dimer_cache (symmetric / triangular):
                            8 N² × {16|8} B  →  128 N² B (default) or 64 N² B (float)
```

The cache (3) is the new piece — it scales as **O(N²)**, the others
as O(N). For modest N it's a small constant; past a few hundred
addresses it dominates.

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

`populate_dimer_cache` precomputes thal results for the
**canonical** half of every `(addr, kind, addr', kind')` pair —
`thal(A, B)` is approximately symmetric in its arguments, so storing
both `(A, B)` and `(B, A)` is wasteful. With `M = 4·N` "sequence
slots" (each address has 4 kinds: f, f_rc, r, r_rc), the cache
holds the `M·(M+1)/2 ≈ 8·N²` canonical pairs `(a, b)` with `a ≤ b`.
Lookups normalize argument order before indexing.

Each entry is `{ThalReal dh, ThalReal ds}` — 16 B default, 8 B with
`-DTHAL_USE_FLOAT`. Total cache bytes:

| flag | per-entry | per N entries | bytes / N² |
|---|---:|---:|---:|
| (default)            | 16 B | 8·N² + 2·N | **128 N²** |
| `-DTHAL_USE_FLOAT`   | 8 B  | 8·N² + 2·N | **64 N²**  |

For specific N:

| N | dimer_cache (double) | dimer_cache (float) |
|---:|---:|---:|
| 100      | 1.28 MB    | 642 kB     |
| 500      | 32 MB      | 16 MB      |
| 1 000    | 128 MB     | 64 MB      |
| **2 000**| **512 MB** | **256 MB** |
| 5 000    | 3.2 GB     | 1.6 GB     |
| 10 000   | 12.8 GB    | 6.4 GB     |
| 100 000  | 1.28 TB ⚠  | 640 GB ⚠   |

The pre-symmetry cache was 2× larger (covered both `(A, B)` and
`(B, A)`) — at N=979 we measured **122 MB** symmetric vs **245 MB**
asymmetric.

## Total memory (single-threaded `sim_pcr` path)

`test_eq.cpp` runs sequentially: one `EQ` is alive at a time,
`addresses` is built once and shared, `dimer_cache` is built once
and shared. Total:

```
  total ≈  ~150 × N    (addresses)
       +  ~1184 × N    (1 EQ × address_k_conc_vec)
       +  128 × N²     (dimer_cache, double mode)
```

| N | addresses | EQ + per-N | dimer_cache | **total (double)** | total (float) |
|---:|---:|---:|---:|---:|---:|
| 100      | 14 kB   | 119 kB   | 1.28 MB   | **1.41 MB**  | 775 kB |
| 500      | 68 kB   | 593 kB   | 32 MB     | **32.7 MB**  | 16.7 MB |
| 1 000    | 136 kB  | 1.18 MB  | 128 MB    | **129 MB**   | 65 MB |
| 2 000    | 272 kB  | 2.37 MB  | 512 MB    | **515 MB**   | 259 MB |
| 5 000    | 680 kB  | 5.92 MB  | 3.2 GB    | **3.21 GB**  | 1.61 GB |
| 10 000   | 1.36 MB | 11.8 MB  | 12.8 GB   | **12.8 GB**  | 6.4 GB |
| 100 000  | 13.6 MB | 118 MB   | 1.28 TB ⚠ | **1.28 TB** ⚠ | 640 GB ⚠ |

The cache exceeds everything else combined past **N ≈ 100** (already
~90 % of total at N = 1000).

Practical ceilings on a typical workstation:
- **8 GB RAM**: comfortable up to N ≈ 1500 (double) or 3000 (float).
- **16 GB RAM**: up to N ≈ 3500 (double) or 5000 (float).
- **64 GB RAM**: up to N ≈ 7000 (double) or 10000 (float).
- **256 GB RAM** (workstation / server): up to N ≈ 14000 (double) or
  20000 (float).
- **N = 100 000**: not feasible on any reasonable machine. Need to
  drop the cache.

`-DTHAL_USE_FLOAT` halves the cache memory and buys roughly √2 ≈ 1.4×
headroom on N (since the cache is O(N²), halving memory doubles the
N² budget, which is √2 in N). On the 979-address workload tested it
also produced no measurable accuracy degradation in the regression
output (max relative error 2.6e-4, all within 1e-3 rel; same outlier
diagnosis as double).

If you want to scale past these ceilings without dropping the cache,
see the alternatives in [docs/performance.md](performance.md):
on-demand sparse caching, parallel `sim_pcr`-without-cache, or
restricting the cache to the kinds `sim_pcr` actually reads (saves
~25 % more).

## Multi-threaded path (`eval_addresses_thread`, num_cpu = 32)

This path isn't currently exercised by `test_eq.cpp`, but
`Primeanneal::eval_addresses` exists and a parallel `sim_pcr` loop
(suggestion #1 in performance.md) would use the same pattern. Each
worker thread holds its own `EQ` (with its own `address_k_conc_vec`),
but the `addresses` and `dimer_cache` are shared. Total:

```
  shared addresses     ≈   150 × N
  shared dimer_cache   ≈ 128 × N²  (double) or 64 × N² (float)
  num_cpu × per-EQ     ≈ 32 × (304 + 1184 × N)
```

| N | shared (cache + addrs) | 32× EQ | total (double) |
|---:|---:|---:|---:|
| 100      | 1.30 MB   | 3.81 MB   | **5.1 MB** |
| 1 000    | 128 MB    | 38 MB     | **166 MB** |
| 2 000    | 513 MB    | 76 MB     | **589 MB** |
| 5 000    | 3.2 GB    | 190 MB    | **3.4 GB** |
| 10 000   | 12.8 GB   | 380 MB    | **13.2 GB** |
| 100 000  | 1.28 TB ⚠ | 3.8 GB    | **1.28 TB** ⚠ |

The 32×-EQ multiplier matters less than it used to; the cache
dominates above N ≈ 200. Useful corollary: if you're spinning up
parallel `sim_pcr` workers on a machine that can already hold the
cache, the per-thread EQ overhead is essentially noise.

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

For small N (< 500), memory is dominated by the per-N pieces and is
trivial. **Above N ≈ 500, the symmetric `dimer_cache` becomes the
dominant storage** and grows quadratically — N = 2 k uses ~515 MB,
N = 5 k uses ~3.2 GB, N = 10 k uses ~12.8 GB. Above N ≈ 14 k on a
256 GB workstation you have to drop the cache.

`-DTHAL_USE_FLOAT` cuts the cache memory in half (and shrinks
`address_k_conc` by 11 %) for free, with no measurable accuracy
loss in the regression at 9 print digits. Buys ~1.4× headroom on N.

If you ever need to scale past these ceilings, the cache strategy
needs to change — see the alternatives discussed in
[docs/performance.md](performance.md).
