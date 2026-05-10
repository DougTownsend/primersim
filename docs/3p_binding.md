# 3' End Binding Cache + Amplifying / Non-Amplifying Split

This is a design doc for what I think `todo.md` item 2 is asking for.
Read it first; correct anything that's wrong, then I'll implement.

## What's being added

A second thal cache `dimer_3p_cache` alongside the existing
`dimer_cache`. Both are keyed identically — by `(primer_index,
kind ∈ {seq, rc})` over `M = 2 * primer_pool.size()` slots — but
they store different thal results:

| cache | thal type | storage |
|---|---|---|
| `dimer_cache` (existing) | `thal_any` (default) | symmetric, triangular `M*(M+1)/2` |
| `dimer_3p_cache` (new) | `thal_end1` | asymmetric, square `M * M` |

Both are populated once after `read_primer_pool` and survive any
pairing reshuffle.

### Why end1, not end2

`thal_end1(seq1, seq2)` anchors the duplex so that **seq1's 3' end**
is at the duplex end. `thal_end2(seq1, seq2)` anchors seq2's 3' end
instead. Whenever the cache is queried, the primer (the one whose 3'
end can extend if anchored) goes in the seq1 position and the strand
goes in seq2; calling end1 mode gives the "primer's 3' end is
anchored" alignment energy directly.

### Why square instead of triangular

`thal_any(A, B) ≈ thal_any(B, A)` (symmetric), so `dimer_cache`
stores only canonical `a ≤ b`. By contrast,
`thal_end1(A, B) ≠ thal_end1(B, A)` in general — different sequence
plays the "3'-end-anchored" role. The new cache must store both
orderings, so it's `M × M`.

Storage cost (16 bytes/entry, ThalReal = double):

| pool size P | M = 2P | square cache size |
|---|---|---|
| 100  | 200  | ~640 KB |
| 979  | 1958 | ~61 MB |
| 1958 | 3916 | ~245 MB |

### Maximum primer pool size that fits in a given memory budget

Both caches combined (triangular + square ≈ `(3/2) * M² * 16` bytes
≈ `96 * P²` bytes with `ThalReal = double`):

| memory budget | max P (double) | max P (float, `-DTHAL_USE_FLOAT`) |
|---:|---:|---:|
| 1 GB   | 3,228   | 4,565   |
| 4 GB   | 6,455   | 9,129   |
| 16 GB  | 12,910  | 18,257  |
| 64 GB  | 25,820  | 36,514  |
| 128 GB | 36,514  | 51,640  |
| 256 GB | 51,640  | 73,028  |
| 500 GB | 72,168  | 102,062 |

Float halves per-entry storage to 8 bytes (`12 * M²` total ≈ `48 * P²`),
so the practical primer cap rises by a factor of √2 with `ThalReal =
float`. At 500 GB and float precision the cap is around 100 K primers.

## Where it's used in `sim_pcr`

`calc_strand_bindings` / `update_strand_concs` walk over
`end5 × end3` binding configurations. The 25 states encode where on
the strand the primer is hybridized:

- 5'-end-of-strand bindings: primer's 3' end is **not** at the strand's
  3' end → polymerase has nothing to extend onto the address → these
  are **always non-amplifying**. Behavior is unchanged: just look up
  the total interaction energy from `dimer_cache` (already happens
  today).

- 3'-end-of-strand bindings: primer's 3' end **is** at the strand's
  3' end → polymerase can extend the primer using the strand as
  template → these can be amplifying *or* not.
  - **Amplifying portion** = `dimer_3p_cache[primer, strand]`
    (thal_end1 result with primer as seq1).
  - **Non-amplifying portion** = `dimer_cache[primer, strand] −
    dimer_3p_cache[primer, strand]` (subtract the amplifying alignment
    from the total ensemble to get what's left).

## Implementation sketch

1. **Header** — add to `Primeanneal`:
   ```cpp
   std::vector<dh_ds_cache_entry> dimer_3p_cache;
   ```
   Layout: `dimer_3p_cache[a * M + b]`, `a` and `b` both `(primer_idx
   * 2 + kind)`. No canonical-ordering swap.

2. **populate_dimer_cache** — extend (or add a sibling
   `populate_dimer_3p_cache`) to fill the square cache via
   `thal(seq_a, seq_b, type = thal_end1)` for every ordered pair.
   Same threading pattern as the existing populate (cyclic by `a`,
   `num_cpu` threads).

3. **sim_pcr cached() lambda** — keep the existing one for
   `dimer_cache`. Add a parallel lookup for `dimer_3p_cache`:
   - takes `(a_pair, ka, b_pair, kb)`
   - resolves to `(primer_idx, kind)` for each side using the same
     pairing/RC logic as today
   - indexes `dimer_3p_cache[a * M + b]` directly (no swap)

4. **3'-end binding split** — in the section of
   `calc_strand_bindings` / `update_strand_concs` that handles the
   "primer at strand's 3' end" cases:
   - Get total dh/ds from the existing cache lookup.
   - Get 3'-end dh/ds from the new cache.
   - Compute non-amplifying dh/ds = total − 3'-end (element-wise on dh
     and ds — see "things I want to confirm" below).
   - Use the amplifying part for amplification accounting (target
     and possibly mispriming product growth) and the non-amplifying
     part the way today's code uses the total interaction.

## Design decisions (confirmed)

1. **Primer goes in seq1**: the primer (whichever the simulator is
   treating as F or R for the active pair) is seq1 / index `a` in the
   cache; the strand is seq2 / index `b`. `thal_end1` anchors the
   primer's 3' end.

2. **Subtraction is element-wise on dh and ds**:
   ```
   dh_nonamp = dh_total - dh_3p
   ds_nonamp = ds_total - ds_3p
   ```
   Subtracting dh and ds separately is exactly equivalent to
   subtracting the corresponding ΔG values at any single temperature
   (ΔG = ΔH - TΔS is linear), so feeding `(dh_nonamp, ds_nonamp)`
   through `dhds_to_eq_const` gives the same K as the equivalent ΔG
   subtraction. No information is lost vs. doing the math in ΔG
   space; element-wise on (dh, ds) is the natural form because the
   downstream code already wants those.

3. **Special case when amplifying ≈ total**: the amplifying
   alignment is a subset of the total ensemble (not strictly so —
   when the most-stable alignment IS the 3'-end-anchored one, total
   = amp). When `(dh_total, ds_total) == (dh_3p, ds_3p)` (within
   floating-point tolerance), treat the binding as 100% amplifying:
   **do not** add a non-amplifying contribution at all. Don't model
   it as a 0-ΔG (`K = 1`) binding — that's a real interaction with
   zero preference, which would corrupt the simulation. Concretely:
   skip the non-amplifying term when both deltas vanish.

4. **Memory cost is acceptable.** See the budget table above.

Implementing now.
