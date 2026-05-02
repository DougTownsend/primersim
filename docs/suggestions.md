# Code review: eq.cpp / eq.hpp / test_eq.cpp

Suggestions for readability, performance, and accuracy. Each item cites
file:line and notes whether it's a bug, a cleanup, or a tradeoff.

---

## Correctness bugs (fix first)

### 1. `read_addresses` allocates `r` and `r_rc` with the length of `f`
[src/eq.cpp:406-409](src/eq.cpp#L406-L409)

```cpp
a.f    = (char *)malloc(strlen(tmp_f)+1);
a.r    = (char *)malloc(strlen(tmp_f)+1);  // ← uses tmp_f length
a.f_rc = (char *)malloc(strlen(tmp_f)+1);
a.r_rc = (char *)malloc(strlen(tmp_f)+1);
```

If `tmp_r` is longer than `tmp_f` (different primer lengths), the next
`strcpy(a.r, tmp_r)` writes past the buffer. Today both inputs happen to be
20 nt so it's latent. Fix: size each buffer with its own `strlen`.

### 2. `evaluate_addresses` passes `double` as `bool`
[src/eq.cpp:677](src/eq.cpp#L677)

```cpp
if(addresses.size() == 0)
    read_addresses(in_filename, dna_conc);   // bool with_temp_c ← gets dna_conc
```

`dna_conc` (3e-15) is non-zero so it converts to `true`, but the file
format being parsed has no temp column. Almost certainly meant `false`.

### 3. `eval_thread` reads uninitialized `eq.tmp[1]`
[src/eq.cpp:645-650](src/eq.cpp#L645-L650)

```cpp
eq.tmp[0].add(eq.c[FX], eq.c[RX]);
//eq.tmp[1].add(eq.c[FY], eq.c[RY]);   ← only assignment, commented out
eq.tmp[2].min(eq.tmp[0], eq.tmp[1]);   ← uses tmp[1]
...
eq.tmp[0].max(eq.tmp[0], eq.tmp[1]);
```

The Y-branch was commented out when Y/B/A were dropped, but `tmp[1]` is
still used by the `min`/`max` calls. The result is whatever lingered in
`tmp[1]` from a previous iteration. Either set `tmp[1] = 0` (degrades to
`min(tmp[0], 0) = 0`, which makes the lin/exp split degenerate) or rewrite
without `min`/`max` now that there's only one nonspec pool.

### 4. `read_addresses` leaks on second call
[src/eq.cpp:393](src/eq.cpp#L393)

`addresses.clear()` drops the vector but doesn't `free` the four `char *`
inside each `address`. The file's TODO acknowledges this. Add a per-entry
free loop before `clear()`, or switch the `char *` quartet to
`std::string`.

### 5. Outer-loop convergence check is bit-exact equality
[src/eq.cpp:293](src/eq.cpp#L293)

```cpp
if (!c[F].cmp(last_val[F]) && !c[R].cmp(last_val[R])) break;
```

`mpfr_cmp` returns 0 only when the values are bitwise identical. With
arithmetic in the inner Newton step, that's rare — outer convergence
relies on the inner solver's "residual stops decreasing" exit producing
a stable fixed point. Today this works (regression matches), but it's
fragile. Replace with `|c[F]_new - c[F]_old| <= eps * c[F]_new`, where
`eps` is something like 2^-(FLOAT_PREC-8).

---

## Dead code (delete)

After replacing the binary-search solver with Newton, the following are
unreferenced:

- `EQ::solve_cf_cr` and its internals: [src/eq.cpp:133-219](src/eq.cpp#L133-L219)
- `EQ::calc_cf` / `EQ::calc_cr`: [src/eq.cpp:103-110](src/eq.cpp#L103-L110)
- `EQ::bounds[3]` / `EQ::solutions[3]`: [include/eq.hpp:328-329](include/eq.hpp#L328-L329)
- `solve.cpp` at the repo root — it was a scratch/handoff file and has
  been merged into `src/eq.cpp`. Either delete or move under `docs/`
  with a clear "reference, not built" comment.
- The commented-out A/B/Y blocks in `eval_thread`
  ([eq.cpp:577-582, 591-594, 622-628, 632, 638-639, 645](src/eq.cpp))
  and the parallel commented enum block in
  [include/eq.hpp:65-84](include/eq.hpp#L65-L84) and
  [include/eq.hpp:92-106](include/eq.hpp#L92-L106).
- Unused EQ members: `spec_fwd_amp`, `spec_rev_amp`, `nonspec_exp_amp`,
  `nonspec_lin_amp`, `best_*` family, `saved_frc_conc`, `saved_rrc_conc`,
  `avg_nonspec_amp`, `last_nonspec_*` — most only set by `eval_thread`,
  which is itself the obsolete address-evaluator. Same for
  `address_k_conc::dhds[2][4]` (only `tmp_dhds` is read).
- Doubled `address_index = 0` in
  [src/eq.cpp:681-682](src/eq.cpp#L681-L682) (also at 1019-1020).

---

## Index-array hygiene

### 6. Right-size the c, c0, k arrays
[include/eq.hpp:331-333](include/eq.hpp#L331-L333)

`Psim_f c[19]; c0[6]; k[13];` — the active enum only uses 10 / 3 / 7 of
those slots. With FLOAT_PREC=1024 each unused Psim_f is ~144 bytes plus
an mpfr_init/clear pair per `EQ` instance. Resize to match the live
enums and the sized `address_k_conc_vec` (one EQ per thread, but tens of
thousands of `address_k_conc` per simulation — that's where the savings
add up).

### 7. Replace bare ints with a typed enum
[include/eq.hpp:54-91](include/eq.hpp#L54-L91)

`const int F = 0;` etc. work but provide no type-checking. A fwd-decl
`enum Species : int { F=0, R=1, X=2, ... };` and `enum Rate : int { K_FH=0, ... };`
would catch `c[K_FF]`-style index mixups at compile time and let you
`#include` the header without polluting the namespace.

---

## Psim_f wrapper hardening

The wrapper has a long list of footguns documented in `CONTEXT.md` but
not in the header itself. Most are easy to fix.

### 8. Delete the implicit copy/move ctors explicitly
[include/eq.hpp:122-271](include/eq.hpp#L122-L271)

The compiler-generated copy ctor shallow-copies `mpfr_t` — both objects
end up aliasing the same limb pointers, so the second destructor
double-frees. CONTEXT documents this. Make it a hard compile error:

```cpp
Psim_f(const Psim_f&) = delete;
Psim_f(Psim_f&&) = delete;
Psim_f& operator=(const Psim_f&) = delete;  // existing op= takes non-const
```

You'd then have to fix `Psim_f(){...}` initializations like
`Psim_f x = some_psim_f;` — but those already produce silent bugs, so
better to find them.

### 9. Add the missing operators
- `double op Psim_f` (so users don't have to reorder operands)
- unary `operator-` (instead of `*(-1.0)`)
- `Psim_f sqrt() const`, `Psim_f abs() const` (or free functions) so the
  Newton code doesn't have to drop into `mpfr_sqrt(x.val, ...)` directly.

### 10. Make the comparison operators `const`
[include/eq.hpp:258-269](include/eq.hpp#L258-L269)

`cmp(Psim_f &o)` and `operator=(Psim_f &o)` take non-const refs, which
is why `const Psim_f` is unusable. Add `const` overloads where it
doesn't hurt the underlying mpfr API (cmp is genuinely read-only;
operator= is genuinely a write to LHS but should accept const RHS).

### 11. Document the temp-pool wrap-around
[include/eq.hpp:124-139](include/eq.hpp#L124-L139)

The 100-slot thread-local pool is a hard limit on intermediate-result
lifetime. Add a comment, or assert when `idx` wraps that no live ref
predates the wrap. Today nothing trips it, but a refactor that builds
larger expressions could.

---

## Performance — likely wins

### 12. Lower FLOAT_PREC
[include/eq.hpp:11](include/eq.hpp#L11)

1024 bits ≈ 308 decimal digits. The output is printed at `%.9Re` (9
significant digits), and the underlying physics has at most ~10 digits
of meaningful precision (DG values are double, exp() of double is
~15 digits). 256 bits gets you ~77 decimals — still wildly more than
needed, and 4× faster mpfr arithmetic on most operations. **Test:** drop
to 256, rerun `test_eq`, diff `regression_new.csv` vs `regression.csv` —
if the columns still match at 9 digits, ship the change.

### 13. Cache `dhds_to_eq_const` results
[src/eq.cpp:646, 716-728](src/eq.cpp)

`calc_strand_bindings` is called 25× per address per cycle and computes
`dhds_to_eq_const` for the same `(address, primer-end, primer-end)` pair
many times. Hoist the per-cycle dhds→eq_const conversions out of the
end5/end3 loops and into a small precomputed table indexed by
`[primer_f_or_r][addr_seq_kind]` once per cycle.

### 14. Avoid expression-temp explosion in `solve_eq`
[src/eq.cpp:257-258, 269-271](src/eq.cpp)

Each `+`, `*`, `/` allocates a new pool slot. `F1 = k[K_RR]*r*r*2.0 + (B1
+ k[K_FR]*f)*r + k[K_RH]*f - c0[R]` is ~8 mpfr_mul/add per evaluation
and ~6 temp slots. Rewriting with `mpfr_fma` (fused multiply-add)
halves the operations and eliminates intermediate rounding:

```cpp
// F1 = 2*k_RR*r*r + B1*r + k_FR*f*r + k_RH*f - c0[R]
mpfr_mul(t.val, k[K_RR].val, r.val, MPFR_RNDN);     // t = k_RR*r
mpfr_fma(t.val, t.val, r.val, c0[R].val, MPFR_RNDN); // not quite — work it out
// or: keep wrapper for readability and add a Psim_f::fma() method
```

Likely 1.3-2× speedup on the inner Newton loop.

### 15. Warm-start `solve_eq` from previous cycle
[src/eq.cpp:238-240](src/eq.cpp#L238-L240)

```cpp
c[F] = c0[F] / 2.0;   // discards the previous solution
c[R] = c0[R] / 2.0;
calc_cx();
```

Across PCR cycles `c[F]` and `c[R]` change slowly (they're 99%-bound to
target strands). Skipping the reset and reusing the prior solution
should cut outer iterations to 1-2 in steady state. Caveat: in the
first cycle there's no prior, and `c0[F]` itself changes between
cycles, so guard with a "first call" flag.

### 16. Cap the outer loop at something realistic
[src/eq.cpp:226](src/eq.cpp#L226)

`max_iter_outer = FLOAT_PREC * 2 = 2048`. In practice it should converge
in ≤5. Cap at 30 with an assert/warn on hitting the cap so divergence
becomes visible instead of silently chewing CPU.

---

## Performance — speculative

### 17. Multi-thread inside `sim_pcr`
The current threading parallelizes across addresses
([eval_addresses_thread, src/eq.cpp:997](src/eq.cpp#L997)) but a single
`sim_pcr` call is fully serial. The 5×5 loops over end5/end3 in
`update_strand_concs` and `calc_strand_bindings` are independent across
addresses (`i`) — could be parallelized for large address sets. Only
worth doing if a single-address simulation is the bottleneck.

### 18. Profile before tuning further
Most of the suggestions above are at-most-2× wins. Run `perf record
./test_eq` on a representative workload first; the hotspot may be
elsewhere (mpfr arithmetic in `calc_strand_bindings`, or `thal` calls
during the per-cycle dhds setup).

---

## Readability

### 19. Function-length budget
- `Primeanneal::sim_pcr` is ~240 lines
  ([src/eq.cpp:755-994](src/eq.cpp#L755-L994))
- `Primeanneal::eval_thread` is ~190 lines
  ([src/eq.cpp:418-611](src/eq.cpp#L418-L611))
- `Primeanneal::assign_addresses` is ~75 lines

Each one has clear sub-phases (per-address dhds precompute / per-cycle
update / output). Extract those into named helpers — even just splitting
the inner `for(temp_c = 40 ... 80)` body into a free function makes
`eval_thread` skimmable.

### 20. Document output column layout
The CSVs written by `sim_pcr` ([src/eq.cpp:956](src/eq.cpp#L956)),
`eval_thread` ([src/eq.cpp:669](src/eq.cpp#L669)), and the `print_state`
helper have wide column lists with no header row. Either:
- emit a header on first write, or
- add a comment block above the `mpfr_fprintf` enumerating the columns.

### 21. Inconsistent style in arithmetic
Code mixes `eq.tmp[0].add(eq.c[FX], eq.c[RX])` (method form) with
`eq.tmp[0] = eq.c[FX] + eq.c[RX]` (operator form). The method form
avoids the temp-pool allocation but is harder to read. Pick one
convention per function — reserve method form for the truly hot inner
loops where allocations measure.

### 22. Move `Psim_f` to its own header
It's a 150-line class definition embedded in [include/eq.hpp](include/eq.hpp).
Splitting `psim_f.hpp` from `eq.hpp` lets you reuse the wrapper outside
this project and shortens the eq header for casual readers.

### 23. Comment what `tmp[0..3]` mean per call site
[src/eq.cpp:704-708](src/eq.cpp#L704-L708) has the only good per-tmp
documentation in the codebase. Other functions reuse `tmp[0..3]` with
no commentary; one-line comments at the start of each function
(`tmp[0] = nonspec total, tmp[1] = ...`) make the body legible.

---

## Quick wins (cheap, no risk)

1. Remove dead `solve_cf_cr`, `calc_cf`, `calc_cr`, `bounds`,
   `solutions`, and unused EQ members. Mechanical and locally testable
   (`test_eq` regression check).
2. Delete the commented-out enum block in `eq.hpp`.
3. Fix the `read_addresses(in_filename, dna_conc)` bool-cast bug.
4. Replace `read_primers_individual`'s `printf` debug spew
   ([src/eq.cpp:309, 322-324](src/eq.cpp#L309)) with conditional logging
   or remove.
5. Add a header row to `regression*.csv` so future diffs are
   self-documenting.
