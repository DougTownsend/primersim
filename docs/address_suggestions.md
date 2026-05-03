# Code review — address evaluation path

This file collects suggestions specific to
[src/address_eval.cpp](../src/address_eval.cpp): the older
primer-design / address-assignment / per-temperature evaluation code
(`eval_thread`, `evaluate_addresses`, `assign_addresses`,
`assign_addresses_nosort`, `read_primers_individual`,
`shuffle_addresses`).

Items that affect `sim_pcr` / `solve_eq` / shared infrastructure live
in [docs/suggestions.md](suggestions.md) instead.

> **Status note.** The active workload (`test_eq.cpp`) doesn't go
> through this code today — it calls `sim_pcr` directly. These
> functions still build and link against the new `Real`-typedef
> tree, but they aren't exercised by the regression run.

---

## Cleanups

### 1. Function-length budget
- `Primeanneal::eval_thread` — ~170 lines [src/address_eval.cpp:154-323](../src/address_eval.cpp#L154-L323)
- `Primeanneal::assign_addresses` — ~75 lines [src/address_eval.cpp:68-141](../src/address_eval.cpp#L68-L141)

Each has clear sub-phases (per-address dhds precompute, per-cycle /
per-temperature update, output). Extract those into named helpers —
even just splitting the inner `for(temp_c = 40 ... 80)` body of
`eval_thread` into a free function makes the outer skimmable.

### 2. Output CSV column documentation
The CSV emitted by `eval_thread`
([src/address_eval.cpp:316](../src/address_eval.cpp#L316)) has 10
columns and no header row:

```cpp
fprintf(outfile, "%s,%s,%lf,%u,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le\n", ...)
```

Either:
- emit a header on first write, or
- add a comment block above the `fprintf` enumerating the columns.

Same idea for the 8-column file written by `assign_addresses` /
`assign_addresses_nosort` ([src/address_eval.cpp:134-138](../src/address_eval.cpp#L134-L138)).

### 3. Strip `read_primers_individual`'s `printf` debug spew
[src/address_eval.cpp:41, 47, 51, 53](../src/address_eval.cpp)

`read_primers_individual` `printf`s "primer is not 20 nt long"
warnings every time a non-20-nt primer is seen, plus dumps each
primer that is 20 nt to stdout. The stride of these prints leaks
into anything downstream that captures stdout. Either remove or
guard behind a `verbose` flag.

`assign_addresses` similarly `printf`s a per-primer counter and
post-iteration "iter: num_swaps" lines
([src/address_eval.cpp:83, 130](../src/address_eval.cpp)).

---

## Open questions

### 4. Is this code still wanted?
After the migration to native FP, this file was kept but the
`test_eq` workload only calls `sim_pcr`. If `eval_thread` /
`evaluate_addresses` / `assign_addresses` aren't part of the
near-term plan, this file is mostly dead weight.

If they are wanted: extract helpers (#1) before further work, since
the long single-function bodies make it hard to reason about. Also
note that the `eval_thread` change in the changelog below took the
"degenerate but defined" branch — `nonspec_exp_amp` will now be 0
for every address, which makes the F30/R30 selection loop pick
`first` as `best_temp` always. If that's not the intended
semantics, this needs a real fix.

---

## Done since previous revision (changelog)

- The "uninitialized `eq.tmp[1]`" bug (was #1 here) is now resolved
  by the `EQ::tmp[4]` deletion — `tmp[1]` no longer exists. Inside
  `eval_thread` it's been replaced with a 0-initialized local
  (`other_pool`) for the min/max split. This takes the
  "degenerate but defined" branch from the original suggestion: the
  min collapses to 0, so `nonspec_exp_amp` accumulates 0 every
  iteration and `nonspec_lin_amp` gets the full bound mass. The
  function compiles, runs, and produces deterministic output, but
  the science is still degenerate (see #4 above).
