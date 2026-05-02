# How `EQ::solve_eq` works

A walkthrough of the equilibrium solver in
[src/eq.cpp:solve_eq](src/eq.cpp), aimed at someone who hasn't seen
Newton's method or Jacobians before.

---

## 1. The problem

In each PCR cycle we want to know how much *free* primer is floating
around in solution after the reactions reach equilibrium. We track
three concentrations:

- `c[F]` — free forward primer
- `c[R]` — free reverse primer
- `c[X]` — free nonspecific binding sites

The total amount of each species is fixed (`c0[F]`, `c0[R]`, `c0[X]`).
Some of each total is "free", and the rest is bound to other species.
Mass balance says

```
total = free + sum of all bound forms
```

Writing that out for the forward primer:

```
c0[F] = c[F]                       ← free
      + k_FH·c[F]                  ← in a hairpin with itself
      + 2·k_FF·c[F]²               ← in an F-F dimer (factor 2 because it
                                     uses two F molecules)
      + k_FR·c[F]·c[R]             ← in an F-R dimer
      + k_FX·c[F]·c[X]             ← bound to a nonspecific site
```

Each `k_*` is an equilibrium constant: at equilibrium, the bound
concentration equals `k * (concentrations of the binders)`. Bigger `k`
means tighter binding.

Doing the same for reverse primer:

```
c0[R] = c[R]
      + k_RH·c[R]                  ← wait, this is k_RH·c[F] in the code,
                                     because the "RH" hairpin actually uses
                                     the reverse primer in a way that scales
                                     with c[F] in this simplified model.
                                     See the comment at line 98 of the
                                     pre-cleanup eq.cpp for the original.
      + 2·k_RR·c[R]²
      + k_FR·c[F]·c[R]
      + k_RX·c[R]·c[X]
```

And for the nonspecific pool:

```
c0[X] = c[X] + k_FX·c[F]·c[X] + k_RX·c[R]·c[X]
```

That last one rearranges cleanly:

```
c[X] = c0[X] / (1 + k_FX·c[F] + k_RX·c[R])
```

So **once we know `c[F]` and `c[R]`, we get `c[X]` for free.** That's
what `EQ::calc_cx()` does.

So the real puzzle is finding `c[F]` and `c[R]` — two unknowns linked
by two equations.

## 2. Why this is hard

The two equations look linear at first glance, but they aren't —
`c[F]²`, `c[R]²`, and `c[F]·c[R]` are quadratic and bilinear terms.
And the equations are *coupled*: you can't solve for `c[F]` without
knowing `c[R]`, and vice versa.

There's no general formula for this kind of system. We have to
**iterate**: guess a solution, see how wrong it is, improve the guess,
repeat until the wrongness shrinks below precision noise.

The classic algorithm for that is **Newton's method**.

## 3. Newton's method in 1D (warm-up)

Imagine just one equation, one unknown: `g(x) = 0`. We want the `x`
that makes `g` zero.

Newton's idea: at any point `x`, you can approximate `g` by its tangent
line. The tangent line crosses zero at a (usually better) guess. So:

1. Pick a starting `x₀`.
2. Compute the value `g(x₀)` and the slope `g'(x₀)`.
3. The tangent line through `(x₀, g(x₀))` with slope `g'(x₀)` hits
   zero at:

   ```
   x₁ = x₀ - g(x₀) / g'(x₀)
   ```

4. Repeat with `x₁`, then `x₂`, ... Each step roughly **doubles** the
   number of correct digits, so it converges very fast once you're
   close.

The formula `x_new = x_old - g/g'` is the heart of Newton's method.
Everything below is just generalizing that to two unknowns.

## 4. Newton's method in 2D — the Jacobian

Now we have **two** functions of **two** variables. Let `f` be our
guess for `c[F]` and `r` be our guess for `c[R]`. Define:

```
F1(f, r) = (the c0[R] equation, rearranged so the right side is zero)
F2(f, r) = (the c0[F] equation, rearranged so the right side is zero)
```

We want the `(f, r)` that makes `F1 = 0` AND `F2 = 0` at the same
time.

In 1D, we needed one slope `g'(x)`. In 2D, each function has *two*
slopes — one in each direction. So we get four slopes total:

```
∂F1/∂f   ∂F1/∂r
∂F2/∂f   ∂F2/∂r
```

`∂F1/∂f` (read "partial F1 partial f") means: "if I bump `f` a little
and hold `r` fixed, how much does `F1` change?" It's just a slope, but
restricted to one direction.

Stacking those four slopes into a 2-by-2 grid gives the **Jacobian
matrix**:

```
J = [ ∂F1/∂f   ∂F1/∂r ]   = [ J11  J12 ]
    [ ∂F2/∂f   ∂F2/∂r ]     [ J21  J22 ]
```

That's it. A Jacobian is just the multi-variable version of a
derivative — a table of partial slopes.

In 1D, Newton's update was `x_new = x_old - g/g'`. In 2D it's
analogous: the "correction" `(df, dr)` we add to our guess satisfies

```
J · [df]   = - [F1]
    [dr]       [F2]
```

— a small linear system. Solve it, add the corrections to `(f, r)`,
repeat.

## 5. Solving the 2-by-2 by Cramer's rule

We need `[df; dr]` from the equation `J · [df; dr] = -[F1; F2]`. For
a 2-by-2 matrix there's a clean closed form:

First, the **determinant** — a single number that summarizes how
"invertible" the matrix is:

```
det = J11·J22 - J12·J21
```

If `det = 0`, the matrix doesn't have a unique inverse and Newton can't
take a step. (We don't currently guard against this; in our problem
the decoupled initial guess is exact when `det → 0`, so the residual
check exits before we'd hit it. But it's worth knowing.)

Cramer's rule then gives:

```
df = (J12·F2 - J22·F1) / det
dr = (J21·F1 - J11·F2) / det
```

Those two lines in the code:

```cpp
det = J11*J22 - J12*J21;
df  = (J12*F2 - J22*F1) / det;
dr  = (J21*F1 - J11*F2) / det;
```

are doing exactly that.

## 6. Computing `F1`, `F2`, and the four `J` entries

Pull out the parts of each equation that don't depend on `f` or `r` —
they're constant within one solve:

```
B1 = 1 + k_RX · c[X]                     (the coefficient of r in F1)
B2 = 1 + k_FH + k_FX · c[X]              (the coefficient of f in F2)
```

Then

```
F1 = 2·k_RR·r²  +  (B1 + k_FR·f) · r  +  k_RH·f  -  c0[R]
F2 = 2·k_FF·f²  +  (B2 + k_FR·r) · f             -  c0[F]
```

Take partial derivatives one variable at a time. For example,
`∂F1/∂f`: treat `r` as a constant and differentiate term-by-term:

- `2·k_RR·r²`     → 0   (no `f`)
- `(B1 + k_FR·f)·r` → `k_FR · r`
- `k_RH · f`      → `k_RH`
- `c0[R]`        → 0

Sum: `J11 = ∂F1/∂f = k_RH + k_FR·r`. Match it to the code:

```cpp
J11 = k[K_RH] + k[K_FR]*r;
```

Same procedure for the others:

```
J11 = ∂F1/∂f = k_RH + k_FR·r
J12 = ∂F1/∂r = 4·k_RR·r + B1 + k_FR·f
J21 = ∂F2/∂f = 4·k_FF·f + B2 + k_FR·r
J22 = ∂F2/∂r = k_FR·f
```

The `4·k_RR·r` comes from `d/dr (2·k_RR·r²) = 4·k_RR·r`. (The
factor-of-2 in the original equation, doubled by the chain rule.) Same
for `4·k_FF·f`.

## 7. The full inner-loop step

Putting it all together, one Newton iteration is:

```
1. Compute F1, F2 from the current (f, r).
2. Compute residual = |F1| + |F2|. If it didn't decrease this step,
   we've hit precision noise — stop iterating.
3. Compute J11, J12, J21, J22.
4. Compute det.
5. df = (J12·F2 - J22·F1) / det
   dr = (J21·F1 - J11·F2) / det
6. f_new = f + df
   r_new = r + dr
7. Clamp f_new and r_new into the physically valid box [0, c0[F]] ×
   [0, c0[R]]:
     - if step would push below 0 → halve the current value
     - if step would push above c0 → take the midpoint with c0
8. Repeat.
```

We give it up to 15 iterations. Newton-on-a-quadratic typically
converges in 4-6 from a decent guess.

## 8. The initial guess (line 251-254)

Newton converges *quadratically* once you're close, but it can
diverge if you start too far away. We help it by giving it a smart
starting point.

The trick: drop the cross-coupling terms (`k_FR` and `k_RH`) from each
equation. That decouples them — each becomes a standalone quadratic
in one variable. Quadratics have a closed-form solution.

For `f` with the coupling dropped:

```
2·k_FF·f² + B2·f - c0[F] = 0
```

Standard quadratic `ax² + bx + c = 0` with `a = 2·k_FF`, `b = B2`,
`c = -c0[F]`. The positive root is

```
f = (-B2 + √(B2² + 8·k_FF·c0[F])) / (4·k_FF)
```

— the `8` and `4` come from `4ac = 8·k_FF·c0[F]` and `2a = 4·k_FF`.
Same form for `r`. The `mpfr_sqrt` call in the code is computing
`√(...)` since the `Psim_f` wrapper doesn't have its own `sqrt`.

When `k_FR` and `k_RH` happen to be very small, this closed form is
already the answer — Newton's first iteration computes `F1` and `F2`,
sees the residual is at the precision floor, and exits immediately.

## 9. The convergence test (line 262-266)

After each Newton step, we check `|F1| + |F2|`. That's the **residual** —
how badly the equations are violated. We stop when the residual stops
*decreasing*:

```cpp
if (tmp1.cmp(prev_res) >= 0) break;
mpfr_set(prev_res.val, tmp1.val, MPFR_RNDN);
```

Why "stops decreasing" instead of "is below 1e-12"?

With mpfr's high precision (1024 bits ≈ 308 decimal digits), the
arithmetic itself has noise at the last bit. Once Newton converges to
that floor, further iterations bounce around without improving — the
residual stops shrinking and may even tick up. Detecting that
automatically gives the tightest answer the precision allows, no magic
constant required.

The downside: this stops *as soon as* the residual stops decreasing,
even if it's a brief plateau. In practice with quadratic convergence
that's fine — by the time we're plateauing, we're at the floor.

## 10. The outer loop (line 244-297)

Everything above solves the F/R system **for a fixed `c[X]`**. But
`c[X]` itself depends on `c[F]` and `c[R]`! Recall:

```
c[X] = c0[X] / (1 + k_FX·c[F] + k_RX·c[R])
```

The outer loop handles this:

```
1. Initial guess: c[F] = c0[F]/2, c[R] = c0[R]/2, then compute c[X].
2. Run the Newton inner loop with that c[X].
3. Recompute c[X] from the new c[F], c[R].
4. If c[F] and c[R] didn't change (bit-exact comparison), we're done.
5. Otherwise, save them to last_val and repeat.
```

For the workloads we run (where `c0[X]` is small relative to
`c0[F]/c0[R]`), this typically converges in 2-4 outer iterations.

The bound `max_iter_outer = FLOAT_PREC * 2 = 2048` is much higher than
needed — leftover from the old binary-search solver. Worth lowering
once we trust the Newton path.

## 11. After the inner-and-outer dance: `calc_bound_concs()`

Once `c[F]`, `c[R]`, `c[X]` are settled, every "bound" species can be
computed as a simple product. `calc_bound_concs` populates the rest of
the `c[]` array (the `FH`, `RH`, `FF`, `RR`, `FR`, `FX`, `RX` slots)
so downstream code can read them.

## 12. Where it could go wrong

In rough order of likelihood:

1. **`det = 0`** — not currently guarded. If both off-diagonal entries
   collapse it would divide by zero. In practice the decoupled initial
   guess is already very close, so the residual check exits before we
   ever evaluate `det`. But a stiff parameter set could trip it.
2. **Newton diverges**. If the initial guess is too far from the root,
   `df` or `dr` can be huge and shoot the iterate out of the valid
   box. The clamping `f_new = f*0.5` etc. handles the "shot below 0
   or above c0" case, but it's heuristic — there's no theoretical
   guarantee of progress in pathological cases.
3. **Outer loop does too few iterations**. The bit-exact `cmp` check
   may not trigger if the last digit of `c[F]` keeps flipping; the
   loop will then run all `FLOAT_PREC * 2` outer iterations. Slow, but
   produces a correct answer.
4. **`mpfr_set_inf(prev_res.val, 1)`** at the start of each outer
   iteration resets the residual history. That means a Newton solve
   that almost converged in outer iter 1 has to "rediscover" its
   precision floor in outer iter 2. Not a bug, but a tiny
   inefficiency.

## TL;DR

- Three unknowns (`c[F]`, `c[R]`, `c[X]`), three equations.
- `c[X]` rearranges in closed form, leaving 2 unknowns and 2
  equations.
- A 2-by-2 nonlinear system is solved by **Newton's method**: at each
  step, build the Jacobian (table of partial slopes), invert it via
  Cramer's rule, take a corrective step.
- Inner loop: Newton on `(f, r)` with `c[X]` held constant.
- Outer loop: refresh `c[X]` from the updated `(f, r)`, repeat.
- Stop when the residual stops shrinking (precision floor) or when
  `(f, r)` stop changing.
