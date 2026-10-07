# What the viscous update actually is

Short version: the repository calls the `Fv` update **"backward Euler"** in
several places. The flow increment is evaluated at the **old** state, so the
scheme is explicit in the flow direction. The label matters because "backward
Euler" implies unconditional stability, and this update has a step limit —
which is exactly the constraint recorded in `biofilm_material_v01.f` and handed
to the partner group.

## The update

Both implementations of the growth law do, per Gauss point:

```
  Fe_trial = F·Fg⁻¹·Fv_n⁻¹          (trial elastic state, from Fv at t_n)
  τ        = τ_dev(Fe_trial)         (flow driver, evaluated at the OLD state)
  Fv_n+1   = ( I + Δt/(2·η·J_e)·τ )·Fv_n
  σ        = σ(F·Fg⁻¹·Fv_n+1⁻¹)      (stress re-evaluated at the NEW state)
```

The re-evaluation of the stress at `Fv_n+1` is presumably what earned it the
name: the *stress* is taken at the end of the step. But `Fv_n+1` itself is
obtained from a driver computed entirely at `t_n`, with no iteration and no
residual. That is a forward-Euler flow increment, not a backward-Euler solve.

A genuinely implicit update would require `τ` at `Fv_n+1`, i.e. solving

```
  Fv_n+1 = ( I + Δt/(2·η·J_e(Fv_n+1))·τ(Fv_n+1) )·Fv_n
```

which nothing here does.

## What follows from it

The step must resolve the relaxation time `τ_relax = η/(2·C10)`. Measured on
the ANSYS core at `η = 5`, `C10 ≈ 167` (`τ_relax ≈ 0.015 s`):

| Δt | Δt/τ_relax | σ₁₁ |
|---|---|---|
| 1e−4 | 0.007 | −513 |
| 5e−3 | 0.33 | −134 |
| 1e−2 | 0.67 | **+347** |
| 2e−2 | 1.34 | +3288 |

The stress crosses zero near `Δt/τ ≈ 0.5` and diverges past `Δt/τ ≈ 1` — the
signature of an explicit increment taken past the timescale it is resolving, and
not something an unconditionally stable scheme does. Pinned by
`tests/test_material_wrapper.py::test_the_viscous_step_must_resolve_the_relaxation_time`.

## Where the label appears

**Corrected in place** (the growth law on the thesis critical path):

- `umat_biofilm_visco.f` — the Abaqus UMAT
- `ansys_usermat/usermat_biofilm.f` — the ANSYS USERMAT, 0 ULP against it
- `ch5_flow/flow_python_material_hook.tex` — feeds a thesis figure

**Deliberately left alone:**

- `rigor_audit_growth_2026-06-26.md` — a dated audit record. Editing what a
  past audit said would be rewriting history; this note is the correction.
- Everything about the **φ/c diffusion solve** (`tooth_pde3d.py`,
  `fem_2d_extension.py`, `fem_report.tex`, `JAXFEM/`). Those genuinely are
  backward Euler — an implicit linear solve per step — and are unaffected.

**Unverified, and not on this thesis's critical path:**

- `umat_biofilm_visco_phase2.f` and `test_umat_2ch.py` implement a *different*
  two-channel Prony model. Its `test_large_dt_stability` asserts a coarse step
  (`Δt/τ₁ = 10`) stays finite and reaches the same long-time value as a fine
  one. That is consistent with the transient being wrong while the long-time
  limit is set by the equilibrium branch and therefore path-independent — but
  this has not been checked, so nothing here should be read as a claim about
  that model either way.

## For the write-up

Say what the scheme is rather than naming it: *a single-step update whose flow
increment is evaluated at the previous viscous state, with the stress
re-evaluated at the updated one*. That is accurate, and it makes the step
restriction follow rather than contradict.

## If it is ever replaced: the exponential map, specified

Attempted 2026-09-30 and **reverted the same day**. This section is what the
attempt bought: the specification, so a later pass does not rediscover the
four complications from scratch. Nothing in the repository changed.

### The replacement

Same flow driver, integrated as a tensor exponential rather than a first-order
increment:

```
  A        = Δt/(2·η·J_e)·τ_dev(Fe_trial)     (symmetric, deviatoric)
  Fv_n+1   = exp(A)·Fv_n                       instead of ( I + A )·Fv_n
```

`A` is symmetric, so `exp(A)` is symmetric positive definite by construction:
`det Fv` cannot cross zero however long the step, and the sign flip in the
table above cannot happen. It is the standard cure and it removes the step
restriction rather than guarding it.

### What the switch is and is not worth — measured, not assumed

The case for it is **accuracy and the removal of a guard**, not the rescue of a
broken step. Two measurements, 300 random symmetric deviatoric `A` each:

| ‖A‖_F | cases with `det(I+A) ≤ 0` | ‖(I+A) − exp A‖ / ‖exp A‖ (median) |
|---|---|---|
| 0.10 | 0/300 | 0.2 % |
| 0.25 | 0/300 | 1.3 % |
| 0.50 | 0/300 | 4.8 % |
| 1.00 | 0/300 | 17 % |
| 2.00 | 290/300 | 43 % |

`‖A‖ ≈ (Δt/τ_relax)·‖dev strain‖/J_e`, so `DTMAX_RATIO = 0.5` in
`biofilm_material_v01.f` keeps the solve in the top rows. **Inside the guarded
range `I + A` does not lose positive determinacy at all** — the loss begins
near `‖A‖ ≈ 2`, which the guard already refuses. An earlier version of this
analysis claimed the failure started at `‖A‖ ≈ 1`; that was wrong, and it was
the main argument for switching, so it is worth stating plainly: the guard is
doing its job, and what remains is a 1–5 % truncation error on the viscous
increment inside the range that is allowed.

That is why the change was reverted rather than pushed through: it is a real
improvement to a part of the code that is **currently exercised at `η = 0` in
every reported run**, it touches the most heavily verified code in the
repository (0 ULP ANSYS↔Abaqus, closed-form growth reference, tangent to
2.4e−8), and it changes no reported number. Two months before submission that
is the wrong trade. It is a good post-submission change, or a good one the
moment a run actually needs `η > 0` at a step the guard refuses.

### The five files that must move together

The core is mirrored four times over, and a partial change would break the
equivalences that make the mirroring worth having:

- `ansys_usermat/usermat_biofilm.f` — the ANSYS USERMAT core
- `umat_biofilm_visco.f` — the Abaqus UMAT, held 0 ULP against it
- `ansys_usermat/coupling/material_server.py` — the NumPy reference
  (`_sigma_and_fv`)
- `ansys_usermat/coupling/material_jax.py` — the JAX mirror
  (`_sigma_and_fv_jax`), which the exact AD tangent inherits from
- `phase2_patch_test.py`

A shared helper (the attempt used `ansys_usermat/coupling/_viscous_exp.py`)
keeps the two Python copies from drifting; the two Fortran copies have to be
kept in step by hand, as they already are.

### The four complications, and their fixes

Each of these cost a debugging cycle. They are not hypothetical.

1. **The degenerate path overflows.** `J_e` is clamped at `1e-15`, so on a
   collapsed configuration `A` reaches `‖A‖ ~ 1e15` and the exponential
   overflows. A norm cap keeps it finite, but it must scale `A` *before*
   exponentiating (`A ← A·(cap/‖A‖)`), not clamp the result — clamping
   afterwards is a different operator and gives different answers in the two
   languages.
2. **`ABA_PARAM.INC` is `IMPLICIT REAL*8(A-H,O-Z)`.** Anything named I–N is
   an *integer* unless declared. A norm called `NRM` silently truncates; the
   symptom was a 2.95e−8 discrepancy against a 1e−11 tolerance, which reads
   like a tolerance problem and is not. Declare every new scalar explicitly.
3. **JAX needs a static squaring count.** Scaling-and-squaring with a
   compile-time `_NSQ`: 11 works. 60 does not — `0.5^60 ≈ 8.7e-19`, so
   `I + B` rounds back to the identity and the function silently returns
   `Fv_n`. Six Taylor terms with a `1/32` scaling threshold was the working
   combination.
4. **The degenerate path still does not reconcile.** With all of the above
   fixed, the degenerate battery cases differed between the Python and Fortran
   implementations by factors of 1e31–1e56: at `‖A‖ ~ cap` the exponential
   amplifies the difference between two orderings of the same arithmetic
   without bound. **This is the open design question, and it should be settled
   before the code is written again:** the right answer is probably to skip the
   viscous update entirely when `J_e` is at the clamp (the increment is being
   cut anyway — `detFe ≤ 1e-12` already sets `sKeyCut = 1` immediately after)
   rather than to integrate a meaningless `A` and then try to make two
   languages agree about the result.

### What else the change owes

- **`ansys_usermat/apdl/reference_values.json`** holds two viscous cases,
  `viscous_a005` and `viscous_a020` (`η = 0.008`, `Δt = 5.0`, `Δt/τ ≈ 0.25`).
  Both would move. Regenerating them is a reviewed step, not a refresh: they
  are the ANSYS↔Abaqus anchor, and a 2026-08-20 incident where
  `reference_values.json` was edited on an unverified basis silently broke
  five tests for a while.
- **`tests/test_material_wrapper.py::test_the_viscous_step_must_resolve_the_relaxation_time`**
  pins the sign flip that the exponential map removes. Its own docstring says
  what to do: *"if the integrator was changed on purpose, retire this test"*.
  Retiring it is correct, and it should be replaced by the assertion the new
  scheme can make and the old one cannot — that `det Fv > 0` and the stress
  sign hold at `Δt/τ = 1, 5, 50`.
- **`DTMAX_RATIO = 0.5` and the `sKeyCut` guard** in
  `biofilm_material_v01.f` become unnecessary for the viscous term. Whether to
  remove them is a separate call: the guard is cheap, and the partner group was
  handed the constraint explicitly.

## The scheme measured against closed form

Added 2026-09-30, after the specification above. The specification argued about
the change without ever measuring the scheme it proposed; this is that
measurement. Code: `ansys_usermat/viscous_integrator_verification.py`, pinned by
`tests/test_viscous_integrator_verification.py`. It runs in about a second, with
no ANSYS and without touching the verified UMAT.

The target is the ODE Soleimani 2019 Eq. 32 is derived from — the 1D Maxwell
element of his Eq. 31, `Q' + Q/τ = S'` — so the scalar case here is not a
convenience, it is the equation. `τ = 0.01 s`, his Table 2 value.

**Relaxation** (`S' = 0`, exact answer `Q0 exp(−t/τ)`):

| Δt/τ | Eq. 32, relative error | forward Euler, relative error |
|---|---|---|
| 0.1 | 3.0e−15 | **0.235** |
| 0.5 | 0 | 0.855 |
| 1.0 | 3.9e−16 | 1.00 |
| 2.0 | 1.9e−16 | 53.6 |
| 10 | 0 | 2.0e+5 |
| 100 | 0 | 2.7e+45 |

Eq. 32 reduces here to `Q_{n+1} = exp(−Δt/τ) Q_n`, whose n-th iterate *is* the
exact solution sampled at `t_n` — so it is exact at every step size, and the
paper's "cannot go unstable however coarse the step" is not merely stability
but exactness on this problem. Forward Euler gives `(1 − Δt/τ)^n`: already 23 %
wrong at a step one tenth of `τ`, no information left at `Δt = τ`, divergent
past 2.

**Ramp** (`S' = r`, exact answer `r τ (1 − exp(−t/τ))`) separates them on
accuracy instead of stability:

| | observed order |
|---|---|
| Eq. 32 | **1.99, 2.00, 2.00, 2.00** |
| forward Euler | 0.69, 0.87, 0.94, 0.97 |

**One honest qualification, pinned by a test so it cannot quietly disappear.**
At the coarsest ramp step measured (`Δt/τ = 0.5`) forward Euler is the *more*
accurate of the two — 5.8e−5 against 1.0e−4. The exponential update's advantage
is order and unconditional stability, not a smaller error at every step size,
and the case for the change should be made on those terms.

This does not alter the decision recorded above: the repository's own runs use
`η = 0`, `DTMAX_RATIO = 0.5` already refuses the steps where the explicit
update misbehaves, and the change moves five mirrored files in the most
verified part of the code. It does mean the specification now rests on a
measurement of the proposed scheme rather than on an argument about it.
