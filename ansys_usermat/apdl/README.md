# ANSYS verification decks

The USERMAT **runs and converges in ANSYS MAPDL 2022 R2 (v222)** — interface
arguments, `keycut`/`cutFactor` substepping and the `dsdePl` Jacobian were
validated on a `SOLID185` uniaxial-tension benchmark with `NLGEOM,ON`.

That benchmark ran with `alpha = 0`, so `Fg = I` and **the growth branch was
never entered**. The mechanical path is confirmed; the growth path is not. This
directory closes that gap.

> **Status, 2026-08-19 (IKMHIWI03): all four rows PASS, plus the `KEYOPT`
> sweep.** Ran inside a real custom-built ANSYS 2022 R2, matching the closed
> form to displayed precision for elastic and viscous, α=0.05 and α=0.20, and
> across `KEYOPT(1,2)` ∈ {0,1,3} (B-bar / enhanced / simplified enhanced
> strain — no volumetric locking detected). Evidence:
> [`out.txt`](out.txt), [`growth_result.txt`](growth_result.txt); full table
> and build/run procedure in [`RUNBOOK.md`](RUNBOOK.md).
>
> Getting the viscous rows to match required fixing a real deck bug first:
> `TBDATA` only accepts 6 values per call, and the deck's original
> `TBDATA,1,<9 values>` for `Fv` silently dropped values 7–9, leaving `Fv`
> singular rather than identity. Harmless for `η=0` (elastic doesn't use
> `Fv`), but for `η>0` it made the material silently return the elastic
> answer — a genuinely deceptive failure signature (bit-identical to the
> elastic row, not just "off"). Fixed in the committed deck by splitting
> across two `TBDATA` calls; see `RUNBOOK.md` Step 3 for the full story. This
> is a general APDL gotcha, not specific to this material — any state/property
> table needing more than 6 values needs multiple `TBDATA` calls.
>
> **Added 2026-08-20, UNRESOLVED: `t_growth_cylinder_shell.dat`, a two-layer
> curved-shell smoke test** (bonded substrate + growth layer, tooth-surface
> proxy geometry — a step toward realism beyond the unit cube; see
> `slides_1005.tex`'s caveat frame). Two real input bugs found and fixed on
> the way (VGLUE drops MAT attributes across its own renumbering; a
> per-volume element size coarser than the thin layer's own thickness
> produces a degenerate mesh) — but the solve still hits "element highly
> distorted" regardless of α or substep count. Not a closed-form case to
> begin with (no simple analytic answer for two bonded curved layers), so
> this is left as an open, documented next step rather than a pass/fail
> result. Full story and next ideas in the deck's own header comment.
>
> **Added 2026-08-19: a second, complementary closed-form check.**
> `t_growth_free.dat` removes only the 6 rigid-body modes (minimal 3-2-1
> constraint) instead of fixing every node, so the element is free to grow
> under zero traction. Closed form: stress ≡ 0 for any α, η. **PASS** —
> ANSYS returned stress ~1.9e−10 (shear ~1e−14), nine orders of magnitude
> below the constrained case's ~1e−4 scale, i.e. zero to solver tolerance.
> This catches sign/transpose errors in `Fg`/`Fg⁻¹` that full kinematic
> constraint can mask. Evidence: [`out_free.txt`](out_free.txt),
> [`growth_free_result.txt`](growth_free_result.txt).

## Algorithm flow

[`growth_verify_flow.tex`](growth_verify_flow.tex) /
[`growth_verify_flow_standalone.tex`](growth_verify_flow_standalone.tex) —
TikZ diagram of the check below (closed-form prediction branch vs. the actual
ANSYS solve branch, converging on a pass/fail comparison with the diagnostic
table). Same style/build convention as [`umat_flow/`](../../umat_flow/README.md):

```bash
cd ansys_usermat/apdl
pdflatex growth_verify_flow_standalone.tex
```

```latex
\begin{figure}[t]\centering
  \resizebox{0.7\linewidth}{!}{\input{ansys_usermat/apdl/growth_verify_flow.tex}}
  \caption{Closed-form verification of the ANSYS growth branch.}
\end{figure}
```

## The check

`t_growth_constrained.dat` fixes **every** node of a single hex. The deformation
gradient is then `F = I` exactly at every integration point, so with
`Fg = (1+α)I` the elastic part is `Fe = Fg⁻¹` and the Cauchy stress follows in
closed form. **No FE solve is needed to predict the answer** — which is what
makes this a verification rather than a comparison.

Expected: pure hydrostatic **compression** (the element wants to grow and
cannot), all three shear components exactly zero.

| Case | α | η | `Je` | `SX = SY = SZ` |
|---|---|---|---|---|
| elastic | 0.05 | 0 | 0.863837599 | −1.019275856e−04 |
| viscous | 0.05 | 8e−3 | 0.903444141 | −6.963159875e−05 |
| elastic | 0.20 | 0 | 0.578703704 | −4.726465185e−04 |
| viscous | 0.20 | 8e−3 | 0.780850893 | −1.795039360e−04 |

Properties `C10 = 2e−4`, `C01 = 0`, `D1 = 5e3`, `mtype = 0`; the viscous rows use
`dt = 5.0`, which the deck's `TIME` must match. Full precision:
[`reference_values.json`](reference_values.json).

### Complementary check: free (traction-free) growth

`t_growth_free.dat` is the opposite extreme: only a minimal 3-2-1 rigid-body
constraint (6 DOF, by node *location* so it's mesher-numbering-agnostic), not
a full fix. With zero traction, equilibrium has no elastic stretch at all —
the element just grows into the shape `Fg` already wants, so `F = Fg`,
`Fe = I`, and stress is **exactly zero**, independent of `α` and `η` (no
stress ⇒ no driving force for viscous flow either, so `Fv` stays at `I`).

| Case | α | η | `Je` | stress |
|---|---|---|---|---|
| free | 0.05 | 0 | 1.0 | 0 (ANSYS: ~1.9e−10, shear ~1e−14) |

Why bother, given `t_growth_constrained.dat` already passes: full kinematic
constraint can *mask* a sign or transpose error in how `Fg`/`Fg⁻¹` enters
`Fe = F·Fg⁻¹` — cancellation under `F = I` doesn't prove the sign is right,
only that the two errors (if any) cancel. Letting the element actually move
removes that cover.

## Reading a failure

The point of a closed-form case is that each failure mode is diagnostic:

| Symptom | Almost certainly |
|---|---|
| **Nonzero shear** (`SXY`/`SYZ`/`SXZ` ≠ 0) | the `VI/VJ` Voigt map is mis-wired — the Abaqus↔ANSYS 5↔6 shear swap |
| **Tensile** hydrostatic stress | `Fg` applied inverted (`Fg` where `Fg⁻¹` belongs) |
| `Je ≈ 1`, stress ≈ 0 | `alpha` never reached the material — check `TB,STATE` initialisation and `nStatev ≥ 10` |
| Elastic case matches, viscous does not | `dTime` mismatch — the deck's `TIME` must equal the `dt` in the reference. **Or**, if the viscous answer lands *exactly* on the elastic one (not just close): a `TBDATA` call set more than 6 values — APDL silently drops the tail, so `Fv`/whatever state came after value 6 never got set. Split into multiple `TBDATA` calls. |
| Right magnitude, wrong sign throughout | stress-sign convention on the `dsdePl`/`stress` return |

## Running

**Step-by-step for the ANSYS machine (Windows, v222): [`RUNBOOK.md`](RUNBOOK.md).**
It covers licence check-out, the Intel Fortran toolchain, building the custom
executable, and what each failure mode means.

```bat
REM Windows, from the directory holding the custom ANSYS.exe
"%AWP_ROOT222%\ANSYS\bin\winx64\ANSYS222.exe" -b -custom .\ANSYS.exe ^
    -i t_growth_constrained.dat -o out.txt
```

```bash
# Linux
ansys222 -b -custom ./ansys.e -i t_growth_constrained.dat -o out.txt
```

`-custom` must point at the executable built with `usermat_biofilm.f`. Without
it ANSYS runs its own stock material and the check silently means nothing.

Results land in `growth_result.txt`. Change `α` at the `TBDATA,10` line and keep
`TIME` equal to the `dt` in `reference_values.json` to pick another row.

## Regenerating the reference

```bash
python3 make_reference.py           # rewrite reference_values.json
python3 make_reference.py --check   # verify it still matches the core
```

The reference is computed from the **same** Python core that is verified
bit-identical to the Fortran (`tests/test_coupling_vs_fortran.py`, worst 6.8e-14
relative), so these numbers are not an independent hand-derivation — they are the
committed constitutive law evaluated at `F = I`. `tests/test_apdl_reference.py`
guards them in CI.

**Every entry in `reference_values.json` must be the `F = I` fully-constrained
scenario** — `test_apdl_reference.py` assumes it (e.g. asserts `stress[0] < 0`,
and keys dicts by `(alpha, eta)` alone, so a second case sharing an
`(alpha, eta)` pair silently overwrites the first in those dict comprehensions).
The complementary free/traction-free case above (stress ≡ 0 for any α, η) does
**not** belong in this file — it was added there once (2026-08-20, since
reverted) and broke 5 tests by colliding with `elastic_a005` on
`(alpha=0.05, eta=0.0)`. Keep that data in `t_growth_free.dat`'s own header
and in this README instead.

## Element formulation

`KEYOPT(1,2)` selects the `SOLID185` formulation. The deck uses the default
B-bar. Because constrained growth is a **volumetric** load on a near-incompressible
law, this is exactly the case where formulation choice bites.

**Swept 2026-08-19 on IKMHIWI03: all three formulations (`0` B-bar, `2`
enhanced strain, `3` simplified enhanced strain) recover the closed-form
stress exactly** — no volumetric locking detected for this load case on
`SOLID185`. See `RUNBOOK.md` Step 4.
