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

## Element formulation

`KEYOPT(1,2)` selects the `SOLID185` formulation. The deck uses the default
B-bar. Because constrained growth is a **volumetric** load on a near-incompressible
law, this is exactly the case where formulation choice bites.

**Swept 2026-08-19 on IKMHIWI03: all three formulations (`0` B-bar, `2`
enhanced strain, `3` simplified enhanced strain) recover the closed-form
stress exactly** — no volumetric locking detected for this load case on
`SOLID185`. See `RUNBOOK.md` Step 4.
