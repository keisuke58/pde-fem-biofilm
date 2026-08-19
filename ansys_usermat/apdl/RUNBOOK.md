# Runbook — verifying the growth path on IKMHIWI03 (ANSYS 2022 R2, Windows)

For whoever is sitting at the ANSYS machine. The repo is already cloned.

**What you are checking and why.** The USERMAT already runs and converges in
ANSYS 2022 R2 — that was established with a `SOLID185` uniaxial-tension
benchmark, which validated the interface arguments, the `keycut`/`cutFactor`
substepping and the `dsdePl` Jacobian. But that benchmark ran with `α = 0`, so
`Fg = I` and **the growth branch was never executed**. Growth is the entire
point of this material. This runbook exercises it, against an answer known in
closed form.

Budget roughly an hour, most of it in step 0 if the build environment needs
setting up.

---

## Step 0 — Pre-flight (do these before anything else)

Run each and note the result. If one fails, stop there — the later steps depend
on it.

**0a. License reachable.** Licences come from the RRZN floating server
(`1055@ansys-lic.rrzn.uni-hannover.de`), so you need to be on the university
network or VPN.

```bat
"C:\Program Files\ANSYS Inc\Shared Files\Licensing\winx64\ansysli_util.exe" -checkout ansys
```

**0b. Fortran toolchain present.** This is the step most likely to be missing —
the environment report lists no Intel Fortran or Visual Studio. ANSYS 2022 R2 on
Windows needs **Intel Fortran (oneAPI) plus a matching Visual Studio** to build a
custom executable; the version pairing is specified in the ANSYS 2022 R2
Installation Guide's platform-support table, and a mismatched pair fails at link
time with unhelpful errors.

Open the **"Intel oneAPI command prompt for Intel 64 for Visual Studio"** from
the Start Menu and run:

```bat
ifort --version
```

> **If you already built the custom executable for the uniaxial benchmark, you
> have a working toolchain — reuse it, and skip to step 1.** Please also write
> down what you did, so step 1 below can be corrected to match reality; the exact
> procedure is not currently recorded anywhere in this repo.

**0c. Working directory.** Do **not** build inside `C:\Program Files\...` — it
needs administrator rights and pollutes the install. Copy the UPF template
directory somewhere writable:

```bat
xcopy /E /I "C:\Program Files\ANSYS Inc\v222\ansys\custom\user\winx64" C:\work\biofilm_upf
```

---

## Step 1 — Build the custom executable

Copy the USERMAT source in, replacing the stock stub:

```bat
copy <repo>\ansys_usermat\usermat_biofilm.f C:\work\biofilm_upf\
cd /d C:\work\biofilm_upf
```

Then run the ANSYS customisation script from the **Intel oneAPI command prompt**
(not a plain `cmd`, or `ifort` will not be found):

```bat
ANSCUST.bat
```

This should produce `ANSYS.exe` in `C:\work\biofilm_upf`.

> ⚠️ **This build command is written from the documented UPF procedure, not from
> a run on this machine** — nobody has recorded the Windows build here yet. If
> the script has a different name in v222, or needs arguments, check
> `C:\Program Files\ANSYS Inc\v222\ansys\custom\user\winx64\` for what is
> actually there and **tell me what worked** so this file can be corrected.

**Sanity check before going further:** confirm the build picked up *your* source
and not the stock stub. The simplest evidence is that step 2 produces nonzero
stress at all — the stock `usermat` stub returns zero.

---

## Step 2 — Run the growth deck

```bat
cd /d C:\work\biofilm_upf
copy <repo>\ansys_usermat\apdl\t_growth_constrained.dat .
"%AWP_ROOT222%\ANSYS\bin\winx64\ANSYS222.exe" -b -custom .\ANSYS.exe ^
    -i t_growth_constrained.dat -o out.txt
```

The `-custom` argument is essential. Without it ANSYS runs its **own** material
and the whole exercise silently means nothing.

Results land in `growth_result.txt` (`PRESOL,S,COMP` and `PRESOL,SVAR`).

### If the solve refuses to run

Every node is fixed and no external load is applied, so ANSYS may object that
there is nothing to solve. That is a deck problem, not a material problem. Two
fallbacks, in order of preference:

1. Free a single corner node's `UX` only. `F` stays `I` at the integration
   points to within solver tolerance, and the reference values still apply.
2. Apply a tiny nonzero prescribed displacement (`1e-9`) on one node instead of
   `0`. Same reasoning.

Tell me if you need either — the reference values would then want regenerating
for the exact `F` rather than `F = I`.

---

## Step 3 — Compare against the closed-form answer

The element is fully constrained, so `F = I` exactly, `Fe = Fg⁻¹`, and the stress
follows from the constitutive law with no FE solve needed to predict it.

With the deck as committed (`α = 0.05`, `η = 0`, `TIME = 5.0`):

| Quantity | Expected |
|---|---|
| `SX = SY = SZ` | **−1.019275856e−04** |
| `SXY = SYZ = SXZ` | **exactly 0** |
| `SVAR(10)` (α) | 0.05 |
| `SVAR(1..9)` (Fv) | identity, since `η = 0` |

Full precision and the other three cases (α = 0.20, and the viscous rows with
`η = 8e−3`) are in [`reference_values.json`](reference_values.json). To run
another row, edit `TBDATA,10,<α>` and keep `TIME` equal to the `dt` in that file.

### What each failure mode means

This is why a closed-form case is worth the setup — the failures are diagnostic,
not mysterious:

| What you see | Almost certainly |
|---|---|
| **Nonzero shear** | the `VI/VJ` Voigt map is mis-wired — the Abaqus↔ANSYS 5↔6 shear swap |
| **Tensile** (positive) hydrostatic stress | `Fg` applied inverted — `Fg` used where `Fg⁻¹` belongs |
| `Je ≈ 1`, stress ≈ 0 | `α` never reached the material — check `TB,STATE` and that `nStatev ≥ 10` |
| Stress exactly 0 everywhere | the build did not pick up `usermat_biofilm.f` — stock stub still linked |
| Elastic row matches, viscous does not | `TIME` ≠ the reference `dt` |

---

## Step 4 — Element formulation sweep (worth the extra 10 minutes)

Constrained growth is a **volumetric** load on a near-incompressible law, which
is exactly where element formulation bites. Re-run step 2 three times, changing
only the `KEYOPT` line in the deck:

```apdl
KEYOPT,1,2,0     ! B-bar (what the deck ships with)
KEYOPT,1,2,2     ! enhanced strain
KEYOPT,1,2,3     ! simplified enhanced strain
```

Record the hydrostatic stress for each. A formulation that misses the
closed-form value here will also distort results on the real tooth-shell mesh —
this is the cheapest possible way to find that out, and it bears directly on a
volumetric-locking question that is currently open on the Abaqus side too
(linear tets, `C3D4`, no hybrid formulation).

---

## Step 5 — Report back

Paste or send:

1. `out.txt` and `growth_result.txt` (or just the `SX/SY/SZ/SXY/SYZ/SXZ` block
   and the `SVAR` block).
2. The three hydrostatic stresses from step 4.
3. **What the build actually took** — the exact script/command, the `ifort` and
   Visual Studio versions, anything that had to be worked around. That is the
   part not currently captured anywhere, and the next person (including you in
   six months) will need it.

Items 1–2 go into the repo as verification evidence; item 3 becomes the
corrected version of step 1.
