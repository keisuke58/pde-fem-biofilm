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

> **2026-08-20: on IKMHIWI03 specifically, use `F:\biofilm_upf`, not
> `C:\work\biofilm_upf`.** `C:` on this machine is chronically near-full
> (likely VSS retaining deleted blocks — see the disk-space memory / draft
> email to Timo) and has hit 0 bytes free before. `F:\` is a separate local
> data drive with ~3.7 TB free, confirmed to run the built executable
> identically. The already-built environment was copied there; new builds
> on this machine should go straight to `F:\biofilm_upf` instead of
> following the `C:\work\biofilm_upf` paths below literally.
> `ansys_usermat/apdl/run_apdl.ps1`'s default already points at `F:`.

---

## Step 0 — Pre-flight (do these before anything else)

Run each and note the result. If one fails, stop there — the later steps depend
on it.

**0a. License reachable.** Licences come from the RRZN floating server
(`1055@ansys-lic.rrzn.uni-hannover.de`), so you need to be on the university
network or VPN.

```bat
"C:\Program Files\ANSYS Inc\v222\licensingclient\winx64\ansysli_util.exe" -checkout ansys
```

> The path above is corrected from an earlier version of this doc — on
> IKMHIWI03 there is no `Shared Files\Licensing\winx64\ansysli_util.exe`; the
> real one lives under `v222\licensingclient\winx64\`. Verified 2026-08-19:
> checkout succeeds and resolves `Ansys Mechanical Enterprise` — but locally
> (`server=55206@ikmhiwi03...`), not via the RRZN address configured in
> `ansyslmd.ini`. Worth investigating if a checkout ever fails unexpectedly.

**0b. Fortran toolchain present.** Verified 2026-08-19 on IKMHIWI03 — both are
present, contrary to what the original environment report said:
`ifort.exe` (Intel Fortran 2025.3.3) is at
`C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin\ifort.exe`, and
Visual Studio **18** (2026 Developer) is at
`C:\Program Files\Microsoft Visual Studio\18\Community\`. VS18 is newer than
the VS2019-era pairing ANSYS 2022 R2's Installation Guide expects for this
compiler version — a link-time version mismatch is still possible and has not
been ruled out, only the toolchain's mere presence.

Plain `ifort --version` from a bare shell, and even the Start Menu's "Intel
oneAPI command prompt for Intel 64 for Visual Studio" via `setvars.bat`, do
**not** reliably work here: `setvars.bat` shells out to a bare `vswhere.exe`
that isn't itself on `PATH`
(`C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe`), so
VS init warns/fails, and the per-component `compiler`/`mpi`/`umf`
`env\vars.bat` calls then fail with "command not found" — `ifort` never lands
on `PATH`. The init sequence that actually works calls the two relevant env
scripts directly, skipping `setvars.bat` entirely:

```bat
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\env\vars.bat"
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

Then, with the env from step 0b active in the same shell (`vcvars64.bat` +
Intel `compiler\...\env\vars.bat`, not `setvars.bat`), run the ANSYS
customisation script — it is named `ANSCUST.BAT` (uppercase) in v222:

```bat
ANSCUST.BAT
```

This should produce `ANSYS.exe` in `C:\work\biofilm_upf`.

> ✅ **Resolved 2026-08-19 on IKMHIWI03 — link error fixed without installing
> classic ifort, no admin/UAC needed.** The first attempt hit:
> ```
> ifmodintr.lib(iso_c_binding.obj) : error LNK2001: unresolved external symbol for_deallocate_handle
> ifmodintr.lib(iso_c_binding.obj) : error LNK2001: unresolved external symbol for_alloc_allocatable_handle
> ANSYS.exe : fatal error LNK1120: 2 unresolved externals
> ```
> The instinct was "this needs the old classic ifort" (oneAPI 2024.2.1, last
> version to ship it) — but `winget install --id Intel.FortranCompiler
> --version 2024.2.1` needs elevation with no `--scope user` option, and the
> UAC consent prompt can't be approved from a non-interactive session. That
> path is a dead end without a human physically at the machine.
>
> **The actual fix needs neither.** `dumpbin /symbols` on oneAPI 2025.3's own
> `ifmodintr.lib` confirms `for_alloc_allocatable_handle`/
> `for_deallocate_handle` are referenced there as `UNDEF` (i.e. `ifmodintr.lib`
> expects something else to supply them) — but they're **not** missing from
> 2025.3 as a whole: `dumpbin /symbols` on `libifcoremt.lib` (same 2025.3
> install) shows both as real, defined (`SECTB`/`SECTD`) symbols. ANSYS's
> shipped `ansys.lrf` link-response-file simply never lists
> `libifcoremt.lib` among its ~150 `-defaultlib:` entries. Adding one line
> fixes it — edit the **copy** of `ansys.lrf` in your working folder (never
> the template under `Program Files`), immediately before the `ansysexe.res`
> line:
> ```
> -defaultlib:libifcoremt.lib
> ```
> Two more snags surfaced re-running just `link @ansys.lrf` directly (to
> avoid re-doing `ANSCUST.BAT`'s interactive prompt each time) — both are
> artifacts of that shortcut, not real problems, and don't come up if you run
> the full `ANSCUST.BAT`:
> - `LINK1181: cannot open ansysexe.res` — `ANSCUST.BAT` sets `LIB` to
>   include `%AWP_ROOT222%\ansys\Custom\Lib\winx64` (where `ansysexe.res` and
>   `WinAnsys.res` actually live) before linking; a bare `link @ansys.lrf`
>   run outside that context doesn't have it on `LIB`. Set
>   `LIB=%AWP_ROOT222%\ansys\Custom\Lib\winx64;%LIB%` first if you do this.
> - `LNK1149: output filename identical to input` — a stale `ANSYS.lib`/
>   `ANSYS.exp` left over from a prior failed link attempt gets swept up by
>   the response file's `*.lib` wildcard and collides with the new output of
>   the same name. Delete `ANSYS.lib`/`.exp`/`.exe`/`.map` from the working
>   folder before relinking (`ANSCUST.BAT` only auto-deletes `ANSYS.exe`,
>   not the `.lib`/`.exp`).
>
> With those three things in place, `link @ansys.lrf` (or the full
> `ANSCUST.BAT`, answering `Y`/`N` to Wind Turbine Aeroelastic as you
> prefer — it's unrelated) exits 0 and produces a 370 MB `ANSYS.exe`, only
> `LNK4286`/`LNK4199` warnings (duplicate-symbol-import and unused-delayload
> noise, harmless).

> ⚠️ **Confirmed 2026-08-19 that this script runs and is genuinely
> interactive** — it prompts `Do you want to link the Wind Turbine Aeroelastic
> library with Mechanical APDL? (Y or N):` and likely further Y/N questions
> after that. The prompts are read by a bundled `ASK.EXE` that reads the
> console directly, **not stdin** — piping answers into a redirected
> `cmd.exe` does not work (it just loops "Please answer Y or N" until it gives
> up). This step must be run interactively by a human at the machine; it
> cannot be scripted or driven by an agent. Answer `N` to Wind Turbine
> Aeroelastic unless you specifically need it — it's unrelated to this
> USERMAT. **Tell me what the remaining prompts were and what you answered**
> so this file can be corrected with the full prompt sequence.

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

> ✅ **Ran successfully 2026-08-19 on IKMHIWI03**, once past two more Step-1-
> adjacent snags:
> - `ANSYS222.exe -custom .\ANSYS.exe` first failed with
>   `EXIT STATUS: -1073741515 (0xc0000135)` — `STATUS_DLL_NOT_FOUND`.
>   `ANSCUST.BAT` has its own late prompt ("Do you want to copy the runtime
>   DLLs?") that a bare `link @ansys.lrf` shortcut skips; copy them manually:
>   `copy "%AWP_ROOT222%\ansys\Bin\winx64\*.dll" .` and
>   `copy "%AWP_ROOT222%\commonfiles\AAS\bin\winx64\*.dll" .`
> - `PRESOL,SVAR` (bare, as originally committed in this deck) failed with
>   `*** WARNING *** No components are specified for the SVAR item.` — it
>   needs an explicit component number per call
>   (`PRESOL,SVAR,1` … `PRESOL,SVAR,10`, one call per state variable — the
>   deck below is already fixed). Even with that fixed, the state variables
>   weren't in the results file at all (`The requested SVAR data is not
>   available`) until `OUTRES,SVAR,ALL` was added before `SOLVE` —
>   `OUTRES,ALL,ALL` does **not** imply it.
>
> With those two fixed (both now baked into the committed
> [`t_growth_constrained.dat`](t_growth_constrained.dat)), the run exits 0,
> `RUN COMPLETED`, and produces exactly the closed-form answer — see Step 3.
> Evidence committed: [`out.txt`](out.txt), [`growth_result.txt`](growth_result.txt).

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

> ✅ **PASS, 2026-08-19 on IKMHIWI03** — real ANSYS 2022 R2 output, custom
> `usermat_biofilm.f` build, `elastic_a005` case (α=0.05, η=0):
>
> | Quantity | Expected | Got |
> |---|---|---|
> | `SX=SY=SZ` | −1.019275856e−04 | **−0.10193E−003** |
> | `SXY=SYZ=SXZ` | exactly 0 | **0** (to ~1e−19–1e−37, i.e. machine noise) |
> | `SVAR(1..9)` (Fv) | identity | **[1,0,0, 0,1,0, 0,0,1]** exactly |
> | `SVAR(10)` (α) | 0.05 | **0.050000** exactly |
>
> No nonzero shear (Voigt map correct), compressive not tensile (`Fg` applied
> the right way round), nonzero stress (real build, not the stock stub) — all
> four failure-mode checks below pass by construction. This is the first time
> the growth branch has been exercised end-to-end inside a real ANSYS solve,
> not just the crosscheck harness's standalone driver. Full output in
> [`out.txt`](out.txt) / [`growth_result.txt`](growth_result.txt).

> 🐛 **Bug found and fixed, 2026-08-19: the viscous rows silently didn't match
> until this deck bug was fixed.** `TBDATA,1,1.0,0.0,0.0, 0.0,1.0,0.0,
> 0.0,0.0,1.0` tries to set all 9 `Fv` components in **one** `TBDATA` call —
> but APDL's `TBDATA` only accepts **6 data values per call**
> (`TBDATA,STLOC,C1,...,C6`). Values 7–9 (the last row of `Fv`, including the
> `Fv(3,3)=1` diagonal entry) silently never got set, leaving `ustatev(9)=0`
> by default — i.e. the material actually received a **singular** prior `Fv`
> (`[[1,0,0],[0,1,0],[0,0,0]]`), not the identity `TBDATA` claimed to set.
>
> Confirmed by instrumenting `usermat_biofilm.f` with debug `WRITE` statements
> (temporary, not committed) showing `FV_OLD` diag `= (1, 1, 0)` at the actual
> call site. For `η=0` (elastic) this happened to not matter — the elastic
> branch's answer doesn't depend on `Fv`. For `η>0` (viscous) it made the
> material silently return the **same answer as the elastic case**, which is
> what the two viscous rows below showed before the fix (bit-identical to
> their elastic counterparts — a genuinely deceptive failure signature, now
> added to the table below). **Fixed in the committed deck** by splitting
> into two `TBDATA` calls (`TBDATA,1,...` for components 1–6,
> `TBDATA,7,0.0,0.0,1.0` for 7–9) — this is a general APDL gotcha, not
> specific to this material: **any `TB,STATE`/`TB,USER` table needing more
> than 6 values must be split across multiple `TBDATA` calls with the right
> `STLOC`.**
>
> **Full result after the fix — all four `reference_values.json` cases and the
> Step 4 `KEYOPT` sweep, all PASS:**
>
> | Case | α | η | `KEYOPT(1,2)` | Expected `SX` | Got |
> |---|---|---|---|---|---|
> | elastic | 0.05 | 0 | 0 (B-bar) | −1.019275856e−04 | **−0.10193E−003** |
> | elastic | 0.05 | 0 | 2 (enh. strain) | −1.019275856e−04 | **−0.10193E−003** |
> | elastic | 0.05 | 0 | 3 (simpl. enh.) | −1.019275856e−04 | **−0.10193E−003** |
> | viscous | 0.05 | 8e−3 | 0 | −6.963159875e−05 | **−0.69632E−004** |
> | elastic | 0.20 | 0 | 0 | −4.726465185e−04 | **−0.47265E−003** |
> | viscous | 0.20 | 8e−3 | 0 | −1.795039360e−04 | **−0.17950E−003** |
>
> Element formulation is a non-issue here — B-bar, enhanced strain, and
> simplified enhanced strain all recover the closed-form hydrostatic stress
> exactly (residual shear noise ~1e−20, i.e. zero). `out.txt`/
> `growth_result.txt` reflect the `elastic_a005`/`KEYOPT=0` row, regenerated
> against the fixed deck.

### What each failure mode means

This is why a closed-form case is worth the setup — the failures are diagnostic,
not mysterious:

| What you see | Almost certainly |
|---|---|
| **Nonzero shear** | the `VI/VJ` Voigt map is mis-wired — the Abaqus↔ANSYS 5↔6 shear swap |
| **Tensile** (positive) hydrostatic stress | `Fg` applied inverted — `Fg` used where `Fg⁻¹` belongs |
| `Je ≈ 1`, stress ≈ 0 | `α` never reached the material — check `TB,STATE` and that `nStatev ≥ 10` |
| Stress exactly 0 everywhere | the build did not pick up `usermat_biofilm.f` — stock stub still linked |
| Elastic row matches, viscous does not | `TIME` ≠ the reference `dt` — **or** check every `TBDATA` call sets ≤6 values; see the 2026-08-19 bug above (viscous answer lands exactly on the elastic one, not just "close") |

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

> ✅ **Done, 2026-08-19 on IKMHIWI03** — all three formulations (`0`, `2`, `3`)
> recover `SX=SY=SZ=-1.0193e-04` exactly, on the `elastic_a005` case. No
> volumetric locking detected for this constrained-growth load on `SOLID185`.
> See the results table in Step 3 above.

---

## Step 4.5 — Free-growth complementary check (added 2026-08-19)

`t_growth_constrained.dat` fixes every node, so `F = I` is forced regardless
of whether `Fg`/`Fg⁻¹` are wired with the right sign — a sign error and its
inverse could in principle cancel under that much kinematic constraint.
`t_growth_free.dat` removes only the 6 rigid-body modes (minimal 3-2-1, by
node location) and lets the element actually grow under zero traction.
Closed form: `Fe = I`, stress ≡ 0, for **any** `α`, `η` — there is no
external resistance to relax against, so this isn't a new numeric target,
just a second independent way to be wrong.

```bat
"%AWP_ROOT222%\ANSYS\bin\winx64\ANSYS222.exe" -b -custom .\ANSYS.exe ^
    -i t_growth_free.dat -o out_free.txt
```

> ✅ **Done, 2026-08-19 on IKMHIWI03** — `SX=SY=SZ` ≈ **−1.9234e−10**, shear
> ≈ **1.8e−14**, `Fv` diagonal = 1 (identity). Nine orders of magnitude below
> the constrained case's ~1e−4 stress scale — zero to solver tolerance.
> Evidence: `out_free.txt`, `growth_free_result.txt`.

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
