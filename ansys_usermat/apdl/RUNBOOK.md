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

> ❌ **Attempted 2026-08-19 on IKMHIWI03, blocked at link time.** Only one
> interactive prompt appeared (Wind Turbine Aeroelastic — answered `Y`) and
> the compile of `usermat_biofilm.f` succeeded. The link step then failed:
> ```
> ifmodintr.lib(iso_c_binding.obj) : error LNK2001: unresolved external symbol for_deallocate_handle
> ifmodintr.lib(iso_c_binding.obj) : error LNK2001: unresolved external symbol for_alloc_allocatable_handle
> ANSYS.exe : fatal error LNK1120: 2 unresolved externals
> ```
> Root cause, consistent with the version-mismatch concern flagged in step
> 0b: the installed oneAPI (2025.3) ships the newer LLVM-based `ifx` runtime
> under the `ifort` name; ANSYS 2022 R2's prebuilt libraries were linked
> against the **classic** Intel Fortran runtime, and the two unresolved
> symbols are classic-runtime-only entry points. This oneAPI install has no
> classic `ifort` fallback. Note that this is an ANSYS-Windows platform
> requirement, not a preference — `ANSCUST.BAT` links against ANSYS's
> prebuilt `.lib`s via MSVC `link.exe` with Intel Fortran's ABI; a
> non-Intel Fortran compiler (e.g. `gfortran`) would fail with a different
> set of unresolved symbols, not fewer.
>
> **Attempted fix 2026-08-19, blocked on UAC, not yet completed.** Classic
> `ifort` was last shipped in oneAPI **2024.2.1** (Intel Fortran Compiler
> Classic 2021.13; 2025.0 dropped it entirely). `winget install --id
> Intel.FortranCompiler --version 2024.2.1` finds and verifies the package,
> but the installer requires elevation (per-machine only — `--scope user`
> returns "No applicable installer found") and the resulting UAC consent
> prompt cannot be approved from a non-interactive/remote session. Needs a
> human physically at IKMHIWI03 to click through UAC, then:
> ```powershell
> winget install --id Intel.FortranCompiler --version 2024.2.1 --source winget --accept-source-agreements --accept-package-agreements --silent
> ```
> Once installed, retry the env-init sequence above pointing at the 2024.2.1
> `compiler\...\env\vars.bat` instead of `2025.3`'s, and re-run `ANSCUST.BAT`.

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
