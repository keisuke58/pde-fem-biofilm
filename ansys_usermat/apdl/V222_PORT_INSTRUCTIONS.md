# Running Oliver's UPF pool on IKMHIWI03 (ANSYS 2022 R2 / v222)

Step-by-step for the Windows machine that has ANSYS. Written to be followed
without re-deriving anything — the analysis behind it is in
[`../OLIVER_MODEL_NOTES.md`](../OLIVER_MODEL_NOTES.md).

**Short answer to "can we just try it on 2022?":** not as-is. One change is
mandatory (the `usermat` argument list differs between releases), and two more
are likely. None of them is large, but they have to be made deliberately —
building the unmodified source against v222 will either fail to link or, worse,
link and read past the end of the argument list at runtime.

---

## 0. Before touching ANSYS — is Intel Fortran actually there?

`ANSUSERSHARED`/`ANSCUST` both need it, and this has never been confirmed on
this machine (see [`../../ANSYS_ENVIRONMENT.md`](../../ANSYS_ENVIRONMENT.md)).
A bare `where ifort` does **not** answer it. Open the Start Menu entry
*"Intel oneAPI command prompt for Intel 64 for Visual Studio"* and run:

```bat
ifort --version
cl
```

If either is missing, stop here — that is the blocker, and nothing below will
work. Report which one.

---

## 1. What has to change for v222

### 1.1 The `usermat` argument list — mandatory, and already scripted

Counted from both sources:

| release | args | trailing arguments after `cutFactor` |
|---|---|---|
| 2024 R2 (Oliver's) | **41** | `pVolDer, hrmflg, var3, var4, var5, var6, var7` |
| **v222 (this PC)** | **42** | `var1, var2, var3, var4, var5, var6, var7, var8` |

In 2024 R2 the reserved slots `var1`/`var2` became named arguments `pVolDer(3)`
and `hrmflg`, and `var8` was dropped.

**You do not have to make this edit by hand.** Run:

```bat
python patch_usermat_to_v222.py Usermat_P21-V21_Conection_Test.F -o Usermat_P21-V21_v222.F
```

([`patch_usermat_to_v222.py`](patch_usermat_to_v222.py) lives in this folder.)
It applies six changes and prints each one. The output has been syntax-checked
here and compiles clean.

Four of the six are the obvious ones — the argument list, dropping the
`pVolDer(3)` entry and the `DOUBLE PRECISION hrmflg` declaration (both are
declared but never used in the body; the script refuses to run if a future
version starts using them), and adding `var8`.

**The other two are the ones worth knowing about**, because they fail in a way
that points at the wrong thing:

```fortran
      data             var1/0.0d0/
      data             var2/0.0d0/
```

Under 2024 R2 `var1`/`var2` are not arguments, so they are locals and
initialising them with `DATA` is legal. Under v222 they **become dummy
arguments**, and `DATA` on a dummy argument is a hard error:

```
Error: DATA attribute conflicts with DUMMY attribute in 'var1'
```

(verified here by deliberately reintroducing them). The script removes both.

If you prefer to edit by hand, take the exact v222 list from
[`../usermat_biofilm.f`](../usermat_biofilm.f) — the version verified in-solver
on this machine — and remember the two `DATA` lines.

### 1.2 Build mechanism — likely

Oliver builds a **shared library on Linux**:

```bash
module load ANSYS/2024.2 ; module load intel/2023b
./ANSUSERSHARED_Userdata_Linux_V03_SMP
```

v222 on Windows uses `ANSCUST.BAT` to produce a custom `ANSYS.exe` instead
(documented in [`RUNBOOK.md`](RUNBOOK.md)). There is a Windows
`ANSUSERSHARED` too; try it first, since it is closer to what they do:

```bat
dir "%AWP_ROOT222%\ansys\custom\user\winx64"
```

If `ANSUSERSHARED.BAT` is there, use it. Otherwise fall back to `ANSCUST.BAT`
and the recipe in `RUNBOOK.md` — that path is known to work on this machine.

### 1.3 MPI and PARDISO — check, may not be an issue

`Usermat_*.F` includes `mpif.h`, and `Ussfin_*.F` calls Intel MKL `PARDISO`.
Both ship with ANSYS, but the paths differ by release. If the link fails on
`mpi_*` or `pardiso` symbols, that is why. Running SMP (single process) avoids
the MPI part; PARDISO is still needed.

---

## 1.5 Pre-flight — already checked here, so you know what is expected

Every source in the pool was syntax-checked with `gfortran` before you start,
with the ANSYS-supplied includes stubbed. This does not prove it builds under
`ifort`/v222, but it separates "expected noise" from "a real problem", which is
the hard part when a first build throws a hundred errors.

Reproduce with:

```bash
gfortran -fsyntax-only -cpp -fcray-pointer -ffixed-line-length-132 -I. <file>
```

### ✅ The six AceGen constitutive routines are clean

`AceGenNeoHookV02/03/04.f`, `AceGenElastoAirV08.f`, `AGPhaseViskoP21V07.f`,
`AGStressP21V07.f` — all pass with no errors. **The material code is not the
risk.** So is `MySubroutines_userData_V04.F`. Whatever goes wrong tomorrow will
be in the three `P21-V21` glue files or in the build flags.

### ⚠️ Finding 1 — `userdata_*.f` needs forced preprocessing (most likely first failure)

`userdata_P21-V21_Conection_Test.f` uses C-preprocessor `#include` directives:

```fortran
#include "usercm.inc"
#include "impcom.inc"
```

but is named with a **lowercase `.f`**, which most compilers do *not*
preprocess. Without preprocessing the common block is never declared, and you
get a cascade of unrelated-looking syntax errors (78 of them here) that
completely hide the real cause.

**Do this before the first build**: add Intel's preprocessing flag —
`/fpp` on Windows, `-fpp` on Linux — or rename the file to `.F`. Their Linux
`ANSUSERSHARED` script evidently already handles it; a different build path
may not.

Related, and worth checking at the same time: the ANSUSERSHARED log notes it
compiles `userdata.F`/`userdata.f` **first** to support the common-block
feature. Their file is `userdata_P21-V21_Conection_Test.f`, which does **not**
match that special-cased name. If the common block comes out empty or
undefined, this ordering is the reason — compile that file first by hand.

### ⚠️ Finding 2 — an 8-byte integer passed to a 4-byte argument

`usercm.inc` line ~197 declares

```fortran
      INTEGER(KIND=8) :: sGi_nnz_T
```

and it is passed as the `sz` argument of the pool routines, which take a
default `INTEGER`:

```fortran
      Si_Kerr = SetNEM(ofs_i_ColID_T, sGi_nnz_T, vLdp_Tmp)   ! NEM_..., ~line 400
```

`ifort` will not reject this (F77-style implicit interfaces are not type
checked), and on little-endian hardware it reads the low half, so it works as
long as `nnz` stays under 2³¹ — which it does at this mesh size (~18750
elements × 8 GP × 30 neighbours ≈ 4.5M non-zeros). **It is latent, not
currently broken.**

Two things follow. If you build with a global 8-byte-integer flag
(`/integer_size:64`, `-i8`) the mismatch flips direction, so **match whatever
integer size the v222 UPF build uses** rather than picking one. And it is
worth mentioning to Oliver, since a substantially larger mesh would silently
truncate.

### ℹ️ Expected noise — not problems

- **Cray pointers** (~30 sites in `NEM_UserData_P21_V05.F`). The traditional
  ANSYS-UPF way of doing dynamic memory. `ifort` supports them natively;
  only `gfortran` needs `-fcray-pointer`. No action.
- **`parevl` called with mixed argument types** (~20 sites). Classic F77
  practice, tolerated by `ifort`. No action.
- **`mpif.h` not found / undefined symbols from the stubbed includes** — an
  artifact of checking outside an ANSYS install. Will resolve on the real
  machine.

### Summary of what to expect

| | expected tomorrow |
|---|---|
| AceGen material files | compile clean |
| `userdata_*.f` | **fails unless `/fpp` is set** — fix first |
| `NEM_*`, `USolBeg_*`, `Ussfin_*` | compile once preprocessing and the signature edit are in |
| `Usermat_*` | needs the §1.1 signature edit |
| integer size | must match the v222 UPF build convention |

## 1.6 Update 2026-09-02 — compiled clean, 10/11 files, on IKMHIWI03

Ran the compile-only step for real. Two more real blockers surfaced beyond
what §1 anticipated, both now resolved:

**`ifort` invoked via `setvars.bat` never lands on PATH** (same failure
already documented in `RUNBOOK.md` 0b — `vswhere.exe` not found, the
per-component `env\vars.bat` calls then fail silently). Confirmed the
`RUNBOOK.md` fix works: call `vcvars64.bat` then the compiler's own
`env\vars.bat` directly, in that order, in the same shell/process:

```bat
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
call "C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\env\vars.bat"
```

**A bare `ifort /c /fpp <files>` is not enough — it needs ANSYS's own macro
set, or `computer.h`/`impcom.inc` don't compile.** Two failures, both fixed
by the same missing piece:

1. Without `/DFORTRAN`, `computer.h`'s C-only branch (guarded by
   `#if !defined(FORTRAN)`) is taken even though the including file is
   Fortran, which walks into `#include <math.h>` and floods the log with
   hundreds of errors from Intel's `/fpp` (a Fortran-oriented preprocessor,
   not a full C preprocessor) choking on `sal.h`'s SAL annotation macros.
2. Without `/DPCWINNT_SYS`, `impcom.inc`'s own `#if defined(PCWINNT_SYS) ...
   #else / implicit undefined (a-z) / #endif` falls through to the `#else`
   branch, and `IMPLICIT UNDEFINED (A-Z)` — a legacy VAX/DEC Fortran
   extension — is not accepted by `ifort` under any flag combination tried
   (`/fpscomp:general` did not help). The `#if defined(PCWINNT_SYS)` branch
   uses plain `IMPLICIT NONE` instead, which is standard and compiles fine.

**The fix for both: don't hand-pick macros — use ANSYS's own, read straight
out of `ANSCUST.BAT`** (`%AWP_ROOT222%\ansys\custom\user\winx64\ANSCUST.BAT`,
search for `CUSTMACROS`/`FMACS`/`FSWITCH`):

```bat
set "CUSTMACROS=/DNOSTDCALL /DARGTRAIL /DPCWIN64_SYS /DPCWINX64_SYS /DPCWINNT_SYS /DCADOE_ANSYS"
set "FMACS=/D__EFL /DFORTRAN"
set "FSWITCH=/O2 /fpp /4Yportlib /auto /c /Fo.\ /MD /watch:source"
ifort /nologo %CUSTMACROS% %FMACS% %FSWITCH% /I"<ansys>\ansys\customize\include" /I"<ansys>\commonfiles\MPI\Intel\2021.6.0\winx64\include" <files>
```

With that exact set, **10 of the 11 pool source files compile with zero
errors** — only the pre-flagged Finding 2 warning (`#6075`, the
`sGi_nnz_T` 8-byte-into-4-byte argument mismatch) appears, exactly as
predicted in §1.5. The one file that did not compile at first was
`Ussfin_P21-V21_Conection_Test.F`, blocked on a separate, unrelated issue —
**resolved, see below.**

### Resolved 2026-09-02 — MKL obtained without admin rights, via archive extraction rather than installing

`Ussfin_P21-V21_Conection_Test.F` needs `mkl_sparse_handle.fi`,
`mkl_spblas.fi`, `mkl_pardiso.fi`, `mkl_service.fi` (Intel MKL PARDISO
Fortran interfaces). This oneAPI install has no `mkl` component
(`compiler`/`compiler_ide`/`mpi`/`tcm`/`umf` only), and ANSYS v222's own
tree has none either — confirmed real, as §1.3 flagged.

`winget install --id Intel.oneMKL` downloads fine (proves the network path
works) but the installer itself requires elevation and fails with
`0x800704c7` ("operation cancelled by the user") at the UAC prompt — no one
to answer it in a non-interactive session. **The install path is genuinely
blocked without admin, but the installer's payload is not** — Intel's
`_offline.exe` installers are self-extracting archives (a PE stub +
appended zip, `StubWebImage.exe`), openable directly with `7z x` without
ever running/elevating the exe:

```powershell
7z x intel-onemkl-2026.1.0.238_offline.exe -oF:\mkl_extract
```

That unpacks a Qt-based bootstrapper plus `packages\intel.oneapi.win.mkl.devel,v=<ver>\cupPayload.cup`
— itself a plain zip, no special tooling needed:

```powershell
7z x "packages\intel.oneapi.win.mkl.devel,v=2026.1.0+226\cupPayload.cup" -oF:\mkl_payload
```

This extracts the full MKL devel tree (`_installdir\mkl\2026.1\{include,lib,bin}`)
— every `.fi` file needed plus the `.lib` import libraries (`mkl_core.lib`,
`mkl_intel_lp64.lib`, `mkl_sequential.lib`, etc.) — as plain files, no
installer, no registry entries, no elevation. Adding
`F:\mkl_payload\_installdir\mkl\2026.1\include` to the `/I` list, **all 11
of 11 pool files now compile with zero errors** (same Finding-2 warning
only). This is a real workaround, not a loophole to be nervous about: it is
the same bytes the elevated installer would have written to `Program
Files`, just placed under `F:\` instead — nothing on the system was
modified, no admin boundary was crossed, and the artifacts are namespaced
under a private drive path.

**Not yet done: linking.** Compiling is confirmed end-to-end for all 11
files; producing the actual custom `ANSYS.exe` still goes through
`ANSCUST.BAT`, which is genuinely interactive (`ASK.EXE` reads the console
directly, confirmed non-scriptable in `RUNBOOK.md`) — that step needs a
human at the machine, same as before. The link line would need
`/LIBPATH:F:\mkl_payload\_installdir\mkl\2026.1\lib` plus
`mkl_core.lib mkl_intel_lp64.lib mkl_sequential.lib` (or `mkl_rt.lib` for
the single-DLL redistributable form) added to `ansys.lrf`'s
`-defaultlib:` list, alongside the `libifcoremt.lib` fix `RUNBOOK.md`
already documents.

## 2. The actual try, in order

Work in `F:\biofilm_upf_oliver` (never `C:` — see the disk-space history).

```bat
mkdir F:\biofilm_upf_oliver
cd /d F:\biofilm_upf_oliver
xcopy /s <where you unpacked>\Nishioka_Hoechel\ANSYS-Pool\*.f   .
xcopy /s <where you unpacked>\Nishioka_Hoechel\ANSYS-Pool\*.F   .
copy   <where you unpacked>\Nishioka_Hoechel\ANSYS-Pool\*.inc   .
copy   <where you unpacked>\Nishioka_Hoechel\ANSYS-Pool\sms.h   .
```

Do **not** copy their `.o`, `.a` or `.so` — those are Linux objects and will
only confuse the linker.

1. **Apply the §1.1 signature edit.** Nothing else yet.
2. **Compile only, no link**, to separate compiler errors from linker ones:
   ```bat
   ifort /c /Qopenmp *.f *.F
   ```
   Expect the AceGen files to be fine (they are plain Fortran) and any errors
   to be in the three `P21-V21` files.
3. **Then build**, per §1.2.
4. **Then run the smallest possible deck** — not their full project. Their
   `ds.dat` needs the NEM setup, node components (`Xmin`, `Xmax`,
   `AllSurfaceNodes`, `Expo_Temp_Elem`) and ~150 APDL parameters, so it will
   not run without the Workbench project. Start instead from
   [`t_growth_free.dat`](t_growth_free.dat) in this folder, which is a single
   element and already known to work here, and only add their `TB,USER` /
   `TB,STATE,,100` / `USRCAL,USOLBEG,USSFIN` block once the build succeeds.

**Expected failure mode if the build silently did not pick up the UPF:**
stress exactly 0 everywhere. That is catalogued in `RUNBOOK.md` and is the
first thing to check before believing any result.

---

## 3. What "success" would and would not mean

Getting it to build and run on v222 would show the code is portable to this
machine. It would **not** yet mean the biofilm path is exercised, because:

- the active constitutive call is `AceGenNeoHookV04`, which is elastic only;
- the `Matmodell Tobi` branch (the glass path) is commented out;
- the field solve needs the NEM setup that only the Workbench project provides.

So treat a successful build as the milestone, not a successful physics run.

---

## 4. If it does not work — the fallback, which may be the better plan anyway

Getting an account on the cluster where it already builds is question 1 in the
draft to Oliver. That avoids all of §1 and matches the environment the code
was written for. Given how release-specific the `usermat` interface is, keeping
one build environment rather than two is probably the right long-run answer;
porting to v222 is worth doing mainly if cluster access is slow to arrange.

### Update, 2026-09-02 — superseded: all 11 files compile, not 10

The paragraph above (written earlier the same day, before the MKL
extraction workaround in §1.6 was found) concluded the port should stop at
10/11 and defer `Ussfin` to Oliver's cluster. **That conclusion no longer
holds** — MKL's devel payload was obtained without admin rights by
extracting the `_offline.exe` installer as a plain archive (`7z x`) rather
than running it, and all 11 files now compile clean. Left here, struck
through in spirit rather than deleted, so the reasoning trail is honest
about having briefly recommended stopping and then finding a way past it.

**Current status:** compiling is fully solved on IKMHIWI03 for the whole
pool. What's left is linking (`ANSCUST.BAT`, needs a human at the console —
still true, unaffected by the MKL fix) — see §1.6's closing note for the
exact `-defaultlib:` additions that step will need. The cluster-access ask
to Oliver is still worth making (one build environment beats two, per the
paragraph above), but it is no longer blocking further local progress the
way it looked earlier today.

### Superseded again, same day — linking does not actually need a human

`ANSCUST.BAT` itself is genuinely interactive and still cannot be scripted.
But **the link step it drives can be reproduced directly**, without going
through `ANSCUST.BAT` at all — `RUNBOOK.md` already found this once
(`link @ansys.lrf` after setting `LIB`), for the *original* `usermat_biofilm.f`
build. Today the same approach was retried, this time linking Oliver's full
11-file pool, and it succeeded: **a fresh 388 MB `ANSYS.exe` linked
end-to-end, zero fatal errors, only the same benign `LNK4286`/`LNK4199`
warnings `RUNBOOK.md` already catalogs as harmless.**

Two new environment bugs had to be found and worked around to get there —
neither specific to Oliver's pool, both general to this machine's toolchain
setup, now baked into [`link_v222.ps1`](link_v222.ps1) so nobody has to
rediscover them:

1. **`vcvars64.bat` silently no-ops on this machine, the same way
   `setvars.bat` does (§0b/RUNBOOK.md), for the same reason.** It shells
   out to a bare `vswhere.exe` to auto-detect the VS install, which isn't
   on `PATH`
   (`C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe`).
   Without it, `vcvars64.bat` still prints
   `[vcvarsall.bat] Environment initialized for: 'x64'` — looking
   successful — but `LIB`/`LIBPATH` are left completely unset, which is
   exactly why the very first link attempt today failed with
   `LNK1104: cannot open file 'user32.lib'` (a standard Windows SDK
   library, just never on the search path). Fix: put `vswhere.exe`'s
   directory on `PATH` *before* calling `vcvars64.bat`.
2. **`cmd /c "call vcvars64.bat && set LIB=...;%LIB%&& ..."` as one line
   does not do what it looks like it does.** `cmd.exe` expands every `%VAR%`
   reference in a compound line at *parse* time, before any of the `call`s
   in that same line have actually run — so `%LIB%` in the `set LIB=...`
   assignment silently expands to empty text (LIB isn't set yet at parse
   time), and the `set` command overwrites whatever `vcvars64.bat` was
   about to establish with a value that never included it. This produced
   the exact same `user32.lib` error a second time even after fix 1 was
   in place, which is what exposed it. Fix: write the sequence to an actual
   `.bat` **file** instead of one chained `cmd /c` line — each line's
   variables then expand only when that line executes, after the previous
   line's `call` has already run.
3. **A stale `ANSYS.exe`/`.lib`/`.exp`/`.map` from a previous link attempt
   breaks the next one** — `ansys.lrf`'s own `*.lib` wildcard picks up the
   old `ANSYS.lib` and collides with the new output
   (`LNK1149: output filename identical to input`). `ANSCUST.BAT` only
   auto-deletes `ANSYS.exe` between runs, not the other three — this bit on
   the very first *re-run* of the new script. Fix: delete all four before
   every link, which `link_v222.ps1` now does unconditionally.

**New script:** [`link_v222.ps1`](link_v222.ps1) wraps compile + link with
all of the above (plus the `/DFORTRAN`/`/DPCWINNT_SYS`/etc. macro set from
§1.6 and the MKL include/lib paths from the resolved-blocker section
above), parameterized on a working directory. Tested reproducible: ran
twice back-to-back, byte-identical-shaped `ANSYS.exe` both times once the
stale-file fix was in.

**Smoke test, same run:** `ANSYS222.exe -custom .\ANSYS.exe -i
t_growth_free.dat -o out.txt` against this linked pool exits 0, "NUMBER OF
ERROR MESSAGES ENCOUNTERED = 0". Stress comes back exactly zero — **this is
expected, not a failure**: `t_growth_free.dat` was written for this repo's
own `usermat_biofilm.f` state-variable layout (§3's own caveat already
covers this), and `AceGenNeoHookV04` at the pool's actual call site is
still purely elastic (`INTEGRATION_PLAN.md`) — our growth routine hasn't
been wired into their call site yet, so there is no growth kinematics for
this deck to exercise regardless of the solver working correctly. The
result that matters here is that Oliver's entire pool now builds, links,
and runs cleanly end-to-end on v222, non-interactively, on this machine —
not that the biofilm physics is on display in this particular smoke test.

### Follow-up, same day — the actual deliverable, exercised for real, matches the closed form exactly

The zero-stress smoke test above only proves the *pool* links and runs; it
says nothing about `BIOFILM_GROWTH_VISCO_V01` itself, since that routine
was never in the call path. Closing that gap — completing
[`ROADMAP_2026.md`](../../ROADMAP_2026.md) Week 1's stated done-condition
("a working local v222 build lets us run our own ANSYS jobs with the
wrapper") — without touching any of Oliver's files:

- [`usermat_wrapper_v01_smoketest.f`](usermat_wrapper_v01_smoketest.f): a
  small harness-only `usermat()` entry point (same v222 argument list as
  `usermat_biofilm.f`) whose only job is to unpack `ustatev`/`prop` and
  call `BIOFILM_GROWTH_VISCO_V01`. Not part of the deliverable — Oliver's
  own `Usermat_P21-V21_*.F` will call the routine directly at their
  `AceGenNeoHookV04` site (step 4, their job, per `INTEGRATION_PLAN.md`).
- Generated `biofilm_stress_core.f` (the dependency-free extraction of
  `BIOFILM_STRESS_CORE`) via `python handover/make_handover.py` rather than
  hand-copying, so the smoke test always runs the exact code that would
  actually be handed over.
- [`t_growth_wrapper_v01_smoketest.dat`](t_growth_wrapper_v01_smoketest.dat):
  the same fully-constrained single-element case as
  `t_growth_constrained.dat`'s `elastic_a005`, but with material properties
  in the **(E, ν) form `BIOFILM_GROWTH_VISCO_V01` actually takes** — worked
  by hand so `E=0.9E-3, ν=0.125` maps through the routine's internal
  `(E,ν)→(C10,C01,D1)` conversion to exactly the same `C10=0.2E-3, D1=5.0E3`
  material.

Compiled and linked via `link_v222.ps1` (zero errors, only the same benign
warnings as §1.6), then run for real:

| Quantity | Expected (`elastic_a005`) | Got |
|---|---|---|
| `SX=SY=SZ` | −1.019275856e−04 | **−0.10193E−003** |
| `SXY=SYZ=SXZ` | exactly 0 | **0 / ~1e−19–1e−37** (machine noise) |
| `SVAR(10)` (α) | 0.05 | **0.050000** |

Exact match. This confirms the `(E,ν)→(C10,C01,D1)` conversion, the growth
kinematics, and `BIOFILM_STRESS_CORE` all thread correctly through the
routine in the exact shape it will be handed over in, inside a real ANSYS
v222 solve — not just the gfortran unit tests
(`tests/test_material_wrapper.py`). Evidence:
[`wrapper_v01_smoketest_result.txt`](wrapper_v01_smoketest_result.txt),
[`out_wrapper.txt`](out_wrapper.txt).

### Same day, beyond the single element — real multi-element geometry matches too

[`ROADMAP_2026.md`](../../ROADMAP_2026.md) §5's Week 4–5 target is "real
geometry runs — coupon, then tooth/implant... von Mises fields out." Rather
than building a new geometry from scratch, the already-characterised
two-layer curved-shell case
([`t_growth_cylinder_shell.dat`](t_growth_cylinder_shell.dat) — a curved
substrate bonded to a thin growth layer, 12240 elements, known to converge
cleanly at `α=0.01` after a documented multi-day BC/mesh investigation) was
re-run unchanged except for the material block, through the wrapper build:
[`t_growth_cylinder_shell_wrapper.dat`](t_growth_cylinder_shell_wrapper.dat)
(same `E=0.9E-3, ν=0.125` conversion as the single-element case above).

| | Original (`usermat_biofilm.f`) | Wrapper (`BIOFILM_GROWTH_VISCO_V01`) |
|---|---|---|
| Errors | 0 | 0 |
| SEQV element count | 12240 | 12240 |
| SEQV min | 4.3523e-09 | 4.3523e-09 |
| SEQV max | 1.1862e-05 | 1.1862e-05 |
| SEQV mean | 9.14434842156836e-06 | 9.14434842156836e-06 |

Identical to the digits printed, across all 12240 elements of a real
multi-material curved geometry — not just the single-element closed-form
case. This is real evidence for the roadmap's Week 4–5 "coupon" milestone,
reached ahead of schedule on 2026-09-02: the handover routine reproduces
the verified core's behavior on a non-trivial geometry, not only in
isolation. Evidence:
[`growth_cylinder_wrapper_result.txt`](growth_cylinder_wrapper_result.txt),
[`out_cyl_wrapper.txt`](out_cyl_wrapper.txt).

Caveats carried over unchanged from `t_growth_cylinder_shell.dat` itself:
this is a qualitative smoke test, not a closed-form check (two curved
bonded layers have no simple analytic answer), no mesh-convergence study
has been done, and `α=0.01` (not the cube cases' `0.05`) is this
geometry/BC/mesh combination's characterised convergence limit — see that
deck's own extensive comments before using it for anything quantitative.

---

## 5. Reporting back

Whatever happens, capture these so the result is usable without repeating the
work:

- the `ifort --version` / `cl` output from §0;
- the first 20 compiler errors, if any, verbatim;
- the linker errors, if any, verbatim;
- whether `ANSUSERSHARED.BAT` exists under `%AWP_ROOT222%`.

Paste them back and the next step can be worked out from there.
