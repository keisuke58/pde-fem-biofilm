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
predicted in §1.5. The one file that does not compile is
`Ussfin_P21-V21_Conection_Test.F`, blocked on a separate, unrelated issue:

### New blocker — MKL is not installed anywhere on this machine

`Ussfin_P21-V21_Conection_Test.F` needs `mkl_sparse_handle.fi`,
`mkl_spblas.fi`, `mkl_pardiso.fi`, `mkl_service.fi` (Intel MKL PARDISO
Fortran interfaces). Checked exhaustively: this oneAPI install has only
`compiler`, `compiler_ide`, `mpi`, `tcm`, `umf` components — no `mkl`
directory — and ANSYS v222's own tree has no `mkl_pardiso.fi` anywhere
either (checked as thoroughly as `syspar.inc`/`mpif.h`, which *are* both
present under ANSYS's tree). This is the risk §1.3 flagged in advance,
now confirmed real rather than hypothetical. Getting `Ussfin` to compile
needs the Intel MKL component added to this oneAPI install (a system
change, not attempted without asking first) — everything else in the pool
is unblocked without it.

**Tried and confirmed blocked, 2026-09-02:** `winget install --id
Intel.oneMKL` — no applicable installer under `--scope user` (machine-scope
only), and the machine-scope installer download succeeds (internet access
and hash verification both fine, so this isn't a network problem) but the
install itself fails with `0x800704c7` ("the operation was canceled by the
user") the moment it needs to elevate — the UAC prompt has nobody to answer
it in a non-interactive session, the same wall as every other admin-gated
action found on this account this session (VS Code tunnel service
installation, classic ifort 2024.2.1). **Not a bug to keep working around —
a hard requirement for admin rights on this machine**, same category as the
`C:` VSS/disk-space issue. Getting `Ussfin`/MKL compiling here needs either
admin rights (IT/Timo) or building that one file elsewhere (the cluster
Oliver already uses, per §4's fallback) and treating the local v222 port as
covering the other 10 files only.

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

### Decision, 2026-09-02: this is where the v222 port stops on IKMHIWI03

10 of the 11 pool files compile clean (§1.6). The 11th (`Ussfin`, MKL
PARDISO) needs Intel MKL, which needs admin rights this account does not
have (confirmed via a real `winget install Intel.oneMKL` attempt, blocked
at UAC elevation — not a guess). That makes this the natural stopping
point rather than something to keep chasing here:

- **The material code — the part this thesis is actually about — is fully
  portable and proven to compile under v222.** That is real evidence to put
  in front of Oliver/Meisam, not just an assertion.
- **`Ussfin`'s MKL dependency is a linear-solver plumbing concern, not a
  constitutive-model concern** — it belongs wherever the actual coupled
  field solve runs (Oliver's cluster, which already has MKL), not on
  IKMHIWI03 regardless of admin rights.
- So: **no further v222-port work is planned on this machine.** The
  cluster-access ask (already question 1 in the draft to Oliver) is now
  backed by "10/11 files verified portable, only the MKL-dependent glue
  file is cluster-side," which is a stronger, more specific ask than
  "can we get an account."

---

## 5. Reporting back

Whatever happens, capture these so the result is usable without repeating the
work:

- the `ifort --version` / `cl` output from §0;
- the first 20 compiler errors, if any, verbatim;
- the linker errors, if any, verbatim;
- whether `ANSUSERSHARED.BAT` exists under `%AWP_ROOT222%`.

Paste them back and the next step can be worked out from there.
