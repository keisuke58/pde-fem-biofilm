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

## 1.9 Before any solve whose numbers you will report

Run this first. It takes a second and it catches something the routine's own
guard deliberately does not:

```
python ansys_usermat/apdl/check_deck.py ansys_usermat/apdl/*.dat
```

It checks only things that are **silent in ANSYS** — no error, no warning, a
stress field that looks entirely normal. Anything the solver already complains
about does not need a script.

Besides the time step it catches: too few `TB,STATE` slots, which leaves α
unread so `Fg = I` and **the solve runs purely elastic while reporting as a
growth run**; α declared but never written, same outcome; and an over-long
`TBDATA` whose tail APDL drops without saying so. The repository's own decks
pass all of these — but they are the decks that will be copied to make new
ones.

The guard in `biofilm_material_v01.f` refuses only `dt/tau > 0.5`, where the
stress changes sign. Below that the answer stays qualitatively right, and how
much accuracy to buy with step size is the analyst's call rather than
something a material routine should force. But "qualitatively right" leaves a
lot of room, and **`t_growth_baseclamped.dat` is currently inside it**:

| step | dt/tau | von Mises error |
|---|---|---|
| initial (`NSUBST` 1st arg = 4) | 0.06 | **−7.5 %** |
| coarsest AUTOTS may take (3rd arg = 1) | 0.25 | **−30 %** |

Nothing warns about either. The third `NSUBST` argument is the *fewest*
substeps, so it sets the **largest** step AUTOTS is allowed to take — the
opposite of the intuitive reading, which is why this is easy to miss.

The deck has not been changed here: its results are compared against stored
output, so tightening it is a decision to take with the comparison in front of
you rather than a silent edit. **For any run whose stress goes into the
thesis, tighten it first** — `NSUBST,200,500,200` puts every permitted step
under 1 % error on this material.

Worth holding next to the deviator-split finding, which does *not* touch von
Mises at all: this one does, and by far more.

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

---

## 5. Reporting back

Whatever happens, capture these so the result is usable without repeating the
work:

- the `ifort --version` / `cl` output from §0;
- the first 20 compiler errors, if any, verbatim;
- the linker errors, if any, verbatim;
- whether `ANSUSERSHARED.BAT` exists under `%AWP_ROOT222%`.

Paste them back and the next step can be worked out from there.
