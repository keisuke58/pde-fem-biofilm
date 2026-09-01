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

### 1.1 The `usermat` argument list — mandatory

Counted from both sources:

| release | args | trailing arguments after `cutFactor` |
|---|---|---|
| 2024 R2 (Oliver's) | **41** | `pVolDer, hrmflg, var3, var4, var5, var6, var7` |
| **v222 (this PC)** | **42** | `var1, var2, var3, var4, var5, var6, var7, var8` |

In 2024 R2 the reserved slots `var1`/`var2` became named arguments `pVolDer(3)`
and `hrmflg`, and `var8` was dropped. So in
`Usermat_P21-V21_Conection_Test.F`, change the subroutine statement to the
v222 form:

```fortran
      subroutine usermat(
     &   matId, elemId, kDomIntPt, kLayer, kSectPt,
     &   ldstep, isubst, keycut,
     &   nDirect, nShear, ncomp, nStatev, nProp,
     &   Time, dTime, Temp, dTemp,
     &   stress, ustatev, dsdePl, sedEl, sedPl, epseq,
     &   Strain, dStrain, epsPl, prop, coords,
     &   var0, defGrad_t, defGrad, tsstif, epsZZ, cutFactor,
     &   var1, var2, var3, var4, var5, var6, var7, var8)
```

and declare `var1, var2, var8` as unused `double precision` scalars alongside
the existing `var3..var7`. **Do not** keep `pVolDer`/`hrmflg`: v222 does not
pass them, and `pVolDer` is an array where `var1` is a scalar.

Take the exact v222 list from
[`../usermat_biofilm.f`](../usermat_biofilm.f), which is the version verified
against this machine's ANSYS.

Check the routine body does not *use* `pVolDer` or `hrmflg` before deleting
them:

```bat
findstr /n "pVolDer hrmflg" Usermat_P21-V21_Conection_Test.F
```

If they are only in the signature, the edit is purely mechanical.

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
