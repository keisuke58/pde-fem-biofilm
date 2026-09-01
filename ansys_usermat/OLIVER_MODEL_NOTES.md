# Oliver's ANSYS model — what is in it

Notes from inspecting two deliveries received from Oliver on 2026-09-01: the
Workbench project archive `BiofilmImplementation.wbpz`, and then the UPF source
pool `Nishioka_Hoechel.zip`. **Neither is committed** — together they are ~26 MB
of binaries and of another group's source, not ours to redistribute. This file
records what they contain so the integration question can be discussed without
re-opening them.

Everything here was read directly from those files — the project deck
(`ds.dat`, `solve.out`) and the Fortran. Where a statement is an interpretation
rather than something the code states outright, it says so.

---

> **Update 2026-09-01 (later the same day): the source arrived.** A second
> delivery (`Nishioka_Hoechel.zip`) contains the full UPF source pool plus the
> same `.wbpz`. The blocking gap below is closed, and everything previously
> marked "Inferred" is now confirmed from the code. See
> **[Source pool](#source-pool-nishioka_hoechelzip)** onwards — that section
> also records the finding that matters most for us: **the ANSYS 2024 R2
> `usermat` signature is not the v222 one this repo is built against.**

## ~~The blocking gap: the USERMAT source is not included~~ (resolved)

The `.wbpz` alone has an **empty** `user_files/` and no `.f`, `.for`, `.f90`,
`.c`, `.obj`, `.dll` or `.bat` anywhere, so the project could not be run from
that archive by itself: without the UPF binary, ANSYS falls back to its own
stub and the run is meaningless (the "stress exactly 0 everywhere" failure
catalogued in [`apdl/RUNBOOK.md`](apdl/RUNBOOK.md)). The follow-up zip supplies
it.

---

## Verified — read directly from the archive

| Item | Value |
|---|---|
| ANSYS release | **2024 R2** (build 24.2, UP20240603) |
| Element | `SOLID185`, 8 integration points, `NLGEOM,ON` |
| Mesh | 18,750 elements / 23,556 nodes |
| Geometry | `PRJ11_TestCube` — a **test cube**, not tooth/implant geometry |
| Material properties | `TB,USER,1,1,1,NONLINEAR` — **one** property, and **no `TBDATA` follows it** |
| State variables | `TB,STATE,1,,100` — **100** |
| User hooks | `USRCAL,USOLBEG,USSFIN` (solution-begin / solution-finish) |
| Workbench material | "Glass, soda lime (common glass)", E = 69 930 MPa, ν = 0.2149 |
| Solve | 1 load step, 11 substeps, dt = 0.1, **converged in 2 equilibrium iterations per substep, 0 errors** |
| Results | `file.rst` present (~9.9 MB) — a completed run |
| Systems | two (`SYS-16`, `SYS-31`); the solved one is `SYS-31` under design point `dp16` |

Model parameters are passed as **APDL parameters** (`*SET`), not through
`TBDATA` — which is why `TB,USER` declares only one property.

### Biofilm parameters already present

```
YOUNG_BIO   = 1000        ! Pa
YOUNG_VOID  = 1.0         ! void/empty-region stiffness
POISSON_BIO = 0.3
MY_BIOSTART1 = 1.0   MY_BIOSTART2 = 0.0      ! biofilm initial condition, two regions
MY_NUTSTART1 = 1.0   MY_NUTSTART2 = 0.0      ! nutrient initial condition
MY_BETA1 = 1.0e-4    MY_BETA2 = 5.0e-2       ! growth rates
MY_DIFF1 = 1.0       MY_DIFF2 = 1.0          ! diffusion coefficients
```

**The scales agree with this repo.** `YOUNG_BIO = 1000` Pa matches
`material_models.E_MAX_PA = 1000`, and `POISSON_BIO = 0.3` matches the ν
default in [`coupling/composition_to_material.py`](coupling/composition_to_material.py).
A two-field (biofilm + nutrient) reaction–diffusion structure is also what
`JAXFEM/` implements for the Klempt model.

### Nonlocal / meshless machinery

```
SEARCH_RADIUS   = -0.1     ! negative => estimated radius
NEIGHBOR_CNT    = 30       ! max neighbours
BETA_STAR       = 2e-1     ! weighting in the W-matrix exponential
WMAT_THRESHOLD  = 1e-8
ASSEMBLE_KEY    = 1        ! 1: Laplace, 2: Dx+Dy+Dz, else: all
MY_NEM_ACTIVE   = 1        ! 1 = CREATE, 0 = READ
```

plus per-node arrays handed from APDL to the USERMAT —
`VGLBNODENORMX/Y/Z` (surface normals), `VGBOUNDARYNODE` (boundary flags),
`VGLNODEICVALUES` — and a built-in verification suite of analytic test fields
selected by `DEBUG_KEY` (`phi = (x²+y²+z²)/6`, `cos(2πx)cos(2πy)cos(2πz)`,
`sin(πx)`, `exp(2x)`, and a boundary-gradient check).

### Glass / laser-processing parameters also present

```
MY_TEMP_MELT              = 2003.15 K
MY_TEMP_GLASS_TRANSITION  = 1683.15 K
MY_LATENTHEAT             = 1.5e11
MY_A_AMORPH / MY_B_AMORPH / MY_T0_AMORPH
MY_A_CRYSTAL / MY_B_CRYSTAL / MY_TRANGE_CRYSTAL
MY_MAX_ETA = 1e15,  MY_MIN_ETA = 1e-7,  MY_LAMBDA_ETA,  MY_ETA_MECH_1/2
```

with an export folder named `120W_Mesh_Dependence`.

---

## Source pool (`Nishioka_Hoechel.zip`)

Received 2026-09-01. Contents: the same `.wbpz`, plus an `ANSYS-Pool/`
directory with the complete UPF source, the objects, the built
`libansuser.so` / `userlib.a`, and the build script. **Not committed** — it is
Oliver's group's code, not ours to redistribute.

### Build mechanism differs from ours

```bash
module load ANSYS/2024.2
module load intel/2023b
export ANS_PATH=.../ANSYS/2024.2/v242/ansys/
./ANSUSERSHARED_Userdata_Linux_V03_SMP
```

**`ANSUSERSHARED` on Linux, producing a shared library** (`libansuser.so`) —
not `ANSCUST.BAT` producing a custom `ANSYS.exe` on Windows, which is what
[`apdl/RUNBOOK.md`](apdl/RUNBOOK.md) documents for this repo. Built cleanly
(`ansusershared.log`: no compiler or linker errors) with `ifx` / `icc` from
Intel 2023.2.1, on what the module paths show to be an HPC cluster.

Compile order (the log notes ANSUSERSHARED compiles `userdata.F`/`userdata.f`
first to support the common-block feature — worth remembering, since their
file is named `userdata_P21-V21_Conection_Test.f` and so does **not** match
that special-cased name):

```
AceGenElastoAirV08.f  AceGenNeoHookV02.f  AceGenNeoHookV03.f  AceGenNeoHookV04.f
AGPhaseViskoP21V07.f  AGStressP21V07.f    MySubroutines_userData_V04.F
NEM_UserData_P21_V05.F  userdata_P21-V21_Conection_Test.f
Usermat_P21-V21_Conection_Test.F  USolBeg_P21-V21_Conection_Test.F
Ussfin_P21-V21_Conection_Test.F
```

Constitutive routines are **AceGen-generated** (Mathematica symbolic → Fortran;
`sms.h` is the AceGen runtime header).

### ⚠️ The `usermat` signature is release-specific — and it changed

Counted directly from both sources:

| | args | trailing arguments after `cutFactor` |
|---|---|---|
| **2024 R2** (Oliver) | **41** | `pVolDer, hrmflg, var3, var4, var5, var6, var7` |
| **v222** (this repo) | **42** | `var1, var2, var3, var4, var5, var6, var7, var8` |

In 2024 R2 the two reserved slots `var1`/`var2` became named arguments —
`pVolDer(3)` (derivatives of the volumetric potential w.r.t. J: dU/dJ,
d²U/dJ², d³U/dJ³) and `hrmflg` (harmonic-analysis flag) — and `var8` was
**dropped entirely**.

So `usermat_biofilm.f` **cannot be built against 2024 R2 as it stands**: it
declares one argument more than the solver passes, and `var8` would read past
the end of the actual argument list. This is exactly the hazard
[`README.md`](README.md) flags ("the argument list is release-specific — recheck
it first when moving to another ANSYS version"), now confirmed concretely
rather than in principle. Adapting it is a small, mechanical edit, but it must
be done deliberately and the result cannot then run on v222 unchanged.

### Confirmed: the routine is nonlocal, with a parallel data pool

What was inferred from the APDL deck is visible in the source:

- `NEM_UserData_P21_V05.F` (72 kB) — the meshless neighbour/Laplacian machinery.
- `Usermat_*.F` includes `mpif.h` and declares an interface to a global data
  pool: `GetVals` / `SetVals` / `GetTMP` / `SetTMP` / `SetNEM`, all indexed by
  a location `iloc`. State therefore lives partly in `ustatev` and partly in
  this **MPI-shared pool** — hence the directory name `ANSYS-Pool`.
- `USolBeg_*.F` (40 kB) and `Ussfin_*.F` (166 kB) do the per-solution-step
  field work around the per-Gauss-point material call.

### `prop` is unused; state layout

`prop(...)` never appears in the executable body — consistent with
`TB,USER,1,1,1` carrying no `TBDATA`. **All model parameters arrive as APDL
parameters** through the pool, not through the standard material-constant
array.

The state vector in this build:

| slot | meaning |
|---|---|
| `ustatev(2)` | `Sdp_Phi` — **the φ field** |
| `ustatev(3)` | temperature |
| `ustatev(4:6)` | `Vdp_L_n(1:3)` |
| `ustatev(7:9)` | viscosities of the crystalline / amorphous / liquid phases |
| `ustatev(10)` | `Vdp_Cv_n(1)` — right Cauchy–Green of the **viscous part** |
| `ustatev(11)` | density |
| `ustatev(12)` | mechanical viscosity |

A viscous kinematic split is therefore already anticipated in the layout —
`Vdp_Cv_n` plays the role this repo's `Fv` (`ustatev(1:9)`) plays.

### The active material is a biofilm material, not a placeholder

*(This corrects an earlier reading in this file, which had it backwards.)*

In `Usermat_P21-V21_Conection_Test.F` the only **active** constitutive call is

```fortran
CALL AceGenNeoHookV04(Vdp_AceGen, defGrad, stress, dsdePl,
     &   sGdp_YoungBio,    sGdp_YoungVoid,
     &   sGdp_PoissonBio,  sGdp_PoissonVoid,
     &   Sdp_sumBio, Sdp_sumLocal, sedEl, ID)            ! line ~554
```

against the signature

```fortran
SUBROUTINE AceGenNeoHookV04(v, mDefGrad, vCauchy, mTangCC,
     &   sYoung, sYoungL, sNu, sNuL, sBiofilm, sAlpha, sElasticWork, sID)
```

so it is a Neo-Hookean that **blends biofilm against void stiffness**
(`YOUNG_BIO` vs `YOUNG_VOID`, `POISSON_BIO` vs `POISSON_VOID` — the `L`
suffix reads as *leer*) under the local biofilm content. That is a
purpose-built biofilm material, and at 2026-08-05 it is the **newest** routine
in the pool (V02 8 Jan, V03 29 Jul, V04 5 Aug, all with the same signature —
so the interface settled in January and only the internals have been moving).

What is commented out is the *other* model, and the block says whose:

```fortran
!------Matmodell Tobi Start
!      IF(Sdp_Phi .EQ. 0.0D0)THEN
!      !Air Phase
!      CALL AceGenElastoAirV08(...)
!      ELSE !Solid Phase
!      CALL AGStressP21V07(... Sdp_T_n, sGdp_T_Ref, ... vGdp_Th_Expans ...)
!------Matmodell Tobi End
```

— an air/solid split with temperature and thermal expansion, i.e. the **glass**
path. So the disabled branch is the glass model and the live one is the
biofilm model, not the other way round.

### ⚠️ `sAlpha` in that routine is NOT the growth α

Worth stating explicitly, because the name invites exactly the wrong
assumption. The two biofilm arguments are fed from:

```fortran
Sdp_sumBio   = Sdp_bio1_n + Sdp_bio1_n              ! -> sBiofilm
Sdp_sumLocal = (Sdp_locbio1_n + Sdp_locbio2_n)/2    ! -> sAlpha
```

with the declaration commented `!Summe biofilm/local Biofilm`. So `sAlpha`
receives a **local biofilm average**, not a growth variable. There is still no
growth kinematics anywhere in the pool: neither this routine nor
`AGPhaseViskoP21V07` carries an `Fg`, so `Fg = (1+α)I` remains ours to add.

### Possible typo in `Sdp_sumBio`

```fortran
Sdp_sumBio = Sdp_bio1_n + Sdp_bio1_n
```

adds `bio1` to itself, where the very next line correctly averages `locbio1`
and `locbio2`. It reads like `bio1 + bio2` was intended. Flagged as a question
for Oliver rather than a conclusion — it may be deliberate, and we cannot run
the build here to tell.

### The phase/viscous material itself (`AGPhaseViskoP21V07.f`)

Tobi's glass model — the branch commented out of the call chain (see above).
161 kB of Fortran generated by **AceGen 8.103** from a Mathematica notebook of
the same name (2398 formulae, generated 2026-06-09). Worth reading anyway,
because it is the closest thing in the pool to a viscous law and shows the
house style for one. Its interface:

```fortran
SUBROUTINE AGPhaseViskoP21V07(v, vChiN(3), vChiN1(3), vLambdaInit(3),
     vRCGviscoN(6), vRCGviscoN1(6), sPhi, TempN, TempRef, TempMelt, TempG,
     mDefGradN(3,3), mDefGradN1(3,3), mEtaLambda(3,3), mEtaVisco(3,3),
     vYoung(4), vNu(4), vAlphaTh(3), vHeatCapacity(4), vRhoInit(4),
     sKthres, sKpen, sDt, OffsetCal, LatentHeat,
     sNRTolAbs, sNRTolRel, sNRTolStep, sNRResidual, sNRLoops, sNRConverged,
     sRHSTemp, vDebug(100), sNRNmax, sNLSMax, sDebug)
```

Reading it against `BIOFILM_STRESS_CORE`:

| | `AGPhaseViskoP21V07` | `BIOFILM_STRESS_CORE` (this repo) |
|---|---|---|
| Elastic law | E, ν per phase (`vYoung(4)`, `vNu(4)`) | Mooney–Rivlin `C10`, `C01` + `D1` |
| Viscous variable | `vRCGvisco(6)` — right Cauchy–Green of the viscous part, Voigt | `Fv(3,3)` — viscous deformation gradient |
| Viscous update | **local Newton with line search** (`sNRTol*`, `sNRLoops`, `sNRConverged`, `sNLSMax`) | **backward Euler, closed form** — no local iteration |
| Bounded internal vars | `χ` unbounded, mapped through a logistic sigmoid `1/(1+exp(-χ))` | φ, ψ kept in (0,1) by a log-barrier (`hamilton_ode_jax.py`) |
| Volumetric driver | thermal expansion `vAlphaTh` + phase change (melt/glass-transition) | **growth** `Fg = (1+α)I` |
| φ | enters as scalar `sPhi` | not part of the local law |

Two things follow, and they set the size of the porting job:

1. **There is no growth kinematics in his material.** Volume change comes from
   thermal expansion and phase transition, not from a growth tensor. The
   thesis' central mechanism, `F = Fe·Fv·Fg` with `Fg=(1+α)I`, is simply not
   there and would have to be added — this is the substantive part of option
   (A), not a rename.
2. **It is machine-generated, so the edit belongs in the notebook.** Adding
   `Fg` means changing the AceGen source and regenerating; hand-patching
   161 kB of generated Fortran would be overwritten on the next generation and
   is not a maintainable route. Whether we can get the notebook is therefore a
   real prerequisite, not a nicety.

Encouragingly, the *shape* matches: a multiplicative viscous split with the
viscous state carried as a Cauchy–Green tensor is the same kinematic family
this repo verified, and `vRCGvisco` maps to `Fv` through `C_v = Fv^T Fv`.

---

### How the whole thing runs — the solution loop

Traced through `USolBeg`, `Usermat` and `Ussfin`. This is the part worth
understanding first, because it decides where our work can attach.

```
ANSYS solve  (SOLID185, NLGEOM,ON)
 │
 ├─ USolBeg ......... once, at solution start
 │     · ~150 parevl calls: pull every APDL parameter into the
 │       /usercm/ common block   (this is why prop() is unused)
 │     · InitVals            — allocate the shared data arena
 │     · NEM_CreateData_Init — build the meshless neighbour operator
 │     · Mapping_Node_ID, Initial_Values
 │
 ├─ usermat ......... per Gauss point, per equilibrium iteration
 │     · GetVals / GetTMP    — read this point's state from the pool
 │     · CALL AceGenNeoHookV04 → stress, dsdePl     <<< our law would go here
 │     · SetVals             — write state back
 │
 └─ USSFin .......... after each substep
       · assemble + PARDISO solve  →  temperature field
       · assemble + PARDISO solve  →  "Nut1" field
       · assemble + PARDISO solve  →  "Nut2" field
       · (AGPhaseViskoP21V07 — commented out here too)
       · CalcLaserIntegralOMP / CALCPYRO — laser & pyrometer (glass process)
```

So it is a **staggered / operator-split scheme**: ANSYS solves the mechanics,
`USSFin` solves the transport fields on the NEM operator with Intel MKL's
PARDISO sparse direct solver, and the two alternate substep by substep. The
transport fields are *not* ANSYS degrees of freedom — they live entirely in
the UPF's own data pool.

*Interpretation, not something the code states:* the two "Nut" systems are
most likely the biofilm and nutrient fields rather than two nutrients — the
deck carries both `MY_BIOSTART1/2` and `MY_NUTSTART1/2`, with two diffusion
coefficients and two rates, and there are exactly two such solves besides
temperature. The naming looks inherited from the glass code. Worth confirming
rather than assuming.

**Why this matters for us.** `JAXFEM/` solves φ–c–α by finite differences in
JAX; this solves its fields by NEM + PARDISO inside ANSYS. Same operator-split
idea, different discretisation and solver. The mechanical answer is still
produced one Gauss point at a time in `usermat`, which is exactly the shape our
verified law has — so the attachment point is clean even though the surrounding
field machinery is completely different from ours.

---

## How this relates to `usermat_biofilm.f`

The two are **not drop-in compatible** — they cover different scopes.

| | Oliver's USERMAT | `usermat_biofilm.f` (this repo) |
|---|---|---|
| Scope | **Nonlocal** — computes the φ / nutrient fields itself | **Local** constitutive law only |
| Growth driver | Solved internally | Arrives **prescribed** as `ustatev(10) = α` |
| State variables | 100 | 10 (14 with `kStateMat=1`) |
| Properties | APDL parameters; `TB,USER` carries 1 | `TB,USER` carries 6–7 (`prop(1:7)`) |
| Gradients | Meshless neighbour operator (NEM) | none needed — point-local |
| ANSYS release | 2024 R2 | verified on 2022 R2 (v222) |

Two ways they could meet. Now that the source is in hand, (A) has a concrete
address — the `CALL AceGenNeoHookV04(...)` site — but which one is right is
still a supervisor-level decision, not a coding one:

- **(A) Port the constitutive law into Oliver's framework.** His NEM machinery
  supplies φ; this repo contributes the verified `Fg=(1+α)I` growth +
  Mooney-Rivlin/`D1` + viscous response as the mechanical answer, replacing the
  Neo-Hookean stand-in. The 0-ULP Abaqus↔ANSYS equivalence (`crosscheck/`)
  travels with it, and `Vdp_Cv_n` already reserves the viscous state. Under
  this option `usermat_biofilm.f` stops being the deliverable and becomes the
  reference implementation the ported law is checked against.
- **(B) Keep them separate.** This repo's USERMAT stays a local law taking a
  prescribed α from `JAXFEM/`; Oliver's computes its own. They then answer
  different questions and should not be compared directly.

## Open questions for Oliver

1. ~~**The USERMAT Fortran source**~~ — **received** (`Nishioka_Hoechel.zip`).
2. **Which ANSYS release do we target?** His pool is built for **2024 R2 on
   Linux** via `ANSUSERSHARED`; IKMHIWI03 has **v222 on Windows** via
   `ANSCUST.BAT`. The `usermat` signatures differ (41 vs 42 arguments — see
   above), so the two cannot share one source file unguarded. Either we adapt
   to 2024 R2 and work on the cluster, or he confirms a v222 build is viable.
3. **Who computes φ?** His NEM solves the field, which makes this repo's
   α-field mapping (`ustatev(10)`) redundant under option (A) and points at
   (A) as the real integration path — but that is his call, not an inference
   we should act on unilaterally.
4. **How far is the biofilm material meant to go?** `AceGenNeoHookV04` is
   live, recent and biofilm-specific, but purely elastic — no viscosity, no
   growth. Is a viscous biofilm law planned (which is what we would bring), or
   is elastic the intended scope? Related: is the `Sdp_bio1_n + Sdp_bio1_n`
   above a typo?
5. **Is AceGen required to modify the material?** The constitutive routines are
   machine-generated from Mathematica. If the `.nb` source is the real master,
   hand-editing the generated `.f` would be overwritten — we would need either
   the notebook or an agreed hand-written call-out point.

---

*Both deliveries inspected 2026-09-01. Internal file paths and cluster
usernames from the originals are deliberately not reproduced here — this
repository is public, and the source pool is another group's code.*
