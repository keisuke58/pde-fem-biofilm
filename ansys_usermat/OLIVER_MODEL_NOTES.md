# Oliver's ANSYS model (`BiofilmImplementation.wbpz`) — what is in it

Notes from inspecting the Workbench project archive received from Oliver
(2026-09-01). The archive itself is **not committed** — it is a ~10 MB binary
that is not ours to redistribute. This file records what it contains so the
integration question can be discussed without re-opening it.

Everything under "Verified" was read directly out of the archive
(`ds.dat`, `solve.out`, the file listing). Everything under "Inferred" is
read *from* those facts and is flagged as such — the Fortran that would settle
it is not in the archive.

---

## ⚠️ The blocking gap: the USERMAT source is not included

`Biofilm-Implementation_files/user_files/` is **empty**, and the archive
contains no `.f`, `.for`, `.f90`, `.c`, `.obj`, `.dll` or `.bat` anywhere.

The project calls a user material (`TB,USER` + `USRCAL,USOLBEG,USSFIN`), so
**it cannot be run as shipped**: without the custom `ANSYS.exe` built from that
Fortran, ANSYS falls back to its own stub and the run is meaningless. (The
stock-stub failure mode is already catalogued in
[`apdl/RUNBOOK.md`](apdl/RUNBOOK.md): "Stress exactly 0 everywhere".)

**This is the first thing to ask Oliver for.**

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

## Inferred — consistent with the above, but not confirmable without the source

1. **The USERMAT is nonlocal: it solves field PDEs at the integration points.**
   A neighbour search with a Laplace assemble key, a weighted W-matrix, node
   normals, boundary flags, and analytic test fields for a gradient/Laplacian
   operator only make sense if a diffusion-type equation is being solved
   *inside* the material routine. 100 state variables and the
   `USOLBEG`/`USSFIN` global hooks fit the same reading.

2. **It is a glass/laser-processing framework with biofilm work grafted on.**
   Soda-lime glass as the Workbench material, melt and glass-transition
   temperatures, latent heat, amorphous/crystalline kinetics, a viscosity
   ladder up to 1e15, and a `120W` (laser power) export folder are not
   biofilm quantities. The biofilm parameters look like a newer layer on an
   existing codebase.

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

Two ways they could meet, both of which are a supervisor-level decision rather
than a coding one:

- **(A) Port the constitutive law into Oliver's framework.** His NEM machinery
  supplies φ; this repo contributes the verified `Fg=(1+α)I` growth +
  Mooney-Rivlin/`D1` + viscous response as the mechanical answer. The
  0-ULP Abaqus↔ANSYS equivalence (`crosscheck/`) travels with it.
- **(B) Keep them separate.** This repo's USERMAT stays a local law taking a
  prescribed α from `JAXFEM/`; Oliver's computes its own. They then answer
  different questions and should not be compared directly.

## Open questions for Oliver

1. **The USERMAT Fortran source** — nothing runs without it.
2. **Will it build under 2022 R2 (v222)?** IKMHIWI03 has v222 only, and the
   `usermat` argument list is release-specific — see the warning in
   [`README.md`](README.md), which this repo already had to pin down for v222.
3. **Who computes φ?** If his NEM solves the field, this repo's α-field mapping
   (`ustatev(10)`) is redundant and option (A) above is the real integration
   path. If not, (B) applies.
4. **Is the biofilm layer finished or in progress?** The parameters are present
   and the run converges, but the surrounding glass/crystallisation machinery
   suggests the biofilm path may not yet be the primary code path.

---

*Archive inspected 2026-09-01. Internal file paths and cluster usernames in
the original are deliberately not reproduced here — this repository is public.*
