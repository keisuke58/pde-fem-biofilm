# Composition coupling in Abaqus (Keio continuation)

The ANSYS runs put the point model (Klempt et al. 2026) and Eq. 36 (Klempt et al.
2024) into the partner's element through the call-site fragments in
`ansys_usermat/apdl/callsite/`. This folder puts the **same fragments, unchanged**,
into an Abaqus UMAT, so that the model can move to Abaqus at Keio without being
re-implemented.

| file | what |
|---|---|
| `make_umat.py` | writes one user-subroutine file: the Python bridge module, the per-point cache, Eq. 36, the stress routine shared with ANSYS (`biofilm_material_v01.f` + `BIOFILM_STRESS_CORE`), and a UMAT that includes the two fragments |
| `make_inp.py` | one-element check deck (C3D8, phi and c as field variables 1 and 2, clamped or free) |
| `run_comp.ps1` | IKMHIWI03: compiles the C shim, links it through a job-local `abaqus_v6.env`, starts the material server on port 8766 (ANSYS uses 8765), runs the job in `F:\abaqus_work\comp_<job>` |

Inputs: field variable 1 = phi (amount of biofilm), field variable 2 = nutrient c;
constants 1-42 = the ANSYS prop layout, 43-46 = E, E_void, nu, nu_void;
`*DEPVAR` 100 (72-83 point model, 84 alpha - 1, 85-93 F_v). In the ANSYS runs
phi comes from the partner's field; in Abaqus it has to be supplied (a field
variable here; a phi solver or a coupled field is the next step at Keio).

Differences from the ANSYS build, all in `make_umat.py`/`run_comp.ps1`: the
fragments' trace units 95-99 become 195-199 (Abaqus keeps units below 100) and are
opened in the job's output directory (Abaqus runs in a scratch directory); the
Fortran is compiled with `/libs:dll /threads` (the trace I/O otherwise pulls in the
static Intel runtime next to the DLL one Abaqus links); `dev-env.ps1` is not used
(its `C:\msys64\usr\bin` shadows `find`/`link` and breaks `vcvars64.bat`).

## Checked on IKMHIWI03, 5 Oct 2026 (Abaqus 2024)

Case 6, phi = 1 (0.9 in the point model), s = 0.15, T* = 1, dt = 0.025, one element:

| check | result |
|---|---|
| `judge_comp_trace.py` on the Abaqus trace (once per increment, carried, amount, one call, Eq. 36, server replay bit for bit, stand-alone scheme) | PASS, difference 0.0 |
| share phi_1/(phi_1+phi_2) at T* = 1 | 0.0209, the ANSYS value (elem 220) |
| c = 0.3, c_ref = 1 | 0.352; Abaqus and the gfortran driver agree to 3e-16 |
| alpha - 1 | 1.0e-3 = k_alpha phi T* |
| clamped: stress | -5.0624e-7 MPa = closed form + the spherical term of `DEVIATOR_SCALING_FINDING.md` (core shared with ANSYS, not fixed by decision) |
| free: stress | ~1e-24 (none) |

## phi solved by Abaqus: the partner's cube (5 Oct 2026)

`make_cube_inp.py` builds the partner's 2 mm cube as a coupled
temperature-displacement model (C3D8T): phi is the temperature, `UMATHT` (in
`make_umat.py`) solves Eq. 34 as the partner's element does (front term inactive,
its penalty keeping phi in [0, 1]), the UMAT takes phi from the temperature
(constant 47 = 1) and integrates Eq. 36. Same seed elements, beta = 0.02,
k_alpha = 1e-3, dt = 0.025, T* = 1.1, Klempt 2024 stiffness. `compare_ansys.py`
prints the measures of `summarize_runs_json.py` next to the ANSYS run.

The partner's phi lives at the integration points: at t = 0 it is 1 at the seed
elements' points and 0 elsewhere, a jump an FE nodal field cannot hold. Two
nodal versions bracket it: 1 at every node of a seed element (`--ic nodes`, too
much phi: seed p 2.0x ANSYS on 8^3 and 16^3) and the share of seed elements
around each node (`--ic fraction`). With `fraction` the two programs approach
each other as the mesh is refined (ratio Abaqus / ANSYS):

| measure | 8^3 | 16^3 | 24^3 |
|---|---|---|---|
| seed mean stress p | 0.64 | 0.83 | 0.91 |
| seed von Mises | 0.51 | 0.84 | 0.92 |
| layer-1 p | 0.57 | 0.81 | 0.90 |
| alpha - 1, seed surface | 0.74 | 0.91 | 0.95 |
| layer 1 in tension | 100 % / 100 % | 100 % / 100 % | 100 % / 100 % |

So the seed values on coarse meshes are set largely by how the initial seed is
represented; the ring of tension is the same in both programs on every mesh.
24^3 in Abaqus: 6 minutes on one core.

## Two species and the nutrient field (5 Oct 2026)

`make_cube_inp.py --case 2sp_case6` runs the composition mode (prop(28) = 7) with
phi from the temperature. `--cons g` overlays a user element (`UEL` in
`make_umat.py`, DOF 12 = c, DOF 11 = phi read only) that solves Eq. 35 as the
partner's element does (quasi-static, zero order: d lap c = g phi), c = 1 held
on the NUTRIENT1 layer (y <= -0.75 mm); the UMAT reads c at its integration
point from a shared module (one iteration later) and keeps it in SDV51.
Step `UNSYMM=YES` (the c-phi coupling block).

8^3, case 6, beta = 0.02, s = 0.15, against ANSYS:

| | Abaqus | ANSYS |
|---|---|---|
| no nutrient: seed share mean | 0.061 | 0.049 |
| consumption 6: c in the seed | 0.33-0.63 | 0.26-0.56 |
| consumption 6: mean c | 0.63 | 0.59 |
| consumption 6: seed share | 0.08-0.34 | 0.12-0.46 |

The same behaviour (the local nutrient lifts species 1 from ~0.05 to 0.1-0.4);
the 8^3 differences are of the size of the one-species ones (seed
representation). 16^3: `F:\abaqus_work\_abq16_1005.log`.

## In progress (5 Oct 2026, evening; logs in F:\abaqus_work)

- **Four species**: one element, `4sp_case1`, prop(37) = 4: the trace equals the
  stand-alone scheme (n = 4) exactly (difference 0.0); end shares
  0.415 / 0.258 / 0.176 / 0.151, sum phi = 0.9 = phi_cap. Cube 8^3: `_abq4sp_1005.log`.
- **Front term of Eq. 34 as printed** (UMATHT constants 4-6: r, k, h; c and
  grad c from the nutrient UEL; streamline diffusion v h / 2):
  `make_klempt_inp.py` builds Klempt 2024's own fig7_high case (20 um cube,
  20^3), `compare_klempt.py` compares the domain means with the
  finite-difference reproduction (`JAXFEM/klempt2024_quantitative.py`,
  growth="printed": mean phi 0.030, mean c 0.972 at T* = 0.2). Running.
- **Parallel**: one species, 16^3, cpus = 4 gives the same printed values as one
  CPU (largest relative difference 0.0), in 0.6 min. `run_comp.ps1 -Cpus`.
  With the nutrient UEL the first version crashed on 4 CPUs (the store was
  grown on demand) and then differed from 1 CPU by 0.6 % (UMAT read c from
  whichever iteration the UEL had reached). Now the store is allocated once in
  `UEXTERNALDB` (LOP = 0; size `BIOFILM_NUT_NEL`, default 200000 elements) and
  UMAT/UMATHT read the c and grad c committed at the start of each increment
  (LOP = 1, the converged previous increment, as the partner's element hands
  the material the c of the previous sub-step): two species with consumption
  6 on 8^3, 4 CPUs = 1 CPU to all printed digits.
- **Front term, well-posed check** (`--case advect`, `check_advect.py`): c held
  1 / 0 on bottom / top, no consumption, ball of phi = 1 (radius 3 um) at the
  centre. The phi centroid moves down as the ODE dz/dt = -r c/(k+c) says:
  by T* = 0.1, 3.41 um in Abaqus (20^3), 3.65 um finite differences, 3.70 um
  ODE; Abaqus starts one increment late (grad c committed) and then runs ~7 %
  slower on this coarse mesh (ball 6 elements across). 40^3 running.
- **Abaqus-only mesh series** (done; one species, beta = 0.02, dt = 0.025,
  T* = 1.1, the ANSYS 8^3 seed region, `--ic fraction`, 6 CPUs), Pa:

  | mesh | seed von Mises | seed p | layer-1 p | minutes |
  |---|---|---|---|---|
  | 8^3 | 2.48e-4 | -1.61e-4 | 3.43e-5 | 0.4 |
  | 16^3 | 4.91e-4 | -2.11e-4 | 4.87e-5 | 0.6 |
  | 24^3 | 5.55e-4 | -2.24e-4 | 5.24e-5 | 1.9 |
  | 32^3 | 5.79e-4 | -2.28e-4 | 5.38e-5 | 5.4 |
  | 40^3 | 5.90e-4 | -2.30e-4 | 5.44e-5 | 16.1 |
  | h -> 0 (second order, 24/32/40) | 6.09e-4 | -2.34e-4 | 5.55e-5 | |
  | ANSYS 24^3 | 6.06e-4 (-0.5 %) | -2.45e-4 (+5 %) | 5.85e-5 (+5 %) | ~120 |

  The differences 24->32->40 shrink by 2.1-2.2, as second-order convergence
  predicts. ANSYS stops at 24^3 (32^3 out of memory); against the extrapolated
  values its 24^3 seed von Mises is within 0.5 % and its mean stresses ~5 % high
  (its seed is a jump at the integration points). 48^3 (470,596 equations,
  unsymmetric direct solver) stopped in the first solve, most likely out of
  memory; `*SOLUTION TECHNIQUE, TYPE=SEPARATED` would split the system.
- **Front term, Klempt 2024 fig7_high** (done, 20^3, dt 0.002, T* = 0.2): the
  domain mean phi agrees with the finite-difference reproduction at T* = 0.04
  (0.0061 / 0.0068) and then falls behind (0.012 / 0.030 at 0.2). In both the
  seed is eroded (max phi 0.11 / 0.17) and phi spreads thinly above it. Not a
  usable check of the implementation: Eq. 34 as printed moves phi with the full
  speed r c/(k+c) along n_c = grad c / |grad c|, and above the seed, where the
  nutrient is hardly consumed, grad c ~ 0 and its direction is set by
  discretisation noise (finite differences there, shape functions here). The
  same ill-posedness as in the partner element (direction undefined where the
  nutrient is uniform). A check needs a case with a well-defined grad c
  everywhere (e.g. c held 1 / 0 on two faces, no consumption) or a regularised
  n_c.
- **16^3 two species against ANSYS** (done, ratio Abaqus / ANSYS, 8^3 -> 16^3):

  | measure | 8^3 | 16^3 |
  |---|---|---|
  | no nutrient: seed share mean | 1.25 | 1.05 |
  | no nutrient: layer-1 share | 0.81 | 0.96 |
  | consumption 6: seed share mean | 0.67 | 0.85 |
  | consumption 6: seed share min | 0.70 | 1.01 |
  | consumption 6: layer-1 share | 0.77 | 0.97 |
  | consumption 6: c seed min / mean | 1.30 / 1.06 | 1.13 / 1.02 |

  Composition and nutrient approach ANSYS under refinement as the one-species
  measures do.
- **Four species, 8^3 cube** (done): seed phi sum 0.447, shares
  0.395 / 0.268 / 0.187 / 0.150.

## Note

Abaqus shortens the last increment by ~5e-16 to land on T* = 1; the judge now
takes the point model's step per increment from the trace (it assumed a constant
step), which is what made the stand-alone comparison exact.
`tests/test_abaqus_composition_umat.py` repeats the share/alpha/stress checks
with gfortran (no Abaqus needed).
