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
| `run_comp.sh` | Linux (Keio server) counterpart of `run_comp.ps1 -Native`, with `mp_mode=threads`; not yet run on Linux (`KEIO_SERVER_HANDOFF.ja.md`) |

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
  by T* = 0.1 (um moved, share of the ODE's):

  | | 20^3 | 40^3 |
  |---|---|---|
  | Abaqus | 3.41 (92 %) | 3.59 (97 %) |
  | finite differences (20^3) | 3.65 (99 %) | |
  | ODE | 3.70 | 3.70 |

  Abaqus starts one increment late (grad c is the committed one) and then
  moves at 93 % (20^3) / 97 % (40^3) of the ODE speed (T* = 0.05-0.1): the gap
  halves with h, the first-order streamline diffusion. 40^3: 130 min on 6 CPUs.
- **Abaqus-only mesh series** (done; one species, beta = 0.02, dt = 0.025,
  T* = 1.1, the ANSYS 8^3 seed region, `--ic fraction`, 6 CPUs), Pa:

  | mesh | seed von Mises | seed p | layer-1 p | minutes |
  |---|---|---|---|---|
  | 8^3 | 2.48e-4 | -1.61e-4 | 3.43e-5 | 0.4 |
  | 16^3 | 4.91e-4 | -2.11e-4 | 4.87e-5 | 0.6 |
  | 24^3 | 5.55e-4 | -2.24e-4 | 5.24e-5 | 1.9 |
  | 32^3 | 5.79e-4 | -2.28e-4 | 5.38e-5 | 5.4 |
  | 40^3 | 5.90e-4 | -2.30e-4 | 5.44e-5 | 16.1 (6.7 separated) |
  | 48^3 (separated) | 5.96e-4 | -2.31e-4 | 5.47e-5 | 16.7 |
  | h -> 0 (second order, 40/48) | 6.10e-4 | -2.34e-4 | 5.55e-5 | |
  | ANSYS 24^3 | 6.06e-4 (-0.7 %) | -2.45e-4 (+5 %) | 5.85e-5 (+5 %) | ~120 |

  The differences 24->32->40 shrink by 2.1-2.2 and 32->40->48 by 1.7 (1.8 for
  second order). ANSYS stops at 24^3 (32^3 out of memory); against the
  extrapolated values its 24^3 seed von Mises is within 1 % and its mean
  stresses ~5 % high (its seed is a jump at the integration points). 48^3 ran
  out of memory with the default unsymmetric coupled solve;
  `make_cube_inp.py --separated` (`*SOLUTION TECHNIQUE, TYPE=SEPARATED` and a
  symmetric solve, possible because phi does not depend on the displacement)
  runs it. On 40^3 it gives the same values as the coupled solve to all printed
  digits, in 6.7 instead of 16.1 minutes.
- **Front term, Klempt 2024 fig7_high** (done, 20^3, dt 0.002, T* = 0.2): the
  domain mean phi agrees with the finite-difference reproduction at T* = 0.04
  (0.0061 / 0.0068) and then falls behind (0.012 / 0.030 at 0.2). In both the
  seed is eroded (max phi 0.11 / 0.17) and phi spreads thinly above it.
  Where the extra phi comes from (finite differences, T* = 0.2): source
  k_alpha alpha 0.0002, clipping to [0, 1] 0.0000, the front term 0.027. The
  front term -v . grad(phi) is not conservative: it adds phi div v, and
  v = r c/(k+c) grad c/|grad c| has a large divergence around the seed, where
  the c isolines bend round the small consuming colony. That is a
  grid-scale quantity on 20^3.
  Regularised direction, n_c = grad c / sqrt(|grad c|^2 + eps^2)
  (`--eps`, also in the reproduction's `run(eps=)`; an assumption of this
  work), mean phi at T* = 0.2 (Abaqus / finite differences): eps = 0
  0.012 / 0.030, 1e-4 0.014 / 0.030, 1e-3 0.011 / 0.020 (|grad c| at t = 0:
  median 1e-5, max 3e-3 per um). eps = 1e-4 switches the term off where c is
  uniform and changes nothing in the reproduction, so the noise direction
  where grad c ~ 0 is not what separates the two. The reproduction from the
  same continuous seed (`--seed-field`, `--fd-n`) gives 0.030 on 21^3 and
  0.056 on 41^3 nodes: the printed term does not converge in the grid, so
  this case cannot check the implementation (the Abaqus 40^3 run ran out of
  memory next to other jobs and was not repeated).
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

## Klempt 2024 test cases 4.1 and 4.2 in Abaqus (5-6 Oct 2026)

The closest reproduction of the paper's figures (KLEMPT2024_REPRODUCTION.md
sec. 15-16) in Abaqus, to check the finite-difference reproduction with an
independent FE solver: UMATHT constant 8 = w (growth on every face,
r c/(k+c)(w |grad phi| + (1 - w)|n_c . grad phi|), w = 0.5), UEL property
4 = 1 (consumption g phi c), Table 2 otherwise with a time scale s per run
(beta, k_alpha, r times s), 4.1 with the nutrient in a 2 um corner block.
All three are departures from the printed equations, stated as a hypothesis.
`make_klempt_inp.py --case fig4_corner|fig7_high|fig7_low --blend 0.5
--first-order --scale s`, `compare_blend.py`. 20^3, RMS difference of the
mean phi / c curves (averaged) to the digitised figures:

| | Abaqus, artificial diffusion v h/2 | Abaqus, none (`--hstab 0`) | finite differences |
|---|---|---|---|
| 4.1 (s = 1) | 0.019 | 0.025 | 0.013 |
| 4.2 high (s = 5) | 0.102 | 0.021 | 0.022 |
| 4.2 low (s = 3) | 0.041 | 0.030 | 0.030 |

The artificial diffusion (first-order upwinding's v h / 2, up to
125 um^2/T* at s = 5 against beta s = 10) smeared the one-layer seed of 4.2 before the
front took it: phi then levels off at the smeared maximum (0.886), since a
growth term proportional to |grad phi| cannot raise a maximum. Without it
the Galerkin solution is stable (beta s keeps the cell Peclet number at ~12)
and all three cases agree with the reproduction to RMS 0.01-0.03.

Table 3 (stress): with the mechanics in the same run (`--E 10 --s-every 25`,
`table3_pressure.py`, `klempt_table3_fig.py`,
`assets/fig_abaqus_klempt2024_table3_nostab.png`) the hydrostatic stress has
the paper's pattern (compressed interior, ring in tension, about zero outside)
and, with E = 10 read as MPa, its size: p = -3.2e-4 ... +1.5e-4 at T* = 0.1
and -1.2e-3 ... +4.9e-4 at 0.25, against the legend's +3e-4 ... -6.5e-4 MPa.
In Pa it would be 1e6 times smaller (one colour in the table), which supports
sec. 13's reading of mu in MPa with the FE mechanics instead of a
small-strain estimate. (p carries the spherical term of
DEVIATOR_SCALING_FINDING.md, ~1.3 % at nu = 0.49.) The shapes agree in kind;
the thin tip towards the corner is missing for every w (sec. 16).

## The point model in Fortran: no material server (6 Oct 2026)

The material server (Python, one socket round trip per integration point and
increment) was ~100x the rest of a run: one species 40^3 (no server) 16 min,
two species on 2880 elements 1.5 h. `ansys_usermat/coupling/ecology_native.f`
is the point model in Fortran, a drop-in for `usermat_py_hook.f` (same module
`biofilm_py_bridge`, same hooks; the fragments are unchanged): the scheme of
`ecology_jax.ecology_substeps` (6 Newton iterations per sub-step on the 12
residuals, clip_state before and after each, the exact Jacobian written out,
LU with partial pivoting). Case constants (n, c*, alpha*, eta_i, theta,
checked against the deck as the server does) from the file in
`BIOFILM_ECO_CASE` (`write_eco_cfg.py`, taken from `material_server.set_case`).

| check | result |
|---|---|
| `tests/test_ecology_native.py`: 16 states x cases (2sp_case6, 2sp_case3, 4sp_case1, none) x sub-steps x c_rel, against ecology_jax | largest relative difference 8e-15 |
| Abaqus, 8^3 two species, consumption 6 (`run_comp.ps1 -Native`) against the server run of the same deck | printed values identical; point-model trace (8213 point-steps) within 6.7e-12; 0.6 min on 1 CPU (with compiling) against 5.2 min on 4 CPUs |
| ANSYS, partner element, week deck `w8_c6_g6_s015_b001` (exe built in `F:\biofilm_upf_native`, ecology_native.f as usermat_py_hook.f, no C shim) against the week chain's server run | all_stress and nut_field identical; trace within 1.6e-12 (absolute); 36 s against 156 s |

The Newton iteration count is fixed at 6 as in the reference: where a step is
too long for 6 iterations to converge, both versions end at the clip bounds
in round-off-dependent places (seen with THETA_DEMO, c* = 25, one step of
0.01 from the default state); the composition runs' steps converge.

## Front term with two species; composition mesh series (6 Oct 2026)

Runs one at a time on 4 CPUs (`F:\abaqus_work\_req3_1006.ps1`; since 11:15
with the Fortran point model, `-Native`), each result in
`results_1006/<job>.txt` (`compare_ansys.py` measures). All: case 6,
beta = 0.02 mm^2/T*, s = 0.15, the seed of the ANSYS 8^3 run, T* = 1.0.

1. Front term with two species and the nutrient, **consumption 1** (Klempt
   2024 Table 2 scaled to the 2 mm cube): `fr8_c6_g1`, `fr16_c6_g1`
   (`make_cube_inp.py --front 10 --blend 0.5`, dt 0.01) and `nf8_c6_g1`
   without the front term. ANSYS reference: the week run `w8_c6_g1_s015`
   (`../ansys_usermat/apdl/results/2026-10-05_week/`, consumption 1, Klempt
   2024 stiffness, no front term; the partner element's front term is not
   active), like for like with `nf8_c6_g1`. `fr8_c6_g6` / `nf8_c6_g6`
   (consumption 6, server build) were the first try, superseded.
2. Composition mesh series, no front term, consumption 6: `cs16_c6_g6`,
   `cs24_c6_g6`, `cs32_c6_g6` (ANSYS column: the 8^3 run
   `ds_c6_nut_g6_b002`, a reference only for the trend).

Results, T* = 1.0 (stress in MPa with E = 10 Pa; seed = the ANSYS seed elements):

| run | seed vM | seed p | alpha-1 interior | phi seed | share seed | share layer 1 | c seed min | c mean |
|---|---|---|---|---|---|---|---|---|
| ANSYS `w8_c6_g1_s015` (8^3) | 4.68e-4 | -2.38e-4 | 8.62e-4 | | 0.051 | 0.334 | 0.876 | 0.932 |
| `nf8_c6_g1` (no front term) | 2.39e-4 | -1.51e-4 | 7.15e-4 | 0.447 | 0.063 | 0.268 | 0.889 | 0.938 |
| `fr8_c6_g1` (front, w = 0.5) | 6.35e-4 | -3.08e-3 | 1.16e-3 | 1.070 | 0.324 | 0.267 | -0.021 | 0.323 |
| `fr16_c6_g1` (front, w = 0.5) | 7.88e-4 | -3.29e-3 | 1.23e-3 | 1.089 | 0.361 | 0.279 | -0.048 | 0.323 |

- Without the front term, Abaqus and ANSYS agree on the nutrient (c within
  1.5 %) and on the composition far from the seed (share in layer 2: 0.4994
  vs 0.4995). The seed stress on 8^3 is 0.51 (vM) and 0.64 (p) of ANSYS:
  the same ratio as the one-species 8^3 run (0.64), which approached ANSYS
  under refinement (0.83 / 0.91 on 16^3 / 24^3). `nf16_c6_g1` and
  `nf24_c6_g1` against the ANSYS week runs at 16^3 (and 24^3 when it is
  done) follow.
- With the front term (r = 10 mm/T*, w = 0.5) the biofilm fills the whole
  cube by T* = 1 (phi about 1.07-1.09 everywhere; Eq. 34 has no upper
  bound on phi). The nutrient is consumed down to c mean 0.32 (minimum just
  below 0), the share at the seed rises from 0.06 to 0.32-0.36, and the
  pressure is about 20 times the run without the front term, because the
  whole cube grows against the fixed base. 8^3 to 16^3 changes seed vM by
  24 % and p by 7 %.

Composition mesh series (consumption 6, no front term):

| | 16^3 | 24^3 | 32^3 | extrapolated |
|---|---|---|---|---|
| seed vM | 4.77e-4 | 5.41e-4 | 5.65e-4 | 6.0e-4 (order 1.8) |
| seed p | -2.00e-4 | -2.13e-4 | -2.17e-4 | -2.23e-4 (order 2.0) |
| share seed mean | 0.250 | 0.256 | 0.258 | |
| share layer 1 | 0.400 | 0.408 | 0.409 | |
| c seed min | 0.278 | 0.266 | 0.262 | |

The composition and the nutrient change by about 1 % from 24^3 to 32^3; the
seed stress by 4 %, 6 % below the extrapolated value at 32^3. (alpha - 1
at the surface elements keeps falling, 5.05/4.88/4.71e-4, because the
surface element's centroid moves closer to the free surface as the mesh is
refined; it is a measure of the element, not of a point.)

40^3 (`cs40_c6_g6`, `--separated`, 20 min on 4 CPUs; the first unsymmetric
solve ran out of memory): seed vM 5.76e-4, seed p -2.19e-4, share seed mean
0.260, share layer 1 0.411, c seed min 0.260. From 32^3 the seed vM rises by
2.0 % and p by 1.1 %; composition and nutrient change by less than 1 %. The
values extrapolated from 16/24/32 lie 4 % (vM) and 1.6 % (p) beyond 40^3.
The ANSYS column in `results_1006/cs40_c6_g6.txt` is the 8^3 run with the
example stiffness (1000 MPa), so its stress ratios mean nothing.

**Felix Klempt's answers (mail of 5 Oct 2026)**, which change how the
front-term runs are labelled:
- Eq. 34 **as printed** is the form in his implementation: where the
  nutrient gradient behind the biofilm points towards it, the biofilm does
  not grow there (more consumption or more nutrient diffusion changes
  that). The growth on every face with w = 0.5 is therefore my
  modification, chosen because it reproduces the paper's figures in my
  code; it is not his model. Runs with the printed form (`pr8_c6_g1`,
  `pr16_c6_g1`) follow.
- grad c/(|grad c| + eps) with eps about 1e-8 to 1e-12 (here
  grad c/sqrt(|grad c|^2 + eps^2), the same for these values).
- Consumption: "I think g phi" (zero order, Eq. 35), to be checked in the
  AceGen file (with Dr. Soleimani). Simulation length 1 in every test
  case. Test case 1: nutrient only at the corner.
- beta: chosen to give sensible results, dependent on the mesh; beta =
  0.02 mm^2/T* stays an assumption.
- Stiffness: the paper's units may be wrong; he suggests E = 10 kPa (mu
  about 3.3 kPa). With growth as the only load the stress is proportional
  to E, so all stresses here scale by 1000 for that value.
- Composition at points the biofilm reaches later: not part of the model;
  his idea is the neighbour's composition or the neighbours' average.

Assumptions, to be stated wherever these results are shown:
- **Front term:** growth on every face, w = 0.5, without artificial diffusion. This is the form that reproduces Klempt 2024's figures (sections above; KLEMPT2024_REPRODUCTION.md sec. 15-16), **not Eq. 34 as printed**. The printed form does not converge in the grid in my FD code. Felix Klempt (5 Oct) confirms the printed form as his implementation, so w = 0.5 is my modification. r = 10 mm/T* and k = 1 are Table 2 converted to the 2 mm cube (KLEMPT2024_REPRODUCTION.md sec. 8).
- **Consumption:** 1 in the front-term runs is Table 2 converted to the 2 mm cube. 6 in the mesh series is an example input value (the partner's deck), about six times Table 2. With the paper's zero-order form, c can drop below 0 where the biofilm is dense (-0.19 on 8^3 at T* = 0.1). The point model then takes c_rel = min(max(c/c_ref, 0), 1) = 0, i.e. no nutrient: `phi_mode_exec.inc`, the same fragment in ANSYS and Abaqus.

## Klempt 2024 test cases with the set-up the first author described (6 Oct 2026)

`make_klempt_inp.py` with Felix Klempt's answers of 5 Oct: Eq. 34 as printed,
n_c = grad c/(|grad c| + eps) with eps = 1e-8, Table 2 with T* = 0..1 (no time
scale per run), the nutrient at the corner only (4.1), no artificial
diffusion (`--hstab 0`), the paper's 1 um mesh (20^3), dt = 1e-3, E = 10 kPa;
phi unbounded (`--pen 0`) or held in [0, 1] (`--pen 100`).
`compare_paper.py`, RMS difference of the mean phi / c curves (averaged) to
the digitised figures (`results_1006/fx*.txt`):

| | 4.1 | 4.2 high | 4.2 low |
|---|---|---|---|
| consumption g phi (as he described), phi unbounded | 0.32 | 0.46 | 0.22 |
| the same, phi held in [0, 1] | 0.32 | 0.54 | 0.24 |
| consumption g phi c, phi unbounded | 0.045 | 0.45 | 0.27 |
| growth on every face (w = 0.5), g phi c, time scale per run (above) | 0.025 | 0.021 | 0.030 |

- With g phi and Table 2 the nutrient is used up at once and goes below 0
  (4.1: mean c 0.004 at T* = 0.05 against the paper's 0.33); the front term
  r c/(k+c) then vanishes and phi hardly grows (0.065 to 0.071).
- With g phi c, 4.1 comes within 0.045 (mean phi 0.61 against 0.74 at
  T* = 1); 4.2 is not reproduced by either form.
- The first author's current implementation (`../ansys_usermat/FELIX_VS_PARTNER_DIFF.md`,
  compared in words only) adds front growth only where it is positive and
  lets it fade as phi approaches 1. None of the runs above has that form; it
  is the next one to test, stated as "as in the first author's current
  implementation". (Done on the paper's mesh, next section.)

## His front form on the paper's mesh: g phi against g phi c (6 Oct 2026, evening)

`make_klempt_inp.py --felix` (his form of the front term as described in his
mail, penalty 100), the paper's 1 um mesh (20^3), Table 2 with T* = 0..1 (no
time scale per run), dt = 1e-3, eps = 1e-8, no artificial diffusion,
E = 10 kPa. Consumption g phi (`ff41`, `ff42h`, `ff42l`) and g phi c
(`--first-order`, `ff*_fo`; the form Dr. Soleimani confirmed for the 2024
figures on 6 Oct). 4 CPUs, 18-27 min each. RMS difference of the mean phi / c
curves to the digitised figures (`compare_paper.py`, `results_1006/ff*.txt`):

| case | g phi: RMS phi, c | g phi c: RMS phi, c | mean phi at T* = 1: g phi / g phi c / paper | lowest element c at T* = 1: g phi / g phi c |
|---|---|---|---|---|
| 4.1 | 0.443, 0.188 | 0.252, 0.094 | 0.071 / 0.339 / 0.740 | -0.143 / 0.126 |
| 4.2 high | 0.709, 0.374 | 0.709, 0.376 | 0.099 / 0.099 / 1.000 | 0.830 / 0.850 |
| 4.2 low | 0.183, 0.306 | 0.174, 0.384 | 0.023 / 0.051 / 0.301 | -0.589 / 0.054 |

- With g phi c the nutrient stays positive in all three cases, and 4.1 moves
  towards the paper (RMS phi 0.44 to 0.25; mean phi 0.34 against 0.74 at
  T* = 1).
- 4.2 high does not change: c stays near 0.9, so the consumption form hardly
  acts, and phi reaches 0.10 against the paper's 1.0.
- 4.2 low: phi a little closer, c further away (mean c 0.21 against the
  paper's 0.08).
- 4.2 is not reproduced with his front form and Table 2 under either
  consumption form.

`ff41_g5` (4.1 with five times the consumption) also ran: c falls to about
-4 at once and phi stays at 0.065. It is not part of the comparison.

## The same set-up against the first author's own code (6 Oct 2026)

The three test cases were also run in Felix Klempt's current code (ANSYS, kept
outside git; `../ansys_usermat/FELIX_VS_PARTNER_DIFF.md`). The ANSYS runs use
8^3 elements and the regions of the partner-style Fig. 7 deck. The Abaqus runs
use the same mesh and regions (`partner_elem_sets.py` writes `pset_41.json` and
`pset_42.json`; `make_klempt_inp.py --elem-sets`):
- c = 1 on every node of the nutrient elements;
- initial phi = the share of seed elements around each node.

Common settings: his form of the front term (`--felix`), consumption g phi,
dt = 0.005, T = 1, E = 10 kPa, penalty 100. Jobs fb41, fb42h, fb42l.

| case | mean phi at T* = 0.5: Abaqus / his code / paper | mean c at T* = 0.5: Abaqus / his code / paper | RMS Abaqus - his code (phi, c) | RMS Abaqus - paper (phi) |
|---|---|---|---|---|
| 4.1 | 0.097 / 0.069 / 0.52 | -0.033 / -0.037 / 0.11 | 0.026, 0.048 | 0.42 |
| 4.2 high | 0.115 / 0.185 / 1.0 | 0.93 / 0.85 / 0.49 | 0.033, 0.028 | 0.68 |
| 4.2 low | 0.032 / 0.037 / 0.29 | 0.007 / 0.15 / 0.08 | 0.005, 0.070 | 0.17 |

The 4.2 rows compare up to T* = 0.78 (high) and 0.95 (low), where the ANSYS
runs ended.

- **The two codes agree with each other far better than either agrees with
  the paper**, although the discretisations differ: NEM at the integration
  points with an explicit update, against nodal Galerkin with an implicit one.
- **The largest element phi differs.** It is 0.99 in his code and 0.29-0.85
  in Abaqus. This is the known difference of the seed: in his code the seed
  is a jump at the integration points, in Abaqus a nodal share.
- **What this shows.** With Table 2 and consumption g phi, the 2024 figures
  are not reproduced by this implementation, nor by his. The gap to the
  paper lies in the inputs or the consumption form of the 2024 runs, not in
  the implementation.

## Biofilm on an implant collar (6 Oct 2026)

`make_implant_inp.py`, `summarize_implant.py`: a 0.25 mm biofilm layer on a
4.1 mm implant's collar, 2 mm high (a sulcus), 90 degree sector, 6x24x20
C3D8T + nutrient UEL, bonded to rigid titanium, nutrient from the gingival
margin, phi = 1 in the ring on the titanium at t = 0; case 6. Geometry
values are mine, not from a paper.

- Consumption 6 (zero order, as the cube) on 2 mm depth drives c to -0.92 at
  the base: too strong for this depth. With consumption 1 (Table 2) c stays
  0.68 (base) to 0.98 (margin).
- Consumption 1 with species-weighted growth prop(36) = 1/3 (as the ANSYS
  week runs), T* = 1.1, on the titanium: share phi_1/(phi_1+phi_2) 0.227 at
  the base to 0.199 at the margin, alpha - 1 2.63e-4 to 2.53e-4 (the
  composition gives a 4 % depth gradient of growth), von Mises 1.5-1.6e-4 Pa,
  sigma_tt ~ sigma_zz ~ -1.4e-4 Pa (the layer is held by the titanium), and
  the shear on the bond sigma_rz peaks at the two ends of the layer
  (6.7e-5 Pa at the base, 6.5e-5 at the margin, ~40 % of von Mises), where
  detachment would start. phi spreads through the 0.25 mm layer by T* ~ 1
  (beta = 0.02), so phi itself is almost uniform; the depth dependence comes
  from the nutrient through the composition only.
- 2 CPUs, about 1.5 hours (the material server call per integration point
  dominates).
- Tooth with the front term (`--ri 4.0 --bulge 0.5 --nut outer --front 10`,
  consumption 1, prop(36) = 1/3, T* = 1.1, dt = 0.002; dt = 0.01 diverged, the
  front crossed about 1.2 elements per step), Fortran point model, 4 CPUs,
  7-8 min each: `tooth_fr_c1` (growth on every face, w = 0.5) and
  `tooth_pr_c1` (Eq. 34 as printed). Against `tooth_c1` (no front term),
  inner ring (on the enamel):

  | | no front term | front, w = 0.5 | front, as printed |
  |---|---|---|---|
  | phi | 0.17-0.20 | 1.01 | 0.73 (base) to 1.14 |
  | share phi_1/(phi_1+phi_2) | 0.17-0.21 | 0.023 | 0.023-0.031 |
  | alpha - 1 | 2.5e-4 | 1.0e-3 | 0.8-1.0e-3 |
  | von Mises | 1.4-2.3e-4 | 1.9-2.5e-2 | 1.5-2.4e-2 |

  The layer fills with biofilm in both front forms; above z = 0.25 mm they
  differ by less than 2 % in von Mises, at the base the printed form gives
  40 % less (phi 0.73 there). The stress rises about 100 times for 4 times the
  growth; I have not checked yet what makes up that factor, so these
  stresses are not to be used before that.

## Next: species carried in space (design, not implemented)

Now only the amount phi is a field. The composition lives in the point model at
each integration point and does not move: where phi spreads into a point, that
point's model starts from the start share (prop(35)), not from the composition
of the biofilm that arrived. With the same beta for every species, conservation
of each species gives for n species n - 1 more fields, e.g. for two
species the amount of species 1, phi_1 (phi_2 = phi - phi_1):

    d phi_1/dt = beta lap(phi_1) + R_1,

R_1 = the point model's rate of species 1. In Abaqus this is one more DOF (13)
on the nutrient UEL, with R_1 committed per increment as c is. The open
question is the point model's state: it would have to take phi_1 from the field
at the start of each increment instead of carrying it, which changes the
call-site fragments (now unchanged from ANSYS) and the stand-alone check.

Does it matter? The 16^3 two-species runs (case 6, start share 0.5), share
phi_1/(phi_1+phi_2) at T* = 1.1, Abaqus (ANSYS the same within 5 %):

| | seed | layer 1 | layer 2 |
|---|---|---|---|
| no nutrient | 0.055 | 0.34 | 0.50 |
| consumption 6 | 0.25 | 0.40 | 0.50 |

One element outside the seed the share is 0.34-0.40 instead of the seed's
0.05-0.25, two elements out it is the start value: the composition of the
spreading biofilm is set by the start share, not by the biofilm it came
from. The stress does not see this as long as the growth uses the total
(prop(36) = 0); with species-weighted growth it would.

## Note

Abaqus shortens the last increment by ~5e-16 to land on T* = 1; the judge now
takes the point model's step per increment from the trace (it assumed a constant
step), which is what made the stand-alone comparison exact.
`tests/test_abaqus_composition_umat.py` repeats the share/alpha/stress checks
with gfortran (no Abaqus needed).
