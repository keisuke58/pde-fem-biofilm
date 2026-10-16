# The unified model as an Abaqus user element (step 3)

`unified_uel.f90` is the two-level model of `p3/unified_model.tex` as an Abaqus UEL: 8-node hexahedron,
dofs per node u (1-3), phi (11), c (12); composition and growth at the 8 Gauss points (the point model
of `ansys_usermat/coupling/ecology_native.f`, damped Newton). Plan: `P3_STEP3_FEM_PLAN.ja.md`.

`driver.py` tests it without Abaqus: it calls the element through C entry points and assembles a
column of elements with Newton's method, as Abaqus would.

```
gfortran -O2 -fPIC -c ../../ansys_usermat/coupling/ecology_native.f -o ecology_native.o
gfortran -O2 -fPIC -ffree-line-length-none -c unified_uel.f90 -o unified_uel.o
gfortran -shared -o libunified_uel.so unified_uel.o ecology_native.o
python driver.py
```

Results (15 Oct 2026):
- one element, uniform: Gauss-point composition = converged point model, 7.8e-16;
- column along x, seed at x <= 0.25, nutrient supplied at x = 1 (20 and 40 elements): alpha - 1 near the
  substrate 0.1371 / 0.1370, phi(0) 0.776 / 0.772, c(0) 0.936 / 0.938, living 0.526; the 1D finite-difference
  model `JAXFEM/unified_two_level_fuel_1d.py` gives 0.137, 0.77, 0.938, 0.53. Newton residual < 1e-10.

The Abaqus entry point `SUBROUTINE UEL` is at the end of `unified_uel.f90`. On every call it advances the
composition from the start-of-increment SVARS (so repeated Newton iterations of one increment give the
same result), returns RHS = -R and, for LFLAGS(3) = 1 or 2, AMATRX = dR/dU; on a point-model failure it
sets PNEWDT = 0.5. NSVARS = 8 x 14, NPROPS = 16 (order in the header of `unified_uel.f90`).
`test_uel_wrapper.f90` calls it twice from the same state:

```
gfortran -O2 -ffree-line-length-none test_uel_wrapper.f90 unified_uel.o ecology_native.o -o test_uel_wrapper
./test_uel_wrapper
```

Result: second call = first call (RHS, AMATRX, SVARS identical), RHS + R = 0, all 8 Gauss points started;
the closed-form tangent (`tangent`) equals central differences of the residual (`tangent_fd`) to 8.7e-10
relative, with a displacement field and growth active. Same results as before in `driver.py` and the
emulation; 8^3 seed 70 s -> 53 s.

Abaqus input: `make_column_inp.py [n] [outdir]` writes `uel_column.inp` (check 2 as an Abaqus job: a 1e-8
step that prescribes phi and c at every node, the only way to give these dofs initial values, then 60
increments) and `uel_eco.txt` (the point model, read through `BIOFILM_ECO_CASE`). The build and run
commands are in its header. `emulate_abaqus.py` runs the same two steps by calling `UEL` itself with
Abaqus's argument list and the point model from that file: 20 elements give alpha - 1 0.1371,
phi(0) 0.776, c(0) 0.936, living 0.526, the numbers of `driver.py` (15 Oct 2026).

3D seed: `make_seed_inp.py [n] [outdir] [--lam 1 --psiref 0.32] [--antibiotic] [--small]` writes the seed deck, the seed of
`JAXFEM/fuel_front_3d_coupled.py` (unit cube, hemisphere of radius 0.25 on z = 0, nutrient at x = 1). Optional
PROPS(17) = total time at which the antibiotic starts and PROPS(18) = alpha* (the element switches
`eco_alstar` from TIME(2)). `emulate_abaqus.py --seed n [--antibiotic]` runs it (8^3: 1 min, 16^3: 45 min, most of it
the dense-fill sparse direct solve in Python). Against `JAXFEM/fuel_front_3d_coupled.py` (lam = 0, 16^3, finite differences,
explicit splitting), 15 Oct 2026:

| 16^3 | UEL (emulated) | Python 3D |
|---|---|---|
| phi at the centre | 0.778 | 0.779 |
| alpha - 1, nutrient side | 0.080 (element) | 0.085 (node) |
| edge 30 deg, nutrient / far side | 0.0135 / 0.0036 | 0.0126 / 0.0034 |
| same, antibiotic from T* = 0.5 | 0.0104 / 0.0027 (-23 %) | 0.0096 / 0.0025 (-24 %) |

The two discretisations agree within 8 %, and the antibiotic effect within one percentage point.

3D seed decks are named `uel_seed<n>[_small][_ab].inp` (the Abaqus job name); `--small` writes JPROPS(1) = 0.

Finite strain (default since 15 Oct 2026): F = F_e F_g, F_g = (1 + xi) I, total Lagrangian, energy
(1+xi)^3 (phi^2 + f) [mu/2 (J_e^(-2/3) tr C_e - 3) + K/2 (J_e - 1)^2] (Klempt et al. 2024 Eq. 20 with a volumetric
penalty in place of the incompressibility constraint, as the UMAT core); `nh_stress`. Small strain with JPROPS(1) = 0
(`--small` in the deck generators and the emulation), as the Python 3D model. Checks (`test_uel_wrapper`): tangent
against central differences, mechanics rows, 1.0e-11 (small displacement) and 7.2e-11 (displacement gradients
about 0.3); free growth u = xi X gives a mechanics residual of 2e-21. 8^3 seed: edge 0.0141 (finite) against
0.0143 (small strain), alpha - 1 about 0.07.

Near incompressibility (nu = 0.49; Klempt et al. 2024 use an incompressible material) the fully integrated
hexahedron locks. JPROPS(2) = 1 switches on the mean dilatation (Q1/P0, Simo, Taylor and Pister 1985): the
volumetric energy sees the element mean of J_e, one pressure per element, as the element pressure that Klempt
et al. 2024 condense. Its residual equals the gradient of the element energy computed independently in Python
(3.9e-11); its tangent is by central differences for now. `bilayer_test.py`: strip clamped at one end, top half
grown, tip deflection against beam theory (0.75 xi L^2 / h):

| nu | mesh | fully integrated | mean dilatation |
|---|---|---|---|
| 0.3 | 40 x 4 | 0.994 | 1.048 |
| 0.3 | 80 x 8 | 1.020 | 1.034 |
| 0.49 | 40 x 4 | **0.674** | 1.066 |
| 0.49 | 80 x 8 | 0.910 | 1.051 |

Both converge to about 1.04; with nu = 0.49 the fully integrated element is a third too stiff on the coarse
mesh. Started from u = 0 with the full growth, a nearly incompressible grown layer is in hydrostatic
compression and the tangent is indefinite (geometric stiffness), so the test raises the growth in steps, as a
run does. F-bar (de Souza Neto et al. 1996) was tried first and dropped: same residual check, but erratic
Newton convergence in this test.

Stress output: with VARIABLES = 128 the element writes von Mises stress and pressure (Cauchy, pressure positive
in compression) per Gauss point to SVARS(113:128). Check: fully constrained uniform growth gives vM = 0 and
p = 3 K xi phi^2 (small strain exactly, finite strain 2 % lower at xi = 0.01). The mean-dilatation tangent is now
in closed form (against central differences 5.8e-10, 9.5e-11 at displacement gradients about 0.3).

3D seed with the stresses (`emulate_abaqus.py --seed n [--nu 0.49] [--md]`, finite strain, E = 1; biofilm =
elements with phi > 0.5):

| case | edge 30 deg, 8^3 | vM mean, 8^3 | edge 30 deg, 16^3 | vM mean, 16^3 |
|---|---|---|---|---|
| nu = 0.3 | 0.0141 | 0.0130 | 0.0132 | 0.0144 |
| nu = 0.49, fully integrated | 0.0154 | **0.0232** | 0.0150 | 0.0201 |
| nu = 0.49, mean dilatation | 0.0184 | 0.0134 | 0.0156 | 0.0163 |

With nu = 0.49 the fully integrated element puts spurious pressure into the stress: the von Mises mean is
73 % (8^3) and 23 % (16^3) above the mean dilatation. The mean-dilatation vM rises with refinement as the edge
is resolved, as the thesis's 8^3 to 24^3 study of the element does. For the Keio Abaqus runs with a UMAT at
nu = 0.49 this means hybrid elements (C3D8H), not C3D8.

Large growth (`--L0 10`, alpha - 1 = 0.28 on the nutrient side, 8^3, nu = 0.3): small strain overstates the edge
displacement by 10 % (0.0560 against 0.0507 finite strain) and the von Mises mean by 3 % (0.0413 against 0.0402).
At the default growth (alpha - 1 = 0.07) the difference is 1.4 %.

Nutrient in the current configuration (JPROPS(3) = 1, finite strain): amount J c per reference volume and the
pulled-back flux J D C^-1 Grad c. `current_config_test.py`: a uniformly expanded element (J = 1.331, no reaction,
no flux) keeps the amount, c = 0.7513 = 1/1.331; a bar stretched to F = diag(2, 1, 1) carries half the reference
flux, as D J / F11^2 says. Tangent by differences. 8^3 seed: edge 0.0136 against 0.0141 (reference configuration),
von Mises mean 0.0123 against 0.0130 (the grown biofilm dilutes the nutrient a little).

Composition acting on the phase (`--lam 1 --psiref 0.32`, 16^3, small strain, as the Python 3D model of
`JAXFEM/fuel_front_3d_coupled.py`):

| 16^3, lam = 1, Psi_ref = 0.32 | UEL (emulated) | Python 3D |
|---|---|---|
| biofilm volume (reference), no antibiotic | 0.0346 | 0.0351 |
| same, antibiotic from T* = 0.5 | 0.0304 (-12 %) | 0.0309 (-12 %) |
| phi at the centre, no antibiotic / antibiotic | 0.783 / 0.676 | 0.78 / 0.67 |
| edge 30 deg, no antibiotic / antibiotic | 0.0124 / 0.0090 | 0.0115 / 0.0083 |
| von Mises mean (phi > 0.5), no antibiotic / antibiotic | 0.0136 / 0.0082 (-39 %) | (not computed) |

The antibiotic shrinks the biofilm itself by the same 12 % in both discretisations.

Abaqus source: `make_abaqus_source.py` writes `unified_uel_abaqus.for`, all in fixed form: the point model as it
is and the element converted from free to fixed form (column 7, continuation in column 6, comments from column 1,
132 columns for Abaqus's /extend-source). The first version joined a free-form part with `!DEC$ FREEFORM`; Abaqus's
`ifort /fpp` read it as fixed form and broke the module (IKMHIWI03, A-11, 15 Oct 2026). Check without Abaqus:
`gfortran -cpp -ffixed-form -ffixed-line-length-132 -c unified_uel_abaqus.for`, then `test_uel_wrapper` and the
emulation with that object give the same results to the last digit (column 0.1371 / 0.776 / 0.936; 8^3 seed edge
0.0141, vM mean 0.012997). Instructions for IKMHIWI03: `ansys_usermat/apdl/CLOUD_TO_IKMHIWI03_2026-10-14_C.ja.md` A-11.

Abaqus cannot print a user element's state variables, so the element keeps the latest Gauss-point record (x, y,
z, phi, c, alpha - 1, living biomass, von Mises, pressure) and `UEXTERNALDB` appends it to `uel_gp.txt` in the job
directory (or at the absolute path in BIOFILM_UEL_GP) every 10th increment and at the end of the analysis (it also reads the point model's constants at
LOP = 0, before any element call). Through the emulation (which now passes each element its own JELEM) the file
equals the state variables: 8^3, vM mean 0.0130, max 0.0164 over the elements with phi > 0.5.
Until 16 Oct 2026 the record was taken before the growth update, so its alpha - 1 was one increment behind (vM and
pressure were not); now it is taken after `store_growth`. Implant sector 2x6x8 (`make_implant_uel_inp.py`): the file
agrees with the emulation to 4e-7 in alpha - 1 and 5e-7 in vM at the last increment (before: 2 % in alpha - 1).

Abaqus procedure: Abaqus does not allow boundary conditions on DOF 11 and 12 in *STATIC (IKMHIWI03, second try of
A-11), so the decks use a transient *COUPLED TEMPERATURE-DISPLACEMENT with fixed increments (no DELTMX) and declare the
element UNSYMM. The steady-state form chose its own increments and gave NT12 = 7.9e6 at x = 0. The transient form
needs a specific heat (third try): dummy C3D8T elements on the same nodes with stiffness, conductivity and specific
heat 1e-9 (`make_column_inp.overlay`). Linking needs the DLL
runtime (`run_abaqus.ps1 -LibsDll`) because of the OPEN / WRITE of the output file.

Abaqus runs (IKMHIWI03, 15 Oct 2026, fourth try): all five jobs pass. Column: NT11 = 0.77532, NT12 = 0.93565 at x = 0.
Seed 8^3 against the emulation, Gauss point by Gauss point: phi within 2.6e-7 and the living biomass identical, but c up
to 0.011 and alpha - 1 up to 0.009 lower (von Mises mean 0.01135 against 0.01300). Abaqus accepted most increments after
one iteration (70 iterations in 61 increments): its default check for temperature-type DOFs, 0.5 % of the average
flux, is too loose for phi and c. The decks now set `*CONTROLS, PARAMETERS=FIELD, FIELD=GLOBAL` to 1e-8 and allow more
iterations, and state AMPLITUDE=STEP on the growth step; the fifth try reruns the seeds.
Fifth try: same values, so the tolerance was not the cause. The cause was the start: Abaqus switched the supply at
x = 1 from 0.1 to 1 inside the first growth increment (start-of-increment value 0.1), the emulation had it at 1 from
the start. With the same switch in the emulation, Abaqus and the emulation agree Gauss point by Gauss point to
c 1.8e-6, alpha - 1 5e-6, von Mises 2e-6: the element works in Abaqus. The jump enters the consistent mass term and
lowered the growth by 12 % on 8^3; the seed deck now prescribes the supply already in the INIT step.
Sixth try (supply from the INIT step): all four seeds pass, and Abaqus gives the emulation's numbers (von Mises mean / max
over the elements with phi > 0.5):

| job | Abaqus | emulation |
|---|---|---|
| 8^3, finite strain | 0.012998 / 0.016366 | 0.012997 / 0.016365 |
| 8^3, small strain | 0.013120 / 0.016585 | 0.013119 / 0.016584 |
| 8^3, nu = 0.49, mean dilatation | 0.013419 / 0.015773 | 0.013419 / 0.015773 |

The unified model runs in Abaqus 2024 as a user element (results in `ansys_usermat/apdl/results/2026-10-uel/`).

## Toward real numbers (15 Oct 2026, `calibrate_growth.py`)

The growth rate is chosen so that the seed's mean alpha - 1 at Day 21 equals the measured CH volume increase,
V(21)/V(3) = 1.80, read as alpha^3 = 1.80 from Day 1 (alpha - 1 = 0.216; assumptions: all of the increase is growth,
little growth from Day 1 to 3). 8^3, finite strain, nu = 0.49, mean dilatation.
- With the yield Y = 1 the growth has a ceiling, alpha - 1 = 0.107 near L0 = 13: faster growth eats the nutrient and
  starves the seed (L0 = 40 gives 0.052). The measured increase is out of reach.
- With less nutrient per growth it is reached: Y = 0.25, L0 = 53.1 gives alpha - 1 = 0.2154 and von Mises / E = 0.0794;
  Y = 0.1, L0 = 113.9 gives 0.2165 and 0.0773. The stress does not depend on the choice of Y within 3 %.
- von Mises / (E (alpha - 1)) = 0.36-0.37 in the seed, about a fifth of a laterally fixed film's 1/(1 - nu).

Seed von Mises mean (max) for measured stiffness values:

| E | von Mises mean | max | pressure |
|---|---|---|---|
| 10 Pa (Klempt et al. 2024) | 0.77-0.79 Pa | 0.89-0.92 Pa | 0.36-0.39 Pa |
| 600 Pa | 46-48 Pa | 53-55 Pa | 22-23 Pa |
| 41 kPa (Pattem et al. 2018, upper end) | 3.2-3.3 kPa | 3.6-3.8 kPa | 1.5-1.6 kPa |

All far below the pressures reported to stop growth (0.1-1 MPa). The stiffness spans three decades and is the
largest uncertainty. Not yet: the tooth and implant geometry, detachment, a stiffness for the measured composition.

## Limit of Klempt et al. 2024 (15 Oct 2026, `klempt2024_limit_test.py`)

With one species, the living fraction frozen near 1, the nutrient fixed at c = 1, the phase fixed at
phi = 1 and lambda = 0, the growth reaction is Eq. 36 of Klempt et al. 2024, alpha_dot = k_alpha phi with
k_alpha = L0 liv ln 10. One element (one Gauss-point state, finite strain, E = 10 Pa, nu = 0.49, f = 1e-3),
compared with references written without the element's source:

| case | quantity | k_alpha = 1e-3, T = 1.1 (Table 2) | k_alpha = 0.5, T = 1 (alpha = 1.5) |
|---|---|---|---|
| free (3-2-1) | alpha - 1 vs Eq. 36 | 1.097802e-3 vs 1.097801e-3 | 0.4990009 vs 0.4990005 |
| free | edge stretch vs alpha, stress | equal to 8 digits, sigma 1e-12 | equal to 8 digits, sigma 4e-14 |
| all u = 0 | sigma vs (phi^2+f) K (alpha^-3 - 1) | rel. 3e-14 | rel. 2e-16 |
| u_x = 0 at both ends | lateral stretch, sigma_xx vs numpy minimisation of W | rel. 8e-9 | rel. 2e-10 |

The small difference in alpha - 1 (1e-6 relative) is the drift of the frozen point model over the run.
This checks the growth law and the finite-strain mechanics of Klempt et al. 2024 (Eq. 20, 36). The phase
equation of Klempt et al. 2024 (the k_alpha alpha term and the front term of Eq. 34) is not a limit of
this model, so their Section 4 test cases cannot be reproduced with the element; the Python reproduction
(`JAXFEM/`, Appendix D of the 19 Oct deck) covers them.

## Cross-code check with ANSYS: seed with the phase frozen (A-16, 15 Oct 2026)

Same problem in the UEL (`seed_frozen_reference.py`, mean dilatation) and in ANSYS (SOLID185 B-bar with the
unchanged core of the thesis USERMAT, phi as nodal temperature; `ansys_usermat/apdl/frozen_seed/`). The two
codes share no source. Seed = elements with mean nodal phi > 0.5.

| | 8^3 UEL | 8^3 ANSYS | diff. | 16^3 UEL | 16^3 ANSYS | diff. |
|---|---|---|---|---|---|---|
| seed alpha - 1 | 6.76058e-4 | 6.76057e-4 | 1.2e-6 | 8.81599e-4 | 8.81598e-4 | 8e-7 |
| seed von Mises | 3.8496e-3 | 3.6559e-3 | -5.0 % | 5.4222e-3 | 5.2881e-3 | -2.5 % |
| largest von Mises | 7.2210e-3 | 7.0229e-3 | -2.7 % | 1.18150e-2 | 1.17799e-2 | -0.3 % |

At 24^3 (UEL run on fifa, `results/seed_frozen_24.json`): seed alpha - 1 8.70129e-4 against 8.70128e-4 (8e-7),
seed von Mises 5.0930e-3 against 5.0100e-3 (-1.6 %), largest 1.34523e-2 against 1.34248e-2 (-0.2 %). In the seed
interior (elements within 0.2 of (0.5, 0.5, 0)) the difference is -2.7 % (8^3), -1.2 % (16^3), -0.2 % (24^3):
6.5492e-3 (UEL) against 6.5360e-3 (ANSYS), and ANSYS at 32^3 gives 6.510e-3.

The growth agrees to 1e-6; the stress difference halves with the mesh, as expected from the element
formulation (mean dilatation against B-bar; in ANSYS alone, enhanced strain against B-bar differs by 12.7 % at
8^3 and 2.5 % at 16^3). The largest value sits at mirror-image elements in the two codes (equal values).

## The composition of a newly colonised Gauss point (WP4, 16 Oct 2026, `init_rule_study.py`)

When phi of the field first exceeds phi_min at a Gauss point, the point model there has no composition.
Until now the element started it from the file's `phi_init` line (the Day-1 state of the calibration, psi = 0.999,
gamma = 0) wherever and whenever that happens, the same state for the seed and for every point the biofilm
spreads to; the ANSYS USERMAT (`INIT_ECO_IF_ZERO`) does the same with a fixed default composition. Felix Klempt's
answer (KEIO_PLAN.ja.md §1) was that this is outside the model, and suggested the neighbour's composition or a
neighbourhood mean. `unified_uel.f90` now has `init_rule` (Abaqus JPROPS(4), `make_seed_inp.py --init-rule k`;
emulation `uel_set_init_rule`): 0 the default (unchanged: column 20, seed 8^3 and implant 2x6x8 reproduce the
earlier runs and `results/implant_2x6x8.json` bit for bit), 1 the composition of the nearest colonised Gauss point,
2 the phi-weighted mean of the colonised Gauss points within 1.5 element sizes (PROPS(19) sets the radius), 3 the
point model's steady state under the local nutrient. A rule applies only once a colonised point exists, so the seed
itself keeps the Day-1 state; the state variable 14 records the start (1 default, 1 + rule). Rules 1 and 2 need the
other elements, which a UEL call does not see: every element writes its Gauss points (position, phi, composition,
started) into a module-level registry at every call, two slots per element keyed by the start-of-increment total
time, and a query reads each element's newest slot older than the current increment, i.e. the state at the end of
the previous increment, written at that increment's last (converged) call. Abaqus guarantees no call order among
elements, only that all are called every iteration, so nothing read is written in the same increment and the
result does not depend on the order; a cutback overwrites the slot. With threads the element calls are safe
(distinct elements, distinct memory, reads of the other slot), the growth of the registry is not: set
`BIOFILM_UEL_NEL` so UEXTERNALDB allocates it at LOP = 0; MPI domains would lose the neighbours across a domain
boundary, so one domain. The decks here run with cpus=1 anyway. `tests/test_uel_init_rule.py` checks the default
against the pre-WP4 build (1e-12) and each rule against the same start advanced by `uel_advance` with the
neighbours found in Python (1e-14).

Study with the point model seeing the local nutrient (c_ref = 1; with c_ref = 0 every point sees c* and a colonised
point's composition depends only on its age, so copying a neighbour copies a clock and rules 0-2 differ by 1e-5).
Column of driver.py check 2 (mechanics fixed) and the seed of `make_seed_inp.py` (finite strain, nu = 0.3, E = 1).
"Late" = the Gauss points colonised after the first growth increment (20 / 44 / 88 of 68-72 / 128-140 / 252-280 in
the column, 184 of 508 at 8^3 and 748 of 1612 at 12^3), the same set for every rule; living = sum phi_i psi_i there.

| column | rule | front x | alpha - 1 at the substrate | living, late points | share of species 1, late | alpha - 1, late |
|---|---|---|---|---|---|---|
| 20 | 0 default | 0.2351 | 0.1353 | 0.495 | 0.541 | 0.0149 |
| 20 | 1 nearest | 0.2359 | 0.1353 | 0.513 | 0.377 | 0.0153 |
| 20 | 2 mean | 0.2359 | 0.1353 | 0.513 | 0.377 | 0.0153 |
| 20 | 3 steady state | 0.2525 | 0.1352 | 0.891 | 0.031 | 0.0444 |
| 40 | 0 / 1 / 2 / 3 | 0.2240 / 0.2251 / 0.2251 / 0.2597 | 0.1353 | 0.433 / 0.513 / 0.513 / 0.891 | 0.472 / 0.377 / 0.377 / 0.031 | 0.0169 / 0.0173 / 0.0173 / 0.0521 |
| 80 | 0 / 1 / 2 / 3 | 0.2178 / 0.2190 / 0.2190 / 0.2767 | 0.1354 | 0.458 / 0.513 / 0.513 / 0.891 | 0.493 / 0.377 / 0.377 / 0.031 | 0.0191 / 0.0196 / 0.0196 / 0.0594 |

| seed | rule | alpha - 1 (phi > 0.5) | alpha - 1, all colonised | vM/E mean / max (phi > 0.5) | edge 30 deg | living, late | share 1, late | alpha - 1, late |
|---|---|---|---|---|---|---|---|---|
| 8^3 | 0 | 0.04155 | 0.01803 | 0.01511 / 0.01906 | 0.01627 | 0.557 | 0.667 | 0.0065 |
| 8^3 | 1 | 0.04157 | 0.01772 | 0.01511 / 0.01906 | 0.01627 | 0.479 | 0.633 | 0.0057 |
| 8^3 | 2 | 0.04157 | 0.01774 | 0.01511 / 0.01906 | 0.01627 | 0.480 | 0.634 | 0.0057 |
| 8^3 | 3 | 0.04152 | 0.01896 | 0.01511 / 0.01907 | 0.01631 | 0.765 | 0.060 | 0.0092 |
| 12^3 | 0 | 0.04168 | 0.01837 | 0.01718 / 0.03057 | 0.01446 | 0.556 | 0.666 | 0.0071 |
| 12^3 | 1 | 0.04170 | 0.01803 | 0.01719 / 0.03057 | 0.01447 | 0.481 | 0.634 | 0.0064 |
| 12^3 | 2 | 0.04170 | 0.01804 | 0.01719 / 0.03057 | 0.01447 | 0.482 | 0.635 | 0.0064 |
| 12^3 | 3 | 0.04162 | 0.01963 | 0.01717 / 0.03056 | 0.01445 | 0.761 | 0.060 | 0.0100 |

What the numbers say. The stress and the growth of the seed do not depend on the rule (vM/E to four digits, the
seed's alpha - 1 to 0.1 %): the newly colonised points carry little biomass and grow little within T, and the seed
keeps its own state under every rule. The rule shows in the composition and growth of the new points only: nearest
and mean are the same to 1e-3 (the colonised neighbours of a new point differ little among themselves), give the
composition of the front a value that does not change with the mesh (column 0.513 / 0.377 at 20, 40 and 80
elements; seed 0.479 / 0.481 at 8^3 and 12^3), and change their growth by -2 to +3 %. The default's late composition
moves with the mesh (column living 0.495 / 0.433 / 0.458) because the age of the Day-1 state at a point then depends on
when the mesh lets phi_min be crossed. The steady state (rule 3) is a different biofilm at the front (living 0.89 and
species 1 at 3 % instead of 38-54 %, three times the growth of the new points) and it moves the column's front by 7 to
27 % with the mesh, because a point that starts at the attractor skips the transient the seed went through.
Recommendation: rule 1 (nearest colonised Gauss point) as the physical choice for runs where the biofilm spreads, rule 2
when a smoother start is wanted (same result here), rule 3 only as a bound. The thesis runs keep the default (the
seed-dominated quantities are the same to four digits, so nothing there changes).

## Rigid substrate (A-18, 15 Oct 2026, ANSYS)

The frozen seed of A-16 (8^3) on a plate 0.25 thick with shared interface nodes (81), plate bottom fixed,
compared with u = 0 on the base (`results/2026-10-uel/rigid_plate_8.json`). Relative difference of the von Mises
stress (seed / seed interior / layer above the interface): plate 10^9 times stiffer than the biofilm +4e-12 /
+8e-12 / +4e-12; 10^3 times -0.045 / -0.058 / -0.045 %; 10 times -4.2 / -5.4 / -4.3 %. alpha - 1 is identical
(phi frozen). A tooth or implant (GPa, over 10^8 times the biofilm) can therefore be replaced by u = 0 on the
surface the biofilm sits on, for the growth-induced stress; external loads such as occlusion need their own check.

## Implant sector: mesh convergence (16 Oct 2026, cloud, Python emulation)

`make_implant_uel_inp.py --run` at nr x nt x nz = f (2 x 6 x 8), f = 1, 2, 3 (4 on fifa); `implant_convergence.py
[--bulge 1] [--fixed-seed]` tabulates `results/implant_*.json`: seed means, the interior without the edge rows, the
seed core (a fixed region: seed elements with 0.25 <= z <= 1.75 and 15 <= theta <= 75 degrees), the largest vM/E and
the observed order. Times are with four runs sharing four cores. The .json now also holds the interior means, the
largest value without the edge rows and the wall time; `compare_implant_gp.py <uel_gp.txt> <json>` compares an
Abaqus run element by element (exit 1 above `--tol`).

Seed = the inner element layer (as `make_implant_inp.py`), bulge 0:

| level | elements | time | seed alpha - 1 | seed vM/E | core alpha - 1 | core vM/E | max vM/E, at (r, theta, z) |
|---|---|---|---|---|---|---|---|
| 2x6x8 | 96 | 12 s | 0.04634 | 0.04775 | 0.01795 | 0.03544 | 0.2728 (2.094, 67.5, 1.875) |
| 4x12x16 | 768 | 160 s | 0.03523 | 0.02358 | 0.01180 | 0.01349 | 0.2375 (2.077, 71.2, 1.938) |
| 6x18x24 | 2592 | 25 min | 0.02962 | 0.01482 | 0.00960 | 0.00729 | 0.2004 (2.069, 77.5, 1.958) |

With bulge 1 at 2x6x8 / 4x12x16: seed alpha - 1 0.04427 / 0.03302, seed vM/E 0.04620 / 0.02383, largest vM/E 0.2927 /
0.2912 at (3.08, 67.5, 1.875) / (3.07, 71.3, 1.938). This seed thins with the mesh (0.125, 0.0625, 0.042 mm), so the
problem changes with the level and nothing converges (seed vM/E -51 %, -37 %). For a convergence study the seed must be
fixed: `--seed-layers f` seeds the inner f layers, 0.125 mm at every level (`_sl<f>` in the file name). With that seed,
bulge 0:

| level | elements | time | seed alpha - 1 | seed vM/E | core alpha - 1 | core vM/E | max vM/E, at (r, theta, z) |
|---|---|---|---|---|---|---|---|
| 2x6x8 | 96 | 12 s | 0.04634 | 0.04775 | 0.01795 | 0.03544 | 0.2728 (2.094, 67.5, 1.875) |
| 4x12x16 | 768 | 201 s | 0.04481 | 0.06206 | 0.01467 | 0.03637 | 0.7791 (2.077, 63.7, 1.938) |
| 6x18x24 | 2592 | 33 min | 0.04440 | 0.06125 | 0.01432 | 0.03606 | 1.2245 (2.069, 62.5, 1.958) |

Bulge 1, same seed (the third level with one cutback, see below; 36 min):

| level | seed alpha - 1 | seed vM/E | core alpha - 1 | core vM/E | max vM/E, at (r, theta, z) |
|---|---|---|---|---|---|
| 2x6x8 | 0.04427 | 0.04620 | 0.01490 | 0.03253 | 0.2927 (3.076, 67.5, 1.875) |
| 4x12x16 | 0.04278 | 0.06100 | 0.01206 | 0.03206 | 0.9808 (3.072, 63.7, 1.938) |
| 6x18x24 | 0.04173 | 0.05913 | 0.01166 | 0.03113 | 1.5766 (3.067, 62.5, 1.958) |

Changes from the level before and the observed order p (three levels, q = q_inf + C h^p; "-" where the differences
do not shrink with one sign):

| quantity | bulge 0: 1 -> 2, 2 -> 3, p (q_inf) | bulge 1: 1 -> 2, 2 -> 3, p (q_inf) |
|---|---|---|
| core alpha - 1 | -18.3 %, -2.4 %, 2.9 (0.01416) | -19.1 %, -3.3 %, 2.5 (0.01142) |
| core vM/E | +2.6 %, -0.9 %, - | -1.4 %, -2.9 %, - |
| seed alpha - 1 | -3.3 %, -0.9 %, 1.4 (0.04384) | -3.4 %, -2.5 %, - |
| seed vM/E | +30 %, -1.3 %, - | +32 %, -3.1 %, - |
| max vM/E (all) | +186 %, +57 %, - | +235 %, +61 %, - |

The core, a fixed region away from the edges, settles: alpha - 1 to 2-3 % at the third level with an order between 2
and 3, vM/E within 1-3 % at every level (its changes are below the level where an order can be read). Bulge 1 lowers
the core by 17-19 % in alpha - 1 and 8-14 % in vM/E (the wider collar at the top takes more nutrient). The seed means
include the top row and the sector cuts, so they move more. The largest vM/E does not converge: it sits in the top row
of elements next to the fixed inner face, where the seed (phi = 0.9) meets the empty layer (phi = 1e-4) and the
displacement is fixed on one side and free on the other, and grows 0.27 -> 0.78 -> 1.22 (bulge 0) and 0.29 -> 0.98 ->
1.58 (bulge 1) as the element reaches closer to the corner, as a singular corner does. It stays at theta = 62-68
degrees, in the inner element column, half an element below the top face. It is reported but left out of the compared quantities, as
the cube's corners are. The 8x24x32 level (sl4) is for fifa (`CLOUD_TO_FIFA.ja.md`, 16 Oct evening and night).

Newton: at 6x18x24 with the 0.125 mm seed the full Newton step of the first growth increment diverged in the mechanics
(phi and c already converged, a 0.235 mm displacement step on 0.042 mm elements, then NaN from an inverted element; the
element itself returned PNEWDT = 1). `emulate_abaqus.increment` now halves the step when the residual grows by more than
half or is not finite (up to 6 times; Abaqus: `*CONTROLS, PARAMETERS=LINE SEARCH`): bulge 0 then needs five halvings
in increment 1 and plain Newton after. The converged state is unchanged (`tests/test_uel_implant.py`, 2x6x8 to 1e-6).
With bulge 1 the line search was not enough: 6x18x24 stopped after 30 iterations in increment 1 (the residual crept
from 2.6e-2 to 3.5e-3 with steps of 1/8 to 1/16).

Cause (16 Oct, evening; first growth increment of 6x18x24 sl3 with bulge 1, own Newton loop with the residual per dof
kind and the corner Jacobians of the deformed elements): the geometry is sound before the step (smallest corner det J
7.8e-5 > 0, the sheared elements at the bulge included) and the residual is ordinary (u 3.8e-3, phi 5.3e-3, c 2.6e-2),
but the first Newton update is 0.97 mm at the top face in the empty layer (phi = 1e-4, r = 3.22), which inverts 94
elements at once. The update depends on the increment length: 6.3e-5, 2.5e-4, 1.0e-3, 2.4e-3 and 5.8e-3 mm at dt/64,
dt/16, dt/4, dt/2 and 3 dt/4 (linear in dt, as a response to the growth of one increment should be), 0.97 mm at dt and
2.1e-2 mm at 2 dt. So the coupled tangent of the start state (alpha = 1, u = 0, the empty layer at its floor stiffness
f = 1e-3) is nearly singular at the base increment dt = T/60: an eigenvalue of the 1/dt-reaction part crosses zero near
that dt and the near-null mode sits in the soft layer. It is not an inverted element and not the bulge geometry; the
bulge only moves the crossing onto the base increment (bulge 0 at the same level: 2.7e-2 mm, five halvings). At dt/4
the increment converges quadratically in three iterations, 4x12x16 sl2 in five at the full dt.

Remedy, as Abaqus would do it: automatic time incrementation (`emulate_abaqus.march`, used by
`make_implant_uel_inp.run`). dt starts at T/60, which is also the largest increment; an increment whose Newton fails
(residual not finite, 30 iterations, or PNEWDT < 1 from the element) is retried from the same start state with dt/4,
at most five times in one increment; after two consecutive increments converged in at most four iterations dt grows by
1.5 up to T/60; the last increment ends at T. Without a failure the increments are the 60 fixed ones, so converged
runs are unchanged (2x6x8, bulge 0 and 1: every per-element value identical, scalars to 2e-15). The .json records
`incrementation` (increments, cutbacks, iterations, smallest dt and the increment log). The deck's GROW step carries
the same settings: `*COUPLED TEMPERATURE-DISPLACEMENT, DELTMX=10` (never reached, phi and c are in [0, 1]; it switches
Abaqus to automatic incrementation), data line `dt0 = T/60, T, DTMIN = dt0/4^5, DTMAX = dt0`, and
`*CONTROLS, PARAMETERS=TIME INCREMENTATION` with I_C = 30 iterations and I_A = 5 cutbacks per increment; the cutback
factor 0.25 and the growth factor 1.5 are Abaqus's defaults. `tests/test_uel_implant.py` forces a failure in
increment 3 of a 2x3x4 sector and checks the retry, the recovery of dt and the end time.

With it, 6x18x24 sl3 with bulge 1 converges: the first attempt at dt = T/60 uses its 30 iterations, the retry at dt/4
converges in three, dt is back at T/60 from increment 8 on, 64 increments and 253 iterations in all (36 min, three
other runs on the four cores). Bulge 0 at the same level ran before the cutback existed (five line-search halvings in
increment 1, no cutback), so its .json has no `incrementation` record.
