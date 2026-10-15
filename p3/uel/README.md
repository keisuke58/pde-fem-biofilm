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

The growth agrees to 1e-6; the stress difference halves with the mesh, as expected from the element
formulation (mean dilatation against B-bar; in ANSYS alone, enhanced strain against B-bar differs by 12.7 % at
8^3 and 2.5 % at 16^3). The largest value sits at mirror-image elements in the two codes (equal values).

## Rigid substrate (A-18, 15 Oct 2026, ANSYS)

The frozen seed of A-16 (8^3) on a plate 0.25 thick with shared interface nodes (81), plate bottom fixed,
compared with u = 0 on the base (`results/2026-10-uel/rigid_plate_8.json`). Relative difference of the von Mises
stress (seed / seed interior / layer above the interface): plate 10^9 times stiffer than the biofilm +4e-12 /
+8e-12 / +4e-12; 10^3 times -0.045 / -0.058 / -0.045 %; 10 times -4.2 / -5.4 / -4.3 %. alpha - 1 is identical
(phi frozen). A tooth or implant (GPa, over 10^8 times the biofilm) can therefore be replaced by u = 0 on the
surface the biofilm sits on, for the growth-induced stress; external loads such as occlusion need their own check.
