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

Abaqus source: `make_abaqus_source.py` writes `unified_uel_abaqus.for`, the fixed-form point model and the
free-form element in one file, joined by `!DEC$ FREEFORM` (Intel Fortran). Instructions for IKMHIWI03:
`ansys_usermat/apdl/CLOUD_TO_IKMHIWI03_2026-10-14_C.ja.md` A-11.

Abaqus cannot print a user element's state variables, so the element keeps the latest Gauss-point record (x, y,
z, phi, c, alpha - 1, living biomass, von Mises, pressure) and `UEXTERNALDB` appends it to `uel_gp.txt` in the job
directory every 10th increment and at the end of the analysis (it also reads the point model's constants at
LOP = 0, before any element call). Through the emulation (which now passes each element its own JELEM) the file
equals the state variables: 8^3, vM mean 0.0130, max 0.0164 over the elements with phi > 0.5.

Not yet: the run inside Abaqus (expect phi(0) = NT11 0.776, c(0) = NT12 0.936 for the column).
