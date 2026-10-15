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

3D seed: `make_seed_inp.py [n] [outdir] [--lam 1 --psiref 0.32] [--antibiotic]` writes `uel_seed.inp`, the seed of
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

Not yet: the run inside Abaqus (Keio or IKMHIWI03; expect phi(0) = NT11 0.776, c(0) = NT12 0.936), finite-strain mechanics
(small strain now), a closed-form tangent (central differences now).
