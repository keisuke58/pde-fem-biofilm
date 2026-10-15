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

Result: second call = first call (RHS, AMATRX, SVARS identical), RHS + R = 0, all 8 Gauss points started.

Not yet: a run inside Abaqus (input deck with `*USER ELEMENT`; Keio or IKMHIWI03), finite-strain mechanics
(small strain now), a closed-form tangent (central differences now).
