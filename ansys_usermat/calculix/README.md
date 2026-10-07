# Growth law in CalculiX: an independent solver check

ANSYS and Abaqus cannot run in a cloud session (licences, installers).
CalculiX 2.21 is open source and calls Abaqus-style UMATs, so the growth law
`F = Fe Fv Fg`, `Fg = alpha I` of `umat_biofilm_visco.f` can be run there
unchanged. This is a check of the material routine in a third solver, not
Abaqus work for the thesis.

## How

- `build_ccx.sh`: builds `ccx_2.21` with the UMAT linked in (adds
  `ABA_PARAM.INC`, system SPOOLES/ARPACK, gfortran `-fallow-argument-mismatch`).
- The material name starts with `ABAQUSNL`, so CalculiX passes the
  deformation gradient (`DFGRD1`) to the UMAT.
- The growth `alpha - 1` enters through the temperature (`TEMP + DTEMP`), ramped
  from 0 to 0.05 over 4 increments with NLGEOM.
- Material as `apdl/t_growth_constrained.dat`: C10 = 0.2e-3 MPa, C01 = 0,
  D1 = 5e3 1/MPa, eta = 0 (elastic), neo-Hooke.

## Result (2026-10-04)

| case | CalculiX | closed form | ANSYS (t_growth_constrained.dat) |
|---|---|---|---|
| `g_con.inp`, all nodes fixed | SX = SY = SZ = −1.019276e−04 MPa, shear 0, all 8 Gauss points | −1.019275856e−04 | −1.0193e−04 |
| `g_free.inp`, symmetry faces only | u = 0.05 on the far faces, stresses < 1e−19 | u = alpha − 1 = 0.05, σ = 0 | – |

Run: `ccx_2.21 g_con` and read `g_con.dat`.

## Seed stress against mesh and element type (2026-10-04)

`seed_mesh_study.py` solves the problem of `apdl/mesh_study_seed.py` (the 32-element
seed of the partner's deck, growth alpha − 1 = 1.1e-3 as an eigenstrain, E (phi² + f),
E = 10 Pa, nu = 0.49, f = 1e-3, three corner nodes) in CalculiX, linear elastic.
Seed von Mises mean / seed mean stress, Pa:

| element | 8³ | 16³ | 32³ |
|---|---|---|---|
| Python hex8, selective reduced (B-bar-like) | 5.70e-5 / −3.77e-5 | 3.59e-5 / −2.52e-5 | 3.26e-5 / −2.26e-5 |
| C3D8, full 2×2×2 | 5.96e-5 / −1.41e-4 | 3.41e-5 / −5.69e-5 | 3.14e-5 / −3.34e-5 |
| C3D8I, incompatible modes | 6.09e-5 / −5.51e-5 | 3.69e-5 / −2.76e-5 | 3.27e-5 / −2.31e-5 |
| C3D20R, quadratic | 3.11e-5 / −4.43e-5 | 2.91e-5 / −2.76e-5 | 3.05e-5 / −2.29e-5 |

- The Python solver with full integration (`solve(..., full=True)`) gives the
  C3D8 numbers to four digits at 8³ and 16³ (5.9627e-5 / −1.4096e-4 at 8³): the
  Python solver is verified by an independent code.
- At 32³ the three element types without locking agree: seed von Mises mean
  3.05–3.27e-5 Pa, seed mean stress −2.26 to −2.31e-5 Pa.
- The 8³ mesh with a locking-free element is 1.75× (von Mises) and 1.6× (mean
  stress) above these values, as found in Python.
- Full integration at 8³ overstates the mean stress six-fold (−1.41e-4 against
  −2.3e-5): at nu = 0.49 the 8³ mesh locks without B-bar. So whether the ANSYS
  runs use B-bar (SOLID185 KEYOPT(2) = 0, the ETLIST check) decides between a
  factor 1.6 and a factor 6 on the pressure.
- Maxima grow with refinement for the linear elements (re-entrant corners of the
  staircase seed); quote averages, not peaks.

## Does step 2 reach the stress? (2026-10-04)

`two_way_stress.py` on the fields of `two_way_step2.py --dt 0.025 --d 0.3333333333 --save`
(Klempt 2024 test case 4.1, alpha − 1 = k_alpha ∫phi dt up to 1.5e-3, C3D8I on the
20³ grid, E (phi² + f), rigid-body constraints only). Biofilm = elements with phi ≥ 0.5;
stresses in Pa:

| case, d | biofilm volume | von Mises mean | von Mises 99 % | mean stress, mean |
|---|---|---|---|---|
| one-way (d = 0) | 0.735 | 3.43e-3 | 6.21e-3 | −1.72e-4 |
| case 3, d = 1/3 | 0.768 (+4.5 %) | 3.57e-3 (+4 %) | 6.34e-3 (+2 %) | −1.46e-4 |
| case 6, d = 1/3 | 0.657 (−11 %) | 3.01e-3 (−12 %) | 5.70e-3 (−8 %) | −2.43e-4 |

The change in spreading carries over to the stress at about the same relative
size: the composition reaches the stress through the amount of biofilm, without
any species-dependent stiffness. Recorded as a cloud check; not in the thesis
(decided 2026-10-04).
