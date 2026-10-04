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
