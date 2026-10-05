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

Abaqus shortens the last increment by ~5e-16 to land on T* = 1; the judge now
takes the point model's step per increment from the trace (it assumed a constant
step), which is what made the stand-alone comparison exact.
`tests/test_abaqus_composition_umat.py` repeats the share/alpha/stress checks
with gfortran (no Abaqus needed).
