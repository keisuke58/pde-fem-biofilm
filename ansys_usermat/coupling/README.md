# Gauss-point material bridge — skeleton

Scaffolding for the thesis' core deliverable: calling the calibrated **Python**
material model from the Fortran **UMAT/USERMAT** at each Gauss point, with the
verified inline Fortran core kept as the reference/fallback. This lets the
IKM-verified constitutive law and the Keio-side computational-mechanics work
meet at a single, well-defined interface.

> **Status: skeleton.** The Python side runs and is tested end-to-end here
> (no Abaqus/ANSYS needed); the Fortran side is a syntax-checked stub with the
> integration marked. Wiring it into a live solver is the continuation work.

## Interface contract

One Gauss-point evaluation, per increment:

| direction | quantity | shape | note |
|---|---|---|---|
| in | `F` (deformation gradient) | 3×3 | row-major |
| in | `Fv` (prior viscous state) | 3×3 | row-major |
| in | `alpha, C10, C01, D1, eta, mtype, dt` | scalars | growth driver + material props + increment |
| out | `stress` (Cauchy) | 6 | Voigt **11,22,33,12,13,23** (Abaqus order) |
| out | `Fv_new` | 3×3 | updated viscous state |
| out | `dsdePl` (material Jacobian) | 6×6 | ∂σ/∂ε, F-perturbation |

## Files

| file | role |
|---|---|
| `material_server.py` | Python side — NumPy reference core (mirrors the verified Fortran `BIOFILM_STRESS_CORE`) + F-perturbation tangent + a socket server. In production, swap the core for the JAX model (`material_models.py` / `JAXFEM/`) behind the same interface. |
| `protocol.py` | wire schema (newline-delimited JSON; one request→one response). Swap for a binary frame later without touching the physics. |
| `usermat_py_hook.f` | Fortran side — `ISO_C_BINDING` interface to a C shim `biofilm_py_eval` + array marshalling; the `PYTHON MATERIAL HOOK` call site. Falls back to the inline core on failure. |
| `../../tests/test_coupling.py` | end-to-end test: marshals a state through the protocol to the server and checks the response equals the in-process evaluation (runs in CI). |

## Two integration mechanisms

1. **ISO_C_BINDING (in-process).** A small C shim embeds/loads Python (or links a
   compiled model) and is called directly from `biofilm_py_hook`. Lowest latency;
   best once the model is stable.
2. **Local socket (out-of-process).** The C shim opens a TCP client to
   `material_server.py` (default `127.0.0.1:8765`) and exchanges one JSON frame
   per call. Easiest to develop and debug — the Python model runs in its own
   process and can be restarted independently.

Both use the *same* array layout and Voigt order, so switching is local to the C shim.

## Run / test

```bash
# start the Python material server
python ansys_usermat/coupling/material_server.py            # 127.0.0.1:8765
# end-to-end round-trip (server started in-process by the test)
python -m pytest tests/test_coupling.py
# syntax-check the Fortran hook
gfortran -c -fsyntax-only -ffixed-line-length-132 ansys_usermat/coupling/usermat_py_hook.f
```

## Verification note

`material_server.stress_core` is a line-by-line NumPy mirror of the verified
Abaqus core `umat_biofilm_visco.f` (`BIOFILM_STRESS_CORE`) — the same core the
`ansys_usermat/crosscheck/` harness proves bit-identical to the ANSYS port. The
Python core can be cross-validated against the compiled Fortran by driving both
with the same states (as `crosscheck.py` does). The included test checks the
end-to-end plumbing plus tangent shape/symmetry.

## Next steps (continuation)

1. Provide the C shim `biofilm_py_eval` (socket client first; ISO_C_BINDING embed later).
2. Single-element ANSYS/Abaqus smoke test with `kUsePy=1`; compare against the inline core.
3. Replace `stress_core` with the calibrated JAX model; keep the inline Fortran core as fallback.
