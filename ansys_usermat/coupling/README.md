# Gauss-point material bridge — skeleton

Scaffolding for the thesis' core deliverable: calling the calibrated **Python**
material model from the Fortran **UMAT/USERMAT** at each Gauss point, with the
verified inline Fortran core kept as the reference/fallback. This lets the
IKM-verified constitutive law and the Keio-side computational-mechanics work
meet at a single, well-defined interface.

> **Status: the bridge is complete and verified end-to-end** (no Abaqus/ANSYS
> needed): the Python material core is proven equivalent to the verified Fortran
> core, the C shim talks to it over a socket, and the Fortran hook links against
> the shim and falls back cleanly when the server is absent. What remains is
> wiring it into a live solver (§ Next steps).

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
| `usermat_py_hook.f` | Fortran side — `ISO_C_BINDING` interface to the C shim `biofilm_py_eval` + array marshalling; the `PYTHON MATERIAL HOOK` call site. Falls back to the inline core on failure. |
| `biofilm_py_eval.c` | **C shim** — persistent local TCP connection to the material server, one JSON frame per Gauss point, one reconnect retry, NaN guard. Returns nonzero on any failure so Fortran falls back. Host/port via `BIOFILM_PY_HOST` / `BIOFILM_PY_PORT`. |
| `test_shim_main.c` | tiny C driver used by the shim test |
| `../../tests/test_coupling.py` | Python-side round-trip through the protocol (CI) |
| `../../tests/test_coupling_vs_fortran.py` | **equivalence proof** — compiles the real Fortran core and compares it against the Python core over 28 states (CI) |
| `../../tests/test_coupling_shim.py` | **C-shim end-to-end** — compiles the shim, drives it against a live server, and checks it fails cleanly with no server (CI) |

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

## Verification status

The chain is closed — swapping the Fortran law for the Python model at the Gauss
point provably does not change the physics:

```
Abaqus UMAT  ==  ANSYS USERMAT  ==  Python material core
   (0 ULP, crosscheck/)        (6.8e-14 relative, this dir)
```

- ✅ **Python core ≡ Fortran core.** `test_coupling_vs_fortran.py` compiles the
  real `BIOFILM_STRESS_CORE` and drives both over 28 states (named corner cases +
  random finite strains). **Worst relative discrepancy 6.8e-14** on stress, `Fv`
  and `detFe`. The battery deliberately spans both regimes: 18 well-conditioned
  states and 10 degenerate ones where a large viscous step drives `Fv` singular
  and the `detFe` clamp fires (stress ~1e30 — non-physical but a real code path).
  Comparison is *relative*, since an absolute tolerance is meaningless at 1e30.
- ✅ **C shim ↔ Python server.** `test_coupling_shim.py` compiles the shim, drives
  it against a live server, and confirms the returned stress/`Fv`/`dsdePl` equal
  the in-process evaluation — and that it exits nonzero (no hang) with no server.
- ✅ **Fortran hook ↔ C shim link.** `usermat_py_hook.f` compiles and links
  against `biofilm_py_eval.c` (ABI symbol resolves); with no server running the
  hook returns `ok = .false.`, so the solver falls back to the inline core.

## Next steps (continuation)

1. **Single-element ANSYS/Abaqus smoke test** with `kUsePy=1`; compare against the
   inline core (should match to machine precision — the equivalence is proven above,
   so any mismatch is interface, not physics).
2. **Replace `stress_core` with the calibrated JAX model**, keeping the inline
   Fortran core as the fallback; re-run `test_coupling_vs_fortran.py` to quantify
   the intended physical difference.
3. Optional: switch the shim from socket to in-process `ISO_C_BINDING` embedding
   once the model is stable (lower per-Gauss-point latency).
