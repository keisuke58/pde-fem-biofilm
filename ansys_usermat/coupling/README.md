# Gauss-point material bridge — skeleton

Scaffolding for the thesis' core deliverable: calling the calibrated **Python**
material model from the Fortran **UMAT/USERMAT** at each Gauss point, with the
verified inline Fortran core kept as the reference/fallback. This lets the
IKM-verified constitutive law and the Keio-side computational-mechanics work
meet at a single, well-defined interface.

> **Status: the bridge is wired and verified end-to-end through the real
> `usermat()` entry point** (no Abaqus/ANSYS needed): the Python material core
> is proven equivalent to the verified Fortran core, the C shim talks to it
> over a socket, the Fortran hook links against the shim and falls back
> cleanly when the server is absent, and `kUsePy=1` matches `kUsePy=0` when
> run through `usermat_biofilm.f`'s actual `usermat()` subroutine (not a
> bypass driver) — see `../../tests/test_usermat_kusepy_e2e.py`. What remains
> is a live ANSYS single-element smoke test and swapping in the calibrated
> JAX model (§ Next steps).

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
| `material_server.py` | Python side — NumPy reference core (mirrors the verified Fortran `BIOFILM_STRESS_CORE`) + F-perturbation tangent + a socket server. `--tangent jax` swaps the tangent for the exact AD one below; the default stays `fd` so `kUsePy=1` vs `kUsePy=0` remains an exact equivalence check. |
| `material_jax.py` | JAX mirror of that core plus an **exact tangent** (`jax.jacfwd`, no finite-difference step) and `dsigma_dparams` — ∂σ/∂θ for posterior/UQ propagation. See the note below on what this is and is not worth. |
| `composition_to_material.py` | CLSM composition φ → per-integration-point `C10, C01, D1, eta` (the `kStateMat=1` path), plus the `TB,USER`/`TB,STATE` block that delivers them. No per-increment Python call — see `../README.md`. |
| `protocol.py` | wire schema (newline-delimited JSON; one request→one response). Swap for a binary frame later without touching the physics. |
| `usermat_py_hook.f` | Fortran side — `ISO_C_BINDING` interface to the C shim `biofilm_py_eval` + array marshalling; the `PYTHON MATERIAL HOOK` call site. Falls back to the inline core on failure. `dsde` is built as `transpose(reshape(d36,[6,6]))` — the wire carries a row-major (NumPy/C-order) flatten, Fortran's `RESHAPE` fills column-major, so the plain (untransposed) reshape silently returned the tangent's transpose. |
| `biofilm_py_eval.c` | **C shim** — persistent local TCP connection to the material server, one JSON frame per Gauss point, one reconnect retry, NaN guard. Returns nonzero on any failure so Fortran falls back. Host/port via `BIOFILM_PY_HOST` / `BIOFILM_PY_PORT`. |
| `test_shim_main.c` | tiny C driver used by the shim test |
| `usermat_endtoend_driver.f` | standalone driver calling the **real `usermat()`** subroutine (not `BIOFILM_STRESS_CORE`), toggling `kUsePy` via `prop(6)` — used by `test_usermat_kusepy_e2e.py` |
| `../../tests/test_coupling.py` | Python-side round-trip through the protocol (CI) |
| `../../tests/test_coupling_vs_fortran.py` | **equivalence proof** — compiles the real Fortran core and compares it against the Python core over 28 states (CI) |
| `../../tests/test_coupling_shim.py` | **C-shim end-to-end** — compiles the shim, drives it against a live server, and checks it fails cleanly with no server (CI) |
| `../../tests/test_usermat_kusepy_e2e.py` | **full-chain end-to-end** — compiles `usermat_biofilm.f` + `usermat_py_hook.f` + `biofilm_py_eval.c` + this driver, and checks `kUsePy=1` matches `kUsePy=0` (stress/`Fv`/`dsdePl`) through the actual `usermat()` entry point, across elastic/viscous/Mooney-Rivlin cases, plus the no-server fallback (CI) |

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

Reproduced live in [`python_core_vs_fortran_verification.ipynb`](python_core_vs_fortran_verification.ipynb)
(bilingual EN/JA) — same 28-case battery as `test_coupling_vs_fortran.py`, run
against a freshly-compiled Fortran core rather than restated from this file.

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
- ✅ **`kUsePy` branch of the real `usermat()` entry point.** `usermat_biofilm.f`
  now has a live (not commented-out) `PYTHON MATERIAL HOOK` body: it calls
  `biofilm_py_hook`, reindexes the Abaqus-order response to ANSYS order via
  `MAP6`, and writes `stress`/`ustatev`/`dsdePl`. `test_usermat_kusepy_e2e.py`
  compiles the actual `usermat()` subroutine and drives it with `kUsePy=0` and
  `kUsePy=1` — stress and the updated viscous state match to numerical
  precision, `dsdePl` matches to floating-point-noise precision (both sides
  use the same F-perturbation scheme, `PERT=1e-7`), and the no-server fallback
  path is confirmed too. Found and fixed a real bug along the way — see the
  `usermat_py_hook.f` row above.

## Exact (AD) tangent — what it is and is not worth

`material_jax.py` computes the same 6×6 the FD path does, but by forward-mode
AD, so there is no step to pick. Measured rather than assumed
(`tests/test_material_jax.py` pins each of these):

- ❌ **Not a convergence fix.** At the USERMAT's `PERT=1e-7` the FD tangent is
  already within **~3e-8 relative** of exact — far tighter than anything
  Newton's convergence *rate* responds to. This specifically **rules the
  tangent out** as the explanation for the cylinder-shell case that stops
  converging at α=0.015; that has to be looked for elsewhere.
- ❌ **Step-size optimum is not material-dependent** — worth recording because
  the opposite is the natural guess. σ scales with `C10`, so `C10` cancels out
  of the relative error and the optimum sits at `h=1e-8` for both the stiffest
  (CH, `C10`=166 Pa) and softest (DS, `C10`=5.4 Pa) condition, despite the ~31×
  spread the composition path introduces. One global `PERT` is fine.
- ✅ **No tuned magic number, no truncation error.**
- ✅ **Independent validation of the FD tangent** — and so of the Fortran one,
  which is 0-ULP identical to the NumPy core.
- ✅ **∂σ/∂θ (`dsigma_dparams`)** — the one thing the FD path cannot practically
  match. Propagating the TMCMC posterior needs the stress response to the
  calibrated parameters; by differencing that would mean extra solves (or extra
  socket round-trips) per Gauss point. Validated column-by-column against
  central differences on the independent NumPy core.

One implementation note worth carrying forward: at `eta = 0` the viscous arm
divides by zero, and a plain `jnp.where` gives the **correct value but a NaN
derivative** — the classic JAX where-trap. The guard is a second `jnp.where` on
the denominator, and the test asserts on the *derivative*, since a value-only
check passes either way.

## Next steps (continuation)

1. **Single-element ANSYS smoke test with `kUsePy=1`**, on real hardware
   (Windows/`IKMHIWI03`) rather than the gfortran standalone driver used here —
   needs `usermat_py_hook.f` and `biofilm_py_eval.c` copied into the UPF build
   directory alongside `usermat_biofilm.f` and compiled/linked together
   (`usermat_py_hook.f` must build first, for its `.mod` file). Not yet done;
   the standalone-driver equivalence proven here should carry over directly,
   so any mismatch there would be build/link, not physics.
2. **Replace `stress_core` with the calibrated JAX model**, keeping the inline
   Fortran core as the fallback; re-run `test_coupling_vs_fortran.py` to quantify
   the intended physical difference.
3. Optional: switch the shim from socket to in-process `ISO_C_BINDING` embedding
   once the model is stable (lower per-Gauss-point latency).
