# Gauss-point material bridge — skeleton

Scaffolding for the thesis' core deliverable: calling the calibrated **Python**
material model from the Fortran **UMAT/USERMAT** at each Gauss point, with the
verified inline Fortran core kept as the reference/fallback. This lets the
IKM-verified constitutive law and the Keio-side computational-mechanics work
meet at a single, well-defined interface.

> **Status: verified end-to-end through the real `usermat()` entry point,
> now including a live ANSYS run, not only the gfortran driver.** The Python
> material core is proven equivalent to the verified Fortran core, the C
> shim talks to it over a socket, the Fortran hook links against the shim
> and falls back cleanly when the server is absent, and `kUsePy=1` matches
> `kUsePy=0` when run through `usermat_biofilm.f`'s actual `usermat()`
> subroutine — see `../../tests/test_usermat_kusepy_e2e.py`.
>
> **2026-09-03, on IKMHIWI03: the single-element ANSYS smoke test from
> Next-steps #1 is done, on real hardware.** Two blockers had to be fixed
> first, neither specific to this material law:
> - `biofilm_py_eval.c` used POSIX sockets (`arpa/inet.h`, `sys/socket.h`,
>   `unistd.h`) and does not compile under MSVC. Ported to Winsock2 behind
>   `#ifdef _WIN32` (kept the POSIX path for Linux); `send()`/`recv()`
>   replace `write()`/`read()` on both platforms rather than maintaining
>   two implementations of the framing helpers.
> - `usermat_py_hook.f` needs 132-column fixed form (its own header already
>   said so, for `gfortran -ffixed-line-length-132`) but ifort defaults to
>   72 and misparses the overflow as a syntax error on the `bind(C,
>   name=...)` continuation line rather than reporting a line-length
>   problem. Fixed by adding `/extend-source:132` to `link_v222.ps1` (now
>   permanent there — safe for the files that already fit in 72 columns).
>
> With both fixed, `usermat_biofilm.f` + `usermat_py_hook.f` +
> `biofilm_py_eval.c` compiled and linked into a custom v222 `ANSYS.exe`
> (`link_v222.ps1 -WorkDir F:\biofilm_upf_kusepy`) and run against the
> `elastic_a005` closed-form case (`α=0.05, η=0`) with `prop(6)=kUsePy=1`
> and `material_server.py` listening on `127.0.0.1:8765`:
> - `NUMBER OF ERROR MESSAGES = 0`; `SX=SY=SZ=-0.10193E-003` — exact match
>   to the closed form (`-1.019275856e-04`), `SVAR(10)=0.05000`.
> - **The live socket path was confirmed actually taken, not a silent
>   fallback that happens to agree**: a diagnostic per-request log line
>   added to a scratch copy of `material_server.py` recorded real requests
>   arriving from `127.0.0.1`, `alpha=0.05...`, one per F-perturbation
>   tangent evaluation, across the run's substeps.
> - Re-ran with the server killed: `NUMBER OF ERROR MESSAGES = 0`, bit-
>   identical stress — the fallback-to-inline-core path also confirmed on
>   real ANSYS, not only the gfortran driver.
>
> **Same day, follow-up: the gfortran/MinGW driver path (this file's own
> `test_usermat_kusepy_e2e.py` etc.) was actually still broken on Windows
> at that point** — the real-ANSYS run above used MSVC-built decks and
> didn't exercise it. Verifying it turned up two more pre-existing,
> Windows-specific bugs, both unrelated to sockets: (1) every test in this
> family hardcoded a POSIX-only `PATH` for the child process env, which
> starves a MinGW-built `.exe` of its runtime DLLs on Windows
> (`STATUS_DLL_NOT_FOUND`) — fixed by inheriting the parent environment
> instead; (2) this particular MinGW-w64 build lacks the POSIX socket
> headers other MinGW distributions ship, confirmed empirically, so it
> needs the Winsock2 path too (kept the `_WIN32` guard rather than
> narrowing to `_MSC_VER`). With both fixed, the whole coupling e2e suite
> passes on IKMHIWI03 — see the commit fixing `tests/test_*` for detail.
>
> Not yet done: porting this same bridge to Oliver's `Usermat_P21-V21_*.F`
> (a separate task from wiring `BIOFILM_GROWTH_VISCO_V01` there directly,
> which was tried the same day — see `V222_PORT_INSTRUCTIONS.md` §6) and
> replacing `stress_core` with the calibrated JAX model (§ Next steps
> below, now down to items 2–3).

> **2026-09-07: the second, independent Gauss-point bridge — the 0D
> Hamilton ecology ODE (`ecology_jax.py`) — is wired in and verified
> end-to-end, gfortran driver AND real ANSYS, the same rigor as the
> material bridge above.** This is next-steps item 4 below: instead of
> reading a precomputed alpha field, `prop(8)=kUseEcology=1` makes
> `usermat()` advance a per-Gauss-point 0D Hamilton state
> (`ustatev(15:26)`, the same Newton-per-step integrator as
> `jax_hamilton_0d_5species_demo.py`) by one increment every call, and
> derives alpha's increment from `k_alpha * sum(phi_i*psi_i)`. Same
> connection/server as the material hook, disambiguated on the wire by a
> `"kind":"ecology"` request field (`protocol.py`); independent of
> `kUsePy` — this drives the growth input, not the stress law.
>
> **Bug found and fixed during the gfortran-driver test, before ANSYS was
> even touched:** the `bind(C)` interface for `biofilm_ecology_eval`
> declared `dt_h` (a plain `double`, by-value in the C prototype) without
> the `VALUE` attribute, so Fortran passed its *address* instead — a
> silent ABI mismatch (wrong register/stack slot for a by-value scalar)
> that corrupted the whole call rather than crashing. It surfaced as an
> unexplained clean fallback (`ierr != 0`) with no error message, not an
> obvious bug — isolated by comparing a standalone C-only test (worked)
> against the same call through the Fortran driver (silently took the
> fallback path) until the two diverged only at the interface declaration.
>
> **2026-09-07, real ANSYS on IKMHIWI03** (`link_v222.ps1 -WorkDir
> F:\biofilm_upf_kusepy`, reusing the material bridge's own working
> directory/scaffolding): `t_growth_ecology.dat` — fully constrained
> single element (`F=I` forced, same trick as `t_growth_kusepy.dat`),
> `ustatev(15:26)` arriving all-zero so `INIT_ECO_IF_ZERO` seeds the
> default composition, one step of `dt=1e-5` against `THETA_DEMO` and
> `k_alpha=50`:
> - `NUMBER OF ERROR MESSAGES = 0`; `SX=SY=SZ=-0.29998`, all six
>   `SVAR(15:26)` (the full ecology state) and `SVAR(10)=0.49714E-003`
>   (alpha) — **exact match**, to displayed precision, against the
>   independent Python reference (`ecology_jax.default_initial_state()` →
>   `ecology_jax.ecology_step(..., THETA_DEMO, 1e-5)` →
>   `ecology_jax.living_fraction_total(...)` → `alpha = 1e-5*50*phi_tot` →
>   `material_server.stress_core(...)`). Evidence:
>   [`out_ecology.txt`](../apdl/out_ecology.txt) /
>   [`growth_result_ecology.txt`](../apdl/growth_result_ecology.txt).
> - **Live socket path confirmed actually taken**: `material_server.py`'s
>   stderr shows the (benign) end-of-run `ConnectionResetError` that only
>   happens after a request from `127.0.0.1` was received and answered —
>   the same signature the material bridge's own 2026-09-03 verification
>   used as proof of a live, not silently-skipped, call.
> - Re-ran with the server killed: `NUMBER OF ERROR MESSAGES = 0`,
>   `SX=SY=SZ=0` and every `SVAR` at its input value (alpha stays 0, the
>   ecology state stays all-zero) — the fail-safe fallback confirmed on
>   real ANSYS too, not just the driver. Evidence:
>   [`out_ecology_fallback.txt`](../apdl/out_ecology_fallback.txt) /
>   [`growth_result_ecology_fallback.txt`](../apdl/growth_result_ecology_fallback.txt).
>
> Not yet done for this path: feeding a TMCMC-calibrated theta (currently
> the demo's `THETA_DEMO` placeholder in `prop(10:29)`) and a real
> CLSM-measured initial composition (currently `INIT_ECO_IF_ZERO`'s fixed
> default) instead of both being placeholders — see next-steps below.

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

A second, independent evaluation shares the same server/connection — the 0D
Hamilton ecology ODE step (`ecology_jax.py`), disambiguated by a `"kind"`
field on the wire (absent = the material request above, backward compatible):

| direction | quantity | shape | note |
|---|---|---|---|
| in | `g` (ecology state) | 12 | phi(5), phi0, psi(5), gamma |
| in | `theta` (interaction params) | 20 | 15 independent A_ij + 5 b_i, TMCMC-calibrated |
| in | `dt_h` | scalar | increment |
| out | `g_new` | 12 | one Newton-per-step 0D Hamilton advance |

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
| `ecology_jax.py` | 0D Hamilton ecology ODE bridge — thin wrapper around `jax_hamilton_0d_5species_demo.py`'s Newton-per-step integrator (`ecology_step`), plus `living_fraction_total` (the growth reaction driver) and `default_initial_state`. No physics duplicated: imports the demo's `newton_step_jit`/`theta_to_matrices` directly. |
| `usermat_ecology_e2e_driver.f` | standalone driver calling the real `usermat()` subroutine with `prop(8)=kUseEcology` toggled — used by `test_usermat_ecology_e2e.py` |
| `test_shim_ecology_main.c` | tiny C driver exercising `biofilm_ecology_eval()`, used by `test_ecology_shim.py` |
| `../../tests/test_ecology_coupling.py` | Python-side round-trip (protocol + in-process + socket dispatch by `"kind"`) against the verified 0D reference (CI) |
| `../../tests/test_ecology_shim.py` | **C-shim end-to-end** for the ecology path, mirrors `test_coupling_shim.py` (CI) |
| `../../tests/test_usermat_ecology_e2e.py` | **full-chain end-to-end** — `usermat()` with `kUseEcology=1`: `g_new`/alpha match the Python reference bit-for-bit through the wire, the all-zero-state default seed, `kUseEcology=0` leaves state untouched, and the no-server fallback (CI) |
| `../apdl/t_growth_ecology.dat` | real-ANSYS single-element smoke test for the ecology hook (fully constrained, `F=I`) — see the 2026-09-07 Status note above |

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

1. ~~Single-element ANSYS smoke test with `kUsePy=1`, on real hardware~~ —
   done 2026-09-03 on IKMHIWI03, see the Status note above.
2. **Replace `stress_core` with the calibrated JAX model**, keeping the inline
   Fortran core as the fallback; re-run `test_coupling_vs_fortran.py` to quantify
   the intended physical difference.
3. Optional: switch the shim from socket to in-process `ISO_C_BINDING` embedding
   once the model is stable (lower per-Gauss-point latency).
4. ~~Live per-Gauss-point Python call for the growth driver (e.g. the 0D
   Hamilton ODE) rather than a precomputed field~~ — done 2026-09-07 on
   IKMHIWI03, gfortran driver AND real ANSYS, see the Status note above.
   What's left on this specific path (not blocking, but needed before the
   result means anything physically):
   - **Feed the actual TMCMC-calibrated theta** into `prop(10:29)` instead
     of the demo's `THETA_DEMO` placeholder — needs a generator analogous
     to `composition_to_material.py`'s `apdl_state_block`.
   - **Seed `ustatev(15:26)` from real CLSM-measured composition** per
     Gauss point instead of `INIT_ECO_IF_ZERO`'s fixed default — the same
     kind of `TB,STATE` block `composition_to_material.py` already emits
     for `ustatev(11:14)`.
5. Port this bridge (hook + shim + server) to Oliver's `Usermat_P21-V21_*.F`,
   the same way `BIOFILM_GROWTH_VISCO_V01` was wired in directly on
   2026-09-03 (`V222_PORT_INSTRUCTIONS.md` §6).
