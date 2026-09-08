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

> **2026-09-07, same day: a multi-element extension
> (`t_growth_ecology_multi.dat`, 8 SOLID185 elements instead of 1) found
> and fixed two real concurrency bugs — RESOLVED.** Under ANSYS's default
> parallel execution the single-element deck's success did not reproduce:
> the run hung indefinitely. Two distinct bugs were involved, and the
> first fix found was necessary but not sufficient — the second, deeper
> one turned out to be the actual cause of the hang:
>
> 1. **`biofilm_py_eval.c` kept one global, unsynchronised TCP connection
>    (`g_fd`) shared by every caller.** If ANSYS ever evaluates the
>    material routine from multiple threads within one process, concurrent
>    send()/recv() on the same fd can interleave requests into malformed
>    frames and misdeliver a reply to the wrong caller. Real, and fixed
>    with a mutex (`SRWLOCK` on Windows, `pthread_mutex_t` elsewhere)
>    serialising the whole request/response exchange in both
>    `biofilm_py_eval()` and `biofilm_ecology_eval()` (they share `g_fd`,
>    so one lock covers both) — see `biofilm_py_eval.c`'s own note. This
>    costs no real parallelism (see point 2). Verified: compiles clean
>    under MinGW and MSVC, full pytest regression suite passes.
> 2. **The actual root cause: ANSYS parallelises this solve via MPI, not
>    (only) OpenMP threads.** A multi-element job spawned 4 *separate*
>    `ANSYS.exe` processes — confirmed by killing a hung run and reading
>    `"BAD TERMINATION ... RANK 0/1/2/3"` in its output — each opening its
>    *own* persistent connection to `material_server.py`
>    (`biofilm_py_eval.c`'s "one connection, reused for the whole run"
>    design, working exactly as intended, once per process). But the
>    server used a plain `socketserver.TCPServer`, which services one
>    accepted connection's *entire* lifetime — its handler's
>    `for line in self.rfile:` loop — before ever accepting the next.
>    Rank 0's connection monopolised the server for the rest of the run;
>    ranks 1–3 hung forever waiting for a response that could never come.
>    This produced a symptom nearly identical to bug 1 (both looked like
>    "hangs under default parallel execution"), which is why fixing bug 1
>    alone looked plausible but didn't resolve it — the CPU-bound (not
>    idle) hang that remained afterward was MPI's own barrier/collective
>    wait spinning while ranks 1–3 never returned. **Fixed** by switching
>    the server to `socketserver.ThreadingTCPServer`
>    (`material_server.py`'s new `_Server` class) — GIL still serialises
>    the actual computation, so this loses no real parallelism either, it
>    just stops the server from starving every client after the first.
>
> **With both fixes, `t_growth_ecology_multi.dat` completes in ~12s under
> ANSYS's default thread/rank count** (no `-np 1` needed): 0 errors, all
> 8 elements independently reproduce the single-element closed form
> exactly. Two diagnostic decks used to isolate the two bugs are kept
> alongside it: `t_growth_multi_baseline.dat` (8 elements, no Python hook
> at all — converged in 2.5s even *before* either fix, proving the bug was
> specific to the socket bridge, not multi-element custom-USERMAT
> execution in general) and `t_growth_kusepy_multi.dat` (8 elements,
> material hook only, no ecology — also hung before the
> `ThreadingTCPServer` fix and also converges after it, proving the bug
> was in the shared bridge infrastructure, not ecology-specific code, and
> so affects the *material* bridge's own future multi-element real-ANSYS
> use too). Evidence:
> [`out_ecology_multi.txt`](../apdl/out_ecology_multi.txt) /
> [`growth_result_ecology_multi.txt`](../apdl/growth_result_ecology_multi.txt),
> [`out_kusepy_multi.txt`](../apdl/out_kusepy_multi.txt) /
> [`growth_result_kusepy_multi.txt`](../apdl/growth_result_kusepy_multi.txt),
> [`out_multi_baseline.txt`](../apdl/out_multi_baseline.txt) /
> [`growth_result_multi_baseline.txt`](../apdl/growth_result_multi_baseline.txt).

> **2026-09-07, same day: extended off the synthetic F=I-forced decks onto
> real, unconstrained geometry -- `t_growth_cylinder_ecology.dat`, VERIFIED.**
> Everything above (`t_growth_ecology.dat`, `t_growth_ecology_multi.dat`) used
> a fully-constrained unit cube/cube, closed-form-checkable but geometrically
> nothing like a tooth surface. This deck reuses `t_growth_cylinder_shell.dat`'s
> real two-layer curved-shell geometry (bonded substrate + growth layer,
> genuinely deforming, no forced F=I) and splits the growth layer into two
> materials by initial ecology seed -- one all-zero (the same
> `INIT_ECO_IF_ZERO` default every other test uses), one an explicit
> synthetic composition with more biomass up front -- standing in for a
> spatially-varying CLSM measurement (still a placeholder; see next-steps
> below). No TMCMC calibration either: theta stays `THETA_DEMO`, deliberately
> -- this deck tests the *spatial/FEM* machinery, which is orthogonal to
> calibration accuracy (RESEARCH_MODEL.md sec.7 item 3).
>
> First attempt reused the base deck's full ARC=60/LEN=3.0 mesh (12240
> elements) and was killed after 14 minutes of near-zero CPU: every
> Gauss point of every growth-layer element makes a socket round trip per
> Newton iteration per substep, so element count costs far more here than
> in the pure-Fortran inline core. **ANSYS access on this machine is only
> good through 2026-09-08** (see the project's `ansys_access_window_2026-09`
> note), so the geometry was shrunk (ARC=5, LEN=0.15, same radii/ESIZE/BCs)
> to 54 elements -- still real, unconstrained, curved, two-material geometry,
> just a smaller patch of it.
>
> Result: 0 errors, 157s solve. Both regions' final alpha (`SVAR(10)`) match
> `ecology_cylinder_reference.py`'s independent 10-step chained
> `ecology_step` trajectory exactly to displayed precision (mat 2/region A:
> 4.5687e-3 vs. reference 4.568735e-3; mat 3/region B: 2.8868e-3 vs.
> reference 2.886820e-3) -- and the two regions are genuinely different from
> each other, a real non-uniform alpha(x) produced by live ecology on a
> deforming geometry, not a closed-form identity. A follow-up sensitivity
> check (`ecology_theta_sensitivity.py`) found this specific result --
> region A's default seed growing faster than region B's synthetic seed --
> is **stable across 9 theta variants** (weak/strong/sign-flipped/5 random),
> B/A ratio staying in 0.62-0.63 throughout: at this short integration
> horizon (10 steps of dt=1e-5) the outcome is dominated by the seed
> composition, not by theta, so using the uncalibrated `THETA_DEMO`
> placeholder here is unlikely to be hiding a qualitatively different
> answer -- though it may also mean the model's theta-sensitivity is weak
> at this horizon specifically, which is worth checking again over a longer
> integration once TMCMC calibration exists. See
> [`out_cylinder_ecology.txt`](../apdl/out_cylinder_ecology.txt) /
> [`growth_result_cylinder_ecology.txt`](../apdl/growth_result_cylinder_ecology.txt).

> **2026-09-07, same day: scale-invariance check —
> `t_growth_cylinder_ecology_big.dat`, VERIFIED.** Same two-material patch,
> same seeds, same theta, same fixed dt=1e-5×10 stepping, but ARC/LEN scaled
> up from 5/0.15 to 15/0.4 (54 → 432 elements, run from a robocopy'd working
> directory with the full custom `ANSYS.exe` + ~186 supporting DLLs, since
> the exe alone will not launch). Result: 0 errors, 6 benign warnings,
> both regions' alpha (`SVAR(10)`) uniform across all 432 elements and
> identical to the small-patch result to displayed precision (region A
> 0.45687E-002, region B 0.28868E-002) — exactly the expected outcome for a
> 0D per-Gauss-point reaction with no inter-element coupling, confirming the
> result is not an artifact of the shrunk-mesh smoke test's small size. See
> [`out_cylinder_ecology_big.txt`](../apdl/out_cylinder_ecology_big.txt) /
> [`growth_result_cylinder_ecology_big.txt`](../apdl/growth_result_cylinder_ecology_big.txt).

> **2026-09-07, same day: genuine 2D spatial variation —
> `t_growth_cylinder_ecology_4region.dat`, VERIFIED.** The base ecology
> deck only varies along theta (2 regions); this splits the growth layer
> into a 2x2 (theta x Z) grid of 4 materials instead, seeds A/B reused
> unchanged, C/D new and deliberately different-again. Result: 0 errors, 7
> benign warnings, all four regions' alpha (`SVAR(10)`) match
> `ecology_4region_reference.py`'s independent reference EXACTLY:
> A 0.45687E-002, B 0.28868E-002, C 0.13122E-005, D 0.17898E-002 — spanning
> **~3 orders of magnitude** despite identical theta and identical dt×10
> integration, and non-monotonic in initial biomass (C's `phi` sums to
> 0.20, higher than D's 0.10, yet C's alpha ends up ~1400x *smaller* than
> D's) — a much stronger demonstration that the live per-Gauss-point
> ecology bridge produces a genuinely composition-driven, non-uniform
> alpha(x) field on real curved geometry, not just a 1D theta split. See
> [`out_cylinder_ecology_4region.txt`](../apdl/out_cylinder_ecology_4region.txt) /
> [`growth_result_cylinder_ecology_4region.txt`](../apdl/growth_result_cylinder_ecology_4region.txt).

> **2026-09-07, same day: full 3D composition field —
> `t_growth_cylinder_ecology_8region.dat`, VERIFIED.** One more axis: the
> growth layer's radial half-split (R_MID vs R_OUT side) is not arbitrary —
> depth from the substrate is exactly what a real CLSM z-stack profiles, so
> this is a genuinely motivated third axis, not just "add more regions."
> 2x2x2 (R x theta x Z) grid, 8 materials, seeds A–D reused from the
> 4-region deck (near-substrate half), E–H new (near-surface half). Result:
> 0 errors, 11 benign warnings, all eight regions' alpha (`SVAR(10)`) match
> `ecology_8region_reference.py`'s independent reference exactly: A
> 4.5687e-3, B 2.8868e-3, C 1.3122e-6, D 1.7898e-3, E 1.5979e-3, F
> 1.6777e-3, G 1.4727e-3, H 2.1593e-6 — still spanning >3 orders of
> magnitude, with the two lowest values (C, H) both being single-dominant-
> species compositions, an emergent pattern under `THETA_DEMO` rather than
> something picked for that outcome. This is the most spatially-realistic
> configuration the live ecology bridge has been exercised on to date: a
> genuine 3D (R, theta, Z) composition-driven alpha(x) field on real,
> unconstrained, curved geometry. See
> [`out_cylinder_ecology_8region.txt`](../apdl/out_cylinder_ecology_8region.txt) /
> [`growth_result_cylinder_ecology_8region.txt`](../apdl/growth_result_cylinder_ecology_8region.txt).

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
| `material_server.py` | Python side — NumPy reference core (mirrors the verified Fortran `BIOFILM_STRESS_CORE`) + F-perturbation tangent + a socket server (`_Server`, `ThreadingTCPServer` — services multiple concurrent persistent connections, e.g. one per MPI rank in a parallel ANSYS solve; see the 2026-09-07 Status note). `--tangent jax` swaps the tangent for the exact AD one below; the default stays `fd` so `kUsePy=1` vs `kUsePy=0` remains an exact equivalence check. |
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
| `../apdl/t_growth_ecology_multi.dat` | 8-element extension of the above — resolved (~12s, 0 errors, default MPI thread/rank count) after the two concurrency fixes; see the second 2026-09-07 Status note |
| `../apdl/t_growth_kusepy_multi.dat`, `../apdl/t_growth_multi_baseline.dat` | diagnostic decks that isolated the two concurrency bugs (material hook only; no Python hook at all) — kept as regression evidence |
| `../apdl/t_growth_cylinder_ecology.dat` | real, unconstrained curved two-material geometry (extends `t_growth_cylinder_shell.dat`) with live ecology driving two spatially distinct growth-layer regions — verified; see the third 2026-09-07 Status note |
| `../apdl/ecology_cylinder_reference.py` | independent reference for the above — chains `ecology_jax.ecology_step` 10x (matching the deck's 10 fixed dt=1e-5 substeps) for each region's seed |
| `../apdl/ecology_theta_sensitivity.py` | checks whether the cylinder-ecology deck's region-A-grows-faster result depends on `THETA_DEMO`'s specific numbers (it doesn't, across 9 variants tried) — see the Status note |
| `../apdl/t_growth_cylinder_ecology_big.dat` | same patch scaled ~8x (54→432 elements) — scale-invariance check, identical alpha result — see the fourth 2026-09-07 Status note |
| `../apdl/t_growth_cylinder_ecology_4region.dat` | growth layer split into a 2x2 (theta x Z) grid of 4 materials/seeds instead of 2 — genuine 2D spatial variation, alpha spans ~3 orders of magnitude — see the fifth 2026-09-07 Status note |
| `../apdl/ecology_4region_reference.py` | independent reference for the above — same chained `ecology_step` logic, 4 seeds |
| `../apdl/t_growth_cylinder_ecology_8region.dat` | full 3D (R x theta x Z) 2x2x2 grid, 8 materials/seeds — genuine depth-resolved composition field, alpha spans >3 orders of magnitude — see the sixth 2026-09-07 Status note |
| `../apdl/ecology_8region_reference.py` | independent reference for the above — same chained `ecology_step` logic, 8 seeds |

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
   IKMHIWI03, gfortran driver AND real ANSYS (single-element, multi-element,
   AND real unconstrained curved geometry — `t_growth_cylinder_ecology.dat`),
   see the Status notes above. What's left on this specific path (not
   blocking, but needed before the result means anything physically):
   - **Feed the actual TMCMC-calibrated theta** into `prop(10:29)` instead
     of the demo's `THETA_DEMO` placeholder — needs a generator analogous
     to `composition_to_material.py`'s `apdl_state_block`. **Still blocked**
     as of 2026-09-08: no `theta_MAP.json`/TMCMC posterior run output exists
     on this machine (`tmcmc_to_fem_coupling.py`/`posterior_ci_0d.py` both
     expect one under `../data_5species/_runs/<sweep>/`, which isn't
     present here — the actual run lives elsewhere, e.g. the Keio server).
   - ~~Seed `ustatev(15:26)` from real CLSM-measured composition per Gauss
     point instead of `INIT_ECO_IF_ZERO`'s fixed default~~ — **partially
     done 2026-09-08, real ANSYS, `phi` only**: `apdl/t_growth_ecology_clsm_phi.dat`
     seeds `phi(1:5)`+`phi0` from the real Day-1 (Tag=1) Commensal/Static
     composition in `data/heine_species_distribution_biofilm.xlsx`
     ("all cells" sheet), regenerable via `apdl/clsm_phi_seed_reference.py`.
     Ran with 0 errors; every `SVAR`/stress value matches the independent
     Python reference to displayed precision (`apdl/growth_result_ecology_clsm_phi.txt`).
     `psi(1:5)` is **still the 0.999 placeholder**, deliberately — the
     workbook's "only living cells" sheet gives a ratio
     (`plot_heine_phi_psi.py`'s own "relative viability/enrichment", values
     up to ~30 seen) that is a *different quantity* from the model's own
     `psi_i ∈ [0,1]` (membrane-intact fraction, RESEARCH_MODEL.md sec.1) —
     using it directly would misrepresent measured data as something it
     doesn't measure; see the deck's own header comment for the numbers
     that ruled it out (two species' Day-1 ratio exceeds 1). A real
     bounded per-species viability measurement, if one exists separately
     from this workbook, would close this the rest of the way.
   - **Extended to real, unconstrained, curved geometry, 2026-09-08**:
     `apdl/t_growth_cylinder_ecology_clsm.dat` reuses
     `t_growth_cylinder_ecology.dat`'s two-region curved-shell setup, but
     replaces region B's (mat 3) synthetic "not measured data" seed with
     the same real Day-1 CLSM composition above — region A (mat 2, the
     default `INIT_ECO_IF_ZERO` seed) is left unchanged as a regression
     baseline. Ran with 0 errors; both regions' `alpha` (region A
     4.5687e-3, region B/real-CLSM 4.6145e-3) and region B's full ecology
     state match `ecology_cylinder_reference_clsm.py`'s independent
     10-substep-chained reference exactly — see
     `out_cylinder_ecology_clsm.txt` / `growth_result_cylinder_ecology_clsm.txt`.
     So the real-CLSM seed now has both a fully-constrained closed-form
     check (`t_growth_ecology_clsm_phi.dat`) and a genuinely spatial,
     unconstrained-deformation check on curved geometry.
   - **Extended to three of the four clinical conditions, 2026-09-08**:
     `apdl/t_growth_cylinder_ecology_4region_clsm.dat` reuses
     `t_growth_cylinder_ecology_4region.dat`'s 2×2 (theta×Z) 4-region
     layout, replacing regions B/C/D's (mat 3/4/5) synthetic seeds with
     real Day-1 CLSM compositions for **CS** (Commensal/Static, same seed
     already verified above), **CH** (Commensal/HOBIC), and **DS**
     (Dysbiotic/Static) — region A (mat 2, the default seed) stays the
     unchanged regression baseline. Ran with 0 errors; all four regions'
     `alpha` (A 4.5687e-3, B/CS 4.6145e-3, C/CH 4.6920e-3, D/DS 4.3557e-3)
     match `ecology_4region_reference_clsm.py`'s independent reference
     exactly.
     **This deck's own header originally said DH (Dysbiotic/HOBIC) had "no
     measurement at all" for F. nucleatum/P. gingivalis and was excluded
     for that reason — that claim was WRONG, corrected the same day after
     being asked to double-check ("DHないっけ").** It was an artifact of a
     real bug in `plot_heine_phi_psi.load()`: the function derived the
     species-block column width once from a sheet's first header row and
     reused it for every condition block in that sheet, but the Dysbiotic
     sheet's `HOBIC ...` blocks are 9 columns/species wide while its
     `Static ...` blocks are 18 — reusing width=18 for HOBIC didn't just
     drop 2 species, it read species 2's real data into the "species 1"
     slot, species 3's into "species 2", species 5's into "species 3", and
     landed the last two slots on blank padding. **Fixed** in
     `plot_heine_phi_psi.py` (and the same latent bug in
     `plot_heine_composition.py` / `validate_composition.py`), figures
     regenerated. See the next item for DH included.
   - **All four clinical conditions, real CLSM, 2026-09-08 (same day,
     later)**: `apdl/t_growth_cylinder_ecology_4region_all_real_clsm.dat`
     retires the "region A = default placeholder" convention now that DH
     is available too — mat 2/3/4/5 = CS/CH/DS/DH, all real Day-1 CLSM.
     Ran with 0 errors; all four regions' `alpha` (CS 4.6145e-3, CH
     4.6920e-3, DS 4.3557e-3, DH 4.8186e-3) match
     `ecology_4region_reference_all_real_clsm.py`'s independent reference
     exactly — `out_cylinder_ecology_4region_all_real_clsm.txt` /
     `growth_result_cylinder_ecology_4region_all_real_clsm.txt`. All four
     of this repo's clinical conditions now drive visibly different
     growth, from real measured composition, on one real spatial FEM
     geometry.
5. Port this bridge (hook + shim + server) to Oliver's `Usermat_P21-V21_*.F`,
   the same way `BIOFILM_GROWTH_VISCO_V01` was wired in directly on
   2026-09-03 (`V222_PORT_INSTRUCTIONS.md` §6).
