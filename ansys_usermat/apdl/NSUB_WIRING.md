# Sub-stepping the ecology: what to run on the machine

2026-09-29. The Python side was verified without ANSYS; this is the part that
needs it. Read `V222_PORT_INSTRUCTIONS.md` §8 first if the failure this fixes
is not fresh.

## What it fixes

The ecology ODE is only trusted below `dt = 1e-4` (measured — a Python sweep
turns erratic from ~1.5e-3). The partner's NEM deck runs at `TIME INC = 0.1`,
its stepping set by Workbench-compiled solution controls rather than an
editable `NSUBST`, and ANSYS's own bisection does not reach a 1000x
reduction. The guard therefore refused every step and the run could not
complete.

The increment is now divided **inside one material call**: one socket round
trip, `n_sub` cheap ODE steps within it. Not a new scheme — every reference
under `apdl/` already advances the ODE by chaining `ecology_step`.

## What changed, and why each piece

| File | Change |
|---|---|
| `coupling/protocol.py` | optional `n_sub` in the request; `phi_int` in the response |
| `coupling/material_server.py` | loops `ecology_step` `n_sub` times, accumulating `phi_int` |
| `coupling/biofilm_py_eval.c` | emits `n_sub`, parses `phi_int` |
| `coupling/usermat_py_hook.f` | `biofilm_ecology_hook` gains `dtmax` and `phi_int`; computes `n_sub = ceil(dt/dtmax)` |
| `usermat_biofilm.f` | passes `DTMAX_ECO = 1e-4`; forms alpha from `phi_int` |

**`phi_int` is the part that is easy to get wrong.** Growth accumulates as
`k_alpha * integral(phi_tot dt)`. With one step, `dt * phi_tot(g_new)` is that
integral. With `n_sub > 1` it is not — the end state does not stand for the
whole increment — so the server returns the sum over sub-steps and the
material forms alpha from it. Reusing the old expression would have committed
a wrong growth increment silently, which is exactly the class of failure this
routine is built to refuse.

## Build

`link_v222.ps1` drives ifort only, so the C shim is compiled separately with
`cl.exe`, as when the bridge was first linked:

```
cl /c /O2 ansys_usermat\coupling\biofilm_py_eval.c
```

then relink as usual. Delete the stale `ANSYS.exe` / `.lib` / `.exp` / `.map`
first or `ansys.lrf`'s `*.lib` wildcard collides with the new output.

## Run, in this order

**1. Regression first — the decks that already pass.** Sub-stepping must not
change them: they run at `dt = 1e-5`, below `DTMAX_ECO`, so `n_sub = 1` and
the request is byte-identical to before.

```
t_growth_cylinder_ecology_4region_all_real_clsm.dat
```

Expect: 0 errors, and the four regions' alpha unchanged —
`4.6145e-3 / 4.6920e-3 / 4.3557e-3 / 4.8186e-3`. **If these move, stop.**
Something other than sub-stepping changed.

**2. The deck this was built for.** Read the next section first — this deck
needs an edit on the partner side that is not in this repository.

```
ds_oliver_wired_ecology.dat      (F:\biofilm_upf_wired)
```

Previously: 0 errors but `PROBLEM TERMINATED`, element 8 repeatedly
triggering the guard's cut-back. With the edit below, `dt = 0.1` gives
`n_sub = 1000`, the blanket refusal is gone, and the state check that remains
should not fire.

What to record either way:

- `NUMBER OF ERROR MESSAGES`
- whether it reaches the end of the load step
- `SVAR(10)` (alpha) and `SVAR(84)` on a few elements — a completed run with
  alpha stuck at its seed means the hook is not being reached
- wall time. 1000 ODE steps per Gauss point per increment is the cost being
  traded for the round trip; if it is unusable, that is a real result and
  worth knowing before Thursday

**3. Only if 1 and 2 are clean:** the tooth/implant geometry, which is the
roadmap's actual open item and a larger job.

## The one edit that is not in this repository

`usermat_biofilm.f` is this repo's own routine, and it is done. The deck in
step 2 does not go through it — it goes through the partner's NEM material
routine, whose source lives only on the machine (`F:\biofilm_upf_wired`) and
is deliberately not committed here. That routine has its own copy of the
ecology call site, and per `V222_PORT_INSTRUCTIONS.md` §8 it currently:

1. **refuses outright** when `dt > DTMAX_ECO`, and
2. sanity-checks the returned state, holding alpha at `ustatev(84)` and
   asking for a cut-back (`keycut=1`) on either failure.

Item 1 is the thing that stopped the run, and sub-stepping does not remove it
by itself: the refusal happens *before* the hook is reached, so with it still
in place `dt = 0.1` fails exactly as it did on 09-07 and the new code never
executes. Three lines change there, mirroring `usermat_biofilm.f`:

- drop the `dt > DTMAX_ECO` refusal, and pass `DTMAX_ECO` to the hook as its
  new `dtmax` argument instead — the hook now divides the increment rather
  than declining it;
- take the new `phi_int` argument and form alpha as
  `ustatev(84) = ustatev(84) + k_alpha*phi_int`, **not** from
  `dt*phi_tot(g_new)` — the reason is in the `phi_int` section above, and it
  matters more here than anywhere else, since this is the deck that actually
  runs `n_sub = 1000`;
- **keep item 2 unchanged.** It guards a different failure (a bad state
  returned at a `dt`/`theta` combination the sweep did not probe) and
  sub-stepping does not make it redundant. Only the blanket refusal goes.

Because that source is not here, none of this is covered by the test suite or
by the syntax checks — it is the only part of the change that arrives at the
machine unverified, so if step 2 fails, check this first.

## If it fails

The likely failures, in order:

- **`ierr = 6`** — the server answered without `phi_int`. An old
  `material_server.py` is running; restart it from this checkout.
- **Slow to the point of unusable** — expected risk. Report the wall time
  rather than tuning `DTMAX_ECO` upward to make it finish: that threshold is
  measured, and raising it to buy speed reintroduces exactly the silent
  divergence the guard exists for.
- **alpha differs in run 1** — do not proceed. `n_sub` should be 1 there.
