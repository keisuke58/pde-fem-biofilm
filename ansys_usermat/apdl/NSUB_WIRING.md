# Sub-stepping the ecology: what landed, and what is still open

2026-09-29. The ecology ODE is only trusted below `dt = 1e-4` (measured — a
Python sweep turns erratic from about `1.5e-3`, and not monotonically, which is
why a guard exists at all; `V222_PORT_INSTRUCTIONS.md` §8 has the failure). A
mechanical deck's step is orders of magnitude coarser, and ANSYS's own
bisection does not reach a 1000× reduction, so the guard refused every step and
the wired deck could not complete.

The increment is now divided **inside one material call**: one socket round
trip, `n_sub` cheap ODE steps within it. Not a new scheme — every reference
under `apdl/` already advances the ODE by chaining `ecology_step`; the bridge
does the same thing, so it cannot drift from the references it is verified
against.

## Where it stands

**Done, and confirmed on real hardware** (`bdddf8a`). `usermat_biofilm.f`
splits `dTime` into `ceiling(dTime/DT_ECO_MAX)` sub-steps with
`DT_ECO_MAX = 1e-4` and passes that as `NSUB` to `biofilm_ecology_hook`; the C
shim and the ISO_C_BINDING interface carry `n_sub` (by VALUE) and `phi_int`
(out). The evidence is `t_growth_ecology_substep.dat` with its output kept
beside it in `out_ecology_substep.txt`: ANSYS v222, custom `ANSYS.exe`, 0
errors, `alpha = 1.2908e-2`, `SX = SY = SZ = -7.6033`, and all twelve ecology
components equal to an independent Python reference at printed precision.

That deck earns its place by telling the two builds apart rather than merely
passing: the pre-`n_sub` build gives `alpha = 1.4914e-2` and `gamma = 150.70`
against `295.16`, so a silent collapse to `NSUB = 1` cannot slip through as a
pass.

`phi_int` is the part that was easy to get wrong. Growth accumulates as
`k_alpha * integral(phi_tot dt)`. With one step, `dTime * phi_tot(g_new)` *is*
that integral; past one sub-step it is not, because the end state does not
stand for the whole increment. The server returns the sum over sub-steps and
the material forms alpha from that. At `dTime <= DT_ECO_MAX` the two agree
exactly, so every existing deck and test gives the same answer as before.

## Build order — the trap that cost a run

**`usermat_py_hook.f` must compile before `usermat_biofilm.f`.** In the other
order ifort reads the stale `biofilm_py_bridge.mod` left by the previous build,
fails on the argument count, and `link_v222.ps1` still links the *old*
`usermat_biofilm.obj` into a mismatched `ANSYS.exe` — a binary that runs and
produces plausible numbers from the previous wiring. Delete the stale
`ANSYS.exe` / `.lib` / `.exp` / `.map` before every relink as well, or
`ansys.lrf`'s `*.lib` wildcard collides with the new output.

`link_v222.ps1` drives ifort only, so the C shim is compiled separately:

```
cl /c /O2 ansys_usermat\coupling\biofilm_py_eval.c
```

## Still open

**1. The partner's own deck, and the call-site edit it needs.** This
repository's `usermat_biofilm.f` is done. The partner's NEM material routine
has a *second* copy of the ecology call site, its source lives only on the
machine (`F:\biofilm_upf_wired`) and is deliberately not committed here, and
per `V222_PORT_INSTRUCTIONS.md` §8 it still:

1. **refuses outright** when `dt > DT_ECO_MAX`, and
2. sanity-checks the returned state, holding alpha at `ustatev(84)` and asking
   for a cut-back (`keycut=1`) on either failure.

Item 1 is what stopped that deck on 09-07, and sub-stepping does not remove it:
the refusal happens *before* the hook is reached, so with it in place the deck
fails exactly as it did and none of the new code runs. Three lines change
there, mirroring `usermat_biofilm.f` — compute `NSUB = ceiling(dt/DT_ECO_MAX)`
and pass it instead of declining; form alpha as
`ustatev(84) + k_alpha*phi_int`, **not** from `dt*phi_tot(g_new)`, which
matters more here than anywhere else since this is the deck that actually runs
`NSUB = 1000`; and **keep item 2 unchanged**, because it guards a different
failure (a bad state at a `dt`/`theta` combination the sweep never probed) that
sub-stepping does not make redundant.

Because that source is not here, none of it is covered by the test suite or the
syntax checks. It is the one part that arrives at the machine unverified.

What to record when that deck runs: `NUMBER OF ERROR MESSAGES`; whether it
reaches the end of the load step; and `SVAR(10)` and `SVAR(84)` on a few
elements, since a completed run with alpha stuck at its seed means the hook was
never reached.

Speed is no longer the thing to watch. The server runs the chain as one
compiled scan, bit-identical to chaining `ecology_step`, which took 1000
sub-steps from ~10 s to ~0.07 s per call: the single-element deck at
`dTime = 0.1` went 357 s to 5.6 s, and the 54-element cylinder 131 s to 13.6 s.

**2. The scale of `k_alpha`, which sub-stepping exposed rather than fixed.**
This is now the real blocker for any deck stepping at `TIME INC = 0.1`, and it
is a modelling decision, not a numerical one. With 1000 sub-steps the ecology
state comes back finite and sane (`gamma` 111–326, phi summing to 1) and alpha
comes back **large**: 4.72 / 4.46 / 0.39 / 4.64 for CS / CH / DS / DH. With
`Fg = (1+alpha)I` that is a 5.7× stretch per direction, so the "element highly
distorted" error on that run is **the growth itself, not a numerical failure**.

`k_alpha = 50` was chosen for decks that run 0.1–1 ms of ecology time in total.
It has no meaning on a deck whose pseudo-time unit is something else, so
`k_alpha` — or the ecology/mechanical time-unit map — has to be set on that
deck's own time axis before any result from it means anything. That needs the
partner's answer to what one unit of their `TIME` physically stands for; it is
question 4 on the 10/1 agenda. Re-running before then only reproduces the same
5.7× stretch.

Either way, do not raise `DT_ECO_MAX` to make a run finish: that threshold is
measured, and raising it to buy speed reintroduces exactly the silent
divergence the guard exists for.

**3. Regression, if the custom exe is rebuilt.**
`t_growth_cylinder_ecology_4region_all_real_clsm.dat` runs at `dt = 1e-5`,
below `DT_ECO_MAX`, so `NSUB = 1` and it must be unaffected: 0 errors and the
four regions' alpha unchanged at
`4.6145e-3 / 4.6920e-3 / 4.3557e-3 / 4.8186e-3`. **If those move, stop** —
something other than sub-stepping changed, and the build-order trap above is
the first thing to check.

**4. The tooth/implant geometry**, the roadmap's actual open item and a larger
job.

## If something fails

- **`ierr = 6`** — the server answered without `phi_int`. An old
  `material_server.py` is running; restart it from this checkout.
- **Plausible-but-wrong numbers** — suspect the build order above before
  suspecting the physics. A mismatched `ANSYS.exe` runs cleanly.
- **Alpha far too large, elements distorting** — not a numerical failure; see
  the `k_alpha` item above. Do not raise `DT_ECO_MAX` to work around it.
