# Biofilm growth + viscoelastic material law — for the ANSYS/UPF framework

This is the constitutive law from the LUH/IKM biofilm thesis work, packaged to
be callable from your framework at the point where `Usermat_P21-V21_*.F`
currently calls `AceGenNeoHookV04`.

Two files, no dependencies:

| File | |
|---|---|
| `biofilm_material_v01.f` | `BIOFILM_GROWTH_VISCO_V01` — the routine to call |
| `biofilm_stress_core.f` | the verified core it calls |

No ANSYS includes, no `usercm` common block, no UPF calls — the same
self-contained shape as `AceGenNeoHookV04`, which means it does not care which
ANSYS release it is linked into. Only the surrounding `usermat` glue is
release-specific, and that is your file, not this one.

`wrapper_driver.f` is a small stdin/stdout driver for testing the routine
outside ANSYS. `solution_loop.png` shows where it attaches.

## The call

```fortran
      subroutine BIOFILM_GROWTH_VISCO_V01(
     &   mDefGrad, vCauchy, mTangCC,
     &   sYoung, sYoungL, sNu, sNuL, sBiofilm,
     &   sGrowth, mFvN, mFvN1,
     &   sEta, sDt, sC01Ratio, sMtype,
     &   sElasticWork, sKeyCut, sID)
```

Deformation gradient and material constants in; Cauchy stress and consistent
tangent out. `(E, ν)` per phase blended by `sBiofilm`, as your deck already
carries them.

What it adds over `AceGenNeoHookV04`, which is purely elastic:

- growth kinematics `F = Fe·Fv·Fg` with `Fg = (1+α)I`
- a viscous branch carrying `Fv`

Voigt order is ANSYS (`11,22,33,12,23,13`).

## Three things to know before wiring it in

**1. The viscous state needs nine slots, not six.** `Vdp_Cv_n` reserves six,
which would be right if `Fv` stayed symmetric. It does not: measured at
~6e-5 relative asymmetry after 20 steps over 200 random states, so it cannot
be reconstructed from a six-component Cauchy–Green tensor without losing
information. `mFvN` / `mFvN1` are full 3×3.

**2. `sGrowth` is the growth variable, and is deliberately not called
`sAlpha`.** In `AceGenNeoHookV04` the argument named `sAlpha` is fed
`Sdp_sumLocal`, a local biofilm average, and is not a growth variable at all.
Reusing that name here would invite exactly the wrong wiring.

**3. `sDt` must resolve the viscous relaxation time `η/(2·C10)` — and the
routine enforces this rather than trusting you to read it.** The flow increment
is evaluated at the old state, so accuracy degrades as `Δt` approaches that
time; past `Δt/τ ≈ 0.5` the stress changes sign and past `Δt/τ ≈ 1` it
diverges. Because a wrong step would otherwise return a plausible-looking wrong
stress, `Δt/τ > 0.5` sets `sKeyCut = 1` and returns without computing one, so
the solver cuts the increment and retries. **You may therefore see cut-backs
you did not expect** — that is this check, not an instability.

Two things follow. Growth shrinks `τ` as `C10` rises, so a step that is fine
early in a solve can stop being fine later; it is checked every call for that
reason. And `sEta = 0` selects the purely elastic path, which has no relaxation
time and is never restricted.

The threshold is the `DTMAX_RATIO` parameter at the top of
`biofilm_material_v01.f`. It is set at the sign-flip rather than tighter,
deliberately: below it the answer loses accuracy but stays qualitatively right,
and trading accuracy against step size is your engineering choice, not ours to
force. Lower it if you want the routine to insist on more resolution.

On a cut-back (`sKeyCut = 1`) every output is defined: stress and tangent are
zeroed and `mFvN1` is returned unchanged, so reading them before checking the
flag gives zeros rather than whatever was in the work array.

## What has been verified

- The core is **0 ULP identical** to an independently written Abaqus UMAT over
  **8017 deformation states** — identity and growth, uniaxial, shear, elastic
  (`η=0`), frozen (`η→∞`), large growth, neo-Hookean and Mooney–Rivlin, and
  random finite-strain states. Those states were found by adversarial search
  for a disagreement, not by sampling.
- `BIOFILM_GROWTH_VISCO_V01` is deliberately a thin adapter around that core
  rather than a reimplementation, and this is tested at zero tolerance: fed
  `(E, ν)`, it must reproduce the core fed the `(C10, C01, D1)` those map to.
  It does, so the verification above applies to the routine as delivered.
- **This exact routine has since been exercised in a real ANSYS solve, not
  only gfortran unit tests (2026-09-02).** Built and linked into a custom
  ANSYS 2022 R2 (v222) `ANSYS.exe` and run through a small harness-only
  `usermat()` entry point that does nothing but unpack `ustatev`/`prop` and
  call `BIOFILM_GROWTH_VISCO_V01` — i.e. the same call shape your
  `Usermat_P21-V21_*.F` will use at your `AceGenNeoHookV04` site. Two cases:
  - a fully-constrained single element, matching the closed-form reference
    to display precision (`SX=SY=SZ` exact to 5 significant figures, shear
    exactly zero to machine noise);
  - the same two-layer curved-shell geometry already characterised for this
    law (12240 elements), re-run unchanged except for the material block —
    SEQV min/max/mean identical to the digits printed against the original,
    already-verified build across all 12240 elements.

  Both runs: 0 errors. This is independent of, and additional to, the
  gfortran-level testing above — it confirms the `(E,ν)→(C10,C01,D1)`
  conversion, the growth kinematics, and the core all thread correctly
  through this exact routine inside a real ANSYS solver process. Full
  narrative and raw output:
  [`../ansys_usermat/apdl/V222_PORT_INSTRUCTIONS.md`](../ansys_usermat/apdl/V222_PORT_INSTRUCTIONS.md)
  §1.6, and the executable walkthrough
  [`../ansys_usermat/apdl/v222_wrapper_verification.ipynb`](../ansys_usermat/apdl/v222_wrapper_verification.ipynb).

One honest caveat, since it is better heard from us than found later.
Cross-implementation agreement establishes that two ports compute the same
expression identically, not that the expression is right. We have since found
that the isochoric split in both applies `J^(-2/3)` to the subtracted trace but
not to the tensor beside it. The resulting error is **purely spherical** — a
pressure error; deviatoric and von Mises stress are unaffected to machine
precision — and it reaches the deviator only weakly through the viscous flow
driver (≤0.01% in von Mises over 40 steps). It is documented rather than
corrected for now, because correcting it invalidates the reference values and
it is not yet settled whether the published formulation defines the split as
implemented. Happy to discuss which way you would prefer it.

## Testing it outside ANSYS

`wrapper_driver.f` reads a state on stdin and prints stress, `Fv`, cut-back
flag and tangent:

```
gfortran -ffixed-line-length-132 wrapper_driver.f biofilm_material_v01.f \
    biofilm_stress_core.f -o wrap

# F (3x3), then Fv_n (3x3), then
# Young YoungL Nu NuL Biofilm Growth Eta Dt C01Ratio Mtype
echo "1.06 0.02 0.0  0.0 0.97 0.01  0.0 0.0 0.98
1 0 0  0 1 0  0 0 1
1000 1 0.30 0.30 1.0 0.2 5.0 0.0001 0.15 1" | ./wrap
```

## Provenance

Generated from the `pde-fem-biofilm` repository by `handover/make_handover.py`
— extracted from the sources under test rather than copied by hand, so it
cannot drift from them. Please send changes back there rather than editing
these files, so the verification still applies.
