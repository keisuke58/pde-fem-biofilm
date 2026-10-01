# One-species coupling into the partner's element

Decided 2026-10-01 in the meeting with Oliver: couple **one species fully**
first, then extend to two. Main results follow Klempt et al. (2024); nothing
is added to their formulation.

## Why one species first

Every obstacle that made the five-species merge "future work" goes away:

| obstacle with five species | with one species |
|---|---|
| no shared variable space (their 2 species / 2 nutrients, our 5 species) | one `φ` on both sides |
| who computes `φ` — two ecology models would both evolve it | their NEM solves `φ` and `c`; we only grow and stress |
| nothing to verify against | the chain reduces to closed form, and `JAXFEM/klempt2024_quantitative.py` is single-species Klempt 2024 |
| a Python server per Gauss point | none: Eq. 36 is one line of Fortran |

## The law

Klempt et al. (2024) Eq. 36, read off the bundled PDF:

    α̇ − k_α φ = 0

Their `Fg = α_K I` with `α_K(0) = 1`; ours `Fg = (1+α) I` with `α(0) = 0`
(Soleimani, Haverich & Wriggers 2021, Eq. 15). `α = α_K − 1`, a constant
shift, so it is the same equation.

## What is in this repository (verified, no ANSYS needed)

- `ansys_usermat/growth_from_phi.f` — `BIOFILM_ALPHA_FROM_PHI`, explicit:
  `α_{n+1} = α_n + k_α φ Δt`. No cap, no clamp, no gate.
- `ansys_usermat/one_species_reference.py` — the same, in Python, written
  independently; also predicts `SVAR(84)` from a printed `φ` history.
- `ansys_usermat/crosscheck/one_species_driver.f` — the chain at one Gauss
  point: `φ → BIOFILM_ALPHA_FROM_PHI → BIOFILM_GROWTH_VISCO_V01 → σ`, with
  `η = 0`, `C01 = 0`, neo-Hookean.
- `tests/test_one_species_coupling.py` — 7 tests: Eq. 36 for constant `φ`;
  bit-identical to the Python reference on a varying `φ`; no biofilm, no
  growth; nothing added to Eq. 36; free growth gives zero stress; constrained
  growth matches closed form (plus the documented spherical term,
  `DEVIATOR_SCALING_FINDING.md`).

## The call-site change (apply on IKMHIWI03, `F:\biofilm_upf_wired`)

The partner's file stays off this public repository; this is the recipe. It
sits inside the existing `prop(1) /= 0` branch (§7 of
`V222_PORT_INSTRUCTIONS.md`), next to the `prop(6) = 1` ecology mode, as a new
`prop(6) = 2` mode:

```fortran
      else if (prop(6) .gt. 1.5d0) then
C       one-species coupling: the partner's local phi drives growth,
C       Klempt 2024 Eq. 36. ustatev(84) is the CONVERGED alpha -- ANSYS
C       hands it back unchanged on every equilibrium iteration, so growth
C       is counted once per increment, not once per iteration.
        PHI_LOC = <local biofilm phi from the pool -- see question 1>
        call BIOFILM_ALPHA_FROM_PHI(ustatev(84), PHI_LOC, prop(7),
     &                              dTime, ALPHA_NEW)
        sGrowth     = ALPHA_NEW
        ustatev(84) = ALPHA_NEW
```

Declare `double precision PHI_LOC, ALPHA_NEW` explicitly. Material constants
for the Klempt setting: `prop(2) = 0` (η), `prop(3) = 0` (C01 ratio),
`prop(4) = 0` (neo-Hookean). Add `growth_from_phi.f` to the compile list; it
has no module dependency, so build order does not matter for it. No material
server is needed in this mode.

## Questions that must be answered first (one answered)

1. **Which pool variable is the local `φ` at the Gauss point** —
   `Sdp_bio1_n` or `Sdp_locbio1_n`? Their file uses both (`Sdp_sumBio` from
   the first, `Sdp_sumLocal` from the second), and the code does not say which
   is the Gauss-point value.
2. ~~**`Sdp_sumBio = Sdp_bio1_n + Sdp_bio1_n`.**~~ **Answered 2026-10-01:
   Oliver thinks it is a typo** for `Sdp_bio1_n + Sdp_bio2_n`. Until his
   source is corrected, in this mode pass `PHI_LOC` as `sBiofilm` as well,
   so the stiffness blend sees `φ` and not `2 φ`. For the two-species step
   the corrected sum is exactly the `φ₁ + φ₂` that drives growth.
3. **The physical time unit of the deck**, which sets `k_α` (`prop(7)`).
   `k_α = 50` was chosen for millisecond decks; at `TIME INC = 0.1` it takes
   `α` to about 4.7 per increment.

## Verifying the first run

1. Record `NUMBER OF ERROR MESSAGES`, whether the load step reaches its end,
   and `SVAR(84)` plus the local `φ` at a few Gauss points over time.
2. Feed that `φ` history to `one_species_reference.alpha_history`: it must
   reproduce `SVAR(84)` to printed precision. A mismatch means the call site
   reads a different `φ` than it prints, or counts growth per iteration.
3. Then, on a matching problem, compare the `α` field against
   `JAXFEM/klempt2024_quantitative.py`. That is the first cross-check between
   the two codes that has a shared variable to compare.

## Then: two species

The same routine with `PHI_LOC = φ₁ + φ₂` — total biofilm drives growth. That
needs question 2 settled, because the sum is exactly what `Sdp_sumBio` was
meant to be.
