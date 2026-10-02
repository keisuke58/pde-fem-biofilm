# Reproducing Klempt et al. (2024): what was found, and what it means

Status 2026-09-30. The 2024 PDE is recorded as **not independently reproduced**,
and this is the account of why. The paper is bundled at the repository root
(`references/Klempt2024_Hamilton_biofilm_growth_BMMB.pdf`, CC BY 4.0 — see [`THIRD_PARTY.md`](THIRD_PARTY.md)),
so everything below is read off it rather than inferred.

Code: [`JAXFEM/klempt2024_quantitative.py`](JAXFEM/klempt2024_quantitative.py)
(the reproduction) and
[`JAXFEM/klempt2024_sensitivity.py`](JAXFEM/klempt2024_sensitivity.py)
(the one-change-at-a-time study). Numbers in
`JAXFEM/klempt2024_results/`.

---

## 1. The headline

**The two variants that fit the paper's figures best are the two the paper
does not use, and the combination the paper does use is the worst fit.**

That is the finding. It is not a small numerical discrepancy to be tuned away.

| | growth term | consumption | Fig. 7 "high", φ at T* = 0.20 |
|---|---|---|---|
| the paper's own equations | `∇φ · n_∇c` | `g φ` | **0.105** |
| what fits best here | `\|∇φ · n_∇c\|` | `g φ c` | 0.886 after also moving `k` |
| the paper's figure | | | **1.000** |

## 2. What the paper actually says

**Growth — Eq. 34.**

```
φ̇ − β ∇²φ − k_α α + ‖∇φ‖ · (r c)/(k + c) · n_∇φ · n_∇c = 0
```

A plain dot product, no absolute value. Section 4.1 depends on it being signed:

> "the gradient of biofilm and the gradient of nutrients are almost
> perpendicular to each other resulting in a small value for the vector product
> and consequently in minimal to no growth"

That suppression on the perpendicular faces is what produces the egg shape in
Fig. 3. So taking the absolute value does not read the paper differently — it
removes the mechanism the paper describes. It is worth 3.2× here, and it is not
available.

**Consumption — Eq. 35, fixed by Eq. 24.**

```
ċ − d ∇²c + g φ = 0        g*_c = ḡ φ
```

Eq. 24 calls this "the simplest possible functional dependency, a linear
relation". Zeroth order in `c`. So first-order consumption is not the paper
either.

**Consequence.** The earlier report that Fig. 4 came "within 0.17 / 0.10"
is the agreement of a variant the paper does not use, and must not be quoted as
agreement with Klempt 2024.

## 3. Two suspicions cleared

**The constants are right.** Table 2's values are *post-division*. The paper
sets `η_φ = η_c = 10⁻¹⁰` and states that "parameters which have been divided by
their respective η will lose their bar", so `d = 10¹⁰`, `β = 2`, `k_α = 10⁻³`,
`k = 1`, `g = 10⁸`, `r = 100` are already the coefficients of Eqs. 34–36. The
`η_φ` visible in Table 1's weak form is the same quantity before division, not a
factor missing from the implementation.

**The seed and the nutrient are not the cause**, measured rather than assumed
(`klempt2024_sensitivity.py`): a six-node initial disk reaches only 0.238 against
the paper's 1.000 and scales with the seed rather than changing the rate, and
weakening consumption a hundredfold moves 0.105 to 0.108. The colony is not
short of nutrient; it is short of growth.

## 4. An inconsistency in the paper

Table 1 step 3a solves

```
(α_{n+1} − α_n) / ((1 + α_{n+1}) Δt)  −  (k_α / α_n)(φ_{n+1} − φ_n) / Δt  =  0
```

driven by `φ̇`, while Eq. 36 is `α̇ = k_α φ`, driven by `φ`. Which was run is not
recoverable from the text. It cannot be the factor of ten — `k_α α` is of order
10⁻³ against a growth term of order 50 — but "the paper's α equation" is
ambiguous, and any future comparison has to say which one it used.

## 5. What is still unseparated

- **Scheme.** Explicit upwind finite differences here against the paper's
  implicit Galerkin FEM, which bisects to as many as 10⁶ substeps. The growth
  term is proportional to `‖∇φ‖`, so whatever under-resolves the interface
  under-states the growth directly.
- **Geometry.** A fixed grid here; there, `Fg` swells the domain the averages
  are taken over (Fig. 13 shows it, at 500× exaggeration).
- **The clip.** `φ` is clipped to [0, 1] here. The paper has no such clip; its
  energy is non-convex and localisation of `φ` is intended.

## 6. What this is worth

Two things, and they are not the same.

**For the thesis.** The 2024 PDE is a *reference*, not a dependency: this
repository's own ecology follows the later multi-species paper, whose fourteen
numerical examples do reproduce (`JAXFEM/klempt2026_reproduction.py`). The
honest statement is that the 2024 single-species PDE was not independently
reproduced, that the discrepancy was localised to two named terms, and that
those terms are the paper's own — which is a stronger position than a silent
omission.

**For the meeting.** Meisam Soleimani is a co-author of this paper. The question
"which growth and consumption terms produced Figs. 4 and 7, given Eq. 34 and
Eq. 24 as printed, and which α equation — Eq. 36 or Table 1 step 3a?" is one
only the authors can answer, and it is worth more asked than worked around.

## 7. Also worth knowing, from the same reading

- **Klempt 2024 reports hydrostatic stress `p`, not von Mises.** This
  repository's known, unfixed isochoric-split defect
  ([`DEVIATOR_SCALING_FINDING.md`](DEVIATOR_SCALING_FINDING.md)) is *purely
  spherical*: harmless for the von Mises this repository reports, and directly
  in the way the moment anything is compared against the paper's stress figures.
- **The paper declares its own absolute stresses uncalibrated.** Fig. 8's
  caption: "The absolute values are dependent on Young's modulus which needs to
  be adjusted to fit real biofilm in future works." That is a precedent for this
  work's `E_SPEC` limitation, from the model's own authors.
- **The paper's validation is qualitative morphology**, not measured stress —
  mushroom then droplet shapes, "in accordance with observations in in vitro
  experiments" (Toyofuku et al. 2016).
- **It contains five numerical examples**, of which two were attempted here:
  §4.1 directional growth, §4.2 biofilm on nutrients, §4.3 growth against rigid
  obstacles (geometry from Albero et al. 2014, later used by Soleimani 2019 and
  Soleimani et al. 2020), §4.4 maze, §4.5 grate.

## 8. The partner element's front term, |∇²φ| instead of |∇φ| (2026-10-02)

The partner's element multiplies the front-growth term by `|∇²φ|` where
Eq. 34 has `|∇φ|`. In the partner element, the Fig. 7 setup with Table 2
values scaled to its deck (IKMHIWI03, commit 4a54c09) does not grow: the seed's
peak φ falls from 1.0 to 0.1. Was the `|∇²φ|` the cause? I tested this here,
in the same 21³ reproduction: `klempt2024_quantitative.run(...,
growth="lap")`, Fig. 7 high nutrient, consumption as printed, dt = 1e−3.

| T* | paper (mean φ) | Eq. 34 as printed, `\|∇φ\|` (mean / max φ) | `\|∇²φ\|` (mean / max φ) |
|---|---|---|---|
| 0.1 | 0.55 | 0.027 / 0.93 | 0.011 / 1.00 |
| 0.2 | 1.00 | 0.068 / 0.60 | 0.032 / 1.00 |
| 0.5 | 1.00 | 0.160 / **0.36** | 0.127 / **1.00** |

- **It is the printed `|∇φ|` form that erodes the peak** (1.0 → 0.36): it is a
  transport of φ towards the nutrient, as §2 explains, so the back of the
  colony empties. With `|∇²φ|` the peak stays at 1.0.
- **Neither form fills the cube** as Fig. 7 does by T* = 0.2. That is the
  negative result of §1, not something new.

So this test does not support `|∇²φ|` as the cause of the partner run's
decay. The decay looks more like what Eq. 34 as printed does on its own. Two
caveats:
- `|∇²φ|` scales with 1/h², so its weight depends on the mesh (here h = 1 µm;
  the partner mesh is 0.25 mm with deck-scaled parameters);
- the partner element also has its own discretisation.

Not a proof either way; it says which test to run next on the machine: the
same deck with the front term switched to `|∇φ|`.

**Correction, same day (IKMHIWI03, f278885): the reason for the decay was
neither form.** In the partner element the front term was **switched off**.
`USolBeg` reads `ORI_WEIGHT12` and `ORI_WEIGHT22` into species 1's weight
`sGdp_OriWeight11`, so the last value read, `ORI_WEIGHT22 = 0`, zeroes it. The
front term never acted, and the seed only diffused. My reading above, that the
decay looked like Eq. 34 as printed, was wrong: it was diffusion with no front
term at all.

On a test copy with that slip fixed:
- with `|∇²φ|` the explicit update is unstable at dt = 0.005;
- with `|∇φ|` (Eq. 34) at dt = 0.001 the seed spreads, from 4 to 94 traced
  points with φ ≥ 0.5 by T* = 0.18. That is the paper's pace (cube filled by
  T* = 0.2). But φ exceeds 1 (up to about 3) and the run breaks down near
  T* = 0.22: the element does not bound φ.

This reproduction clips φ to [0, 1] and stays far slower than the paper, so
the two now differ in a way worth separating later.

**Why φ exceeds 1 there: the cell Péclet number, not the time step.**
Estimated from `klempt2024_cases.deck_values` (2 mm cube, 0.25 mm elements):
- front speed v = r·c/(K+c) ≤ 10·½ = 5 mm/T*, β = 0.02 mm²/T*;
- CFL v·dt/h = 0.02 at dt = 0.001, far inside any limit;
- **cell Péclet number Pe = v·h/(2β) ≈ 31**.

The front term is an advection of φ up the nutrient gradient (§2). Advection
by central differences at Pe ≫ 1 gives over- and undershoots, which matches
φ reaching about 3. The `|∇²φ|` form scales with 1/h², which matches the
explicit instability at a front rate of about 560/T*.

Remedies, each a change to the partner's code, so for Oliver to decide:
- upwind the front term, as this reproduction does (bounded);
- clip φ to [0, 1] (simple, but not mass-conserving);
- refine to Pe < 1, which needs h < 0.008 mm (impractical).



## 9. Test case 1 (Fig. 3/4) revisited with the paper's set-up (2026-10-02)

`JAXFEM/klempt2024_case1_bc.py` (results in `klempt2024_results/case1_bc.json`;
the full table is in its docstring). The earlier runs held the nutrient on one
edge line; the paper's Fig. 2 draws a strip along the edge. With the strip the
early part of Fig. 4 follows. The late part does not with Eq. 34/35 as
printed, under any numerical variant (seed held, no clipping, 2nd-order ENO):
Eq. 34 as printed advects the colony rather than growing it. With growth on
both faces (|∇φ·n_∇c|, as the paper's text describes) and first-order
consumption, mean φ matches Fig. 4 to 0.05 over T* ∈ [0, 1], but mean c is
about twice the paper's. φ and c together need a stronger consumption and a
smaller Monod constant than Table 2 (diagnostic: g ×4, k = 0.25 gives 0.07 /
0.09), the same direction as the Fig. 7 finding above. So: φ reproduced with
two stated departures from the printed equations; φ and c together not
reproduced with Table 2. Questions for the authors: which form of Eq. 34 and
35 was run (Table 1 differs from both), the width of the nutrient strip, and
whether c in Fig. 4 is normalised as in Eq. 35.

## 10. The figures read by colour, and test case 4.2 (Fig. 7) (2026-10-02)

`JAXFEM/klempt2024_digitize.py` reads Fig. 4 and Fig. 7 off the PDF by colour,
calibrated on the tick marks (`klempt2024_results/paper_curves_digitized.json`).
Fig. 4 and three of Fig. 7's curves agree with the earlier eye readings within
0.02; Fig. 7 "high" φ was read far too low early (0.12 / 0.55 at T* = 0.05 /
0.10; the figure gives 0.35 / 0.73). Results quoted against `PAPER["fig7_high"]`
in `klempt2024_quantitative.py` / `klempt2024_sensitivity.py` used the eye
reading.

`JAXFEM/klempt2024_case2.py` runs 4.2 with the paper's set-up (bottom face
c = 1, 5 µm disk above it) against the digitised curves. Not reproduced under
any reading tried, and two numbers say why: the "high" case fills the cube at a
front speed near 150 µm/T*, beyond Table 2's r·c/(k+c) ≤ 50 µm/T*; and Eq. 34
as printed moves the biofilm towards the nutrient (down), while Table 4 shows
it growing away from it (up). The nutrient plateau of the "high" case does
hold, and is a clean check: first-order consumption with Table 2's d and g
gives 0.482 (tanh 2 / 2), the paper 0.489, the printed zero-order form 0.139.
With §9 this points to first-order consumption in the paper's computations.
Questions for the authors (add to §9): the growth term actually run (sign,
both faces?), r and k, and the size of the initial biofilm (the paper's mean
φ at T* = 0.01 is ten times a one-layer disk).
