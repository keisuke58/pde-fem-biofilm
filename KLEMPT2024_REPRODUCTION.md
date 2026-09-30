# Reproducing Klempt et al. (2024): what was found, and what it means

Status 2026-09-30. The 2024 PDE is recorded as **not independently reproduced**,
and this is the account of why. The paper is bundled at the repository root
(`felix_s10237-024-01883-x.pdf`, CC BY 4.0 — see [`THIRD_PARTY.md`](THIRD_PARTY.md)),
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
