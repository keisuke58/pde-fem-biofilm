# Soleimani et al. (2023), multi-species co-aggregation — what it settles

> Soleimani M, Szafranski SP, Qu T, Mukherjee R, Stiesch M, Wriggers P, Junker P.
> **Numerical and experimental investigation of multi-species bacterial
> co-aggregation.** *Scientific Reports* **13**:11839 (2023).
> doi:[10.1038/s41598-023-38806-2](https://doi.org/10.1038/s41598-023-38806-2)

Read 2026-09-30. This is the supervisor's own multi-species paper, with the
dental clinic side (Stiesch, MHH) and both Klempt co-authors on it, and it was
not cited anywhere in this repository. It answers several questions that were
open here.

---

## 1. It is a different formulation, and that matters for our headline result

| | Soleimani 2023 | this repository |
|---|---|---|
| species variables | φ₁, φ₂ **independent phase fields**, each ∈ [0,1] | φᵢ **volume fractions on a simplex**, Σφᵢ + φ₀ = 1 |
| equation | Allen–Cahn, `□̇ = −f′(□) + ε²∇²□ + S(□,c)` (Eq. 1) | Hamilton-principle residual, Newton-solved |
| both species at one point | **yes** — that is the co-aggregation zone | no: they share one budget |

This bears directly on the finding that the four clinical conditions come out
within 0.4%. That collapse happens because `φ_tot = Σφᵢψᵢ` is one total and the
CLSM compositions are fractions summing to 1, so the differences cancel **in
the sum**. Soleimani 2023 has no such constraint: φ₁ and φ₂ are independent
fields that can both be 1 in the same place.

**So the degeneracy is a consequence of the simplex, not of multi-species
modelling.** That is a much more precise statement than "composition does not
reach the stress", and it is the one the thesis should make.

## 2. A route to make composition matter, from the same group

Eqs. 3 and 5 are

```
S₁ = R_s1 · H(φ₁ − φ_cri) · (1 + τ α) · φ₁ · c
S₂ = R_s2 · H(φ₂ − φ_cri) · (1 + τ α) · φ₂ · c
```

with `α` the co-aggregation field and `τ = 1` (Table 1). Growth is **boosted
where two species co-locate**. That is a pairwise interaction that survives
into the growth rate, which is exactly what summing into `φ_tot` destroys here.
It is the mechanism this repository's model is missing, described by the
supervisor.

## 3. Consumption: the same group uses first order in c

Eq. 7:

```
S_c(φ₁, φ₂, c) = −R_c (φ₁ + φ₂) c
```

**First order in `c`.** Klempt 2024's Eq. 35 is `ċ − d∇²c + gφ = 0`, zeroth
order, and Eq. 24 there calls the linear-in-φ form "the simplest possible
functional dependency".

`klempt2024_quantitative.py` found that first-order consumption fits Fig. 4
far better than the printed zeroth-order form, and that was recorded here as a
deviation from the paper. It now looks less like a deviation and more like the
form the group actually uses — see
[`KLEMPT2024_REPRODUCTION.md`](KLEMPT2024_REPRODUCTION.md), which should be read
with this in mind. It does not rescue the reproduction (the growth term is a
separate and larger problem), but it changes which variant is the odd one.

## 4. The α source is the same advection idiom as Klempt 2024's growth term

Eq. 6:

```
S_α = R_α · H(φ₁ − φ_cri) H(φ₂ − φ_cri) · (∇α · ∇c)/|∇c|
```

described in the text as `v · ∇α` with `v = R_α ∇c/|∇c|`. Structurally the
same as Klempt 2024's `r c/(k+c) ∇φ · n_∇c`: a signed advection along the
nutrient-gradient direction. Signed, again — no absolute value.

## 5. Mesh: their own statement that phase-field needs a fine grid

Table 1 gives a mesh size of **0.5 μm** on a 100×100 μm 2D domain, and the text
says why:

> "Due to the nature of phase-field equations, one needs to use sufficiently
> fine mesh. Even in two-dimensional examples, a mesh size of 0.5 μm leads to
> approximately 2×10⁴ nodes ... for a three-dimensional model with the same mesh
> size, the total number of DoFs exceeds 10⁷ and is computationally
> prohibitive."

**Klempt 2024 runs 1 μm in 3D** — coarser than this group's own stated
requirement, and the quote says why: 3D at their preferred size is not
affordable. That is direct support for the resolution hypothesis in
`JAXFEM/klempt2024_resolution.py`, and it is a statement by the authors rather
than an inference.

## 6. It answers the UserElement question THESIS_ASSIGNMENT.md left open

`THESIS_ASSIGNMENT.md` §4.1 asks whether the spatial field should be (a) extra
DOFs inside a UserElement or (b) precomputed and passed in per integration
point, and says not to implement until it is settled.

This paper does **(a)**: "a new multi-field user element is developed and
implemented in ANSYS", with **four scalar DOFs per node** (φ₁, φ₂, α, c),
weak forms in Eqs. 8–11, backward Euler in Eqs. 12–15, linearised with AceGen
and tailored to a FORTRAN user element. That is the supervisor's own answer.

## 7. The validation bar in this group is qualitative

Worth knowing precisely, because this repository keeps asking what it can claim
without stress measurements.

- "in qualitatively good agreement with the experimental observations"
- "Experimental investigations also comply qualitatively with the mathematical model"
- "the numerical predictions show plausible qualitative compliance"
- For test case 3 the experiment was **ruled out entirely**: a directional
  nutrient supply cannot be done without convection destroying the gradient.

What is compared is **morphology** — dendritic branching under diffusion-limited
growth against compact colonies under diffusion-unlimited growth — not numbers.
No quantitative error metric appears.

## 8. A measurement for ψ exists, and this group does it

The experiments stain with **LIVE/DEAD BacLight** (SYTO9 green for
non-permeable cells, propidium iodide red for permeable, yellow for
co-localised) and image by CLSM.

This repository declines to substitute ψ because the workbook's "living cells"
ratio exceeds 1 and is not the model's ψ ∈ [0,1]. A LIVE/DEAD partition is a
different quantity and is bounded. Together with
`JAXFEM/psi_spread_sensitivity.py` — which puts the required per-species spread
at about 17% to recover the 9.1% once reported — that turns "we have no data
for ψ" into a specific, answerable request to the experimental side.

---

## 9. The basis for two-way step 2 (added 2026-10-04)

The PDF is now in `references/Soleimani2023_coaggregation_SciRep.pdf` (CC BY 4.0,
`THIRD_PARTY.md`). Read again for `ansys_usermat/two_way_step2.py`:

- **Summing Eqs. 3 and 5 gives a composition-weighted rate.** With both
  Heaviside factors equal to 1 and no co-aggregation (α = 0),
  `S₁ + S₂ = (χ₁ R_s1 + χ₂ R_s2) φ c` with `φ = φ₁ + φ₂` and `χᵢ = φᵢ/φ`. This is
  the same weighting as step 2's `r = Σ χᵢ rᵢ`. The mechanism differs: here the
  growth is a local reaction, in step 2 it is the front term of the Klempt 2024
  field. So the weighting itself has a published precedent from the
  supervisor's group; the front-term form is this work's.
- **No species-specific values.** Table 1 sets `R_s1 = R_s2 = 500`; the paper
  gives no rates that differ between species. Step 2's `d` stays a sensitivity
  parameter.
- **Eq. 4 as printed** gives `H = 1` for `φ ≤ φ_cri` and 0 above, which would
  switch growth off once a colony is established; the text ("the minimum level
  of bacterial density that prevents the bacterial colonies from vanishing")
  reads as the opposite. Noted, not used here.
- **Validation is qualitative** (colony shapes, aggregate size; Figs. 3, 6), and
  the conclusion says so. The averaged curves of Figs. 9-10 come from randomly
  seeded colonies, so they cannot be reproduced number for number.

## What to do with this

1. **Cite it.** It is the nearest precedent and it was missing.
2. **Sharpen the headline claim** to the simplex, per §1.
3. **Re-read `KLEMPT2024_REPRODUCTION.md` §2** in light of §3 here.
4. **Settle the UserElement question** with §6 rather than leaving it open.
5. **Ask the experimental side for LIVE/DEAD per species**, with the precision
   §8 and the ψ sweep together specify.
