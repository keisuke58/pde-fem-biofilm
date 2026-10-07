# Soleimani, Haverich & Wriggers (2021) — the same machinery, and two fixes we need

> Soleimani M, Haverich A, Wriggers P. **Mathematical Modeling and Numerical
> Simulation of Atherosclerosis Based on a Novel Surgeon's View.** *Archives of
> Computational Methods in Engineering* **28**:4263–4282 (2021).
> doi:[10.1007/s11831-021-09623-5](https://doi.org/10.1007/s11831-021-09623-5)

Read 2026-09-30. Atherosclerosis, not biofilm — but the same author applying
the same framework: `F = Fe·Fg`, a nutrient diffusion–reaction, and a
phase-field with a double well separating healthy from overgrowing tissue. Two
of its choices are direct fixes for problems open in this repository.

---

## 1. It settles `CITATION_AUDIT.md` F1c

Eq. 15:

> "If `Fg` is assumed to be isotropic, it can be characterized by a scalar `α`
> and the identity tensor `I`, as follows: **`Fg = (1 + α)I`**. Here, the scalar
> `α` is introduced to capture the overgrowth. **The initial value of `α` is set
> to be zero.**"

That is **exactly this repository's convention**, verbatim, from the same
author. F1c recorded `Fg = (1+α)I` as an undocumented deviation from Klempt
2024's `Fg = αI` and offered two explanations; explanation (a) is right and it
is not a deviation at all. Both conventions live in this group's own work:

| | | `α` at zero growth |
|---|---|---|
| Soleimani 2019 Eq. 3, Klempt 2024 | `Fg = α I`, `α := Jg^{1/3}` | `α = 1` |
| **Soleimani 2021 Eq. 15, this repository** | `Fg = (1 + α) I` | **`α = 0`** |

Identical, with `α_ours = α_Klempt − 1`. The thesis can cite Eq. 15 for its own
convention instead of explaining a shift.

## 2. It caps `α` — which is the fix for the `k_alpha` blocker

Eqs. 16 and 17:

```
Lg = Ḟg Fg⁻¹ = α̇/(1+α) · I                                   (16)
α̇/(1+α) = k_g · H(α − α_cri) · φ                              (17)
```

with `H` the Heaviside step, and the paper says plainly why it is there:

> "The Heaviside function `H(α − α_cri)` is introduced to **prevent the variable
> `α` from growing unboundedly**. ... Without the imposition of limiting
> constraints on the growth function, **it can literally approach infinity that
> is physically inadmissible**."

**That is our blocker.** On the partner's deck at `TIME INC = 0.1`, `α` reaches
≈ 4.7 — a 5.7× stretch per direction — and the element distorts. We recorded
that the distortion "is the growth itself, not a numerical failure", which is
right, and concluded the time unit had to be settled first. This paper says the
growth law should carry a cap regardless, and that an uncapped one is
physically inadmissible.

`usermat_biofilm.f` integrates `ALPHA = ALPHA + KALPHA*PHIINT`, uncapped, and
the only guard is `ALPHA < 0 → 0`. → `ansys_usermat/apdl/NSUB_WIRING.md`.

Note also the **rate form**: growth is written on `α̇/(1+α)`, i.e. on `tr(Lg)/3`,
not on `α̇`. That is the principled measure, and it explains Klempt 2024's
Table 1 step 3a, which `KLEMPT2024_REPRODUCTION.md` §4 flagged as inconsistent
with its own Eq. 36 — the `α̇/(1+α)` there is this, not a typo. Which driver
belongs on the right-hand side (`φ` here, `φ̇` there, `φ` in Eq. 36) remains the
open question.

## 3. It warns about the advection term we are struggling with

Eq. 13's phase-field driver is the group's recurring idiom, now seen in three
papers:

```
S(φ,c) = R_s · H(c − c_cri) · (∇φ · ∇c)/|∇c|
```

and the Remark after it is the useful part:

> "The authors are aware of the **torturous complexity** arising naturally from
> the advection term in form of `v·∇φ` ... An unsymmetric contribution to the
> stiffness matrix is one of the cumbersome consequences. Furthermore, **if the
> divergence of `v` is positive, it leads to the loss of coercivity ... and
> ultimately to stability issues**. If the advection term becomes dominant, the
> resulting instability becomes so tenacious that it **necessitates employing
> particular numerical remedies**. Since the inflammatory process is slow, we
> can use **small values for `R_s`** and hence the numerical code does not fail
> even in the absence of stabilizing terms."

Klempt 2024 runs `r = 100`, a *large* coefficient on the same kind of term, in
a regime this Remark calls advection-dominant and says needs remedies. Our
reproduction uses explicit upwind finite differences, a crude stabilisation,
against their Galerkin FEM.

**This is the strongest support yet for the scheme hypothesis** in
`JAXFEM/klempt2024_resolution.py`, and it is the group's own statement rather
than our inference.

## 4. Boundedness by penalty, not by clipping

Eqs. 18–19 enforce `φ ∈ [0,1]` with a penalty inside the variational form,
`P(φ) = ⟨φ − 1⟩² + ⟨−φ⟩²`, rather than by clipping after the update.
`klempt2024_quantitative.py` clips, which is listed there as an unseparated
difference from the paper; this is what the lineage does instead.

## 5. Smaller confirmations

- **Quasi-static nutrient, third time**: "The time scale of the diffusion
  process is substantially smaller than that of the inflammation process ...
  This is why the time dependent term is eliminated from the nutrient transport
  equation." (Klempt 2024 and Soleimani 2019 say the same.)
- `D = φ D_min + (1−φ) D_max` (Eq. 10) — diffusivity varies with the phase
  field, which we do not do.
- Double well `f(φ) = 16Mφ²(1−φ)²` (Eq. 12), against `2M□²(1−□)²` in Soleimani
  2023 — the coefficient is a convention, worth not copying blindly.
- The mechanics is an **anisotropic HGO** free energy with collagen fibres
  (Eq. 5); ours is isotropic neo-Hookean. Not a defect, but it is the shape the
  group reaches for when the tissue has fibres.
