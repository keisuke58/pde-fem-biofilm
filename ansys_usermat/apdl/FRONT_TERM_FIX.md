# Front-growth term in the partner element: fix and acceptance checks (2 Oct)

Approved on 2 Oct: correct the reading of the orientation weights, and use
`|∇φ|` as in Klempt et al. 2024 Eq. 34 with an upwind discretisation.
The partner's source stays on IKMHIWI03 (never committed); this sheet says
what to change and how to judge it.

Background:
- `KLEMPT2024_REPRODUCTION.md` §8 (the slip, the |∇²φ| instability, the
  cell Péclet number ≈ 31);
- chapter 5, framework property 4.

## Changes (in a copy first, e.g. `F:\biofilm_upf_diag`)

1. **Reading slip, USolBeg.** `ORI_WEIGHT12` and `ORI_WEIGHT22` are both read
   into `sGdp_OriWeight11`. Read each into its own variable (`ORI_WEIGHT22` →
   the species-2 weight), so that species 1 keeps its `ORI_WEIGHT11`.
2. **`|∇φ|` instead of `|∇²φ|`, USSFin.** In the front term
   `−Ori·Growth`, replace the factor `|Lap(bio_i)|` by `|∇bio_i|`. The routine
   already computes the gradient for the orientation.
3. **Upwind the front term.** The term is an advection of φ with velocity
   `v = r·c/(K+c)·∇c/|∇c|`, and the paper's form is `−v·∇φ`. Evaluate `∇φ`
   from the **upwind side** of v: of the NEM neighbours, use those lying
   against the direction of v (`(x_j − x_i)·v < 0`) when forming the gradient
   used in this term. Keep the central gradient for everything else.
   - If the NEM operator does not allow a one-sided gradient, fall back to
     clipping φ to [0, 1] after the update, as `klempt2024_quantitative.py`
     does. That is bounded but not mass-conserving, so report it as such.
4. Leave the β-term, the local source `K_LOCAL·locbio` and the penalty
   unchanged.

## Acceptance checks, Klempt 2024 Fig. 7 setup (`ds_fig7h.dat`, seed 92 93 100 101)

| check | pass |
|---|---|
| bounds | 0 ≤ φ ≤ 1 everywhere, every substep (allow 1e−6) |
| spreading | the seed spreads towards NUTRIENT1 (y = −1). Paper Fig. 7 (high nutrient): mean φ ≈ 0.55 at T* = 0.1, ≈ 1.0 at T* = 0.2. Report the curve; matching the paper is a result, not a pass condition |
| direction | growth on the side facing the nutrient, none on the far side |
| time step | same result at dt and dt/2 (difference in mean φ < 1 %) |
| mesh | if time allows, one refinement: the front position should move by less than one element |
| unchanged parts | checks 1–3 of the 5 Oct slides still pass (gradient-free region and whole-model α ratio do not involve the front term) |

Also record the global mass `Σ φ·V` against time. Upwinding conserves it up
to the source terms; clipping does not.

## After it passes

- Re-run modes 7 and 8 on the Fig. 7 setup. With a spreading front, the
  composition should vary in space: by local amount in mode 7, by arrival time
  in mode 8. This is the result chapter 5 still lacks (§5.7, limitation "no
  spatial composition yet").
- Tell me the numbers; chapter 5, property 4 and the limitation will be
  rewritten from them.
