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

## Result on IKMHIWI03, 2 Oct (test copy `F:\biofilm_upf_diag`)

**Implementation:**
- the reading slip is fixed;
- the front term is `−w·r·c/(K+c)·n_c·∇φ`;
- ∇φ comes from a weighted least-squares fit (w = 1/|d|²) over the NEM
  neighbours upstream of `a = w·r·c/(K+c)·n_c`, with the central gradient
  where fewer than 3 neighbours lie upstream.

Fig. 7 setup (`ds_fig7h`, scaled Table 2), T* = 0.3, 0 errors.

| check | result | |
|---|---|---|
| upper bound | max φ 0.999 | pass |
| lower bound | min φ −0.011 | fail (small) |
| spreading | mean φ 0.031 / 0.045 / 0.052 / 0.053 / 0.051 / 0.047 at T* = 0.05 … 0.30 (paper: 0.55 at 0.1, 1.0 at 0.2) | far slower than the paper |
| time step | dt 1e−3 vs 5e−4: mean φ differs 1.1–2.6 % | fail (1 %) |
| direction | slab means along y at T* = 0.3 (slab 1 = nutrient face): front on 0.039 0.077 0.066 0.056 0.047 0.040 0.032 0.019; front off 0.010 0.044 0.010 0 0 0 0 0 | fail: growth towards the nutrient **and** away from it |

**Cause (diagnosis on IKMHIWI03, agreed):**
- Beyond the seed the nutrient is almost uniform (c by slab 1.000 0.996
  0.993 0.993 …), so ∇c ≈ 0.
- Eq. 34 normalises the direction, `n_c = ∇c/|∇c|`, which stays a unit
  vector there while the speed `r·c/(K+c) ≈ 0.5 r` does not vanish, so
  round-off sets the direction.
- This is a property of Eq. 34's normalised `n_c` (0/0 where c is uniform),
  not of the upwinding.

**Options:**
- (a) switch the term off where |∇c| < 1 % of its maximum: a
  regularisation, with the threshold as a new assumption;
- (b) unnormalised ∇c: a different model, so not recommended;
- (c) stop and record it as a limitation.

## Front build for WP3, 7 Oct (`F:\biofilm_upf_front`)

Felix Klempt's answer of 5 Oct (n_c = grad c/(|grad c| + eps), eps about
1e-8 to 1e-12) settles the 0/0 above, so the fix is carried into a build for
the Keio WP3 runs (`RUN_WP3_IKMHIWI03.md`, block B):

- base: the sources of the native build (`F:\biofilm_upf_native`: the 5 Oct
  wired sources with the point model in Fortran, no material server needed);
- USolBeg: the reading slip fixed as on 2 Oct (`ORI_WEIGHT12` and
  `ORI_WEIGHT22` into their own variables). In the old builds the last
  assignment set the species-1 weight to `ORI_WEIGHT22` (0 in `ds_fig7h`),
  which is why the front term did not act there;
- USSFin: the front term of species 1 is the 2 Oct upwind version, with
  n_c = grad c/(|grad c| + 1e-8) in place of the cut-off at |grad c| > 1e-14,
  so the term goes to 0 where c is uniform. Eq. 34 as printed otherwise (no
  w-blend, no factor in phi);
- one output line per substep, `FRONTPHI time mean min max` of phi, for the
  bounds and mass checks; the 2 Oct per-point and slab prints are not kept.

Built 7 Oct with `link_v222.ps1` (0 errors). Acceptance runs on `ds_fig7h`
(8^3, T* = 1): dt 0.005, dt 0.0025, and dt 0.005 without the front term
(`MAX_GROWTH11 = 0`) for the direction check; queued behind the week chain
(`_chain_front_accept.log`). Results follow here.

**Recommendation (cloud session): (c) for the thesis.** No result of chapter 5
depends on the front term. (a) would add an assumption of my own. How Klempt
et al. handle `∇c/|∇c|` where `∇c → 0` is a question for Felix: Table 1 of
the 2024 paper uses the normalised form in the weak form. Decide (a) or its
equivalent at Keio from his answer.
