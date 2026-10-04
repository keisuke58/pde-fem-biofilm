# Towards a two-way coupling (written 2 Oct 2026)

**State at the end of the thesis:** the coupling is one-way, macro → micro.
- The element's field φ (Klempt et al. 2024) drives growth, α, and the
  stress.
- The point model (Klempt et al. 2026) runs at every Gauss point and gives
  the species composition (two schemes, `prop(28) = 7` and `8`; see
  `ansys_usermat/apdl/SPLIT_COUPLING.md`).
- Nothing flows back from the point model to the field or the mechanics.

This is consistent with Klempt 2024, where no material parameter depends on
the species. Chapter 5 states it as a limitation.

## Four steps back, ordered by how well the literature supports them

| step | what flows back | literature basis | prerequisite |
|---|---|---|---|
| 1 | **field nutrient → point model**: replace the constant c* by the local c | Klempt 2026: c* is "a given variable … possibly time-dependent" | implementation only: a fourth input to the bridge (`g`, `theta`, `dt`, **`c`**) |
| 2 | **composition → spreading**: front factor `R_s = Σ χᵢ R_s,i`, consumption `g = Σ χᵢ gᵢ` | Klempt, Soleimani, Junker 2023 (PAMM): both "depend on the type of microorganism" | a working front term in the element (`FRONT_TERM_FIX.md`: reading slip, \|∇φ\|, upwinding, and the authors' treatment of ∇c/\|∇c\| where ∇c → 0) |
| 3 | **composition → stiffness**: `E = Σ χᵢ Eᵢ` | none for species-specific moduli; Böl et al. 2013 explain why no reliable values exist | measured moduli per species |
| 4 | **stress → growth** | Klempt 2024 Eq. 30 keeps the mechanical term `φμ(I:C_e − 3)` that Eq. 34 drops; Soleimani et al. 2020 (stress-induced anisotropic growth) | the largest change: the field equation itself, and a stiffness-consistent time scale |

**Suggested order: 1 → 2 → 3 → 4.**
- Step 1 needs no new physics.
- Step 2 is the first real two-way coupling, and the PAMM paper supports it.
  Together they are a natural first goal for the continuation at Keio.
- Steps 3 and 4 need data or a new formulation.

## Step 1 tried in Python (4 Oct 2026)

`ansys_usermat/composition_local_nutrient.py` puts the local nutrient into the
point model on the partner's 8^3 model (steady nutrient with first-order
consumption in the seed, c held at y = -1 mm; Thiele number Lambda scanned
because the partner's d and g are example inputs). Case 3 hardly changes
(share 0.63-0.66). In case 6 the takeover is slower where c is low: for
Lambda = 4 the share in the seed ranges from 0.03 near the nutrient face to
0.44 in the interior, converged in grid (16^3, 32^3) and coupling step (0.025)
(`assets/fig_composition_local_nutrient.png`). So step 1
alone gives a composition that varies in space in ANSYS, without the front
term and without waiting for the authors.

**Bridge side done (4 Oct 2026), tested off-machine** (`tests/test_local_nutrient.py`):
- `material_server.py`: an optional `c_rel` in an ecology request scales c* for
  that call (c* = c*_0 c_rel; absent or 1: bit-identical to before; negative:
  refused);
- `biofilm_py_eval.c`: `biofilm_ecology_eval_c(g, theta, dt, n_sub, c_rel, g_new,
  phi_int)`; a negative `c_rel` sends nothing;
- `usermat_py_hook.f`: `biofilm_ecology_hook_c(g, theta, dt, n_sub, c_rel, g_new,
  phi_int, ok)`.

Left for IKMHIWI03: read the element's nutrient at the Gauss point in the
material call (the partner's `Nut1`; check that it is available there),
normalise it by the held value, call `biofilm_ecology_hook_c` instead of
`biofilm_ecology_hook` in the composition fragment, relink, and re-run case 6
with a consumption strong enough for c to drop across the seed.

## Steps 2-4 prototyped in Python (4 Oct 2026)

- **Step 3, composition -> stiffness** (`ansys_usermat/two_way_step3.py`): with the
  seed composition of step 1 (case 6, chi_1 = 0.12-0.44) and E_2/E_1 up to 5 at
  the same mean modulus, the seed's average stresses change by at most 0.2 %.
  The seed is about 1000 times stiffer than the void, so the soft surroundings
  set its stress. The step matters only where the biofilm itself carries the
  load (a contiguous biofilm on a substrate), not for a seed in a void.
- **Step 4, stress -> growth**: in this model |p|/mu in the seed is about
  4e-5 Pa / 3.36 Pa = 1e-5, so a stress-dependent growth law changes nothing
  with the Klempt 2024 parameters, while Eq. 30's mechanical term as written
  (divided by eta_phi = 1e-10) would stop all growth (KLEMPT2024_REPRODUCTION.md
  sec. 12). Step 4 needs a different growth law (e.g. Soleimani et al. 2020) and
  realistic stress levels.
- **Step 2, composition -> spreading** (`ansys_usermat/two_way_step2.py`, on the
  reproduced Klempt 2024 test case 4.1 field): front growth rate
  r = sum chi_i r_i with r_{1,2} = (1 +- d) r. d = 1/3 follows from eta_1 = 1,
  eta_2 = 2 of Klempt 2026 cases 3 and 6 if the front rate scales like 1/eta_i
  (assumption): case 3 0.704 -> 0.735, case 6 0.704 -> 0.648 (dt 0.025).
  Sensitivity d = 0.5 (r_1 = 1.5 r, r_2 = 0.5 r): at T* = 1
  the mean phi goes from 0.704 to 0.750 in case 3 (the faster species is the
  majority) and to 0.618 in case 6 (the slower species takes over); phi changes
  by up to 0.41 / 0.55 locally, the composition hardly at all. Converged in the
  coupling step (0.1 / 0.05 / 0.025: case 6 0.599 / 0.611 / 0.618). The feedback is
  real and large; in ANSYS it needs the element's front term.
  Grid (dt 0.05, d = 0.5): n = 21 gives +6 % / −13 %, n = 41 (0.5 um) +5 % / −8 %,
  with the one-way mean phi itself moving 0.704 -> 0.674. Sign and order hold;
  the size depends on the grid of the reproduced field.
  It reaches the stress: with d = 1/3 the mean von Mises stress in the biofilm
  changes by +4 % (case 3) and −12 % (case 6), about as much as the biofilm volume
  (CalculiX, `ansys_usermat/calculix/README.md`).

## Checks to carry over

Each step must keep the checks that pass today:
- the gradient-free exact solution (sinh/cosh);
- reduction to the point model when the amounts agree
  (`tests/test_composition_reduces_to_point_model.py`);
- bit-exact replay of every point-model call.

Step 2 adds one check of its own: with all `R_s,i` and `gᵢ` equal, it must
reduce exactly to the one-way scheme.

## Literature for the two-way steps (searched 2026-10-04)

Aim: replace assumptions of this work by published precedents. Feng et al.
2021 read in full (PDF in `references/`, 4 Oct); the others from abstracts.

- **Step 2, composition into the spreading: Feng, Neuweiler, Nogueira,
  Nackenhorst (2021)**, Bull. Math. Biol., doi:10.1007/s11538-021-00888-2
  (open access, PMC7990864; Hannover, Nackenhorst group).
  - Two-species oral biofilm (S. gordonii, Veillonella), continuum, FEM.
  - The biomass spreads with a potential flow driven by the **sum of the
    species' growth**: ∇²Φ = g₁/ρ + g₂/ρ in the biofilm, u = ∇Φ. Same idea as
    step 2 (composition sets how fast the biofilm spreads), so the mechanism
    would be a published one rather than mine.
  - Species rates: μ₁ = 3e-5 1/s (S. gordonii, measured, Rath et al. 2017,
    FEMS Microbiol. Ecol.), μ₂ = 8e-5 1/s (Veillonella, estimated). Ratio 2.7,
    i.e. d = (2.7 − 1)/(2.7 + 1) ≈ 0.45, between d = 1/3 (from η) and the
    sensitivity case d = 0.5. Run (dt 0.025): case 3 0.704 -> 0.746, case 6
    0.704 -> 0.626. With d in [1/3, 0.5] case 6 gives mean φ
    0.648 to 0.618 at T* = 1.
  - No mechanics, no species-specific stiffness ("fluid-structure interaction
    is not considered").
  - Checked in the PDF (Eq. 1, 3, 4, 14, 15): g_i = ϑ_i ρ μ_i × Monod factors,
    and with Σϑ_i = 1 the divergence of the growth velocity is
    ∇·u = Σ_i ϑ_i μ_i (Monod)_i, the **share-weighted rate** of step 2 exactly.
    Both species move with the same velocity (their assumption, after Alpkvist
    and Klapper 2007). The difference to step 2: there the weighted rate drives
    a potential flow of the whole biofilm, here it scales the front term of
    Klempt 2024.
  - μ₂ is marked "Estimated" in their Table 1, only μ₁ is measured; so the
    ratio 2.7 is a literature value for one species and an estimate for the
    other. Usable as a second, independent choice of d, not as a measurement.
  - CC BY 4.0, bundled in `references/` with a text extraction.
- **Rath, Feng, Neuweiler, Stumpp, Nackenhorst, Stiesch (2017)**, FEMS
  Microbiol. Ecol.: measured S. gordonii biofilm growth, the source of μ₁.
- **Soleimani et al. 2023** (already cited): summing S₁ + S₂ gives the
  share-weighted rate; same rate for both species (R_s = 500).
- **Klempt, Soleimani, Junker (2025), arXiv:2509.01274** = the point model
  (Klempt et al. 2026). Material point model, no mechanics, no spatial
  extension in the outlook; "a higher viscosity leads to a slower reaction"
  (Sec. 3.1) supports the direction of the 1/η_i assumption, not its form.
- **Step 3, species stiffness:** no moduli for oral species found. Species
  differences in biofilm rheology exist for other species (nonlinear rheology
  of single-species biofilms, PMC7156450; Peterson et al., viscoelasticity
  review). Stays a sensitivity parameter.
- **Step 4, stress into growth:** growth inhibited by compressive stress,
  coarse-grained to a continuum law without free parameters (arXiv:2603.28630,
  not yet read); Soleimani 2020 (already cited).

### Measured species parameters (searched 2026-10-04)

- **Growth rates, planktonic:** P. gingivalis doubling time about 3 h to 9 h
  depending on the medium (PMC9788703, PMC4083621, Frontiers fcimb 2023.1193198);
  oral streptococci and Actinomyces on glucose-supplemented saliva 1.6 h to 4 h
  (PMC240009-series, Springer bf00393856); S. gordonii μ = 3e-5 1/s, i.e. about
  6.4 h (Rath 2017, via Feng 2021). A. naeslundii does not grow on saliva alone,
  S. oralis only weakly and not reproducibly as a biofilm (Palmer group,
  saliva-grown communities). So the ratio between a fast early coloniser and
  P. gingivalis is roughly 2 to 5 (d ≈ 0.33 to 0.67), and it depends on the
  medium at least as much as on the species.
- **Stiffness:** no Young's modulus per oral species found. Multi-species
  saliva microcosm biofilms (Pattem et al. 2018, Sci. Rep., AFM): 14 to 41 kPa
  nutrient-poor, 0.55 to 2.6 kPa with 5 % sucrose, i.e. a factor 10 to 70 from
  the nutrient condition alone. The condition changes the modulus more than any
  species ratio used in step 3 (up to 5).
- **Consequence:** literature values bound d to about 1/3 to 2/3, but do not
  fix it, because they are condition dependent and planktonic. A value for this
  system could come from the Heine time courses themselves (initial growth of
  each species per condition), i.e. from the data already used for TMCMC.
