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
(share 0.596-0.608). In case 6 the takeover is slower where c is low: for
Lambda = 4 the share in the seed ranges from 0.12 near the nutrient face to
0.44 in the interior (`assets/fig_composition_local_nutrient.png`). So step 1
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
- **Step 2, composition -> spreading** (`ansys_usermat/two_way_step2.py`): result
  to follow.

## Checks to carry over

Each step must keep the checks that pass today:
- the gradient-free exact solution (sinh/cosh);
- reduction to the point model when the amounts agree
  (`tests/test_composition_reduces_to_point_model.py`);
- bit-exact replay of every point-model call.

Step 2 adds one check of its own: with all `R_s,i` and `gᵢ` equal, it must
reduce exactly to the one-way scheme.
