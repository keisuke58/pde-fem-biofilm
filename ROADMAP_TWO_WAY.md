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

## Checks to carry over

Each step must keep the checks that pass today:
- the gradient-free exact solution (sinh/cosh);
- reduction to the point model when the amounts agree
  (`tests/test_composition_reduces_to_point_model.py`);
- bit-exact replay of every point-model call.

Step 2 adds one check of its own: with all `R_s,i` and `gᵢ` equal, it must
reduce exactly to the one-way scheme.
