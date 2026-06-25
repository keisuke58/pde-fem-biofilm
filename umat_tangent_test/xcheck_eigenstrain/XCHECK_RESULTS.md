# Growth cross-check: eigenstrain vs Klempt UMAT (2026-06-26)

Goal: elevate the eigenstrain growth path to an independent cross-check of the
Klempt UMAT (Fg=(1+α)I), after the rigor audit found the legacy `eps=α/3`
mis-scaling. Single C3D8, bottom face clamped, NLGEOM, E=1e-3 MPa.

## Results (growth-induced Mises, α=0.1, ν=0.49 unless noted)

| Implementation | growth | elastic law | Mises [MPa] | vs UMAT |
|---|---|---|---|---|
| `umat_klempt_voigt.f` (production) | Fg=(1+α)I | Klempt Hencky | 6.375e-5 | — |
| `umat_biofilm_visco.f` (η=0) | Fg=(1+α)I | Neo-Hooke C10/D1 | 6.382e-5 | **+0.1%** |
| built-in `*Expansion` eps=α | thermal | Neo-Hooke C10/D1 | 7.375e-5 | +16% |
| built-in `*Expansion` eps=α/3 | thermal | Neo-Hooke C10/D1 | 2.642e-5 | **−59%** |

## Findings

1. **Growth magnitude — `eps=α` is correct, `eps=α/3` is wrong.**
   Free-expansion diagnostic: `*Expansion` eps=0.1 gives nodal U=0.1 per
   direction ⇒ stretch exactly 1.1 = Fg(0.1). So the thermal eigenstrain
   reproduces Klempt's Fg kinematics *exactly* at eps=α; eps=α/3 would give
   stretch 1+α/3 and under-predicts stress 2.4× (−59% Mises). The legacy /3 in
   `compute_alpha_eigenstrain.py` is decisively confirmed wrong.

2. **Dual-UMAT agreement = the strong cross-check.**
   Two independent UMATs — `umat_klempt_voigt.f` (Klempt-Hencky elastic) and
   `umat_biofilm_visco.f` (Neo-Hooke elastic), both explicit Fg=(1+α)I — agree
   on the constrained growth-induced stress to **0.1%**. Different elastic
   potentials, different code, same growth kinematics → same stress. This
   independently validates the headline's growth-stress computation. (Two
   different W agreeing to 0.1% here because α=0.1 strains are still moderate.)

3. **Built-in `*Expansion` is NOT a clean equivalent for *constrained* stress.**
   Despite identical free stretch (1.1) and matched C10/D1, it sits ~16% high
   at ν=0.49 and ~10% at ν=0.30 (in-plane S11 differs ~60% throughout). The gap
   is partly volumetric locking (ν-dependent) and partly a genuine difference in
   how Abaqus couples thermal eigenstrain to hyperelasticity vs an explicit-Fg
   UMAT. → `*Expansion` carries a systematic offset and should not be presented
   as an equivalent reproduction of the UMAT.

## Recommendation

- **Cross-check vehicle = the dual-UMAT agreement (0.1%)**, not the eigenstrain.
  This already delivers the goal ("independent implementations agree ⇒ headline
  validated"), and more cleanly than the eigenstrain could.
- **Apply the `eps=α` fix** in the eigenstrain helper (the magnitude is correct;
  /3 was wrong) — done as an annotated `eps_growth_klempt`. If the eigenstrain
  path is kept for any figure, use eps=α and report its ~15% formulation offset
  vs the UMAT; otherwise drop it in favour of the UMAT.
- Headline σ_CH/σ_DH=6.44× is unaffected (UMAT-based) and now has an independent
  dual-UMAT validation of its growth-stress kernel.

Artifacts: `umat_tangent_test/xcheck_eigenstrain/` (inps + extract scripts;
re-run any with `abaqus job=<f> input=<f>.inp [user=...]`).
