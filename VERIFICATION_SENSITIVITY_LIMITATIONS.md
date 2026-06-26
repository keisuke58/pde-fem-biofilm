# Verification, Sensitivity, and Limitations (rigor audit 2026-06-26)

Pipeline:  CLSM species fractions -> TMCMC -> composition phi -> stiffness E(phi),
growth alpha -> FEM (Klempt UMAT) -> max von Mises sigma -> stress ratio
sigma_CH/sigma_DH.  This section reports what is verified, what the result is
sensitive to, and the disclosed limitations.  (Audit: JAXFEM/audit_all.py = ALL
CLEAR; eq/fig/reg/thesis/ci all pass.)

## 1. Verification (what is rigorous / verified correct)

V1. **Constitutive implementation (Klempt-faithful).** The growth kinematics
    F = Fe.Fg with Fg = (1+alpha)I (Klempt 2024, Sec. 2.1) are implemented
    exactly in the production UMAT (umat_klempt_voigt.f / umat_klempt2025.f,
    s_iso = 1+alpha). Compressible neo-Hookean (Hencky) elastic part.

V2. **Consistent tangent (DDSDDE).**
    - Production UMATs carry a closed-form analytic Jaumann tangent; verified by
      (i) the small-strain limit reducing to the exact isotropic tensor
      [lam+2mu, lam; mu], and (ii) a real single-element Abaqus job that
      converges quadratically (<=2 equilibrium iters, 0 cut-backs, auto-growing
      increment).
    - The viscoelastic extension UMAT (F=Fe.Fv.Fg) was given an exact algorithmic
      tangent by deformation-gradient perturbation (Sun et al. 2008); verified
      vs central-difference to 2.4-2.9e-8 (elastic, viscoelastic, Mooney-Rivlin)
      and converges in a real Abaqus increment loop.

V3. **Independent-implementation cross-check.** Two independent UMATs (Klempt-
    Hencky elastic and neo-Hookean elastic), both Fg=(1+alpha)I, agree on the
    constrained growth-induced stress to 0.1%.

V4. **Numerical convergence.** The 2D growth PDE is mesh-converged (domain-
    averaged fields < 0.03%); the through-thickness alpha profile is resolved
    (Nz 5->101 changes < 3%); the tooth-stress mesh has a convergence study.

V5. **Direct pipeline validation.** A surrogate-free pipeline (phi -> PDE ->
    alpha field -> tooth UMAT -> max Mises) reproduces the four calibration
    (MAP) stresses to +-0.0%.

V6. **Posterior credible interval.** The original CI surrogate
    (log sigma = c + a*log E_voigt + b*log k_eff, 4 points / 3 params) was
    statistically invalid (E_voigt and k_eff collinear, cond=2890; real leave-
    one-out CV: median 46%, max 7258%). It was replaced by the single-predictor,
    well-conditioned sigma = C*k_eff^b (cond=3.7; LOO median 2.3%, max 10.8%),
    with leave-one-out CV now computed and reported. The MAP ratio 6.44x is
    unchanged; the 90% CI tightened from [0.97, 16.02]x to [0.98, 10.29]x.

## 2. Sensitivity (what the headline result depends on)

S1. **Depth / nutrient model -> ROBUST.** The headline imposes a linear depth
    ramp alpha = alpha_max*(1-depth). Solving the growth reaction-diffusion PDE
    directly in 3D on the real 0.2 mm biofilm shell (established phi=1, nutrient
    from the oral surface) gives a different, outer-peaked profile. Sweeping the
    nutrient penetration depth across the literature range (O2 penetration
    50-200 um; de Beer 1994, Stewart 2003) and beyond (Lp/L = 0.15..10),
    sigma_CH/sigma_DH = 5.3-6.6x, bracketing the headline 6.44x. The RATIO is
    robust to the depth/nutrient model; absolute sigma is not.

S2. **Per-species stiffness E_SPEC -> SENSITIVE.** sigma scales linearly with
    E_voigt (verified: 2x-E -> 2.00x sigma), so sigma_CH/sigma_DH =
    (E_voigt_CH/E_voigt_DH) x G, with G = 3.37 (E-independent growth factor).
    Varying the (assumed) per-species moduli: uniform E -> 3.7x; current
    assumption -> 6.4x; strong contrast -> 12.4x. ~half of the 6.44x comes from
    the assumed So-stiff / Pg-soft contrast. The robust floor (growth alone) is
    ~3.4x.

S3. **Biofilm maturation (composition timepoint) -> STRONG.** Composition is
    strongly time-varying. Using the raw CLSM composition at each day:
    HOBIC sigma_CH/sigma_DH peaks early (~6x at D3, So-dominant commensal) and
    falls to ~2x at maturity (D10-21); Static sigma_CS/sigma_DS is opposite
    (~1x early -> ~5.7x at D10 -> ~4x). The single number 6.44x corresponds to
    early biofilm; the full sigma(t) trajectory is reported instead of one day.

S4. **Growth magnitude alpha -> magnitude data-anchored.** Measured biofilm
    thickness growth (CLSM z-profiles: ~1.5-3.5x) implies alpha ~ 0.5-2.5; the
    model alpha ~ 1-2 lies within this measured range.

## 3. Limitations (honest disclosure)

L1. **E_SPEC is an assumption, not a measurement.** The per-species elastic
    moduli have no direct literature source for all five species; only the
    order of magnitude (biofilm E ~ kPa) is literature-supported (AFM). The
    per-species CONTRAST, which materially affects the ratio (S2), requires
    species-specific AFM/rheology. Reported with the S2 sensitivity band.

L2. **Growth alpha cannot be precisely per-condition calibrated.** The CLSM
    thickness data are noisy / non-monotonic (e.g. CH: 43->64->28->53->42 um),
    incomplete (no Static), and inconsistent with the volume data. Only the
    magnitude is anchored (S4); precise per-condition alpha is not claimed to
    avoid over-fitting noise.

L3. **Composition is experimentally anchored, not TMCMC-inferred.** With at most
    25 observations for 15 parameters the inverse problem is under-identified;
    the posterior is strongly multimodal (>=4 modes confirmed; the paper
    discloses this). The 15-D posterior does not uniquely recover the measured
    composition, so the stress composition is taken from the CLSM measurement.
    TMCMC is therefore framed as calibrating the interaction parameters
    (A matrix), not the composition.

L4. **Corrected data bug (now fixed).** ref_0d_dysbiotic_static had a commensal-
    like So-dominant composition (copy error); the raw CLSM dysbiotic-static is
    V.dispar-dominant. This made sigma_DS ~8x too high (13.6 vs ~1.6 kPa). Fixed
    end-to-end (ref_0d, tooth/implant MAP, samples_0d). The headline ratio
    sigma_CH/sigma_DH (CH/DH) is independent of DS and unaffected.

L5. **Surrogate vs direct solve.** The CI uses the validated k_eff^b surrogate;
    a surrogate-free direct UMAT propagation pipeline exists (V5) and can replace
    it for the final tables if desired.

## Bottom line
The constitutive/FEM core is verified to continuum-mechanics (IKM) standard. The
qualitative claim -- commensal biofilm carries several-fold higher growth-induced
stress than dysbiotic -- is robust (>~3-6x) to the depth/nutrient model and to a
3D treatment. The precise number and the absolute stresses depend on assumptions
(E_SPEC contrast, biofilm maturation stage, growth magnitude) that are anchored
to data where possible and otherwise reported as sensitivity bands; the
remaining gaps require new measurements (species AFM, clean thickness) and are
disclosed rather than hidden.
