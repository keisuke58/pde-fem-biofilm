# Dysbiotic-static (DS) composition bug — fix status (2026-06-26)

## Bug (confirmed, pervasive)
_multiscale_2d_results/ref_0d_dysbiotic_static.json had phi_final =
[So 0.944, 0.011, 0.011, 0.011, 0.011] -- a commensal-like So-dominant
composition (copy error) for a DYSBIOTIC condition. Raw CLSM dysbiotic-static is
V.dispar-dominant (So 0.03-0.15) at ALL days; the real 15-D TMCMC DS is Fn/An-
dominant. The So-dominant value propagated to gen_tooth, the Abaqus DS run
(E_voigt=9.6e-4, alpha_max=1.9), MAP_SIGMA[DS], MAP_PHI[DS], and samples_0d.

## Corrected sigma (real Abaqus, raw CLSM DS composition)
  buggy So-dominant : 13.92 kPa  (~ headline 13.63)
  corrected D6      : 1.72 kPa
  corrected D10     : 1.55 kPa   <- used (robust ~1.6 across days)
  corrected D21     : 1.56 kPa
=> correct DS sigma ~ 1.6 kPa; the headline 13.63 was ~8x too high.

## Fixed
 - ref_0d_dysbiotic_static.json: phi_final -> raw CLSM D10 (Vd-dom), di_0d
   recomputed (1-H/ln5): 0.845 -> 0.257.
 - posterior_klempt_stress_ci.py (tooth): MAP_PHI[DS], MAP_SIGMA[DS] 13.63->1.55.

## STILL TO FIX (incomplete)
 - samples_0d.json (dysbiotic_static): still So-dominant -> the DS posterior
   distribution is still ~13 kPa (now inconsistent with MAP 1.55). Needs the 0D
   sample set regenerated with the correct composition (or replaced by the
   experimental composition + replicate-based CI).
 - MAP_SIGMA_IMPLANT[DS] = 6.10e-3: still from the So-dominant comp; needs the
   implant DS Abaqus re-run.

## Impact
 - Headline sigma_CH/sigma_DH = 6.44x: UNAFFECTED (CH/DH, DS-independent).
 - Any DS-based number (sigma_DS, CS/DS ratio, 4-condition figure) was wrong
   (DS ~8x too high) and is now partially corrected (tooth MAP only).
