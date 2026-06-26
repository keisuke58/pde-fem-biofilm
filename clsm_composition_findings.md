# Composition provenance & raw-CLSM recomputation of sigma_CH/sigma_DH (2026-06-26)

The downstream stress pipeline uses MAP_PHI compositions that are NOT the 15-D
TMCMC posterior (multimodal, sparse data: 25 obs / 15 params; the real posterior
has commensal So<=0.36). Composition is directly MEASURED by CLSM, so it should
be anchored to the measurement, not inferred from the under-identified TMCMC.

## Canonical source data
data_5species/experiment_data/expected_species_volumes.csv  (all 4 conditions,
species_fraction_median, time-resolved D1..D21). NOTE: fig3_species_distribution_
summary.csv in the same dir disagrees on V.dispar (Vd~=0 vs ~0.3) -- a data-
cleaning inconsistency to resolve. expected_species_volumes matches the dysbiotic
MAP_PHI (Vd-dominant) and is the better source.

## Composition is strongly time-varying -> the ratio is timepoint-dependent
Validated framework sigma_cond = E_voigt(phi) * A * k_eff(phi)^b (sigma ∝ E_voigt
from the 2x-E test; A,b fit to the 4 MAP, reproduces them <1% and 6.43x).
Raw CLSM composition (HOBIC) -> sigma_CH/sigma_DH per day:
  D1 5.1x | D3 6.2x | D6 3.0x | D10 2.2x | D15 1.9x | D21 2.7x
The headline 6.44x corresponds to ~D3 (early, So-peak commensal So~0.89). At
mature timepoints (D10-21) the commensal community shifts away from So-dominance
(So->0.36-0.46, An/Vd rise) and the ratio falls to ~2x.

## Conclusion / recommendation
- Composition must be anchored to expected_species_volumes.csv at a JUSTIFIED
  timepoint (the headline implicitly uses ~D3 / early biofilm).
- "Commensal stress ~6x dysbiotic" is an EARLY-biofilm result; mature biofilm
  is ~2x. State the timepoint, or report the range across maturation.
- Reframe: composition = experimentally measured (CLSM); TMCMC calibrates the
  interaction parameters (A matrix), not the composition.

Tools: experiment_data/clsm_composition.py, experiment_data/recompute_with_clsm.py
(in Tmcmc202601/data_5species/).
