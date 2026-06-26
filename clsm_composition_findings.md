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

## CORRECTION (provenance, verified in TMCMC code)
The TMCMC estimation code (Tmcmc202601/data_5species/main/estimate_reduced_
nishioka.py L1688) loads **fig3_species_distribution_summary.csv** (or
species_distribution_data.csv) -- NOT expected_species_volumes.csv (no code in
Tmcmc202601 reads expected_species_volumes). So:
 - TMCMC INPUT = fig3_summary  (dysbiotic HOBIC: V.dispar ~= 0, An/Fn-dominant)
 - expected_species_volumes.csv (V.dispar-dominant) is an UNUSED auxiliary file.
 - The downstream stress MAP_PHI (V.dispar-dominant, Vd=0.47) matches the UNUSED
   expected_species_volumes, NOT the TMCMC input (fig3, Vd~=0). The two Heine
   CSVs disagree on V.dispar.
=> The stress composition is consistent with neither (a) the 15-D TMCMC posterior
   nor (b) the data the TMCMC was fit to (fig3). Which Heine CSV is authoritative
   for composition must be resolved by the author; it changes E_voigt(phi) and
   hence the ratio. The recompute_with_clsm.py timepoint result (6.2x@D3 -> ~2x
   mature) used expected_species_volumes (the MAP_PHI-consistent convention); a
   fig3-based recompute would differ for dysbiotic (no V.dispar).
