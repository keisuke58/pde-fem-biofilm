# E_SPEC sensitivity of sigma_CH/sigma_DH (2026-06-26)

E_SPEC (per-species elastic moduli) is flagged in code as an assumed scaling with
no literature source (requires AFM validation). How much does the headline ratio
depend on it?

## Method (validated)
2x-E UMAT test: scaling all E_i by 2 scales sigma by exactly 2.00x (13.72->27.43
kPa). So for growth-induced stress the strain is set by constrained growth
(E-independent) and sigma_cond = E_voigt_cond * strain, with E_voigt=sum_i phi_i*E_i.
Hence  sigma_CH/sigma_DH = (E_voigt_CH/E_voigt_DH) * G,  G = E-independent growth
factor = 3.365 (calibrated from the headline).

## Result: the ratio IS sensitive to E_SPEC (unlike the depth model)
| E_SPEC scenario          | E_voigt CH/DH | sigma_CH/sigma_DH |
| current assumed (100x)   | 1.91          | 6.4x  (headline)  |
| uniform E (no contrast)  | 1.09          | 3.7x              |
| mild contrast (2x)       | 1.43          | 4.8x              |
| strong contrast (1000x)  | 3.69          | 12.4x             |
| reversed (Pg stiff)      | 0.08          | 0.3x              |

## Conclusion
- Robust floor ~3.4x from growth alone (G, E-independent, also depth-model-robust):
  the QUALITATIVE claim (commensal stress > dysbiotic by several-fold) is robust.
- ~half of the 6.44x (a factor ~1.9) comes from the assumed So-stiff/Pg-soft
  species stiffness contrast, which has NO literature/AFM basis. Uniform E -> 3.7x.
- Thesis should report the robust growth contribution (~3.4x) as the core finding
  and disclose the E_SPEC-dependent amplification as an assumption pending AFM data.
