"""E_SPEC sensitivity of sigma_CH/sigma_DH.
Validated (2x-E UMAT test: sigma scales 2.00x): for growth-induced stress the
strain is set by constrained growth (E-independent), so sigma_cond ∝ E_voigt_cond
= sum_i phi_i*E_i. Hence
    sigma_CH/sigma_DH = (E_voigt_CH/E_voigt_DH) * G,   G = E-independent growth factor.
Calibrate G from the headline, then sweep the (literature-unjustified) E_SPEC."""
import numpy as np

phi_CH = np.array([0.942, 0.012, 0.012, 0.011, 0.011])   # So-dominant
phi_DH = np.array([0.097, 0.119, 0.474, 0.123, 0.093])   # Vd-dominant
SP = ["So", "An", "Vd", "Fn", "Pg"]

E_assumed = np.array([1e-3, 8e-4, 6e-4, 2e-4, 1e-5])      # current assumption

def evr(E):  # E_voigt ratio CH/DH
    return float(phi_CH @ E) / float(phi_DH @ E)

# Calibrate G from headline (sigma ratio 6.44 at E_assumed)
G = 6.44 / evr(E_assumed)
print("Calibrated growth factor G (E-independent) = %.3f" % G)
print("(headline sigma_CH/sigma_DH=6.44x = E_voigt-ratio %.3f x G %.3f)\n" % (evr(E_assumed), G))

scenarios = {
    "current assumed         ": E_assumed,
    "uniform E (no contrast) ": np.full(5, 6e-4),
    "mild contrast (2x range)": np.array([8e-4, 7e-4, 6e-4, 5e-4, 4e-4]),
    "strong contrast (1000x) ": np.array([1e-3, 5e-4, 2e-4, 5e-5, 1e-6]),
    "reversed (Pg stiff)     ": np.array([1e-5, 2e-4, 6e-4, 8e-4, 1e-3]),
    "So/Pg only 10x          ": np.array([1e-3, 1e-3, 1e-3, 1e-3, 1e-4]),
}
print("%-26s E_voigt_CH/DH   sigma_CH/sigma_DH" % "E_SPEC scenario")
for name, E in scenarios.items():
    r = evr(E)
    print("  %-26s %6.2f          %5.1fx" % (name, r, r * G))

print("\nNote: with the calibrated growth factor G=%.2f, the LOWER bound of the" % G)
print("ratio (uniform E, no species stiffness contrast) is %.1fx; the headline" % (evr(np.ones(5)) * G))
print("6.44x requires the assumed So-stiff / Pg-soft contrast (no literature/AFM).")
