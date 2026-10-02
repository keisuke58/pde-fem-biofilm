"""klempt2024_cases.py -- the parameters and numerical examples of Klempt,
Soleimani, Wriggers & Junker, "A Hamilton principle-based model for
diffusion-driven biofilm growth", Biomech Model Mechanobiol 23:2091-2113
(2024), doi:10.1007/s10237-024-01883-x (PDF in references/,
references/Klempt2024_Hamilton_biofilm_growth_BMMB.pdf), as data -- the same idea as
klempt2026_cases.py / fritsch2025_cases.py -- plus how they map onto the
parameters of the partner's (Oliver's) ANSYS deck.

    from klempt2024_cases import TABLE2, CASES, PAPER_CURVES, OLIVER_DECK, deck_values
    deck_values("fig7_high")     # deck cube 2 mm -> {"MY_DIFF1": ..., "CONSUMPTION11": ..., ...}

Equations (Table 2 values are already divided by eta, so they are the
coefficients of these equations; time is the normalised T* in [0, 1]):
    Eq. 34  phi_dot - beta lap(phi) - k_a alpha_K
                    + |grad phi| r c/(K_M+c) n_gradphi . n_gradc = 0
    Eq. 35  c_dot - d lap(c) + g phi = 0      (quasi-static, d >> growth)
    Eq. 36  alpha_K_dot = k_a phi,  Fg = alpha_K I,  alpha_K(0) = 1
    Eq. 20  elastic energy weighted with phi^2; the void has no stiffness

CAUTION (klempt2024_quantitative.py, 2026-09-29/30): an independent
finite-difference solution of Eq. 34-36 exactly as printed, with these
values, does NOT reproduce the paper's Fig. 4 / Fig. 7 (variants the paper
does not use fit better). Treat the curves below as the paper's plotted
values, not as something the printed equations are known to give.
"""
from __future__ import annotations

# ---- Table 2 (post-division by eta; lengths in um, time in T*) ---------------
TABLE2 = {
    "d":     1e10,    # nutrient diffusivity          [um^2 / T*]   Eq. 35
    "beta":  2.0,     # phase-field regularisation     [um^2 / T*]   Eq. 34
    "k_a":   1e-3,    # growth factor                  [1 / T*]      Eq. 34, 36
    "K_M":   1.0,     # half-velocity (Monod) constant [-]           Eq. 34
    "g":     1e8,     # consumption                    [1 / T*]      Eq. 35 (Fig. 7 "low": 1e10)
    "r":     100.0,   # front growth rate              [um / T*]     Eq. 34
    "E":     10.0,    # Young's modulus                [Pa]          Eq. 20
    "nu":    0.49,    # Poisson's ratio (incompressible in the paper, H1P0)
    "mu":    3.3557,  # shear modulus E / (2 (1 + nu)) [Pa]
}

# ---- numerical examples (3D, 20 um cube, 1 um mesh, T* in [0, 1]) -------------
CASES = {
    "fig4_edge": {
        "figure": "Fig. 4 (test case 1)",
        "domain_um": 20.0, "mesh_um": 1.0,
        "seed": "sphere, radius 5 um, at the cube centre",
        "nutrient": "c = 1 supplied at one edge of the cube",
        "g": 1e8,
    },
    "fig7_high": {
        "figure": "Fig. 7, 'biofilm on nutrients', high consumption",
        "domain_um": 20.0, "mesh_um": 1.0,
        "seed": "disk, diameter 5 um, just above the bottom face",
        "nutrient": "c = 1 on the bottom face",
        "g": 1e8,
    },
    "fig7_low": {
        "figure": "Fig. 7, 'biofilm on nutrients', low consumption",
        "domain_um": 20.0, "mesh_um": 1.0,
        "seed": "disk, diameter 5 um, just above the bottom face",
        "nutrient": "c = 1 on the bottom face",
        "g": 1e10,
    },
}

# domain averages read off the paper's figures by eye (+-0.01..0.02);
# same numbers as klempt2024_quantitative.PAPER
PAPER_CURVES = {
    "fig4_edge": {
        "t":   [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00],
        "phi": [0.13, 0.21, 0.28, 0.325, 0.41, 0.47, 0.53, 0.58, 0.63, 0.67, 0.705, 0.74],
        "c":   [0.33, 0.24, 0.19, 0.165, 0.135, 0.12, 0.11, 0.105, 0.098, 0.093, 0.089, 0.085],
    },
    "fig7_high": {
        "t":   [0.01, 0.05, 0.10, 0.15, 0.20, 0.50, 1.00],
        "phi": [0.02, 0.12, 0.55, 0.90, 1.00, 1.00, 1.00],
        "c":   [1.00, 0.80, 0.55, 0.50, 0.49, 0.49, 0.49],
    },
    "fig7_low": {
        "t":   [0.01, 0.05, 0.10, 0.15, 0.20, 0.50, 1.00],
        "phi": [0.02, 0.07, 0.17, 0.22, 0.25, 0.29, 0.30],
        "c":   [0.80, 0.42, 0.15, 0.10, 0.09, 0.08, 0.075],
    },
}

# ---- the partner's deck (F:\biofilm_upf_wired, not in the repo) ---------------
# Read off Ussfin_P21-V21_Conection_Test.F on 2026-10-02:
#   bio1 update: dt * ( - Ori * Growth + K_LOCAL1 * locbio1 + MY_BETA1 * lap(bio1) + penalty )
#   Growth  = |lap(bio1)| * MAX_GROWTH11 * c / (HALF_VELO11 + c) (+ nutrient 2 term),
#   Ori = n_gradc . n_gradphi
#   !! DIFFERS FROM Eq. 34: the front term uses |lap(phi)| (Sdp_LapBio1, the
#      same Laplacian as the beta term), where Eq. 34 has |grad phi|; |grad phi|
#      (Sdp_NormBio1) is computed in the same routine but not used there.
#      So MAX_GROWTH11 has units L^2/T* there (r: L/T*); it is NOT r one-to-one,
#      and the r scaling in deck_values() only holds for a |grad phi| front term.
#   !! AND THE FRONT TERM IS SWITCHED OFF: USolBeg reads ORI_WEIGHT12 and
#      ORI_WEIGHT22 into sGdp_OriWeight11 (copy-paste slip), so the last one
#      read, ORI_WEIGHT22 = 0, sets species 1's weight to 0 and Ori = 0 always
#      (found 2026-10-02 by printing the factors in a test copy). In every
#      partner-element run so far the biofilm cannot spread.
#      Test copy F:\biofilm_upf_diag, Fig. 7 setup: slip fixed + |grad phi| ->
#      the seed spreads (phi >= 0.5 at 4 -> 94 traced points by T* = 0.18), but
#      phi exceeds 1 (up to ~3; the penalty does not bound it) and the run
#      breaks down near T* = 0.22. Slip fixed + |lap phi| -> unstable at
#      dt = 0.005 (front rate ~560 / T*). Neither is used for thesis results.
#   nutrient: quasi-static  MY_DIFF1 * lap(c) = CONSUMPTION11 * bio1 (+ CONSUMPTION12 * bio2),
#             c held at MY_NUTSTART1 on the Gauss points of component NUTRIENT1
#   index order of MAX_GROWTH / HALF_VELO / CONSUMPTION: %nutrient%%species%
OLIVER_DECK = {
    "map": {                        # deck parameter -> Klempt 2024 symbol
        "MY_DIFF1": "d", "CONSUMPTION11": "g", "MY_BETA1": "beta",
        "K_LOCAL1": "k_a", "MAX_GROWTH11": "r (but multiplies |lap phi|, not |grad phi|)",
        "HALF_VELO11": "K_M",
        "YOUNG_BIO": "E (deck units MPa: 1e-5)", "POISSON_BIO": "nu",
        "YOUNG_VOID": "negative = Klempt phi^2 stiffness, value = -floor (e.g. -1e-3)",
    },
    "example_input": {              # what the example deck carries (NOT the paper)
        "MY_DIFF1": 1.0, "CONSUMPTION11": 0.5, "MY_BETA1": 1e-4, "K_LOCAL1": 0.1,
        "MAX_GROWTH11": 100.0, "HALF_VELO11": 0.1, "MY_NUTSTART1": 1.0,
        "MY_BIOSTART1": 1.0, "MY_BIOSTART2": 1.0, "YOUNG_BIO": 1000.0, "POISSON_BIO": 0.3,
    },
    "geometry": {
        "units": "mm (/units,MPA)", "cube_mm": 2.0, "range_mm": (-1.0, 1.0),
        "elements_per_side": 8, "element_mm": 0.25,
        "NUTRIENT1": "elements 1-64 = the layer on the y = -1 face -> c held at MY_NUTSTART1 there "
                     "(like Fig. 7's nutrient face)",
        "NUTRIENT2": "the layer on the z = -1 face (nutrient 2, start value 0 -> unused)",
        "BIOFILM1": "32 elements, a 4x4x2-ish block around the centre (incl. 220)",
        "BIOFILM2": "8 elements: 220 221 228 229 284 285 292 293",
        "fig7_seed": "elements 92 93 100 101: 2x2 at the centre of the 2nd layer from the "
                     "nutrient face, 0.5 mm wide = 1/4 of the cube, like the paper's 5 um disk "
                     "in 20 um (make_wired_deck.py --cmblock BIOFILM1=92,93,100,101)",
    },
}


def deck_values(case: str, L_deck: float = 2.0) -> dict:
    """Klempt 2024 coefficients for `case`, rescaled to a deck cube of side
    L_deck (deck length unit) with T* kept as the time unit, so that every
    dimensionless group of the paper is preserved:
        beta T*/L^2, r T*/L, g L^2/d, k_a T*, K_M
    d itself only has to stay >> growth (quasi-static nutrient), so the deck's
    MY_DIFF1 = 1 is kept and g is chosen to give the paper's g L^2/d."""
    c = CASES[case]
    L = c["domain_um"]
    s = L_deck / L                       # deck length per um
    d_deck = 1.0
    return {
        "MY_BETA1": TABLE2["beta"] * s ** 2,
        "MAX_GROWTH11": TABLE2["r"] * s,
        "HALF_VELO11": TABLE2["K_M"],
        "K_LOCAL1": TABLE2["k_a"],
        "MY_DIFF1": d_deck,
        "CONSUMPTION11": c["g"] * L ** 2 / TABLE2["d"] * d_deck / L_deck ** 2,
        "YOUNG_BIO": TABLE2["E"] * 1e-6,  # Pa -> MPa
        "POISSON_BIO": TABLE2["nu"],
        "_dimensionless": {
            "beta_T/L2": TABLE2["beta"] / L ** 2,
            "r_T/L": TABLE2["r"] / L,
            "gL2/d": c["g"] * L ** 2 / TABLE2["d"],
        },
    }


if __name__ == "__main__":
    import json
    for k in CASES:
        print(k, json.dumps(deck_values(k), indent=1))
