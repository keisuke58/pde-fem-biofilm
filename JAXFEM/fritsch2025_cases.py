"""fritsch2025_cases.py -- the case studies of Fritsch, Geisler, Grashorn, Klempt,
Soleimani, Broggi, Junker, Beer, "Bayesian Updating of constitutive parameters
under hybrid uncertainties with a novel surrogate model applied to biofilms"
(arXiv:2512.15145v1; bib key Fritsch2025BayesianMicrofilms), transcribed as data
in the same layout as klempt2026_cases.py, so the repo's 0D Hamilton integrator
(klempt2026_reproduction.build / simulate) can run them directly.

Read from the arXiv HTML version (arxiv.org/html/2512.15145v1) on 2026-10-01;
the PDF itself was not opened. Same multi-species Hamilton model as Klempt et
al. 2026 (A = nutrient/interaction matrix, b = antibiotic susceptibility).

  Each entry carries its own "dt"; klempt2026_reproduction.build reads its
  module-level DT, so set klempt2026_reproduction.DT = case["dt"] before
  simulate() for the 1e-5 submodels.

  Case I  (two species, Table 1-2): synthetic data generated from the mean
          theta* = [a11, a12, a22, b1, b2] = [1, 0.1, 1, 1, 2]; 20 data points
          sampled on t in [0.05, 1]; aleatory CoV 0.5 % and 2.0 %.
  Case II (four species, Table 3): hierarchical updating in three submodels
          M1 (species 1-2), M2 (species 3-4), M3 (cross terms a13, a14, a23, a24),
          plus a validation run. theta* for Case II is given ONLY as the yellow
          "True Mean" bars of Fig. 16 (p. 41), read off by eye from the PDF on
          2026-10-01; every bar sits on a round value, listed in A4_TRUE / B4_TRUE.
          Posterior means from the same figure, for reference: a11 0.72, a12 2.09,
          a13 1.98, a14 1.02, a22 0.92, a23 2.01, a24 0.98, a33 1.50, a34 1.01,
          a44 2.00, b1 0.10, b2 0.20, b3 0.30, b4 0.57 (b4 poorly identified).

Note: the paper's Eq. 19 AND Eq. 20 both end in "+gamma"; the psi-equation
(Eq. 20) carrying gamma is the same misprint klempt2026_cases.py documents for
Klempt 2026 Eq. 17 (the constraint does not depend on psi).
Table 3 gives the validation antibiotics as 50*1[t > 500], but Fig. 18 switches
at normalised t = 0.5 of a 1500-step run (= step 750) and says t < 0.5 equals
Fig. 15 (the 750-step M3 run) -- the figure and the table disagree. Settled by
running (2026-10-01): switch at step 750 reproduces Fig. 18 (max |diff| 0.011),
step 500 does not (0.31), so Table 3's "500" is the misprint.

Checked 2026-10-01 with klempt2026_reproduction.simulate (DT set per case)
against the paper's figures read by eye (phibar = phi*psi per species):
Case I vs Fig. 3/5a max |diff| 0.013; M1 vs Fig. 11 0.003; M2 vs Fig. 13 0.007;
M3 vs Fig. 15 0.002; validation (step 750) vs Fig. 18 0.011. This also confirms
the Fig. 16 reading of theta*.

Unlike Klempt 2026, this paper states psi0 = 0.999 and Kp = 1e-4 explicitly.
No species-specific mechanical (stiffness / growth-to-deformation) values are
given anywhere in the paper.
"""
from __future__ import annotations

PSI0 = 0.999        # Table 1, stated in the paper
KP = 1.0e-4         # Table 1, penalty term


# Case II true parameters, Fig. 16 (read off the bar chart; see docstring)
A4_TRUE = [[0.8, 2.0, 2.0, 1.0],
           [2.0, 1.0, 2.0, 1.0],
           [2.0, 2.0, 1.5, 1.0],
           [1.0, 1.0, 1.0, 2.0]]
B4_TRUE = [0.1, 0.2, 0.3, 0.4]


def _alpha_validation(step):
    return 50.0 if step > 500 else 0.0                     # Table 3 as printed: alpha* = 50 * 1[t > 500]


def _alpha_validation_fig18(step):
    return 50.0 if step > 750 else 0.0                     # Fig. 18: switch at normalised t = 0.5 of 1500 steps


CASES = {
    # ---- Case I: two species (Table 1, Table 2) --------------------------------
    "fritsch_2sp": {
        "n": 2, "A": [[1.0, 0.1], [0.1, 1.0]], "b": [1.0, 2.0],
        "eta": [1.0, 2.0], "phi_init": [0.25, 0.30],
        "c_star": 100.0, "alpha_star": 10.0, "dt": 1.0e-4, "steps": 1000,
        "prior": {"a11": (0, 3), "a12": (0, 0.5), "a22": (0, 3), "b1": (0, 3), "b2": (0, 3)},
        "data": {"n_points": 20, "t_range": (0.05, 1.0), "cov": [0.005, 0.02]},
        "source": "Table 1-2, theta* = [1, 0.1, 1, 1, 2]",
    },
    # ---- Case II: four species, submodels (Table 3); theta* from Fig. 16 --------
    "fritsch_4sp_M1": {
        "n": 2, "species": [1, 2], "A": [[0.8, 2.0], [2.0, 1.0]], "b": [0.1, 0.2],
        "eta": [1.0, 1.0], "phi_init": [0.2, 0.2],
        "c_star": 100.0, "alpha_star": 100.0, "dt": 1.0e-5, "steps": 2500,
        "source": "Table 3, M1 (a11, a12, a22, b1, b2)",
    },
    "fritsch_4sp_M2": {
        "n": 2, "species": [3, 4], "A": [[1.5, 1.0], [1.0, 2.0]], "b": [0.3, 0.4],
        "eta": [1.0, 1.0], "phi_init": [0.2, 0.2],
        "c_star": 100.0, "alpha_star": 10.0, "dt": 1.0e-5, "steps": 5000,
        "source": "Table 3, M2 (a33, a34, a44, b3, b4)",
    },
    "fritsch_4sp_M3": {
        "n": 4, "species": [1, 2, 3, 4], "A": A4_TRUE, "b": B4_TRUE,
        "eta": [1.0, 1.0, 1.0, 1.0], "phi_init": [0.02, 0.02, 0.02, 0.02],
        "c_star": 25.0, "alpha_star": 0.0, "dt": 1.0e-4, "steps": 750,
        "source": "Table 3, M3 (a13, a14, a23, a24; diagonal blocks from M1/M2). "
                  "Likely origin of the repo's c* = 25.",
    },
    "fritsch_4sp_val": {
        "n": 4, "species": [1, 2, 3, 4], "A": A4_TRUE, "b": B4_TRUE,
        "eta": [1.0, 1.0, 1.0, 1.0], "phi_init": [0.02, 0.02, 0.02, 0.02],
        "c_star": 25.0, "alpha_star": _alpha_validation, "dt": 1.0e-4, "steps": 1500,
        "source": "Table 3, M_val^3 (validation), switch as printed (step 500)",
    },
    "fritsch_4sp_val_fig18": {
        "n": 4, "species": [1, 2, 3, 4], "A": A4_TRUE, "b": B4_TRUE,
        "eta": [1.0, 1.0, 1.0, 1.0], "phi_init": [0.02, 0.02, 0.02, 0.02],
        "c_star": 25.0, "alpha_star": _alpha_validation_fig18, "dt": 1.0e-4, "steps": 1500,
        "source": "Fig. 18 reading of M_val^3: switch at step 750",
    },
}
# prior for every Case II parameter: U(0, 3); CoV = 0.5 % throughout (Table 3 text)
