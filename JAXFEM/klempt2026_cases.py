"""klempt2026_cases.py -- the numerical examples of Klempt, Geisler, Soleimani,
Junker, "A continuum multi-species biofilm model with a novel interaction
scheme" (arXiv:2509.01274v1; published as Klempt2026ContinuumBacterialGrowth),
transcribed as data so this repo's own 0D Hamilton integrator can be run on
exactly the paper's inputs and compared with the paper's figures.

Transcribed from the PDF directly (not a summary): Table 1 (two species,
cases 1-6), Table 3 / Appendix A (two species, cases 4*, 5*, 4**, 5**),
Eq. 19-20 and Table 2 (four species, cases 1-4). Common to all:
dt = 1e-4 Time, Kp_i = eta_i * 1e-4, material point, implicit Newton.
The paper does not state the initial psi or gamma; PSI0 below is this
repo's own default and is a free choice here, not a paper value.

Two places where the printed equations and the paper's own figures disagree,
both settled by running this repo's integrator on every case (2026-09-29;
the published version, Arch Appl Mech 96:164, prints both unchanged):
  * Eq. 17 ends in "+gamma". With it, even the 2-species figures are not
    reproduced (case 3 phibar 0.22/0.19 vs 0.585/0.385; case 6 the wrong
    species wins). Without it -- the form a variation w.r.t. psi gives, and
    the form this repo's residual has -- all ten 2-species cases match.
  * Eq. 19 prints A with a prefactor 1/2. The 4-species figures are
    reproduced with the integer matrix itself (A_FIGURE_SCALE = 2 below):
    e.g. case 1 phibar 0.410/0.248/0.168/0.151 vs ~0.41/0.25/0.17/0.15,
    phi0 reaching 0 at step 287 vs ~290. With the 1/2 the small-phi cases
    never take off within the plotted window. Table 1's 2-species A are
    integers with no prefactor, consistent with this.

`paper` holds values read off the paper's figures by eye -- the paper gives
no tables of results. They are approximate (about +-0.02 in value, +-30 in
time step) and are there to be compared against, not to be fitted to.
Keys: phibar_end = living fraction phi_i*psi_i at the last plotted step;
phi_end / psi_end likewise; phi0_zero_step = step at which the empty space
phi0 first reaches ~0 (the kink in the curves).
"""
from __future__ import annotations

import math

DT = 1.0e-4
PSI0 = 0.999

_A4 = [[0.5, 2.5, 2.5, 2.5],
       [2.5, 0.5, 1.5, 1.5],
       [2.5, 1.5, 0.5, 1.0],
       [2.5, 1.5, 1.0, 0.5]]           # Eq. 19: (1/2) * [[1,5,5,5],[5,1,3,3],[5,3,1,2],[5,3,2,1]]
_ETA4 = [0.8, 1.0, 1.5, 2.0]           # Eq. 20
A_FIGURE_SCALE_4SP = 2.0               # figures were computed without Eq. 19's 1/2 -- see module docstring


def _two(a11, a12, a22, b1, b2, eta1, eta2, phi1, phi2, steps, paper, fig):
    return {
        "n": 2, "A": [[a11, a12], [a12, a22]], "b": [b1, b2],
        "eta": [eta1, eta2], "phi_init": [phi1, phi2], "steps": steps,
        "c_star": 100.0, "alpha_star": 10.0, "figure": fig, "paper": paper,
    }


def _four(b, phi, steps, c_star, alpha_star, paper, fig):
    return {
        "n": 4, "A": _A4, "A_figure_scale": A_FIGURE_SCALE_4SP, "b": b, "eta": _ETA4, "phi_init": phi,
        "steps": steps, "c_star": c_star, "alpha_star": alpha_star,
        "figure": fig, "paper": paper,
    }


def _c_sine(step):
    return 50.0 + 50.0 * math.sin(500.0 * step * DT)      # Table 2 case 3; t in Time (Fig. 9a: ~126-step period)


def _alpha_step500(step):
    return 100.0 if step > 500 else 0.0                    # Table 2 case 4; switch at step 500 (Fig. 9b)


CASES = {
    # ---- two species, Table 1 -------------------------------------------------
    "2sp_case1": _two(2, 0, 1, 0, 0, 1, 1, 0.2, 0.2, 500, fig="Fig. 1", paper={
        "phi_end": [0.97, 0.02], "psi_end": [1.0, 0.70], "phi0_zero_step": 240}),
    "2sp_case2": _two(1, 0, 1, 0, 0, 1, 2, 0.2, 0.2, 1500, fig="Fig. 2", paper={
        "phi_end": [0.97, 0.02], "psi_end": [1.0, 0.70], "phi0_zero_step": 420}),
    "2sp_case3": _two(1, 1, 1, 0, 0, 1, 2, 0.2, 0.2, 500, fig="Fig. 3", paper={
        "phi_end": [0.60, 0.40], "phibar_end": [0.585, 0.385], "phi0_zero_step": 220}),
    "2sp_case4": _two(1, 0, 1, 1, 2, 1, 2, 0.2, 0.3, 1500, fig="Fig. 4", paper={
        "phibar_end": [0.0, 0.95], "psi_end": [0.07, 0.97], "phi0_zero_step": 1080}),
    "2sp_case5": _two(1, 0, 1, 1, 2, 1, 2, 0.25, 0.3, 1500, fig="Fig. 5", paper={
        "phibar_end": [0.95, 0.0], "psi_end": [0.97, 0.07], "phi0_zero_step": 450}),
    "2sp_case6": _two(1, -1, 1, 0, 0, 1, 2, 0.2, 0.2, 1500, fig="Fig. 6", paper={
        "phibar_end": [0.0, 0.95], "psi_end": [0.05, 0.97], "phi0_zero_step": 1050}),
    # ---- two species, Table 3 / Appendix A ----------------------------------------
    "2sp_case4s": _two(1, 0, 1, 1, 2, 2, 1, 0.2, 0.3, 1500, fig="Fig. 12", paper={
        "phibar_end": [0.0, 0.95], "phi0_zero_step": 1100}),
    # 5s/5ss phi0 re-read on the published figures: phi0 reaches 0 where phibar_1 plateaus
    "2sp_case5s": _two(1, 0, 1, 1, 2, 2, 1, 0.25, 0.3, 1500, fig="Fig. 13", paper={
        "phibar_end": [0.95, 0.0], "phi0_zero_step": 900}),
    "2sp_case4ss": _two(1, 0, 1, 1, 2, 1, 1, 0.2, 0.3, 1500, fig="Fig. 14", paper={
        "phibar_end": [0.0, 0.95], "phi0_zero_step": 800}),
    "2sp_case5ss": _two(1, 0, 1, 1, 2, 1, 1, 0.25, 0.3, 1500, fig="Fig. 15", paper={
        "phibar_end": [0.95, 0.0], "phi0_zero_step": 580}),
    # ---- four species, Eq. 19-20 + Table 2 ---------------------------------------
    "4sp_case1": _four([0.4, 0.3, 0.2, 0.1], [0.02, 0.02, 0.02, 0.02], 500, 100.0, 10.0,
                       fig="Fig. 7", paper={"phibar_end": [0.41, 0.25, 0.17, 0.15],
                                            "phi0_zero_step": 290}),
    "4sp_case2": _four([0.4, 0.3, 0.2, 0.1], [0.02, 0.02, 0.02, 0.2], 500, 100.0, 10.0,
                       fig="Fig. 8", paper={"phibar_end": [0.41, 0.24, 0.13, 0.18],
                                            "phi0_zero_step": 70}),
    "4sp_case3": _four([0.4, 0.3, 0.2, 0.1], [0.02, 0.02, 0.02, 0.02], 1500, _c_sine, 10.0,
                       fig="Fig. 10", paper={"phibar_end": [0.40, 0.24, 0.16, 0.15],
                                             "phi0_zero_step": 580}),
    "4sp_case4": _four([10.0, 2.0, 1.0, 0.01], [0.02, 0.02, 0.02, 0.02], 1500, 100.0, _alpha_step500,
                       fig="Fig. 11", paper={"phibar_end": [0.0, 0.0, 0.0, 0.93],
                                             "phi0_zero_step": 90}),   # steep drop, readable only as ~60-120
}
