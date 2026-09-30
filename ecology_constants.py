"""ecology_constants.py -- the constants the Hamilton ecology model is run with,
in one place, for every path that solves it: the 0D integrator
(jax_hamilton_0d_5species_demo.py, ecology_jax.py -> the ANSYS bridge) and the
1D / 2D Hamilton+nutrient PDEs (JAXFEM/core_hamilton_*).

C_STAR: the nutrient constant c* multiplying A (Klempt et al. 2026, Eq. 10).
    25 is what the TMCMC calibration of A used (Tmcmc202601
    data_5species/model_config/model_constants.json and the estimation
    scripts; the manuscript, Sec. 2.3). c* enters only through c*A, so any
    other value rescales every calibrated interaction. In the PDEs the
    nutrient field c(x) is normalised to 1 at the supply boundary, and the
    reaction sees C_STAR * c(x).
    (Klempt et al.'s own numerical examples use c* = 100; the reproduction
    of those, JAXFEM/klempt2026_reproduction.py, sets c* per case and does
    not use this constant.)

K_HILL, N_HILL: the Hill gate on species 5 (a repo extension, not part of
    Klempt et al.). Off (K_HILL = 0), by decision 2026-09-29. N_HILL only
    matters if the gate is switched on; 2 is the TMCMC estimation value.

ALPHA_STAR: the antibiotic concentration alpha* (Klempt et al. 2026, Eq. 17).
    0, by decision 2026-09-29: the TMCMC calibration ran with alpha* = 0 and
    its 15-parameter A has no antibiotic coefficients b, so the b entries of
    the repo's 20-parameter theta vectors have no effect.
"""

C_STAR = 25.0
K_HILL = 0.0
N_HILL = 2.0
ALPHA_STAR = 0.0
