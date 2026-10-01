#!/usr/bin/env python3
"""Reference for the one-species coupling, independent of the Fortran.

Klempt et al. (2024) Eq. 36, ``alpha_dot = k_alpha * phi``, integrated the way
``growth_from_phi.f`` integrates it (explicit in phi). Used by
``tests/test_one_species_coupling.py`` and, once the partner-side deck runs on
IKMHIWI03, to predict ``SVAR(84)`` from the phi history that deck prints.

Kept deliberately small: if this file grows a cap, a clamp or a gate, it has
stopped being Klempt 2024.
"""
from __future__ import annotations

import numpy as np


def alpha_history(phi, k_alpha: float, dt: float) -> np.ndarray:
    """alpha after each step, from alpha(0) = 0, for a sequence of phi."""
    alpha = 0.0
    out = []
    for p in np.asarray(phi, dtype=float):
        alpha = alpha + k_alpha * p * dt
        out.append(alpha)
    return np.array(out)


def alpha_klempt(alpha_repo: float) -> float:
    """This repository's alpha in Klempt's convention, Fg = alpha_K * I."""
    return 1.0 + alpha_repo
