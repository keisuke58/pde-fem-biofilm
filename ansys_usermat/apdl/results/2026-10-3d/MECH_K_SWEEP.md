# Mechanical term of Klempt 2024 Eq. 30: K sweep and checks (Nov 2026)

Term added to Eq. 34: $-K\,\phi_{3D}\,(\mathbf I:\mathbf C_e-3)$ with $K=\mu/\eta_\phi$ ($\eta_\phi$ is not given in the paper, so K is
swept); `coupled_living_consumption.py --mech-feedback K`. The stress is solved on the 16^3 stress mesh at the start of each
outer step. Not part of the LUH thesis (the model there is frozen); Keio continuation, `KEIO_PLAN.ja.md` 0b'.

## 1. Sweep on 64^3, 128^3, 192^3 (fifa, part 12; table in `KEIO_CHAIN_LOG.md`, 8-10 Nov)

Seed means at T* = 1 against K = 0 on the same grid, superbee, `--transient-c 1e8`, outer step 0.025, with antibiotic.
Figure: `assets/fig_mech_K_sweep.png` (`fig_mech_K_sweep.py`).

| K [1/T*] | separate: alpha - 1, von Mises | mixed: alpha - 1, von Mises |
|---|---|---|
| 1e3 | -0.01 %, -0.04 % | -0.1 %, +0.1 % |
| 1e4 | -0.2 %, -0.5 % | -1.1 %, +1.3 % |
| 3e4 | -0.7 %, -2.4 % | -3.4 %, +0.4 % |
| 1e5 | -3.7 %, -20 % | -9.0 %, -19 % |
| 3e5 | diverged | -19 %, -68 % |
| 1e6 | diverged | diverged |

(finest grid available: 192^3 for 1e4 and 1e5, 128^3 otherwise.)

- The term changes alpha - 1 or the stress by 1 % or more from K = 1e4 (mixed) and 3e4 (separate). Below 1e3 it does nothing,
  which is the range in which Klempt et al. dropped it.
- At K = 1e5 the change of alpha - 1 is converged on the grids (GCI 0.1 % separate, 0.8 % mixed); the stress has GCI 4.6 %
  (mixed) and 24 % (separate, observed order 0.25).
- Also in this free example phi_3D exceeds 1 where the biofilm is compressed: maximum 1.02 at K = 1e4, 1.33 at 1e5, 3.0 in the
  mixed run at 3e5 (`phi_max` in the JSON files).

## 2. Checks on 32^3 (cloud, 13 Nov), separate example with antibiotic

Changes against K = 0 at the same outer step (K = 0 itself: von Mises +2.1 % and +4.3 % at half and quarter step, alpha - 1
-0.4 % and -0.5 %).

| K | variant | outer step | alpha - 1 | von Mises | max phi_3D |
|---|---|---|---|---|---|
| 1e5 | as on fifa | 0.025 | -4.6 % | -25.6 % | 1.32 |
| 1e5 | as on fifa | 0.0125 | -4.6 % | -23.1 % | 1.33 |
| 1e5 | as on fifa | 0.00625 | -4.5 % | -20.9 % | 1.33 |
| 1e5 | `--mech-exp` (term integrated exactly) | 0.025, 0.0125 | as above | as above (4 digits) | 1.32, 1.33 |
| 1e5 | `--mech-refresh 4` (stress solved 4 times per outer step) | 0.025 | -4.8 % | -27.6 % | 1.33 |
| 1e5 | `--mech-refresh 4` | 0.0125 | -4.7 % | -24.0 % | 1.33 |
| 1e5 | `--phi-clip` (phi_3D kept in [0, 1]; a change of the model) | 0.025 | -6.0 % | -34.5 % | 1.00 |
| 1e5 | `--phi-clip` | 0.0125 | -5.9 % | -32.0 % | 1.00 |
| 1e5 | `--phi-clip` | 0.00625 | -5.9 % | -30.0 % | 1.00 |
| 3e5 | as on fifa | 0.025 | NaN | NaN | - |
| 3e5 | as on fifa | 0.00625 | -18.7 % | -94 % | 477 |
| 3e5 | `--mech-exp` | 0.025 | -13.2 % | -59 % | 6.0 |
| 3e5 | `--mech-exp` | 0.00625 | -18.6 % | -93 % | 11 |
| 3e5 | `--mech-exp --phi-clip` | 0.025 | -13.2 % | -20 % | 1.00 |
| 3e5 | `--mech-exp --phi-clip` | 0.0125 | -14.7 % | -49 % | 1.00 |

What this shows:

1. **alpha - 1 at K = 1e5 is robust**: -4.5 to -4.6 % at every outer step (32^3), -3.7 % on 192^3.
2. **The stress at K = 1e5 is not converged in the outer step**: the reduction shrinks by about 2.3 points per halving
   (-25.6, -23.1, -20.9 %), with no sign of a limit yet. It is not the lag of the term: solving the stress four times per outer
   step moves the value the other way by 2 points, and the exact integration of the term changes nothing. So the "-20 %" of the
   sweep is an estimate with an uncertainty of several points from the time step alone.
3. **The divergence at K >= 3e5 is the term, not the scheme.** With exact integration or a small step the runs finish, but
   phi_3D grows to 6 to 477 where the biofilm is compressed, and the results depend strongly on the step. The elastic energy
   $(\mu/2)\,\phi^2\,(\mathbf I:\mathbf C_e-3)$ has no lower bound in $\phi$ under compression; the penalty that keeps phi_3D in [0, 1]
   (5 /T*) is weaker than the term. This is the same mechanism as in the rigid-bars example (`KEIO_PLAN.ja.md` 0b').
4. **Clipping phi_3D to [0, 1] removes the runaway but changes the answer** (K = 1e5: stress -30 to -35 % instead of -21 to
   -26 %) and at 3e5 the result still depends strongly on the step. It is a change of the model, not a fix.

## 3. Reading

- Up to K of about 3e4 the term is a small, well-behaved correction (below 3.5 % in alpha - 1, below 3 % in the stress).
- Around K = 1e5 the growth result is usable (alpha - 1 lower by 4 to 9 %), the stress result only as an order (lower by about a
  fifth), and phi_3D already leaves [0, 1].
- Beyond that the term as printed is not usable without a bound on phi_3D or a different energy. This and the sign of the
  effect under confinement (stress up by a factor of 4 between rigid bars, against experiments in which compression slows
  growth) are the points to discuss with Felix Klempt before building on the term.

New option (13 Nov): `--mech-refresh N` solves the stress N times per outer step (tag `rfN`); N = 1 is the old behaviour.
Files: `coupled_living_consumption_front1_[dtc..._]..._mech<K>[exp][rfN][_clip]_transc1e+08_tvd_superbee_wgmirror_fast_separate_ab.json`.
