# CLAUDE.md

Companion code for a LUH/IKM master's thesis: 3D FEM stress analysis of oral
biofilms on tooth/implant geometry, built on the Klempt (2024) continuum
growth model. Pipeline: CLSM composition → TMCMC-calibrated 5-species ecology
→ JAXFEM PDE growth field α(x) → Abaqus/ANSYS UMAT (`Fg=(1+α)I`) → von Mises
stress comparison across conditions. See `README.md` and `REPO_MAP.md` for
the full guided tour; `RESEARCH_MODEL.md` for the modeling details.

This machine (host `IKMHIWI03`) is the primary **ANSYS** work environment —
prefer ANSYS/APDL workflows here over Abaqus when both are viable for a task.

## Key directories

- `ansys_usermat/` — ANSYS USERMAT (Fortran) port of the Klempt growth model,
  cross-checked against the Abaqus UMAT. `usermat_biofilm.f` is the core;
  `crosscheck/` compares ANSYS vs Abaqus outputs; `coupling/` is a
  Python-material-server coupling shim. `apdl/` (closed-form growth
  verification deck + RUNBOOK) lands via PR #29 — not yet pulled into this
  working tree as of 2026-08-19.
- `JAXFEM/` — JAX PDE solver for the growth/composition field, TMCMC
  calibration, posterior propagation.
- `ch5_flow/`, `umat_flow/` — thesis chapter LaTeX + associated flow/UMAT docs.
- `tier2b_real/`, `configs/`, `runs/` — Abaqus coupon/implant job generation,
  configs, and run logs.
- `tests/` — pytest unit tests (`pytest tests/`).

## ANSYS environment on this PC

Full hardware/license/product inventory: `ANSYS_ENVIRONMENT.md`. Summary:

- **ANSYS 2022 R2 (v222)** only, at `C:\Program Files\ANSYS Inc\v222`.
  Env vars `AWP_ROOT222`, `ANSYS222_DIR` already set system-wide.
- License: floating, via RRZN Uni Hannover server
  (`1055@ansys-lic.rrzn.uni-hannover.de` / `2325@...` for ANSYSLI) — needs
  campus network or VPN to check out.
- Custom UMAT build location: `C:\Program Files\ANSYS Inc\v222\ansys\custom\user\winx64\` —
  **not yet verified**; the RUNBOOK's `ANSCUST.bat` build steps are
  doc-derived, not confirmed against what's actually on this machine.
- **Intel Fortran (ifort) / Visual Studio presence: unconfirmed.** Not found
  via plain `where ifort`; must check from the "Intel oneAPI command prompt
  for Intel 64 for Visual Studio" Start Menu entry, not a bare cmd/PowerShell.
  Do this before assuming the USERMAT build will work.

## Git on this machine — important quirks

- No `git` on PATH. Use the MSYS64 install directly:
  `C:\msys64\usr\bin\git.exe`, and prepend `C:\msys64\usr\bin` to `$env:Path`
  for that PowerShell call — otherwise git's https/credential helper
  subprocesses fail with "shared libraries" errors (they need MSYS DLLs on
  PATH, not just the git.exe path).
- **This repo's working tree has a massive line-ending mismatch** — `git
  status` shows ~480+ files as modified with equal insertions/deletions
  (pure CRLF↔LF churn, zero real content change). **Never `git add -A` or
  `git commit -a`.** Always stage specific files by name.
- `origin` is `https://github.com/keisuke58/pde-fem-biofilm.git`. No
  credential helper was configured as of 2026-08-19 — pushes prompt for
  username/password (PAT) on the console. If a push hangs, it's waiting on
  that prompt.

## Working style for this repo

- Keep changes scoped to named files; don't touch the pre-existing
  line-ending noise even incidentally.
- Prefer direct edits over spawning subagents for small, well-scoped tasks —
  this user is cost-conscious about agent/token usage.
