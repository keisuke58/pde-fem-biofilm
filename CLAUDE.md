# CLAUDE.md

Companion code for a LUH/IKM master's thesis: 3D FEM stress analysis of oral
biofilms on tooth/implant geometry, built on the Klempt (2024) continuum
growth model. Pipeline: CLSM composition → TMCMC-calibrated 5-species ecology
→ JAXFEM PDE growth field α(x) → Abaqus/ANSYS UMAT (`Fg=(1+α)I`) → von Mises
stress comparison across conditions. See `README.md` and `REPO_MAP.md` for
the full guided tour; `RESEARCH_MODEL.md` for the modeling details.

This machine (host `IKMHIWI03`) is the primary **ANSYS** work environment —
prefer ANSYS/APDL workflows here over Abaqus when both are viable for a task.
**Abaqus 2024 is also actually installed and licensed here** (confirmed
2026-08-20: `C:\SIMULIA\Commands\abaqus.bat`, `abaqus information=release`
completes with a valid Site ID) — earlier session notes assumed otherwise
and were wrong. What genuinely is NOT on this machine is any *prior Abaqus
run output* (no `.odb`/`.sta`/`.msg`/`.dat` anywhere on `C:`), so a fresh
Abaqus run is possible here but nothing has been run here yet.

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
- **Never add a `Co-Authored-By: Claude` (or similar AI-attribution) trailer
  to commit messages, and don't otherwise mark commits as AI-assisted.**
  Commits should read as the user's own work (git identity is already set
  locally to Keisuke Nishioka <kei128608@gmail.com> — see `.git/config`).
  Verified 2026-08-20: none of the ~20 commits made this session carry any
  such trailer, and there's no commit template/hook in this repo that would
  add one — keep it that way.
- **Incident, 2026-08-20: 46 commits (2026-07-02 to 2026-08-19) had author
  `Claude <noreply@anthropic.com>`**, visible on GitHub's Contributors page
  — traced to a Claude Code environment other than this machine (no global
  `.gitconfig` exists here, so it wasn't this machine's default; likely a
  cloud/remote session). Fixed via `git filter-branch` (author/committer
  rewritten to Keisuke Nishioka, content byte-identical, verified) + force
  push — see the `git-history-rewrite-2026-08-20` memory for the full
  incident and what it means for other clones (e.g. the Keio server one
  needs re-cloning, not pulling). A local `pre-commit` hook now blocks any
  commit with an `anthropic.com`/`Claude`-looking identity on this clone;
  install it in any fresh clone with:
  `cp scripts/pre-commit-no-ai-identity.sh .git/hooks/pre-commit && chmod +x .git/hooks/pre-commit`
  — this only protects clones where it's installed, not every environment.

## This PC vs. claude.ai (web) — don't mix them up

The user also discusses this repo with Claude on claude.ai (browser, no file/
tool access). Ground rules so nothing said there gets mistaken for verified
fact here, or vice versa:

- **claude.ai has no access to this repo, ANSYS, or git.** Anything it says
  about specific files, line numbers, current test/build status, or "what the
  code currently does" is inference from whatever was pasted into that chat —
  not a live read of the repository. Treat it as a source of ideas/drafts to
  bring back here and verify, never as a substitute for actually checking.
- **This machine (Claude Code) is the only place that can confirm anything** —
  build success, test results, ANSYS output, git state. If a claude.ai
  conversation concluded something works, re-verify it here before relying on
  it (see the 2026-08-20 incident: an unverified change to
  `reference_values.json` silently broke 5 tests for a while).
- **Never paste `.env`, the GitHub PAT, or other secrets into claude.ai.**
  This PC's push workflow already isolates the token to local PowerShell
  calls with redacted output — keep it that way.
- If the user brings a plan or code snippet over from a claude.ai chat,
  treat file paths/API shapes/current-state claims in it as unverified until
  checked against the actual repo, the same as any other secondhand claim.
