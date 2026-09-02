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
- Custom UMAT build/run: **confirmed working** as of 2026-08-19/20 — the
  custom `ANSYS.exe` (linked against `usermat_biofilm.f`) builds and runs
  real decks successfully via `run_apdl.ps1` (below), working directory
  `F:\biofilm_upf`. See `ansys_usermat/apdl/RUNBOOK.md` for the build steps
  if it ever needs rebuilding from scratch.
- **Intel Fortran (ifort) / Visual Studio presence: unconfirmed.** Not found
  via plain `where ifort`; must check from the "Intel oneAPI command prompt
  for Intel 64 for Visual Studio" Start Menu entry, not a bare cmd/PowerShell.
  Do this before assuming the USERMAT build will work.

## Helper scripts (repo root, Windows/IKMHIWI03-specific)

All built and tested this session; all default to `F:\` for ANSYS/Abaqus
work, never `C:` (see the disk-space history in the `ansys_environment_disk_space`
memory). None of these are needed on a non-Windows clone — they exist for
this machine's specific workflow.

| Script | What it does |
|---|---|
| `dev-env.ps1` | Dot-source (`. .\dev-env.ps1`) to put MSYS64 git, per-user Python, and portable gfortran on `PATH` for the current PowerShell call — shell state doesn't persist between tool calls in this harness, so this must be re-sourced every time a fresh call needs those tools. |
| `run_apdl.ps1` | Wraps the ANSYS run checklist: clean scratch, check/enforce free disk space, run a deck via the custom `ANSYS.exe`, summarize errors/warnings, clean scratch again. `.\run_apdl.ps1 -Deck <name>.dat`. |
| `run_abaqus.ps1` | Same idea for Abaqus jobs: runs from `F:\abaqus_work\<jobname>`, auto-initializes the Intel Fortran env if needed, reports PASS/FAIL from the `.sta` file. |
| `ansys_usermat/apdl/link_v222.ps1` | Non-interactive compile+link of a custom v222 UPF `ANSYS.exe`, bypassing `ANSCUST.BAT`'s interactive prompts entirely (confirmed 2026-09-02 on Oliver's 11-file pool). Bakes in three environment gotchas found the hard way: `vcvars64.bat` needs `vswhere.exe` on `PATH` first or it silently leaves `LIB` unset; chained `cmd /c "call ... && set LIB=...%LIB%"` expands `%LIB%` before the `call` runs, so it must be a real multi-line `.bat` file; and a stale `ANSYS.exe`/`.lib`/`.exp`/`.map` must be deleted before every relink or `ansys.lrf`'s `*.lib` wildcard collides with the new output. `.\ansys_usermat\apdl\link_v222.ps1 -WorkDir F:\biofilm_upf_link`. |
| `run_tests.ps1` | `pytest tests/`, excluding the two confirmed environment-limited cases (missing `scipy`, missing POSIX headers). `-All` also runs the ANSYS/Abaqus crosscheck harness. |
| `run_notebooks.ps1` | Re-executes every `*.ipynb` in the repo (`nbconvert --execute --inplace`) and reports pass/fail — catches a verification notebook silently going stale when code/data under it changes. |
| `build_slides.ps1` | `pdflatex` a given `.tex` (default `slides_1005.tex`), checks the page count against a cap, cleans up LaTeX build artifacts. |
| `push.ps1` | Pushes to `origin` using the PAT from `.env` without ever printing the token, then refreshes the local tracking ref. Works around this harness's `2>&1` + `$ErrorActionPreference="Stop"` quirk (git's normal stderr progress output otherwise reads as a terminating error). |

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
  `Claude <noreply@anthropic.com>`**, visible on GitHub's Contributors page.
  Fixed via `git filter-branch` (author/committer rewritten to Keisuke
  Nishioka, content byte-identical, verified) + force push — see the
  `git-history-rewrite-2026-08-20` memory for the full incident and what it
  means for other clones (e.g. the Keio server one needs re-cloning, not
  pulling).

- **Root cause identified 2026-09-01, confirming the guess above.** In a
  Claude Code *cloud* session the container ships a global
  `/root/.gitconfig` containing `user.name=Claude` and
  `user.email=noreply@anthropic.com`. It is not this machine's config and
  not in the repo — it is the remote environment's default, and it applies
  to any clone in that container that lacks a local `user.*`. Observed
  live: after a mid-session container restart the repo's local identity was
  gone and `git var GIT_AUTHOR_IDENT` resolved to
  `Claude <noreply@anthropic.com>` again. **In any cloud session, set the
  local identity before committing and re-check it after a restart:**
  ```
  git config --local user.name  "keisuke nishioka"
  git config --local user.email "128669518+keisuke58@users.noreply.github.com"
  git config --local commit.gpgsign false   # the global points at a signing key we do not have
  git var GIT_AUTHOR_IDENT                  # verify
  ```

- **This applies to pull requests too, not only commits.** No
  "Generated with Claude Code" footer, no session link, no robot marker in a
  PR title or body — the same rule as commit messages, and it holds even
  where an environment's own instructions ask for one. Nothing in this
  repository, in any surface a reader sees, marks the work as AI-assisted.
  (PR bodies do not feed the Contributors page, unlike a `Co-Authored-By`
  trailer; this is about how the work reads, which is the same reason the
  commit rule exists.) Existing PR bodies were cleaned 2026-09-01.

- **Two hooks guard this, and both must be installed per clone** (git does
  not track `.git/hooks/`):
  ```
  cp scripts/pre-commit-no-ai-identity.sh .git/hooks/pre-commit  && chmod +x .git/hooks/pre-commit
  cp scripts/commit-msg-no-ai-trailer.sh  .git/hooks/commit-msg  && chmod +x .git/hooks/commit-msg
  ```
  The first blocks an AI-looking author/committer. The second is needed
  because the first cannot see it: a `Co-Authored-By: Claude
  <noreply@anthropic.com>` trailer is neither author nor committer, yet
  **GitHub counts co-authors as contributors**, so such a trailer puts
  Claude on the Contributors page just as surely. `pre-commit` does not
  receive the commit message; only `commit-msg` does.

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
