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

## Which machine is this?

- **IKMHIWI03** (Windows, LUH): everything below about ANSYS, `F:\`, PowerShell
  scripts and Portable Git applies. Available until mid-December 2026.
- **The Keio Linux server (fifa)**: Abaqus and Python only, no ANSYS (for about
  half a year). After cloning, run `bash scripts/setup_keio_server.sh` (git
  identity, hooks, checks), then read `KEIO_SERVER_HANDOFF.ja.md` and
  `KEIO_PLAN.ja.md` §0. The Windows sections below do not apply there.
- **A Claude Code cloud session**: see the git-identity note under "Working
  style"; no ANSYS or Abaqus.

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
- `references/` — Klempt et al. 2024 (BMMB), the paper this work follows, and
  Soleimani et al. 2023 (Sci. Rep., two-species co-aggregation in an ANSYS user
  element) and Feng et al. 2021 (Bull. Math. Biol., two-species oral biofilm,
  spreading driven by the summed species growth), each with
  a searchable text extraction. Licence and credit in `THIRD_PARTY.md`.
- `ecology_constants.py` — the one place the Hamilton ecology model's c*
  (25, the TMCMC calibration value) and Hill gate (off) are set. Every path
  (0D / ANSYS bridge `ecology_jax`, 1D and 2D PDEs) imports it; never
  hard-code `c`, `K_hill`, `n_hill` again (decision 2026-09-29). Exception:
  reproductions of Klempt et al.'s own examples
  (`JAXFEM/klempt2026_reproduction.py`, c* = 100 per case).

## ANSYS environment on this PC

Full hardware/license/product inventory: `ANSYS_ENVIRONMENT.md`. Summary:

- **ANSYS 2022 R2 (v222)** only, at `C:\Program Files\ANSYS Inc\v222`.
  Env vars `AWP_ROOT222`, `ANSYS222_DIR` already set system-wide.
- License: floating, via RRZN Uni Hannover server
  (`1055@ansys-lic.rrzn.uni-hannover.de` / `2325@...` for ANSYSLI) — needs
  campus network or VPN to check out.
- Custom UMAT build/run: **confirmed working**, re-verified 2026-09-03 — the
  custom `ANSYS.exe` (linked against `usermat_biofilm.f`) at `F:\biofilm_upf`
  runs real decks and reproduces the closed-form reference exactly
  (`t_growth_constrained.dat`: SX=SY=SZ=−1.0193e−04, SVAR(10)=0.05, 0 errors).
  **There is no `run_apdl.ps1` in this repo** — despite an earlier version of
  this doc claiming one — `git log --all -- run_apdl.ps1` finds no commit
  that ever added it. Invoke ANSYS directly, from a writable working
  directory (`F:\biofilm_upf`, never `C:`):
  ```bat
  "%AWP_ROOT222%\ANSYS\bin\winx64\ANSYS222.exe" -b -custom .\ANSYS.exe ^
      -i <deck>.dat -o out.txt
  ```
  `-custom` is essential — without it ANSYS runs its own stock material and
  the run means nothing. See `ansys_usermat/apdl/RUNBOOK.md` for the full
  build/run procedure and `link_v222.ps1` (below) for rebuilding from
  scratch.
- **Which `ANSYS.exe` (2026-09-29):** `F:\biofilm_upf\ANSYS.exe` is the
  8/19 build (inline core only). It gives **non-deterministic, broken
  solves under `-smp -np` > 1** — use it only with `-np 1` or the default
  DMP launch. The current build is `F:\biofilm_upf_kusepy\ANSYS.exe`
  (usermat_biofilm.f + Python/ecology bridge with n_sub, via
  `link_v222.ps1`), safe with `-smp -np 4/8`. When relinking,
  `usermat_py_hook.f` must be compiled **before** `usermat_biofilm.f`, or
  ifort reads a stale `biofilm_py_bridge.mod` and link_v222.ps1 still
  links the old object into a mixed exe.
- **Native point model (6 Oct 2026):** `F:\biofilm_upf_native\ANSYS.exe` is
  the partner-element build with `ansys_usermat/coupling/ecology_native.f`
  (the point model in Fortran) in place of `usermat_py_hook.f`: no material
  server, no port 8765; set `BIOFILM_ECO_CASE` to a file from
  `abaqus_composition/write_eco_cfg.py <case>` before starting ANSYS (-np 1).
  Same results as the server build (abaqus_composition/README.md), ~4x faster
  on 8^3. Abaqus: `run_comp.ps1 -Native`. The week chain still uses
  `F:\biofilm_upf_wired` with the server.
- **Intel Fortran (ifort) / Visual Studio: confirmed working** through
  `link_v222.ps1` (ifort 2025.3 + VS 18, 2026-09-02/09-29). A bare
  `where ifort` still finds nothing — the script sets up the environment
  itself.

## Helper scripts (repo root, Windows/IKMHIWI03-specific)

All built and tested this session; all default to `F:\` for ANSYS/Abaqus
work, never `C:` (see the disk-space history in the `ansys_environment_disk_space`
memory). None of these are needed on a non-Windows clone — they exist for
this machine's specific workflow.

| Script | What it does |
|---|---|
| `dev-env.ps1` | Dot-source (`. .\dev-env.ps1`) to put MSYS64 git, per-user Python, and portable gfortran on `PATH` for the current PowerShell call — shell state doesn't persist between tool calls in this harness, so this must be re-sourced every time a fresh call needs those tools. |
| `run_abaqus.ps1` | Runs Abaqus jobs from `F:\abaqus_work\<jobname>`, auto-initializes the Intel Fortran env if needed, reports PASS/FAIL from the `.sta` file. |
| `ansys_usermat/apdl/link_v222.ps1` | Non-interactive compile+link of a custom v222 UPF `ANSYS.exe`, bypassing `ANSCUST.BAT`'s interactive prompts entirely (confirmed 2026-09-02 on Oliver's 11-file pool). Bakes in three environment gotchas found the hard way: `vcvars64.bat` needs `vswhere.exe` on `PATH` first or it silently leaves `LIB` unset; chained `cmd /c "call ... && set LIB=...%LIB%"` expands `%LIB%` before the `call` runs, so it must be a real multi-line `.bat` file; and a stale `ANSYS.exe`/`.lib`/`.exp`/`.map` must be deleted before every relink or `ansys.lrf`'s `*.lib` wildcard collides with the new output. `.\ansys_usermat\apdl\link_v222.ps1 -WorkDir F:\biofilm_upf_link`. |
| `run_tests.ps1` | `pytest tests/`, excluding the two confirmed environment-limited cases (missing `scipy`, missing POSIX headers). `-All` also runs the ANSYS/Abaqus crosscheck harness. |
| `run_notebooks.ps1` | Re-executes every `*.ipynb` in the repo (`nbconvert --execute --inplace`) and reports pass/fail — catches a verification notebook silently going stale when code/data under it changes. |
| `build_slides.ps1` | Builds decks twice (pdflatex, or lualatex for the Japanese deck) and reads the log for `!` errors, frames that run off the page (`Overfull \vbox`), undefined references and the page count; cleans up. `-All` = the three meeting decks, `-Chapter` = `thesis_ch5/_build_check.tex`. Exit 1 on any problem. |
| `commit.ps1` | `-Files a,b -Message $m [-Push]`: stages exactly the named files, refuses if anything else is staged or the message has an AI co-author trailer, writes the message as UTF-8 without BOM. |
| `push.ps1` | Pushes HEAD to `origin/master` with the PAT from `.env` (token redacted from all output), then checks the GitHub API that the remote `master` equals local HEAD — the local tracking ref is not trusted (rename-lock). Rewritten 2026-09-29: the old version used the msys64 path, which has no `git.exe`. |
| `ansys_usermat/apdl/run_chain.ps1` | `-Name n -Runs deck[:case[:minutes]], ... [-WaitFor other.log] [-Export results\<dir>] [-Push] [-DryRun]`: runs partner-element decks one after another through `run_wired.ps1` in a process **detached from the session** (WMI `Win32_Process.Create`), so a Claude session restart does not stop them (5 Oct: a background-shell chain died with the session). Log `F:\biofilm_upf_wired\_chain_<n>.log`; at the end optionally `export_runs_json.py` + commit + push of the JSON. Use it for any ANSYS work longer than a few minutes. |
| `ansys_usermat/apdl/run_apdl.ps1` | `-Deck x.dat [-Ecology] [-Np 4] [-WorkDir ...]`: runs a deck through the custom exe (default `F:\biofilm_upf_kusepy`, `-smp -np 4`; refuses `-Np>1` with the thread-unsafe `F:\biofilm_upf`), starts/stops the material server for ecology decks, clears solver scratch before and after, prints the error/warning counts. |

## Git on this machine — important quirks

- No `git` on PATH by default. **The working git is Portable Git at
  `C:\Users\nishioka\git\cmd\git.exe`** (and `...\git\mingw64\bin\git.exe`) —
  corrected 2026-09-03; an earlier version of this doc pointed at
  `C:\msys64\usr\bin\git.exe`, which only has `git-shell.exe`/
  `git-cvsserver`, no actual `git.exe`. Prepend
  `C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin` to
  `$env:Path` for that PowerShell call.
- `.git\refs\remotes\origin\master` has a persistent rename-lock on this
  machine — `git fetch`/`pull` reliably fails with `error: couldn't set
  'refs/remotes/origin/master'` even though the fetch itself succeeded
  (`FETCH_HEAD` is correct). This only corrupts the *local* tracking ref's
  freshness, not real push/fetch success — `git status`'s "diverged" count
  can be stale/wrong. Verify actual remote state via
  `https://api.github.com/repos/keisuke58/pde-fem-biofilm/commits/master`
  rather than trusting local `git status`, and merge/rebase against
  `FETCH_HEAD` directly when the tracking ref won't update.
- **This repo's working tree has a massive line-ending mismatch** — `git
  status` shows ~480+ files as modified with equal insertions/deletions
  (pure CRLF↔LF churn, zero real content change). **Never `git add -A` or
  `git commit -a`.** Always stage specific files by name.
- `origin` is `https://github.com/keisuke58/pde-fem-biofilm.git`. No
  credential helper was configured as of 2026-08-19 — pushes prompt for
  username/password (PAT) on the console. If a push hangs, it's waiting on
  that prompt.

## Working style for this repo

- **Always reply to the user in Japanese, and keep replies short**
  (decided 2026-10-02). Documents keep their own language (thesis and decks
  in English unless asked otherwise).
- **No PhD** (the user, 6 Oct 2026: never). Plans and strategy aim at the
  Keio master's (final presentation December 2027, graduation March 2028;
  papers should be published before the presentation), papers and a job, not a doctorate or a
  DFG proposal for him.
- **In replies, call the partner "Oliver" (オリバー), not 相方** (asked
  2026-10-06). His element is an earlier version of Felix Klempt's
  implementation (same NEM and AceGen files); Felix's is the later version
  with fixes. Documents keep their own wording.
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
  commit rule exists.)

  **The footer is appended by tooling after the body is written, so it comes
  back on every new PR and a one-off cleanup does not hold.** A note here
  claimed existing PR bodies were cleaned on 2026-09-01; on 2026-09-03,
  nineteen of them still carried one, #31 among them — the PR whose own
  subject is forbidding exactly this. Check the body after opening a PR, and
  run `python scripts/strip_pr_ai_footer.py` (dry run; `--apply` to write) to
  sweep the repository. There is no `gh` CLI and no direct API access from a
  cloud session, so fixing them by hand means retyping each body in full,
  which is how a historical record gets damaged — use the script.

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

## People and how to address them (2026-10-04)

- **Meisam Soleimani is not a professor** (his own clarification, mail of
  6 Oct 2026: he teaches as a lecturer). Write "Dr. Soleimani" (title page:
  Dr.-Ing. Meisam Soleimani); in mail he is fine with "Meisam". This
  replaces the 2026-10-04 rule "Prof. Soleimani"; documents already handed
  over keep their wording.
- **Mayu Muramatsu (Keio) is an associate professor (准教授)**: "Assoc. Prof.
  Muramatsu" / 村松准教授 (村松先生). Address: muramatsu@mech.keio.ac.jp (lab
  contact page). Not part of the December colloquium.
- **Keita Ando (Keio) is an associate professor (准教授)**, Department of
  Mechanical Engineering: "安藤先生". Address kando@mech.keio.ac.jp (from the
  department's 2021 lab brochure; not yet confirmed in use). Contact for the
  Keio 課題研究報告 (the Keio-side thesis presentation). **Mail sent 6 Oct
  2026 from keisuke58@keio.jp** (cc Assoc. Prof. Muramatsu), asking for its
  date and whether it can be online as for earlier double-degree students;
  stated: LUH thesis submitted in November, oral exam in December, return
  to Japan about January 2027. Waiting for his reply. **Settled by the K-LMS
  announcement (7 Oct 2026): no presentation meeting.** Submit the report
  (A4, about 6 pages or more, Japanese or English) and a 5-minute recorded
  PowerPoint talk by **12 Mar 2027 (Fri) 16:00 JST**; in Google Calendar with
  reminders. A report based on a submitted or published single-author paper
  must state the paper's details (e.g. a footnote). Keio student number
  82519093, 開放環境科学専攻. The report itself:
  `luh_summer_2026/1050_Keio/kadaikenkyu2026/kadaikenkyu_nishioka.tex`
  (updated to the thesis content on 6 Oct, branch claude/kadaikenkyu-update).
- Examiners as registered: Prof. Junker (IKM) first, Dr.-Ing. Matthias
  Wangenheim (IDS, wangenheim@ids.uni-hannover.de) second; supervisors Prof.
  Soleimani and Dr.-Ing. Hendrik Geisler.
- Addresses (from earlier mail): junker@ikm.uni-hannover.de,
  soleimani@ikm.uni-hannover.de, geisler@ikm.uni-hannover.de,
  klempt@ikm.uni-hannover.de (Felix Klempt). The colloquium draft in Gmail
  (4 Oct) already has To: Junker, Wangenheim; Cc: Soleimani, Geisler.

## Felix Klempt's code and dissertation chapter (6 Oct 2026): confidential

Felix offered his USERMAT implementation (two species; one species set to
zero is the 2024 paper's model) and the theory of the two-species model,
which is a chapter of his dissertation, still under development. He asked
that it is not shared with anyone without asking him first.

- Never commit his code, his theory notes or excerpts of them to this
  repository (it is public), never upload them to claude.ai, artifacts,
  Drive shares or any other service, and never paste them into mail to
  others. Keep them outside git (IKMHIWI03: `F:\felix_private\`; a cloud
  session: the scratchpad only).
- Notes and results in the repository may say what was compared and what
  agreed or differed (e.g. "the consumption term is g phi c in his
  implementation"), but must not reproduce his code or his derivations.
- Before anything based on his two-species theory goes into the thesis, a
  paper or slides, ask him. If a mistake is found, tell him.

## Slides, notes and other documents for supervisors (decided 2026-10-02)

These apply to every deck, speaker script, email draft or report written for
the supervisors (Dr. Soleimani, Oliver, Assoc. Prof. Muramatsu, the
examiners). Check each one before the document is handed over.

- **Abaqus is not part of this thesis.** Abaqus work is the Keio
  continuation. Do not list Abaqus runs or ANSYS-vs-Abaqus comparisons as
  done work in thesis material. Mentioning Abaqus as *future work at Keio*
  is fine.
- **First person singular.** Write "I", not "we" or "our". Use neutral
  wording such as "this work" or "added in this work" where "I" reads
  badly. In Japanese, write 私, not 私たち.
- **Avoid wording that reads as AI-generated.**
  - Avoid emphatic slogans: "strictly", "exactly as published", "nothing
    else is tuned", "a property of X, not of Y", "The reason is simple",
    "This matters:".
  - Do not use dashes (---) as the main punctuation; use commas, colons or
    a new sentence.
  - Keep bold to a few key numbers.
  - State the result plainly and let the numbers carry it.
- **Klempt et al. 2024 is the reference this work follows.** The paper is in
  the repository: `references/Klempt2024_Hamilton_biofilm_growth_BMMB.pdf`
  (CC BY 4.0; searchable text in the `.txt` next to it). Check the model,
  parameters and notation against it, not against memory or older notes.
- **Notation as in the papers** (decided 2026-10-02), in every document and
  figure:
  - $\phi$ (`\phi`, not `\varphi`) for volume fractions, as printed in
    Klempt 2024, Klempt et al. 2026 and PAMM 2023;
  - $\alpha$ with $\mathbf F_g=\alpha\mathbf I$ and $\alpha(0)=1$ (Klempt 2024);
    the growth is written $\alpha-1$. Do not use $\alpha_K$ or
    $\mathbf F_g=(1+\alpha)\mathbf I$ in documents (the UMAT's internal
    variable is $\alpha-1$; say so where code values are quoted);
  - $k_\alpha$ for the growth rate (Klempt 2024 Eq. 34/36, Table 2);
  - point model (Klempt et al. 2026): $\phi_i$, $\psi_i$,
    $\bar\phi_i=\phi_i\psi_i$, $\phi_0$, $\gamma$, $\eta_i$, $c^*$, $\alpha^*$.
    The papers have no symbol for the share of a species: write
    $\phi_1/(\phi_1+\phi_2)$, not a new symbol such as $\chi_i$.
- **Follow Klempt et al. 2024 (Felix) for the model and its parameters.**
  - Any value not taken from a paper must be marked as such on the slide,
    e.g. the bring-up growth rates or the partner's example-input
    stiffness.
  - Modelling assumptions (e.g. the time link s, φ_cap) are stated as
    assumptions, each with its sensitivity study.
  - **β (diffusion of φ, Eq. 34) = 0.02 mm²/T*** in the ANSYS runs
    (decided 2026-10-05): Klempt 2024 Table 2 (β = 2) converted to the
    partner's 2 mm cube (KLEMPT2024_REPRODUCTION.md §8). The partner's
    example value 1e−4 is shown only as a sensitivity case: its diffusion
    length (~0.01 mm) is below every mesh, so the seed stress does not
    converge. The conversion is an assumption and is stated as such.
- **Units.** The partner's decks are `/units,MPA`, so `YOUNG_BIO = 1000`
  means 1000 MPa. Give units for every material constant. The Klempt 2024
  values are μ = 3.3557 Pa (E = 10 Pa, ν = 0.49).
- **Background appendix.** Every deck gets an appendix with the notation
  (basic variables of the growth field, the point model and the coupling)
  and the background equations, so questions can be answered from the
  slides. The 5 Oct deck (`slides_1005.tex`, appendix A–C) is the template.
- **No internal labels** (stage numbers, `prop(28)` modes, run names) in
  anything the partner or the supervisors see. Describe what a run does.
- **Figures in Times New Roman** (decided 2026-10-02). Every figure script
  calls `figstyle.apply()` from `ansys_usermat/figstyle.py` (Times New Roman,
  Liberation Serif where it is not installed, STIX mathematics). New figure
  scripts use it too; do not set fonts per script.
- **Build and check before handing over.** Keep the deck at 20 pages or
  fewer, with no LaTeX errors and no overfull frames, and look at the
  rendered pages. The build is `build_slides.ps1` on IKMHIWI03, or
  pdflatex/lualatex in a cloud session.

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
