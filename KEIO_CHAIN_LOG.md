# Keio chain log

Record of the decisions made on `keio/fifa-first-run` about what fifa should
run next, and why. Newest at the bottom. Written for a reader who was not
there when the decision was made.

---

## 2026-10-07: first decision, chain `20261008b`

**What I checked.** Both files currently in `abaqus_composition/results_keio_1008/`
against the reference named in their own `# vs ...` / comparison line, after
stripping `#` header lines and blank lines from both sides:

- `fx41_fo_p0.txt` vs `results_1006/fx41_fo_p0.txt`: **identical**, every
  printed digit, including the RMS line (phi 0.066, c 0.025).
- `wp2_imp_free_ph01.txt` vs `results_1007/wp2_imp_free_ph01.txt`:
  **identical**, every printed digit.

Both results' own header already said so; the independent diff confirms it.
No disagreement to chase before starting anything new, so step 5 of
hand-off section 4 ("if a comparison disagreed, find out why before
starting anything new") does not apply this time.

**What I chose.** Chain `20261008b`, manifest `scripts/keio_runs/NEXT.json`,
is exactly `scripts/keio_runs/1008b_fig4_corner.json` (committed by the fifa
session itself in 3c73b8f, described there as "the next chain ... not
started") plus the top-level `"chain"` key. Five runs: `fx41_p100`,
`fx41_p0`, `ff41`, `ff41_fo`, `ff41_g5` — the other five `fig4_corner`
variants already in `results_1006/` (0th/1st-order consumption, with and
without the stiffness penalty, and Felix's front-term variant), each
compared against its own `results_1006/*.txt` reference. 4 CPUs each (the
cap), nice 10, sequential — nothing changed from the manifest fifa proposed.

**Why this and not something else.**
- It is already-reviewed, already-committed work that fifa itself flagged
  as the obvious next step, not something I invented. It extends the one
  check that just passed (`fx41_fo_p0`) to the rest of the same reference
  family, so it directly finishes hand-off §4 step 3 ("the other 5 runs
  with `results_1006/`, to see whether the digit-exact agreement holds
  across the consumption, front-term and phi-cap variants") and KEIO_PLAN
  WP1 ("Klempt 2024 の再現を仕上げる。消費 0次/1次 × φ の上限あり/なし",
  budgeted "1 week, runs automatically").
- It is bounded: 5 jobs (under the ~6-job guideline), 4 CPUs each (at, not
  above, the cap), and based on how long `fx41_fo_p0` took (37 min 55 s,
  with pytest running alongside, so likely an overestimate of the clean
  time) the whole chain should land somewhere around 2.5-3.5 hours of wall
  clock — under the ~4-hour guideline, though I note the uncertainty here
  explicitly rather than treat it as settled.
- I did not reach for WP2 (stress-to-growth) or WP3 (NEM vs Galerkin)
  instead: WP2 is still at the Python-prototype stage and needs a growth
  law decided first (not yet, per KEIO_PLAN §4); WP3's ANSYS side is the
  one with a hard IKMHIWI03 deadline (mid-December), but WP3's Abaqus/Linux
  side is not blocking on fifa capacity this week, and finishing the WP1
  reproduction first is what the hand-off itself asked to do before moving
  on. Running Klempt reproduction variants is also the lowest-risk thing to
  ask of a shared server I have not yet seen fully idle (`uptime`/`nproc`
  are not available to me here; the chain script's own load/disk guards are
  what actually protect the other ~15 users).

**Open questions from PR #58, answered.** The PR is already closed and
merged (into master, 2026-10-07 21:57 JST); two of its four questions were
already resolved by the fifa session's own follow-up work (the scratch
move to `/home` is done and documented in `KEIO_SERVER_ENV.ja.md` §2; the
hand-off §4 step-3 reference choice was settled by actually running
`fx41_fo_p0`, which is the one that matches the paper). The other two are
still open:

1. **Pin numpy/matplotlib/pytest to `requirements.txt`, or relax the pins
   to match fifa?** Recommend **relaxing the pins to match what is on
   fifa** (numpy 2.4.2, matplotlib 3.10.8, pytest 9.0.2), not raising fifa
   to meet `requirements.txt` (numpy 2.4.6, matplotlib 3.11.0, pytest
   9.1.1). The digit-exact agreement just confirmed above was produced with
   numpy 2.4.2; the Abaqus `.dat` output itself doesn't depend on numpy
   (Abaqus compiles Fortran with a fixed, reproducible `ifort` flag set
   regardless of the Python environment — `KEIO_SERVER_ENV.ja.md` §4), but
   the summarise/compare scripts (`summarize_implant.py`,
   `compare_paper.py`, the RMS figures quoted above) do run through numpy,
   and a numpy point release can occasionally change rounding in corner
   cases. There is nothing to gain from upgrading fifa's packages right
   now and a small, silent risk of breaking a comparison that currently
   matches to the last printed digit. If a real reason to upgrade shows up
   later (a CVE, a feature a script needs), re-run the fig4_corner family
   after upgrading and check the agreement still holds before relying on
   it again.
2. **Hand-off §1 names `claude/plan-next-hxjjve` as the working branch, but
   that work merged into master as #56; only one commit of automatically
   exported ANSYS JSON remains on it.** Recommend **correcting it**: the
   statement is simply out of date now and would send a reader of the
   hand-off looking at the wrong branch for anything except that one
   JSON-export commit. I have not made this edit myself this round — the
   task that reached me framed this as a question to answer with a
   recommendation here, not an edit to make, and it touches a document
   outside the result/manifest/log files this decision is scoped to. Flagging
   it here so the correction can be made (by a future decision on this
   branch, or by whoever reads this log) rather than silently dropped.

**`pytest tests/` (473 tests expected per `KEIO_SERVER_ENV.ja.md` §8, this
is a cloud session, not fifa).** One failure:
`tests/test_docs_index.py::test_every_top_level_document_is_reachable_from_the_index`,
which checks that every top-level `*.md` is mentioned in `REPO_MAP.md`,
`DOCS.md` or `README.md`. It named `KEIO_CHAIN_LOG.md` (this file, new) and
`KEIO_SERVER_ENV.ja.md` (added by fifa in 3c73b8f, already committed, never
indexed) as missing. Both are genuine, this decision's own doing in one
case and a small pre-existing gap in the other, not an environment
limitation, so I added one line for each to `REPO_MAP.md` (next to the
existing `KEIO_SERVER_HANDOFF.ja.md`/`KEIO_PLAN.ja.md` entries) and
re-ran `test_docs_index.py` on its own to confirm the fix: 6 passed. I did
not re-run the full 473-test suite a second time (it takes about 14
minutes on this machine); the one failure found was isolated to this test
and the fix touches only `REPO_MAP.md`, so I am not re-running the whole
suite to confirm nothing else regressed from a documentation-only edit.

**Not touched.** `KEIO_SERVER_HANDOFF.ja.md`, `KEIO_PLAN.ja.md` and
`KEIO_SERVER_ENV.ja.md` are read-only this round (see above on the
hand-off). No `.inp`, Abaqus output, or mesh-sidecar `.json` is committed.
Nothing from `F:\felix_private`, Felix's code/theory, or Oliver's ANSYS
sources appears anywhere in what I read or write here.

**Files committed this round.** `scripts/keio_runs/NEXT.json` (the
decision), `KEIO_CHAIN_LOG.md` (this entry), `REPO_MAP.md` (two index
lines, see above).

---

## 2026-10-08: fifa notes, following CLOUD_TO_FIFA.ja.md

**Item 5, the sentence of the hand-off that is out of date.** In
`KEIO_SERVER_HANDOFF.ja.md` section 1, the third bullet:

> - 作業ブランチ：`claude/plan-next-hxjjve`（10月7日の時点。master にはまだ入れていない）。

That work reached master as #56, and later #58, #60 and #61 as well. What is
left on `claude/plan-next-hxjjve` that master does not have is the automatic
ANSYS JSON export from IKMHIWI03, which keeps arriving there. A reader
following the bullet as written would look at the wrong branch for anything
else. The cloud session makes the correction (item 5).

**Item 4, done.** `requirements.txt` now pins numpy 2.4.2, matplotlib 3.10.8
and pytest 9.0.2, the versions on fifa that produced the digit-exact
agreement, with a comment saying what has to be re-run before raising them.

**Item 1, done.** `scripts/run_chain_keio.py` no longer comments on a pull
request and no longer calls `gh`: it appends its report to this file and
commits it with the summaries. Note for the cloud session: `gh` *is* already
authenticated on fifa under `~/.config/gh/hosts.yml` (found while checking how
to reach GitHub from here, before item 1 was written). The chain no longer
uses it, but the token is still on the shared machine; removing it is a
decision for whoever owns that login, not something this session did on its
own.

**Item 3, noted.** `scripts/keio_dispatch.py` was written but never installed
on fifa - the transfer was refused by the sandbox before item 3 arrived, and
item 3 then said not to build it. Nothing of it is committed. What starts a
chain on fifa today is a person running `run_chain_keio.py`; the manifest to
run comes from `scripts/keio_runs/` on master, as item 3 describes.

**A bug the first real chain run found.** `run_comp.sh` copied only the `.inp`
into the work directory, while `compare_ansys.py` reads `<job>.seed.txt` next
to the `.dat`, so every cube comparison on Linux failed after the Abaqus run
had already finished. Fixed in 8200dde, with `scripts/keio_runs/1008s_smoke.json`
(nf8_c6_g1, 25 s, has a reference) as the cheap end-to-end check to run before
trusting the chain with a long one. `run_comp.ps1` copies only the `.inp` too;
on IKMHIWI03 the inputs were generated inside the work directory, so the gap
never showed there.

**Still to merge.** This branch's `run_comp.sh` and master's differ and both
changes are needed: master adds the `ABQ_SCRATCH`/`scratch=`/`TMPDIR`
handling, this branch adds the seed copy. They are complementary. The merge is
deliberately deferred until chain `20261008b` finishes, so that the summarise
scripts cannot change between runs of one chain.

---

## 2026-10-08: chain `20261008b` ran on fifa

The five remaining `fig4_corner` variants, as decided in the entry above. All
five completed and all five match their IKMHIWI03 reference in
`results_1006/` in every printed digit. Comparison re-verified here from the
committed summaries, not copied from their headers.

| job | status | wall clock | vs reference |
|---|---|---|---|
| `fx41_p100` | PASS | 25 min 42 s | identical |
| `fx41_p0` | PASS | 25 min 44 s | identical |
| `ff41` | PASS | 25 min 13 s | identical |
| `ff41_fo` | PASS | 30 min 24 s | identical |
| `ff41_g5` | PASS | 24 min 45 s | identical |

Total 2 h 11 min, sequential, 4 CPUs each, nice 10. Wall clock is
comparable with IKMHIWI03 only if the machine was otherwise idle: the numbers
are deterministic, the timing is not.

**What this settles.** Hand-off section 4 is finished: the user subroutines
compile on Linux (step 1), `wp2_imp_free_ph01` matches (step 2), and the whole
`fig4_corner` family matches (step 3). The agreement holds across every axis
the family varies - consumption g phi (0th order) and g phi c (1st order), the
phi cap on (`--pen 100`) and off (`--pen 0`), the front term as printed and in
Felix's current form (`--felix`), and `--gscale 5`. Together with the
`compile_fortran` flags recorded in `KEIO_SERVER_ENV.ja.md` section 4, a future
disagreement between fifa and IKMHIWI03 should be read as an implementation
difference, not as rounding.

**Note on this entry.** The chain's worker process was started before
`scripts/run_chain_keio.py` was changed to append its own report (ac94d34);
it still had the pull-request path in memory, so it posted one last comment on
PR #58 and did not write here. Chains started after that commit append their
own entry. Nothing else about the run is affected: the summaries, the commit
and the push all came from the chain itself (d10cc81).

<!-- keio-chain -->
```json
{
 "chain": "20261008b",
 "runs": [
  {
   "job": "fx41_p100",
   "cpus": 4,
   "status": "PASS",
   "seconds": 1542,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/fx41_p100.txt",
   "reference": "abaqus_composition/results_1006/fx41_p100.txt"
  },
  {
   "job": "fx41_p0",
   "cpus": 4,
   "status": "PASS",
   "seconds": 1544,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/fx41_p0.txt",
   "reference": "abaqus_composition/results_1006/fx41_p0.txt"
  },
  {
   "job": "ff41",
   "cpus": 4,
   "status": "PASS",
   "seconds": 1513,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/ff41.txt",
   "reference": "abaqus_composition/results_1006/ff41.txt"
  },
  {
   "job": "ff41_fo",
   "cpus": 4,
   "status": "PASS",
   "seconds": 1824,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/ff41_fo.txt",
   "reference": "abaqus_composition/results_1006/ff41_fo.txt"
  },
  {
   "job": "ff41_g5",
   "cpus": 4,
   "status": "PASS",
   "seconds": 1485,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/ff41_g5.txt",
   "reference": "abaqus_composition/results_1006/ff41_g5.txt"
  }
 ]
}
```

---

## 2026-10-08: lab is behind master

`mmc-research-group/nishioka-biofilm-fem` needs a push. fifa cannot reach that
organisation (every attempt is refused by the sandbox), so this is only a note,
per CLOUD_TO_FIFA.ja.md of 8 October, evening.

```
origin/master  ac91c8a
lab/main       0186a2c   (5 commit(s) behind)
```

The push the user has to run, which goes around the stale local master on fifa
and is a fast-forward, so no force:

```
git fetch origin master && git push lab origin/master:main
```

What lab is missing is today's Keio work: the Linux verification of hand-off
section 4 (#58, #60, #66), the chain and its guards, `CLOUD_TO_FIFA.ja.md`
itself (#61, #63), the scratch move (#59), the pins (#62) and the README
update (#64).

---

## 2026-10-08: second decision, chain `20261008c`

**What I checked.** All eight files in `abaqus_composition/results_keio_1008/`,
against the reference named in each one's own header line, stripping `#`
lines and blank lines from both sides:

- `ff41.txt`, `ff41_fo.txt`, `ff41_g5.txt`, `fx41_fo_p0.txt`, `fx41_p0.txt`,
  `fx41_p100.txt`, `nf8_c6_g1.txt`: all **identical** to their
  `results_1006/` reference, independently re-diffed, not just read off the
  `# vs ...` header (these are the five jobs chain `20261008b` ran plus the
  two from the smoke/first-run chains before it).
- `wp2_imp_free_ph01.txt`: no `# vs` header, but its body matches
  `results_1007/wp2_imp_free_ph01.txt` digit for digit, confirming
  `KEIO_SERVER_ENV.ja.md` section 4's claim.

So every comparison fifa has produced so far agrees with its IKMHIWI03
reference in every printed digit. Nothing disagreed, so there is nothing to
chase before starting something new.

**The two open questions from PR #58.**

1. *Package pins.* Already done: `requirements.txt` now reads
   `numpy==2.4.2`, `matplotlib==3.10.8`, `pytest==9.0.2` (commit `ac94d34`),
   i.e. pinned to what fifa actually has, not the newer versions the pins
   originally asked for. This is the right call: the digit-exact agreement
   above was produced with these exact versions running through
   `compare_paper.py`/`compare_ansys.py`/`summarize_implant.py`, all of
   which go through numpy, and there is no numeric reason on the table to
   prefer the newer, unverified versions over an agreement that already
   holds. Recommendation: keep the pins as fifa's versions; only raise them
   again if a specific newer feature is needed, and re-run at least the
   `fig4_corner` family afterwards to confirm the agreement still holds.
2. *Stale branch name in the hand-off.* Already corrected: current
   `KEIO_SERVER_HANDOFF.ja.md` section 1 says the working branch is
   `master` ("10月7日に PR #56 をマージした"), with no mention of
   `claude/plan-next-hxjjve`. For the record, that branch still exists on
   `origin` (`42229b3`, a merge of `keio/fifa-first-run` into it) but is not
   where fifa's instructions point any more. No further edit needed.

**What I chose and why.** `KEIO_PLAN.ja.md` section 3 lists the Klempt 2024
reproduction (WP1) as needing both the 4.1 (`fig4_corner`) and 4.2
(`fig7_high`/`fig7_low`) families across consumption order and the phi cap,
and section 4 gives WP1 a one-week target "実行は自動" (the chain does the
work unattended). Chain `20261008b` finished the 4.1 family; `results_1006/`
already has IKMHIWI03 references for the matching 4.2 family
(`fx42h_p0`, `fx42h_p100`, `fx42h_fo_p0`, `ff42h`, `ff42h_fo`, and the
equivalent `*l*` "low" variants), so verifying them on fifa is the direct
continuation of the same, already-justified check, not a new kind of work.
I queued only the five "high" variants as chain `20261008c`
(`scripts/keio_runs/20261008c_fig7_high.json`, copied into `NEXT.json`): on
IKMHIWI03 these five took 40.5 + 39.2 + 41.3 + 27 + 27.2 = 175.2 minutes, and
the `fig4_corner` five just run on fifa took about as long there as on
IKMHIWI03 (2 h 11 min vs. roughly comparable), so a sixth or more job would
risk pushing past the "about four hours of wall clock" guidance in one go.
The five "low" variants (`fx42l_p0`, `fx42l_p100`, `fx42l_fo_p0`, `ff42l`,
`ff42l_fo`) are the natural next chain after this one finishes; I did not
queue them now so as not to exceed the chain-length guidance.

I considered, and did not queue, three alternatives:

- **Bigger meshes (40^3/48^3, `cs16/24/32/40_c6_g6`, `fr`/`pr`/`nf` mesh
  families).** `KEIO_SERVER_ENV.ja.md` section 9 says these are
  specifically unconfirmed on Linux yet, but each one is a bigger,
  slower job than the ones just proven to work (section 2's disk note: a
  single 20^3 run was 132 MB of scratch; 40^3/48^3 is "nearly 2 GB" each),
  and WP1's own plan puts the 4.1/4.2 family reproduction before the mesh
  study. Finishing the smaller, already-planned family first is lower risk
  per hour of shared-machine time.
- **The implant/tooth stress jobs (`imp_c1`, `tooth_c1`, etc. in
  `results_1006/`, and `results_1007/`'s remaining `wp2_*` variants beyond
  `wp2_imp_free_ph01`).** These belong to WP2 (the application paper
  candidate), which `KEIO_PLAN.ja.md` section 0 says is not decided as the
  main target until December, after the ANSYS WP3 results are in. Nothing
  about them is wrong to run, but WP1 (the reproduction this fifa loop
  exists to finish first, per the hand-off) is the one with an open,
  dated task.
- **Starting WP2 application work (new growth-law runs) or WP4/5/6.**
  These need a growth-law and alternating-scheme decision
  (`KEIO_PLAN.ja.md` section 4, WP2) that has not been made yet; queuing
  compute for them now would be inventing work ahead of a decision that is
  not mine to make.

**Validation before writing `NEXT.json`.** Ran
`python3 scripts/run_chain_keio.py --name 20261008c --runs
scripts/keio_runs/20261008c_fig7_high.json --dry-run --worker`: all five
`make_klempt_inp.py` calls succeeded (rc 0, each `.inp` generated) before
committing the manifest. The dry run also appended a throwaway entry to this
file (`seconds: 0`, `verdict: skipped`); that entry was reverted and is not
part of this one. `python3 -m pytest tests/ -q` was also started as the
general Python-side check this session is allowed to do, but did not finish
within this session's available time (it was still running when this entry
was committed); that is a cloud-session timing limit, not a test failure, and
does not block this decision since the dry run above is the actual
validation `run_chain_keio.py`'s own docstring asks for before trusting a
manifest.

---

## 2026-10-08 23:32: chain `20261008c` ran on fifa

Automatic entry from `scripts/run_chain_keio.py`. No analysis: what runs next is for the cloud session to decide.

| job | status | wall clock | vs reference |
|---|---|---|---|
| `fx42h_p0` | PASS | 52 min 19 s | identical |
| `fx42h_p100` | PASS | 50 min 18 s | identical |
| `fx42h_fo_p0` | PASS | 52 min 38 s | identical |
| `ff42h` | PASS | 33 min 36 s | identical |
| `ff42h_fo` | PASS | 33 min 48 s | identical |



Wall clock is comparable with IKMHIWI03 only if the machine was otherwise idle: the numbers are deterministic, the timing is not.

<!-- keio-chain -->
```json
{
 "chain": "20261008c",
 "runs": [
  {
   "job": "fx42h_p0",
   "cpus": 4,
   "status": "PASS",
   "seconds": 3139,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/fx42h_p0.txt"
  },
  {
   "job": "fx42h_p100",
   "cpus": 4,
   "status": "PASS",
   "seconds": 3018,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/fx42h_p100.txt"
  },
  {
   "job": "fx42h_fo_p0",
   "cpus": 4,
   "status": "PASS",
   "seconds": 3158,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/fx42h_fo_p0.txt"
  },
  {
   "job": "ff42h",
   "cpus": 4,
   "status": "PASS",
   "seconds": 2016,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/ff42h.txt"
  },
  {
   "job": "ff42h_fo",
   "cpus": 4,
   "status": "PASS",
   "seconds": 2028,
   "verdict": "identical",
   "detail": "15 lines, every digit equal",
   "file": "abaqus_composition/results_keio_1008/ff42h_fo.txt"
  }
 ]
}
```

## 2026-10-08: third decision, chain `20261008d`

**What I checked.** All thirteen files currently in
`abaqus_composition/results_keio_1008/` (the eight from before plus the five
`fig7_high` jobs chain `20261008c` just ran: `fx42h_p0`, `fx42h_p100`,
`fx42h_fo_p0`, `ff42h`, `ff42h_fo`), independently re-diffed against the
reference named in each file's own header, stripping `#` lines and blank
lines from both sides rather than trusting the `# vs ...: identical` line
printed in the file. All thirteen agree with their `results_1006`/`results_1007`
reference in every printed digit; `wp2_imp_free_ph01.txt` again has no
`# vs` header but matches `results_1007/wp2_imp_free_ph01.txt` body for body.
Nothing disagreed, so there is nothing to chase before starting something new.

**The two open questions from PR #58.** Both were already answered and acted
on in the previous decision (2026-10-08, chain `20261008c`): `requirements.txt`
still pins `numpy==2.4.2`, `matplotlib==3.10.8`, `pytest==9.0.2` (fifa's actual
versions, not the original newer pins), and `KEIO_SERVER_HANDOFF.ja.md`
section 1 still names `master` as the working branch with no mention of
`claude/plan-next-hxjjve`. Re-checked on disk, nothing has regressed; no
further edit needed.

**What I chose and why.** `scripts/keio_runs/20261008d_fig7_low.json` was
already sitting prepared (commit `76e8ed1`, from this session's own earlier
turn): the five `fig7_low` variants (`fx42l_p0`, `fx42l_p100`, `fx42l_fo_p0`,
`ff42l`, `ff42l_fo`) that close out `KEIO_PLAN.ja.md` section 4's WP1 task
("消費 0次/1次 × φ の上限あり/なし × 4.1/4.2"). `fig4_corner` (4.1) finished in
chain `20261008b`, `fig7_high` (4.2 high) just finished in chain `20261008c`;
`fig7_low` (4.2 low) is the only remaining piece of that family, and
`results_1006/` already holds IKMHIWI03 references for all five variants, so
this is a direct continuation of the same already-justified check, not new
work. I re-ran `scripts/check_manifest.py`'s comparison by hand (the script's
own `ROOT` is hardcoded to fifa's path, so I reproduced its logic against
this checkout instead of invoking it directly): all five `make_klempt_inp.py`
command lines in the manifest match the flags recorded in their
`results_1006` reference headers word for word, and each `summarise --case`
matches too. `python3 -c "import json;json.load(...)"` parses the manifest
cleanly. I also ran the klempt/compare_paper/keio-related slice of
`tests/` (`pytest tests/ -k "klempt or compare_paper or keio" -q`, bounded to
90 s rather than the full ~14 min suite): 26 passed, 0 failed.

By the IKMHIWI03 times (23 + 25.7 + 34.1 + 22.5 + 23.9 = 129.2 min) this chain
is shorter than the `fig7_high` chain that just ran on fifa (175.2 min on
IKMHIWI03, 222.6 min actual on fifa); even at the same roughly 1.3x slowdown
fifa showed for `fig7_high`, `fig7_low` should land under 3 hours, well inside
the four-hour/six-job guidance, and it is the last chain needed to close WP1's
reproduction family — after this, no planned-but-unqueued WP1 work remains and
the next decision has to look at WP2/WP3 instead.

I considered, and did not queue, the same three alternatives as last time,
for the same reasons: bigger meshes (`KEIO_SERVER_ENV.ja.md` section 9 still
marks 40^3/48^3 as unconfirmed on Linux and they are the slower, riskier next
step); the implant/tooth stress jobs (WP2, not decided as the main paper
target until December per `KEIO_PLAN.ja.md` section 0); and new WP2/WP4/5/6
compute, which needs a growth-law or scope decision that is not mine to make.
Promoted `scripts/keio_runs/20261008d_fig7_low.json` to `NEXT.json` unchanged
(same chain name, already new: `git log --all` shows no commit with
"chain 20261008d" or any chain-log entry for it).

## 16 Oct 2026: tooth data provenance and licence (CLOUD_TO_FIFA, 15 Oct late night, item 1), UEL on the implant sector (item 2)

**Tooth data on fifa** (not in the repository; nothing was copied into it):
- `/home/nishioka/IKM_Hiwi/FEM/external_tooth_models/OpenJaw_Dataset/` (92 MB): Patient_1 to 3, per-tooth STL files
  (`Patient_1/Teeth/P1_Tooth_18.stl` and so on), `Dataset_Structure.pdf`. It has no LICENSE or README of its own.
- `/home/nishioka/IKM_Hiwi/Open-Full-Jaw` (clone of github.com/diku-dk/Open-Full-Jaw, 95 MB) and a copy of its code in
  `FEM/external_tooth_models/Open-Full-Jaw-main`. Its LICENSE is **CC BY-NC-SA 4.0**.
- Tooth models and inp files for the 3-tooth assembly: `/home/nishioka/IKM_Hiwi/FEM/*.inp` (38 files, e.g.
  `biofilm_3tooth.inp` 26 MB, `p23_*`, `p30_*`, `p31_*`), outside the repository.
- Source: Gholamalizadeh et al., Open-Full-Jaw, Comput. Methods Programs Biomed. 224 (2022) 107009. Web search found
  the OpenJaw archive at Univ. Copenhagen (erda.ku.dk) described as "research purposes only, cite the work".
- `tier2b_real/_crown_real_watertight.stl` (already in the public repository) is derived from `P1_Tooth_30.stl`
  (`tier2b_real/mesh_crown_real.py`: clinical crown cut out, cervical loop closed).

**Licence reading (not legal advice, the licence text of the dataset itself was not found):** the repository code is
CC BY-NC-SA 4.0 (attribution, non-commercial, share-alike); the dataset archive says research use and citation.
Neither text was found to allow redistribution of the STL files. Suggestion for the cloud session: keep the tooth data
out of the repository; for `_crown_real_watertight.stl` either remove it from the public repository or state the
source, the cut and CC BY-NC-SA 4.0 in `THIRD_PARTY.md` after the owner has decided; ask the dataset authors
if redistribution is wanted.

**UEL on the implant sector** (`p3/uel/make_implant_uel_inp.py`, branch `keio/uel-implant`): Abaqus and the Python
emulation agree at 2x6x8 to 4.1e-7 (alpha - 1) and 4.7e-7 (vM) after a fix of the Gauss-point record in
`unified_uel.f90` (alpha - 1 had been one increment behind; commit 55062c1). With `--bulge 1`: 3.9e-7 and 5.2e-7.
Seed means (inner layer, 48 elements): alpha - 1 0.04634 (bulge 0), 0.04427 (bulge 1); vM/E 0.04775 and 0.04620;
largest vM/E 0.2728 and 0.2927, both in the top row of elements at theta 67.5 deg (r 2.09 and 3.08).

## 16 Oct 2026, fifa: tooth data provenance and licence, UEL on the implant sector

Request: `CLOUD_TO_FIFA.ja.md`, entry of 15 Oct late night, item 1 (no computation). Nothing from the dataset was copied into the repository.

**Where the tooth data is on fifa** (`/home/nishioka/IKM_Hiwi/`, not in this repository)
- `FEM/external_tooth_models/OpenJaw_Dataset/`: Patient_1, Patient_2, Patient_3, `Dataset_Structure.pdf`; per patient a `Teeth/` folder of STL files (`P1_Tooth_18.stl` ...); 92 MB in all. It contains no LICENSE or README of its own.
- `FEM/external_tooth_models/Open-Full-Jaw-main/` and `Open-Full-Jaw/` (clone of `https://github.com/diku-dk/Open-Full-Jaw.git`, 95 MB): the pipeline and docs of the same group; the dataset folder holds only a README.
- `FEM/*.inp`: 38 decks built on these STL files (`biofilm_3tooth.inp` 26 MB, `p23_*`, `p30_*`, `p31_*`).
- `tier2b_real/mesh_crown_real.py` reads `.../OpenJaw_Dataset/Patient_1/Teeth/P1_Tooth_30.stl`; `tier2b_real/_crown_real_watertight.stl` (194 KB, already in this public repository) is the clinical crown of that tooth, uprighted, scaled and capped.

**Source**: Gholamalizadeh et al., Open-Full-Jaw: An open-access dataset and pipeline for finite element models of human jaw, Comput. Methods Programs Biomed. 224 (2022) 107009, University of Copenhagen (cited in `refs_openjaw.bib`, `biofilm_3tooth_refs.bib`). The OpenJaw archive (ERDA, Univ. Copenhagen) states research use and citation; I did not find a redistribution grant.

**Licence**: `Open-Full-Jaw/LICENSE` is CC BY-NC-SA 4.0 (attribution, non-commercial, share-alike). It is the licence of that repository (code, docs); whether it covers the OpenJaw STL files is not stated in the data folder. Unknown: the exact terms of the three-patient STL archive.

**Consequences / proposal for the cloud session**
- Do not put any OpenJaw STL or mesh derived from it into the public repository until the terms are confirmed (ask the authors or read the ERDA archive page).
- `tier2b_real/_crown_real_watertight.stl` is a derivative of P1_Tooth_30. Under CC BY-NC-SA it would need attribution, NC and the same licence; under "research use only" it should probably be removed. Decide after the terms are known.
- Add a provenance entry to `THIRD_PARTY.md` (citation above, licence as found, "STL files not redistributed").

## 16 Oct 2026, fifa: UEL checks (branches keio/uel-check, keio/uel-implant)

- Frozen seed 24^3 (`p3/uel/results/seed_frozen_24.json`): vM mean over elements within 0.2 (3D) of (0.5, 0.5, 0) is 6.549e-3, ANSYS 24^3 6.536e-3 (+0.2 %). Seed (phi > 0.5, 420 elements): alpha - 1 8.70e-4, vM 5.093e-3, largest vM 1.345e-2 at (0.479, 0.646, 0.021).
- Implant sector (`p3/uel/make_implant_uel_inp.py`, 2x6x8 elements, 90 degrees): Abaqus and the Python emulation agree to 4e-7 (alpha - 1) and 5e-7 (vM/E), also with `--bulge 1` (3.9e-7, 5.2e-7). Inner layer: alpha - 1 0.04634, vM/E 0.04775, largest vM/E 0.2728 at (r, theta, z) = (2.09, 67.5, 1.875). With `--bulge 1`: 0.04427, 0.04620, 0.2927 at (3.08, 67.5, 1.875).
- Bug found and fixed in `unified_uel.f90`: the Gauss-point record (`uel_gp.txt`) was taken before `store_growth`, so its alpha - 1 was one increment behind (vM and pressure were not). Now taken after. The 8^3 seed in Abaqus still gives vM mean 0.01342, max 0.01577 (emulation 0.013419 / 0.015773).

## 16 Oct 2026, fifa: Abaqus column and seed 8^3 (entry of 15 Oct, `KEIO_SERVER_HANDOFF.ja.md` section 7)

Built with the corrected `unified_uel.f90` (Gauss-point record after the growth update). `abaqus job=... user=p3/uel/unified_uel_abaqus.for cpus=1 interactive`; linked and ran without error, no `/libs:dll` needed on Linux.
- `uel_column` (20 elements): x = 0, NT11 0.77532, NT12 0.93565 (IKMHIWI03 0.7753, 0.9357).
- `uel_seed8` (8^3, default options): elements with phi > 0.5: vM mean 0.01300, max 0.01637 (expected 0.01300 / 0.01637). alpha - 1 mean 0.03571.
Output: `p3/uel/results/abaqus_fifa/` (`.sta` and the Gauss-point files, gzip).

## 16 Oct 2026, fifa (overnight): implant sector mesh convergence, Abaqus checks, UEL seed 32^3
```
bulge 0
mesh      elems  seed a-1   seed vM/E  int a-1    int vM/E   max vM/E  (r, theta, z)            layer row_from_top  max w/o edge rows
2x6x8        96  0.04634   0.04775   0.01795   0.03544   0.2728   (2.094, 67.5, 1.875)    0     0   0.1285
4x12x16     768  0.03523   0.02358   0.02176   0.01825   0.2375   (2.077, 71.2, 1.938)    0     0   0.1608
6x18x24    2592  0.02962   0.01482   0.02137   0.01199   0.2004   (2.069, 77.5, 1.958)    0     0   0.1586
bulge 1
mesh      elems  seed a-1   seed vM/E  int a-1    int vM/E   max vM/E  (r, theta, z)            layer row_from_top  max w/o edge rows
2x6x8        96  0.04427   0.04620   0.01490   0.03253   0.2927   (3.076, 67.5, 1.875)    0     0   0.1262
4x12x16     768  0.03302   0.02383   0.01901   0.01782   0.2912   (3.072, 71.2, 1.938)    0     0   0.1809
6x18x24    2592  0.02774   0.01599   0.01908   0.01262   0.2641   (3.067, 72.5, 1.958)    0     0   0.1956
```
Abaqus against the emulation (largest relative difference per element, alpha - 1 / vM), last increment:
- 4x12x16 (768 elements): 1.8e-07 / 2.9e-07
- 4x12x16_b1 (768 elements): 1.5e-07 / 2.7e-07
UEL frozen seed 32^3: seed alpha - 1 8.8340e-04, seed vM 5.0865e-03, largest vM 1.4426e-02; vM mean within 0.2 of (0.5, 0.5, 0): 6.5117e-03 (544 elements); ANSYS 32^3: 6.5e-3

---

## 2026-10-28: UEL implant sector after the shape-gradient fix (CLOUD_TO_FIFA, 17 Oct)

**Done.** Merged `origin/master` (to 27 Oct) into `keio/uel-implant`. Five files conflicted
(`p3/uel/implant_convergence.py`, `results/implant_{4x12x16,4x12x16_b1,6x18x24,6x18x24_b1}.json`);
I took the master versions, since the fifa results of 16 Oct were computed with the wrong shape
gradient. Rebuilt `libunified_uel.so` and `unified_uel_abaqus.for`.

**Abaqus 2x6x8 (bulge 0, free outer face): PASS.** Seed alpha - 1 0.04412878 in both, seed vM/E
0.02501274 / 0.02501273, max vM/E 0.1061752 in both (element 85); largest element-wise difference
4.1e-8 (alpha - 1) and 3.2e-8 (vM/E). Record: `results/abaqus_fifa/uel_implant2x6x8_gp.txt.gz`.

**Abaqus 6x18x24, fixed seed (`--seed-layers 3`), bulge 1: FAIL, all NaN.** The job reports
"completed successfully", but every Gauss-point record (7 blocks, increments 10..60) has NaN in
phi, c, von Mises and pressure, and alpha - 1 = 0. `compare_implant_gp.py` gives
"Abaqus nan" for the seed vM/E and the maximum.
The `.sta` shows the first growth increment (dt = T/60 = 3.958e-3) accepted after 4 iterations, with
no cutback; then 59 further increments of 1 iteration each, all at dt = 3.958e-3. The expected
sequence was a cutback (attempt 2, dt 9.9e-4) and a return to 3.96e-3 by factors of 1.5. The `.msg`
has 20357 warnings, all "system matrix has N negative eigenvalues" (5 to 73), which the
coupled temperature-displacement step also gives in the 2x6x8 job. Not investigated further.
Likely cause (not tested): the near-singular first-increment tangent that the emulation cuts
back is not caught in Abaqus, and NaN enters the state variables at that increment.
Files: `results/abaqus_fifa/uel_implant6x18x24_b1_sl3.sta`, `..._gp.txt.gz`.
Deck: `make_implant_uel_inp.py abaqus_in --nr 6 --nt 18 --nz 24 --seed-layers 3 --bulge 1`.

**Running, not finished.** 8x24x32 fixed seed, bulge 0 and 1 (Python emulation, nohup, started
17:44 JST). At increment 3 both show residual halving in the line search. If one stops with
"5 cutbacks did not converge", the increment number and the last 30 lines of stderr go here.

**Not done.** Homeostatic-pressure Abaqus check (2x6x8, p_h/E = 0.05) and the Lean/Mathlib check.

**Finished (29 Oct).** 8x24x32 fixed seed (`--seed-layers 4`, Python emulation), bulge 0 and 1:
seed alpha - 1 0.04127 / 0.03979, seed vM/E 0.02904 / 0.02608, interior vM/E 0.02553 / 0.02254,
max vM/E 0.3840 / 0.3645. The seed has 4 layers here and 3 in the 6x18x24 runs, so the two
meshes are not compared directly; a 6x18x24 run with 4 layers is needed for the convergence check.
Files: `p3/uel/results/implant_8x24x32_sl4.json`, `implant_8x24x32_b1_sl4.json`.

## 2026-10-29: Abaqus 6x18x24 bulge 1, NaN found and guarded

The NaN of 28 Oct came from the residual, not from the state: at the first attempt of increment 1
(dt = 3.958e-3) 72 UEL elements returned a non-finite R (outer two radial layers, top three z rows,
all theta; the seed layer at the top edge). Abaqus read the NaN residual as converged, so there was no
cutback. `unified_uel.f90` now checks sv after `advance_composition` and R after `residual`; a non-finite
value gives PNEWDT = 0.25 and a line "UEL non-finite ..." on stdout (element, time). The job then
cut back (attempt 2, dt 9.896e-4) and ran 64 increments to the end, no NaN.
Against the emulation (`compare_implant_gp.py`): seed alpha - 1 0.04467319 / 0.04469628, max vM/E
0.3017095 / 0.3018817 (same element 2485); largest element-wise difference 1.9e-4 (alpha - 1), 1.7e-4
(vM/E), tolerance 1e-5, so FAIL by the tolerance. The step sequences differ slightly (Abaqus three
increments of 9.896e-4, emulation two), which is a time-discretisation difference, not shown to be the
whole cause. Not done: why the Newton iterate of the first attempt becomes non-finite; a comparison
with the same increment sequence (fixed dt = 9.896e-4 in both).
Files: `results/abaqus_fifa/uel_implant6x18x24_b1_sl3.sta`, `..._gp.txt.gz`.

## 2026-10-29 (morning): CLOUD_TO_FIFA 28 Oct, steps 1-3 (6x18x24, bulge 1 NaN isolation)

Merged `origin/master` (to `c839bcc`) into `keio/uel-implant`, no conflicts. The UEL source in the repository is
unchanged. Steps 1-3 were run with the UEL of `9776d59`, which is **before** the guard of `3ab9c00`, because with the
guard the NaN no longer appears. The only change in that copy is the Gauss-point output, written at KINC <= 5 as well as
every 10th increment. The source and decks are on fifa only, in `~/abaqus_work/nan_probe_20261029/`
(`probe_unguarded.for`, `run_chain2.sh`). Each job was `abaqus job=... user=probe_unguarded.for cpus=1 interactive` with
`BIOFILM_ECO_CASE` and `BIOFILM_UEL_GP` set.

Two problems with the procedure, found on the way (not fixed, for the cloud):
- **`make_implant_uel_inp.py --ninc n` does not reach the deck.** It changes only the `--run` emulation. The GROW data
  line (line 151) always uses the global `NINC = 60`, so the deck of step 2 was the same as step 1's. For step 2 the data
  line was set by hand in the copy of the deck: `0.0009895833333, 0.2375, 9.6639e-07, 0.0009895833333`, which is
  dt0 = DTMAX = T/240 and DTMIN = T/240/4^5, as the generator writes them for T/60.
- **Without `BIOFILM_UEL_GP`, `uel_gp.txt` is written to the solver's scratch directory and deleted with it.** This
  happens on fifa because `abaqus_v6.env` sets `scratch`. The first runs of steps 1 and 2 lost their records this way
  and were run again.

**1. Where the NaN enters: GROW increment 1, attempt 1, equilibrium iteration 4.** dt = 3.958e-3. The `.sta` line is
`2 1 1 0 4 4 0.00396 0.00396 0.003958`: 4 iterations, no cutback, then 59 increments of 1 iteration each. The `.msg`
of increment 1:
- iteration 1: largest scaled residual force 1.27e-2, correction to displacement -4.1e-2 (node 3197, DOF 3);
  73 negative eigenvalues in the next matrix
- iteration 2: residual 1.27e-3 (node 3181, DOF 2), correction 1.96e-2
- iteration 3: residual 1.13e-2 (node 3098, DOF 3); **correction 0.297 (node 3233, DOF 3), increment of displacement
  0.302 (node 3198, DOF 1)**
- iteration 4: "AVERAGE HEAT FLUX NaN", "ALL FORCE RESIDUALS ARE ZERO", "ALL HEAT FLUX RESIDUALS ARE ZERO". Abaqus
  takes this as converged.

The first record after GROW starts (end of increment 1) already has every element non-finite: all 2592 elements, all
8 Gauss points. phi, c, von Mises and pressure are NaN, and alpha - 1 = 0. So the record alone cannot narrow it down to
an element. The 29 Oct guarded run localised it to the residual of 72 elements: the outer two radial layers and the top
three z rows. Values just before, in the INIT record (step 1, increment 1), over all elements: phi 1e-4..0.9001,
c 0.1..0.8098, alpha - 1 = 0 everywhere.
Files: `p3/uel/results/abaqus_fifa/uel_implant6x18x24_b1_sl3_unguarded.sta`, `..._unguarded_gp.txt.gz`.

**2. dt = T/240 (9.896e-4, data line set by hand, unguarded UEL): runs, no NaN.** 240 increments, 1 attempt each, no
cutback. Against the emulation `results/implant_6x18x24_b1_sl3.json` (`compare_implant_gp.py`):

| | Abaqus T/240 | emulation |
| --- | --- | --- |
| seed alpha - 1 | 0.0443474 | 0.04469628 |
| seed vM/E | 0.04267654 | 0.04313718 |
| max vM/E | 0.3002227 (element 2485) | 0.3018817 (element 2485) |

The largest element-wise difference is 1.80e-3 (alpha - 1) and 1.94e-3 (vM/E): FAIL at tol 1e-5. The two do not use
the same steps. The emulation cuts back once to 9.9e-4 and then grows dt by 1.5 back to T/60; Abaqus stays at T/240.
The difference is therefore a time-step difference at least in part. An emulation `--run --ninc 240` was not run.
Files: `p3/uel/results/abaqus_fifa/uel_implant6x18x24_b1_sl3_n240.sta`, `..._n240_gp.txt.gz`.

**3. bulge 0, same 6x18x24, `--seed-layers 3`, no `--ninc` (T/60, unguarded UEL): PASS.** 60 increments, 1 attempt
each, no NaN. Against `results/implant_6x18x24_sl3.json`: seed alpha - 1 0.04787072 in both, seed vM/E 0.0490564 in
both, max vM/E 0.3236237 in both (element 2485). The largest element-wise difference is 5.08e-8 (alpha - 1) and
5.04e-8 (vM/E).
So with dt = T/60 the 6x18x24 mesh passes with bulge 0 and gives NaN with bulge 1: **the bulge matters, not the mesh**.
With bulge 1 a smaller first step (T/240) avoids the NaN.
Files: `p3/uel/results/abaqus_fifa/uel_implant6x18x24_sl3.sta`, `..._sl3_gp.txt.gz`.

**4.** The 8x24x32 emulation was already finished and recorded above (29 Oct). **5.** (homeostatic-pressure Abaqus,
Lean/Mathlib) not done.

## 2026-10-29: CLOUD_TO_FIFA 28 Oct, step 5 (homeostatic-pressure Abaqus check, 16 Oct): PASS, 4 of 4

The UEL is the one in the repository at `1f34cfb`, with the guard. `unified_uel_abaqus.for` was regenerated with
`make_abaqus_source.py`; the file is identical to the committed one. The decks are 2x6x8, bulge 0, p_h/E = 0.05:
`make_implant_uel_inp.py abaqus_in --p-ref 0.05 [--gate-implicit] [--outer fixed]`. The jobs were run as in the
16 Oct instructions, with `BIOFILM_UEL_GP` set (script `p3/uel/abaqus_in/run_ph.sh`). Each job took about 25 s.

All four jobs ran 60 increments, 1 attempt each, no cutback, and wrote no "UEL non-finite" line.
`compare_implant_gp.py --tol 1e-6`:

| job | seed alpha - 1 (Abaqus / Python) | seed vM/E | max vM/E (element 85) | largest difference alpha - 1 / vM/E | iterations per increment |
| --- | --- | --- | --- | --- | --- |
| `_ph0.05` (explicit, outer free) | 0.0277612 / 0.02776119 | 0.0151908 / 0.0151908 | 0.05699419 / 0.05699419 | 2.7e-8 / 4.9e-9 | 2: 51, 3: 8, 4: 1 |
| `_ph0.05_gi` (implicit, outer free) | 0.02758512 / 0.02758512 | 0.01509242 / 0.01509242 | 0.05657986 / 0.05657986 | 2.0e-8 / 2.5e-9 | 2: 53, 3: 6, 4: 1 |
| `_ph0.05_of` (explicit, outer fixed) | 0.008987526 / 0.008987525 | 0.00974648 / 0.00974648 | 0.04093257 / 0.04093258 | 2.3e-9 / 5.6e-9 | 2: 52, 3: 7, 4: 1 |
| `_ph0.05_gi_of` (implicit, outer fixed) | 0.00888262 / 0.00888262 | 0.009664487 / 0.009664486 | 0.04059057 / 0.04059057 | 3.6e-9 / 4.1e-9 | 2: 57, 3: 3 |

Files: `p3/uel/results/abaqus_fifa/uel_implant2x6x8_ph0.05{,_gi,_of,_gi_of}.sta` and `..._gp.txt.gz`.
Step 4 of the 16 Oct instructions (6x18x24 emulation with the gate, about 35 min each) was not run.

**Lean / Mathlib (16 Oct night): skipped.** Installing elan with `curl ... elan-init.sh | sh`, as step 0 says, was
refused by the permission settings of the Claude session on fifa. The user decided (29 Oct) to skip the Lean / Mathlib
check on fifa: nothing was installed and `keio/lean-check` was not created.

## 2026-10-29: CLOUD_TO_FIFA 29 Oct, step 1 (same dt = T/240 in both): PASS

Merged `origin/master` (with the `--ninc` fix of `make_implant_uel_inp.py`); the deck now has the GROW data line
`0.0009895833333, 0.2375, 9.6639e-07, 0.0009895833333` by itself. The UEL is the guarded one in the repository.
Deck and emulation: `make_implant_uel_inp.py [abaqus_in | --run] --nr 6 --nt 18 --nz 24 --seed-layers 3 --bulge 1 --ninc 240`
(script `p3/uel/abaqus_in/run_n240.sh`).

- Abaqus: 240 increments, 1 attempt each, no cutback, no "UEL non-finite" line, 12 min.
- Emulation: 240 increments, no cutback, 3 Newton iterations each (720 in all), 1 h 42 min (6147 s).
- `compare_implant_gp.py`: seed alpha - 1 0.0443474 in both, seed vM/E 0.04267654 in both, max vM/E 0.3002227 in
  both (element 2485); largest element-wise difference 4.90e-8 (alpha - 1), 4.24e-8 (vM/E). **PASS** (tol 1e-5).

So the FAIL of 28/29 Oct (1.8e-3) came only from the different increment sequences. Step 2 is not needed.

Step 3 (recorded only, no conclusion): det J at the 8 Gauss points from the deck geometry. The elements are hardly
distorted: the min/max ratio of det J per element is at least 0.971 with bulge 1 (median 0.983) and 0.9885 with bulge 0.
The elements around the nodes of the jump in iteration 3 (3233, 3198, 3197; top outer edge, z = 2) have 0.992, so they
are not the worst. Bulge 1 widens the outer radius at the top edge from 2.26 to 3.26: det J ranges 7.8e-5..1.24e-4 with
bulge 1 and 7.8e-5..8.7e-5 with bulge 0, so the top outer elements are up to about 1.4 times larger.

Files: `p3/uel/results/implant_6x18x24_b1_sl3_n240.json`, `p3/uel/results/abaqus_fifa/uel_implant6x18x24_b1_sl3_n240.sta`
and `..._gp.txt.gz` (this replaces the unguarded T/240 record of the morning, same name).

## 2026-10-29 (evening, in progress): growth-living (B) and growth-living-energy (A), 4 cases x 32/64/128

CLOUD_TO_FIFA 29 Oct evening. `pyamg` 5.3.0 was installed into the venv on fifa. `ecology_native_ctypes` built, so the
runs keep `--fortran`. The options are the common P of the instructions (`--dtc 0.0025`). The runs go one at a time
(nohup, `~/runs/growliv_20261029/run.sh`), in the order 32 B, 32 A, 64 B, 64 A, 128 B, 128 A.
Interim push (29 Oct, 23:47 JST): the 32^3 runs (A and B) and the 64^3 B runs are done, about 40 s each at 32^3. The JSONs
have `_dtc0.0025_` in the name; the 32^3 `_growliv` files of master have no dtc tag (outer step 0.025), so they are a
different setting. The table, ratios and GCI follow when 128^3 is done.

## 2026-10-30: growth-living B/A, 4 cases x 32/64/128, limit and dc/dt (CLOUD_TO_FIFA 29 Oct evening, night, late night)

All 26 dc/dt runs (`--transient-c 1e8`) and all limit runs finished with exit 0. The options are the common P of
29 Oct evening (`--dtc 0.0025`, `--fortran` built). Both series ran one at a time. The dc/dt chain ran beside the limit chain,
because the runs use 1 core and about 1.1 GB.
Wall time per run (limit / dc/dt): 32^3 38-47 s / 32-40 s; 64^3 269-355 s / 225-288 s; 128^3 54-66 min / 45-55 min.
Scripts: `~/runs/growliv_20261029/run.sh`, `run_base64.sh`, `~/runs/transc_20261030/run.sh` (fifa only); table script
`~/runs/analysis/growliv_table.py`.

The 64^3 dc/dt baseline `mixed_ab` (`..._nf64_..._mixeq_..._transc1e+08_..._mixed_ab.json`) is also on master (cloud).
The fifa run agrees with it to round-off: seed alpha - 1 relative difference 7e-15. The master file was kept and the
fifa copy set aside.

**64^3 limit baselines (29 Oct night) against `tab:ch5-gci`:** separate_ab alpha - 1 8.9310e-4, vM 5.9998e-4 (table
8.931 / 6.000); mixed_ab 8.0402e-4 / 2.8151e-4 (table 8.04 / 2.82). Both agree within 0.1 % of the table's last digit.

Quantities are seed means at the end: sum phibar_i (living), share A, c, alpha - 1, von Mises.
"/ none" is the ratio to the no-growth-law run of the same grid and the same dc/dt setting. It is left out ("-") where that
run is not in `results/2026-10-3d/`: limit 128^3 separate_ab and mixed_ab (they are in the thesis table and A-89, but not as
JSON), and dc/dt separate / mixed at any grid.

### limit (no dc/dt)
| law | case | grid | sum phibar | share A | c | alpha-1 (1e-4) | vM/Pa (1e-4) | alpha-1 / none | vM / none |
|---|---|---|---|---|---|---|---|---|---|
| none | separate | 32^3 | 0.9645 | 0.9779 | 0.5293 | 8.7289 | 7.4404 | - | - |
| none | separate | 64^3 | 0.9646 | 0.9780 | 0.5358 | 9.0140 | 6.7297 | - | - |
| none | separate | 128^3 | 0.9647 | 0.9780 | 0.5394 | 9.1105 | 6.5853 | - | - |
| none | separate_ab | 32^3 | 0.1556 | 0.9517 | 0.9149 | 8.4529 | 5.5229 | - | - |
| none | separate_ab | 64^3 | 0.1562 | 0.9517 | 0.9133 | 8.9310 | 5.9998 | - | - |
| none | mixed | 32^3 | 0.9622 | 0.8776 | 0.6424 | 7.4612 | 2.1290 | - | - |
| none | mixed | 64^3 | 0.9613 | 0.8621 | 0.6194 | 8.0799 | 3.0794 | - | - |
| none | mixed | 128^3 | 0.9610 | 0.8567 | 0.6099 | 8.3500 | 3.6493 | - | - |
| none | mixed_ab | 32^3 | 0.0378 | 0.6078 | 0.9769 | 7.4215 | 1.9977 | - | - |
| none | mixed_ab | 64^3 | 0.0378 | 0.6061 | 0.9743 | 8.0402 | 2.8151 | - | - |
| B | separate | 32^3 | 0.9645 | 0.9779 | 0.5293 | 7.9614 | 6.9490 | 0.9121 | 0.9340 |
| B | separate | 64^3 | 0.9646 | 0.9780 | 0.5358 | 8.2239 | 6.2574 | 0.9123 | 0.9298 |
| B | separate | 128^3 | 0.9647 | 0.9780 | 0.5394 | 8.3124 | 6.1180 | 0.9124 | 0.9290 |
| B | separate_ab | 32^3 | 0.1556 | 0.9517 | 0.9149 | 4.3510 | 4.4139 | 0.5147 | 0.7992 |
| B | separate_ab | 64^3 | 0.1562 | 0.9517 | 0.9133 | 4.5950 | 4.9178 | 0.5145 | 0.8197 |
| B | separate_ab | 128^3 | 0.1572 | 0.9518 | 0.9132 | 4.6848 | 5.2408 | - | - |
| B | mixed | 32^3 | 0.9622 | 0.8776 | 0.6424 | 6.4120 | 1.8713 | 0.8594 | 0.8789 |
| B | mixed | 64^3 | 0.9613 | 0.8621 | 0.6194 | 6.9899 | 2.7331 | 0.8651 | 0.8875 |
| B | mixed | 128^3 | 0.9610 | 0.8567 | 0.6099 | 7.2428 | 3.2675 | 0.8674 | 0.8954 |
| B | mixed_ab | 32^3 | 0.0378 | 0.6078 | 0.9769 | 2.7591 | 1.3172 | 0.3718 | 0.6594 |
| B | mixed_ab | 64^3 | 0.0378 | 0.6061 | 0.9743 | 2.8388 | 1.9320 | 0.3531 | 0.6863 |
| B | mixed_ab | 128^3 | 0.0378 | 0.6056 | 0.9732 | 2.8652 | 2.3120 | - | - |
| A | separate | 32^3 | 0.9645 | 0.9779 | 0.5293 | 7.9610 | 6.9462 | 0.9120 | 0.9336 |
| A | separate | 64^3 | 0.9646 | 0.9780 | 0.5359 | 8.2235 | 6.2567 | 0.9123 | 0.9297 |
| A | separate | 128^3 | 0.9647 | 0.9780 | 0.5394 | 8.3120 | 6.1180 | 0.9124 | 0.9290 |
| A | separate_ab | 32^3 | 0.1557 | 0.9517 | 0.9150 | 4.3514 | 4.4103 | 0.5148 | 0.7985 |
| A | separate_ab | 64^3 | 0.1564 | 0.9517 | 0.9134 | 4.5955 | 4.9203 | 0.5146 | 0.8201 |
| A | separate_ab | 128^3 | 0.1574 | 0.9518 | 0.9133 | 4.6855 | 5.2505 | - | - |
| A | mixed | 32^3 | 0.9622 | 0.8776 | 0.6424 | 6.4110 | 1.8701 | 0.8592 | 0.8784 |
| A | mixed | 64^3 | 0.9613 | 0.8622 | 0.6194 | 6.9895 | 2.7327 | 0.8651 | 0.8874 |
| A | mixed | 128^3 | 0.9610 | 0.8568 | 0.6100 | 7.2411 | 3.2648 | 0.8672 | 0.8946 |
| A | mixed_ab | 32^3 | 0.0378 | 0.6078 | 0.9769 | 2.7589 | 1.3140 | 0.3717 | 0.6578 |
| A | mixed_ab | 64^3 | 0.0378 | 0.6061 | 0.9743 | 2.8386 | 1.9281 | 0.3530 | 0.6849 |
| A | mixed_ab | 128^3 | 0.0378 | 0.6056 | 0.9732 | 2.8649 | 2.3072 | - | - |

### dc/dt (--transient-c 1e8)
| law | case | grid | sum phibar | share A | c | alpha-1 (1e-4) | vM/Pa (1e-4) | alpha-1 / none | vM / none |
|---|---|---|---|---|---|---|---|---|---|
| none | separate_ab | 32^3 | 0.1557 | 0.9517 | 0.9149 | 8.4552 | 5.5287 | - | - |
| none | separate_ab | 64^3 | 0.1563 | 0.9517 | 0.9133 | 8.9327 | 6.0028 | - | - |
| none | mixed_ab | 32^3 | 0.0378 | 0.6081 | 0.9769 | 7.4217 | 1.9978 | - | - |
| none | mixed_ab | 64^3 | 0.0378 | 0.6065 | 0.9743 | 8.0406 | 2.8156 | - | - |
| B | separate | 32^3 | 0.9645 | 0.9779 | 0.5293 | 7.9637 | 6.9465 | - | - |
| B | separate | 64^3 | 0.9646 | 0.9780 | 0.5358 | 8.2263 | 6.2544 | - | - |
| B | separate | 128^3 | 0.9647 | 0.9780 | 0.5394 | 8.3136 | 6.1167 | - | - |
| B | separate_ab | 32^3 | 0.1557 | 0.9517 | 0.9149 | 4.3549 | 4.4213 | 0.5151 | 0.7997 |
| B | separate_ab | 64^3 | 0.1563 | 0.9517 | 0.9133 | 4.5992 | 4.9226 | 0.5149 | 0.8200 |
| B | separate_ab | 128^3 | 0.1573 | 0.9518 | 0.9132 | 4.6870 | 5.2416 | - | - |
| B | mixed | 32^3 | 0.9623 | 0.8781 | 0.6424 | 6.4145 | 1.8723 | - | - |
| B | mixed | 64^3 | 0.9613 | 0.8627 | 0.6194 | 6.9929 | 2.7348 | - | - |
| B | mixed | 128^3 | 0.9610 | 0.8570 | 0.6099 | 7.2442 | 3.2682 | - | - |
| B | mixed_ab | 32^3 | 0.0378 | 0.6081 | 0.9769 | 2.7617 | 1.3186 | 0.3721 | 0.6600 |
| B | mixed_ab | 64^3 | 0.0378 | 0.6065 | 0.9743 | 2.8417 | 1.9343 | 0.3534 | 0.6870 |
| B | mixed_ab | 128^3 | 0.0378 | 0.6058 | 0.9732 | 2.8667 | 2.3132 | - | - |
| A | separate | 32^3 | 0.9645 | 0.9779 | 0.5293 | 7.9633 | 6.9437 | - | - |
| A | separate | 64^3 | 0.9646 | 0.9780 | 0.5358 | 8.2259 | 6.2538 | - | - |
| A | separate | 128^3 | 0.9647 | 0.9780 | 0.5394 | 8.3132 | 6.1167 | - | - |
| A | separate_ab | 32^3 | 0.1558 | 0.9517 | 0.9149 | 4.3553 | 4.4173 | 0.5151 | 0.7990 |
| A | separate_ab | 64^3 | 0.1566 | 0.9518 | 0.9134 | 4.5997 | 4.9249 | 0.5149 | 0.8204 |
| A | separate_ab | 128^3 | 0.1575 | 0.9518 | 0.9133 | 4.6877 | 5.2512 | - | - |
| A | mixed | 32^3 | 0.9623 | 0.8781 | 0.6424 | 6.4135 | 1.8712 | - | - |
| A | mixed | 64^3 | 0.9613 | 0.8627 | 0.6194 | 6.9925 | 2.7343 | - | - |
| A | mixed | 128^3 | 0.9610 | 0.8570 | 0.6100 | 7.2427 | 3.2657 | - | - |
| A | mixed_ab | 32^3 | 0.0378 | 0.6081 | 0.9769 | 2.7615 | 1.3154 | 0.3721 | 0.6584 |
| A | mixed_ab | 64^3 | 0.0378 | 0.6065 | 0.9743 | 2.8415 | 1.9304 | 0.3534 | 0.6856 |
| A | mixed_ab | 128^3 | 0.0378 | 0.6058 | 0.9732 | 2.8664 | 2.3085 | - | - |

### GCI (fine grid 128^3, Fs 1.25, r 2): order p and GCI %
| dc/dt | law | case | alpha-1: p, GCI % | vM: p, GCI % |
|---|---|---|---|---|
| no | B | separate | 1.57, 0.677 | 2.31, 0.718 |
| no | B | separate_ab | 1.44, 1.397 | 0.64, 13.748 |
| no | B | mixed | 1.19, 3.393 | 0.69, 33.349 |
| no | B | mixed_ab | 1.59, 0.570 | 0.69, 33.267 |
| no | A | separate | 1.57, 0.677 | 2.31, 0.713 |
| no | A | separate_ab | 1.44, 1.402 | 0.63, 14.428 |
| no | A | mixed | 1.20, 3.341 | 0.70, 32.807 |
| no | A | mixed_ab | 1.60, 0.564 | 0.70, 33.114 |
| yes | B | separate | 1.59, 0.654 | 2.33, 0.699 |
| yes | B | separate_ab | 1.48, 1.314 | 0.65, 13.327 |
| yes | B | mixed | 1.20, 3.333 | 0.69, 33.063 |
| yes | B | mixed_ab | 1.68, 0.496 | 0.70, 32.794 |
| yes | A | separate | 1.59, 0.653 | 2.33, 0.694 |
| yes | A | separate_ab | 1.48, 1.317 | 0.64, 13.968 |
| yes | A | mixed | 1.21, 3.286 | 0.70, 32.584 |
| yes | A | mixed_ab | 1.68, 0.491 | 0.70, 32.695 |

### dc/dt minus limit, relative (same grid)
| law | case | grid | sum phibar | c | alpha-1 | vM |
|---|---|---|---|---|---|---|
| B | separate | 32^3 | -1.66e-07 | -1.19e-05 | +2.86e-04 | -3.57e-04 |
| B | separate | 64^3 | -2.98e-07 | -2.12e-05 | +2.99e-04 | -4.76e-04 |
| B | separate | 128^3 | -1.48e-07 | -1.04e-05 | +1.51e-04 | -2.18e-04 |
| B | separate_ab | 32^3 | +6.09e-04 | -8.44e-05 | +9.05e-04 | +1.67e-03 |
| B | separate_ab | 64^3 | +7.91e-04 | -7.71e-05 | +9.20e-04 | +9.66e-04 |
| B | separate_ab | 128^3 | +4.62e-04 | -3.37e-05 | +4.63e-04 | +1.69e-04 |
| B | mixed | 32^3 | +1.38e-05 | +1.07e-05 | +3.89e-04 | +5.53e-04 |
| B | mixed | 64^3 | +1.42e-05 | -1.06e-05 | +4.19e-04 | +6.04e-04 |
| B | mixed | 128^3 | +7.57e-06 | +2.30e-06 | +2.00e-04 | +2.12e-04 |
| B | mixed_ab | 32^3 | +1.03e-04 | -1.74e-06 | +9.45e-04 | +1.09e-03 |
| B | mixed_ab | 64^3 | +1.12e-04 | -3.50e-06 | +1.02e-03 | +1.17e-03 |
| B | mixed_ab | 128^3 | +5.80e-05 | -1.09e-06 | +5.27e-04 | +5.12e-04 |
| A | separate | 32^3 | -1.75e-07 | -1.25e-05 | +2.86e-04 | -3.52e-04 |
| A | separate | 64^3 | -2.97e-07 | -2.10e-05 | +3.00e-04 | -4.72e-04 |
| A | separate | 128^3 | -1.52e-07 | -1.07e-05 | +1.51e-04 | -2.19e-04 |
| A | separate_ab | 32^3 | +5.74e-04 | -8.43e-05 | +8.98e-04 | +1.60e-03 |
| A | separate_ab | 64^3 | +7.75e-04 | -7.75e-05 | +9.17e-04 | +9.36e-04 |
| A | separate_ab | 128^3 | +4.52e-04 | -3.31e-05 | +4.60e-04 | +1.24e-04 |
| A | mixed | 32^3 | +1.38e-05 | +1.04e-05 | +3.90e-04 | +5.55e-04 |
| A | mixed | 64^3 | +1.42e-05 | -1.12e-05 | +4.20e-04 | +6.07e-04 |
| A | mixed | 128^3 | +7.28e-06 | -1.03e-05 | +2.19e-04 | +3.00e-04 |
| A | mixed_ab | 32^3 | +1.03e-04 | -1.87e-06 | +9.45e-04 | +1.09e-03 |
| A | mixed_ab | 64^3 | +1.12e-04 | -3.77e-06 | +1.02e-03 | +1.17e-03 |
| A | mixed_ab | 128^3 | +5.79e-05 | -2.27e-06 | +5.30e-04 | +5.84e-04 |
| none | separate_ab | 32^3 | +6.09e-04 | -8.44e-05 | +2.78e-04 | +1.04e-03 |
| none | separate_ab | 64^3 | +7.91e-04 | -7.71e-05 | +1.88e-04 | +5.10e-04 |
| none | mixed_ab | 32^3 | +1.03e-04 | -1.74e-06 | +2.16e-05 | +7.14e-05 |
| none | mixed_ab | 64^3 | +1.12e-04 | -3.50e-06 | +4.28e-05 | +1.62e-04 |

Reading (no conclusion beyond the numbers):
- **A and B agree** within 0.2 % in alpha - 1 and 0.4 % in vM at every grid and case.
- **dc/dt against the limit:** the largest relative difference is 1.7e-3 (vM, separate_ab, 32^3). At 128^3 every
  difference is below 6e-4.
- **Growth by the living amount** lowers alpha - 1 by about 9 % without antibiotic (separate 0.912, mixed 0.86-0.87). With
  antibiotic it lowers alpha - 1 by 49 % (separate_ab, 0.515) and 65 % (mixed_ab, 0.35-0.37).
- **GCI of alpha - 1** is 0.5-3.4 % (p 1.2-1.7).
- **GCI of vM is large** for separate_ab (13-14 %, p 0.64) and for both mixed cases (33 %, p 0.7). There vM still rises by
  20-30 % per refinement, so it is not in the asymptotic range on 32/64/128. Only separate is converged (0.7 %, p 2.3).

## 2026-10-30: 256^3 mixed example with antibiotic, dc/dt (CLOUD_TO_FIFA 30 Oct morning), running

Options as in the instructions (`--nf 256 ... --transient-c 1e8 --only mixed_ab`), run 1 with `--mix-equal`, run 2 with the
Fritsch shares. Run 1 started 30 Oct 16:33 JST (nohup).

**First steps of run 1** (400 outer steps, dt 0.0025; read from the checkpoint `k_done`, which is written every 120 s):
- step 0 done after 2 min 50 s, including setup
- steps 0 to 8 took 1171 s, about **146 s per step**
- RSS 3.8-4.7 GB, **peak 4.7 GB**
- So one run takes about 16 h, alone on fifa. The instructions estimated 22-29 h, as 3.2 times IKMHIWI03's 192^3.

**Run 2 started beside run 1** at 16:57 JST (the user, via the cloud: run both at once if >= 8 cores and >= 10 GB free).
fifa had 12 cores and 43 GB available (free 10 GB). The chain script of run 1 was stopped so that it does not start run 2
a second time; run 1 goes on as its own process. The two checkpoints have different names (`_mixeq`), so they do not clash.
Logs: `~/runs/n256_20261030/` (`chain.log`, `mixeq_ab.out`, `fritsch_ab.out`, `monitor.log`).

## 2026-10-31: 256^3 mixed example with antibiotic, dc/dt: done

Both runs ran in parallel on fifa and ended with no error. Run 1 finished 31 Oct 08:29 JST (57276 s, 15.9 h); its exit code
was not captured, because its chain script had been stopped (see 30 Oct), but the `.out` has the final line and `.err` is
empty. Run 2 finished 31 Oct 09:46 JST (60521 s, 16.8 h) with exit 0.

| shares | sum phibar | share A | c | alpha - 1 | vM [Pa] |
| --- | --- | --- | --- | --- | --- |
| equal (`--mix-equal`) | 0.0378 | 0.6055 | 0.9726 | 8.4728e-4 | 3.6327e-4 |
| Fritsch | 0.0370 | 0.5065 | 0.9729 | 8.4486e-4 | 4.5012e-4 |

The vM of the equal-share run still rises with the grid. With the cloud's 32/64/128/192^3 values (2.00/2.82/3.30/3.52e-4,
CLOUD_TO_FIFA 30 Oct morning), the steps are +0.22e-4 (128 to 192) and +0.11e-4 (192 to 256).
The GCI is computed with Celik 2008 (Fs 1.25), with the iteration for unequal refinement ratios:

| grids | p | GCI | extrapolated vM |
| --- | --- | --- | --- |
| 64/128/256 | 0.53 | 25.9 % | 4.38e-4 |
| 128/192/256 | 0.93 | 12.6 % | 4.00e-4 |

So at 256^3, where the cell Peclet number is 0.96 by the instructions' estimate, the vM is still not in the asymptotic
range, although the order rises toward 1 on the finest three grids. No lower-grid run of the Fritsch shares with dc/dt
is in `results/2026-10-3d/`, so that run has no grid series. The JSONs are in `results/2026-10-3d/` (`_nf256_`).

## 2026-11-02: surrogate pilot, the 12 remaining points at 128^3 (CLOUD_TO_FIFA 1 Nov night)

Command: `python3 papers/surrogate_pilot/run_pilot.py --tag g --nf 128 --reuse papers/surrogate_pilot/design_g_in.json --jobs 11`.
`--jobs 11` is 12 cores minus the 256^3 run that was running beside it. The run started 1 Nov 23:47 JST, ended 2 Nov 01:28 JST,
exit 0; all 12 points report "done".
- **Time per run** (from the finish times of the JSONs): 55-59 min for the 11 that ran together (12 processes on 12 cores).
  The 12th started when the first one finished and took 46 min.
- **Memory per run**: RSS about 1.1 GB, peak 1.23 GB (sampled once a minute).

Files: `papers/surrogate_pilot/runs/g000..g011_*.json` and `papers/surrogate_pilot/design_g.json`, which the script wrote.
The analysis of the ratios is left to the cloud.
The 256^3 mixed run without antibiotic (1 Nov) is still running; its report follows.

## 2026-11-02: 256^3 mixed (equal shares) without antibiotic, dc/dt (CLOUD_TO_FIFA 1 Nov): done

Started 1 Nov 22:45 JST and ended 2 Nov 15:10 JST, exit 0, `.err` empty, 59039 s (16.4 h). The 128^3 surrogate pilot runs
(11 jobs) ran beside it from 23:47 to 01:28.

Seed means at the end: sum phibar 0.9608, share A 0.8540, c 0.6047, alpha - 1 8.5038e-4, **vM 4.0167e-4 Pa**.

vM against the grid. The coarser values are the limit values of `tab:ch5-gci` (2.13/3.08/3.65/3.89e-4 at 32/64/128/192^3);
dc/dt and the limit differ by less than 0.1 %. The steps are +0.24e-4 (128 to 192) and +0.13e-4 (192 to 256).
GCI by Celik 2008 (Fs 1.25), with the iteration for unequal ratios:

| grids | p | GCI | extrapolated vM |
| --- | --- | --- | --- |
| 128/192/256 | 0.85 | 14.3 % | 4.48e-4 |
| 64/128/256 | 0.64 | 20.6 % | 4.68e-4 |

Against the antibiotic run at 256^3 (3.6327e-4, 31 Oct), the run without antibiotic is **10.6 % higher in vM**
(ratio 1.106). The alpha - 1 is 8.504e-4 here and 8.473e-4 with antibiotic.
File: `ansys_usermat/apdl/results/2026-10-3d/coupled_living_consumption_nf256_..._mixeq_..._transc1e+08_..._mixed.json`.

## 2026-11-03: CLOUD_TO_FIFA 2 Nov night, late night, late night part 2, 3 Nov early morning (superbee, asymmetric 64^3), fifa

Branch `claude/plan-next-hxjjve`, run in a separate worktree (`~/IKM_Hiwi/wt_cloud1102`, branch `cloud-1102`), so
the main checkout stayed on `keio/uel-implant`. The superbee 128^3 run and the first asymmetric 64^3 run went in parallel.
The first start of the superbee 128^3 run failed after 2 s: both runs built `libecology_native.so` in the fresh worktree
at the same moment. It was restarted alone.
**First asymmetric 64^3 run, without `--amg`: stopped.** Started 2 Nov 23:44 and stopped by hand 3 Nov 06:41 JST (the user
agreed), after 6 h 57 min at 2 of 400 outer steps (100 % of one core, RSS 7.3 GB). Without `--fast` and without `--amg` the
fields are solved directly (spsolve) on the full domain. The cloud's 3 Nov early-morning correction added `--amg`.
Order after that: superbee 32/64/192^3, then the edge run and the corner run with `--amg`, one at a time. Each `--amg` run
had a speed guard: stop if fewer than 15 of 400 steps were done after 15 min.
The worktree was fast-forwarded to `42234db` and then `6ac609e`; the solver changes add output for more than 2 species and
the option `--b-corner` (off by default), so they do not touch the other runs.
This section was written by a script (`~/runs/cloud1102/report.py`) when the runs had finished.

**2) Asymmetric example, B along one edge** (`--b-one-edge --seed-shift -1 --amg`). B is kept only in the block along one vertical edge of the box (0.25 x 0.25 x 1.0 mm at x > 0.75, y < -0.75 mm): one edge block, not a corner. The cloud chose this example; it was not Oliver's request.
Both runs: `--nf 64`, full domain, `--transient-c 1e8 --only separate_ab --snap-every 0.025 --snap-stress`.
Wall time 17:22.38, peak RSS 1.35 GB. b1edge64amg speed check: 32 of 400 steps after 900 s.
- Seed means at the end: sum phibar 0.1559, share A 0.9518, c 0.9181, alpha - 1 8.2161e-04, vM 3.9396e-04 Pa.
- B block (the JSON's `edges`) at the end: sum phibar 0.0412, share A 0.0530. The JSON has no c, alpha - 1 or vM for the B block.
JSON: `ansys_usermat/apdl/results/2026-10-3d/coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_snapsvm_transc1e+08_b1edge_seedx-1_tvd_wgmirror_amg_separate_ab.json`.
Fields for the video (not pushed, on fifa): `/home/nishioka/IKM_Hiwi/wt_cloud1102/figdata/coupled_living_separate_ab_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_snapsvm_transc1e+08_b1edge_seedx-1_tvd_wgmirror_amg_separate_ab.npz` (261 MB).
No video was made on fifa (left to the cloud).

**Part 2: B in one corner** (`--b-corner --amg`). B is only in one cube of 0.25 mm at the corner x > 0.75, y < -0.75, z > 0.75 mm of the box, next to the nutrient layer.
Both runs: `--nf 64`, full domain, `--transient-c 1e8 --only separate_ab --snap-every 0.025 --snap-stress`.
Wall time 17:25.14, peak RSS 1.36 GB. bcorner64amg speed check: 31 of 400 steps after 900 s.
- Seed means at the end: sum phibar 0.1579, share A 0.9520, c 0.9170, alpha - 1 8.3385e-04, vM 3.4117e-04 Pa.
- B block (the JSON's `edges`) at the end: sum phibar 0.0412, share A 0.0530. The JSON has no c, alpha - 1 or vM for the B block.
JSON: `ansys_usermat/apdl/results/2026-10-3d/coupled_living_consumption_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_snapsvm_transc1e+08_bcorner_tvd_wgmirror_amg_separate_ab.json`.
Fields for the video (not pushed, on fifa): `/home/nishioka/IKM_Hiwi/wt_cloud1102/figdata/coupled_living_separate_ab_nf64_front1_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_snapsvm_transc1e+08_bcorner_tvd_wgmirror_amg_separate_ab.npz` (317 MB).
No video was made on fifa (left to the cloud).

**1) and late night: superbee grid series** (mixed, equal shares, antibiotic, `--dtc 0.0025`, `--fast`). minmod: `tab:ch5-gci`.

| grid | vM superbee | vM minmod | alpha - 1 superbee | alpha - 1 minmod | wall time, peak RSS (superbee) |
| --- | --- | --- | --- | --- | --- |
| 32^3 | 2.3875e-04 | 2.00e-04 | 7.7879e-04 | 7.42e-04 | 0:42.70, 0.55 GB |
| 64^3 | 3.2841e-04 | 2.82e-04 | 8.3666e-04 | 8.04e-04 | 5:05.19, 0.63 GB |
| 128^3 | 3.5866e-04 | 3.30e-04 | 8.4946e-04 | 8.32e-04 | 1:00:45, 1.30 GB |
| 192^3 | 3.6745e-04 | 3.52e-04 | 8.5125e-04 | 8.42e-04 | 5:09:23, 3.22 GB |

GCI by Celik 2008 (Fs 1.25) on 64/128/192^3 (r21 = 1.5, r32 = 2):

| quantity | limiter | p | GCI | extrapolated |
| --- | --- | --- | --- | --- |
| vM | superbee | 1.24 | 4.6 % | 3.809e-04 |
| vM | minmod | 0.44 | 40.1 % | 4.648e-04 |
| alpha - 1 | superbee | 2.46 | 0.2 % | 8.523e-04 |
| alpha - 1 | minmod | 0.88 | 3.5 % | 8.653e-04 |

(The cloud's minmod figures, p 0.43, extrapolated 4.68 and GCI 41 %, come from a different grid triple; the row above uses 64/128/192^3 for both limiters.)

**Videos (fifa, 2026-11-03).** Rendered off screen with pyvista 0.49.2 / VTK 9.7.1 under xvfb-run, in a separate
venv (`~/IKM_Hiwi/.venv_render`), so the solver venv was not changed:
- Video `b1edge seedx-1`: `assets/anim_coupled_tvd_separate_ab_inhliv_b1edge_seedx-1_transc_alpha_5panel.mp4` (4.4M) and `_end.png`, 70 s.
- Video `bcorner`: `assets/anim_coupled_tvd_separate_ab_inhliv_bcorner_transc_alpha_5panel.mp4` (4.2M) and `_end.png`, 68 s.

## 2026-11-04: CLOUD_TO_FIFA 4 Nov, superbee runs (fifa, written by ~/runs/sb1104/watch.py)

All runs as separate PBS jobs (1 core each) on the compute nodes, started 4 Nov 08:52 JST, in parallel.

**2) separate, antibiotic** (`--limiter superbee --dtc 0.0025 --fast`, no `--transient-c`). Compare: minmod (thesis table): vM 6.000 / 6.284e-4, alpha - 1 8.931 / 9.112e-4 on 64 / 128^3.

| run | grid | seed vM [Pa] | seed alpha - 1 | wall, peak RSS | JSON |
| --- | --- | --- | --- | --- | --- |
| `sep_ab_64` | 64^3 | 6.1223e-04 | 9.1614e-04 | 25:10.61, 0.61 GB | `coupled_living_consumption_nf64_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_separate_ab.json` |
| `sep_ab_128` | 128^3 | 6.3382e-04 | 9.2092e-04 | 4:51:07, 1.27 GB | `coupled_living_consumption_nf128_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_separate_ab.json` |



## 2026-11-04 (evening): CLOUD_TO_FIFA 4 Nov part 2 c), part 3 and part 5 (fifa)

**Part 2 c) maze** (`maze_coupled_2d.py --tend 60 --maze-file ../../abaqus_living/maze_k24.txt --transient-c`,
with and without `--ab`). Values at T* = 20 and 60; area = colonised cells x h^2.

| h [um] | antibiotic | T* | reach [um] | living mean | c mean | area [um^2] |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | no | 20 | 71 | 0.9383 | 0.2078 | 46928 |
| 4 | no | 60 | 42 | 0.9194 | 0.1119 | 72384 |
| 2 | no | 20 | 68 | 0.9359 | 0.2012 | 48320 |
| 2 | no | 60 | 40 | 0.9176 | 0.1101 | 74448 |
| 1 | no | 20 | 67 | 0.9350 | 0.1983 | 47924 |
| 1 | no | 60 | 39 | 0.9179 | 0.1099 | 75793 |
| 4 | yes | 20 | 71 | 0.9383 | 0.2078 | 46928 |
| 4 | yes | 60 | 133 | 0.3467 | 0.4070 | 73408 |
| 2 | yes | 20 | 68 | 0.9359 | 0.2012 | 48320 |
| 2 | yes | 60 | 43 | 0.3408 | 0.3902 | 74928 |
| 1 | yes | 20 | 67 | 0.9350 | 0.1983 | 47924 |
| 1 | yes | 60 | 42 | 0.3259 | 0.3827 | 78812 |

Reading (maze): at T* = 20 the reach is 71 / 68 / 67 um for h = 4 / 2 / 1 um, the living mean 0.938 / 0.936 / 0.935 and the
area 46 900 / 48 300 / 47 900 um^2, so here the three cells agree to within 6 %. The cloud's check of 4 Nov had 65, 195 and
71 um at T* = 20; fifa's runs do not show the 195 um. The one outlier is the reach with antibiotic at T* = 60 for h = 4 um
(133 um against 43 and 42 um): 'reach' is the smallest path distance of a cell with phi >= 0.5, so one far cell is
enough to change it; the area and the living mean agree (73 400 / 74 900 / 78 800 um^2, 0.347 / 0.341 / 0.326).

**Part 3: four limiters** (mixed, equal shares, antibiotic, `--dtc 0.0025 --fast`, no `--transient-c`; minmod = the
thesis series, superbee = 3 Nov, van Leer and MC = this run).

| limiter | vM 32 / 64 / 128 / 192^3 [1e-4 Pa] | alpha - 1 32 / 64 / 128 / 192^3 [1e-4] | vM: p, extrapolated, GCI (64/128/192) | alpha - 1: p, extrapolated, GCI |
| --- | --- | --- | --- | --- |
| minmod | 1.998 / 2.815 / 3.300* / 3.520 | 7.422 / 8.040 / 8.320* / 8.421 | 0.46, 4.5883e-04, 38.0 % | 0.86, 8.6613e-04, 3.6 % |
| vanleer | 2.156 / 2.985 / 3.404 / 3.579 | 7.585 / 8.168 / 8.383 / 8.458 | 0.60, 4.2143e-04, 22.2 % | 0.92, 8.6246e-04, 2.5 % |
| mc | 2.264 / 3.071 / 3.439 / 3.609 | 7.687 / 8.230 / 8.405 / 8.475 | 0.43, 4.5018e-04, 30.9 % | 0.69, 8.6912e-04, 3.2 % |
| superbee | 2.387 / 3.284 / 3.587 / 3.675 | 7.788 / 8.367 / 8.495 / 8.513 | 1.24, 3.8095e-04, 4.6 % | 2.46, 8.5230e-04, 0.2 % |

* thesis table value (tab:ch5-gci); fifa has no minmod 128^3 JSON with step 0.0025.
Extrapolated seed vM of the four limiters: 3.809e-04 to 4.588e-04 Pa (spread 20.4 % of the lowest); GCIs minmod 38.0 %, vanleer 22.2 %, mc 30.9 %, superbee 4.6 %.

Reading (limiters): the four extrapolated values do not agree (3.81 to 4.59e-4 Pa, 20 % apart), but only superbee is in
the asymptotic range on 64/128/192^3 (p 1.24, GCI 4.6 %); minmod, van Leer and MC have p 0.43-0.60 and GCI 22-38 %, so their
extrapolations are not reliable. On each grid the order is minmod < van Leer < MC < superbee, and the gaps shrink with
refinement (192^3: 3.52 / 3.58 / 3.61 / 3.68e-4). alpha - 1 agrees within 1 % at 192^3 for all four.

**Part 5: AMR** (`coupled_amr.py`, mixed, equal shares, antibiotic; `ECO_THREADS=4`, 4 cores each).

| run | seed alpha - 1 | seed vM [Pa] | against | wall [s] | time split [s] (base / fine / nutrient / point / retag / restrict) | finest cells, mean share of (2N)^3 |
| --- | --- | --- | --- | --- | --- | --- |
| L3 minmod (32^3 base, 3 levels = 128^3) | 8.3450e-04 | 3.3417e-04 | uniform 128^3 minmod: alpha - 1 +0.30 %, vM +1.26 % | 4970 | 8 / 1089 / 2127 / 1593 / 55 / 57 | 54.7 % |
| L3 superbee | 8.5246e-04 | 3.6321e-04 | uniform 128^3 superbee: alpha - 1 +0.35 %, vM +1.26 % | 2341 | 3 / 594 / 978 / 693 / 26 / 22 | 54.7 % |
| 64^3 base, c every 1 step | 8.3520e-04 | 3.3509e-04 | (reference) | 1702 | 20 / 270 / 630 / 741 / 17 / 11 | 54.7 % |
| 64^3 base, c every 2 steps | 8.4225e-04 | 3.4588e-04 | k = 1: alpha - 1 +0.84 %, vM +3.22 % | 852 | 12 / 182 / 246 / 385 / 12 / 5 | 54.7 % |
| 64^3 base, c every 4 steps | 8.5560e-04 | 3.6660e-04 | k = 1: alpha - 1 +2.44 %, vM +9.40 % | 787 | 12 / 181 / 182 / 384 / 12 / 5 | 54.6 % |

Reading (AMR): 3 levels on a 32^3 base give the uniform 128^3 within 0.3 % (alpha - 1) and 1.3 % (vM) for both limiters.
Updating the nutrient every 2 / 4 steps halves the wall time (1702 -> 852 / 787 s) but moves vM by +3.2 / +9.4 %.
For the 3-level runs the 'share of (2N)^3' is the share of the finest level against a uniform grid of twice the base.
Part 4 a) (128^3 base), the obstacle runs of part 4 b) and WENO5 192^3 are still running; their tables follow.
**1) mixed, equal shares, no antibiotic** (`--limiter superbee --dtc 0.0025 --fast`, no `--transient-c`). Compare: superbee with antibiotic (3 Nov): vM 2.39/3.28/3.59/3.67e-4, extrapolated 3.81e-4, GCI 4.6 %; minmod without antibiotic (thesis table): 3.08/3.65/3.89e-4 on 64/128/192^3.

| run | grid | seed vM [Pa] | seed alpha - 1 | wall, peak RSS | JSON |
| --- | --- | --- | --- | --- | --- |
| `mix_32` | 32^3 | 2.5729e-04 | 7.8406e-04 | 1:03.91, 0.53 GB | `coupled_living_consumption_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_mixed.json` |
| `mix_64` | 64^3 | 3.6361e-04 | 8.4145e-04 | 22:56.76, 0.61 GB | `coupled_living_consumption_nf64_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_mixed.json` |
| `mix_128` | 128^3 | 3.9919e-04 | 8.5317e-04 | 4:39:29, 1.22 GB | `coupled_living_consumption_nf128_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_mixed.json` |
| `mix_192` | 192^3 | 4.0721e-04 | 8.5458e-04 | 19:22:53, 3.00 GB | `coupled_living_consumption_nf192_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_mixed.json` |

- 64/128/192^3 vM: p 1.67, extrapolated 4.1551e-04, GCI 2.5 %
- 64/128/192^3 alpha - 1: p 2.71, extrapolated 8.5529e-04, GCI 0.1 %
- Antibiotic effect on vM (superbee, extrapolated): without 4.155e-04, with 3.81e-4 (3 Nov); difference 3.451e-05 = 9.1 % of the value with antibiotic; sum of the two GCIs 2.5 + 4.6 = 7.1 %. The difference is larger than the sum: the antibiotic effect is resolved.

**3) mixed, Fritsch shares, antibiotic** (`--limiter superbee --dtc 0.0025 --fast`, no `--transient-c`). Compare: minmod (thesis table): vM 2.01/2.96/4.01/4.35e-4, extrapolated 5.01e-4, GCI 19 %.

| run | grid | seed vM [Pa] | seed alpha - 1 | wall, peak RSS | JSON |
| --- | --- | --- | --- | --- | --- |
| `frab_32` | 32^3 | 2.3930e-04 | 7.7838e-04 | 1:02.32, 0.53 GB | `coupled_living_consumption_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_mixed_ab.json` |
| `frab_64` | 64^3 | 3.4551e-04 | 8.3647e-04 | 21:45.40, 0.61 GB | `coupled_living_consumption_nf64_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_mixed_ab.json` |
| `frab_128` | 128^3 | 4.4056e-04 | 8.4851e-04 | 4:37:37, 1.29 GB | `coupled_living_consumption_nf128_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_mixed_ab.json` |
| `frab_192` | 192^3 | 4.5174e-04 | 8.4887e-04 | 19:39:23, 3.23 GB | `coupled_living_consumption_nf192_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_mixed_ab.json` |

- 64/128/192^3 vM: p 2.75, extrapolated 4.5721e-04, GCI 1.5 %
- 64/128/192^3 alpha - 1: p 4.90, extrapolated 8.4893e-04, GCI 0.0 %

**4) mixed, equal shares, antibiotic, 256^3** (`--limiter superbee --dtc 0.0025 --fast`, no `--transient-c`). Compare: superbee equal shares with antibiotic (3 Nov): 64/128/192^3 extrapolated 3.81e-4, GCI 4.6 %; 192^3 vM 3.67e-4.

| run | grid | seed vM [Pa] | seed alpha - 1 | wall, peak RSS | JSON |
| --- | --- | --- | --- | --- | --- |
| `mixab_256` | 256^3 | 3.7574e-04 | 8.5427e-04 | 31:07:48, 6.72 GB | `coupled_living_consumption_nf256_front1_dtc0.0025_phimin0.001_v2_fortran_ansysbc_inhliv_mixeq_a100_gab1_gc1_tab0.333333_tvd_superbee_wgmirror_fast_mixed_ab.json` |

- 128/192/256^3 vM (128 and 192 from 3 Nov): p 1.27, extrapolated 3.9457e-04, GCI 6.3 %; 64/128/192^3 vM: extrapolated 3.81e-4, GCI 4.6 %
- 128/192/256^3 alpha - 1 (128 and 192 from 3 Nov): p 4.02, extrapolated 8.5566e-04, GCI 0.2 %; 64/128/192^3 vM: extrapolated 3.81e-4, GCI 4.6 %


## 2026-11-06: CLOUD_TO_FIFA 4 Nov part 2 a) four species, part 2 b) and part 4 b) obstacle (fifa)

All with `--transient-c 1e8 --dtc 0.0025 --fast`, the common `$B`; one PBS job per run.

**Part 2 a) four species** (`--case k26c4`, Klempt 2026 case 4), antibiotic; separate and mixed (equal shares).

*separate_ab, minmod*

| grid | sum phibar | share A | phibar_i (species 1-4) | c | alpha - 1 | vM [Pa] |
| --- | --- | --- | --- | --- | --- | --- |
| 32^3 | 0.7677 | 0.1708 | 0.0022 / 0.0072 / 0.0017 / 0.7565 | 0.5888 | 8.8456e-04 | 7.0312e-04 |
| 64^3 | 0.7715 | 0.1671 | 0.0021 / 0.0071 / 0.0017 / 0.7606 | 0.5947 | 9.1065e-04 | 6.5833e-04 |
| 128^3 | 0.7740 | 0.1648 | 0.0021 / 0.0070 / 0.0017 / 0.7631 | 0.5981 | 9.2028e-04 | 6.5202e-04 |
| 192^3 | 0.7749 | 0.1639 | 0.0020 / 0.0070 / 0.0017 / 0.7641 | 0.5992 | 9.2289e-04 | 6.5204e-04 |

- 64/128/192^3 vM: not monotone; alpha - 1: p 1.36, extrapolated 9.2644e-04, GCI 0.5 %

*separate_ab, superbee*

| grid | sum phibar | share A | phibar_i (species 1-4) | c | alpha - 1 | vM [Pa] |
| --- | --- | --- | --- | --- | --- | --- |
| 32^3 | 0.7666 | 0.1719 | 0.0023 / 0.0072 / 0.0017 / 0.7554 | 0.5890 | 9.0982e-04 | 6.9826e-04 |
| 64^3 | 0.7710 | 0.1676 | 0.0021 / 0.0071 / 0.0017 / 0.7601 | 0.5949 | 9.2847e-04 | 6.5507e-04 |
| 128^3 | 0.7738 | 0.1649 | 0.0021 / 0.0070 / 0.0017 / 0.7630 | 0.5983 | 9.2907e-04 | 6.5172e-04 |
| 192^3 | 0.7748 | 0.1639 | 0.0020 / 0.0070 / 0.0017 / 0.7641 | 0.5994 | 9.2779e-04 | 6.5192e-04 |

- 64/128/192^3 vM: not monotone; alpha - 1: not monotone

*mixed_ab, minmod*

| grid | sum phibar | share A | phibar_i (species 1-4) | c | alpha - 1 | vM [Pa] |
| --- | --- | --- | --- | --- | --- | --- |
| 32^3 | 0.8358 | 0.1066 | 0.0006 / 0.0061 / 0.0015 / 0.8277 | 0.6676 | 7.6446e-04 | 2.2861e-04 |
| 64^3 | 0.8208 | 0.1203 | 0.0008 / 0.0066 / 0.0015 / 0.8119 | 0.6459 | 8.3839e-04 | 3.4324e-04 |
| 128^3 | 0.8122 | 0.1281 | 0.0010 / 0.0067 / 0.0016 / 0.8029 | 0.6345 | 8.7960e-04 | 4.2853e-04 |
| 192^3 | 0.8100 | 0.1301 | 0.0010 / 0.0067 / 0.0016 / 0.8006 | 0.6316 | 8.9150e-04 | 4.5716e-04 |

- 64/128/192^3 vM: p 0.99, extrapolated 5.1530e-04, GCI 15.9 %; alpha - 1: p 1.25, extrapolated 9.0959e-04, GCI 2.5 %

*mixed_ab, superbee*

| grid | sum phibar | share A | phibar_i (species 1-4) | c | alpha - 1 | vM [Pa] |
| --- | --- | --- | --- | --- | --- | --- |
| 32^3 | 0.8279 | 0.1141 | 0.0007 / 0.0065 / 0.0015 / 0.8193 | 0.6578 | 7.9898e-04 | 2.6904e-04 |
| 64^3 | 0.8132 | 0.1273 | 0.0010 / 0.0067 / 0.0016 / 0.8039 | 0.6372 | 8.6803e-04 | 3.9313e-04 |
| 128^3 | 0.8078 | 0.1322 | 0.0011 / 0.0068 / 0.0016 / 0.7983 | 0.6293 | 8.9714e-04 | 4.6397e-04 |
| 192^3 | 0.8074 | 0.1326 | 0.0011 / 0.0067 / 0.0016 / 0.7979 | 0.6285 | 9.0199e-04 | 4.7885e-04 |

- 64/128/192^3 vM: p 1.79, extrapolated 4.9281e-04, GCI 3.6 %; alpha - 1: p 2.18, extrapolated 9.0541e-04, GCI 0.5 %

Reading (four species): the separate example is converged: vM 6.520 / 6.520e-4 at 128 / 192^3 with either limiter
("not monotone" only because the change is below the last digits), alpha - 1 within 0.5 %. The mixed example still
rises: minmod 4.29 -> 4.57e-4 (GCI 15.9 %), superbee 4.64 -> 4.79e-4 (p 1.79, extrapolated 4.93e-4, GCI 3.6 %), as in
the two-species mixed example; species 4 holds 97-99 % of the living amount in both.

**Part 2 b) / part 4 b) obstacle** (`--obstacle columns`, four rigid columns), antibiotic. `nm` = stress mesh (`--nm`, default 16).

| case | limiter | nm | vM 32 / 64 / 128 / 192^3 [1e-4 Pa] | alpha - 1 32 / 64 / 128 / 192^3 [1e-4] | vM 64/128/192^3 | alpha - 1 64/128/192^3 |
| --- | --- | --- | --- | --- | --- | --- |
| separate_ab | minmod | 16 | 25.228 / 37.309 / 50.363 / 53.965 | 6.916 / 7.294 / 7.641 / 7.732 | p 1.33, extrapolated 5.9023e-03, GCI 11.7 % | p 1.42, extrapolated 7.8488e-04, GCI 1.9 % |
| separate_ab | minmod | 32 | - / 36.121 / 48.738 / - | - / 7.294 / 7.641 / - | - | - |
| separate_ab | superbee | 16 | 26.890 / 42.120 / 54.557 / - | 7.068 / 7.515 / 7.789 / - | - | - |
| mixed_ab | minmod | 16 | 9.334 / 16.136 / 22.115 / 24.579 | 5.854 / 6.445 / 6.828 / 6.967 | p 0.63, extrapolated 3.3082e-03, GCI 43.2 % | p 0.84, extrapolated 7.3103e-04, GCI 6.2 % |
| mixed_ab | superbee | 16 | 11.316 / 19.806 / 24.930 / - | 6.135 / 6.759 / 7.016 / - | - | - |

- Stress mesh 32^3 against 16^3 (separate, minmod): vM lower by 3 % at both 64^3 and 128^3 fields, so the rise with the
  field grid comes from the fields near the columns, not from the stress mesh (as the cloud expected).
- `--nm 64` (separate, minmod, fields 64^3 and 128^3): both runs finished the field part (1297 s, 14 710 s) and were killed
  (signal 9) in the stress solve at 36.5 GB and 43.5 GB peak RSS; no JSON. The 64^3 one died at 4 Nov about 20:46 JST, just
  before frontale01's pbs_mom logged a shutdown (20:47:21); the memory pressure of this job may have been part of it.
- The obstacle 192^3 mixed run (minmod) was dropped from the instructions later but had already started; it is reported here.


## 2026-11-06: CLOUD_TO_FIFA 6 Nov part 9, asymmetric starts at 128^3 (fifa)

Common `$B`, `--transient-c 1e8 --only separate_ab`, the default outer step (no `--dtc`, as the 64^3 asymmetric runs),
`--nm 16`. One PBS job per run, 4 cores, `ECO_THREADS=4`.

| run | grid | sum phibar | share A | c | alpha - 1 | vM [Pa] | vM vs symmetric | wall, peak RSS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B along one edge (`--b-one-edge --seed-shift -1`, full domain, `--amg`) | 64^3 | 0.1559 | 0.9518 | 0.9181 | 8.2161e-04 | 3.9396e-04 | -30.7 % | (earlier run) |
| B along one edge (`--b-one-edge --seed-shift -1`, full domain, `--amg`) | 128^3 | 0.1509 | 0.9511 | 0.9154 | 8.6208e-04 | 4.8719e-04 | -17.4 % | 3:44:18, 3.57 GB |
| B in one corner (`--b-corner`, full domain, `--amg`) | 64^3 | 0.1579 | 0.9520 | 0.9170 | 8.3385e-04 | 3.4117e-04 | -40.0 % | (earlier run) |
| B in one corner (`--b-corner`, full domain, `--amg`) | 128^3 | 0.1530 | 0.9514 | 0.9143 | 8.7502e-04 | 4.2911e-04 | -27.3 % | 3:41:36, 3.55 GB |
| symmetric (`--fast`, quarter domain) | 64^3 | 0.1516 | 0.9511 | 0.9114 | 8.9517e-04 | 5.6871e-04 | - | (earlier run) |
| symmetric (`--fast`, quarter domain) | 128^3 | 0.1520 | 0.9512 | 0.9110 | 9.1594e-04 | 5.9009e-04 | - | 54:02.01, 1.28 GB |

The 64^3 asymmetric rows are the 3 Nov `--amg` runs (with `--snap-every 0.025 --snap-stress`); the instructions quote
3.94e-4 and 3.41e-4 against the symmetric 5.69e-4 (-31 % and -40 %).

## 2026-11-07: CLOUD_TO_FIFA 6 Nov part 9 addition (fifa proposals 1 and 2)

**1) Obstacle, superbee 192^3** (`--obstacle columns --transient-c 1e8 --dtc 0.0025 --fast`, stress mesh 16^3), against minmod.

| case | limiter | vM 32 / 64 / 128 / 192^3 [1e-4 Pa] | alpha - 1 32 / 64 / 128 / 192^3 [1e-4] | vM 64/128/192^3 | alpha - 1 64/128/192^3 |
| --- | --- | --- | --- | --- | --- |
| separate_ab | minmod | 25.228 / 37.309 / 50.363 / 53.965 | 6.916 / 7.294 / 7.641 / 7.732 | p 1.33, extrapolated 5.9023e-03, GCI 11.7 % | p 1.42, extrapolated 7.8488e-04, GCI 1.9 % |
| separate_ab | superbee | 26.890 / 42.120 / 54.557 / 56.967 | 7.068 / 7.515 / 7.789 / 7.827 | p 1.92, extrapolated 5.9007e-03, GCI 4.5 % | p 2.49, extrapolated 7.8485e-04, GCI 0.3 % |
| mixed_ab | minmod | 9.334 / 16.136 / 22.115 / 24.579 | 5.854 / 6.445 / 6.828 / 6.967 | p 0.63, extrapolated 3.3082e-03, GCI 43.2 % | p 0.84, extrapolated 7.3103e-04, GCI 6.2 % |
| mixed_ab | superbee | 11.316 / 19.806 / 24.930 / 26.052 | 6.135 / 6.759 / 7.016 / 7.062 | p 1.72, extrapolated 2.7167e-03, GCI 5.3 % | p 2.06, extrapolated 7.0971e-04, GCI 0.6 % |

Wall time, peak RSS: `ob_sb_sep_192` 8:29:24, 3.14 GB; `ob_sb_mix_192` 7:20:24, 3.15 GB

**2) Four limiters, mixed, equal shares, antibiotic**, all from fifa JSONs (minmod 128^3 now run on fifa;
no `--transient-c`, `--dtc 0.0025`).

| limiter | vM 32 / 64 / 128 / 192^3 [1e-4 Pa] | alpha - 1 [1e-4] | vM 64/128/192^3 | alpha - 1 64/128/192^3 |
| --- | --- | --- | --- | --- |
| minmod | 1.998 / 2.815 / 3.298 / 3.520 | 7.422 / 8.040 / 8.316 / 8.421 | p 0.43, extrapolated 4.6746e-04, GCI 41.0 % | p 0.77, extrapolated 8.7094e-04, GCI 4.3 % |
| vanleer | 2.156 / 2.985 / 3.404 / 3.579 | 7.585 / 8.168 / 8.383 / 8.458 | p 0.60, extrapolated 4.2143e-04, GCI 22.2 % | p 0.92, extrapolated 8.6246e-04, GCI 2.5 % |
| mc | 2.264 / 3.071 / 3.439 / 3.609 | 7.687 / 8.230 / 8.405 / 8.475 | p 0.43, extrapolated 4.5018e-04, GCI 30.9 % | p 0.69, extrapolated 8.6912e-04, GCI 3.2 % |
| superbee | 2.387 / 3.284 / 3.587 / 3.675 | 7.788 / 8.367 / 8.495 / 8.513 | p 1.24, extrapolated 3.8095e-04, GCI 4.6 % | p 2.46, extrapolated 8.5230e-04, GCI 0.2 % |

minmod 128^3 on fifa: 1:27:44, 1.27 GB; the thesis table had vM 3.30e-4, alpha - 1 8.32e-4 for this cell.

## 2026-11-07/08: part 9 at 192^3 (fifa addition with the user's OK, 7 Nov; not in the instructions)

Same commands as part 9 with `--nf 192` (asymmetric: full domain, `--amg`; symmetric: `--fast`); 4 cores each, `ECO_THREADS=4`.

| run | vM 64 / 128 / 192^3 [1e-4 Pa] | vs symmetric 64 / 128 / 192^3 | vM 64/128/192^3: p, extrapolated, GCI | wall, peak RSS (192^3) |
| --- | --- | --- | --- | --- |
| B along one edge | 3.940 / 4.872 / - | -30.7 % / -17.4 % / - | - | ? |
| B in one corner | 3.412 / 4.291 / - | -40.0 % / -27.3 % / - | - | ? |
| symmetric | 5.687 / 5.901 / 5.967 | - | p 1.14, 6.0784e-04, GCI 2.3 % | 3:32:05, 3.22 GB |

Difference of the extrapolated values to the symmetric one: .
