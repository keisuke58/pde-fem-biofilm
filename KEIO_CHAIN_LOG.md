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
