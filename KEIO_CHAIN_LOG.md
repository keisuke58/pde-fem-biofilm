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
