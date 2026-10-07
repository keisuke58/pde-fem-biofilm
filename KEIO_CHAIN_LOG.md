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
