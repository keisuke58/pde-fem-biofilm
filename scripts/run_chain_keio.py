"""run_chain_keio.py -- the Linux/Abaqus counterpart of
ansys_usermat/apdl/run_chain.ps1, for the Keio server (fifa): run a list of
jobs one after another in a process detached from the calling shell (a closed
terminal or a restarted Claude session does not stop it), summarise each one
against its reference, commit and push the summaries, and report to a pull
request so the decision about what to run next is made there.

    python3 scripts/run_chain_keio.py --name 1008 --runs scripts/keio_runs/1008.json \
        [--pr 58] [--push] [--cpus 4] [--dry-run]

fifa is a shared lab server (about 15 users), so the chain is deliberately
conservative and never runs two jobs at once:

- it waits until the one-minute load average is below --max-load before each
  job, and gives up on the chain if it stays high for --load-timeout minutes;
- it refuses to start a job unless "/" has --min-free-root GiB and the work
  disk has --min-free-work GiB left (Abaqus scratch is redirected to /home by
  ~/abaqus_v6.env, but a full "/" would still stop every user of the machine);
- every Abaqus process runs under `nice -n <--nice>`;
- --cpus is capped at --max-cpus, a third of the machine by default;
- the .odb, .sim and scratch of a finished job are deleted once its summary
  has been written, since the summary is what the repository keeps.

The manifest (--runs) is JSON:

    {"out": "abaqus_composition/results_keio_1008",
     "runs": [{"job": "fx41_fo_p0",
               "make": ["abaqus_composition/make_klempt_inp.py", "{inp}",
                        "--case", "fig4_corner", "--n", "20"],
               "summarise": ["abaqus_composition/compare_paper.py", "{dat}",
                             "--case", "fig4_corner"],
               "reference": "abaqus_composition/results_1006/fx41_fo_p0.txt",
               "case": "2sp_case6", "cpus": 4, "timeout_min": 180}]}

{inp}, {dat} and {json} are replaced by the input file, the .dat of the
finished job and the mesh sidecar make_*_inp.py writes next to the .inp.
"reference" is optional; without it the summary is only recorded, not compared.
"""
import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO = Path(__file__).resolve().parents[1]
logger = logging.getLogger("keio_chain")


def load1() -> float:
    """One-minute load average."""
    return os.getloadavg()[0]


def free_gib(path: Path) -> float:
    """Free space at `path`, in GiB."""
    return shutil.disk_usage(path).free / 2**30


def wait_for_load(max_load: float, timeout_min: float) -> bool:
    """Block until the one-minute load average drops below `max_load`.

    Parameters
    ----------
    max_load : float
        Threshold on the one-minute load average.
    timeout_min : float
        Give up after this many minutes.

    Returns
    -------
    bool
        True if the load dropped in time, False on timeout.
    """
    deadline = time.monotonic() + timeout_min * 60
    while load1() >= max_load:
        if time.monotonic() > deadline:
            logger.error("load stayed at or above %.1f for %.0f min; giving up",
                         max_load, timeout_min)
            return False
        logger.info("load %.2f >= %.1f (shared server): waiting 2 min",
                    load1(), max_load)
        time.sleep(120)
    return True


def check_disk(workroot: Path, min_root: float, min_work: float) -> bool:
    """True if "/" and the work disk both have enough free space."""
    root, work = free_gib(Path("/")), free_gib(workroot)
    if root < min_root:
        logger.error("/ has %.1f GiB free, below %.1f: refusing to run "
                     "(a full / stops every user of this server)",
                     root, min_root)
        return False
    if work < min_work:
        logger.error("%s has %.1f GiB free, below %.1f: refusing to run",
                     workroot, work, min_work)
        return False
    logger.info("disk ok: / %.1f GiB, %s %.1f GiB free", root, workroot, work)
    return True


def sh(argv: List[str], cwd: Path, timeout: Optional[float] = None,
       nice: int = 0, stdin_text: Optional[str] = None) -> Tuple[int, str]:
    """Run `argv` and return (returncode, combined output)."""
    cmd = (["nice", "-n", str(nice)] + argv) if nice else argv
    logger.info("$ %s", " ".join(cmd))
    try:
        p = subprocess.run(cmd, cwd=cwd, timeout=timeout, text=True,
                           input=stdin_text,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except subprocess.TimeoutExpired as exc:
        return 124, "timeout after %s s: %s" % (timeout, exc)
    except OSError as exc:
        return 127, "could not run %s: %s" % (cmd[0], exc)
    return p.returncode, p.stdout or ""


def strip_header(text: str) -> List[str]:
    """Body lines of a summary: no '#' header line, no blank line."""
    return [l for l in text.splitlines() if not l.startswith("#") and l.strip()]


def compare_to_reference(summary: str, reference: Path) -> Tuple[str, str]:
    """Compare a fresh summary with a committed reference.

    Parameters
    ----------
    summary : str
        The summary script's stdout.
    reference : Path
        A committed `results_*/<job>.txt`.

    Returns
    -------
    (verdict, detail)
        verdict is "identical", "differs" or "no reference"; detail names the
        first differing lines, at most ten of them.
    """
    if not reference.is_file():
        return "no reference", "%s not found" % reference
    want, got = strip_header(reference.read_text()), strip_header(summary)
    if want == got:
        return "identical", "%d lines, every digit equal" % len(got)
    diff = ["  line %d:\n    ref %r\n    new %r" % (i + 1, a, b)
            for i, (a, b) in enumerate(zip(want, got)) if a != b]
    if len(want) != len(got):
        diff.insert(0, "  line count %d -> %d" % (len(want), len(got)))
    return "differs", "\n".join(diff[:10]) or "(only the line count differs)"


def write_summary(path: Path, job: str, make: List[str], case: str, cpus: int,
                  seconds: int, name: str, reference: Optional[str],
                  verdict: str, table: str) -> None:
    """Write `results_*/<job>.txt` in the house format of results_1006/1007."""
    mins, secs = divmod(seconds, 60)
    head = [
        "# %s  (%s, %d min %d s, %d CPUs, Fortran point model)"
        % (job, datetime.now().strftime("%Y-%m-%d %H:%M"), mins, secs, cpus),
        "# Keio server fifa (Linux, Abaqus 2024), automatic chain %s." % name,
        "# " + " ".join(make),
        "# abaqus_composition/run_comp.sh abaqus_composition/%s.inp %s %d"
        % (job, case, cpus),
        "# PASS: %s completed (native point model)" % job,
    ]
    if reference:
        head.append("# vs %s: %s" % (reference, verdict))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(head) + "\n" + table)


def cleanup(workroot: Path, job: str) -> None:
    """Delete what the repository does not keep: .odb, .sim and the scratch."""
    wd = workroot / ("comp_%s" % job)
    for f in (wd / ("%s.odb" % job), wd / ("%s.sim" % job)):
        f.unlink(missing_ok=True)
    shutil.rmtree(wd / job, ignore_errors=True)
    for d in Path("/tmp").glob("nishioka_%s_*" % job):
        shutil.rmtree(d, ignore_errors=True)
    scratch = workroot / "scratch"
    shutil.rmtree(scratch, ignore_errors=True)
    scratch.mkdir(exist_ok=True)


def run_one(run: Dict, name: str, out_dir: Path, workroot: Path,
            args: argparse.Namespace) -> Dict:
    """Generate, run, summarise and compare one job.

    Returns
    -------
    dict
        job, cpus, status, seconds, verdict, detail and (on success) file.
    """
    job = run["job"]
    cpus = min(int(run.get("cpus", args.cpus)), args.max_cpus)
    case = run.get("case", "2sp_case6")
    inp = REPO / "abaqus_composition" / ("%s.inp" % job)
    mesh_json = inp.with_suffix(".json")
    dat = workroot / ("comp_%s" % job) / ("%s.dat" % job)
    rec = {"job": job, "cpus": cpus, "status": "", "seconds": 0,
           "verdict": "", "detail": ""}
    fmt = {"inp": str(inp), "dat": str(dat), "json": str(mesh_json)}

    make = [a.format(**fmt) for a in run["make"]]
    rc, log = sh([sys.executable] + make, REPO, timeout=600)
    if rc != 0:
        rec.update(status="make failed (rc=%d)" % rc, detail=log[-2000:])
        return rec
    if args.dry_run:
        rec.update(status="dry run", verdict="skipped")
        return rec
    if not check_disk(workroot, args.min_free_root, args.min_free_work):
        rec["status"] = "skipped: not enough free disk"
        return rec
    if not wait_for_load(args.max_load, args.load_timeout):
        rec["status"] = "skipped: machine too busy"
        return rec

    t0 = time.monotonic()
    rc, log = sh([str(REPO / "abaqus_composition" / "run_comp.sh"),
                  "abaqus_composition/%s.inp" % job, case, str(cpus)],
                 REPO, timeout=float(run.get("timeout_min", 180)) * 60,
                 nice=args.nice)
    rec["seconds"] = int(time.monotonic() - t0)
    sta = workroot / ("comp_%s" % job) / ("%s.sta" % job)
    if not (sta.is_file() and "COMPLETED SUCCESSFULLY" in sta.read_text()):
        rec.update(status="run failed (rc=%d)" % rc, detail=log[-2000:])
        return rec
    rec["status"] = "PASS"

    rc, table = sh([sys.executable] + [a.format(**fmt) for a in run["summarise"]],
                   REPO, timeout=600)
    if rc != 0:
        rec.update(status="summarise failed (rc=%d)" % rc, detail=table[-2000:])
        return rec

    reference = run.get("reference")
    verdict, detail = (compare_to_reference(table, REPO / reference) if reference
                       else ("no reference", "none given"))
    rec.update(verdict=verdict, detail=detail)
    out = out_dir / ("%s.txt" % job)
    write_summary(out, job, make, case, cpus, rec["seconds"], name,
                  reference, verdict, table)
    rec["file"] = str(out.relative_to(REPO))
    if not args.keep:
        cleanup(workroot, job)
    return rec


def report(name: str, recs: List[Dict]) -> str:
    """The pull-request comment: a table for people, JSON for the next session."""
    rows = "\n".join(
        "| `%s` | %s | %d min %d s | %s |"
        % (r["job"], r["status"], r["seconds"] // 60, r["seconds"] % 60,
           r["verdict"] or "-") for r in recs)
    details = "\n\n".join(
        "<details><summary><code>%s</code>: %s</summary>\n\n```\n%s\n```\n</details>"
        % (r["job"], r["verdict"] or r["status"], r["detail"][:3000])
        for r in recs if r["verdict"] == "differs" or "failed" in r["status"])
    payload = json.dumps({"chain": name, "runs": recs}, indent=1)
    return ("## chain `%s` finished\n\n"
            "Automatic report from `scripts/run_chain_keio.py` on fifa.\n\n"
            "| job | status | wall clock | vs reference |\n|---|---|---|---|\n%s\n\n"
            "%s\n\nWall clock is comparable with IKMHIWI03 only if the machine was "
            "otherwise idle: the numbers are deterministic, the timing is not.\n\n"
            "<!-- keio-chain -->\n```json\n%s\n```\n"
            % (name, rows, details, payload))


def commit_push(files: List[str], message: str, branch: str) -> str:
    """Stage exactly `files`, commit, rebase on origin/<branch> and push."""
    if not files:
        return "nothing to commit"
    out = []
    for argv in (["git", "add"] + files,
                 ["git", "commit", "-m", message],
                 ["git", "fetch", "origin", branch],
                 ["git", "rebase", "--autostash", "FETCH_HEAD"],
                 ["git", "push", "origin", "HEAD:%s" % branch]):
        rc, log = sh(argv, REPO, timeout=300)
        out.append("%s: rc=%d %s" % (argv[1], rc, log.strip()[:300]))
        if rc != 0 and argv[1] in ("add", "commit"):
            break
    return "\n".join(out)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Command line of the chain."""
    cpus = os.cpu_count() or 4
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True, help="chain name, used in the log and the report")
    ap.add_argument("--runs", required=True, help="manifest JSON (see the module docstring)")
    ap.add_argument("--pr", default="", help="pull request to comment on (number or URL)")
    ap.add_argument("--push", action="store_true", help="commit and push the summaries")
    ap.add_argument("--cpus", type=int, default=4, help="default CPUs per job")
    ap.add_argument("--max-cpus", type=int, default=max(1, cpus // 3),
                    help="hard cap on CPUs per job (shared server)")
    ap.add_argument("--nice", type=int, default=10, help="nice level for Abaqus")
    ap.add_argument("--max-load", type=float, default=max(2.0, cpus / 2),
                    help="wait while the one-minute load average is at or above this")
    ap.add_argument("--load-timeout", type=float, default=120,
                    help="minutes to wait for the load to drop before giving up")
    ap.add_argument("--min-free-root", type=float, default=8.0,
                    help='GiB that must be free on "/" before a job starts')
    ap.add_argument("--min-free-work", type=float, default=100.0,
                    help="GiB that must be free on the work disk")
    ap.add_argument("--keep", action="store_true", help="keep .odb/.sim and the scratch")
    ap.add_argument("--dry-run", action="store_true", help="generate the inputs only")
    ap.add_argument("--log", default="", help="log file (default $WORKROOT/_chain_<name>.log)")
    ap.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    return ap.parse_args(argv)


def main() -> int:
    args = parse_args()
    workroot = Path(os.environ.get("WORKROOT", Path.home() / "abaqus_work"))
    log_path = Path(args.log) if args.log else workroot / ("_chain_%s.log" % args.name)

    if not args.worker:
        workroot.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(Path(__file__).resolve()), "--worker"] + sys.argv[1:]
        with open(log_path, "a") as fh:
            p = subprocess.Popen(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL, start_new_session=True)
        print("chain %s started, PID %d; log %s" % (args.name, p.pid, log_path))
        return 0

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    manifest = json.loads(Path(args.runs).read_text())
    out_dir = REPO / manifest["out"]
    runs = manifest["runs"]
    logger.info("chain %s start: %d run(s), %s", args.name, len(runs),
                ", ".join(r["job"] for r in runs))
    logger.info("guards: max-cpus %d, nice %d, max-load %.1f, / >= %.0f GiB, "
                "work >= %.0f GiB", args.max_cpus, args.nice, args.max_load,
                args.min_free_root, args.min_free_work)

    recs = []
    for i, run in enumerate(runs, 1):
        logger.info("--- run %d/%d: %s", i, len(runs), run["job"])
        try:
            rec = run_one(run, args.name, out_dir, workroot, args)
        except (OSError, ValueError, KeyError) as exc:
            logger.exception("run %s raised", run["job"])
            rec = {"job": run["job"], "cpus": 0, "status": "exception: %s" % exc,
                   "seconds": 0, "verdict": "", "detail": ""}
        recs.append(rec)
        logger.info("--- %s: %s %s", rec["job"], rec["status"], rec["verdict"])

    files = [r["file"] for r in recs if r.get("file")]
    if args.push and files:
        _, branch = sh(["git", "rev-parse", "--abbrev-ref", "HEAD"], REPO)
        verdicts = ", ".join("%s %s" % (r["job"], r["verdict"])
                             for r in recs if r["verdict"])
        msg = ("Keio server (fifa): chain %s summaries (%d run(s), automatic)\n\n"
               "scripts/run_chain_keio.py on fifa. %s.\n\n"
               "No analysis yet: what to run next is asked for on the pull "
               "request.\n" % (args.name, len(files), verdicts))
        logger.info("commit/push:\n%s", commit_push(files, msg, branch.strip()))

    if args.pr:
        rc, log = sh(["gh", "pr", "comment", args.pr, "--body-file", "-"], REPO,
                     timeout=120, stdin_text=report(args.name, recs))
        logger.info("pr comment: rc=%d %s", rc, log.strip()[:300])
    logger.info("chain %s end", args.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
