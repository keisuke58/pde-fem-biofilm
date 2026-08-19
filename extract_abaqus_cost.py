#!/usr/bin/env python3
"""Extract computational-cost evidence from Abaqus job files.

The thesis needs a "computational cost" section backed by measurement, not
recollection. Abaqus already writes everything required -- increment/iteration
history in ``.sta``, timing and model size in ``.msg``/``.dat`` -- but the repo
only ever kept the stress results, so the cost numbers were lost with the
scratch directories.

This script reads those files and writes a committable artifact. It needs no
Abaqus licence and no ODB access: the inputs are plain text, so it runs
anywhere, including in CI against archived job files.

Usage
-----
    # every job in a directory (matched by <jobname>.sta)
    python3 extract_abaqus_cost.py /path/to/abaqus/scratch -o runs/abaqus_cost.json

    # specific jobs, plus a Markdown table for pasting into the thesis
    python3 extract_abaqus_cost.py p23_klempt_A_commensal_hobic --markdown

    # print to stdout without writing anything
    python3 extract_abaqus_cost.py . --stdout

What is extracted
-----------------
``.sta``  increments, total/severe-discontinuity/equilibrium iterations,
          attempts, cutbacks (an attempt number > 1 means the increment was
          retried at a smaller step), final step time.
``.msg``  wallclock and CPU seconds (``JOB TIME SUMMARY``), plus any
          convergence warnings worth reporting alongside a timing.
``.dat``  model size: number of elements, nodes and total DOF.

Derived: iterations per increment, and seconds per equilibrium iteration --
the figure to quote when comparing UMAT cost against a reference material.

Missing files are not an error. A job with only a ``.sta`` yields iteration
counts with ``wallclock_s: null``; the report says so rather than guessing.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

# ── .msg "JOB TIME SUMMARY" block ────────────────────────────────────────────
# Abaqus prints e.g.  "  WALLCLOCK TIME (SEC) =        456"
_TIME_PATTERNS = {
    "user_cpu_s": re.compile(r"USER\s+TIME\s*\(SEC\)\s*=\s*([0-9.eE+-]+)", re.I),
    "system_cpu_s": re.compile(r"SYSTEM\s+TIME\s*\(SEC\)\s*=\s*([0-9.eE+-]+)", re.I),
    "total_cpu_s": re.compile(r"TOTAL\s+CPU\s+TIME\s*\(SEC\)\s*=\s*([0-9.eE+-]+)", re.I),
    "wallclock_s": re.compile(r"WALLCLOCK\s+TIME\s*\(SEC\)\s*=\s*([0-9.eE+-]+)", re.I),
}

# ── .dat model-size block ────────────────────────────────────────────────────
_SIZE_PATTERNS = {
    "n_elements": re.compile(r"NUMBER\s+OF\s+ELEMENTS\s+IS\s+([0-9]+)", re.I),
    "n_nodes": re.compile(r"NUMBER\s+OF\s+NODES\s+IS\s+([0-9]+)", re.I),
    "n_dof": re.compile(
        r"TOTAL\s+NUMBER\s+OF\s+(?:VARIABLES|DEGREES\s+OF\s+FREEDOM)\s+IN\s+THE\s+MODEL\s*[:=]?\s*([0-9]+)",
        re.I,
    ),
}


def _first_float(pattern: re.Pattern, text: str):
    m = pattern.search(text)
    return float(m.group(1)) if m else None


def _first_int(pattern: re.Pattern, text: str):
    m = pattern.search(text)
    return int(m.group(1)) if m else None


def parse_sta(path: Path) -> dict:
    """Parse an Abaqus ``.sta`` increment history.

    Column layout (Abaqus/Standard):
        STEP  INC  ATT  SEVERE-DISCON-ITERS  EQUIL-ITERS  TOTAL-ITERS
        STEP-TIME/LPF  INC-OF-TIME/LPF ...

    Header and banner lines are skipped by requiring the first six fields to be
    integers; that is the only reliable discriminator across Abaqus versions,
    which vary the header wording and the trailing columns.
    """
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 6:
            continue
        try:
            step, inc, att, disc, equil, total = (int(f) for f in fields[:6])
        except ValueError:
            continue
        row = {
            "step": step,
            "inc": inc,
            "attempt": att,
            "severe_discon_iters": disc,
            "equil_iters": equil,
            "total_iters": total,
        }
        # trailing float columns are version-dependent; take them if parseable
        try:
            row["step_time"] = float(fields[6])
        except (IndexError, ValueError):
            row["step_time"] = None
        rows.append(row)

    if not rows:
        return {"increments": 0, "rows": []}

    # An increment is *converged* when it is the last attempt recorded for that
    # (step, inc) pair. Earlier attempts are cutbacks.
    last_attempt: dict[tuple[int, int], dict] = {}
    for r in rows:
        key = (r["step"], r["inc"])
        if key not in last_attempt or r["attempt"] >= last_attempt[key]["attempt"]:
            last_attempt[key] = r

    converged = list(last_attempt.values())
    cutbacks = sum(1 for r in rows if r["attempt"] > 1)
    step_times = [r["step_time"] for r in converged if r["step_time"] is not None]

    return {
        "increments": len(converged),
        "attempts": len(rows),
        "cutbacks": cutbacks,
        "total_iters": sum(r["total_iters"] for r in rows),
        "equil_iters": sum(r["equil_iters"] for r in rows),
        "severe_discon_iters": sum(r["severe_discon_iters"] for r in rows),
        "steps": sorted({r["step"] for r in rows}),
        "final_step_time": max(step_times) if step_times else None,
        "rows": rows,
    }


def parse_msg(path: Path) -> dict:
    text = path.read_text(errors="replace")
    out = {k: _first_float(p, text) for k, p in _TIME_PATTERNS.items()}
    warn = len(re.findall(r"^\s*\*\*\*WARNING", text, re.M))
    err = len(re.findall(r"^\s*\*\*\*ERROR", text, re.M))
    out["warnings"] = warn
    out["errors"] = err
    return out


def parse_dat(path: Path) -> dict:
    text = path.read_text(errors="replace")
    return {k: _first_int(p, text) for k, p in _SIZE_PATTERNS.items()}


def collect_job(stem: Path) -> dict:
    """Gather every available file for one job, keyed by its path stem."""
    job: dict = {"job": stem.name, "files": {}}
    # Seed every field so the JSON schema is identical across jobs: a value that
    # could not be read is an explicit null, never an absent key. Downstream
    # consumers (and the thesis table) can then distinguish "not measured" from
    # "zero" without guessing.
    job.update(
        {k: None for k in
         ("increments", "attempts", "cutbacks", "total_iters", "equil_iters",
          "severe_discon_iters", "final_step_time", "warnings", "errors")}
    )
    job["steps"] = []
    job.update({k: None for k in _TIME_PATTERNS})
    job.update({k: None for k in _SIZE_PATTERNS})

    for ext, parser in ((".sta", parse_sta), (".msg", parse_msg), (".dat", parse_dat)):
        p = stem.with_suffix(ext)
        if p.exists():
            job["files"][ext] = str(p)
            parsed = parser(p)
            rows = parsed.pop("rows", None)
            if rows is not None:
                job["_sta_rows"] = rows
            job.update(parsed)

    # derived quantities -- only where the inputs actually exist
    incs, iters = job.get("increments"), job.get("total_iters")
    if incs:
        job["iters_per_increment"] = round(iters / incs, 2) if iters else None
    wall, iters = job.get("wallclock_s"), job.get("total_iters")
    job["s_per_iteration"] = round(wall / iters, 3) if wall and iters else None
    wall, nel = job.get("wallclock_s"), job.get("n_elements")
    job["us_per_element_iteration"] = (
        round(1e6 * wall / (nel * job["total_iters"]), 2)
        if wall and nel and job.get("total_iters")
        else None
    )
    return job


def discover(target: Path) -> list[Path]:
    """Return job stems: every ``*.sta`` in a directory, or one explicit job."""
    if target.is_dir():
        return sorted({p.with_suffix("") for p in target.glob("*.sta")})
    # a bare job name, or a path with any of the known extensions
    stem = target.with_suffix("") if target.suffix in (".sta", ".msg", ".dat", ".odb") else target
    return [stem]


_COLUMNS = [
    ("job", "Job"),
    ("n_elements", "Elements"),
    ("n_dof", "DOF"),
    ("increments", "Increments"),
    ("cutbacks", "Cutbacks"),
    ("total_iters", "Iterations"),
    ("iters_per_increment", "Iters/inc"),
    ("wallclock_s", "Wallclock (s)"),
    ("total_cpu_s", "CPU (s)"),
    ("s_per_iteration", "s/iter"),
]


def to_markdown(jobs: list[dict]) -> str:
    head = "| " + " | ".join(label for _, label in _COLUMNS) + " |"
    rule = "|" + "|".join("---" for _ in _COLUMNS) + "|"
    lines = [head, rule]
    for j in jobs:
        cells = []
        for key, _ in _COLUMNS:
            v = j.get(key)
            cells.append("—" if v is None else (f"`{v}`" if key == "job" else str(v)))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Extract Abaqus cost evidence from .sta / .msg / .dat files.",
    )
    ap.add_argument(
        "target",
        type=Path,
        help="directory to scan for *.sta, or a single job name / file path",
    )
    ap.add_argument("-o", "--output", type=Path, help="write JSON here")
    ap.add_argument("--markdown", action="store_true", help="also print a Markdown table")
    ap.add_argument("--stdout", action="store_true", help="print JSON instead of writing it")
    ap.add_argument(
        "--keep-rows",
        action="store_true",
        help="keep the per-increment .sta rows in the JSON (verbose)",
    )
    args = ap.parse_args(argv)

    stems = discover(args.target)
    jobs = [collect_job(s) for s in stems]
    jobs = [j for j in jobs if j["files"]]

    if not jobs:
        print(f"no Abaqus job files (.sta/.msg/.dat) found under {args.target}", file=sys.stderr)
        return 1

    if not args.keep_rows:
        for j in jobs:
            j.pop("_sta_rows", None)

    missing = [j["job"] for j in jobs if j.get("wallclock_s") is None]
    if missing:
        print(
            f"note: no timing found for {len(missing)} job(s) "
            f"(no .msg, or no JOB TIME SUMMARY): {', '.join(missing)}",
            file=sys.stderr,
        )

    payload = {"jobs": jobs}
    text = json.dumps(payload, indent=2)

    if args.stdout or not args.output:
        print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")
        print(f"wrote {args.output}  ({len(jobs)} job(s))", file=sys.stderr)
    if args.markdown:
        print("\n" + to_markdown(jobs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
