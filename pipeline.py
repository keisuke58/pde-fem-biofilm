#!/usr/bin/env python3
"""pipeline.py
=============
T2.1 — single, config-driven entry point for the mechano-ecological pipeline.

    input (condition + geometry)
        -> posterior FEM ensemble        [stage: posterior]
        -> forward PDE alpha -> INP       [stage: forward]
        -> posterior stress CI (sigma)    [stage: stress_ci]
        -> clinical risk  P[sigma>tau]    [stage: risk]

Design
------
This does NOT reimplement the existing stage scripts. It is a thin orchestrator
that (i) reads ONE problem-specific config (a JSON file under ``configs/``),
(ii) separates that config from the generic stage machinery, and (iii) runs the
requested stages, delegating each to the authoritative existing script:

  posterior   -> run_posterior_pipeline.py       (posterior FEM ensemble)
  forward     -> run_end_to_end_pipeline.py       (MAP theta -> PDE -> eigenstrain -> INP -> viz)
  stress_ci   -> JAXFEM/posterior_klempt_stress_ci.py   (posterior sigma CI json)
  risk        -> JAXFEM/risk_metric.py            (P[sigma > threshold]; in-process)

Each stage checks its prerequisites first. Stages whose inputs live in sibling
repositories (``../data_5species/_runs`` TMCMC posteriors) or require Abaqus are
SKIPPED with an explicit reason when those inputs are absent, instead of
crashing -- so the pipeline is runnable in an isolated checkout for the stages
whose inputs are committed (``stress_ci`` reuse + ``risk``).

A manifest of what ran / was skipped is written to the output directory.

Config format (JSON)
--------------------
  {
    "condition": "commensal_hobic",   # CH | CS | DS | DH long name
    "geom": "tooth",                  # tooth | implant
    "stages": ["stress_ci", "risk"],  # subset + order of stages to run
    "risk": {"threshold_kpa": 5.0, "reference_kpa": [2.5, 5.0, 10.0]},
    "posterior": {"n_samples": 20, "nx": 15, "ny": 15, "nz": 15},
    "forward": {"quick": false, "no_inp": true}
  }

Usage
-----
  python pipeline.py --config configs/tooth_commensal_hobic.json
  python pipeline.py --config configs/tooth_commensal_hobic.json --stages risk
  python pipeline.py --list-stages
  python pipeline.py --config configs/tooth_commensal_hobic.json --dry-run
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path

logger = logging.getLogger("pipeline")

_HERE = Path(__file__).resolve().parent
_JAXFEM = _HERE / "JAXFEM"
_DATA_ROOT = _HERE.parent / "data_5species" / "_runs"
_CI_DIR = _JAXFEM / "_posterior_ci"

sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_JAXFEM))

# ── condition -> TMCMC run dir (problem-specific; mirrors the existing scripts) ─
CONDITION_RUNS: dict[str, Path] = {
    "dh_baseline": _DATA_ROOT / "sweep_pg_20260217_081459" / "dh_baseline",
    "dysbiotic_hobic": _DATA_ROOT / "sweep_pg_20260217_081459" / "dh_baseline",
    "commensal_static": _DATA_ROOT / "Commensal_Static_20260208_002100",
    "commensal_hobic": _DATA_ROOT / "Commensal_HOBIC_20260208_002100",
    "dysbiotic_static": _DATA_ROOT / "Dysbiotic_Static_20260207_203752",
}

STAGE_ORDER = ["posterior", "forward", "stress_ci", "risk", "risk_field"]


# ── config ────────────────────────────────────────────────────────────────────
@dataclass
class PipelineConfig:
    condition: str = "commensal_hobic"
    geom: str = "tooth"
    stages: list[str] = field(default_factory=lambda: ["stress_ci", "risk"])
    risk: dict = field(default_factory=dict)
    risk_field: dict = field(default_factory=dict)
    posterior: dict = field(default_factory=dict)
    forward: dict = field(default_factory=dict)
    stress_ci: dict = field(default_factory=dict)
    out_dir: str | None = None

    @staticmethod
    def load(path: Path) -> "PipelineConfig":
        with open(path) as f:
            d = json.load(f)
        known = {f.name for f in PipelineConfig.__dataclass_fields__.values()}  # type: ignore[attr-defined]
        unknown = set(d) - known
        if unknown:
            raise ValueError(f"unknown config keys: {sorted(unknown)} in {path}")
        return PipelineConfig(**d)

    def resolved_out_dir(self) -> Path:
        if self.out_dir:
            return Path(self.out_dir)
        return _HERE / "_pipeline" / f"{self.geom}_{self.condition}"


# ── stage result ────────────────────────────────────────────────────────────────
@dataclass
class StageResult:
    stage: str
    status: str  # "ran" | "skipped" | "failed" | "dry-run"
    detail: str = ""


def _run_cmd(cmd: list[str], dry_run: bool) -> tuple[bool, str]:
    printable = " ".join(str(c) for c in cmd)
    logger.info("$ %s", printable)
    if dry_run:
        return True, f"(dry-run) {printable}"
    ret = subprocess.run([str(c) for c in cmd])
    if ret.returncode != 0:
        return False, f"exit code {ret.returncode}: {printable}"
    return True, printable


# ── stages ────────────────────────────────────────────────────────────────────
def stage_posterior(cfg: PipelineConfig, out_dir: Path, dry_run: bool) -> StageResult:
    run_dir = CONDITION_RUNS.get(cfg.condition)
    if run_dir is None:
        return StageResult("posterior", "skipped", f"no TMCMC run mapping for {cfg.condition!r}")
    if not run_dir.exists():
        return StageResult(
            "posterior", "skipped",
            f"TMCMC run dir absent (sibling repo): {run_dir}",
        )
    p = cfg.posterior
    cmd = [
        sys.executable, _HERE / "run_posterior_pipeline.py",
        "--conditions", cfg.condition,
        "--n-samples", p.get("n_samples", 20),
        "--nx", p.get("nx", 15), "--ny", p.get("ny", 15), "--nz", p.get("nz", 15),
    ]
    ok, detail = _run_cmd(cmd, dry_run)
    return StageResult("posterior", "dry-run" if dry_run else ("ran" if ok else "failed"), detail)


def stage_forward(cfg: PipelineConfig, out_dir: Path, dry_run: bool) -> StageResult:
    run_dir = CONDITION_RUNS.get(cfg.condition)
    theta = (run_dir / "theta_MAP.json") if run_dir else None
    if theta is None or not theta.exists():
        return StageResult(
            "forward", "skipped",
            f"theta_MAP.json absent (sibling repo): {theta}",
        )
    f = cfg.forward
    cmd = [sys.executable, _HERE / "run_end_to_end_pipeline.py", "--condition", cfg.condition]
    if f.get("quick"):
        cmd.append("--quick")
    if f.get("no_inp", True):
        cmd.append("--no-inp")
    ok, detail = _run_cmd(cmd, dry_run)
    return StageResult("forward", "dry-run" if dry_run else ("ran" if ok else "failed"), detail)


_ALL_CONDITIONS = ["commensal_hobic", "dysbiotic_hobic", "commensal_static", "dysbiotic_static"]


def stage_stress_ci(cfg: PipelineConfig, out_dir: Path, dry_run: bool) -> StageResult:
    """Produce the posterior stress CI json.

    Reuse-by-default: the committed CI json already carries the full per-condition
    posterior samples. Recomputing requires the upstream 0D samples for ALL four
    conditions (``_ci_0d_results/<cond>/samples_0d.json``); a partial set would
    silently collapse the missing conditions to MAP-only and clobber the
    committed json. Recompute therefore only runs when the config sets
    ``stress_ci.recompute=true`` AND every condition's 0D samples are present.
    """
    ci_json = _CI_DIR / f"klempt_stress_ci_{cfg.geom}.json"
    recompute = bool(cfg.stress_ci.get("recompute", False))

    if not recompute:
        if ci_json.exists():
            return StageResult("stress_ci", "skipped", f"reusing committed {ci_json.name} (recompute=false)")
        return StageResult("stress_ci", "skipped", f"no committed CI json ({ci_json.name}); set stress_ci.recompute")

    ci0d = _HERE / "_ci_0d_results"
    missing = [c for c in _ALL_CONDITIONS if not (ci0d / c / "samples_0d.json").exists()]
    if missing:
        return StageResult(
            "stress_ci", "skipped",
            f"recompute requested but 0D samples missing for {missing} — "
            f"would clobber committed json; reusing {ci_json.name}",
        )
    cmd = [sys.executable, _JAXFEM / "posterior_klempt_stress_ci.py", "--geom", cfg.geom]
    ok, detail = _run_cmd(cmd, dry_run)
    return StageResult("stress_ci", "dry-run" if dry_run else ("ran" if ok else "failed"), detail)


def stage_risk(cfg: PipelineConfig, out_dir: Path, dry_run: bool) -> StageResult:
    """Clinical risk metric — runs in-process from the (committed) CI json."""
    ci_json = _CI_DIR / f"klempt_stress_ci_{cfg.geom}.json"
    if not ci_json.exists():
        return StageResult(
            "risk", "skipped",
            f"CI json absent ({ci_json}); run stress_ci first",
        )
    if dry_run:
        return StageResult("risk", "dry-run", f"risk_metric on {ci_json.name}")

    import risk_metric as rm

    r = cfg.risk
    ci = rm.load_ci(cfg.geom)
    summary = rm.compute_risk(
        ci,
        threshold_kpa=r.get("threshold_kpa", 5.0),
        reference_kpa=r.get("reference_kpa"),
        seed=r.get("seed", 0),
    )
    summary["geom"] = cfg.geom
    summary["condition_of_interest"] = cfg.condition
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / f"risk_summary_{cfg.geom}.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    rm._print_table(summary, cfg.geom)
    if not r.get("no_plot", False):
        try:
            rm.plot_risk(ci, summary, cfg.geom, out_dir)
        except Exception as exc:
            logger.warning("risk plot skipped: %s", exc)
    return StageResult("risk", "ran", f"summary -> {summary_path}")


def stage_risk_field(cfg: PipelineConfig, out_dir: Path, dry_run: bool) -> StageResult:
    """Fig4 per-location risk field — needs a posterior stress stack (ensemble output)."""
    rfcfg = cfg.risk_field
    stack_dir = rfcfg.get("stack_dir")
    if not stack_dir:
        return StageResult("risk_field", "skipped", "no risk_field.stack_dir in config (needs FEM ensemble)")
    stack_dir = Path(stack_dir)
    if not (stack_dir / "sigma_stack.npy").exists():
        return StageResult("risk_field", "skipped", f"stack absent: {stack_dir}/sigma_stack.npy")
    if dry_run:
        return StageResult("risk_field", "dry-run", f"risk_field on {stack_dir}")

    import risk_field as rfmod

    sigma, coords, line = rfmod.load_stack(stack_dir)
    summary = rfmod.build_fig4(
        sigma, coords,
        tag=rfcfg.get("tag", f"{cfg.geom}_{cfg.condition}"),
        threshold_kpa=rfcfg.get("threshold_kpa", cfg.risk.get("threshold_kpa", 5.0)),
        line_nodes=line,
        out_dir=out_dir,
        make_plots=not rfcfg.get("no_plot", False),
    )
    rfmod._print_summary(summary)
    return StageResult("risk_field", "ran", f"summary -> {out_dir}/risk_field_summary_{summary['tag']}.json")


STAGES = {
    "posterior": stage_posterior,
    "forward": stage_forward,
    "stress_ci": stage_stress_ci,
    "risk": stage_risk,
    "risk_field": stage_risk_field,
}


# ── driver ──────────────────────────────────────────────────────────────────────
def run(cfg: PipelineConfig, stages: list[str], dry_run: bool) -> list[StageResult]:
    out_dir = cfg.resolved_out_dir()
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 64)
    logger.info("mechano-ecological pipeline")
    logger.info("  condition : %s", cfg.condition)
    logger.info("  geometry  : %s", cfg.geom)
    logger.info("  stages    : %s", stages)
    logger.info("  out_dir   : %s", out_dir)
    logger.info("=" * 64)

    results: list[StageResult] = []
    for name in stages:
        if name not in STAGES:
            results.append(StageResult(name, "failed", "unknown stage"))
            continue
        logger.info("── stage: %s ──", name)
        res = STAGES[name](cfg, out_dir, dry_run)
        logger.info("   [%s] %s", res.status.upper(), res.detail)
        results.append(res)

    manifest = {
        "config": asdict(cfg),
        "stages_requested": stages,
        "results": [asdict(r) for r in results],
    }
    with open(out_dir / "pipeline_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)

    logger.info("=" * 64)
    for r in results:
        logger.info("  %-10s %s", r.stage, r.status)
    logger.info("  manifest -> %s", out_dir / "pipeline_manifest.json")
    logger.info("=" * 64)
    return results


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    ap = argparse.ArgumentParser(
        description="Config-driven mechano-ecological pipeline (T2.1)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--config", help="Path to a JSON config (see configs/)")
    ap.add_argument("--stages", nargs="+", default=None, help="Override the stages to run")
    ap.add_argument("--condition", default=None, help="Override the condition")
    ap.add_argument("--geom", default=None, choices=["tooth", "implant"], help="Override geometry")
    ap.add_argument("--dry-run", action="store_true", help="Print the plan without running stages")
    ap.add_argument("--list-stages", action="store_true", help="List available stages and exit")
    args = ap.parse_args(argv)

    if args.list_stages:
        print("Available stages (canonical order):")
        for s in STAGE_ORDER:
            print(f"  {s:<10} {STAGES[s].__doc__ or ''}".rstrip())
        return 0

    if not args.config and not (args.condition and args.geom):
        ap.error("provide --config, or both --condition and --geom")

    if args.config:
        cfg = PipelineConfig.load(Path(args.config))
    else:
        cfg = PipelineConfig()

    if args.condition:
        cfg.condition = args.condition
    if args.geom:
        cfg.geom = args.geom
    stages = args.stages or cfg.stages

    results = run(cfg, stages, args.dry_run)
    # non-zero exit only on a genuine failure (skips are fine)
    return 1 if any(r.status == "failed" for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
