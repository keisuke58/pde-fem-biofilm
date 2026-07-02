#!/usr/bin/env python3
"""
Unit tests for the Fig4 per-location risk-field module (JAXFEM/risk_field.py).

Run:  python -m pytest tests/test_risk_field.py -v
"""

import sys
from pathlib import Path

import numpy as np
import pytest

_FEM_DIR = Path(__file__).resolve().parent.parent
_JAXFEM = _FEM_DIR / "JAXFEM"
for p in (str(_FEM_DIR), str(_JAXFEM)):
    if p not in sys.path:
        sys.path.insert(0, p)

import risk_field as rf  # noqa: E402

KPA = rf.MPA_PER_KPA


def test_risk_field_shapes_and_values():
    # 3 samples, 4 nodes; node 3 always above 6 kPa, node 0 always below.
    stack = np.array([
        [1.0, 4.0, 7.0, 10.0],
        [2.0, 5.0, 8.0, 11.0],
        [3.0, 6.0, 9.0, 12.0],
    ]) * KPA
    fld = rf.risk_field(stack, tau_mpa=6.0 * KPA)
    assert fld["n_samples"] == 3 and fld["n_nodes"] == 4
    # node 0: none > 6 -> 0 ; node 3: all > 6 -> 1
    assert fld["risk"][0] == 0.0
    assert fld["risk"][3] == 1.0
    # mean is per-node column mean
    assert fld["mean"][0] == pytest.approx(2.0 * KPA)
    # quantiles ordered
    assert np.all(fld["p05"] <= fld["p50"]) and np.all(fld["p50"] <= fld["p95"])


def test_risk_field_rejects_1d():
    with pytest.raises(ValueError):
        rf.risk_field(np.array([1.0, 2.0, 3.0]), 1.0)


def test_auto_pocket_line_orders_by_arclength():
    # points strung along x with small y,z noise -> line should order by x
    x = np.linspace(0, 1, 20)
    coords = np.column_stack([x, 0.01 * np.sin(x), 0.01 * np.cos(x)])
    line = rf.auto_pocket_line(coords)
    assert len(line) >= 2
    xs = coords[line, 0]
    assert np.all(np.diff(xs) >= -1e-9)  # monotonic along principal (x) axis


def test_pocket_band_brackets_mean_and_monotone_arclen():
    stack, coords = rf.synthesize(n_samples=30, nx=8, seed=1)
    band = rf.pocket_line_band(stack, coords)
    assert np.all(np.diff(band["arclen"]) >= -1e-12)     # arc length non-decreasing
    assert np.all(band["p05"] <= band["mean"] + 1e-12)   # band brackets the mean
    assert np.all(band["mean"] <= band["p95"] + 1e-12)


def test_synthesize_deterministic():
    a1, c1 = rf.synthesize(n_samples=10, nx=6, seed=42)
    a2, c2 = rf.synthesize(n_samples=10, nx=6, seed=42)
    assert np.array_equal(a1, a2) and np.array_equal(c1, c2)
    assert a1.shape == (10, 6 ** 3)
    assert c1.shape == (6 ** 3, 3)


def test_risk_peaks_at_pocket_centre():
    # synthetic ridge peaks at y=z=0.5; risk there should exceed the corners.
    stack, coords = rf.synthesize(n_samples=60, nx=12, seed=0)
    fld = rf.risk_field(stack, tau_mpa=5.0 * KPA)
    r2 = (coords[:, 1] - 0.5) ** 2 + (coords[:, 2] - 0.5) ** 2
    centre = fld["risk"][r2 < 0.02].mean()
    corner = fld["risk"][r2 > 0.3].mean()
    assert centre > corner
    assert fld["risk"].max() == pytest.approx(1.0, abs=1e-9)


def test_build_fig4_writes_outputs(tmp_path):
    stack, coords = rf.synthesize(n_samples=20, nx=6, seed=3)
    summary = rf.build_fig4(stack, coords, "unit", 5.0, None, tmp_path, make_plots=False)
    assert (tmp_path / "risk_field_summary_unit.json").exists()
    assert (tmp_path / "risk_field_unit.csv").exists()
    assert summary["n_nodes"] == 6 ** 3
    assert 0.0 <= summary["risk_field"]["mean"] <= 1.0
    # csv has header + one row per node
    rows = (tmp_path / "risk_field_unit.csv").read_text().strip().splitlines()
    assert len(rows) == 1 + 6 ** 3


def test_load_stack_roundtrip(tmp_path):
    stack, coords = rf.synthesize(n_samples=12, nx=5, seed=7)
    np.save(tmp_path / "sigma_stack.npy", stack)
    np.save(tmp_path / "coords.npy", coords)
    s2, c2, line = rf.load_stack(tmp_path)
    assert np.array_equal(s2, stack) and np.array_equal(c2, coords)
    assert line is None


def test_load_stack_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        rf.load_stack(tmp_path)
