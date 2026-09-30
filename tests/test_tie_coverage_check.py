"""Tests for tier2b_real/tie_coverage_check.py.

The tool exists to say how much of a *TIE interface is actually tied, after
the cylinder decks turned out to be bonded at four corner nodes. A diagnostic
that is itself wrong is worse than none: over-report the gap and it sends you
hunting a bug that is not there, under-report it and it misses the one that
is. So the distance it rests on is pinned against hand-computable answers in
all three Voronoi regions of a triangle -- face, edge and vertex -- not just
the easy interior case.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
_SRC = _ROOT / "tier2b_real" / "tie_coverage_check.py"

_spec = importlib.util.spec_from_file_location("tie_coverage_check", _SRC)
tcc = importlib.util.module_from_spec(_spec)
sys.modules["tie_coverage_check"] = tcc
_spec.loader.exec_module(tcc)

# The unit triangle in z = 0: A=(0,0,0), B=(1,0,0), C=(0,1,0).
TRI = np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0]]], dtype=float)


@pytest.mark.parametrize("p,expected,region", [
    ((0.2, 0.2, 0.3), 0.3,            "above the face interior"),
    ((0.0, 0.0, 0.0), 0.0,            "exactly on vertex A"),
    ((2.0, 0.0, 0.0), 1.0,            "past vertex B"),
    ((-1.0, -1.0, 0.0), np.sqrt(2.0), "past vertex A"),
    ((0.5, -1.0, 0.0), 1.0,           "past edge AB"),
    ((1.0, 1.0, 0.0), np.sqrt(0.5),   "past edge BC"),
])
def test_point_triangle_distance_is_exact(p, expected, region):
    got = tcc._tri_dist(np.asarray(p, dtype=float), TRI)[0]
    assert got == pytest.approx(expected, abs=1e-12), region


def test_vertex_distance_would_not_do():
    """Why the exact routine is needed at all: a node 0.1 above the middle of
    a large face is 0.1 from the surface and much further from any vertex, so
    a vertex-distance check would call a sound interface unbonded."""
    big = np.array([[[0, 0, 0], [10, 0, 0], [0, 10, 0]]], dtype=float)
    p = np.array([3.0, 3.0, 0.1])
    assert tcc._tri_dist(p, big)[0] == pytest.approx(0.1, abs=1e-12)
    assert np.sqrt(((big[0] - p) ** 2).sum(axis=1)).min() > 4.0


def _write_inp(tmp_path):
    """Master face in z=0; slave nodes at 0.10, 0.50 and 2.00 above it."""
    inp = tmp_path / "tie.inp"
    inp.write_text(
        "*NODE\n"
        "1, 0.0, 0.0, 0.0\n2, 1.0, 0.0, 0.0\n3, 0.0, 1.0, 0.0\n"
        "4, 0.0, 0.0, -1.0\n"
        "11, 0.10, 0.10, 0.10\n12, 0.20, 0.20, 0.50\n"
        "13, 0.10, 0.20, 2.00\n14, 0.10, 0.10, 3.00\n"
        "*ELEMENT, TYPE=C3D4, ELSET=MASTER\n101, 1, 2, 3, 4\n"
        "*ELEMENT, TYPE=C3D4, ELSET=SLAVE\n201, 11, 12, 13, 14\n"
        "*SURFACE, NAME=M_SURF, TYPE=ELEMENT\n 101, S1\n"
        "*SURFACE, NAME=S_SURF, TYPE=ELEMENT\n 201, S1\n"
        "*TIE, NAME=T_TEST, ADJUST=NO, POSITION TOLERANCE=1.0\n S_SURF, M_SURF\n"
    )
    return inp


def test_parses_nodes_elements_surfaces_and_the_tie_tolerance(tmp_path):
    nodes, elems, etype, surfaces, ties = tcc.parse_inp(str(_write_inp(tmp_path)))
    assert len(nodes) == 8 and len(elems) == 2
    assert etype[101] == "C3D4"
    assert surfaces["M_SURF"] == [(101, 1)] and surfaces["S_SURF"] == [(201, 1)]
    assert ties == [["T_TEST", "S_SURF", "M_SURF", 1.0]]


def test_gaps_are_the_heights_the_nodes_were_placed_at(tmp_path):
    """Face S1 of a C3D4 is local nodes 1-2-3, so node 14 is not on the slave
    surface and must not appear; the other three sit directly above the master
    triangle, so their gaps are exactly their z."""
    nodes, elems, etype, surfaces, _ = tcc.parse_inp(str(_write_inp(tmp_path)))
    sn = tcc.face_nodes(surfaces["S_SURF"], elems, etype, corners_only=False)
    assert sn == [11, 12, 13]
    P = np.array([nodes[i] for i in sn], dtype=float)
    T = tcc.face_triangles(surfaces["M_SURF"], elems, etype, nodes)
    d = np.sort(tcc.point_tri_dist(P, T))
    np.testing.assert_allclose(d, [0.10, 0.50, 2.00], atol=1e-12)


def test_exit_code_flags_an_interface_that_is_not_fully_tied(tmp_path, capsys):
    """One of the three slave nodes sits 2.0 from a surface tied at 1.0, so the
    tool must say so and exit non-zero -- that is the whole point of it."""
    rc = tcc.main([str(_write_inp(tmp_path))])
    out = capsys.readouterr().out
    assert "beyond tolerance: 1 of 3" in out
    assert rc == 1
