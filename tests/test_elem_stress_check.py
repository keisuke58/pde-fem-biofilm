"""The stage-5 stress checker, on synthetic result rows.

check_elem_stress judges the point-model element from the CSV that
apdl/callsite/post_elem_stress.mac writes and the point-model trace. Here it
is fed rows built to pass, then each property is broken in turn, to make
sure every check can fail. The macro is checked for what this repository
can check without ANSYS: that it writes the columns the reader expects.
"""
import sys
from pathlib import Path

_AU = Path(__file__).resolve().parents[1] / "ansys_usermat"
sys.path.insert(0, str(_AU))
import one_species_reference as ref   # noqa: E402

E = 220


def _rows(tmp_path, alphas=(0.0, 0.1, 0.2, 0.3), p_sign=-1.0, seqv=None,
          nbr=1.0e-3, a_shift=0.0):
    lines = ["lstep,sbstep,time,sx,sy,sz,sxy,syz,sxz,seqv,alpha,"
             "seqv_nbr_max"]
    pm = []
    for i, a in enumerate(alphas, start=1):
        s = p_sign * a * 1.0e-2
        q = seqv[i - 1] if seqv else a * 5.0e-3
        lines.append(f"    1.,{i:8.0f}.,{0.1 * i:24.16E},"
                     + ",".join(f"{x:24.16E}" for x in
                                (s, s, s, 0, 0, 0, q, a + a_shift,
                                 nbr if a > 0 else 0.0)))
        for it in range(3):                   # three iterations, last wins
            pm.append({"elem": E, "ldstep": 1, "isubst": i,
                       "alpha_new": a if it == 2 else -1.0})
    p = tmp_path / "elem_stress.csv"
    p.write_text("\n".join(lines) + "\n")
    return ref.read_elem_stress(p), pm


def test_a_consistent_run_passes(tmp_path):
    rows, pm = _rows(tmp_path)
    c = ref.check_elem_stress(rows, pm, E)
    assert c["n_sets"] == 4
    assert all(c[k] for k in ("alpha_matches", "compressive", "seqv_grows",
                              "loads_neighbours")), c
    assert abs(c["alpha_end"] - 0.3) < 1e-15


def test_alpha_not_from_the_point_model_fails(tmp_path):
    rows, pm = _rows(tmp_path, a_shift=1e-3)
    assert not ref.check_elem_stress(rows, pm, E)["alpha_matches"]


def test_tensile_growth_fails(tmp_path):
    rows, pm = _rows(tmp_path, p_sign=+1.0)
    assert not ref.check_elem_stress(rows, pm, E)["compressive"]


def test_falling_von_mises_fails(tmp_path):
    rows, pm = _rows(tmp_path, seqv=[0.0, 3e-3, 2e-3, 4e-3])
    assert not ref.check_elem_stress(rows, pm, E)["seqv_grows"]


def test_unloaded_neighbours_fail(tmp_path):
    rows, pm = _rows(tmp_path, nbr=0.0)
    assert not ref.check_elem_stress(rows, pm, E)["loads_neighbours"]


def test_the_macro_writes_the_columns_the_reader_expects():
    mac = (_AU / "apdl" / "callsite" / "post_elem_stress.mac").read_text()
    head = mac.split("*VWRITE\n(", 1)[1].split(")", 1)[0].strip("'")
    assert head.split(",") == ["lstep", "sbstep", "time", "sx", "sy", "sz",
                               "sxy", "syz", "sxz", "seqv", "alpha",
                               "seqv_nbr_max"]
    vw = [l for l in mac.splitlines() if l.strip().startswith("*VWRITE,")]
    assert len(vw) == 1 and len(vw[0].split(",")) - 1 == 12
    assert "SVAR,84" in mac
