"""The trace checker that decides whether the first real run succeeded.
Synthetic traces only: each test builds the trace a correct or a broken call
site would write, and asserts the checker tells them apart."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ansys_usermat"))
import one_species_reference as ref   # noqa: E402

K = 0.5


def _trace(per_iteration_bug=False, wrong_rate=False, n_sub=4, iters=3):
    rows, header = [], "elem,ip,ldstep,isubst,dtime,bio1,locbio1,alpha_n,alpha_new"
    for elem, ip, phi0 in ((1, 1, 0.2), (7, 1, 0.6)):
        alpha = 0.0
        for s in range(1, n_sub + 1):
            phi = phi0 + 0.05 * s                       # grows in time
            dt = 0.1
            a_n = alpha
            for it in range(iters):
                rate = 2 * K if wrong_rate else K
                a_new = a_n + rate * phi * dt
                rows.append(f"{elem},{ip},1,{s},{dt},{phi},{phi/2},{a_n},{a_new}")
                if per_iteration_bug:
                    a_n = a_new                         # counts every iteration
            alpha = a_new
    return [header] + rows


def _write(tmp_path, lines):
    p = tmp_path / "phi_trace.csv"
    p.write_text("\n".join(lines) + "\n")
    return ref.read_trace(p)


def test_a_correct_call_site_passes_every_check(tmp_path):
    c = ref.check_trace(_write(tmp_path, _trace()), K)
    assert c["once_per_increment"] and c["eq36"] and c["carried"]


def test_counting_growth_per_iteration_is_caught(tmp_path):
    c = ref.check_trace(_write(tmp_path, _trace(per_iteration_bug=True)), K)
    assert not c["once_per_increment"]


def test_a_wrong_k_alpha_is_caught(tmp_path):
    c = ref.check_trace(_write(tmp_path, _trace(wrong_rate=True)), K)
    assert not c["eq36"]


def test_phi_statistics_describe_each_candidate(tmp_path):
    c = ref.check_trace(_write(tmp_path, _trace()), K)
    s = c["phi_bio1"]
    assert s["varies_between_points"] and s["changes_in_time"] and s["in_0_1"]


def test_k_alpha_for_a_target_growth():
    k = ref.k_alpha_for_target(0.02, phi_max=1.0, time_total=1.1)
    assert abs(k * 1.0 * 1.1 - 0.02) < 1e-15
