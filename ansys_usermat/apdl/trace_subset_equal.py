"""trace_subset_equal.py -- regression check across the comp_trace threshold change.

    python ansys_usermat/apdl/trace_subset_equal.py OLD.csv NEW.csv [--stride 37] [--phi-min 0.5]

Since 5 Oct the composition mode traces every integration point 1 with
phi >= 0.01 (CM_TRACE_PHI_MIN) instead of phi >= 0.5. A run with the new
executable is the same as before when the rows of its trace that the old rule
would have written (MOD(elem, stride) = 1 or phi_used >= phi_min) equal the old
trace line for line: integers exactly, reals to --rtol (default 1e-10). Not as
text: on 5 Oct the n-species loops changed the compiled arithmetic of the
rescale, and the rerun differed by up to 3e-12 (relative) with identical
stresses; a text comparison called that a failure and the old exe was restored.
Exit 0 if equal, 1 otherwise.
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--stride", type=int, default=37)
    ap.add_argument("--phi-min", type=float, default=0.5)
    ap.add_argument("--rtol", type=float, default=1e-10)
    a = ap.parse_args()
    old = [l.rstrip("\r\n") for l in open(a.old) if l.strip()]
    new = []
    for l in open(a.new):
        if not l.strip():
            continue
        v = l.split(",")
        if int(v[0]) % a.stride == 1 or float(v[8]) >= a.phi_min:
            new.append(l.rstrip("\r\n"))
    def close(x, y):
        a, b = x.split(","), y.split(",")
        if len(a) != len(b) or a[:6] != b[:6]:
            return False
        return all(abs(float(p) - float(q)) <= a_rtol * max(abs(float(p)), abs(float(q)), 1e-300)
                   for p, q in zip(a[6:], b[6:]))
    a_rtol = a.rtol
    bad = [i for i, (x, y) in enumerate(zip(old, new)) if x != y and not close(x, y)]
    same = len(old) == len(new) and not bad
    print(f"old {len(old)} rows, new {len(new)} rows after the old rule: "
          f"{'same to rtol ' + str(a.rtol) if same else 'DIFFERENT'}")
    if bad:
        i = bad[0]
        print(f"first difference at row {i + 1}:\n  old {old[i][:120]}\n  new {new[i][:120]}")
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()
