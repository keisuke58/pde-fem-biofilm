"""trace_subset_equal.py -- regression check across the comp_trace threshold change.

    python ansys_usermat/apdl/trace_subset_equal.py OLD.csv NEW.csv [--stride 37] [--phi-min 0.5]

Since 5 Oct the composition mode traces every integration point 1 with
phi >= 0.01 (CM_TRACE_PHI_MIN) instead of phi >= 0.5. A run with the new
executable is the same as before when the rows of its trace that the old rule
would have written (MOD(elem, stride) = 1 or phi_used >= phi_min) equal the old
trace line for line, as text. Exit 0 if equal, 1 otherwise.
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--stride", type=int, default=37)
    ap.add_argument("--phi-min", type=float, default=0.5)
    a = ap.parse_args()
    old = [l.rstrip("\r\n") for l in open(a.old) if l.strip()]
    new = []
    for l in open(a.new):
        if not l.strip():
            continue
        v = l.split(",")
        if int(v[0]) % a.stride == 1 or float(v[8]) >= a.phi_min:
            new.append(l.rstrip("\r\n"))
    same = old == new
    print(f"old {len(old)} rows, new {len(new)} rows after the old rule: {'identical' if same else 'DIFFERENT'}")
    if not same:
        for i, (x, y) in enumerate(zip(old, new)):
            if x != y:
                print(f"first difference at row {i + 1}:\n  old {x[:120]}\n  new {y[:120]}")
                break
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()
