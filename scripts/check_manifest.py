"""Check a manifest against the IKMHIWI03 reference headers.

Each results_1006/<job>.txt records, on its second line, the make_*_inp.py
command the reference was produced with. A run that is meant to reproduce it
must pass the same flags, so compare the two argument lists directly instead
of trusting that they were copied correctly by hand.
"""
import json
import sys
from pathlib import Path

ROOT = Path("/home/nishioka/IKM_Hiwi/pde-fem-biofilm")
manifest = json.loads((ROOT / sys.argv[1]).read_text())

print("chain:", manifest["chain"], "| out:", manifest["out"])
ok = True
for run in manifest["runs"]:
    job = run["job"]
    ref = ROOT / run["reference"]
    if not ref.is_file():
        print("  %-14s MISSING REFERENCE %s" % (job, run["reference"]))
        ok = False
        continue
    header = [l for l in ref.read_text().splitlines() if l.startswith("#")]
    cmd = next((l for l in header if "make_klempt_inp.py" in l), None)
    if cmd is None:
        print("  %-14s reference has no make_ command line" % job)
        ok = False
        continue
    # the header drops the output path; the manifest carries it as {inp}
    want = cmd.lstrip("# ").split()[1:]
    got = [a for a in run["make"][1:] if a != "{inp}"]
    if want == got:
        print("  %-14s flags match the reference (%d args)" % (job, len(got)))
    else:
        ok = False
        print("  %-14s FLAGS DIFFER" % job)
        print("      reference: %s" % " ".join(want))
        print("      manifest : %s" % " ".join(got))
    case_ref = "fig7_low" if "fig7_low" in cmd else ("fig7_high" if "fig7_high" in cmd else "fig4_corner")
    if run["summarise"][-1] != case_ref:
        ok = False
        print("      summarise --case is %s, reference is %s"
              % (run["summarise"][-1], case_ref))

print("\nall flags match the references" if ok else "\nMISMATCHES ABOVE")
sys.exit(0 if ok else 1)
