"""extract_history.py -- dump the full stress-vs-time history from an odb.

Run with Abaqus's own Python (odbAccess is only importable there), not a
regular Python interpreter:

    abaqus python extract_history.py <path-to-odb> [element_index]

Prints one CSV-formatted line per output frame: step_time,S11,S22,S33,Mises.
"""
import sys

from odbAccess import openOdb

odb_path = sys.argv[1]
elem_idx = int(sys.argv[2]) if len(sys.argv) > 2 else 0

o = openOdb(odb_path)
step = o.steps[o.steps.keys()[-1]]

print("step_time,S11,S22,S33,Mises")
for fr in step.frames:
    v = fr.fieldOutputs["S"].values[elem_idx]
    print("%.6f,%.6e,%.6e,%.6e,%.6e" % (fr.frameValue, v.data[0], v.data[1], v.data[2], v.mises))

o.close()
