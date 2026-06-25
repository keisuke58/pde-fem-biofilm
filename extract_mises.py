"""abaqus python extract_mises.py <odb>  -> prints 'MISES_MAX <value MPa>'"""
import sys
from odbAccess import openOdb
o = openOdb(sys.argv[1])
fr = o.steps[o.steps.keys()[-1]].frames[-1]
mx = 0.0
for v in fr.fieldOutputs["S"].values:
    if v.mises > mx:
        mx = v.mises
print("MISES_MAX %.8e" % mx)
o.close()
