from odbAccess import openOdb
import sys
o=openOdb(sys.argv[1])
fr=o.steps[o.steps.keys()[-1]].frames[-1]
v=fr.fieldOutputs["S"].values[0]
print("  S11,S22,S33 = %.4e %.4e %.4e  Mises=%.4e" % (v.data[0],v.data[1],v.data[2],v.mises))
o.close()
