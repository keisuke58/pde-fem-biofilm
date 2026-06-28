from odbAccess import openOdb
o=openOdb("kv1.odb")
step=o.steps[o.steps.keys()[-1]]
fr=step.frames[-1]
S=fr.fieldOutputs["S"]
v=S.values[0]
print("FRAME time", fr.frameValue)
print("S11 S22 S33 S12 S13 S23 =", [round(x,8) for x in v.data])
print("MISES =", round(v.mises,8))
o.close()
