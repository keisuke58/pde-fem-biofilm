from odbAccess import openOdb
o=openOdb("free_exp.odb")
fr=o.steps[o.steps.keys()[-1]].frames[-1]
U=fr.fieldOutputs["U"]
for v in U.values:
    if v.nodeLabel in (2,7):
        print("  node",v.nodeLabel,"U=",[round(x,6) for x in v.data])
o.close()
