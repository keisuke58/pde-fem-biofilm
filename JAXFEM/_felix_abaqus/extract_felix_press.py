"""Abaqus Python ODB extraction for felix_klempt snapshots.
Run: abaqus python extract_felix_press.py <jobname>
Writes: <jobname>_press.npy  -- shape (N_elem, 4) = [x, y, PRESS, alpha]
"""
from __future__ import print_function
import sys, os
import numpy as np
from odbAccess import openOdb

job = sys.argv[1] if len(sys.argv) > 1 else 'felix_klempt_t0'
odb = openOdb(job + '.odb')
step = list(odb.steps.values())[-1]
fr   = step.frames[-1]

S     = fr.fieldOutputs['S']
COORD = fr.fieldOutputs['COORD']

try:
    SDV3 = fr.fieldOutputs['SDV3']
    sdv3 = {v.elementLabel: v.data for v in SDV3.values}
except KeyError:
    sdv3 = {}

coord_map = {v.elementLabel: v.data for v in COORD.values}

rows = []
for v in S.values:
    el   = v.elementLabel
    sd   = v.data          # (S11, S22, S33, S12, S13, S23) or (S11,S22,S33,S12)
    if len(sd) >= 3:
        press = -(sd[0] + sd[1] + sd[2]) / 3.0
    else:
        press = 0.0
    xy = coord_map.get(el, [0.0, 0.0, 0.0])
    alpha = sdv3.get(el, 0.0)
    if hasattr(alpha, '__len__'):
        alpha = alpha[0] if len(alpha) > 0 else 0.0
    rows.append([xy[0], xy[1], press, alpha])

odb.close()
arr = np.array(rows, dtype=float)
out = job + '_press.npy'
np.save(out, arr)
print('Saved', out, arr.shape)
