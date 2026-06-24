# render_odb_png.py - headless multi-angle contour PNGs for a set of .odb files.
# Config via environment variables (abaqus eats '--', so we avoid argv):
#   ODB_OUT   : output directory
#   ODB_LIST  : newline- or space-separated odb paths
#   ODB_FIELD : 'mises' | 'dspss' | 'both'   (default 'both')
# Run: ODB_OUT=... ODB_LIST=... xvfb-run -a abq2024 viewer noGUI=render_odb_png.py
from abaqus import *
from abaqusConstants import *
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
import os, sys, shutil, tempfile

def out(msg):
    sys.__stdout__.write(msg + "\n"); sys.__stdout__.flush()

OUT   = os.environ['ODB_OUT']
LIST  = os.environ['ODB_LIST'].split()
FIELD = os.environ.get('ODB_FIELD', 'both').lower()
if not os.path.isdir(OUT):
    os.makedirs(OUT)

vp = session.viewports['Viewport: 1']
vp.makeCurrent(); vp.maximize()
session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='#FFFFFF')
session.pngOptions.setValues(imageSize=(2400, 1600))
session.printOptions.setValues(vpDecorations=OFF, vpBackground=OFF, reduceColors=False)

# named views + a couple of rotated isometrics -> "various nice angles"
def set_iso():   vp.view.setValues(session.views['Iso'])
def set_front(): vp.view.setValues(session.views['Front'])
def set_top():   vp.view.setValues(session.views['Top'])
def set_right(): vp.view.setValues(session.views['Right'])
def set_iso2():  vp.view.setValues(session.views['Iso']); vp.view.rotate(xAngle=0, yAngle=130, zAngle=0, mode=MODEL)
def set_iso3():  vp.view.setValues(session.views['Iso']); vp.view.rotate(xAngle=25, yAngle=-110, zAngle=0, mode=MODEL)
VIEWS = [('iso', set_iso), ('iso2', set_iso2), ('iso3', set_iso3),
         ('front', set_front), ('top', set_top), ('right', set_right)]

def add_dspss(odb):
    """Write a scalar field DSPSS = sigma1+sigma2+sigma3 = tr(sigma) into each frame."""
    for sn in odb.steps.keys():
        for fr in odb.steps[sn].frames:
            if 'S' not in fr.fieldOutputs or 'DSPSS' in fr.fieldOutputs:
                continue
            S = fr.fieldOutputs['S']
            sp = (S.getScalarField(invariant=MAX_PRINCIPAL)
                  + S.getScalarField(invariant=MID_PRINCIPAL)
                  + S.getScalarField(invariant=MIN_PRINCIPAL))
            nf = fr.FieldOutput(name='DSPSS', description='sigma1+sigma2+sigma3 (tr sigma)', type=SCALAR)
            nf.addData(field=sp)

def last_valid(odb):
    """(stepIndex, frameIndex, stepName) of the last non-empty frame, or None."""
    keys = odb.steps.keys()
    for si in range(len(keys) - 1, -1, -1):
        if len(odb.steps[keys[si]].frames) > 0:
            return si, len(odb.steps[keys[si]].frames) - 1, keys[si]
    return None

def disp_positions(odb, var, fr):
    """Candidate display positions ordered by preference, based on stored position."""
    p = fr.fieldOutputs[var].locations[0].position
    if p == CENTROID:
        return [ELEMENT_CENTROID]
    if p == INTEGRATION_POINT:
        return [ELEMENT_NODAL, INTEGRATION_POINT]   # ELEMENT_NODAL = smooth averaged contour
    return [ELEMENT_NODAL, p]

def render(odb, tag, var, kind, comp, name):
    lv = last_valid(odb)
    if lv is None:
        raise RuntimeError('no frames with results')
    si, fi, sname = lv
    fr = odb.steps[sname].frames[fi]
    if var not in fr.fieldOutputs:
        raise RuntimeError('%s not in last valid frame' % var)
    vp.setValues(displayedObject=odb)
    vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF,))
    # auto-pick a valid display position (centroid vs integration-point odbs differ)
    ok = False
    for pos in disp_positions(odb, var, fr):
        try:
            if kind == INVARIANT:
                vp.odbDisplay.setPrimaryVariable(variableLabel=var, outputPosition=pos,
                                                 refinement=(INVARIANT, comp))
            else:
                vp.odbDisplay.setPrimaryVariable(variableLabel=var, outputPosition=pos)
            ok = True; break
        except Exception:
            continue
    if not ok:
        raise RuntimeError('no valid display position for %s' % var)
    vp.odbDisplay.setFrame(step=si, frame=fi)
    vp.viewportAnnotationOptions.setValues(triad=OFF, title=ON, state=OFF, compass=OFF)
    for vname, setter in VIEWS:
        setter(); vp.view.fitView()
        session.printToFile(os.path.join(OUT, '%s_%s_%s' % (name, tag, vname)), PNG, (vp,))
    out('  OK %s [%s] x%d views' % (name, tag, len(VIEWS)))

tmpdir = tempfile.mkdtemp(prefix='odbrender_', dir=os.environ.get('ODB_TMP') or None)
for path in LIST:
    name = os.path.splitext(os.path.basename(path))[0]
    try:
        if FIELD in ('dspss', 'both'):
            # operate on a COPY so originals stay pristine: write DSPSS, save, close,
            # then re-open read-only for display (viewer cannot display a writable odb).
            cp = os.path.join(tmpdir, name + '.odb')
            shutil.copy(path, cp)
            odbw = session.openOdb(cp, readOnly=False)
            add_dspss(odbw); odbw.save(); odbw.close()
            odb = session.openOdb(cp, readOnly=True)
            render(odb, 'dspss', 'DSPSS', SCALAR, None, name)
            if FIELD == 'both':
                render(odb, 'mises', 'S', INVARIANT, 'Mises', name)
            odb.close(); os.remove(cp)
        else:
            odb = session.openOdb(path, readOnly=True)
            render(odb, 'mises', 'S', INVARIANT, 'Mises', name)
            odb.close()
    except Exception as e:
        out('  FAIL %s : %s' % (name, e))
shutil.rmtree(tmpdir, ignore_errors=True)
out('DONE -> %s' % OUT)
