#!/usr/bin/env python3
"""Per-part contact sheets (top=Mises x6 angles, bottom=DSPSS x6) for any png_out dir.
Usage: python3 make_contact_sheets_generic.py <png_dir> [png_dir2 ...]
Auto-discovers parts from filenames '<part>_<field>_<angle>.png'."""
import os, re, sys
from PIL import Image, ImageDraw, ImageFont

ANGLES = ["iso", "iso2", "iso3", "front", "top", "right"]
FIELDS = ["mises", "dspss"]
PAT = re.compile(r"^(?P<part>.+)_(?P<field>mises|dspss)_(?P<angle>iso|iso2|iso3|front|top|right)\.png$")
TW, PAD, LBLW, HEAD, SUBH = 560, 10, 90, 46, 26

def font(sz):
    for p in ("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
              "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
              "/usr/share/fonts/gnu-free/FreeSansBold.ttf"):
        if os.path.exists(p):
            return ImageFont.truetype(p, sz)
    return ImageFont.load_default()
F_HEAD, F_LBL = font(30), font(20)

def thumb(src, part, field, angle):
    p = os.path.join(src, "%s_%s_%s.png" % (part, field, angle))
    if not os.path.exists(p):
        return None
    im = Image.open(p).convert("RGB")
    bb = Image.eval(im, lambda x: 255 - x).getbbox()
    if bb:
        im = im.crop(bb)
    return im.resize((TW, int(im.height * TW / im.width)), Image.LANCZOS)

def build(src):
    parts = set()
    for fn in os.listdir(src):
        mt = PAT.match(fn)
        if mt:
            parts.add(mt.group("part"))
    if not parts:
        print("  no part PNGs in", src); return
    out = os.path.join(src, "contact"); os.makedirs(out, exist_ok=True)
    for part in sorted(parts):
        cells = {(f, a): thumb(src, part, f, a) for f in FIELDS for a in ANGLES}
        rows = [f for f in FIELDS if any(cells[(f, a)] for a in ANGLES)]  # only fields present
        have = [c for c in cells.values() if c]
        if not have or not rows:
            continue
        rowh = max(c.height for c in have); cols = len(ANGLES)
        W = LBLW + cols * TW + (cols + 1) * PAD
        H = HEAD + SUBH + len(rows) * (rowh + PAD) + PAD
        sheet = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(sheet)
        d.text((PAD, 8), part, fill="black", font=F_HEAD)
        for j, a in enumerate(ANGLES):
            d.text((LBLW + PAD + j * (TW + PAD) + 4, HEAD), a, fill="#444444", font=F_LBL)
        for i, f in enumerate(rows):
            y = HEAD + SUBH + i * (rowh + PAD) + PAD
            d.text((6, y + rowh // 2 - 10), f.upper(), fill="black", font=F_LBL)
            for j, a in enumerate(ANGLES):
                c = cells[(f, a)]
                if c:
                    sheet.paste(c, (LBLW + PAD + j * (TW + PAD), y))
        outp = os.path.join(out, "%s_contact.png" % part)
        sheet.save(outp); print("  OK", outp)

for src in sys.argv[1:]:
    print("DIR", src); build(src)
print("DONE")
