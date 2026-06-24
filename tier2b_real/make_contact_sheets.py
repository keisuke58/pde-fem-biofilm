#!/usr/bin/env python3
"""Per-part contact sheet: top row Mises x6 angles, bottom row DSPSS x6 angles."""
import os
from PIL import Image, ImageDraw, ImageFont

SRC = "/home/nishioka/IKM_Hiwi/FEM/tier2b_real/png_out"
OUT = os.path.join(SRC, "contact")
os.makedirs(OUT, exist_ok=True)

PARTS = ["tier2b_crown", "tier2b_crown_E10000", "tier2b_crown_E210000",
         "tier2b_crownreal", "tier2b_generic", "tier2b_real",
         "tm_healthy", "tm_moderate", "tm_deep", "pimp_e4"]
ANGLES = ["iso", "iso2", "iso3", "front", "top", "right"]
FIELDS = ["mises", "dspss"]

TW = 560          # thumbnail width
PAD = 10          # gap between cells
LBLW = 90         # left label column width
HEAD = 46         # header strip height
SUBH = 26         # angle-label strip under header

def font(sz):
    for p in ("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
              "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
              "/usr/share/fonts/gnu-free/FreeSansBold.ttf"):
        if os.path.exists(p):
            return ImageFont.truetype(p, sz)
    return ImageFont.load_default()

F_HEAD, F_LBL = font(30), font(20)

def thumb(part, field, angle):
    p = os.path.join(SRC, "%s_%s_%s.png" % (part, field, angle))
    if not os.path.exists(p):
        return None
    im = Image.open(p).convert("RGB")
    bbox = Image.eval(im, lambda x: 255 - x).getbbox()  # trim white border
    if bbox:
        im = im.crop(bbox)
    th = int(im.height * TW / im.width)
    return im.resize((TW, th), Image.LANCZOS)

for part in PARTS:
    cells = {(f, a): thumb(part, f, a) for f in FIELDS for a in ANGLES}
    have = [c for c in cells.values() if c]
    if not have:
        print("skip (no images):", part); continue
    rowh = max(c.height for c in have)
    cols = len(ANGLES)
    grid_w = LBLW + cols * TW + (cols + 1) * PAD
    grid_h = HEAD + SUBH + len(FIELDS) * (rowh + PAD) + PAD
    sheet = Image.new("RGB", (grid_w, grid_h), "white")
    d = ImageDraw.Draw(sheet)
    d.text((PAD, 8), part, fill="black", font=F_HEAD)
    # angle column headers
    for j, a in enumerate(ANGLES):
        x = LBLW + PAD + j * (TW + PAD)
        d.text((x + 4, HEAD), a, fill="#444444", font=F_LBL)
    # rows
    for i, f in enumerate(FIELDS):
        y = HEAD + SUBH + i * (rowh + PAD) + PAD
        d.text((6, y + rowh // 2 - 10), f.upper(), fill="black", font=F_LBL)
        for j, a in enumerate(ANGLES):
            c = cells[(f, a)]
            if c is None:
                continue
            x = LBLW + PAD + j * (TW + PAD)
            sheet.paste(c, (x, y))
    outp = os.path.join(OUT, "%s_contact.png" % part)
    sheet.save(outp)
    print("OK", outp, sheet.size)
print("DONE ->", OUT)
