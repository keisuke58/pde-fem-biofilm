"""make_slides_3d_demo_pptx.py -- PowerPoint version of slides_3d_demo.tex with the videos embedded (they play inside
the slide; the still of the last frame is the poster). Same content as the beamer deck.

    python make_slides_3d_demo_pptx.py   -> slides_3d_demo.pptx
"""
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.util import Emu, Inches, Pt

ROOT = Path(__file__).resolve().parent
A = ROOT / "assets"
BLUE = RGBColor(30, 80, 160)
FONT = "Times New Roman"

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
W, H = prs.slide_width, prs.slide_height
blank = prs.slide_layouts[6]


def text(slide, x, y, w, h, lines, size=18, color=None, bold=False):
    tb = slide.shapes.add_textbox(x, y, w, h)
    tf = tb.text_frame
    tf.word_wrap = True
    for i, ln in enumerate(lines if isinstance(lines, list) else [lines]):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = ln
        p.space_after = Pt(6)
        for r in p.runs:
            r.font.size = Pt(size); r.font.name = FONT; r.font.bold = bold
            if color is not None:
                r.font.color.rgb = color
    return tb


def frame(title, n):
    s = prs.slides.add_slide(blank)
    bar = s.shapes.add_shape(1, 0, 0, W, Inches(0.85))
    bar.fill.solid(); bar.fill.fore_color.rgb = BLUE; bar.line.fill.background()
    text(s, Inches(0.4), Inches(0.12), Inches(12.5), Inches(0.7), title, size=26, color=RGBColor(255, 255, 255))
    text(s, Inches(11.8), Inches(7.05), Inches(1.4), Inches(0.4), f"{n}", size=12, color=RGBColor(120, 120, 120))
    return s


def movie(s, name, x, y, w):
    png = A / f"{name}.png"
    from PIL import Image
    iw, ih = Image.open(png).size
    h = Emu(int(w * ih / iw))
    s.shapes.add_movie(str(A / f"{name}.mp4"), x, y, w, h, poster_frame_image=str(png), mime_type="video/mp4")
    return h


# 1 title
s = prs.slides.add_slide(blank)
text(s, Inches(1), Inches(2.3), Inches(11.3), Inches(1), "Biofilm Growth and Stress in 3D", size=40, color=BLUE, bold=True)
text(s, Inches(1), Inches(3.3), Inches(11.3), Inches(1),
     "Videos of the growth model: verification problem, seed position, growth law, tooth-like collar", size=20)
text(s, Inches(1), Inches(4.6), Inches(11.3), Inches(1.5),
     ["Keisuke Nishioka", "Institute of Continuum Mechanics (IKM), Leibniz Universität Hannover", "20 October 2026"],
     size=18)

# 2 overview
s = frame("What the videos show", 2)
text(s, Inches(0.6), Inches(1.2), Inches(12.1), Inches(5.8), [
    "• Biofilm phase φ, growth α − 1 with F_g = α I, and the von Mises stress, in 3D over time.",
    "• Growth model of Klempt et al. (2024): dφ/dt = β∇²φ + k_α α (plus the bound on φ), dα/dt = k_α φ; "
    "stiffness E(φ² + f).",
    "• Homeostatic-pressure law (added in this work): growth is scaled by g = min(1, max(0, 1 − p/p_h)), "
    "p = −tr σ / 3.",
    "• Cases: (1) the verification problem of the thesis, seed in the centre of the 2 mm cube; (2) the seed in a "
    "corner, on an edge, or two seeds that meet; (3) the corner seed with and without the growth law; "
    "(4) a tooth-like collar with synthetic four-species data (Fritsch et al. 2025).",
    "• Only case 1 is a thesis result. Cases 2–4 are demonstrations with the same model. No experimental data are used.",
], size=20)

n = 3
# 3 verification
s = frame("1. Verification problem: seed in the centre of the cube", n); n += 1
h = movie(s, "anim_wp3_seed_3d", Inches(0.7), Inches(1.0), Inches(11.9))
text(s, Inches(0.7), Inches(1.05) + h, Inches(11.9), Inches(1.6), [
    "β = 0.02 mm²/T*, k_α = 1e−3 /T*, E = 10 Pa, ν = 0.49, T* = 1.1; one octant cut away. Field on 48³ "
    "(reference solution), stress on 16³.",
    "Seed von Mises at T* = 1.1: 6.19e−4 Pa (thesis, 16³: 6.18e−4 Pa). The stress forms a ring at the edge of the seed."],
    size=16)

# 4 corner / edge
s = frame("2. Seed in a corner and on an edge", n); n += 1
h = movie(s, "anim_seed_corner_3d", Inches(0.3), Inches(1.1), Inches(6.3))
movie(s, "anim_seed_edge_3d", Inches(6.75), Inches(1.1), Inches(6.3))
text(s, Inches(0.5), Inches(1.2) + h, Inches(12.3), Inches(2.5), [
    "Same model and seed volume as case 1; bottom face clamped (the substrate). Left: corner (three faces), right: edge "
    "(two faces).",
    "The biofilm spreads along the walls; the stress concentrates in the biofilm next to the clamped face: the mean over "
    "its bottom quarter millimetre is 7.6 mPa (corner) and 6.4 mPa (edge), 10 to 12 times the seed stress of the centre "
    "case (0.61 mPa); these means change by less than 1 % from 24³ to 48³. The pointwise maximum sits where the clamped "
    "face meets a free wall and grows with refinement (about h^−0.4), so the colour scale stops at the 95th percentile, "
    "about 15 mPa, which does not depend on the mesh."], size=16)

# 5 two seeds (if rendered)
if (A / "anim_two_3d.mp4").exists():
    s = frame("2. Two seeds that grow into each other", n); n += 1
    h = movie(s, "anim_two_3d", Inches(0.7), Inches(1.0), Inches(11.9))
    text(s, Inches(0.7), Inches(1.05) + h, Inches(11.9), Inches(1.6), [
        "Two seeds of 0.5 mm on the clamped bottom face, 0.125 mm apart; front half cut away. Same model as case 1."],
        size=16)

# 6 law
s = frame("3. The growth law removes the concentration", n); n += 1
h = movie(s, "anim_law_3d", Inches(0.7), Inches(1.0), Inches(11.9))
text(s, Inches(0.7), Inches(1.05) + h, Inches(11.9), Inches(1.6), [
    "Corner seed, same colour scale in both stress panels. Without the law the stress in the biofilm next to the clamped "
    "face reaches about 15 mPa (95th percentile); with the law it stays well below. p_h = 3.3 mPa is chosen by me (0.4 times the 95th percentile of the "
    "pressure without the law), not taken from a paper."], size=16)

# 7 collar
s = frame("4. Tooth-like collar with synthetic four-species data", n); n += 1
h = movie(s, "anim_cervical_f4_3d", Inches(0.7), Inches(1.0), Inches(11.9))
text(s, Inches(0.7), Inches(1.05) + h, Inches(11.9), Inches(2), [
    "A 0.25 mm biofilm layer on a crown-shaped collar (radius 2.05 to 3.05 mm, 2 mm high), inner face bonded, nutrient "
    "from the outside. Point model: the four-species case of Fritsch et al. (2025) with their true parameters "
    "(synthetic, no experimental data). One quarter cut away; displacement up to 0.31 mm, drawn at true scale."], size=16)

# 8 section
s = frame("4. Section and composition", n); n += 1
s.shapes.add_picture(str(A / "fig_cervical_f4.png"), Inches(1.4), Inches(1.0), width=Inches(10.5))
text(s, Inches(0.7), Inches(5.6), Inches(11.9), Inches(1.6), [
    "(a) Species 3 takes over the seed (share 0.47 at the end), species 4 almost disappears (0.02).",
    "(b), (c) Section at 45°: without the law the von Mises stress reaches 0.68 E at the bonded base; with p_h/E = 0.05 "
    "the largest value is 0.12 E."], size=16)

# 9 summary
s = frame("Summary", n); n += 1
text(s, Inches(0.6), Inches(1.3), Inches(12.1), Inches(5.6), [
    "• Case 1 reproduces the thesis value of the seed stress (0.2 % difference) and shows it in 3D.",
    "• Moving the seed to a wall breaks the symmetry; the stress then concentrates in the biofilm next to the fixed "
    "surface, 10 to 12 times the seed stress of the centre case (band mean, converged to 1 %; the pointwise maximum is "
    "mesh-dependent and is not quoted).",
    "• The homeostatic-pressure law removes this concentration, in the cube and on the tooth-like collar "
    "(0.68 to 0.12 times E).",
    "• Synthetic data only (Fritsch et al. 2025) keep the demonstration independent of the calibration.",
    "• Next: a real tooth or implant surface, and the calibrated parameters of Chapter 3 at Keio."], size=20)

# speaker notes (EN) from speaker_notes_3d_demo.md, in slide order
import re
secs = re.split(r"^## \d+\. ", (ROOT / "speaker_notes_3d_demo.md").read_text(), flags=re.M)[1:]
notes = []
for sec in secs:
    m = re.search(r"EN: (.*?)\nJA:", sec, flags=re.S)
    notes.append(" ".join(m.group(1).split()) if m else "")
if not (A / "anim_two_3d.mp4").exists():
    del notes[4]
for sl, nt in zip(prs.slides, notes):
    sl.notes_slide.notes_text_frame.text = nt

out = ROOT / "slides_3d_demo.pptx"
prs.save(out)
print("wrote", out, "slides", len(prs.slides), f"{out.stat().st_size / 1e6:.1f} MB")
