// English deck: the rigid obstacles and the maze of Klempt et al. 2024 with the coupled model (2 Nov 2026).
// Built by build.sh geom.
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const path = require('path');
const R = path.resolve(__dirname, '../..') + '/';
const B = __dirname + '/build/', EQ = B + 'eq/';
const vid = (c) => R + 'assets/anim_coupled_tvd_' + c + '_ab_inhliv_transc_alpha_5panel.mp4';   // panel (d) = alpha - 1 (2 Nov)
const poster = (c) => R + 'assets/anim_coupled_tvd_' + c + '_ab_inhliv_transc_alpha_5panel_end.png';

const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka';
pres.title = 'Rigid obstacles and a maze in the coupled growth model';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', RULE = '8A99AD', ACC = 'A61E1E', LIGHT = 'F3F5F8';
const TM = 'Times New Roman';
let n = 0;

function frame(title) {
  const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addText(title, { x: 0.55, y: 0.32, w: 12.2, h: 0.6, fontFace: TM, fontSize: 26, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 0.98, w: 12.23, h: 0, line: { color: NAVY, width: 1.25 } });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 7.0, w: 12.23, h: 0, line: { color: RULE, width: 0.5 } });
  s.addText('K. Nishioka  |  Rigid obstacles and a maze in the coupled growth model', { x: 0.55, y: 7.05, w: 10, h: 0.3, fontFace: TM, fontSize: 11, color: GREY, isTextBox: true, margin: 0 });
  s.addText(String(n), { x: 11.78, y: 7.05, w: 1.0, h: 0.3, fontFace: TM, fontSize: 12, color: GREY, align: 'right', isTextBox: true, margin: 0 });
  return s;
}
const cap = (s, t, y, x = 0.55, w = 12.23) => s.addText(t, { x, y, w, h: 0.55, fontFace: TM, fontSize: 13, color: INK, isTextBox: true, margin: 0, valign: 'top' });
const bullets = (s, items, o) => s.addText(items.map((t, i) => ({ text: t, options: { bullet: { indent: 18 }, breakLine: i < items.length - 1 } })),
  Object.assign({ fontFace: TM, fontSize: 18, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 10 }, o));
function pngSize(p) { const b = fs.readFileSync(p); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
const img = (s, p, x, y, w) => { const [pw, ph] = pngSize(p); s.addImage({ path: p, x, y, w, h: w * ph / pw }); return w * ph / pw; };
const H = (t) => ({ text: t, options: { bold: true, color: 'FFFFFF', fill: { color: NAVY }, align: 'center' } });
const C = (t, o = {}) => ({ text: t, options: Object.assign({ align: 'center' }, o) });
const L = (t, o = {}) => C(t, Object.assign({ align: 'left' }, o));


const A = R + 'assets/';
const obstVid = A + 'anim_coupled_tvd_separate_ab_inhliv_obst_transc_alpha_5panel.mp4';
const obstPoster = obstVid.replace('.mp4', '_end.png');
const vidSlide = (s, v, poster, y, w, hmax = 9) => { const [pw, ph] = pngSize(poster); if (w * ph / pw > hmax) w = hmax * pw / ph; const hh = w * ph / pw; const x = (13.33 - w) / 2;
  if (v && fs.existsSync(v)) s.addMedia({ type: 'video', path: v, cover: 'data:image/png;base64,' + fs.readFileSync(poster).toString('base64'), x, y, w, h: hh });
  else s.addImage({ path: poster, x, y, w, h: hh }); return hh; };

const bcTable = (s, rows, y) => { const c_ = (t, o = {}) => ({ text: t, options: Object.assign({ align: 'left', valign: 'middle' }, o) });
  s.addTable([[H('field'), H('where it is held / fixed'), H('value'), H('from'), H('elsewhere')]].concat(rows.map((r) => [c_(r[0], { bold: true }), c_(r[1]), c_(r[2]), c_(r[3]), c_(r[4])])),
    { x: 0.55, y, w: 12.23, colW: [1.9, 5.0, 1.0, 1.9, 2.43], fontFace: TM, fontSize: 13, color: INK, rowH: 0.36, border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } }); };

// 1 title
{ const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 2.0, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText('Rigid obstacles and a maze in the coupled growth model', { x: 0.9, y: 2.2, w: 11.5, h: 0.7, fontFace: TM, fontSize: 32, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addText('The geometries of Klempt et al. 2024, Sections 4.3 and 4.4, with two species and an antibiotic', { x: 0.9, y: 3.0, w: 11.5, h: 0.6, fontFace: TM, fontSize: 20, color: INK, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 3.8, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText([{ text: 'Keisuke Nishioka', options: { fontSize: 20, bold: true, breakLine: true } },
    { text: 'Institute of Continuum Mechanics, Leibniz Universität Hannover', options: { fontSize: 16, breakLine: true } },
    { text: 'November 2026', options: { fontSize: 16 } }],
    { x: 0.9, y: 4.25, w: 11.5, h: 1.5, fontFace: TM, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 6 }); }

// 2 what is repeated
{ const s = frame('1. What is repeated');
  bullets(s, [
    'Klempt et al. 2024 grow the biofilm against rigid obstacles (Section 4.3) and through a maze (Section 4.4). Both are repeated here with the coupled model: the point model of Klempt et al. 2026 at every point, nutrient and antibiotic taken up by the living biomass, g φ₃D (Σᵢφ̄ᵢ) c (Eq. 35).',
    'The equations, the point model (Fritsch et al. 2025, case I, two species) and the constants are those of the cube examples; only the geometry changes.',
    'Obstacles: four rigid bars in the 2 mm cube, nutrient held at y = −1 mm, the stress is computed (frictionless contact with the bars).',
    'Maze: the maze of Fig. 11 of the paper, in two dimensions, constants of its Table 6; no stress (2D). The question here is the composition: what happens to the living biomass after the antibiotic, and which rule for newly reached cells this needs.'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 }); }

// 3 obstacles: set-up and growth
{ const s = frame('2. Rigid obstacles: set-up and growth');
  const h = img(s, B + 'geom_obst3d.png', 3.2, 1.15, 6.9);
  cap(s, 'Fig. 1  Four rigid bars (dark grey) after Fig. 9 of Klempt et al. 2024, positions read off the sketch: two along x and two along z, 0.25 mm square, 0.5 to 0.75 mm from the centre. '
    + 'They hold no biofilm and no nutrient; the biofilm slides along them without friction, and they replace the three corner supports of the cube examples. '
    + '(a) Set-up: start biomass (red), nutrient supply on y = −1 mm (green). (b)-(d) Surface φ₃D = 0.5 (separate start) and 0.2 (mixed start) at T* = 0.05, 0.1 and 0.2. Python, 32³, no antibiotic.', 1.15 + h + 0.1); }

// 3b obstacles: boundary conditions
{ const s = frame('2. Rigid obstacles (cont.): initial and boundary conditions');
  bcTable(s, [
    ['biofilm φ₃D', 'start: species A in the seed (32 elements of 0.25 mm), species B in four blocks along the edges parallel to z (separate); A and B together in the seed (mixed)', 'φ₃D = 1', 'T* = 0', 'φ₃D = 0'],
    ['nutrient c', 'element layer on the face y = −1 mm (0.25 mm thick)', 'c = 1', 'T* = 0', 'zero flux (walls and bar faces)'],
    ['antibiotic a', 'outer element shell on all six faces (0.25 mm thick)', 'a = 1', 'T* = 1/3 (video only)', 'zero flux at the bar faces'],
    ['bars', 'four rigid bars, 0.25 mm square: two along x and two along z, 0.5 to 0.75 mm from the centre', 'no biofilm, no c, a', '—', '—'],
    ['displacement u', 'frictionless unilateral contact with the bars (the biofilm slides along them)', 'u·n ≤ 0', '—', 'free; no corner supports']], 1.3);
  cap(s, 'The bars replace the three corner supports of the cube examples: they are the only mechanical constraint, so the swelling pushes the biofilm against them. '
    + 'The numbers of the next slide are without antibiotic (32³); the video there is the separate start with the antibiotic (64³).', 4.6); }

// 4 obstacles: video and numbers
{ const s = frame('3. Rigid obstacles: time evolution and stress');
  const hh = vidSlide(s, obstVid, obstPoster, 1.15, 11.6, 3.7);
  const y0 = 1.15 + hh + 0.12;
  s.addTable([[H('seed mean at T* = 1 (32³)'), H('without bars'), H('with bars'), H('ratio')],
    [L('separate: von Mises [Pa]'), C('7.4 × 10⁻⁴'), C('7.1 × 10⁻³'), C('≈ 10', { bold: true })],
    [L('mixed: von Mises [Pa]'), C('2.1 × 10⁻⁴'), C('9.5 × 10⁻⁴'), C('≈ 4.5')],
    [L('growth α − 1 (separate / mixed)'), C('—'), C('9 % / 21 % lower'), C('')]],
    { x: 2.2, y: y0, w: 8.9, colW: [3.5, 1.8, 1.8, 1.8], fontFace: TM, fontSize: 13, color: INK, rowH: 0.3, border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  cap(s, 'The bars restrain the swelling, so the seed stress rises about tenfold; the composition of the seed hardly changes. Values on one grid (32³), where the stress without bars is 14 % above its extrapolated value. Click the frame to play.', y0 + 1.35); }

// 4b maze: boundary conditions
{ const s = frame('4. Maze: initial and boundary conditions');
  bcTable(s, [
    ['biofilm φ₃D', 'colony of both species at the entrance in the outer wall (top), shares of Fritsch et al. (0.25, 0.30)', 'φ₃D = 1', 'T* = 0', 'φ₃D = 0'],
    ['nutrient c', 'central disk of the maze', 'c = 1', 'T* = 0', 'zero flux through all walls'],
    ['antibiotic a', 'the entrance (where the colony starts)', 'a = 1', 'T* = 20 (one third of T* = 60)', 'zero flux through all walls'],
    ['walls', 'maze walls read off Fig. 11 of Klempt et al. 2024, 10 µm wide', 'no biofilm, no c, a', '—', '—'],
    ['displacement u', 'not computed (two-dimensional model)', '—', '—', '—']], 1.3);
  cap(s, 'Two-dimensional, 4 µm cells, 560 µm across. Constants of Table 6 of Klempt et al. 2024 for the maze (d = 10¹³ µm²/T*, β = 10 µm²/T*, k_α = 10⁻⁶ T*⁻¹, r = 200 µm/T*), '
    + 'the rest as Table 2 (g = 10⁸ T*⁻¹, k = 1); point model Fritsch et al. 2025, case I, α*₀ = 100, φ_min = 10⁻³. Newly reached cells take a state by one of two rules (slide 6).', 4.6); }

// 5 maze: growth
{ const s = frame('5. Maze: growth through the corridors');
  const h = img(s, B + 'geom_maze.png', 2.4, 1.15, 8.5);
  cap(s, 'Fig. 2  The maze of Fig. 11 of Klempt et al. 2024, read off the figure (560 µm across), 2D with 4 µm cells, constants of their Table 6 (β = 10 µm²/T*, r = 200 µm/T*, k_α = 10⁻⁶ T*⁻¹). '
    + 'Nutrient held in the centre, the colony of both species at the entrance. φ₃D at T* = 12, 24, 36 and 60, (a) without and (b) with the antibiotic held at the entrance from T* = 20. '
    + 'At T* = 60 the nearest point with φ₃D ≥ 0.5 is about 40 µm from the centre with and without the antibiotic (39 to 43 µm on 4, 2 and 1 µm cells; 133 µm only with antibiotic on 4 µm cells, where one cell decides it).', 1.15 + h + 0.1); }

// 6 maze: the rule for new cells
{ const s = frame('6. Maze with antibiotic: the rule for newly reached cells decides');
  const hh = vidSlide(s, A + 'anim_maze_rules.mp4', A + 'anim_maze_rules_end.png', 1.15, 10.2, 4.6);
  cap(s, 'Living volume fraction Σᵢφ̄ᵢ over time when newly reached cells (left) take the state of their colonised neighbours and (right) start from the colony\'s start state. '
    + 'After the antibiotic at T* = 20 it falls to 0.10; it recovers to 0.35 at T* = 60 with the first rule, because new cells start with a large volume fraction and the antibiotic term scales as 1/φᵢ², and stays at 0.04 with the second. Click to play.', 1.15 + hh + 0.12); }

// 7 summary
{ const s = frame('7. Summary');
  bullets(s, [
    'The coupled model runs in the geometries of Klempt et al. 2024 without any change of the equations or constants.',
    'Rigid bars restrain the swelling: the seed stress rises about tenfold (separate start), the composition of the seed hardly changes.',
    'In the maze the biofilm keeps spreading after the antibiotic, because the front speed r c/(k+c) depends on the nutrient only; dead biofilm spreads too.',
    'Whether the living biomass recovers in the maze depends on the state given to newly reached cells (0.35 against 0.04 at T* = 60). The point model has no transport, so this rule is an assumption of the coupling. In this model the biofilm spreads by the front term, not by material flow (the growth-driven velocity is about 10⁻³ of the front speed), so even a species-resolved continuum would need such a rule at the edge.'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 }); }

// 8 references
{ const s = frame('References');
  bullets(s, [
    'F. Klempt, M. Soleimani, P. Wriggers, P. Junker: A Hamilton principle-based model for diffusion-driven biofilm growth. Biomech. Model. Mechanobiol. 23, 2091–2113 (2024).',
    'F. Klempt, H. Geisler, M. Soleimani, P. Junker: A continuum multi-species bacterial growth model with a novel interaction scheme. Arch. Appl. Mech. 96, 164 (2026).',
    'L. Fritsch, H. Geisler, J. Grashorn, F. Klempt, M. Soleimani, M. Broggi, P. Junker, M. Beer: Bayesian updating of constitutive parameters under hybrid uncertainties with a novel surrogate model applied to biofilms. arXiv:2512.15145 (2025).'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4, fontSize: 14 }); }

pres.writeFile({ fileName: B + 'geometries_en.pptx' }).then(() => console.log('ok ' + B + 'geometries_en.pptx ' + n + ' slides'));
