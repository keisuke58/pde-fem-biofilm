// English deck, results only: the four-species antibiotic case of Klempt et al. 2026 (Table 2, case 4) in the coupled cube
// examples (2 Nov 2026; to send in December). Built by build.sh k26c4.
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const path = require('path');
const R = path.resolve(__dirname, '../..') + '/';
const B = __dirname + '/build/', EQ = B + 'eq/';
const vid = (c) => R + 'assets/anim_coupled_tvd_' + c + '_ab_inhliv_transc_alpha_5panel.mp4';   // panel (d) = alpha - 1 (2 Nov)
const poster = (c) => R + 'assets/anim_coupled_tvd_' + c + '_ab_inhliv_transc_alpha_5panel_end.png';

const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka';
pres.title = 'Four species with an antibiotic in the coupled growth model';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', RULE = '8A99AD', ACC = 'A61E1E', LIGHT = 'F3F5F8';
const TM = 'Times New Roman';
let n = 0;

function frame(title) {
  const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addText(title, { x: 0.55, y: 0.32, w: 12.2, h: 0.6, fontFace: TM, fontSize: 26, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 0.98, w: 12.23, h: 0, line: { color: NAVY, width: 1.25 } });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 7.0, w: 12.23, h: 0, line: { color: RULE, width: 0.5 } });
  s.addText('K. Nishioka  |  Four species with an antibiotic in the coupled growth model', { x: 0.55, y: 7.05, w: 10, h: 0.3, fontFace: TM, fontSize: 11, color: GREY, isTextBox: true, margin: 0 });
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

const kvid = (c) => A + 'anim_coupled_tvd_' + c + '_ab_inhliv_k26c4_transc_alpha_5panel.mp4';

// 1 title
{ const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 2.0, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText('Four species with an antibiotic in the coupled growth model', { x: 0.9, y: 2.2, w: 11.5, h: 0.7, fontFace: TM, fontSize: 32, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addText('Case 4 of Klempt et al. 2026 (Table 2) in the cube examples: results', { x: 0.9, y: 3.0, w: 11.5, h: 0.6, fontFace: TM, fontSize: 20, color: INK, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 3.8, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText([{ text: 'Keisuke Nishioka', options: { fontSize: 20, bold: true, breakLine: true } },
    { text: 'Institute of Continuum Mechanics, Leibniz Universität Hannover', options: { fontSize: 16, breakLine: true } },
    { text: 'December 2026', options: { fontSize: 16 } }],
    { x: 0.9, y: 4.25, w: 11.5, h: 1.5, fontFace: TM, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 6 }); }

// 2 set-up
{ const s = frame('1. Set-up');
  bullets(s, [
    'The coupled model and the 2 mm cube of the two-species examples (Klempt et al. 2024, Eqs. 34–36, with the point model of Klempt et al. 2026 at every point); only the point-model case changes.',
    'Point model: Table 2, case 4 of Klempt et al. 2026, four species, used as published. The antibiotic acts on species i through the term −α* bᵢ ψᵢ, with sensitivities b = (10, 2, 1, 0.01): species 4 is almost insensitive.',
    'Antibiotic a = 1 on the outer shell from T* = 1/3 (α*₀ = 100), as in the two-species examples; c*₀ = 100. Time link s = 0.15 (T* = 1 spans the 1500 steps of the case).',
    'Separate start: species 1 and 2 in the seed, species 3 and 4 in the four blocks along the edges. Mixed start: all four species in equal amounts (0.02 each) in the seed.',
    '64³ finite volumes for the fields, 16³ hexahedra for the stress.'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 }); }

// 2b what changes from two to four species
{ const s = frame('2. From two to four species: what changes');
  const r = (a, b, c) => [L(a, { bold: true }), L(b), L(c)];
  s.addTable([[H(''), H('two-species examples'), H('four-species case')],
    r('point-model case', 'Fritsch et al. 2025, case I', 'Klempt et al. 2026, Table 2, case 4'),
    r('interaction matrix A', '[[1, 0.1], [0.1, 1]]', '[[1, 5, 5, 5], [5, 1, 3, 3], [5, 3, 1, 2], [5, 3, 2, 1]]'),
    r('viscosities η', '(1, 2)', '(0.8, 1, 1.5, 2)'),
    r('antibiotic sensitivities b', '(1, 2)', '(10, 2, 1, 0.01)'),
    r('start amounts φᵢ (point model)', '(0.25, 0.30)', '0.02 each'),
    r('time link s', '0.10 (1000 steps)', '0.15 (1500 steps)'),
    r('inner steps per outer step 0.0025', '25', '38'),
    r('separate start', 'A in the seed, B in the edge blocks', 'species 1, 2 in the seed, 3, 4 in the edge blocks'),
    r('mixed start', 'A and B in equal amounts (0.275 each)', 'all four in equal amounts (0.02 each)'),
    r('run time, 64³ (quarter domain)', '170–350 s', '350–470 s (memory about the same)')],
    { x: 0.55, y: 1.2, w: 12.23, colW: [3.6, 3.7, 4.93], fontFace: TM, fontSize: 13, color: INK, rowH: 0.36, border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  bullets(s, [
    'Unchanged: the field equations (Klempt et al. 2024, Eqs. 34–36), the constants, the grid, the boundary conditions, the antibiotic from T* = 1/3 (α*₀ = 100) and c*₀ = 100.',
    'The uptake g φ₃D (Σᵢφ̄ᵢ) c keeps its form; the sum runs over four species instead of two.'],
    { x: 0.75, y: 5.45, w: 11.8, h: 1.45, fontSize: 15 }); }

// 3, 4 videos
for (const [c, lab] of [['separate', 'separate start'], ['mixed', 'mixed start']]) {
  const s = frame((c === 'separate' ? '3. ' : '4. ') + 'Time evolution, ' + lab);
  const p = kvid(c).replace('.mp4', '_end.png');
  const hh = vidSlide(s, kvid(c), p, 1.15, 12.2, 4.7);
  cap(s, '(a) Dominant species (colour of the species with the largest share, lighter where it leads by less), (b) living volume fraction Σᵢφ̄ᵢ, (c) nutrient c, (d) growth α − 1, (e) von Mises stress. '
    + 'One octant removed; the antibiotic from T* = 1/3 is written on the frames. Click to play.', 1.15 + hh + 0.12); }

// 5 shares
{ const s = frame('5. Share of each species');
  const w = 5.35, y = 1.08, x0 = (13.33 - 2 * w - 0.5) / 2, x1 = x0 + w + 0.5;
  const h1 = img(s, A + 'fig_four_species_shares_separate_ab.png', x0, y, w);
  img(s, A + 'fig_four_species_shares_mixed_ab.png', x1, y, w);
  s.addText('(a) separate start', { x: x0, y: y + h1 + 0.02, w, h: 0.3, fontFace: TM, fontSize: 13, bold: true, color: INK, align: 'center', isTextBox: true, margin: 0 });
  s.addText('(b) mixed start', { x: x1, y: y + h1 + 0.02, w, h: 0.3, fontFace: TM, fontSize: 13, bold: true, color: INK, align: 'center', isTextBox: true, margin: 0 });
  cap(s, 'Share φᵢ/Σⱼφⱼ of each species where φ₃D > 0.5, just before the antibiotic (T* = 0.325) and at T* = 1, and its mean over the biofilm. All four species coexist before the antibiotic; afterwards species 4 takes over (mean share 0.78 and 0.83 at T* = 1).', y + h1 + 0.36); }

// 6 against the point model alone
{ const s = frame('6. Against the point model alone');
  const h = img(s, A + 'fig_four_species_living_vs_0d.png', 1.65, 1.1, 10.0);
  cap(s, 'φ̄ᵢ = φᵢψᵢ of each species: the point model alone (c* = 100 throughout) and the coupled cube (mean over the biofilm, solid; over the seed, dashed). Dotted line: antibiotic at T* = 1/3.', 1.1 + h + 0.05);
  s.addTable([[H(''), H('point model alone'), H('cube, separate start'), H('cube, mixed start')],
    [L('φ̄₄ at T* = 1 (biofilm mean)', { bold: true }), C('0.93'), C('0.75'), C('0.79')],
    [L('nutrient c in the seed at T* = 1', { bold: true }), C('1 (c* = 100)'), C('0.60'), C('0.65')],
    [L('seed: α − 1  /  von Mises at T* = 1', { bold: true }), C('—'), C('9.1 × 10⁻⁴  /  6.5 × 10⁻⁴ Pa'), C('8.4 × 10⁻⁴  /  3.5 × 10⁻⁴ Pa')]],
    { x: 0.55, y: 1.1 + h + 0.6, w: 12.23, colW: [3.9, 2.5, 2.9, 2.93], fontFace: TM, fontSize: 13, color: INK, rowH: 0.34, border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } }); }

// 7 summary
{ const s = frame('7. Summary');
  bullets(s, [
    'The four-species case runs in the coupled cube examples without any change of the equations; only the point-model case and its time link differ from the two-species examples.',
    'All four species coexist until the antibiotic; afterwards species 4 (b₄ = 0.01) takes over in both starts, as in the point model alone.',
    'Species 4 reaches φ̄₄ = 0.75 (separate) and 0.79 (mixed) at T* = 1 against 0.93 for the point model alone, because the biofilm consumes the nutrient (c ≈ 0.6 in the seed, so c* = c*₀c is lower than 100).',
    'The mixed start keeps more nutrient in the seed and ends with about half the seed stress of the separate start.'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 }); }

// 8 references
{ const s = frame('References');
  bullets(s, [
    'F. Klempt, M. Soleimani, P. Wriggers, P. Junker: A Hamilton principle-based model for diffusion-driven biofilm growth. Biomech. Model. Mechanobiol. 23, 2091–2113 (2024).',
    'F. Klempt, H. Geisler, M. Soleimani, P. Junker: A continuum multi-species bacterial growth model with a novel interaction scheme. Arch. Appl. Mech. 96, 164 (2026).'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4, fontSize: 14 }); }

pres.writeFile({ fileName: B + 'four_species_en.pptx' }).then(() => console.log('ok ' + B + 'four_species_en.pptx ' + n + ' slides'));
