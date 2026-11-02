// English deck: the two coupled examples (separate and mixed start, with antibiotic) in one set, for Oliver,
// Dr. Soleimani and Felix Klempt. Built by build.sh: node deck_en.js -> build/coupled_examples_en.pptx
// Numbers: thesis Ch. 4, tab:ch5-two-species and tab:ch5-gci (64/128/192^3); check them again before sending.
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const path = require('path');
const R = path.resolve(__dirname, '../..') + '/';
const B = __dirname + '/build/', EQ = B + 'eq/';
const vid = (c) => R + 'assets/anim_coupled_tvd_' + c + '_ab_inhliv_transc_alpha_5panel.mp4';   // panel (d) = alpha - 1 (2 Nov)
const poster = (c) => R + 'assets/anim_coupled_tvd_' + c + '_ab_inhliv_transc_alpha_5panel_end.png';

const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka';
pres.title = 'Two species with an antibiotic in the coupled growth model';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', RULE = '8A99AD', ACC = 'A61E1E', LIGHT = 'F3F5F8';
const TM = 'Times New Roman';
let n = 0;

function frame(title) {
  const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addText(title, { x: 0.55, y: 0.32, w: 12.2, h: 0.6, fontFace: TM, fontSize: 26, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 0.98, w: 12.23, h: 0, line: { color: NAVY, width: 1.25 } });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 7.0, w: 12.23, h: 0, line: { color: RULE, width: 0.5 } });
  s.addText('K. Nishioka  |  Two species with an antibiotic in the coupled growth model', { x: 0.55, y: 7.05, w: 10, h: 0.3, fontFace: TM, fontSize: 11, color: GREY, isTextBox: true, margin: 0 });
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

// 1 title
{ const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 2.0, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText('Two species with an antibiotic in the coupled growth model', { x: 0.9, y: 2.2, w: 11.5, h: 0.7, fontFace: TM, fontSize: 32, bold: true, color: NAVY, isTextBox: true, margin: 0, fit: 'shrink' });
  s.addText('Separate and mixed start: set-up, results and grid convergence', { x: 0.9, y: 3.0, w: 11.5, h: 0.6, fontFace: TM, fontSize: 22, color: INK, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 3.8, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText([{ text: 'Keisuke Nishioka', options: { fontSize: 20, bold: true, breakLine: true } },
    { text: 'Institute of Continuum Mechanics, Leibniz Universität Hannover (master\'s thesis, Chapter 4)', options: { fontSize: 16, breakLine: true } },
    { text: 'November 2026', options: { fontSize: 16 } }],
    { x: 0.9, y: 4.25, w: 11.5, h: 1.5, fontFace: TM, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 6 }); }

// 2 objective
{ const s = frame('1. Objective and set-up');
  bullets(s, [
    'I place the point model of Klempt et al. 2026 at every integration point and couple it to the biofilm field of Klempt et al. 2024 (Eqs. 34–36) through a nutrient and an antibiotic.',
    'Question: how much does a change of the composition, here caused by an antibiotic, change the growth and the stress?',
    'Domain: a 2 mm cube. Species: synthetic data of Fritsch et al. 2025, case I (two species, true parameters).',
    'Two start configurations: separate (species A in the seed, B in four small blocks along the edges parallel to z) and mixed (A and B in equal amounts in the seed).',
    'Each with and without antibiotic; reported are means over the seed. Growth and stress are extrapolated from 64³, 128³ and 192³ cells with the grid convergence index (GCI).'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 }); }

// 3 equations
{ const s = frame('2. Governing equations');
  const h = img(s, EQ + 'eqs.png', 1.2, 1.3, 10.6);
  cap(s, 'φ₃D: biofilm volume fraction; α: growth (F_g = α I, α(0) = 1); c: nutrient; a: antibiotic; φᵢ, ψᵢ: volume and living fraction of species i in the point model; '
    + 'Σᵢφ̄ᵢ: living volume fraction. Equation numbers of Klempt et al. 2024. The uptake in Eq. 35 is weighted by the living volume fraction (the coupling of this work); '
    + 'the antibiotic follows the same equation. The point model is used as published (A, b, η of Fritsch 2025, case I). '
    + 'Constants not shown: ℂ is the isotropic elasticity tensor (Poisson ratio ν); the point model runs only where φ₃D > φ_min and advances s ΔT* per field step ΔT*; '
    + 'g is given as gℓ²/d with the half edge ℓ = 1 mm of the cube.', 1.4 + h + 0.2); }

// 2b point model and where the antibiotic acts
{ const s = frame('2. Governing equations (cont.): the point model and the antibiotic');
  const h = img(s, EQ + 'point.png', 0.55, 1.3, 6.9);
  cap(s, 'Eqs. 16 to 18 of Klempt et al. 2026, rearranged (barrier terms omitted). φᵢ, ψᵢ: volume and living fraction of species i, φ̄ᵢ = φᵢψᵢ, '
    + 'φ₀: void, γ: multiplier of the constraint, A: interactions, bᵢ: sensitivity to the antibiotic, ηᵢ: viscosity. Red: the antibiotic term. '
    + 'Last line: rate of ψᵢ when one species fills its point; a thinly present species dies fast. α* (antibiotic, point model) and α (growth, Eq. 36) are different quantities.',
    1.3 + h + 0.25, 0.55, 6.9);
  const steps = ['antibiotic field a (outer shell, T* ≥ 1/3)', 'α* = α*₀ a passed to the point model', 'the term −α* bᵢ ψᵢ of Eq. 17 lowers ψᵢ',
    'living volume fraction Σᵢφ̄ᵢ = Σᵢφᵢψᵢ falls', 'uptake g φ₃D (Σᵢφ̄ᵢ) c of Eq. 35 falls, nutrient c rises', 'front speed r c/(k+c) rises (Eq. 34)'];
  const X = 7.85, W = 4.93, Hb = 0.55, G = 0.22;
  s.addText('How the antibiotic reaches the field', { x: X, y: 1.25, w: W, h: 0.4, fontFace: TM, fontSize: 17, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  steps.forEach((t, i) => {
    const y = 1.7 + i * (Hb + G);
    s.addShape(pres.shapes.RECTANGLE, { x: X, y, w: W, h: Hb, fill: { color: i === 2 ? 'FBEAEA' : LIGHT }, line: { color: i === 2 ? ACC : RULE, width: 0.75 } });
    s.addText(t, { x: X + 0.12, y, w: W - 0.24, h: Hb, fontFace: TM, fontSize: 14, color: INK, valign: 'middle', isTextBox: true, margin: 0 });
    if (i < steps.length - 1) s.addShape(pres.shapes.LINE, { x: X + W / 2, y: y + Hb, w: 0, h: G, line: { color: NAVY, width: 1.25, endArrowType: 'triangle' } });
  });
  s.addText('φ₃D (living and dead biomass) has no antibiotic term, so the antibiotic does not remove biofilm.',
    { x: X, y: 1.7 + 6 * (Hb + G) - 0.12, w: W, h: 0.6, fontFace: TM, fontSize: 13, color: ACC, isTextBox: true, margin: 0, valign: 'top' });
}

// 4 constants
{ const s = frame('3. Constants and numerics');
  const rows = [[H('symbol'), H('value'), H('meaning'), H('source')],
    ['β', '0.02 mm²/T*', 'diffusion of φ₃D', 'Klempt 2024, Table 2, converted to 2 mm'],
    ['r, k', '10 mm/T*, 1', 'front speed, Monod constant', 'same'],
    ['k_α', '10⁻³ /T*', 'growth rate', 'same'],
    ['d, gℓ²/d', '10⁸ mm²/T*, 1', 'diffusion and uptake of c and a', 'same (ℓ = 1 mm)'],
    ['E, ν', '10 Pa, 0.49', 'stiffness of the biofilm', 'same'],
    ['f', '10⁻³', 'stiffness floor of the void', 'this work (f → 0: 0.13 %)'],
    ['c*₀, s', '100, 0.1', 'nutrient of the point model, time link', 'Fritsch 2025, case I'],
    ['α*₀', '100 from T* = 1/3', 'antibiotic strength and onset', 'Klempt 2026, case 4'],
    ['φ_min', '10⁻³', 'threshold to start a point model', 'this work (converged)']]
    .map((r, i) => i === 0 ? r : r.map((t, j) => ({ text: t, options: { align: j >= 2 ? 'left' : 'center' } })));
  s.addTable(rows, { x: 0.55, y: 1.3, w: 7.9, colW: [1.1, 1.9, 2.4, 2.5], fontFace: TM, fontSize: 13, color: INK, rowH: 0.42,
    border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  s.addText('Numerics', { x: 8.85, y: 1.3, w: 4.0, h: 0.4, fontFace: TM, fontSize: 18, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  bullets(s, ['Finite volumes, cell centres, 64³ (a quarter of the cube by symmetry)',
    'Front term: second-order upwind (MUSCL, minmod limiter)',
    'Time: RK4 for φ₃D and α; backward Euler for Eq. 35 (with ċ); implicit Euler with Newton in the point model',
    'Staggered coupling; Σᵢφ̄ᵢ in the uptake lags one step',
    'Stress: 16³ hexahedra, B-bar, small strain',
    'Wall time, one core, 400 steps: 64³ about 5 min, 128³ about 1 h, 192³ about 3.3 h'],
    { x: 8.85, y: 1.8, w: 4.0, h: 5.0, fontSize: 14, paraSpaceAfter: 7 }); }

// 4-9: example 1 (separate start), then example 2 (mixed start): set-up, results, video
['separate', 'mixed'].forEach((c, k) => {
  const ex = 'Example ' + (k + 1) + ', ' + c + ' start';
  {
  const s = frame((4 + 3 * k) + '. ' + ex + ': initial and boundary conditions');
  const h = img(s, R + 'assets/fig_setup_' + c + '_ab.png', 1.6, 1.08, 10.1);
  const y0 = 1.08 + h + 0.05;
  const c_ = (t, o = {}) => ({ text: t, options: Object.assign({ align: 'left', valign: 'middle' }, o) });
  s.addTable([[H('field'), H('where it is held'), H('value'), H('from'), H('other walls')],
    [c_('nutrient c', { bold: true }), c_('element layer on the face y = −1 mm (0.25 mm thick, 64 elements)'), c_('c = 1'), c_('T* = 0'), c_('zero flux')],
    [c_('antibiotic a', { bold: true }), c_('outer element shell on all six faces (0.25 mm thick, 296 elements)'), c_('a = 1'), c_('T* = 1/3 (a = 0 before)'), c_('— (held on all)')],
    [c_('displacement u', { bold: true }), c_('three corner nodes (▲: ux, uy, uz / uy, uz / uz)'), c_('0'), c_('—'), c_('free swelling')]],
    { x: 0.55, y: y0, w: 12.23, colW: [1.6, 5.2, 0.95, 2.4, 2.08], fontFace: TM, fontSize: 12.5, color: INK, rowH: 0.33,
      border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  cap(s, 'Fig. ' + (3 * k + 1) + '  ' + (c === 'separate'
    ? 'Separate start: species A in the seed (32 elements of 0.25 mm, 0.5 mm³), species B in four blocks of 0.25 × 0.25 × 1.0 mm along the edges parallel to z (z = −0.5 to 0.5 mm, 0.25 mm³ in all).'
    : 'Mixed start: species A and B in equal amounts (φ₁ = φ₂ = 0.275) in the seed (32 elements of 0.25 mm, 0.5 mm³).')
    + ' c and a diffuse inward from the held layers and are taken up by the living biomass (Eq. 35). Drawn on 8³ elements; the runs use 64³ cells.', y0 + 1.42);
  }
  {
  const s = frame((5 + 3 * k) + '. ' + ex + ': results with antibiotic');
  const h = img(s, B + (c === 'separate' ? 's12_coupled_3d_separate_ab.png' : 's12m_coupled_3d_mixed_ab.png'), 2.05, 1.1, 9.2);
  cap(s, 'Fig. ' + (3 * k + 2) + '  ' + (c === 'separate' ? 'Separate' : 'Mixed') + ' start with antibiotic, T* = 0.325 (just before the antibiotic, top) and T* = 1 (bottom), one octant removed. '
    + '(a) share φ₁/(φ₁+φ₂) and (b) living volume fraction Σᵢφ̄ᵢ where φ₃D > 0.5; (c) nutrient c; (d) growth α − 1 where φ₃D > 0.5 (as in the video); (e) von Mises stress on 16³ elements. 64³ cells.', 1.12 + h + 0.08);
  }
  {
  const s = frame((6 + 3 * k) + '. ' + ex + ': time evolution (T* = 0 to 1)');
  const [pw, ph] = pngSize(poster(c)); const w = 12.2, hh = w * ph / pw;
  if (fs.existsSync(vid(c))) s.addMedia({ type: 'video', path: vid(c), cover: 'data:image/png;base64,' + fs.readFileSync(poster(c)).toString('base64'), x: 0.57, y: 1.3, w, h: hh });
  else s.addImage({ path: poster(c), x: 0.57, y: 1.3, w, h: hh });
  cap(s, 'Fig. ' + (3 * k + 3) + '  The run of Fig. ' + (3 * k + 2) + ' over time (40 frames; panel (d) shows the growth α − 1 instead of the antibiotic, whose onset is marked on the frames; the caption is part of the video). Click to play.', 1.3 + hh + 0.15);
  }
});

// 11 time history
{ const s = frame('10. Seed means over time');
  const h = img(s, R + 'assets/fig_coupled_time_history.png', 0.9, 1.25, 11.5);
  cap(s, 'Fig. 7  Means over the seed for both starts, with (dashed) and without (solid) antibiotic: (a) living volume fraction Σᵢφ̄ᵢ, '
    + '(b) share φ₁/(φ₁+φ₂), (c) nutrient c. The antibiotic acts from T* = 1/3 (dotted). 64³ cells, outer step 0.025, the runs of Table 1.', 1.25 + h + 0.1); }

// 10 seed means
{ const s = frame('11. Seed means at T* = 1');
  s.addText('Table 1  Means over the seed. Living fraction, share and nutrient on 64³ cells (outer step 0.025); growth and stress extrapolated from 64³/128³/192³ (outer step 0.0025), GCI in brackets.',
    { x: 0.55, y: 1.2, w: 12.23, h: 0.55, fontFace: TM, fontSize: 13, color: INK, isTextBox: true, margin: 0 });
  const rows = [
    [H('quantity'), H('separate'), H('separate, antibiotic'), H('mixed'), H('mixed, antibiotic')],
    [L('living volume fraction Σᵢφ̄ᵢ'), C('0.965'), C('0.152', { bold: true, color: ACC }), C('0.961'), C('0.038', { bold: true, color: ACC })],
    [L('share φ₁/(φ₁+φ₂)'), C('0.978'), C('0.951'), C('0.859'), C('0.610')],
    [L('nutrient c'), C('0.536'), C('0.911'), C('0.608'), C('0.973')],
    [L('growth α − 1 [10⁻⁴]'), C('9.17 (0.4 %)'), C('9.19 (0.5 %)'), C('8.74 (4.2 %)'), C('8.71 (4.3 %)')],
    [L('von Mises stress [10⁻⁴ Pa]'), C('6.57 (0.1 %)'), C('6.42 (1.4 %)', { bold: true, color: ACC }), C('4.86 (31 %)'), C('4.68 (41 %)')]];
  s.addTable(rows, { x: 0.55, y: 1.85, w: 12.23, colW: [3.43, 2.2, 2.2, 2.2, 2.2], fontFace: TM, fontSize: 16, color: INK, rowH: 0.55,
    border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  s.addShape(pres.shapes.RECTANGLE, { x: 0.55, y: 5.4, w: 12.23, h: 1.3, fill: { color: LIGHT }, line: { color: 'D5DCE5', width: 0.5 } });
  s.addText('The antibiotic lowers the living volume fraction by 84 % (separate) and 96 % (mixed). The growth changes by less than 0.3 %. '
    + 'The stress of the separate start falls by 2 %, more than its GCI; the 4 % difference of the mixed start lies within its GCI of 31–41 %.',
    { x: 0.8, y: 5.45, w: 11.75, h: 1.2, fontFace: TM, fontSize: 17, color: INK, isTextBox: true, margin: 0, valign: 'middle' }); }

// 11 grid convergence
{ const s = frame('12. Grid convergence of growth and stress');
  s.addText('Table 2  Seed means on 32³ to 192³ cells, observed order p, extrapolated value and GCI of the finest grid (Celik et al. 2008, three grids 64³/128³/192³).',
    { x: 0.55, y: 1.2, w: 12.23, h: 0.55, fontFace: TM, fontSize: 13, color: INK, isTextBox: true, margin: 0 });
  const g = [
    ['separate, antibiotic', 'α − 1 [10⁻⁴]', '8.453', '8.931', '9.112', '9.152', '1.7', '9.19', '0.5 %'],
    ['separate, antibiotic', 'von Mises [10⁻⁴ Pa]', '5.523', '6.000', '6.284', '6.349', '1.6', '6.42', '1.4 %'],
    ['separate', 'α − 1 [10⁻⁴]', '8.729', '9.014', '9.110', '9.135', '1.5', '9.17', '0.4 %'],
    ['separate', 'von Mises [10⁻⁴ Pa]', '7.440', '6.730', '6.585', '6.571', '3.1', '6.57', '0.1 %'],
    ['mixed, antibiotic', 'α − 1 [10⁻⁴]', '7.42', '8.04', '8.32', '8.42', '0.77', '8.71', '4.3 %'],
    ['mixed, antibiotic', 'von Mises [10⁻⁴ Pa]', '2.00', '2.82', '3.30', '3.52', '0.43', '4.68', '41 %'],
    ['mixed', 'α − 1 [10⁻⁴]', '7.46', '8.08', '8.35', '8.45', '0.77', '8.74', '4.2 %'],
    ['mixed', 'von Mises [10⁻⁴ Pa]', '2.13', '3.08', '3.65', '3.89', '0.55', '4.86', '31 %']];
  s.addTable([[H('example'), H('quantity'), H('32³'), H('64³'), H('128³'), H('192³'), H('p'), H('extrap.'), H('GCI')],
    ...g.map(r => r.map((t, j) => j < 2 ? L(t) : C(t, j === 8 && parseFloat(t) > 10 ? { bold: true, color: ACC } : {})))],
    { x: 0.55, y: 1.85, w: 12.23, colW: [2.4, 2.4, 1.05, 1.05, 1.05, 1.05, 0.9, 1.1, 1.23], fontFace: TM, fontSize: 14, color: INK, rowH: 0.42,
      border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  s.addText('The mixed stress converges slowly (p ≈ 0.5): its edge, of width β/v ≈ 0.004 mm, is narrower than a cell up to 192³ (cell Péclet number 1.28). On 256³ (Pe_h = 0.96, with ċ) the mixed stress with antibiotic is 3.63 (equal) and 4.50 × 10⁻⁴ Pa (Fritsch shares); on 128³/192³/256³ p rises to 0.93 and 1.3 and the GCI falls to 13 % and 9 % (extrapolated 4.0 and 4.8 × 10⁻⁴ Pa). Without antibiotic (equal) 256³ gives 4.02 × 10⁻⁴ Pa, GCI 14 %, extrapolated 4.48: 12 % above the run with antibiotic, within the two bands.',
    { x: 0.55, y: 5.75, w: 12.23, h: 1.1, fontFace: TM, fontSize: 14, color: INK, isTextBox: true, margin: 0, valign: 'top' }); }

// 12 discussion
{ const s = frame('13. Discussion');
  bullets(s, [
    'The composition reaches the stress only through the nutrient and the front speed r c/(k+c). The growth law, Eq. 36, counts all biomass, living or dead, so a large change of the living fraction leaves growth and stress almost unchanged.',
    'The stress is set not by the amount of growth but by its mismatch in stiff biofilm, (φ₃D²+f)|∇α| (correlation 0.94 over 51 runs on 64³). α records when the front arrived, so a lower nutrient and a slower front leave a larger difference between the seed and its surroundings and a higher stress. The mixed seed carries less because the biofilm around it is softer (intermediate φ₃D).',
    'Species B is removed first: its nutrient term c*A₂₂/η₂ is half that of A, while the antibiotic term per viscosity, α*bᵢ/ηᵢ, is the same for both. Where the nutrient is low, A dies as well.',
    'In the mixed start both species meet the antibiotic in the seed, and the living volume fraction falls further (0.038 against 0.152).',
    'For the composition to act on the mechanics directly, the model needs a second path, for example a growth law driven by the living biomass. I would like to discuss whether this is worth pursuing.'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 }); }

// 13 appendix notation
{ const s = frame('Appendix A. Notation');
  const rows = [[H('symbol'), H('meaning'), H('symbol'), H('meaning')],
    ['φ₃D', 'biofilm volume fraction (field, Klempt 2024)', 'φᵢ', 'volume fraction of species i (point model)'],
    ['α', 'growth, F_g = α I, α(0) = 1', 'ψᵢ', 'living fraction of species i'],
    ['c', 'nutrient (field), held at 1 on its source', 'φ̄ᵢ = φᵢψᵢ', 'living volume fraction of species i'],
    ['a', 'antibiotic (field), held at 1 on its source', 'φ₀, γ', 'auxiliary void and its multiplier'],
    ['n_c', '∇c/|∇c|, direction of the front term', 'Aᵢⱼ, bᵢ, ηᵢ', 'interactions, antibiotic sensitivity, viscosity'],
    ['T*', 'time; T* = 1 is the duration of the point model\'s data', 'c*, α*', 'inputs: c* = c*₀c, α* = α*₀a'],
    ['σ, ε', 'stress and strain (small strain)', 's', 'time link of the point model']]
    .map((r, i) => i === 0 ? r : r.map((t, j) => ({ text: t, options: { align: j % 2 ? 'left' : 'center' } })));
  s.addTable(rows, { x: 0.55, y: 1.35, w: 12.23, colW: [1.3, 4.8, 1.5, 4.63], fontFace: TM, fontSize: 15, color: INK, rowH: 0.55,
    border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } }); }

// 14 appendix definitions
{ const s = frame('Appendix B. Definitions of the reported quantities');
  const col = (x, w, t, eq, eqw, body) => {
    s.addShape(pres.shapes.RECTANGLE, { x, y: 1.25, w, h: 5.55, fill: { color: LIGHT }, line: { color: 'D5DCE5', width: 0.5 } });
    s.addText(t, { x: x + 0.2, y: 1.38, w: w - 0.4, h: 0.45, fontFace: TM, fontSize: 19, bold: true, color: NAVY, isTextBox: true, margin: 0 });
    const h = img(s, EQ + eq, x + (w - eqw) / 2, 1.95, eqw);
    s.addText(body.map((t, i) => ({ text: t, options: { breakLine: i < body.length - 1 } })),
      { x: x + 0.2, y: 2.1 + h, w: w - 0.4, h: 4.6 - h, fontFace: TM, fontSize: 14, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 6 }); };
  col(0.55, 3.75, 'von Mises stress', 'vm.png', 3.35, [
    's is the deviatoric stress; σ_vM vanishes for a uniform swelling.',
    'Tables 1 and 2 give the mean over the seed (the 32 elements of the start).']);
  col(4.45, 4.5, 'Grid convergence index', 'gci.png', 4.1, [
    'Celik et al. 2008, three grids 1, 2, 3 = 192³, 128³, 64³, safety factor 1.25.',
    'p is the observed order (solved iteratively), f_ext the extrapolated value.',
    'A small p gives a small denominator and a large GCI: separate stress p = 3.1, 0.1 %; mixed stress p = 0.55, 31 %.']);
  col(9.05, 3.73, 'Cell Péclet number', 'pe.png', 3.4, [
    'v is the speed of the front term −r c/(k+c) n_c·∇φ₃D (r = 10 mm/T*, k = 1, c ≈ 0.96); N cells per edge.',
    'Edge width of the biofilm about β/v ≈ 0.004 mm.',
    'For Pe_h > 1 the edge is narrower than a cell and the limiter lowers the order. 128³: 1.9, 192³: 1.28, 256³: 0.96.']); }

// 15 references
{ const s = frame('References');
  const refs = [
    'F. Klempt, M. Soleimani, P. Wriggers, P. Junker: A Hamilton principle-based model for diffusion-driven biofilm growth. Biomech. Model. Mechanobiol. 23, 2091–2113 (2024). doi:10.1007/s10237-024-01883-x',
    'F. Klempt, H. Geisler, M. Soleimani, P. Junker: A continuum multi-species bacterial growth model with a novel interaction scheme. Arch. Appl. Mech. 96, 164 (2026). doi:10.1007/s00419-026-03160-y',
    'L. Fritsch, H. Geisler, J. Grashorn, F. Klempt, M. Soleimani, M. Broggi, P. Junker, M. Beer: Bayesian updating of constitutive parameters under hybrid uncertainties with a novel surrogate model applied to biofilms. arXiv:2512.15145 (2025)',
    'I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, P. E. Raad: Procedure for estimation and reporting of uncertainty due to discretization in CFD applications. J. Fluids Eng. 130, 078001 (2008). doi:10.1115/1.2960953'];
  s.addText(refs.map((t, i) => ({ text: '[' + (i + 1) + ']  ' + t, options: { breakLine: i < refs.length - 1 } })),
    { x: 0.75, y: 1.4, w: 11.8, h: 5.2, fontFace: TM, fontSize: 16, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 14 }); }

pres.writeFile({ fileName: B + 'coupled_examples_en.pptx' }).then(f => console.log('ok', f, n, 'slides'));
