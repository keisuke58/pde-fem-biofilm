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
const SL = [];   // every slide in order, for the speaker notes at the end

function frame(title) {
  const s = pres.addSlide(); n += 1; SL.push(s); s.background = { color: 'FFFFFF' };
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
{ const s = pres.addSlide(); n += 1; SL.push(s); s.background = { color: 'FFFFFF' };
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

// 2c the point model alone: Fritsch case I (added 25 Nov)
{ const s = frame('2. Governing equations (cont.): the point model alone, Fritsch et al. 2025, case I');
  const h = img(s, R + 'assets/fig_fritsch_case1.png', 1.55, 1.15, 10.2);
  cap(s, 'The point model alone (no field), with the true parameters of case I: A = [[1, 0.1], [0.1, 1]], b = (1, 2), η = (1, 2), start φ = (0.25, 0.30), '
    + 'c* = 100, 1000 steps of 10⁻⁴ (T* = 1 is the whole run, s = 0.1). Solid: volume fraction φᵢ; dashed: living volume fraction φ̄ᵢ = φᵢψᵢ; dotted: Σᵢφ̄ᵢ. '
    + '(a) no antibiotic; (b) the antibiotic of the coupled examples, α* = 100 from T* = 1/3 (Klempt et al. 2026, case 4); (c) case I as published, α* = 10 from the start.',
    1.15 + h + 0.08);
  bullets(s, [
    'Without antibiotic, species A takes over (end φ₁ = 0.97, φ₂ = 0.02); both stay alive.',
    'With α* = 100 from T* = 1/3 (the setting of all coupled examples), both die within about 0.1 T*: Σᵢφ̄ᵢ falls from 0.96 to 0.04, while φ₁, φ₂ stay (0.51, 0.42).'],
    { x: 0.75, y: 1.15 + h + 0.82, w: 11.8, h: 6.9 - (1.15 + h + 0.82), fontSize: 14, paraSpaceAfter: 4 }); }

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
    'Staggered coupling; Σᵢφ̄ᵢ in the uptake lags one step; c* = c*₀c and α* = α*₀a linear in time over the point model\'s inner steps',
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
  s.addText('Table 1  Means over the seed. Living fraction, share and nutrient on 64³ cells (outer step 0.025); growth and stress extrapolated from 64³/128³/192³ (outer step 0.0025), GCI in brackets. ¹ Mixed stress with the superbee limiter, extrapolated from 128³/192³/256³.',
    { x: 0.55, y: 1.2, w: 12.23, h: 0.55, fontFace: TM, fontSize: 13, color: INK, isTextBox: true, margin: 0 });
  const rows = [
    [H('quantity'), H('separate'), H('separate, antibiotic'), H('mixed'), H('mixed, antibiotic')],
    [L('living volume fraction Σᵢφ̄ᵢ'), C('0.965'), C('0.152', { bold: true, color: ACC }), C('0.961'), C('0.038', { bold: true, color: ACC })],
    [L('share φ₁/(φ₁+φ₂)'), C('0.978'), C('0.951'), C('0.859'), C('0.610')],
    [L('nutrient c'), C('0.536'), C('0.911'), C('0.608'), C('0.973')],
    [L('growth α − 1 [10⁻⁴]'), C('9.17 (0.4 %)'), C('9.19 (0.5 %)'), C('8.74 (4.2 %)'), C('8.71 (4.3 %)')],
    [L('von Mises stress [10⁻⁴ Pa]'), C('6.57 (0.1 %)'), C('6.42 (1.4 %)', { bold: true, color: ACC }), C('4.26 (2.4 %) ¹'), C('3.95 (6.3 %) ¹')]];
  s.addTable(rows, { x: 0.55, y: 1.85, w: 12.23, colW: [3.43, 2.2, 2.2, 2.2, 2.2], fontFace: TM, fontSize: 16, color: INK, rowH: 0.55,
    border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  s.addShape(pres.shapes.RECTANGLE, { x: 0.55, y: 5.4, w: 12.23, h: 1.3, fill: { color: LIGHT }, line: { color: 'D5DCE5', width: 0.5 } });
  s.addText('The antibiotic lowers the living volume fraction by 84 % (separate) and 96 % (mixed). The growth changes by less than 0.3 %. '
    + 'The stress of the separate start falls by 2 %, more than its GCI. In the mixed start (superbee) the 7 % difference is about the two GCIs together (9 %).',
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
    ['mixed, antibiotic', 'von Mises, superbee', '2.39', '3.28', '3.59', '3.67', '1.2', '3.81', '4.6 %'],
    ['mixed', 'α − 1 [10⁻⁴]', '7.46', '8.08', '8.35', '8.45', '0.77', '8.74', '4.2 %']];
  s.addTable([[H('example'), H('quantity'), H('32³'), H('64³'), H('128³'), H('192³'), H('p'), H('extrap.'), H('GCI')],
    ...g.map(r => r.map((t, j) => j < 2 ? L(t) : C(t, j === 8 && parseFloat(t) > 10 ? { bold: true, color: ACC } : {})))],
    { x: 0.55, y: 1.85, w: 12.23, colW: [2.4, 2.4, 1.05, 1.05, 1.05, 1.05, 0.9, 1.1, 1.23], fontFace: TM, fontSize: 14, color: INK, rowH: 0.42,
      border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  s.addText('Stress in von Mises [10⁻⁴ Pa]. The mixed stress needs a sharp edge: its width β/v ≈ 0.004 mm is narrower than a cell up to 192³ (cell Péclet number 1.28), and the minmod limiter smears it. With the less diffusive superbee limiter it converges with p = 1.2 and a GCI of 4.6 % (extrapolated 3.81 × 10⁻⁴ Pa); 256³ gives 3.76, and 128³/192³/256³ extrapolate to 3.95 (GCI 6.3 %). In the separate start, where minmod converges, superbee agrees within about 1 %.',
    { x: 0.55, y: 5.75, w: 12.23, h: 1.1, fontFace: TM, fontSize: 14, color: INK, isTextBox: true, margin: 0, valign: 'top' }); }

// 12b the element against the Python code (added 25 Nov; same run as the thesis, Ch. 4)
{ const s = frame('12b. Your element and the Python code');
  s.addText('Separate start with antibiotic, seed means at T* = 1. Both solve the same equations (Eqs. 34 to 36 with the coupling of this work).',
    { x: 0.55, y: 1.2, w: 12.23, h: 0.5, fontFace: TM, fontSize: 14, color: INK, isTextBox: true, margin: 0 });
  s.addTable([[H(''), H('IKM element (ANSYS)'), H('Python, finite volumes')],
    [L('grids'), C('8³, 16³, 24³ elements'), C('64³, 128³, 192³ cells')],
    [L('front term −r c/(k+c) n_c·∇φ₃D'), C('upwind-biased gradients (MLS)'), C('MUSCL, minmod limiter')],
    [L('growth α − 1 [10⁻⁴]'), C('6.38, 8.01, 8.36  →  8.68 (4.8 %)'), C('8.93, 9.11, 9.15  →  9.19 (0.5 %)')],
    [L('von Mises stress [10⁻⁴ Pa]'), C('2.14, 4.56, 5.15  →  5.8 (16.5 %)'), C('6.00, 6.28, 6.35  →  6.42 (1.4 %)', { bold: true, color: ACC })]],
    { x: 0.55, y: 1.8, w: 12.23, colW: [3.8, 4.2, 4.23], fontFace: TM, fontSize: 15, color: INK, rowH: 0.55,
      border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  cap(s, 'Values on the three grids, then the extrapolated value with its GCI in brackets (Celik et al. 2008).', 4.6);
  bullets(s, [
    'The stress of the element is 10 % below the Python value; the Python value lies inside the element\'s uncertainty of 16.5 %. The growth of the element is 6 % lower.',
    'The element reaches 24³ only, so the grid series of this deck come from the Python code, which solves the same equations on finer grids.',
    'With a constant added diffusion r h/4 in place of the upwind gradients, the element gives about a quarter of this stress on 16³: the front term needs a second-order treatment.'],
    { x: 0.75, y: 5.05, w: 11.8, h: 1.9, fontSize: 15, paraSpaceAfter: 6 }); }

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
    'A small p gives a small denominator and a large GCI: separate stress p = 3.1, 0.1 %; mixed stress (superbee) p = 1.2, 4.6 %.']);
  col(9.05, 3.73, 'Cell Péclet number', 'pe.png', 3.4, [
    'v is the speed of the front term −r c/(k+c) n_c·∇φ₃D (r = 10 mm/T*, k = 1, c ≈ 0.96); N cells per edge.',
    'Edge width of the biofilm about β/v ≈ 0.004 mm.',
    'For Pe_h > 1 the edge is narrower than a cell and the limiter lowers the order. 128³: 1.9, 192³: 1.28, 256³: 0.96.']); }

// appendix C: a nutrient supply that oscillates in time (added 25 Nov; not in the thesis)
{ const s = frame('Appendix C. A nutrient supply that oscillates in time');
  img(s, R + 'assets/fig_nutrient_sine.png', 0.55, 1.25, 7.4);

  s.addTable([[H(''), H('α − 1 [10⁻⁴]'), H('von Mises [10⁻⁴ Pa]')],
    [L('c* = 100'), C('9.06'), C('6.57')],
    [L('c* = 50 + 50 sin(500 t)'), C('8.15  (−10 %)'), C('8.49  (+29 %)')],
    [L('c* = 50'), C('8.24  (−9 %)'), C('8.57  (+31 %)')]],
    { x: 8.2, y: 1.25, w: 4.58, colW: [1.9, 1.3, 1.38], fontFace: TM, fontSize: 12.5, color: INK, rowH: 0.42,
      border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  cap(s, 'Seed means at T* = 1, coupled, 64³.', 3.05, 8.2, 4.58);
  bullets(s, [
    'Half the nutrient: growth −10 %, but stress +30 %. The front r c/(k+c) is slower and leaves a sharper edge.',
    'The sine and its mean give the same result within 1 %: at this period the oscillation acts as its mean.',
    'In the point model alone the sine leaves small ripples, and species B survives in both cases with c* = 50.'],
    { x: 8.2, y: 3.45, w: 4.58, h: 2.3, fontSize: 13, paraSpaceAfter: 5 });
  s.addText('Fritsch case I, no antibiotic. Top: point model alone; bottom: coupled, separate start, 64³, seed means. t is the point-model time (Klempt et al. 2026, case 3). '
    + 'The held value is updated once per outer step of 0.025 T*, five values per period of 0.126 T*.',
    { x: 8.2, y: 5.8, w: 4.58, h: 1.1, fontFace: TM, fontSize: 11, color: GREY, isTextBox: true, margin: 0, valign: 'top' }); }

// 15 references
{ const s = frame('References');
  const refs = [
    'F. Klempt, M. Soleimani, P. Wriggers, P. Junker: A Hamilton principle-based model for diffusion-driven biofilm growth. Biomech. Model. Mechanobiol. 23, 2091–2113 (2024). doi:10.1007/s10237-024-01883-x',
    'F. Klempt, H. Geisler, M. Soleimani, P. Junker: A continuum multi-species bacterial growth model with a novel interaction scheme. Arch. Appl. Mech. 96, 164 (2026). doi:10.1007/s00419-026-03160-y',
    'L. Fritsch, H. Geisler, J. Grashorn, F. Klempt, M. Soleimani, M. Broggi, P. Junker, M. Beer: Bayesian updating of constitutive parameters under hybrid uncertainties with a novel surrogate model applied to biofilms. arXiv:2512.15145 (2025)',
    'I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, P. E. Raad: Procedure for estimation and reporting of uncertainty due to discretization in CFD applications. J. Fluids Eng. 130, 078001 (2008). doi:10.1115/1.2960953'];
  s.addText(refs.map((t, i) => ({ text: '[' + (i + 1) + ']  ' + t, options: { breakLine: i < refs.length - 1 } })),
    { x: 0.75, y: 1.4, w: 11.8, h: 5.2, fontFace: TM, fontSize: 16, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 14 }); }


// Speaker notes: a simple English text to read aloud, then the Japanese meaning (one entry per slide, in order).
const NOTES = [
['Thank you for your time. Today I will show two examples of the coupled model. In both examples there are two bacterial species and an antibiotic. The two examples differ only in where the species start: separate or mixed. For each example I will show the set-up, the results in 3D and over time, and then the grid convergence of growth and stress. At the end I would like to discuss one question with you.',
 'お時間ありがとうございます．今日は連成モデルの例を 2 つ見せます．どちらも 2 種類の菌と抗菌薬の例です．違いは菌の置き方だけで，分けて置くか，混ぜて置くかです．それぞれ，設定，3D の結果と時間変化，そして成長と応力の格子収束を説明します．最後に 1 つ相談したいことがあります．'],
['First, the objective. I put the point model of Klempt and co-workers, 2026, at every integration point. The point model describes the composition at one point: how much of each species there is, and how much of it is alive. I couple it to the biofilm field of the 2024 paper through two fields, a nutrient and an antibiotic. The question is: if the antibiotic changes the composition, how much do the growth and the stress change? The domain is a two millimetre cube. The species data are synthetic, from Fritsch and co-workers, case one, with their true parameters. There are two starts, separate and mixed, and each is run with and without antibiotic. I report means over the seed, and I extrapolate growth and stress from three grids with the grid convergence index.',
 'まず目的です．Klempt ら 2026 の点モデルを各積分点に置きます．点モデルは 1 点での組成，つまり各菌種がどれだけあって，そのうちどれだけ生きているかを表します．これを 2024 年の論文のバイオフィルムの場と，栄養と抗菌薬の 2 つの場を通してつなぎます．問いは，抗菌薬で組成が変わると，成長と応力がどれだけ変わるか，です．領域は 2 mm の立方体，菌のデータは Fritsch らの case I の合成データで，真のパラメータを使います．始め方は分けた例と混ぜた例の 2 つで，それぞれ抗菌薬ありとなしを計算します．種子の平均を示し，成長と応力は 3 つの格子から GCI で外挿します．'],
['These are the governing equations, Equations 34 to 36 of the 2024 paper. Equation 34 is the amount of biofilm, phi 3D. It has diffusion, a growth source, and the front term, which moves the edge of the biofilm towards the nutrient. Equation 35 is the nutrient: it diffuses and is taken up. Equation 36 is the growth alpha, which drives the swelling and so the stress. The only change of this work is in Equation 35: the uptake is weighted by the living volume fraction from the point model, so only living biomass consumes the nutrient. The antibiotic follows the same equation as the nutrient. The point model itself is used as published.',
 '支配方程式です．2024 年の論文の式 34〜36 です．式 34 はバイオフィルムの量 φ₃D で，拡散，成長の源，そしてバイオフィルムの端を栄養の方へ動かす前線の項があります．式 35 は栄養で，拡散して消費されます．式 36 は成長 α で，膨らみ，つまり応力のもとになります．この研究で変えたのは式 35 だけで，消費に点モデルの生きている体積分率を掛けています．つまり生きている菌だけが栄養を使います．抗菌薬も栄養と同じ式に従います．点モデルそのものは論文のまま使っています．'],
['On the left is the point model, Equations 16 to 18 of the 2026 paper. Each species has a volume fraction phi i and a live share psi i. The antibiotic term is shown in red: it lowers the live share of each species. On the right is the chain of effects. The antibiotic field gives alpha star to the point model. The live shares fall, so the living volume fraction falls. Then less nutrient is taken up, so the nutrient rises, and the front speed r c over k plus c becomes larger. One important point: phi 3D, the amount of biofilm, has no antibiotic term. So the antibiotic kills bacteria, but it does not remove biofilm.',
 '左は点モデルで，2026 年の論文の式 16〜18 です．各菌種に体積分率 φᵢ と，生きている割合 ψᵢ があります．赤が抗菌薬の項で，各菌種の生きている割合を下げます．右は効き方の流れです．抗菌薬の場から点モデルに α* が渡ります．生きている割合が下がり，生きている体積分率が下がります．すると栄養の消費が減り，栄養が増え，前線の速さ r c/(k+c) が大きくなります．大事な点として，バイオフィルムの量 φ₃D には抗菌薬の項がありません．抗菌薬は菌を殺しますが，バイオフィルムそのものは減らしません．'],
['Before the coupled examples, this slide shows the point model alone, without any field. These are the species data of the examples: Fritsch and co-workers, case one, with the true parameters. The solid lines are the volume fractions of the two species, the dashed lines are the living parts, and the dotted line is the total living fraction. In panel a, without antibiotic, species A takes over, and both stay alive. In panel b, the antibiotic of the coupled examples starts at one third with a strength of one hundred. Both species die within about a tenth of T star: the living fraction falls from 0.96 to 0.04. But the volume fractions stay, because dead bacteria are still there. Panel c is case one as published, with a weak antibiotic from the start. In the coupled examples, panel b is what happens at each point of the seed.',
 '連成の例の前に，点モデルだけ（場なし）の結果を見せます．例で使う菌のデータ，Fritsch らの case I で，真のパラメータです．実線が 2 つの菌種の体積分率，破線がそのうち生きている部分，点線が生きている体積分率の合計です．(a) 抗菌薬なしでは菌 A が勝ち，両方とも生きたままです．(b) は連成の例と同じ抗菌薬で，T* = 1/3 から強さ 100 です．約 0.1 T* で両方とも死に，生きている量は 0.96 から 0.04 に下がります．ただ，死んだ菌も残るので体積分率はそのままです．(c) は論文どおりの case I で，弱い抗菌薬が最初からかかります．連成の例では，種子の各点で (b) が起きています．'],
['These are the constants. Almost all of them are from Table 2 of the 2024 paper, converted to the two millimetre cube: the diffusion beta, the front speed r, the growth rate, the diffusion and uptake of the nutrient, and the stiffness. The point-model values are from Fritsch, case one, and the antibiotic strength and onset are from case four of the 2026 paper. Only two values are mine: the stiffness floor f of the empty space, and the threshold phi min at which a point model starts. I checked both, and they change the results very little. On the right are the numerics: finite volumes, a second-order upwind scheme for the front term, Runge-Kutta in time, and the stress on sixteen cubed hexahedral elements.',
 '定数の一覧です．ほとんどは 2024 年の論文の表 2 の値を 2 mm の立方体に換算したものです．拡散 β，前線の速さ r，成長の速さ，栄養の拡散と消費，硬さです．点モデルの値は Fritsch の case I，抗菌薬の強さと開始時刻は 2026 年の論文の case 4 から取りました．私が決めた値は 2 つだけで，空いた所の硬さの下限 f と，点モデルを始める閾値 φ_min です．どちらも結果にほとんど影響しないことを確かめました．右は数値解法で，有限体積法，前線の項に 2 次の風上法，時間は Runge-Kutta，応力は 16³ の六面体要素で解きます．'],
['This is the first example, the separate start. Species A is in the seed in the centre, which is thirty-two elements of a quarter millimetre. Species B is in four small blocks along the edges parallel to z. The nutrient comes from the bottom face, the layer at y equals minus one millimetre. The antibiotic comes from all six outer faces, from T star equals one third; before that it is zero. For the mechanics, only three corner nodes are fixed, so the cube can swell freely.',
 '1 つ目の例，分けた例です．菌 A は中央の種子に置きます．種子は 0.25 mm の要素 32 個です．菌 B は z 方向の辺に沿った 4 つの小さなブロックに置きます．栄養は下の面，y = −1 mm の層から入ります．抗菌薬は T* = 1/3 から外側の 6 面すべてから入り，それまでは 0 です．力学では角の 3 つの節点だけを固定するので，立方体は自由に膨らめます．'],
['These are the results of the separate start with antibiotic. The top row is just before the antibiotic, and the bottom row is at the end, T star equals one. One eighth of the cube is cut away so that you can see inside. From left to right: the share of species A, the living volume fraction, the nutrient, the growth alpha minus one, and the von Mises stress. Look at panel b: the living fraction falls strongly after the antibiotic. But look at panels d and e: the growth and the stress look almost the same as before.',
 '分けた例の，抗菌薬ありの結果です．上の段は抗菌薬の直前，下の段は最後の T* = 1 です．中が見えるように立方体の 8 分の 1 を切り取っています．左から，菌 A の割合，生きている体積分率，栄養，成長 α − 1，von Mises 応力です．(b) を見ると，抗菌薬のあと生きている量が大きく下がります．しかし (d) と (e) を見ると，成長と応力はほとんど前と同じです．'],
['This video shows the same run over time, from T star zero to one. The antibiotic starts at one third, and the onset is written on the frames. Watch how the living fraction drops, while the growth keeps increasing at the same rate. You can click on it to play.',
 '同じ計算の T* = 0 から 1 までの動画です．抗菌薬は 1/3 で始まり，その時刻は画面に書いてあります．生きている量は下がるのに，成長は同じ速さで増え続けるところを見てください．クリックすると再生できます．'],
['This is the second example, the mixed start. Here both species start in the seed, in equal amounts, each with a volume fraction of 0.275. There are no separate blocks. All boundary conditions are the same as in the first example: the nutrient from the bottom face, and the antibiotic from all outer faces from one third.',
 '2 つ目の例，混ぜた例です．こちらは 2 つの菌種を種子に同じ量ずつ，体積分率 0.275 ずつ置きます．別のブロックはありません．境界条件は 1 つ目の例と同じで，栄養は下の面から，抗菌薬は 1/3 から外側のすべての面から入ります．'],
['These are the results of the mixed start with antibiotic, in the same layout as before. Here both species are in the seed, so both meet the antibiotic there. The living fraction falls even further than in the separate start, to about four percent. Species B is removed first, so the share of A rises. Again, the growth and the stress change very little.',
 '混ぜた例の抗菌薬ありの結果で，並べ方は前と同じです．こちらは両方の菌種が種子にいるので，両方がそこで抗菌薬に当たります．生きている量は分けた例よりさらに下がり，約 4 % になります．菌 B が先に死ぬので，A の割合が上がります．ここでも成長と応力はほとんど変わりません．'],
['And this is the video of the mixed start. Here you can see the share of species A rise after the antibiotic, because species B dies first.',
 '混ぜた例の動画です．抗菌薬のあと，菌 B が先に死ぬので，菌 A の割合が上がっていくのが見えます．'],
['These curves show the means over the seed over time, for both starts, with antibiotic as dashed lines and without as solid lines. The dotted line marks the onset of the antibiotic. In panel a, the living fraction stays near one without antibiotic, and it drops quickly after one third with antibiotic. In panel b, the share of A rises in the mixed start, because B dies first. In panel c, the nutrient rises after the antibiotic, because fewer living bacteria take it up.',
 '種子の平均の時間変化です．両方の例について，抗菌薬ありを破線，なしを実線で描いています．点線は抗菌薬の開始です．(a) の生きている量は，抗菌薬なしでは 1 近くのまま，抗菌薬ありでは 1/3 のあと急に下がります．(b) では混ぜた例で A の割合が上がります．B が先に死ぬからです．(c) では抗菌薬のあと栄養が上がります．栄養を使う生きた菌が減るからです．'],
['This table is the main result. The first rows are on sixty-four cubed. The antibiotic lowers the living fraction by 84 percent in the separate start and by 96 percent in the mixed start. The last two rows are growth and stress, extrapolated from three grids, with the grid uncertainty in brackets. The growth changes by less than 0.3 percent. The stress of the separate start falls by 2 percent, which is more than its uncertainty. In the mixed start the difference is 7 percent, which is about the size of the two uncertainties together. So a very large change of the composition gives only a small change of the mechanics.',
 'この表が主な結果です．上の行は 64³ の値です．抗菌薬で生きている量は，分けた例で 84 %，混ぜた例で 96 % 下がります．下の 2 行は 3 つの格子から外挿した成長と応力で，括弧の中が格子の不確かさです．成長の変化は 0.3 % 未満です．分けた例の応力は 2 % 下がり，これは不確かさより大きいです．混ぜた例の差は 7 % で，2 つの不確かさを合わせたくらいの大きさです．つまり，組成が大きく変わっても，力学はほとんど変わりません．'],
['This table shows the grid convergence from thirty-two cubed to one hundred and ninety-two cubed. The growth converges well in both examples, with uncertainties below five percent. The stress of the separate start converges very well, with 1.4 percent with antibiotic. The stress of the mixed start is harder: its edge is narrower than a cell, and the minmod limiter smears it. With the less diffusive superbee limiter it converges, with an uncertainty of about five percent. On two hundred and fifty-six cubed both limiters give about 4.0 times ten to the minus four pascal.',
 '32³ から 192³ までの格子収束の表です．成長はどちらの例もよく収束し，不確かさは 5 % 未満です．分けた例の応力はとてもよく収束し，抗菌薬ありで 1.4 % です．混ぜた例の応力は難しく，端がセルより細いため，minmod の制限関数ではなまってしまいます．拡散の少ない superbee を使うと収束し，不確かさは約 5 % です．256³ では，2 つの制限関数がどちらも約 4.0 × 10⁻⁴ Pa を出します．'],
['This slide compares your element with my Python code. Both solve the same equations, for the separate start with antibiotic. Your element runs on eight, sixteen and twenty-four cubed elements, with upwind-biased gradients for the front term. The Python code runs on sixty-four to one hundred and ninety-two cubed cells. The stress of the element extrapolates to 5.8 times ten to the minus four pascal, with a grid uncertainty of 16.5 percent. The Python value is 6.4, ten percent higher, and it lies inside the uncertainty of the element. So the two agree. The growth of the element is six percent lower. Because the element reaches only twenty-four cubed, I took the grid series of this deck from the Python code. One more point: with a constant added diffusion instead of the upwind gradients, the stress is only about a quarter of this value, so the front term needs a second-order treatment.',
 'このスライドは，あなたの要素と私の Python コードの比較です．どちらも同じ式を，分けた例・抗菌薬ありで解いています．要素は 8³，16³，24³ の要素で，前線の項には風上寄りの勾配を使っています．Python は 64³ から 192³ のセルです．要素の応力は外挿で 5.8×10⁻⁴ Pa，格子の不確かさは 16.5 % です．Python の値は 6.4 で 10 % 高く，要素の不確かさの中に入ります．なので両者は合っています．成長は要素のほうが 6 % 低いです．要素は 24³ までしか回らないので，このスライドの格子の列は Python のものを使いました．もう 1 点，風上寄りの勾配の代わりに一定の人工拡散を使うと，応力はこの約 4 分の 1 しか出ません．前線の項には 2 次の扱いが必要です．'],
['Why does the stress hardly change? The growth law, Equation 36, counts all biomass, living or dead. So when bacteria die, the biofilm still grows. The composition reaches the stress only through the nutrient and the front speed. The stress itself is set by the difference in growth between the seed and the stiff biofilm around it. For the composition to act on the mechanics directly, the model needs a second path, for example a growth law driven by the living biomass. I would like to hear your opinion: is this worth pursuing in the next step?',
 'なぜ応力がほとんど変わらないのか．成長の式 36 は，生きている菌も死んだ菌も同じに数えます．なので菌が死んでも，バイオフィルムは成長を続けます．組成は，栄養と前線の速さを通してしか応力に届きません．応力そのものは，種子と周りの硬いバイオフィルムとの成長の差で決まります．組成が力学に直接効くには，生きている量で成長を決めるような 2 本目の道が必要です．次の段階でこれを進める価値があるか，意見を聞きたいです．'],
['This appendix lists the notation of the growth field and of the point model, in case you have questions about a symbol.',
 '付録：成長の場と点モデルの記号の一覧です．記号について質問があったときに使います．'],
['This appendix defines the three quantities I report: the von Mises stress, the grid convergence index of Celik and co-workers, and the cell Péclet number, which tells whether the edge of the biofilm is resolved by the grid.',
 '付録：報告した 3 つの量の定義です．von Mises 応力，Celik らの格子収束指標（GCI），そしてバイオフィルムの端が格子で解像できているかを示すセルのペクレ数です．'],
['This appendix shows a nutrient supply that oscillates in time, in the form of case three of the 2026 paper: c star equals fifty plus fifty times the sine of five hundred t. The species are again Fritsch case one, without antibiotic. The top row is the point model alone. With the sine, the curves have small ripples, and they follow closely the curves with a constant c star of fifty, which is the mean of the sine. The bottom row is the coupled model. The nutrient at the seed oscillates, but the living fraction and the share of species A are almost the same as with the constant mean. In the table, halving the nutrient lowers the growth by ten percent but raises the stress by about thirty percent, because the front becomes slower and the edge stays sharp. The sine and its mean agree within one percent. So at this period, the oscillation acts like its mean.',
 '付録：栄養を時間で振る例です．2026 年の論文の case 3 の形で，c* = 50 + 50 sin(500 t) です．菌は同じ Fritsch case I で，抗菌薬はありません．上の段は点モデルだけです．正弦にすると曲線に小さな波が乗りますが，正弦の平均である c* = 50 の一定のときとほぼ同じ曲線になります．下の段は連成モデルです．種子の栄養は振れますが，生きている量と菌 A の割合は，平均の一定のときとほぼ同じです．表を見ると，栄養を半分にすると成長は 10 % 下がるのに，応力は約 30 % 上がります．前線が遅くなり，端が鋭いまま残るからです．正弦とその平均は 1 % 以内で一致します．つまりこの周期では，振動は平均として効きます．'],
['These are the references. Thank you very much for your attention. I am happy to answer questions.',
 '参考文献です．ご清聴ありがとうございました．質問があればどうぞ．']];
if (NOTES.length !== SL.length) throw new Error('notes ' + NOTES.length + ' for ' + SL.length + ' slides');
SL.forEach((sl, i) => sl.addNotes('EN: ' + NOTES[i][0] + '\n\nJA: ' + NOTES[i][1]));

pres.writeFile({ fileName: B + 'coupled_examples_en.pptx' }).then(f => console.log('ok', f, n, 'slides'));
