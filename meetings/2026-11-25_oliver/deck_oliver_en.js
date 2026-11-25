// Meeting with Oliver, 25 Nov 2026: short status and three questions. Slides carry the sentences to read aloud;
// the notes repeat them with the Japanese meaning.
//   NODE_PATH=../2026-11_coupled_examples/build/node_modules node deck_oliver_en.js   -> build/oliver_1125_en.pptx (not in git)
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const B = __dirname + '/build/';
if (!fs.existsSync(B)) fs.mkdirSync(B, { recursive: true });

const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka';
pres.title = 'Status and questions, 25 Nov 2026';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', RULE = '8A99AD', LIGHT = 'EEF2F6', ACC = '0F6E7A';
const F = 'Calibri';
const FOOT = 'Keisuke Nishioka  |  Meeting with Oliver, 25 Nov 2026';
let n = 0;

function frame(title) {
  const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addText(title, { x: 0.6, y: 0.35, w: 12.1, h: 0.7, fontFace: F, fontSize: 28, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.6, y: 1.1, w: 12.13, h: 0, line: { color: NAVY, width: 1.25 } });
  s.addShape(pres.shapes.LINE, { x: 0.6, y: 7.0, w: 12.13, h: 0, line: { color: RULE, width: 0.5 } });
  s.addText(FOOT, { x: 0.6, y: 7.05, w: 10.5, h: 0.3, fontFace: F, fontSize: 11, color: GREY, isTextBox: true, margin: 0 });
  s.addText(String(n), { x: 11.73, y: 7.05, w: 1.0, h: 0.3, fontFace: F, fontSize: 11, color: GREY, align: 'right', isTextBox: true, margin: 0 });
  return s;
}
const lines = (s, items, o = {}) => s.addText(items.map((t, i) => ({ text: t, options: { bullet: { indent: 20 }, breakLine: i < items.length - 1 } })),
  Object.assign({ x: 0.8, y: 1.45, w: 11.8, h: 5.3, fontFace: F, fontSize: 22, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 18 }, o));
const TB = { fontFace: F, fontSize: 18, color: INK, border: { type: 'solid', color: 'BBBBBB', pt: 0.5 }, align: 'center', valign: 'middle', rowH: 0.55 };
const H = t => ({ text: t, options: { bold: true, fill: { color: LIGHT } } });

// 1 title
{ const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 2.3, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText('Status and three questions', { x: 0.9, y: 2.5, w: 11.5, h: 0.9, fontFace: F, fontSize: 36, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addText('Meeting with Oliver, 25 November 2026', { x: 0.9, y: 3.4, w: 11.5, h: 0.6, fontFace: F, fontSize: 22, color: INK, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 4.2, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText('Keisuke Nishioka', { x: 0.9, y: 4.5, w: 11.5, h: 0.5, fontFace: F, fontSize: 18, color: GREY, isTextBox: true, margin: 0 });
  s.addNotes('EN: Thank you for your time. Today I have a short status update and three questions about your element.\n' +
    'JA: お時間ありがとうございます．今日は短い近況報告と，要素についての質問が 3 つあります．');
}

// 2 thesis
{ const s = frame('1. Status of the thesis');
  lines(s, [
    'I went through the whole thesis again and fixed the inconsistencies I found.',
    'It has 89 pages now.',
    'I will send the PDF to Meisam by 15 December.',
    'The colloquium is on 16 December, 13:00, in room 320.']);
  s.addNotes('EN: I went through the whole thesis again and fixed the inconsistencies I found. It has 89 pages now. ' +
    'I will send the PDF to Meisam by the fifteenth of December. The colloquium is on the sixteenth of December at one o\'clock, in room 320.\n' +
    'JA: 修論を全部読み直して，見つけた食い違いを直しました．今 89 ページです．12 月 15 日までに Meisam に PDF を送ります．コロキウムは 12 月 16 日 13 時，320 号室です．');
}

// 3 element vs python
{ const s = frame('2. Your element and the Python code');
  s.addText('Separate example with antibiotic, seed mean at T* = 1, extrapolated from three grids', { x: 0.8, y: 1.35, w: 11.8, h: 0.45, fontFace: F, fontSize: 16, color: GREY, isTextBox: true, margin: 0 });
  s.addTable([
    [H(''), H('Your element (8³, 16³, 24³)'), H('Python code (up to 192³)')],
    ['Growth α − 1', '8.7 × 10⁻⁴', '9.2 × 10⁻⁴'],
    ['von Mises stress', '5.8 × 10⁻⁴ Pa', '6.4 × 10⁻⁴ Pa'],
    ['Grid uncertainty of the stress', '16.5 %', '1.4 %']],
    Object.assign({}, TB, { x: 0.8, y: 1.95, w: 11.8, colW: [4.0, 3.9, 3.9] }));
  lines(s, [
    'The stress of the element is 10 % lower than the Python value.',
    'The Python value lies inside the uncertainty of the element, so the two agree.',
    'The growth of the element is 6 % lower.'], { y: 4.4, h: 2.5, fontSize: 20, paraSpaceAfter: 12 });
  s.addNotes('EN: This compares your element with my Python code for the same equations. The stress of the element is ten percent lower than the Python value. ' +
    'But the element only reaches twenty-four cubed, so its grid uncertainty is sixteen and a half percent, and the Python value lies inside it. So the two agree. The growth of the element is six percent lower.\n' +
    'JA: 同じ式を解く，あなたの要素と私の Python コードの比較です．要素の応力は Python より 10 % 低いです．ただ，要素は 24³ までなので格子の不確かさが 16.5 % あり，Python の値はその中に入ります．なので両者は合っています．成長は要素のほうが 6 % 低いです．');
}

// 4 code check
{ const s = frame('3. I checked the Python code');
  lines(s, [
    'The quarter domain and the full cube give the same result, to 13 digits.',
    'With a known exact solution, the error falls almost with the square of the cell size.',
    'The point model converged in every one of its 1.9 million steps.',
    'I found no error that changes the results in the thesis.']);
  s.addNotes('EN: I checked my Python code carefully. The quarter domain and the full cube give the same result, to thirteen digits. ' +
    'With a known exact solution, the error falls almost with the square of the cell size. The point model converged in every one of its one point nine million steps. I found no error that changes the results in the thesis.\n' +
    'JA: Python のコードを丁寧に確かめました．4 分の 1 の領域と立方体全体で，結果が 13 桁まで同じです．正解の分かっている解で，誤差はセルの大きさのほぼ 2 乗で減ります．点モデルは 190 万ステップすべてで収束しました．修論の結果を変える誤りは見つかりませんでした．');
}

// 5 questions
{ const s = frame('4. Three questions about your element');
  lines(s, [
    'Question 1. The point model gets the nutrient and the antibiotic one increment late. Is that right?',
    'Question 2. The element integrates the growth from the mean of the two species fields. Is that intended?',
    'Question 3. Does the element limit the time step because of the front term?'], { fontSize: 21, paraSpaceAfter: 22 });
  s.addText('Front term of Eq. 34:  − r · c/(k + c) · n_c · ∇φ₃D', { x: 0.8, y: 6.2, w: 11.8, h: 0.5, fontFace: F, fontSize: 16, color: ACC, isTextBox: true, margin: 0 });
  s.addNotes('EN: I have three questions about your element. ' +
    'Question one. As I understand it, the fields are updated after the material routine, so the point model gets the nutrient and the antibiotic one increment late. Is that right? ' +
    'Question two. The element integrates the growth from the mean of the two species fields. In my thesis I use the growth of my material routine for the stress. Is the mean intended? ' +
    'Question three. Does the element limit the time step because of the front term, minus r c over k plus c, times n c dot grad phi?\n' +
    'JA: 要素について質問が 3 つあります．1 つ目：場は材料ルーチンのあとで更新されるので，点モデルは栄養と抗菌薬を 1 増分遅れて受け取る，という理解で合っていますか．' +
    '2 つ目：要素は 2 つの菌種の場の平均から成長を積分しています．修論では，応力には私の材料ルーチンの成長を使っています．平均を取るのは意図どおりですか．' +
    '3 つ目：前線の項のために，要素は時間の刻みを制限していますか．');
}

// 6 next
{ const s = frame('5. Next steps');
  lines(s, [
    'For the paper, I will put the coupling in LaTeX and send you the parts one by one.',
    'Which part would you like to see first?',
    'I am in Hannover from 28 November to 2 December. Can we meet in person then?']);
  s.addNotes('EN: For the paper, I will put the coupling in LaTeX and send you the parts one by one. Which part would you like to see first? ' +
    'I am in Hannover from the twenty-eighth of November to the second of December. Can we meet in person then? Thank you.\n' +
    'JA: 論文については，連成の部分を LaTeX にまとめて，少しずつ送ります．最初にどの部分を見たいですか．11 月 28 日から 12 月 2 日までハノーファーにいます．そのとき対面で会えますか．ありがとうございました．');
}

// appendix
{ const s = frame('Appendix: notation and equations');
  lines(s, [
    'φ₃D: amount of biofilm (field);  c: nutrient;  a: antibiotic;  α: growth, F_g = α I, growth = α − 1.',
    'Eq. 34:  φ̇₃D = β ∇²φ₃D + k_α α − r · c/(k + c) · n_c · ∇φ₃D,   n_c = ∇c / |∇c|',
    'Eq. 35:  ċ = d ∇²c − g φ₃D (Σᵢ φ̄ᵢ) c;  the antibiotic a follows the same equation.',
    'Eq. 36:  α̇ = k_α φ₃D',
    'Point model (Klempt et al. 2026): φᵢ, ψᵢ (volume fraction and live share of species i), φ̄ᵢ = φᵢ ψᵢ.',
    'Constants: β = 0.02 mm²/T*, r = 10 mm/T*, k_α = 10⁻³ /T*, E = 10 Pa, ν = 0.49.'], { fontSize: 18, paraSpaceAfter: 12 });
  s.addNotes('Appendix for questions: notation and equations of the coupled model.\nJA: 質問に答えるための付録：記号と式．');
}

pres.writeFile({ fileName: B + 'oliver_1125_en.pptx' }).then(f => console.log('wrote', f));
