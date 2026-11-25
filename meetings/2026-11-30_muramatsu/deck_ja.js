// 村松先生との打ち合わせ（30 Nov 2026, 17:00 JST）：前回（2 Nov）にいただいた 2 つの質問への回答．日本語，ノートつき．
//   NODE_PATH=../2026-11_coupled_examples/build/node_modules node deck_ja.js   -> build/muramatsu_1130_ja.pptx (not in git)
// 図：assets/fig_uncut_symmetry.png（fig_uncut_symmetry.py），assets/fig_point_model_survival.png（fig_point_model_survival.py），
// 式：../2026-11_coupled_examples/build/eq/point.png（同じデッキの build.sh で作る）．
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const path = require('path');
const R = path.resolve(__dirname, '../..') + '/';
const B = __dirname + '/build/', AS = R + 'assets/', EQ = R + 'meetings/2026-11_coupled_examples/build/eq/';

const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka';
pres.title = '前回の質問への回答（2 菌種の連成モデル）';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', RULE = '8A99AD', ACC = 'A61E1E', LIGHT = 'F3F5F8';
const JP = 'Yu Mincho', TM = 'Times New Roman';
const FOOT = '西岡 佳祐　｜　前回（11 月 2 日）の質問への回答';
let n = 0;
const SL = [];

function frame(title) {
  const s = pres.addSlide(); n += 1; SL.push(s); s.background = { color: 'FFFFFF' };
  s.addText(title, { x: 0.55, y: 0.32, w: 12.2, h: 0.6, fontFace: JP, fontSize: 24, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 0.98, w: 12.23, h: 0, line: { color: NAVY, width: 1.25 } });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 7.0, w: 12.23, h: 0, line: { color: RULE, width: 0.5 } });
  s.addText(FOOT, { x: 0.55, y: 7.05, w: 10.5, h: 0.3, fontFace: JP, fontSize: 10, color: GREY, isTextBox: true, margin: 0 });
  s.addText(String(n), { x: 11.78, y: 7.05, w: 1.0, h: 0.3, fontFace: TM, fontSize: 11, color: GREY, align: 'right', isTextBox: true, margin: 0 });
  return s;
}
const cap = (s, t, y, x = 0.55, w = 12.23, fs_ = 11) => s.addText(t, { x, y, w, h: 0.5, fontFace: JP, fontSize: fs_, color: GREY, isTextBox: true, margin: 0, valign: 'top' });
const bullets = (s, items, o) => s.addText(items.map((t, i) => ({ text: t, options: { bullet: { indent: 18 }, breakLine: i < items.length - 1 } })),
  Object.assign({ fontFace: JP, fontSize: 16, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 10 }, o));
function pngSize(p) { const b = fs.readFileSync(p); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
const img = (s, p, x, y, w) => { const [pw, ph] = pngSize(p); s.addImage({ path: p, x, y, w, h: w * ph / pw }); return w * ph / pw; };
if (!fs.existsSync(B)) fs.mkdirSync(B, { recursive: true });
const NOTES = [];

// 1 symmetry
{ const s = frame('1. 応力は x と z について鏡映対称です');
  const h = img(s, AS + 'fig_uncut_symmetry.png', 0.55, 1.2, 8.3);
  cap(s, '図 1　分離配置・抗菌薬あり，T* = 1 の von Mises 応力（16³ の応力格子）．(a) 切らずに全体を通常の視点から，(b) 上（+y）から，(c) y = 0 の断面と鏡映面 x = 0，z = 0（青の破線）．', 1.2 + h + 0.08, 0.55, 8.3);
  bullets(s, [
    '設定（種子と 4 本の B の柱，栄養の層，抗菌薬の殻）は x と z について鏡映対称で，解も対称です．鏡映した要素どうしの差は，計算の丸め誤差の大きさでした．',
    '前回の図が非対称に見えたのは，中を見せるために 1/8（x, y, z > 0）を切り取り，斜めから見ていたためです．切った側は内部が，反対側は外面が見えます．',
    'y 方向には対称ではありません．栄養が下の層からだけ入るためです．'],
    { x: 9.15, y: 1.25, w: 3.65, h: 5.6, fontSize: 15 });
  NOTES.push('前回，応力が非対称に見えるという質問をいただきました．設定は x と z について鏡映対称で，解も対称です．切らずに全体を見ると (a)，上から見ると (b)，y = 0 の断面に鏡映面を引くと (c) のとおり左右対称で，鏡映した要素どうしの差は丸め誤差の大きさでした．前回の図は中を見せるために 1/8 を切り取って斜めから見ていたので，切った側だけ内部が見えて非対称に見えました．なお y 方向は，栄養が下からだけ入るので対称ではありません．'); }

// 2 antibiotic term
{ const s = frame('2. 点モデルの抗菌薬の項（Klempt et al. 2026, 式 16〜18）');
  const h = img(s, EQ + 'point.png', 0.55, 1.15, 6.4);
  cap(s, 'φᵢ：菌種 i の体積分率，ψᵢ：生きている割合，φ̄ᵢ = φᵢψᵢ，φ₀：空き，γ：拘束の乗数，A：相互作用，bᵢ：抗菌薬への感度，ηᵢ：粘性．赤：抗菌薬の項．最後の式は，菌種 i が点をほぼ占めるときの ψᵢ の変化率．', 1.15 + h + 0.05, 0.55, 6.4);
  const h2 = img(s, AS + 'fig_point_model_survival.png', 7.15, 1.15, 5.65);
  cap(s, '図 2　点モデルだけ（Fritsch case I）．(a) 抗菌薬 α* = 100（T* = 1/3 から）での生きている量，(b) 抗菌薬なしでの A の割合．', 1.15 + h2 + 0.05, 7.15, 5.65);
  bullets(s, [
    '抗菌薬は式 17 の −α* bᵢ ψᵢ で生きている割合 ψᵢ だけを下げます．菌の量 φᵢ を減らす項はないので，死んだ菌は場所を占めたまま残ります．',
    'ψᵢ の減り方は bᵢ/(ηᵢφᵢ²) に比例します．感度 b の大きい種（B，b = 2）と，点の中で少ない種ほど速く死にます．',
    '栄養の項 c*Aᵢᵢ/ηᵢ がこれに逆らうので，生き残るのは「A が多く，栄養も多い点」だけです（図 2 (a)）．連成計算で栄養の入口に近い A の領域だけが生き残るのは，このためです．'],
    { x: 0.75, y: 4.35, w: 12.0, h: 2.55, fontSize: 14 });
  NOTES.push('二つ目は点モデルの抗菌薬の項です．抗菌薬は式 17 の赤い項，マイナス α* bᵢ ψᵢ で入り，生きている割合 ψ だけを下げます．菌の量 φ を減らす項はないので，死んだ菌は場所を占めたまま残ります．最後の式のとおり，ψ の減り方は b/(η φ²) に比例するので，感度 b の大きい B と，点の中で少ない種ほど速く死にます．栄養の項がこれに逆らうので，右の図のように，A が多く栄養も多い点だけが生き残ります．連成の計算で，栄養の入口に近い A の領域だけが生き残るのはこのためです．'); }

// appendix: notation
{ const s = frame('付録：記号');
  const rows = [['記号', '意味', '記号', '意味'],
    ['φ₃D', 'バイオフィルムの体積分率（場，Klempt 2024）', 'φᵢ', '菌種 i の体積分率（点モデル）'],
    ['c', '栄養（場）', 'ψᵢ', '菌種 i の生きている割合'],
    ['a', '抗菌薬（場）', 'φ̄ᵢ = φᵢψᵢ', '生きている体積分率'],
    ['α', '成長，F_g = αI（式 36）', 'c*，α*', '点モデルが受け取る栄養と抗菌薬（c* = c*₀c，α* = α*₀a）'],
    ['T*', '正規化した時間（点モデルのデータの期間 = 1）', 'A，bᵢ，ηᵢ', '相互作用，抗菌薬への感度，粘性']];
  s.addTable(rows.map((r, i) => r.map(t => ({ text: t, options: i === 0 ? { bold: true, color: 'FFFFFF', fill: { color: NAVY } } : {} }))),
    { x: 0.55, y: 1.3, w: 12.23, colW: [1.4, 4.7, 1.6, 4.53], fontFace: JP, fontSize: 14, color: INK, rowH: 0.55, border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  NOTES.push('記号の一覧です．質問があったときに使います．'); }

if (NOTES.length !== SL.length) throw new Error('notes ' + NOTES.length + ' for ' + SL.length + ' slides');
SL.forEach((sl, i) => sl.addNotes(NOTES[i]));
pres.writeFile({ fileName: B + 'muramatsu_1130_ja.pptx' }).then(f => console.log('ok', f, n, 'slides'));
