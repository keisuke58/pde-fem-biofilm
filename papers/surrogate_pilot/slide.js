// One slide (Japanese): where machine learning helps in this model, with the pilot. node slide.js -> build/surrogate_slide.pptx
const pptxgen = require('pptxgenjs'); const fs = require('fs'); const path = require('path');
const R = path.resolve(__dirname, '../..') + '/'; const OUT = __dirname + '/build/';
fs.mkdirSync(OUT, { recursive: true });
const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka'; pres.title = '代理モデルの使いどころ';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', LIGHT = 'F3F5F8', JP = 'Yu Mincho';
const s = pres.addSlide(); s.background = { color: 'FFFFFF' };
s.addText('機械学習の使いどころ：不確かさを応力まで運ぶ代理モデル', { x: 0.5, y: 0.3, w: 12.3, h: 0.6, fontFace: JP, fontSize: 24, bold: true, color: NAVY, isTextBox: true, margin: 0 });
s.addShape(pres.shapes.LINE, { x: 0.5, y: 0.95, w: 12.33, h: 0, line: { color: NAVY, width: 1.25 } });
const box = (x, y, w, h, title, items) => {
  s.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: LIGHT }, line: { color: 'D5DCE5', width: 0.5 } });
  s.addText(title, { x: x + 0.15, y: y + 0.1, w: w - 0.3, h: 0.35, fontFace: JP, fontSize: 15, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addText(items.map((t, i) => ({ text: t, options: { bullet: { indent: 14 }, breakLine: i < items.length - 1 } })),
    { x: x + 0.15, y: y + 0.5, w: w - 0.3, h: h - 0.6, fontFace: JP, fontSize: 12, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 4 }); };
box(0.5, 1.2, 5.6, 2.15, '使う：パラメータ → 応力の代理モデル', [
  '事後分布を応力まで運ぶには数百〜数千本の計算が要る（1 本：64³ 約 5 分，192³ 約 3 時間）',
  '数百本でガウス過程（または小さな NN）を学習し，応力と生存体積分率を即座に返す（応力は 64³ 以上，生存体積分率は 32³ で可）',
  'Fritsch et al. 2025（代理モデルでベイズ更新）と同じ流れ']);
box(0.5, 3.5, 5.6, 1.75, '使わない', [
  'PINN で場を解く：前線の項 −r c/(k+c) n_c·∇φ₃D で縁が急（ペクレ数大）．有限体積法に勝てない',
  '点モデルの置き換え：Newton で十分速い．時間の大半は場の計算',
  '式そのもの：モデルは固定']);
box(0.5, 5.4, 5.6, 1.45, '順番', [
  '12 月まで修論とコロキウム → 1 月から慶應で小さく試す（32³ で 200 本なら 1 日）',
  '入力を 3 章の事後標本に替え，5 菌種の点モデルで → 連成論文の次の論文']);
const im = R + 'assets/fig_surrogate_pilot.png'; const b = fs.readFileSync(im); const pw = b.readUInt32BE(16), ph = b.readUInt32BE(20);
const w = 6.6, h = w * ph / pw;
s.addText('試し（32³，分離・抗菌薬あり，16 本，約 10 分）', { x: 6.35, y: 1.2, w: 6.5, h: 0.35, fontFace: JP, fontSize: 15, bold: true, color: NAVY, isTextBox: true, margin: 0 });
s.addImage({ path: im, x: 6.3, y: 1.65, w, h });
s.addText([
  { text: 'α*₀（10〜300）と gℓ²/d（0.5〜2）をラテン超方格で 16 本．(a) ガウス過程による生存体積分率と計算点，(b)(c) 1 つ抜きの予測（±2σ）．', options: { breakLine: true } },
  { text: '生存体積分率は α*₀ ≈ 60〜100 で急に切り替わり，16 本では足りない（切り替わり付近に点を足す，ロジット変換）．32³ でも 64³ と 0.001 以内．', options: { breakLine: true } },
  { text: 'わかったこと：32³ の応力は学習に使えない．抗菌薬による応力の差が 32³ では −26 %，収束値では −2 % で，ばらつきの大部分が格子の誤差．応力は 64³ 以上，または 32³ と 64³ の組み合わせ（multi-fidelity）で学習する．', options: { bold: true } }],
  { x: 6.35, y: 1.75 + h, w: 6.45, h: 6.85 - (1.75 + h), fontFace: JP, fontSize: 12, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 5 });
s.addText('西岡 佳祐　｜　慶應での続き（修論には入れない）', { x: 0.5, y: 7.05, w: 10, h: 0.3, fontFace: JP, fontSize: 10, color: GREY, isTextBox: true, margin: 0 });
pres.writeFile({ fileName: OUT + 'surrogate_slide.pptx' }).then(f => console.log('ok', f));
