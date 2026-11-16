// 生きている菌から成長する式（慶應の続き，組成から応力への道）．学術的な日本語のスライド，ノートつき．結果が出るたびに書き足す．
//   NODE_PATH=../2026-11_coupled_examples/build/node_modules node deck_living_ja.js   -> build/living_growth_ja.pptx (not in git)
// 図：assets/fig_growth_living_over_time.png．数値：ansys_usermat/apdl/results/2026-10-3d/GROWTH_LIVING.md．
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const path = require('path');
const R = path.resolve(__dirname, '../..') + '/';
const B = __dirname + '/build/', EQ = B + 'eq/', AS = R + 'assets/';

const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka';
pres.title = '生きている菌から成長する式';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', RULE = '8A99AD', ACC = 'A61E1E', LIGHT = 'F3F5F8';
const JP = 'Yu Mincho', TM = 'Times New Roman';
const FOOT = '西岡 佳祐　｜　生きている菌から成長する式（慶應での続き）';
let n = 0;

function frame(title) {
  const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
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
function fit(s, p, x, y, w, h) {
  if (!fs.existsSync(p)) { cap(s, '（図がまだない：' + path.basename(p) + '）', y + h / 2 - 0.2, x, w); return [w, h]; }
  const [pw, ph] = pngSize(p); let ww = w, hh = w * ph / pw;
  if (hh > h) { hh = h; ww = h * pw / ph; }
  s.addImage({ path: p, x: x + (w - ww) / 2, y: y + (h - hh) / 2, w: ww, h: hh }); return [ww, hh];
}

if (!fs.existsSync(B)) fs.mkdirSync(B, { recursive: true });
const TB = { fontFace: JP, fontSize: 13, color: INK, border: { type: 'solid', color: 'BBBBBB', pt: 0.5 }, align: 'center', valign: 'middle', rowH: 0.4 };
const H = t => ({ text: t, options: { bold: true, fill: { color: LIGHT } } });

// 1 title
{ const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 2.05, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText('生きている菌から成長する式', { x: 0.9, y: 2.25, w: 11.5, h: 0.8, fontFace: JP, fontSize: 32, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addText('抗菌薬が成長と応力に届くようにする（組成から応力への 2 本目の道）', { x: 0.9, y: 3.05, w: 11.5, h: 0.6, fontFace: JP, fontSize: 20, color: INK, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 3.85, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText([{ text: '西岡 佳祐', options: { fontSize: 18, bold: true, breakLine: true } },
    { text: '慶應義塾大学 大学院 理工学研究科 開放環境科学専攻（村松研究室）', options: { fontSize: 14, breakLine: true } },
    { text: '作業中の記録（結果が出るたびに更新）', options: { fontSize: 14 } }],
    { x: 0.9, y: 4.3, w: 11.5, h: 1.5, fontFace: JP, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 6 });
  s.addNotes('慶應での研究の柱の一つ、生きている菌から成長する式についての記録です。今のモデルでは死んだ菌も成長するので、抗菌薬が応力にほとんど効きません。成長を生きている量で決める形に変えると何が起きるかをまとめます。');
}

// 2 why
{ const s = frame('1. 問題：今の式では死んだバイオフィルムも成長する');
  bullets(s, [
    'Klempt et al. (2024) の式 36：α̇ = k_α φ₃D．成長はバイオフィルムの量 φ₃D だけで決まり，菌が生きているかどうかは入らない．',
    '連成モデルでは，抗菌薬で菌が死ぬ（種子の生きている体積分率 Σᵢφ̄ᵢ：0.96 → 0.15，混ぜた例は 0.04）．しかし φ₃D は残るので，成長は変わらない．',
    '結果：抗菌薬による種子の変化は，成長 α − 1 が ±0 %，von Mises 応力が −5 %（分けた例，128³）．抗菌薬は栄養の消費を通してしか力学に届かない．',
    'LUH の修論では限界として書いた点（モデルは凍結）．慶應ではここを変える．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 });
  s.addNotes('今のモデルの成長の式は、バイオフィルムの量だけで決まります。抗菌薬で菌が死んでも量は残るので、成長は変わりません。その結果、抗菌薬を入れても応力は 5 パーセントしか変わりません。修論では限界として書いた点で、慶應ではここを変えます。');
}

// 3 the law
{ const s = frame('2. 式：生きている量で成長を決める');
  s.addTable([
    [H('版'), H('変えるところ'), H('出てくる式')],
    ['A（エネルギーから）', 'エネルギーの結合項 −k_α α φ を\n−k_α α φ Σᵢφ̄ᵢ にする（式 20）', 'α̇ = k_α φ₃D Σᵢφ̄ᵢ（式 36）\n式 34 の生成項も k_α α Σᵢφ̄ᵢ'],
    ['B（式 36 だけ）', 'α の粘性を η_α / Σᵢφ̄ᵢ にする', 'α̇ = k_α φ₃D Σᵢφ̄ᵢ（式 36）\n式 34 はそのまま']],
    Object.assign({}, TB, { x: 0.75, y: 1.25, w: 11.8, colW: [2.6, 4.9, 4.3], rowH: 0.75 }));
  bullets(s, [
    'Σᵢφ̄ᵢ = Σᵢ φᵢψᵢ は点モデルの生きている体積分率．栄養と抗菌薬の消費に使っているのと同じ量で，新しい定数はない（k_α は Klempt 2024 の値のまま）．',
    '活性のある菌だけが成長するという仮定は，バイオフィルムのモデルでは標準（Wanner and Gujer 1986）．',
    'A と B の結果はどの計算でも 0.1 % 以内で一致（式 34 の生成項 k_α α は前線の項 −r c/(k+c) n_c·∇φ₃D に比べて小さい）．導出のある A を本筋にする．',
    '実装：`coupled_living_consumption.py --growth-living-energy`（A），`--growth-living`（B）．'],
    { x: 0.75, y: 3.75, w: 11.8, h: 3.1, fontSize: 15 });
  s.addNotes('式は二つの版を試しました。A はエネルギーの結合項に生きている体積分率を掛ける形で、成長の式と、バイオフィルムの量の式の生成項の両方が変わります。B は成長の式だけを変えます。生きている体積分率は、すでに栄養の消費に使っている量なので、新しい定数は要りません。二つの版の結果は 0.1 パーセント以内で一致したので、導出のある A を本筋にします。');
}

// 4 result table
{ const s = frame('3. 結果：抗菌薬が成長と応力に届く（128³）');
  s.addTable([
    [H('例'), H('成長の式'), H('α − 1 の変化'), H('von Mises の変化')],
    ['分けた例', '式 36：α̇ = k_α φ₃D', '±0 %', '−4.6 %'],
    ['分けた例', '生きている菌から：α̇ = k_α φ₃D Σᵢφ̄ᵢ', '−44 %', '−14 %'],
    ['混ぜた例', '式 36', '−0.4 %', '−9.6 %'],
    ['混ぜた例', '生きている菌から', '−60 %', '−29 %']],
    Object.assign({}, TB, { x: 0.75, y: 1.25, w: 11.8, colW: [2.0, 5.2, 2.3, 2.3], fontSize: 14, rowH: 0.45 }));
  cap(s, '抗菌薬あり（T* = 1/3 から）と抗菌薬なしの比較，種子の平均，T* = 1．2 菌種の合成データ（Fritsch et al. 2025 の case I）．', 3.6, 0.75, 11.8, 11);
  bullets(s, [
    '成長は 4〜6 割下がる．応力の下がりは 14 %（分けた例），29 %（混ぜた例）で，成長より小さい：種子の応力は種子と周りの成長の差で決まり，抗菌薬は両方を遅くするため．',
    '分けた例の応力の変化は格子でまだ動く（32³：−36 %，64³：−21 %，128³：−14 %；式 36 でも −26，−11，−5 %）．2 つの式の差は約 10 ポイントで安定．混ぜた例は −29〜−30 % で安定．',
    '抗菌薬なしでも，同じ k_α では成長が 9〜13 %，応力が 7〜11 % 小さい（生きている体積分率が 0.6 から始まり 1 より小さいため）．'],
    { x: 0.75, y: 4.1, w: 11.8, h: 2.8, fontSize: 14 });
  s.addNotes('結果です。生きている菌から成長する式にすると、抗菌薬で種子の成長が 4 割から 6 割下がります。応力の下がりは分けた例で 14 パーセント、混ぜた例で 29 パーセントで、成長より小さいです。種子の応力は、種子と周りの成長の差で決まり、抗菌薬は両方を遅くするからです。分けた例の応力の変化は格子でまだ動きますが、二つの式の差は約 10 ポイントで安定しています。');
}

// 5 time history
{ const s = frame('4. 時間変化（64³）');
  fit(s, AS + 'fig_growth_living_over_time.png', 0.55, 1.1, 7.4, 5.75);
  bullets(s, [
    '上段：生きている体積分率 Σᵢφ̄ᵢ は，抗菌薬（T* = 1/3，点線）のあと 0.15（分けた例），0.04（混ぜた例）まで下がる．2 つの式で同じ．',
    '中段：式 36（灰色の破線）では α − 1 が同じ速さで増え続ける．生きている菌からの式（赤）では，分けた例でゆるやかに頭打ち，混ぜた例では T* ≈ 0.4 でほぼ止まる．',
    '下段：応力は抗菌薬の前から差がつき（生きている体積分率が 1 より小さいため），そのあともほぼ一定の差で進む．'],
    { x: 8.15, y: 1.25, w: 4.6, h: 5.5, fontSize: 13.5 });
  s.addNotes('時間変化です。上の段は生きている体積分率で、抗菌薬のあと大きく下がります。中の段は成長です。今の式では抗菌薬のあとも同じ速さで増え続けますが、生きている菌からの式では、混ぜた例でほぼ止まります。下の段は応力で、抗菌薬の前から差がつき、そのあともほぼ一定の差で進みます。');
}

// 6 reading and next
{ const s = frame('5. 読み方と次にやること');
  bullets(s, [
    '言えること：成長を生きている量で決めると，抗菌薬が成長（−44〜−60 %）と応力（−14〜−29 %）に届く．新しい定数は要らず，エネルギーの枠組みの中で導ける．',
    '言えないこと：効果の大きさが正しいかどうか．実験との照合はしていない．分けた例の応力の変化は格子に依存する．',
    '照合に使えそうなデータ：クロルヘキシジンでの菌種ごとの生き残り（Gränicher et al. 2021 は 1 濃度・1 時点；時間変化のある Schwendener et al. 2025 が向く可能性，本文は未確認）．',
    '次：(1) 歯頸部の窓の例で同じ比較，(2) 組成で変わる硬さ（死んだバイオフィルムを柔らかくする）と組み合わせる，(3) 最終の成長をそろえる k_α の決め方．',
    '位置づけ：「組成から応力へ」の 2 本目の道（1 本目は栄養の消費，3 本目は硬さ）．変形との連成（式 30 の項）とは別の軸．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.5, fontSize: 15 });
  s.addNotes('まとめです。成長を生きている量で決めると、抗菌薬が成長と応力に届きます。新しい定数は要らず、エネルギーの枠組みの中で導けます。ただし、効果の大きさが正しいかどうかは、実験と照合していないので言えません。次は、歯頸部の窓の例で同じ比較をし、組成で変わる硬さと組み合わせます。');
}

// appendix
{ const s = frame('付録：記号と設定');
  bullets(s, [
    'φ₃D：バイオフィルムの体積分率（場），α：成長（F_g = α I，α(0) = 1，成長は α − 1），k_α：成長の速さ（Klempt 2024：10⁻³ /T*）．',
    'φᵢ，ψᵢ：点モデルの菌種 i の体積分率と生きている割合，φ̄ᵢ = φᵢψᵢ，Σᵢφ̄ᵢ：生きている体積分率．',
    '設定：一辺 2 mm の立方体，2 菌種（Fritsch et al. 2025 case I），抗菌薬は T* = 1/3 から α*₀ = 100，栄養の式は ċ つき，minmod，外側の刻み 0.0025（表），0.025（時間変化の図）．',
    '種子の平均：始めにバイオフィルムを置いた領域（十字形，32 要素）の平均．応力は 16³ の六面体要素で解く（E = 10 Pa，ν = 0.49）．',
    '数値の一覧：`ansys_usermat/apdl/results/2026-10-3d/GROWTH_LIVING.md`．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.5, fontSize: 15 });
  s.addNotes('付録です。記号と計算の設定です。');
}

pres.writeFile({ fileName: B + 'living_growth_ja.pptx' }).then(f => console.log('wrote', f));
