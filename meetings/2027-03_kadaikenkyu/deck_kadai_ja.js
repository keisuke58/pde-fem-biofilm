// 課題研究（慶應）の 5 分の発表動画用スライド（日本語，発表者ノートに原稿）．
//   cd meetings/2027-03_kadaikenkyu && node deck_kadai_ja.js   -> build/kadaikenkyu_ja.pptx (not in git)
// 図：keio_report/ の図（EML の PDF は build/ で PNG にする：pdftoppm -r 220 -png -singlefile），assets/，
// 式：meetings/2026-11_coupled_examples/build/eq/eqs.png（その build.sh で作る）．
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const path = require('path');
const R = path.resolve(__dirname, '../..') + '/';
const B = __dirname + '/build/';
const KR = R + 'keio_report/';
const EQ = R + 'meetings/2026-11_coupled_examples/build/eq/';

const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka';
pres.title = '課題研究：多菌種バイオフィルムモデルのベイズ推定と連続体成長モデルとの連成';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', RULE = '8A99AD', ACC = 'A61E1E', LIGHT = 'F3F5F8';
const JP = 'Yu Mincho', TM = 'Times New Roman';
let n = 0;

function frame(title) {
  const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addText(title, { x: 0.55, y: 0.32, w: 12.2, h: 0.6, fontFace: JP, fontSize: 24, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 0.98, w: 12.23, h: 0, line: { color: NAVY, width: 1.25 } });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 7.0, w: 12.23, h: 0, line: { color: RULE, width: 0.5 } });
  s.addText('西岡 佳祐　｜　課題研究：多菌種バイオフィルムモデルのベイズ推定と連続体成長モデルとの連成', { x: 0.55, y: 7.05, w: 10.5, h: 0.3, fontFace: JP, fontSize: 10, color: GREY, isTextBox: true, margin: 0 });
  s.addText(String(n), { x: 11.78, y: 7.05, w: 1.0, h: 0.3, fontFace: TM, fontSize: 11, color: GREY, align: 'right', isTextBox: true, margin: 0 });
  return s;
}
const cap = (s, t, y, x = 0.55, w = 12.23, fs_ = 11) => s.addText(t, { x, y, w, h: 0.5, fontFace: JP, fontSize: fs_, color: GREY, isTextBox: true, margin: 0, valign: 'top' });
const bullets = (s, items, o) => s.addText(items.map((t, i) => ({ text: t, options: { bullet: { indent: 18 }, breakLine: i < items.length - 1 } })),
  Object.assign({ fontFace: JP, fontSize: 16, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 10 }, o));
function pngSize(p) { const b = fs.readFileSync(p); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }
// picture fitted into the box (x, y, w, h), centred
function fit(s, p, x, y, w, h) {
  const [pw, ph] = pngSize(p); let ww = w, hh = w * ph / pw;
  if (hh > h) { hh = h; ww = h * pw / ph; }
  s.addImage({ path: p, x: x + (w - ww) / 2, y: y + (h - hh) / 2, w: ww, h: hh }); return [ww, hh];
}
function stat(s, x, y, big, small) {
  s.addShape(pres.shapes.RECTANGLE, { x, y, w: 3.7, h: 1.5, fill: { color: LIGHT }, line: { color: LIGHT } });
  s.addText(big, { x, y: y + 0.12, w: 3.7, h: 0.8, fontFace: TM, fontSize: 40, bold: true, color: NAVY, align: 'center', isTextBox: true, margin: 0 });
  s.addText(small, { x: x + 0.15, y: y + 0.92, w: 3.4, h: 0.5, fontFace: JP, fontSize: 13, color: INK, align: 'center', isTextBox: true, margin: 0 });
}

// 1 title
{ const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 1.85, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText([{ text: '拡張ハミルトン原理に基づく多菌種バイオフィルムモデルの', options: { breakLine: true } },
    { text: 'ベイズ推定と連続体成長モデルとの連成' }],
    { x: 0.9, y: 2.05, w: 11.5, h: 1.3, fontFace: JP, fontSize: 30, bold: true, color: NAVY, isTextBox: true, margin: 0, valign: 'top' });
  s.addText('Bayesian Calibration of an Extended-Hamilton-Principle Multi-Species Biofilm Model and its Coupling to a Continuum Growth Model',
    { x: 0.9, y: 3.4, w: 11.5, h: 0.6, fontFace: TM, fontSize: 15, italic: true, color: INK, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 4.1, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText([{ text: '西岡 佳祐（82519093）', options: { fontSize: 18, bold: true, breakLine: true } },
    { text: '慶應義塾大学 大学院 理工学研究科 開放環境科学専攻 応用力学・計算力学専修　指導：村松 眞由 准教授', options: { fontSize: 14, breakLine: true } },
    { text: 'ハノーファー大学（LUH）連続体力学研究所でのダブルディグリーの修士研究に基づく', options: { fontSize: 14, breakLine: true } },
    { text: '課題研究　2027 年 3 月', options: { fontSize: 14 } }],
    { x: 0.9, y: 4.45, w: 11.5, h: 1.9, fontFace: JP, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 6 });
  s.addNotes('慶應義塾大学 開放環境科学専攻の西岡佳祐です。課題研究として、拡張ハミルトン原理に基づく多菌種バイオフィルムモデルのベイズ推定と、連続体の成長モデルとの連成について報告します。この研究は、ハノーファー大学の連続体力学研究所でのダブルディグリーの修士研究に基づいています。');
}

// 2 background
{ const s = frame('1. 背景と目的');
  const box = (x, head, body, col) => {
    s.addShape(pres.shapes.RECTANGLE, { x, y: 1.35, w: 5.6, h: 2.3, fill: { color: LIGHT }, line: { color: col, width: 1.5 } });
    s.addText(head, { x: x + 0.25, y: 1.5, w: 5.1, h: 0.5, fontFace: JP, fontSize: 18, bold: true, color: col, isTextBox: true, margin: 0 });
    s.addText(body, { x: x + 0.25, y: 2.05, w: 5.1, h: 1.5, fontFace: JP, fontSize: 14, color: INK, isTextBox: true, margin: 0, valign: 'top' });
  };
  box(0.55, '点モデル（Klempt et al. 2026）', '1 点での菌種の組成 φᵢ と生存率 ψᵢ の時間変化．相互作用行列 A で菌種がつながる．空間を持たない．', NAVY);
  box(7.18, '成長の場のモデル（Klempt et al. 2024）', 'バイオフィルムの量 φ₃D，栄養 c，成長 α の場と応力．菌種を区別しない．', ACC);
  s.addText('どちらも拡張ハミルトン原理から導かれる', { x: 0.55, y: 3.78, w: 12.23, h: 0.4, fontFace: JP, fontSize: 15, color: GREY, align: 'center', isTextBox: true, margin: 0 });
  bullets(s, [
    '実験（Heine et al. 2025）：口腔内 5 菌種，共生・疾患型 × 静置・流れの 4 条件，21 日間の組成．菌種間の促進と競合の強さは実験だけでは決まらない．',
    '目的 1：5 菌種の点モデルの相互作用 15 成分を実測データから推定し，データが決める成分を示す．',
    '目的 2：点モデルを有限要素の各積分点で動かし，栄養と抗菌薬の場と連成して，組成が成長と応力に及ぼす影響を調べる．'],
    { x: 0.75, y: 4.35, w: 11.8, h: 2.5, fontSize: 15 });
  s.addNotes('バイオフィルムのモデルには二つの流れがあります。点モデルは、一点での菌種の組成と生存率の変化を表しますが、空間を持ちません。成長の場のモデルは、バイオフィルムがどこで成長し、どんな応力を持つかを表しますが、菌種を区別しません。どちらも拡張ハミルトン原理から導かれます。本研究の目的は二つです。一つ目は、Heine らの実験データから、点モデルの相互作用を推定すること。二つ目は、点モデルを有限要素の中で動かして場と連成し、組成が成長と応力にどう効くかを調べることです。');
}

// 3 estimation method
{ const s = frame('2. 相互作用の推定：TMCMC と GPU 計算');
  stat(s, 0.55, 1.3, '15', '推定する成分（対称な A）');
  stat(s, 4.82, 1.3, '4 × 3', '段階 × 独立なシード');
  stat(s, 9.08, 1.3, '147 倍', 'GPU による高速化（CPU 12 コア比）');
  bullets(s, [
    'データ：生きている菌の組成と全菌の組成の 2 チャネル，群集の生存率．重み付きのガウス尤度．pH は尤度に入れず確認にだけ使う．',
    'A は自由エネルギーの二次形式から対称になる：非対称なら 25 成分で，各条件の観測 25 個に対して識別できない．',
    '遷移マルコフ連鎖モンテカルロ法（TMCMC）．事前分布の箱を狭めながら 4 段階で推定し，各段を 3 シードで繰り返して 5 つの基準を満たす段だけを採用．',
    '各時間ステップの陰解法の Newton 法を JAX で 1 つのカーネルにし，全粒子を GPU で並列に評価（1 条件 21,600 s → 147 s）．'],
    { x: 0.75, y: 3.15, w: 11.8, h: 3.7, fontSize: 15 });
  s.addNotes('推定するのは、対称な相互作用行列の 15 成分です。対称性は自由エネルギーの形から自然に出てきます。非対称だと 25 成分になり、データの数に対して決まりません。方法は遷移マルコフ連鎖モンテカルロ法で、事前分布の範囲を狭めながら 4 段階で推定し、各段を 3 つの乱数の種で繰り返して、結果がそろった段だけを使いました。計算はモデルの Newton 法を JAX で GPU 向けにまとめ、CPU より約 150 倍速くしたことで、この繰り返しが現実的になりました。');
}

// 4 estimation results
{ const s = frame('3. 推定結果：データが決める相互作用');
  fit(s, B + 'eml_ident.png', 0.55, 1.15, 7.6, 5.0);
  bullets(s, [
    '識別できた成分：CS 5，CH 6，DS 9，DH 8 個（15 個のうち）．',
    '共生条件では初期定着菌（So，An，Vei）の間だけが決まり，So–Vei は正（乳酸の交差給餌と整合）．',
    '疾患型では後期定着菌の成分も決まる．Fn–Pg は識別されるが符号は決まらない．',
    '推定に使っていない pH：DH の上昇を R² = 0.71 で追った．'],
    { x: 8.4, y: 1.3, w: 4.4, h: 4.9, fontSize: 14 });
  cap(s, '最終段，3 シードを合わせた事後中央値（色と数字），斜線：符号のみ，灰色：識別できない．CS/CH：共生・静置/流れ，DS/DH：疾患型・静置/流れ．'
    + '本節の内容は投稿済みの論文に基づく（Extreme Mechanics Letters，査読中，EML-D-26-00377）．', 6.25);
  s.addNotes('推定結果です。色のついたマスが、データで決まった成分です。15 個のうち、条件によって 5 個から 9 個が決まりました。共生条件では、初期に定着する 3 菌種の間だけが決まり、ストレプトコッカスとベイロネラの間は正、つまり協力で、乳酸を受け渡す関係と合っています。疾患型の条件では、後期に定着する菌の成分も決まりました。推定に使っていない pH の変化も、疾患型の流れ条件で決定係数 0.71 で追えています。この部分は投稿済みの論文に基づいています。');
}

// 5 coupling
{ const s = frame('4. 連成：点モデルを有限要素の各積分点で動かす');
  fit(s, KR + 'fig_coupling_schematic_agreed.png', 0.55, 1.15, 12.23, 4.15);
  bullets(s, [
    '成長 F_g = α I を ANSYS の材料ルーチンとして実装し，閉じた形の解で検証して，LUH のバイオフィルム要素に組み込んだ．',
    '点モデルは公表された式のまま各積分点で動かす．栄養と抗菌薬の場から c*，α* を受け取り，生存体積分率 Σᵢφ̄ᵢ を返す．',
    '栄養と抗菌薬は生きている菌だけが消費する：ċ = d∇²c − g φ₃D (Σᵢφ̄ᵢ) c．組成は栄養と前線の項 −r c/(k+c) n_c·∇φ₃D を通して成長と応力に届く．'],
    { x: 0.75, y: 5.4, w: 11.8, h: 1.5, fontSize: 13, paraSpaceAfter: 4 });
  s.addNotes('次に連成です。まず成長の式を ANSYS の材料ルーチンとして実装し、厳密解で検証してから、研究所のバイオフィルム要素に組み込みました。そのうえで、点モデルを各積分点で、公表された式のまま動かします。点モデルは場から栄養と抗菌薬の強さを受け取り、生きている量を返します。栄養と抗菌薬を消費するのは生きている菌だけです。こうして菌の組成は、栄養と、バイオフィルムを栄養の方へ広げる前線の項を通して、成長と応力に届きます。');
}

// 6 two species with antibiotic
{ const s = frame('5. 抗菌薬を含む 2 菌種の例');
  fit(s, KR + 'fig_coupled_3d_panels_separate_ab.png', 0.55, 1.1, 12.23, 3.75);
  cap(s, '分けて置いた例（菌 A を種子に，菌 B を 4 辺の小領域に），抗菌薬は T* = 1/3 から外側の面に，64³．上：T* = 0.325，下：T* = 1．(a) A の割合，(b) 生存体積分率 Σᵢφ̄ᵢ，(c) 栄養 c，(d) 抗菌薬 a，(e) von Mises 応力．', 4.9);
  stat(s, 0.55, 5.4, '0.96 → 0.15', '種子の生存体積分率（分けて置く）');
  stat(s, 4.82, 5.4, '0.04', '混ぜて置いたとき');
  stat(s, 9.08, 5.4, '≈ 1 %', '抗菌薬による成長 α − 1 の変化');
  s.addNotes('合成データの 2 菌種で、抗菌薬を途中から外側の面に与えた例です。種子の生きている量は、菌種を分けて置くと 0.96 から 0.15 に、混ぜて置くと 0.04 まで下がります。同じ菌でも、相手と同じ点にいると早く死にます。一方、成長はほとんど変わりません。成長の式は死んだ菌も含むバイオフィルムの量で決まるので、抗菌薬は生きている割合を下げても、量は下げないからです。');
}

// 7 numerics and geometry
{ const s = frame('6. 数値の信頼性と形状');
  const H = { fontFace: JP, fontSize: 13, bold: true, color: 'FFFFFF', fill: { color: NAVY }, align: 'center', valign: 'middle' };
  const C = t => ({ text: t, options: { fontFace: TM, fontSize: 13, color: INK, align: 'center', valign: 'middle' } });
  const L = t => ({ text: t, options: { fontFace: JP, fontSize: 13, color: INK, align: 'left', valign: 'middle' } });
  s.addTable([
    [{ text: '前線の項の制限関数', options: H }, { text: 'minmod', options: H }, { text: 'van Leer', options: H }, { text: 'MC', options: H }, { text: 'superbee', options: H }],
    [L('混ぜた例の応力の GCI（64/128/192³）'), C('41 %'), C('22 %'), C('31 %'), { text: '4.6 %', options: { fontFace: TM, fontSize: 13, bold: true, color: ACC, align: 'center' } }]],
    { x: 0.55, y: 1.3, w: 7.3, colW: [2.9, 1.1, 1.1, 1.1, 1.1], rowH: 0.45, border: { type: 'solid', color: RULE, pt: 0.5 } });
  bullets(s, [
    '前線の項はセルのペクレ数約 31 の移流：2 次の風上差分が要る．数値拡散の小さい superbee だけが漸近域に入る（128–256³ で 3.95 × 10⁻⁴ Pa，GCI 6.3 %）．',
    '前線の近くだけ細かくする 2 段の格子：一様な 256³ と 0.1 % で一致．',
    '剛な 4 本の柱（Klempt 2024）：種子の応力は約 8 倍（分けた例）．superbee で GCI 4.5 %．',
    '菌種 B の置き方：組成はほとんど変わらないが，応力は 17–27 % 変わる（128³）．'],
    { x: 0.75, y: 2.45, w: 7.0, h: 4.4, fontSize: 14 });
  fit(s, R + 'thesis/figures/fig_asymmetric_starts.png', 8.05, 1.15, 4.73, 5.75);
  s.addNotes('数値の信頼性です。前線の項は拡散に比べて移流が強く、2 次の風上差分が必要です。傾きの制限関数を 4 つ比べると、混ぜた例の応力は数値拡散の小さい superbee だけが収束しました。前線の近くだけを細かくする格子は、一様な細かい格子と 0.1 % で一致します。形状の例では、剛な柱があると種子の応力が約 8 倍になります。右の図は菌種 B の置き方を変えた例で、組成はほとんど同じでも応力は 2 割前後変わります。');
}

// 8 summary
{ const s = frame('7. まとめと今後');
  bullets(s, [
    '5 菌種の点モデルの相互作用 15 成分を GPU-TMCMC で推定し，データが決める成分を示した（5–9 成分，GPU で約 150 倍）．',
    '成長の運動学を ANSYS の材料ルーチンとして実装・検証し，点モデルを各積分点で動かして栄養と抗菌薬の場と連成した．',
    '抗菌薬は生存量を大きく下げるが，成長と応力はほとんど変えない．応力は GCI 1.4 %（分けた例），混ぜた例は superbee の制限関数で 6.3 % で求まる．',
    '限界：力学は一方向の連成．検証（verification）は行ったが，実験との比較（validation）はこれから．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 3.0, fontSize: 16 });
  s.addShape(pres.shapes.RECTANGLE, { x: 0.55, y: 4.55, w: 12.23, h: 1.75, fill: { color: LIGHT }, line: { color: LIGHT } });
  s.addText([{ text: '今後（慶應での修士研究）', options: { bold: true, color: NAVY, breakLine: true } },
    { text: '組成から応力へ：生きている菌から成長する式，菌種や生死で変わる硬さ，実際の歯・インプラントの形での応力．', options: { breakLine: true } },
    { text: '連成の内容を論文にまとめて投稿する（2027 年）．' }],
    { x: 0.8, y: 4.65, w: 11.8, h: 1.6, fontFace: JP, fontSize: 15, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 4 });
  s.addNotes('まとめます。点モデルの相互作用をデータから推定し、どの成分がデータで決まるかを示しました。次に、成長の式を有限要素に組み込み、点モデルを各積分点で動かして場と連成しました。抗菌薬は生きている量を大きく下げますが、成長と応力はほとんど変えません。今後は慶應での修士研究として、組成が応力に届く道、つまり生きている菌から成長する式や、菌種で変わる硬さ、実際の歯の形での応力に取り組みます。以上です。');
}

// appendix: notation
{ const s = frame('付録：記号と支配方程式');
  const p = EQ + 'eqs.png';
  if (fs.existsSync(p)) fit(s, p, 0.9, 1.2, 11.5, 3.9);
  cap(s, 'φ₃D：バイオフィルムの体積分率，α：成長（F_g = α I，α(0) = 1），c：栄養，a：抗菌薬，φᵢ・ψᵢ：点モデルの体積分率と生存率，'
    + 'φ̄ᵢ = φᵢψᵢ，Σᵢφ̄ᵢ：生存体積分率，A：相互作用行列，c*・α*：点モデルが受け取る栄養と抗菌薬の強さ．式番号は Klempt et al. 2024．'
    + '式 35 の消費を生存体積分率で重み付けしたのが本研究の連成で，抗菌薬も同じ形とする．', 5.3, 0.55, 12.23, 12);
  s.addNotes('付録です。質問に答えるための記号と式の一覧で、発表では使いません。');
}

pres.writeFile({ fileName: B + 'kadaikenkyu_ja.pptx' }).then(f => console.log('wrote', f));
