// Japanese report deck for one coupled example (with speaker notes). Built by build.sh: node deck_ja.js separate|mixed
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const CASE = process.argv[2] || 'separate';
const MIX = CASE === 'mixed';
const path = require('path');
const R = path.resolve(__dirname, '../..') + '/';
const F = __dirname + '/build/', EQ = F + 'eq/';
const VIDEO = R + 'assets/anim_coupled_tvd_' + (MIX ? 'mixed_ab' : 'separate_ab') + '_inhliv_transc_5panel.mp4';
const POSTER = R + 'assets/anim_coupled_tvd_' + (MIX ? 'mixed_ab' : 'separate_ab') + '_inhliv_transc_5panel_end.png';

const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka';
pres.title = '2菌種と抗菌薬の連成計算（' + (MIX ? '混合配置' : '分離配置') + '）';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', RULE = '8A99AD', ACC = 'A61E1E', LIGHT = 'F3F5F8';
const JP = 'Yu Mincho', TM = 'Times New Roman';
const LABEL = MIX ? '混合配置' : '分離配置';
let n = 0;

function frame(title) {
  const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addText(title, { x: 0.55, y: 0.32, w: 12.2, h: 0.6, fontFace: JP, fontSize: 24, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 0.98, w: 12.23, h: 0, line: { color: NAVY, width: 1.25 } });
  s.addShape(pres.shapes.LINE, { x: 0.55, y: 7.0, w: 12.23, h: 0, line: { color: RULE, width: 0.5 } });
  s.addText('西岡 佳祐　｜　抗菌薬を含む 2 菌種バイオフィルムの連成計算（' + LABEL + '）', { x: 0.55, y: 7.05, w: 10, h: 0.3, fontFace: JP, fontSize: 10, color: GREY, isTextBox: true, margin: 0 });
  s.addText(String(n), { x: 11.78, y: 7.05, w: 1.0, h: 0.3, fontFace: TM, fontSize: 11, color: GREY, align: 'right', isTextBox: true, margin: 0 });
  return s;
}
const cap = (s, t, y, x = 0.55, w = 12.23) => s.addText(t, { x, y, w, h: 0.5, fontFace: JP, fontSize: 12, color: INK, isTextBox: true, margin: 0, valign: 'top' });
const bullets = (s, items, o) => s.addText(items.map((t, i) => ({ text: t, options: { bullet: { indent: 18 }, breakLine: i < items.length - 1 } })),
  Object.assign({ fontFace: JP, fontSize: 16, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 10 }, o));
const img = (s, path, x, y, w) => { const [pw, ph] = pngSize(path); s.addImage({ path, x, y, w, h: w * ph / pw }); return w * ph / pw; };
function pngSize(p) { const b = fs.readFileSync(p); return [b.readUInt32BE(16), b.readUInt32BE(20)]; }

// 1 title
{ const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 2.05, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText('抗菌薬を含む 2 菌種バイオフィルムの連成計算', { x: 0.9, y: 2.25, w: 11.5, h: 0.8, fontFace: JP, fontSize: 32, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addText((MIX ? '混合配置' : '分離配置') + 'の例：初期・境界条件，結果，格子収束', { x: 0.9, y: 3.05, w: 11.5, h: 0.6, fontFace: JP, fontSize: 20, color: INK, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 3.85, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText([{ text: '西岡 佳祐', options: { fontSize: 18, bold: true, breakLine: true } },
    { text: 'Leibniz Universität Hannover, Institut für Kontinuumsmechanik（修士論文 4 章）', options: { fontSize: 14, breakLine: true } },
    { text: '2026 年 11 月', options: { fontSize: 14 } }],
    { x: 0.9, y: 4.3, w: 11.5, h: 1.5, fontFace: JP, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 6 }); }

// 2 purpose
{ const s = frame('1. 目的と問題設定'); s.addNotes('この計算の目的です。修論で作った連成モデル、つまり各積分点に Klempt 2026 の点モデルを置き、2024 の膜の場と栄養・抗菌薬でつないだモデルで、菌の組成が変わると成長と応力がどう変わるかを調べます。菌は Fritsch の合成データの 2 菌種で、抗菌薬あり・なしを同じ条件で比べます。数字はシードの平均で、成長と応力は 3 つの格子から外挿して誤差の幅も出しています。（約 1 分）');
  bullets(s, [
    '目的：点モデル（Klempt et al. 2026）を各積分点に置き，膜の場（Klempt et al. 2024，式 34–36）と栄養・抗菌薬の場を介して連成したとき，組成の変化が成長と応力に及ぼす影響を定量化する．',
    '対象：一辺 2 mm の立方体．菌は合成データ（Fritsch et al. 2025，case I，2 菌種，真のパラメータ）．',
    MIX ? '配置：菌 A と B を等量（φ₁ = φ₂ = 0.275）でシードに混ぜる（混合配置）．'
        : '配置：菌 A をシードに，菌 B を z 方向の 4 辺に沿った小領域に置く（分離配置）．',
    '比較：同一条件で抗菌薬の有無を比べる．評価量はシード領域の平均（生存体積分率，組成，栄養，膨張 α − 1，von Mises 応力）．',
    '数値の信頼性：成長と応力は 64³/128³/192³ の 3 格子から外挿し，GCI（Celik et al. 2008）で離散化誤差を示す．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 }); }

// 3 equations
{ const s = frame('2. 支配方程式'); s.addNotes('式は 3 本です。式 34 が膜の量で、拡散、成長、それから前線の項、マイナス r c/(k+c) n_c・∇φ が栄養の方へ膜を広げます。式 35 が栄養で、消費を生きている量 Σφ̄ で重み付けしたのがこの研究のつなぎ方です。抗菌薬も同じ形です。式 36 が膨張 α です。点モデルは場から c* と α* を受け取り、生きている量を返します。応力は α から小ひずみで出します。（約 1.5 分）');
  const h = img(s, EQ + 'eqs.png', 1.2, 1.3, 10.6);
  cap(s, 'φ₃D：膜の体積分率，α：膨張（F_g = α I，α(0) = 1），c：栄養，a：抗菌薬，φᵢ・ψᵢ：点モデルの体積分率と生存率，Σᵢφ̄ᵢ：生存体積分率．'
    + '式番号は Klempt et al. 2024．式 35 の消費は生存体積分率で重み付けし（本研究の連成），抗菌薬も同じ形とする．点モデル（Klempt et al. 2026）は公表された式のまま用いる．', 1.4 + h + 0.25);
}

// 4 constants
{ const s = frame('3. 材料定数と数値条件'); s.addNotes('定数はすべて論文の値です。Klempt 2024 の表 2 を 2 mm の立方体に換算しました。この研究で決めたのは、空隙の剛性の下限 f と、点モデルを始める閾値の 2 つで、どちらも収束を確かめています。数値は有限体積の 64³、前線の項は 2 次の風上差分、栄養は後退オイラーです。（約 1 分）');
  const H = (t) => ({ text: t, options: { bold: true, color: 'FFFFFF', fill: { color: NAVY }, align: 'center' } });
  const rows = [[H('記号'), H('値'), H('意味'), H('出典')],
    ['β', '0.02 mm²/T*', 'φ₃D の拡散', 'Klempt 2024 表 2（2 mm に換算）'],
    ['r, k', '10 mm/T*, 1', '前線の速さ，Monod 定数', '同上'],
    ['k_α', '10⁻³ T*⁻¹', '成長速度', '同上'],
    ['d, gℓ²/d', '10⁸ mm²/T*, 1', '栄養の拡散と消費', '同上（ℓ = 1 mm）'],
    ['E, ν', '10 Pa, 0.49', '膜の剛性', '同上'],
    ['f', '10⁻³', '空隙の剛性の下限', '本研究（f → 0 で 0.13 %）'],
    ['c*₀, s', '100, 0.1', '点モデルの栄養，時間の対応', 'Fritsch 2025 case I'],
    ['α*₀', '100（T* ≥ 1/3）', '抗菌薬の強さと開始', 'Klempt 2026 case 4'],
    ['φ_min', '10⁻³', '点モデルを始める閾値', '本研究（収束値）']]
    .map((r, i) => i === 0 ? r : r.map((t, j) => ({ text: t, options: { align: j === 2 || j === 3 ? 'left' : 'center', fontFace: j < 2 ? TM : JP } })));
  s.addTable(rows, { x: 0.55, y: 1.3, w: 7.6, colW: [1.2, 1.9, 2.2, 2.3], fontFace: JP, fontSize: 12, color: INK, rowH: 0.42,
    border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  s.addText('数値条件', { x: 8.6, y: 1.3, w: 4.2, h: 0.4, fontFace: JP, fontSize: 16, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  bullets(s, ['空間：有限体積，セル中心，64³（対称面で 1/4 領域）',
    '前線の項：2 次の風上（MUSCL＋minmod）',
    '時間：φ₃D, α は RK4，式 35 は後退オイラー（ċ 付き），点モデルは陰的オイラー＋Newton',
    '連成：ずらして解く．消費の Σᵢφ̄ᵢ は 1 ステップ遅れ',
    '応力：16³ の六面体，B̄ 法，小ひずみ',
    '計算時間（1 コア，400 ステップ）：64³ 約 5 分，128³ 約 1 時間，192³ 約 3.3 時間'],
    { x: 8.6, y: 1.8, w: 4.2, h: 5.0, fontSize: 13, paraSpaceAfter: 7 }); }

// 5 setup
{ const s = frame('4. 初期条件と境界条件'); s.addNotes('条件です。(a) 菌 A を中央のシードに、菌 B を z 方向の 4 辺の近くに置きます。三角は剛体の動きを止めるための支点で、それ以外は自由に膨らみます。(b) 栄養は手前の面の層で 1 に保ちます。(c) 抗菌薬は T* = 1/3 から外側の殻全体で 1 に保ちます。（約 1 分）');
  const h = img(s, R + 'assets/fig_setup_' + (MIX ? 'mixed' : 'separate') + '_ab.png', 1.95, 1.08, 9.4);
  const T = (t) => ({ text: t, options: { bold: true, color: 'FFFFFF', fill: { color: NAVY }, align: 'center' } });
  const c_ = (t, o = {}) => ({ text: t, options: Object.assign({ align: 'left', valign: 'middle' }, o) });
  const y0 = 1.08 + h + 0.05;
  s.addTable([[T('量'), T('どこで値を固定するか'), T('値'), T('いつから'), T('それ以外の壁')],
    [c_('栄養 c', { bold: true }), c_('y = −1 mm の面に接する要素の層（厚さ 0.25 mm，64 要素）'), c_('c = 1'), c_('T* = 0'), c_('流束 0（出入りなし）')],
    [c_('抗菌薬 a', { bold: true }), c_('6 面すべてに接する外側の要素の殻（厚さ 0.25 mm，296 要素）'), c_('a = 1'), c_('T* = 1/3（前は a = 0）'), c_('—（全面で固定）')],
    [c_('変位 u', { bold: true }), c_('角の 3 節点（図の▲）'), c_('0'), c_('—'), c_('自由（自由膨張）')]],
    { x: 0.55, y: y0, w: 12.23, colW: [1.25, 5.6, 0.9, 2.35, 2.13], fontFace: JP, fontSize: 11.5, color: INK, rowH: 0.33,
      border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  cap(s, '図 1・表　' + (MIX ? '混合配置（A と B を等量でシード（0.25 mm の要素 32 個，0.5 mm³）に混ぜる）' : '分離配置（菌 A：シード，0.25 mm の要素 32 個で 0.5 mm³．菌 B：0.25 × 0.25 × 1.0 mm のブロック 4 つ，z = −0.5〜0.5 mm，計 0.25 mm³）')
    + '．c と a は固定した層から拡散し，生存体積分率に取り込まれながら内側へ届く（式 35）．図は 8³ 要素で表示，計算は 64³ のセル．', y0 + 1.45, 0.55, 12.23); }

// 6 result figure
{ const s = frame('5. 結果：組成・場・応力の分布'); s.addNotes('結果です。上が抗菌薬を入れる直前、下が最後です。(a) 青が A、赤が B で、B は辺から広がって外側を覆います。(b) 生きている量は、抗菌薬の後に大きく下がり、栄養の来る下の方だけ少し残ります。(e) 応力はシードとその周りに集まります。上と下で応力が大きくなるのは、時間とともに膨張 α が増えるためで、抗菌薬のせいではありません。抗菌薬あり・なしの比較は次の表で見せます。(d1) が空なのは、抗菌薬がまだ入っていないからです。（約 1.5 分）');
  const h = img(s, F + (MIX ? 's12m_coupled_3d_mixed_ab.png' : 's12_coupled_3d_separate_ab.png'), 1.4, 1.15, 10.5);
  cap(s, '図 2　抗菌薬ありの' + LABEL + '．上段 T* = 0.325（抗菌薬の開始直前），下段 T* = 1．1/8 を除いて表示．'
    + '(a) 菌 A の割合 φ₁/(φ₁+φ₂)，(b) 生存体積分率 Σᵢφ̄ᵢ（(a)(b) は φ₃D > 0.5），(c) 栄養 c，(d) 抗菌薬 a，(e) von Mises 応力（16³）．Python，64³．', 1.15 + h + 0.1); }

// 7 video (if rendered)
if (fs.existsSync(VIDEO)) {
  const s = frame('6. 時間発展'); s.addNotes('同じ計算の動画です。T* = 1/3 で抗菌薬が外から入り、そこから緑の生きている量が減り続けます（シードで 0.97 から 0.15 へ）。一方で膜は広がり続け、膨張と応力は増え続けます。抗菌薬があってもなくても応力の増え方はほぼ同じ、というのが次の表です。（約 1 分）');
  const [pw, ph] = pngSize(POSTER); const w = 12.0, h = w * ph / pw;
  s.addMedia({ type: 'video', path: VIDEO, cover: 'data:image/png;base64,' + fs.readFileSync(POSTER).toString('base64'), x: 0.67, y: 1.3, w, h });
  cap(s, '図 3　図 2 と同じ計算の時間発展（T* = 0 → 1，40 コマ）．パネルは図 2 と同じ．クリックで再生．', 1.3 + h + 0.15);
}

// 8 table
{ const s = frame((fs.existsSync(VIDEO) ? '7' : '6') + '. 定量評価：シード領域の平均'); s.addNotes('数字です。生きている量は 0.965 から 0.152 へ、84 % 下がります。一方、膨張はほぼ同じで、応力は約 2 % 下がるだけです。分離の例は格子収束が良く、GCI は 0.1〜1.4 % なので、この 2 % は誤差より大きい差です。（約 1 分）');
  const H = (t) => ({ text: t, options: { bold: true, color: 'FFFFFF', fill: { color: NAVY }, align: 'center' } });
  const C = (t, o = {}) => ({ text: t, options: Object.assign({ align: 'center' }, o) });
  const L = (t) => C(t, { align: 'left' });
  const d = MIX
    ? [['0.961', '0.038', '−96 %'], ['0.859', '0.610', '−29 %'], ['0.608', '0.973', '+60 %'],
       ['8.74（GCI 4.2 %）', '8.71（4.3 %）', '−0.3 %'], ['4.86（GCI 31 %）', '4.68（41 %）', '−4 %']]
    : [['0.965', '0.152', '−84 %'], ['0.978', '0.951', '−3 %'], ['0.536', '0.911', '+70 %'],
       ['9.17（GCI 0.4 %）', '9.19（0.5 %）', '+0.2 %'], ['6.57（GCI 0.1 %）', '6.42（1.4 %）', '−2.3 %']];
  const names = ['生存体積分率 Σᵢφ̄ᵢ', '菌 A の割合 φ₁/(φ₁+φ₂)', '栄養 c', '膨張 α − 1 [10⁻⁴]（外挿）', 'von Mises 応力 [10⁻⁴ Pa]（外挿）'];
  s.addTable([[H('量'), H('抗菌薬なし'), H('抗菌薬あり'), H('相対変化')],
    ...names.map((nm, i) => [L(nm), C(d[i][0]), C(d[i][1]), C(d[i][2], i === 0 || i === 4 ? { bold: true, color: ACC } : {})])],
    { x: 1.2, y: 1.75, w: 10.9, colW: [4.0, 2.6, 2.4, 1.9], fontFace: JP, fontSize: 15, color: INK, rowH: 0.55,
      border: { type: 'solid', pt: 0.5, color: 'B8C2CF' } });
  s.addText('表 1　T* = 1 におけるシード領域の平均．生存体積分率・組成・栄養は 64³（外側の刻み 0.025），膨張と応力は 64³/128³/192³（刻み 0.0025）からの外挿値と GCI（付録 A）．',
    { x: 1.2, y: 1.25, w: 10.9, h: 0.45, fontFace: JP, fontSize: 12, color: INK, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.RECTANGLE, { x: 1.2, y: 5.4, w: 10.9, h: 1.25, fill: { color: LIGHT }, line: { color: 'D5DCE5', width: 0.5 } });
  s.addText(MIX
    ? '抗菌薬により生存体積分率はほぼ 0 まで低下する．膨張の変化は 0.3 % にとどまり，応力の差（約 4 %）は応力の GCI（31–41 %）より小さいため，有意な差とはいえない．'
    : '抗菌薬により生存体積分率は 84 % 低下するが，膜の体積分率 φ₃D は減らないため，膨張はほぼ不変であり，応力の低下は約 2 % にとどまる（GCI 0.1–1.4 % より大きい）．',
    { x: 1.45, y: 5.45, w: 10.4, h: 1.15, fontFace: JP, fontSize: 15, color: INK, isTextBox: true, margin: 0, valign: 'middle' }); }

// 9 discussion
{ const s = frame((fs.existsSync(VIDEO) ? '8' : '7') + '. 考察'); s.addNotes('なぜ応力がほとんど変わらないか。成長則の式 36 は、生きている菌と死んだ菌を区別せず、膜の量から膨張を決めるからです。組成は栄養と前線の速さを通してしか応力に届きません。組成を力学に効かせるには、生きている量で育つ成長則のような、もう一つの道が要ります。これが慶應での続きです。（約 1.5 分）');
  bullets(s, MIX ? [
    '組成は栄養と前線の速さ r c/(k+c) を通してのみ応力に届く．成長則（式 36）は生死を区別せず膜の量 φ₃D から α を決めるため，生存体積分率が大きく変わっても膨張と応力はほぼ変わらない．',
    '混合配置では両菌がシード内で抗菌薬にさらされ，分離配置（0.152）より生存体積分率が大きく下がる（0.038）．',
    '混合配置の応力は格子収束が遅い（p = 0.55，GCI 31 %）．前線の項による膜の縁の幅 β/v ≈ 0.004 mm に対し，192³ でもセルのペクレ数が 1.28 と 1 を超え，制限関数で精度が 1 次以下に落ちるためである．256³（Pe_h = 0.96）を加えると次数は 0.93 に上がり，GCI は 13 % に下がる（外挿 4.0×10⁻⁴ Pa）．',
    '組成が力学に直接効くには，生存体積分率で成長する法則などの第 2 の経路が要る（今後の課題）．']
  : [
    '組成は栄養と前線の速さ r c/(k+c) を通してのみ応力に届く．成長則（式 36）は生死を区別せず膜の量 φ₃D から α を決めるため，生存体積分率が大きく変わっても膨張と応力はほぼ変わらない．',
    '菌 B が先に減る：栄養の項 c*A₂₂/η₂ が A の半分である一方，抗菌薬の項 α*bᵢ/ηᵢ は両菌で等しい．栄養の乏しい所では A も死ぬ．',
    '分離配置の応力は格子収束が速い（p = 3.1，GCI 0.1 %）．抗菌薬による約 2 % の低下は離散化誤差より大きく，有意である．',
    '組成が力学に直接効くには，生存体積分率で成長する法則などの第 2 の経路が要る（今後の課題）．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 }); }

// 10 definitions
{ const s = frame('付録 A. 評価量の定義'); s.addNotes('質問用です。von Mises は偏差応力の大きさ、GCI は Celik の 3 格子の方法で格子誤差の幅、ペクレ数は前線の移流と拡散の比を 1 セルで見たものです。');
  const col = (x, w, t, eq, eqw, body) => {
    s.addShape(pres.shapes.RECTANGLE, { x, y: 1.25, w, h: 5.55, fill: { color: LIGHT }, line: { color: 'D5DCE5', width: 0.5 } });
    s.addText(t, { x: x + 0.2, y: 1.38, w: w - 0.4, h: 0.45, fontFace: JP, fontSize: 17, bold: true, color: NAVY, isTextBox: true, margin: 0 });
    const h = img(s, EQ + eq, x + (w - eqw) / 2, 1.95, eqw);
    s.addText(body.map((t, i) => ({ text: t, options: { breakLine: i < body.length - 1 } })),
      { x: x + 0.2, y: 2.1 + h, w: w - 0.4, h: 4.6 - h, fontFace: JP, fontSize: 12.5, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 6 }); };
  col(0.55, 3.75, 'von Mises 応力', 'vm.png', 3.35, [
    's は偏差応力．静水圧を除いた応力の大きさで，一様な膨張だけなら 0 となる．',
    '表 1 の値はシード領域（初期に菌がいた 32 要素）での平均．']);
  col(4.45, 4.5, 'GCI（格子収束指数）', 'gci.png', 4.1, [
    'Celik et al. 2008 の 3 格子法．格子 1, 2, 3 = 192³, 128³, 64³，安全率 1.25．',
    'p は観測された収束次数（反復で解く），f_ext は外挿値．',
    'GCI は最も細かい格子の値が，格子を無限に細かくした値からどれだけ離れうるかを表す．p が小さいと分母が小さく GCI は大きい（分離の応力 p = 3.1 → 0.1 %，混合の応力 p = 0.55 → 31 %）．']);
  col(9.05, 3.73, 'セルのペクレ数', 'pe.png', 3.4, [
    'v は前線の項 −r c/(k+c) n_c·∇φ₃D の速さ（r = 10 mm/T*，k = 1，c ≈ 0.96），β = 0.02 mm²/T*，N は一辺のセル数．',
    '移流と拡散の比を 1 セルの幅で見た量．膜の縁の幅は約 β/v ≈ 0.004 mm．',
    'Pe_h > 1 では縁がセルより細く，制限関数で精度が落ちる．128³：1.9，192³：1.28，256³：0.96．']); }

// 11 references
{ const s = frame('参考文献'); s.addNotes('参考文献です。');
  const refs = [
    'F. Klempt, M. Soleimani, P. Wriggers, P. Junker: A Hamilton principle-based model for diffusion-driven biofilm growth. Biomech. Model. Mechanobiol. 23, 2091–2113 (2024). doi:10.1007/s10237-024-01883-x',
    'F. Klempt, H. Geisler, M. Soleimani, P. Junker: A continuum multi-species bacterial growth model with a novel interaction scheme. Arch. Appl. Mech. 96, 164 (2026). doi:10.1007/s00419-026-03160-y',
    'L. Fritsch, H. Geisler, J. Grashorn, F. Klempt, M. Soleimani, M. Broggi, P. Junker, M. Beer: Bayesian updating of constitutive parameters under hybrid uncertainties with a novel surrogate model applied to biofilms. arXiv:2512.15145 (2025)',
    'I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, P. E. Raad: Procedure for estimation and reporting of uncertainty due to discretization in CFD applications. J. Fluids Eng. 130, 078001 (2008). doi:10.1115/1.2960953'];
  s.addText(refs.map((t, i) => ({ text: '[' + (i + 1) + ']  ' + t, options: { breakLine: i < refs.length - 1 } })),
    { x: 0.75, y: 1.4, w: 11.8, h: 5.2, fontFace: TM, fontSize: 14, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 14 }); }

pres.writeFile({ fileName: F + 'report_' + CASE + '_ab_academic.pptx' }).then(f => console.log('ok', f, n, 'slides'));
