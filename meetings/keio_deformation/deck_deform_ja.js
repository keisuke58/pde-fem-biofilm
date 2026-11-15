// 変形との連成（慶應の修論の柱）：Klempt 2024 式 30 の力学の項．学術的な日本語のスライド，ノートつき．結果が出るたびに書き足す．
//   cd meetings/keio_deformation && bash build.sh   -> build/deform_ja.pptx (not in git)
// 式：eq/*.tex -> build/eq/*.png（build.sh）．図：assets/fig_mech_feedback_section*.png など．
const pptxgen = require('pptxgenjs');
const fs = require('fs');
const path = require('path');
const R = path.resolve(__dirname, '../..') + '/';
const B = __dirname + '/build/', EQ = B + 'eq/', AS = R + 'assets/';

const pres = new pptxgen(); pres.layout = 'LAYOUT_WIDE'; pres.author = 'Keisuke Nishioka';
pres.title = '変形との連成：Klempt (2024) 式 30 の力学の項';
const NAVY = '1E3A5F', INK = '1F2933', GREY = '5F6B7A', RULE = '8A99AD', ACC = 'A61E1E', LIGHT = 'F3F5F8';
const JP = 'Yu Mincho', TM = 'Times New Roman';
const FOOT = '西岡 佳祐　｜　変形との連成（慶應での続き）';
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

// 1 title
{ const s = pres.addSlide(); n += 1; s.background = { color: 'FFFFFF' };
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 2.05, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText('バイオフィルムの成長と変形の双方向連成', { x: 0.9, y: 2.25, w: 11.5, h: 0.8, fontFace: JP, fontSize: 32, bold: true, color: NAVY, isTextBox: true, margin: 0 });
  s.addText('Klempt et al. (2024) の式 30 で省かれた力学の項を戻す', { x: 0.9, y: 3.05, w: 11.5, h: 0.6, fontFace: JP, fontSize: 20, color: INK, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.LINE, { x: 0.9, y: 3.85, w: 11.5, h: 0, line: { color: NAVY, width: 1.5 } });
  s.addText([{ text: '西岡 佳祐', options: { fontSize: 18, bold: true, breakLine: true } },
    { text: '慶應義塾大学 大学院 理工学研究科 開放環境科学専攻（村松研究室）', options: { fontSize: 14, breakLine: true } },
    { text: '作業中の記録（結果が出るたびに更新）', options: { fontSize: 14 } }],
    { x: 0.9, y: 4.3, w: 11.5, h: 1.5, fontFace: JP, color: INK, isTextBox: true, margin: 0, valign: 'top', paraSpaceAfter: 6 });
  s.addNotes('慶應での修士研究の柱の一つ、成長と変形の双方向の連成についての記録です。Klempt ら 2024 の式 30 で省かれている力学の項を戻すと何が起きるかを調べています。結果が出るたびにこのスライドに書き足します。');
}

// 2 one-way vs two-way
{ const s = frame('1. 一方向の連成と双方向の連成');
  fit(s, EQ + 'flow.png', 0.9, 1.3, 11.5, 1.3);
  bullets(s, [
    '今のモデル（LUH の修論・連成の論文）：成長 α̇ = k_α φ₃D は量 φ₃D だけで決まり，応力 σ は計算の最後に出るだけ．応力は φ₃D にも α にも戻らない．',
    '双方向：変形（弾性ひずみ）が φ₃D の式に戻り，φ₃D → α → σ → φ₃D の輪になる．',
    'たとえ：一方向は「風船をふくらませるとゴムが張る」だけ．双方向は「ゴムが張るほど空気が入りにくくなる」．'],
    { x: 0.75, y: 3.0, w: 11.8, h: 3.8 });
  s.addNotes('今のモデルでは、栄養からバイオフィルムの量、量から成長、成長から応力へと一方向にしか情報が流れません。応力は最後に出てくるだけで、成長には戻りません。双方向の連成では、変形がバイオフィルムの量の式に戻り、輪になります。');
}

// 3 the term
{ const s = frame('2. 式 30 の力学の項');
  fit(s, EQ + 'eq30.png', 0.7, 1.2, 11.9, 2.6);
  bullets(s, [
    '自由エネルギーの弾性の部分 (μ/2) φ² (I:C_e − 3) を φ で微分した項．硬さが φ² に比例するので，変形したバイオフィルムは φ を変えて弾性エネルギーを下げようとする．',
    'Klempt らは「μ が他の定数より何桁も小さい」として省いた（式 33）．μ = 3.36 Pa（E = 10 Pa）．実測の口腔バイオフィルムは 0.6〜41 kPa（Pattem 2018）．',
    'K = μ/η_φ の η_φ は論文にない → K を掃引する．小さいひずみでは I:C_e − 3 ≈ 2 tr ε_e：引っ張られた所で φ₃D が減り，押された所で増える．'],
    { x: 0.75, y: 4.0, w: 11.8, h: 2.9, fontSize: 14 });
  s.addNotes('式 30 の最後の項が対象です。自由エネルギーの弾性の部分を φ で微分すると出てくる項で、ハミルトン原理の枠組みの中の項です。Klempt らは硬さが小さいとして省きましたが、実測の硬さは論文の値の 100 倍から 4000 倍あります。係数 K の値は論文にないので、掃引して調べます。符号は、引っ張られた所で φ が減り、押された所で増える向きです。');
}

// 4 implementation
{ const s = frame('3. 実装（Python の連成コード）');
  bullets(s, [
    '`coupled_living_consumption.py --mech-feedback K`：外側の時間ステップの始めに 16³ の応力格子で応力を解き，各セルに 2 tr ε_e を渡す（1 ステップの遅れ）．',
    'φ₃D の式に −K φ₃D (I:C_e − 3) を足す．他の式（栄養，抗菌薬，成長，点モデル）は変えない．K = 0 で既存の結果と一致．',
    '計算の手間：応力を毎ステップ解くので，32³ で約 10 倍（11 s → 110 s）．',
    '比べるためのモデル：恒常圧の成長則 α̇ = k_α φ₃D · max(0, 1 − p/p_h)（現象論，`--homeostatic`）．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.4 });
  s.addNotes('実装は Python の連成コードに入れました。各時間ステップの始めに応力を解き、弾性ひずみの体積成分を各セルに渡して、バイオフィルムの量の式に項を足します。ほかの式は変えていません。比べるために、圧力で成長を止める現象論の式も入っています。');
}

// 5 result: free growth
{ const s = frame('4. 結果：自由に育つ例（32³，K = 10⁵）');
  fit(s, AS + 'fig_mech_feedback_section_K100000.png', 0.55, 1.15, 12.23, 3.6);
  bullets(s, [
    '引っ張られた外側で φ₃D が最大 0.5 減り，押された種子の下側ではわずかに増える（+0.07）．',
    '外側のバイオフィルムが薄くなり，種子を押さえる周りが柔らかくなる → 外側の von Mises は 12 → 4〜5（×10⁻⁴ Pa）．',
    '種子の平均：α − 1 −4.6 %，von Mises −26 %（32³）．192³ では −3.7 %，−20 %（次のスライド）．'],
    { x: 0.75, y: 4.95, w: 11.8, h: 1.9, fontSize: 14 });
  cap(s, '分けた例・抗菌薬あり，T* = 1，z = 0 断面．(a) 項なしの平均応力（赤：引張，青：圧縮），(b) 項による φ₃D の変化，(c) α − 1 の変化，(d) y = 0 での von Mises．', 4.7, 0.55, 12.23, 10);
  s.addNotes('自由に育つ例です。引っ張られている外側でバイオフィルムの量が減り、押されている種子の下側ではわずかに増えます。外側が薄くなると周りが柔らかくなるので、応力が下がります。種子の平均で、成長は約 5 パーセント、応力は約 4 分の 1 小さくなりました。まだ粗い格子の値です。');
}

// 5b K sweep
{ const s = frame('4b. 結果：K の掃引（64³〜192³，自由に育つ例）');
  fit(s, AS + 'fig_mech_K_sweep.png', 0.55, 1.15, 7.6, 3.2);
  cap(s, '種子の平均の変化（K = 0 の同じ格子に対して，T* = 1，抗菌薬あり）．灰色の帯は ±1 %．点線より右は発散．', 4.4, 0.55, 7.6, 10);
  s.addTable([
    [{ text: 'K [1/T*]', options: { bold: true } }, { text: '分けた例', options: { bold: true } }, { text: '混ぜた例', options: { bold: true } }],
    ['10³', '−0.01 %，−0.04 %', '−0.1 %，+0.1 %'],
    ['10⁴', '−0.2 %，−0.5 %', '−1.1 %，+1.3 %'],
    ['3×10⁴', '−0.7 %，−2.4 %', '−3.4 %，+0.4 %'],
    ['10⁵', '−3.7 %，−20 %', '−9.0 %，−19 %'],
    ['3×10⁵', '発散', '−19 %，−68 %'],
    ['10⁶', '発散', '発散']],
    { x: 8.4, y: 1.25, w: 4.4, colW: [1.1, 1.65, 1.65], fontFace: JP, fontSize: 11, color: INK, border: { type: 'solid', color: 'BBBBBB', pt: 0.5 }, align: 'center', valign: 'middle', rowH: 0.36 });
  cap(s, '各欄：α − 1 の変化，von Mises の変化（最も細かい格子）', 3.85, 8.4, 4.4, 10);
  bullets(s, [
    '境目：K ≈ 10⁴（混ぜた例）〜3×10⁴（分けた例）から 1 % 以上変わる．K ≤ 10³ では項は効かない（Klempt らが省いた範囲）．',
    'K = 10⁵ の α − 1 は格子で収束（GCI 0.1〜0.8 %）．応力は混ぜた例で GCI 4.6 %，分けた例は収束が遅い（GCI 24 %）．',
    '注意：K = 10⁵ の応力は外側の時間刻みに敏感（128³ で刻みを半分にすると +7.9 %）．K ≥ 3×10⁵ は発散．原因は次のスライド．'],
    { x: 0.75, y: 4.85, w: 11.8, h: 2.0, fontSize: 14 });
  s.addNotes('係数 K を千から百万まで変えて、三つの格子で確かめました。K が一万から三万のあたりから、成長や応力が 1 パーセント以上変わり始めます。K が千以下では項は効きません。K が十万では成長が 4 から 9 パーセント、応力が約 2 割小さくなり、成長の変化は格子で収束しています。ただし応力は時間の刻みに敏感で、刻みを半分にすると 8 パーセント変わります。K が三十万以上では計算が発散します。');
}

// 5c checks
{ const s = frame('4c. 確認：時間刻みと発散の原因（32³，分けた例）');
  s.addTable([
    [{ text: 'K', options: { bold: true } }, { text: '計算の仕方', options: { bold: true } }, { text: '外側の刻み', options: { bold: true } }, { text: 'α − 1', options: { bold: true } }, { text: 'von Mises', options: { bold: true } }, { text: 'φ₃D の最大', options: { bold: true } }],
    ['10⁵', 'そのまま', '0.025', '−4.6 %', '−25.6 %', '1.32'],
    ['10⁵', 'そのまま', '0.0125', '−4.6 %', '−23.1 %', '1.33'],
    ['10⁵', 'そのまま', '0.00625', '−4.5 %', '−20.9 %', '1.33'],
    ['10⁵', '応力を 1 ステップに 4 回解き直す', '0.025', '−4.8 %', '−27.6 %', '1.33'],
    ['10⁵', 'φ₃D を [0, 1] に切る（モデルの変更）', '0.025', '−6.0 %', '−34.5 %', '1.00'],
    ['3×10⁵', 'そのまま', '0.025', '発散', '発散', '—'],
    ['3×10⁵', '項を厳密に積分', '0.025', '−13 %', '−59 %', '6.0'],
    ['3×10⁵', '項を厳密に積分', '0.00625', '−19 %', '−93 %', '11']],
    { x: 0.75, y: 1.2, w: 11.8, colW: [1.1, 4.6, 1.5, 1.4, 1.6, 1.6], fontFace: JP, fontSize: 12, color: INK, border: { type: 'solid', color: 'BBBBBB', pt: 0.5 }, align: 'center', valign: 'middle', rowH: 0.34 });
  cap(s, '変化は同じ刻みの K = 0 に対して（T* = 1，種子の平均，抗菌薬あり）．', 4.35, 0.75, 11.8, 10);
  bullets(s, [
    '成長の変化（α − 1：−4.6 %）は刻みによらない．応力の変化は刻みを半分にするたびに約 2 ポイント動き，まだ収束していない．',
    '原因は項の 1 ステップの遅れではない：応力を解き直しても，項を厳密に積分しても直らない．',
    'K ≥ 3×10⁵ の発散は数値ではなく項そのもの：計算は回せるが，圧縮された所で φ₃D が 6〜11 まで増える．弾性エネルギーが圧縮の下で φ に下限を持たないため（柱の例と同じ）．',
    'まとめ：K ≲ 3×10⁴ は小さな補正として使える．K ≈ 10⁵ は成長は使えるが応力は目安．それ以上は φ₃D の上限かエネルギーの形の見直しが必要．'],
    { x: 0.75, y: 4.75, w: 11.8, h: 2.2, fontSize: 13 });
  s.addNotes('時間の刻みと発散の原因を粗い格子で確かめました。成長の変化は刻みによらず約 5 パーセントですが、応力の変化は刻みを半分にするたびに 2 ポイントほど動き、まだ収束していません。項が 1 ステップ遅れることが原因かと考えましたが、応力を解き直しても直りませんでした。K が三十万以上で発散するのは数値の問題ではなく、圧縮された所で項がバイオフィルムの量を増やし続けるためです。体積分率が 1 を大きく超えるので、そのままでは使えません。K が三万くらいまでは小さな補正として使えます。');
}

// 6 result: confinement (bars)
{ const s = frame('5. 結果：剛な柱で閉じ込めた例（K = 10⁵）');
  bullets(s, [
    'そのままでは発散：圧縮された所で項が φ₃D を増やし，罰則（5 /T*）より強く φ₃D > 1 が増え続ける．弾性エネルギー (μ/2) φ² (I:C_e − 3) は圧縮の下で φ に下限がない．',
    'φ₃D を毎ステップ [0, 1] に切る上限（`--phi-clip`，モデルの変更）を入れると最後まで回る．K = 0 では上限の有無で結果は同じ．'],
    { x: 0.75, y: 1.2, w: 11.8, h: 1.4, fontSize: 14 });
  fit(s, AS + 'fig_mech_feedback_section_bars_K100000_clip.png', 0.55, 2.65, 12.23, 3.0);
  bullets(s, [
    '柱のあいだは全体が圧縮 → φ₃D が上限まで増え，α − 1 が +33 %，種子の von Mises は 4.2 倍（30 → 127 ×10⁻⁴ Pa）．自由な例（−26 %）と逆向き．',
    '閉じ込めの下では正の帰還（圧縮 → 量が増える → 成長と圧縮が増える）．圧縮で成長が遅れる実験（Delarue 2016）とも，恒常圧の考え方とも逆．項の符号やエネルギーの形を Klempt 氏と確認する．'],
    { x: 0.75, y: 5.7, w: 11.8, h: 1.25, fontSize: 13 });
  s.addNotes('剛な柱で閉じ込めた例です。そのままでは発散しました。バイオフィルムの量に上限を入れると最後まで回りますが、柱のあいだが全体に圧縮されるので量が増え、成長が 3 割増え、応力は 4 倍になりました。自由に育つ例とは逆向きで、閉じ込めの下では圧縮が量を増やし、さらに圧縮を強める正の帰還になります。圧縮で成長が遅れるという実験とも逆なので、この項の符号やエネルギーの形を Klempt さんと確認します。');
}

// 7 position in the literature
{ const s = frame('6. 先行研究との関係');
  bullets(s, [
    '応力で成長が変わる式は，生体組織や腫瘍の成長力学では標準（恒常応力）．',
    '微生物：閉じ込めた酵母は約 1 MPa の接触圧で成長が遅れる（Delarue et al. 2016）．大腸菌は 34 kPa の外からの圧縮では成長速度がほぼ変わらない（Si et al. 2015）．',
    'バイオフィルムの力学モデルはほとんど一方向（成長 → 応力 → 皺・座屈の形）．圧力で成長が変わる連続体モデルは見当たらない．',
    '多菌種の連続体モデル（Klempt et al. 2025, arXiv:2509.01274）はあるが，組成・成長・応力の帰還を一つにした例は見当たらない．',
    '課題：モデルの応力（mPa，実測の硬さで 0.01〜0.1 Pa）は止まる圧力（MPa）よりずっと小さい → 閉じ込めの下で効くかが問い．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.5, fontSize: 15 });
  s.addNotes('先行研究です。応力で成長が変わる考え方は組織の分野では標準で、微生物でも圧力で成長が遅れる実験があります。一方、バイオフィルムの力学モデルはほとんど一方向です。組成と成長と応力の帰還を一つにしたモデルは見当たりませんでした。ただ、モデルの応力は実際に成長が止まる圧力よりずっと小さいので、どの条件で効くかが問いになります。');
}

// 8 next
{ const s = frame('7. 次にやること');
  bullets(s, [
    'K の掃引と確認は済み（境目 K ≈ 10⁴〜3×10⁴；大きい K の発散は項そのもの）．残り：K = 10⁵ の応力がなぜ時間刻みに敏感かを調べる．',
    '閉じ込めの下の正の帰還をどう扱うか：上限，体積変化の項，項の符号（エネルギーの形）を Klempt 氏と確認．歯頸部の窓でも同じ比較．',
    '実測の硬さ（kPa）で K と応力の大きさを読み替える．η_φ の目安を Klempt 氏に聞く．',
    '組成で変わる硬さ（菌種・生死で E を変える）と組み合わせる．',
    '多菌種への拡張は IKM の多菌種モデルとの切り分けを確認してから．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.5 });
  s.addNotes('次にやることです。係数 K の範囲を広げて、どこから効くかと格子での収束を確かめます。閉じ込めのある例で同じ比較をし、実測の硬さで大きさを読み替えます。');
}

// appendix
{ const s = frame('付録：記号');
  bullets(s, [
    'φ₃D：バイオフィルムの体積分率，α：成長（F_g = α I，α(0) = 1），c：栄養，σ：応力，ε_e：弾性ひずみ，C_e：弾性の右コーシー・グリーン テンソル．',
    'μ：せん断弾性係数（Klempt 2024：3.36 Pa），η_φ：φ の粘性（論文に値なし），K = μ/η_φ．',
    'p：圧力（圧縮が正），p_h：恒常圧（比べるためのモデル）．',
    '格子：一辺 2 mm の立方体を n³ のセル，応力は 16³ の六面体要素（φ₃D と α を要素ごとに平均）．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.5, fontSize: 15 });
  s.addNotes('付録です。記号の一覧です。');
}

pres.writeFile({ fileName: B + 'deform_ja.pptx' }).then(f => console.log('wrote', f));
