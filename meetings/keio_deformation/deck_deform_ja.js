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
    '種子の平均：α − 1 −4.6 %，von Mises −26 %．32³ の目安（64³〜192³ の確認は fifa で実行中）．'],
    { x: 0.75, y: 4.95, w: 11.8, h: 1.9, fontSize: 14 });
  cap(s, '分けた例・抗菌薬あり，T* = 1，z = 0 断面．(a) 項なしの平均応力（赤：引張，青：圧縮），(b) 項による φ₃D の変化，(c) α − 1 の変化，(d) y = 0 での von Mises．', 4.7, 0.55, 12.23, 10);
  s.addNotes('自由に育つ例です。引っ張られている外側でバイオフィルムの量が減り、押されている種子の下側ではわずかに増えます。外側が薄くなると周りが柔らかくなるので、応力が下がります。種子の平均で、成長は約 5 パーセント、応力は約 4 分の 1 小さくなりました。まだ粗い格子の値です。');
}

// 6 result: confinement (bars)
{ const s = frame('5. 結果：剛な柱で閉じ込めた例（K = 10⁵）：発散');
  bullets(s, [
    '計算が発散した（16³・32³）．閉じ込めで圧縮された種子では I:C_e − 3 < 0（16³ で約 −1×10⁻⁴）で，項 −K φ₃D (I:C_e − 3) が φ₃D を増やす．',
    'K = 10⁵ では増える速さが約 10 /T* になり，φ₃D > 1 を押し戻す罰則（係数 5）より強い → φ₃D が 1 を超えて増え続ける（2.1 まで確認）．',
    '理由：弾性エネルギー (μ/2) φ² (I:C_e − 3) は，圧縮の下では φ を増やすほど下がり，下限がない．自由に育つ例では圧縮が弱いので表に出なかった．',
    '数値の対処は済み：項を各ステップで厳密に積分（`--mech-exp`，自由な例で一致を確認），空所（φ₃D < 0.05）では項を切る（`--mech-phi-min 0.05`）．それでも発散は残る．',
    '必要なもの：φ₃D ≤ 1 を守る仕組み（強い罰則または上限），または体積変化の項（ln J など）を含むエネルギー．Klempt 氏に相談する．'],
    { x: 0.75, y: 1.35, w: 11.8, h: 5.5, fontSize: 15 });
  s.addNotes('剛な柱で閉じ込めた例では計算が発散しました。圧縮された種子でこの項がバイオフィルムの量を増やし、上限を押し戻す罰則より強くなって、量が 1 を超えて増え続けます。もとの弾性エネルギーは、圧縮の下では量を増やすほど下がり、下限がないためです。項を戻すには、量の上限を守る仕組みか、体積変化の項を含むエネルギーが必要で、Klempt さんに相談します。');
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
    'K の掃引（10³〜10⁶）を 64³〜192³ で：どの K から α − 1 や応力が 1 % 以上変わるかの境目と格子収束（fifa で実行中）．',
    '閉じ込めで発散しない形にする（φ₃D ≤ 1 の上限，体積変化の項）．そのうえで剛な柱，歯頸部の窓（歯面）で同じ比較．',
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
