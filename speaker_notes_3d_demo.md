# Speaker notes: Biofilm Growth and Stress on a Real Molar (slides_3d_demo, 20 Oct 2026)

About 30 seconds per slide. EN: a text that can be said aloud. JA: the same meaning.
The PowerPoint version (`slides_3d_demo.pptx`) carries the EN text in the notes field of each slide.

## 1. Title
EN: I will show the growth model of my thesis in three dimensions, first on the thesis test problem and then on a real
molar, with stresses that I checked against mesh refinement.
JA: 修論の成長モデルを 3 次元で見せます。最初に修論の検証問題、次に本物の大臼歯の上で、メッシュを細かくして確かめた応力を示します。

## 2. What I show
EN: The growth model is the one of Klempt et al. 2024. I added a nutrient that the saliva supplies and the biofilm uses up,
and a limit on growth by pressure. Only the first case is a thesis result, every stress I quote is mesh-converged, and no
experimental data are used.
JA: 成長モデルは Klempt ら 2024 のものです。唾液から来て膜が消費する栄養と、圧力で成長を止める則を加えました。修論の結果は最初の例
だけで、示す応力はすべてメッシュで収束を確かめたものです。実験データは使っていません。

## 3. Verification
EN: This is the verification problem of Chapter 4: a seed in the centre of the 2 mm cube. The seed stress at the end agrees
with the thesis value to 0.2 percent, and the stress forms a ring at the edge of the seed.
JA: 第 4 章の検証問題で、2 mm の立方体の中心にシードがあります。最後のシードの応力は修論の値と 0.2 % で一致し、応力はシードの縁に
輪のように出ます。

## 4. Stress next to a fixed surface
EN: When the seed touches a clamped face, the stress concentrates in the biofilm next to it. The mean over that layer
converges, to within one percent, and is ten to twelve times the free-seed value. The pointwise maximum at the corner keeps
growing with refinement, so I only quote means and percentiles.
JA: シードが固定面に接すると、接する層で応力が集中します。その層の平均は 1 % 以内に収束し、自由なシードの 10〜12 倍です。角の
点の最大値は細かくするほど増えるので、平均と百分位点だけを言います。

## 5. A real molar
EN: This is a lower first molar at its real size, from a public data set under CC BY 4.0. A half-millimetre film sits on one
flank near the gum line. The growth field converges to two percent; the stress colours come from the film model of the next
slide.
JA: CC BY 4.0 の公開データの下顎第一大臼歯を実寸で使っています。歯茎の近くの片側に 0.5 mm の膜があります。成長の場は 2 % まで
収束していて、応力の色は次のスライドの膜のモデルから来ています。

## 6. Resolving the film
EN: On meshes of the whole box the film is only one or two elements thick and the stress does not converge. I solve the
film alone with the tooth as rigid, down to 28 micrometres. With six elements across the film the mean converges to
3.3 millipascal, and the highest values sit in the tenth of a millimetre next to the tooth.
JA: 箱全体のメッシュでは膜が要素 1〜2 個分しかなく、応力が収束しません。歯を剛体として膜だけを 0.028 mm まで細かく解くと、膜を
6 要素以上で切れば平均が 3.3 mPa に収束します。一番大きいのは歯に接する 0.1 mm です。

## 7. Nutrient and the growth limit
EN: With the nutrient used up in the film, growth and stress drop to about one seventh; at the paper's ratio to two thirds.
The pressure limit alone gives about 0.7. The consumption rate and the threshold are my choices.
JA: 膜の中で栄養が使われると、成長と応力は約 1/7 に下がります。論文の比では約 2/3 です。圧力の則だけでは約 0.7 です。消費の速さと
しきい値は私が選んだ値です。

## 8. Tooth-like collar
EN: On a crown-shaped collar I use the four-species synthetic case of Fritsch et al. 2025 with their true parameters, so no
experimental data enter. Growth is large here, so this run is finite strain.
JA: 歯冠形のカラーでは、Fritsch ら 2025 の 4 菌種の合成ケースを真の値のまま使っているので、実験データは入っていません。成長が
大きいので、この計算は有限ひずみです。

## 9. Section and composition
EN: Species three takes over and species four almost disappears. Without the limit the stress reaches 0.68 times E at the
bonded base; with it, 0.12 times E.
JA: 種 3 が優勢になり、種 4 はほぼ消えます。則なしでは接着した付け根で 0.68E、則ありでは 0.12E です。

## 10. Summary
EN: The 3D view reproduces the thesis value. Next to fixed surfaces the stress is ten to twelve times the free value. On the
real molar the film needs six elements across it, and then the film mean is 3.3 millipascal. Next at Keio: the calibrated
parameters on this tooth, and an implant.
JA: 3 次元でも修論の値を再現します。固定面のそばでは応力が自由な場合の 10〜12 倍です。本物の大臼歯では膜を 6 要素以上で切る必要
があり、そのとき膜の平均は 3.3 mPa です。慶應での次は、較正した値をこの歯に入れることと、インプラントです。
