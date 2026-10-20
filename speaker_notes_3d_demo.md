# Speaker notes: Biofilm Growth and Stress in 3D (slides_3d_demo, 20 Oct 2026)

About 30 seconds per slide. EN: a text that can be said aloud. JA: the same meaning.
The PowerPoint version (`slides_3d_demo.pptx`) carries the EN text in the notes field of each slide.

## 1. Title
EN: I will show the growth model of my thesis in three dimensions, first the verification problem and then three
demonstrations with the same model.
JA: 修論の成長モデルを 3 次元で見せます。最初に検証問題、そのあと同じモデルを使った 3 つの例です。

## 2. What the videos show
EN: The videos show the biofilm phase, the growth and the von Mises stress over time. The growth model is the one of
Klempt et al. 2024. I added a homeostatic-pressure law that slows growth when the pressure approaches p_h. Only the
first case is a thesis result; the others are demonstrations, and no experimental data are used.
JA: 動画は、バイオフィルムの相、成長、von Mises 応力の時間変化です。成長モデルは Klempt ら 2024 のもので、圧力が p_h に近づくと
成長を遅くする恒常圧の則を私が加えました。修論の結果は最初の例だけで、ほかは例示です。実験データは使っていません。

## 3. Verification problem
EN: This is the verification problem of Chapter 4: a seed in the centre of the 2 mm cube. One octant is cut away. The
seed stress at the end is 6.19 times 10 to the minus 4 pascal, which agrees with the thesis value to 0.2 percent. The
stress forms a ring at the edge of the seed.
JA: 第 4 章の検証問題です。2 mm の立方体の中心にシードがあり、1/8 を切り取って見せています。最後のシードの応力は 6.19e−4 Pa
で、修論の値と 0.2 % で一致します。応力はシードの縁に輪のように出ます。

## 4. Seed in a corner and on an edge
EN: Here I move the seed to a corner and to an edge and clamp the bottom face, like a substrate. The biofilm spreads
along the walls, and the stress concentrates where it meets the clamped face, about twenty times the centre value.
The peak sits at the clamped edge and depends on the mesh, so the colour scale stops at the 95th percentile.
JA: シードを角と辺に動かし、底面を基材のように固定しました。バイオフィルムは壁に沿って広がり、固定面と接するところで応力が
集中します。中心の場合の約 20 倍です。最大値は固定面の縁にありメッシュに依存するので、色の上限は 95 % 点にしました。

## 5. Two seeds
EN: Two seeds on the clamped face, a little apart. They grow into each other, and the stress is highest under each
seed, at the clamped face.
JA: 固定面の上に、少し離して 2 つのシードを置きました。広がって互いにつながり、応力はそれぞれのシードの下、固定面のところで
最大になります。

## 6. The growth law
EN: The same corner seed, with and without the homeostatic-pressure law, on one colour scale. Without the law the stress
at the clamped edge reaches about 15 millipascal. With the law the growth stops where the pressure reaches p_h, and the
concentration goes away. The value of p_h is chosen by me.
JA: 同じ角のシードで、恒常圧の則なしとありを同じ色の目盛りで比べています。則なしでは固定面の縁で約 15 mPa になります。
則ありでは圧力が p_h に達したところで成長が止まり、集中が消えます。p_h の値は私が決めたものです。

## 7. Tooth-like collar
EN: Now a biofilm layer on a crown-shaped collar, as around a tooth or an implant. The point model uses the
four-species synthetic case of Fritsch et al. 2025 with their true parameters, so no experimental data enter. One
quarter is cut away so that the layer is visible.
JA: 次は、歯やインプラントのまわりのような、歯冠形のカラーの上のバイオフィルム層です。点モデルには Fritsch ら 2025 の
4 菌種の合成ケースを真の値のまま使っているので、実験データは入っていません。層が見えるように 1/4 を切り取っています。

## 8. Section and composition
EN: On the left, species three takes over the seed and species four almost disappears. On the right, the section at 45
degrees: without the law the stress reaches 0.68 times E at the bonded base; with the law the largest value is 0.12
times E.
JA: 左は、シードで種 3 が優勢になり、種 4 はほぼ消えることを示します。右は 45° の断面で、則なしでは接着した付け根で
0.68E、則ありでは最大でも 0.12E です。

## 9. Summary
EN: The 3D view reproduces the thesis value. When the seed touches a wall, the stress concentrates where growing biofilm
meets a fixed surface, and the homeostatic-pressure law removes this concentration. The next step at Keio is a real
tooth surface with the calibrated parameters.
JA: 3 次元でも修論の値を再現します。シードが壁に触れると、成長するバイオフィルムが固定面と接するところで応力が集中し、
恒常圧の則がそれを消します。慶應での次の段階は、本物の歯の表面と較正したパラメータです。
