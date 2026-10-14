# P3 の骨組み（慶應での続き、10月15日の案）

仮題：*Growth from the living biomass couples species composition to stress in oral biofilms*
（慶應の修論の仮題 *From Species Composition to Stress* の中心の章になる）。
修論（LUH）には入れない（10月15日の決定、`CLAUDE.md`）。文献は `P3_LITERATURE.ja.md`、材料の一覧は下の表。

## 言いたいこと（一文）

成長をバイオフィルムの量ではなく生きている量 Σφ̄_i = Σφ_iψ_i で決めると、組成と溶菌（抗菌薬）が応力に効くようになり、
抗菌薬の実験で二つの成長則を見分けられる。

## 節の案

| 節 | 中身 | 材料（今あるもの） | 足りないもの |
|---|---|---|---|
| 1. はじめに | 組成と応力をつなぐ問い。Klempt 2024 の成長は量だけで決まる | 修論の第 1 章 | 歯・インプラントでの意味づけ（`P3_LITERATURE` §5） |
| 2. モデル | 場（Klempt 2024）、点モデル（Klempt 2026、5 菌種、較正済み）、成長則 (a) α̇ = k_α φ と (c) α̇ = κΣφ̄_i。(c) の根拠：Wanner–Gujer 1986、Alpkvist–Klapper 2007（活性な菌だけが体積を増やす）、Feng 2021 | 修論の第 2、4 章、`P3_LITERATURE` §1 | (c) の変分原理からの位置づけ（ハミルトン原理に入るか） |
| 3. 時間と速さの決め方 | T* と日数の対応、体積データ（Heine et al.）で速さを決める | `RESEARCH_IDEAS` §10、`volume_figs.py` | 剥がれ・縮みを入れるか |
| 4. 検証 | 点モデル単独との一致（4e−12）、メッシュ（8³、16³、24³）、κ への比例 | `IKMHIWI03_TO_CLOUD.ja.md`（V2、A-2、A-4） | 24³（B の列） |
| 5. 結果 1：組成と応力 | 4 条件（CS、CH、DS、DH）で (a) と (c) の応力 | CH、CS の 8³・16³ | DS・DH（EML の MAP 待ち） |
| 6. 結果 2：抗菌薬 | T* = 0 と 0.5 からの抗菌薬、菌種ごとの感受性、体積の予測 | A-1、A-6、`volume_prediction_antibiotic.py` | DH での抗菌薬（MAP 待ち） |
| 7. 結果 3：大きい成長 | 体積に合わせた κ で Day 21 まで、E = 600 Pa での応力の大きさ | A-5（T* 0.35 で停止） | A-8 の結果 |
| 8. 歯とインプラント | 形の上での (a) と (c) の違い | `keio_wp2/`（Abaqus） | (c) を Abaqus の UMAT に入れる |
| 9. 考察 | どちらが正しいかを決める実験、両方向の結合（応力で成長が遅くなる、文献の圧力は 0.1〜1 MPa の桁） | 会議用まとめ、`RESEARCH_IDEAS` §11、§12 | A-7 の結果 |

## 順番

1. 会議（10月19日）で Oliver の意見を聞く（φ₀ の意味、(c) の扱い）。
2. IKMHIWI03 で取れるものを 12 月中旬までに取る（B の列、A-6〜A-8、DS・DH が出たら 4 条件）。そのあと慶應には半年ほど ANSYS がない。
3. 慶應で：(c) を Abaqus の UMAT に入れ、歯・インプラントの形で回す（節 8）。
4. 2027 年春〜夏に投稿（`KEIO_PLAN.ja.md` の予定どおり）。

## 決めること

- 雑誌（Biomech. Model. Mechanobiol.、J. Mech. Behav. Biomed. Mater.、J. R. Soc. Interface など）。Klempt 2024 と同じ BMMB が第一候補。
- 共著者：Meisam（点モデルと較正）、Oliver（要素）、Heine（体積・組成のデータ）。Felix の 2 菌種の理論は使わない。
- 点モデルの記号は Klempt et al. 2026 のまま（φ_i、ψ_i、φ̄_i）。
