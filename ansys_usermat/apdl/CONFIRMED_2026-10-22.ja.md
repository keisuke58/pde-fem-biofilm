# 確定事項（2026-10-22、生きている割合の結合、ANSYS と Python の照合）

まずここを読む。根拠は `IKMHIWI03_TO_CLOUD.ja.md` の A-xx と `CLOUD_TO_IKMHIWI03_2026-10-20.ja.md` の追記 xx。
結果のどれを使うかは `results/2026-10-living/README.md`（ANSYS）と `results/2026-10-3d/README_coupled_living.md`（Python）。
モデルの式は 10 月 21 日に凍結（CLAUDE.md）。下はすべて**数値解法と実装**の話で、式は変えていない。

## 1. 確定したこと

| # | 内容 | 根拠 |
|---|---|---|
| 1 | ANSYS の前線項は −r c/(K + c) n_c·∇φ（Ori × Growth、FRONT_GRAD = 1）で Python の式と同じ | A-38 §1 |
| 2 | 前線項を切ると 8³ の分離例で ANSYS と Python が合う（種 Σφ̄ 0.971 対 0.970）：拡散、消費、点モデルは正しい | A-37 §3 |
| 3 | ANSYS の φ の 1 歩は、ANSYS 自身の NEM の行列から Python で計算し直すと 8 桁で一致する（ラプラシアン、勾配、更新） | A-42 §1 |
| 4 | Python の `--artdiff`（中心差分の前線項）は壁のマスで片側差分の勾配を使っていた：**誤り**。角の塊が φ = 1.1 に張り付く | A-38 §3 |
| 5 | ANSYS は表面点の φ を逆距離二乗の平均で決めていて、壁ぎわの値が低くなり φ が壁から漏れていた：**誤り** | A-39 §2 |
| 6 | NEM のラプラシアンは行の和 0、列の和が 0 でない（壁ぎわ正、内部負）：等しい積分点の重みで質量を保たない | A-42 §2 |
| 7 | 対策 SURF_NEAREST（表面点 = 最も近い積分点）で漏れ 1/4、LAP_FILE（保存形のラプラシアン）で 0 | A-39 §3、A-42 §3 |
| 8 | NEM の勾配は壁のマスで Python の鏡映の差分とほぼ同じ（比 0.98〜1.24） | A-39 §2 |
| 9 | 時間刻み：外側の増分 0.005 → 0.000625 で α−1 は 3 % 以内、0.00125 と 0.000625 は 0.4 %：収束している | A-39 §1 |
| 10 | 栄養場の初期値（供給層だけ 1、ほか 0）は結果を変えない（NUT1_INIT_ALL で 3 桁同じ） | A-38 §4 |
| 11 | dc/dt（d = 1e8 mm²/T\*）と準定常は同じ（ANSYS 3 桁、Python 1e−4 以下） | A-34、A-38、追記 30 |
| 12 | 人工拡散 r h/4 の陽解法は 8³ で外側の増分 0.005 以下、16³・24³ で 0.00125 が要る | A-34、A-35 |
| 13 | 8³ の分離例、両方の修正後：ANSYS の種 α−1 3.46e−4、Python 3.27e−4（壁 0）〜3.80e−4（鏡映）。混合例は 2.92e−4 対 3.09e−4 | A-42 §3 |
| 14 | 残りの差は離散化の違い（ANSYS：8³ 要素 × 8 積分点 = 16³ 点の MLS、30 近傍；Python：マスの 7 点差分）。格子を細かくして比べる | A-42 §1 |
| 15 | 16³ では Python の壁の扱い（0 か鏡映か）は `--artdiff` の結果にほとんど効かない（α−1 6.32 対 6.27e−4） | A-40 前の手元計算 |
| 16 | Python の本番（`--fast --tvd`）は壁の勾配を鏡映にした（クラウドの判断、16³ で応力 10 %、32³ で 7 % 変わる） | 追記 30、クラウドの 22 日の報告 |

## 2. 本番の設定（22 日から）

- **ANSYS**：ビルド `F:\biofilm_upf_living7`（living5 ＋ NUT1_INIT_ALL、FRONT_DUMP、SURF_NEAREST、NEM_DUMP、LAP_FILE；
  差分は `oliver_patches/living6_USolBeg.patch`、`living7_Ussfin.patch`）。デッキ：
  `make_living_deck.py --front 1 --v2 --start-step 0.0025 --dt-ussfin 0.1 --cons 1 --cons-ab 1 --phi-min 0.001 --artdiff N
  --deltim 0.00125 --transient-c 1e8 --nut1-init-all --surf-nearest --lap-file`（＋例ごとの `--layout`、`--ab --t-ab 0.333333`、`--eco`）。
  作業フォルダに**同じ格子の** `lap_cons.txt`（`checks_a38_a39/build_lap_cons_a42.py`、短い `--nem-dump` の本から）。名前の頭 `u7_`。
- **Python**：本番は `--fast --tvd --wall-grad mirror`（タグ `_tvd_wgmirror_fast`、128³）。ANSYS との照合は
  `--artdiff --dtc <ANSYS の増分>`、壁は `--wall-grad zero`（`--artdiff` の既定）または `mirror` を明記。

## 3. 使わないもの（間違えやすい順）

| もの | なぜ | 代わり |
|---|---|---|
| ANSYS `u4_`（A-35 の 16³、A-36 の 24³） | 表面点から φ が漏れる | `u7_`（16³ は A-43、24³ は A-41）。`u4_16` は「表面点の感度」としてだけ |
| ANSYS `u3a_`、`u3a5_`、`u3t_`、`u6_`、`u6t_` | 8³ の照合用、漏れあり | 照合の経過としてだけ |
| ANSYS `u6s_`、`u5_` | SURF_NEAREST だけ（列の和が残る） | `u7_`（`u5_16_sep_ab` は追記 30 の u4s_ にあたる中間） |
| ANSYS `u3_` | 安定化なしで 8³・16³ が壊れる（A-32） | — |
| ANSYS `u1_` | 論文の値になる前の設定（g l²/d = 4 など、README の表）、漏れあり | `u7_` |
| ANSYS `lv8_`、`v2_`、`v2s_`、`sens_`、`i8*`、`tsens_`、`v2f_`、`fsens_`、`fthr_`、`v2t_`、`dbg_` | L が凍っていた（A-22）か、φ が 10 倍速く進んでいた（A-25） | — |
| ANSYS `sym_*`、`symt_` | 結果ではない：近傍と対称性の試験（A-24）としてだけ有効 | — |
| Python `_artdiff`（`_wg` なし） | 壁が片側差分（誤り） | `_artdiff_wgzero` / `_artdiff_wgmirror` |
| Python `_tvd_fast`（`_wg` なし） | 栄養の向きの勾配が壁で片側差分 | `_tvd_wgmirror_fast` |
| 「Python は種を 1 近くに保つ」（追記 29 の文） | 誤り：Python も種は痩せる。1.1 は壁の角の塊（片側差分のせい） | A-37、A-38 |
| 「違いは時間刻み」（クラウドの仮説） | 否定：0.4〜3 % | A-39 §1 |

## 4. まだ確かめていないこと

- 16³ と 24³ で両方の修正後に ANSYS と Python がどこまで合うか（A-43、A-41、自動）。
- 漏れのもとの列の和が SURF_NEAREST で ±22 残る理由（辺と角の表面点か、NEM の近傍の非対称）。保存形にすれば効かない。
- オリバーへの連絡（表面点の扱いと保存形のラプラシアン）。まだ送っていない。

## 5. ビルドの一覧（IKMHIWI03、`F:\`）

| ビルド | 中身 | 使い道 |
|---|---|---|
| `biofilm_upf_living4` | 生きている割合を毎回読む修正（A-22） | `u1_`、`u3_`、`u3a_` |
| `biofilm_upf_living5`（＋`5s`） | ＋ NUT_TRANSIENT（dc/dt） | `u4_`、`u3t_` |
| `biofilm_upf_living6`（＋`6s`） | ＋ NUT1_INIT_ALL、FRONT_DUMP、SURF_NEAREST | `u6*`、`u5_` |
| **`biofilm_upf_living7`**（＋`7s`、`7_16a〜e`、`7_16n`） | ＋ NEM_DUMP、LAP_FILE | **`u7_`（本番）**、`u7c_`、`u7n_` |

スイッチはどれも既定 0（従来どおり）。新しいフォルダを作るときは `.f/.F/.inc/.h/.mod/.obj/.lrf/.def/.manifest/.dll/.mac`、`ANSYS.exe`、
`eco_*.txt`、`base_*.dat` だけ写す（living6 には 7 GB の `comp_trace.csv` がある）。
