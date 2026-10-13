# EML → 修論: Heine の体積・生きた割合（実測値）を FEM の検証に（2026-10-13、クラウドから）

Heine（MHH）から 10/13 に LIVE/DEAD の画像ごとの値（生きた割合・バイオフィルム体積）が届いた。
EML の推定（0D の点モデル）では体積を使えない（モデルの Σφ は 0.91〜0.98 でほぼ一定、実測の体積は 2〜3 倍変わる）。
体積の変化は **成長の変数 α（F_g = (1+α)I）** の側に対応するので、**修論の FEM（成長・応力）の検証**に使う。

## 置いたもの

| ファイル | 中身 | git |
|---|---|---|
| `data/heine_volume/volume_summary.csv` | 条件 × 日: 画像数・バイオフィルム数、体積の中央値・四分位（画像）、バイオフィルム平均の平均と SD | 入れる（集計） |
| `data/heine_volume/live_fraction_summary.csv` | 条件 × 日: 生きた割合のバイオフィルム平均の平均と SD | 入れる（集計） |
| `data/heine_volume/live_volume_species.csv` | 条件 × 日 × 種: 種ごとの生きた体積（体積 × 生きた割合 × 生きた菌の組成）の中央値と 90% 区間 | 入れる（集計） |
| `data/raw_heine/biofilm_level.csv` | 条件 × 日 × バイオフィルム: 体積と生きた割合のバイオフィルム平均 | **入れない**（gitignore） |

- 体積の単位は **10^6 µm³**（画像 1 枚の視野の体積）。**5 枚 = 1 バイオフィルム**（N. Heine）、1 時点あたり 3〜9 個。
- CH の Day 1 は 40 倍対物の値を 63 倍に換算した値。
- 作り方: Tmcmc202601 の `tools/export_volume_summary.py`（`tools/raw_data_report.py` のブートストラップ、バイオフィルム単位）。

## 公開の決まり

- データは「**請求に応じて提供**」（N. Heine の希望、10/13）。**このリポジトリは public** なので、反復ごと・画像ごとの値は入れない。
  上の 3 つの集計（論文の Fig. 2A の箱ひげ図と同程度）は入れてよい。
- `data/raw_heine/` は gitignore。中身はこのコンテナにしか無い。別のマシンで要るときは、Tmcmc202601 の
  `docs/handoff/2026-10-13i_raw_data_fetch.md` の手順で生データを置き、`tools/convert_raw_heine.py` → `tools/export_volume_summary.py`。
- **注意**: `data/heine_species_distribution_biofilm.xlsx`（qPCR の生データ、10/12 から入っている）も同じ扱いのはず。
  private にするか履歴から消すかはユーザーの判断待ち。

## FEM での使い方（案）

- FEM は 1 つの代表的なバイオフィルムを解くので、使う量は **日ごとの平均とバイオフィルム間の SD**（画像ごとのばらつきは使わない）。
- 成長の検証: 較正済みの組成・生存率で駆動した α(t)（生きた体積 Σφψ で駆動、EML 補足 S_growth と同じ）から出る体積の変化
  V(t)/V(1) = (1+α(t))^3 と、実測の `mean_biofilm` の比を比べる。k_α を合わせるなら SD を重みに。
- 条件ごとの違い: 実測では DS・DH は Day 15〜21 で増え、CH は Day 3 で一度減ってから戻る。静置（CS・DS）と流れ（CH・DH）で
  形が違う → 成長則が条件差を出せるかの確認になる。
- 種ごとの生きた体積（`live_volume_species.csv`）は、FEM の材料点に種ごとの体積を与えるときの目安に使える（推定には使っていない）。
