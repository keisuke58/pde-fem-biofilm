# 修論用に描き直した図 3.6（2026-10-19、EML 側セッション）

図 3.2（`eml_fit.pdf`）はユーザーの判断で**前の版のまま**（描き直し版は使わない）。

幅 16 cm（本文幅）で作ってある。`\includegraphics[width=\textwidth]` なら拡大縮小なし（文字は印刷時に目盛り 8 pt・軸ラベル 9 pt・列見出し 10 pt、Times = Nimbus Roman を埋め込み）。
上のタイトルとパネル内の凡例は外し、凡例は図の下に 1 つ。縦軸は 0〜1、目盛り 0・0.5・1。
図 3.2 と同じく、データは反復の箱ひげ（四分位）＋中央値の点。Day 1 の点は初期値（φ_init を正規化）で、線と帯は Day 1 から（横軸も図 3.2 と同じ 1・3・6・10・15・21）。

| ファイル | 大きさ | 元（Tmcmc202601） |
|---|---|---|
| `figures/eml_seed_maps.pdf`（図 3.6） | 160 × 180 mm | `tools/plot_seed_maps.py --thesis` |

## 図 3.6 のキャプション用: seed ごとの max logL と組成 RMSE

- max logL は最終段の粒子の最大値（polish 前）。RMSE は seed の MAP の生きた菌の組成（φψ/Σφψ）と、データの中央値の差。
- 段は CS・CH・DS が段 (iv)、DH は段 (iii)（段 (iv) は seed の一致の基準 (e) を満たさないため）。

| 条件 | seed 42 | seed 123 | seed 7 | 事後の中央値の RMSE |
|---|---|---|---|---|
| CS (stage iv) | −32.91 / 0.092 | −32.86 / 0.090 | −33.08 / 0.093 | 0.094 |
| CH (stage iv) | −26.84 / 0.060 | −26.52 / 0.060 | −27.11 / 0.066 | 0.086 |
| DS (stage iv) | −96.80 / 0.057 | −96.84 / 0.057 | −96.94 / 0.057 | 0.057 |
| DH (stage iii) | −51.41 / 0.073 | −51.54 / 0.076 | −51.77 / 0.073 | 0.083 |

（各セル = max logL / RMSE）

キャプション例（英語）:
> Viable species fractions of the MAP estimate of each seed (coloured lines; 42, 123, 7), 10–90% band (grey) and median (dashed) of the posterior prediction pooled over the three seeds, and measured viable composition (dots: median, boxes: replicates, interquartile range; Day 1 is the initial condition). Stage (iv) for CS, CH and DS; stage (iii) for DH. The seeds' maximum log-likelihoods and composition RMSEs are listed in Table X / given in the text.


## 第 3 章の残りの図（2026-10-19 追加、同じ基準）

すべて幅 160 mm、Times、文字は印刷時に 8 pt 以上（目盛り 8・軸ラベル 9・列の見出し 10 pt）、上のタイトルなし、凡例は図の下に 1 つ。
元: Tmcmc202601 の `tools/thesis_figs_ch3.py`（入力は `make_eml_figs56.py` が書く `_thesis_inputs.npz` と `figures/raw/recovery.json`）、
パイプラインは `docs/revision/EML_submission/figures/fig1_pipeline_thesis.tex`。

| ファイル | 大きさ | 変えたところ | キャプションで直すこと |
|---|---|---|---|
| `eml_pipeline.pdf` | 160 × 86 mm | 3 段に並べ替え（モデル・データ → 4 段の推定 → 受け入れの条件）。中身は同じ | なし |
| `eml_ident.pdf` | 160 × 125 mm | 4 条件を 2 × 2 に。各条件の左が行列 (a)、右がネットワーク (b)。色の目盛りは下に横向き | 「(a) 各条件の左の行列、(b) 右のネットワーク」と分かるように。条件の見出しは「CS (stage iv)」など |
| `eml_posteriors.pdf` | 160 × 110 mm | 見出しは「CS (stage iv)」＋その下に「N identified」 | なし |
| `eml_recovery.pdf` | 160 × 112 mm | 条件の順を CS・CH・DS・DH に。凡例に「† not covered」 | 各パネルの下の「covered N/15」はそのまま |
| `eml_checks.pdf` | 160 × 70 mm | 見出しを短く（「DH: knockout」「CH: pH (stage iv)」「DH: pH (stage iii)」）。凡例は 3 パネル共通で下に 1 つ。pH の点を小さく | なし |
| `eml_cs_tradeoff.pdf` | 160 × 68 mm | 幅 16 cm に。凡例は下 | **`\includegraphics[width=0.5\textwidth]` → `width=\textwidth` に** |

図 3.2（`eml_fit.pdf`）: **前の版の配置のまま**（見出しに条件名と RMSE、凡例は上、箱ひげ）、幅 160 mm × 195 mm・Times・8 pt 以上にした（2026-10-19、ユーザー決定）。
見出しの RMSE（CS 0.090・CH 0.060・DS 0.057・DH 0.073）は各 seed の最良の粒子の値で、EML の本文の「0.090–0.093」などと同じ。磨いた MAP の値（CS 0.087・CH 0.059）とは少し違うので、キャプション・本文の値をそろえること。


## Latin Modern 版（`*_lm.pdf`、2026-10-19 夜）

修論の本文（`main.tex` は `lmodern`）とほかの図 27 枚は **Latin Modern** なので、EML の図 8 枚を Latin Modern でも描いた。
Times 版（今までのファイル名）も残す（ユーザー決定: どちらも必要）。中身・大きさ・配置・文字の大きさは同じで、字体だけ違う。

| Times 版 | Latin Modern 版 |
|---|---|
| `eml_pipeline.pdf` | `eml_pipeline_lm.pdf` |
| `eml_fit.pdf` | `eml_fit_lm.pdf` |
| `eml_seed_maps.pdf` | `eml_seed_maps_lm.pdf` |
| `eml_ident.pdf` | `eml_ident_lm.pdf` |
| `eml_posteriors.pdf` | `eml_posteriors_lm.pdf` |
| `eml_recovery.pdf` | `eml_recovery_lm.pdf` |
| `eml_checks.pdf` | `eml_checks_lm.pdf` |
| `eml_cs_tradeoff.pdf` | `eml_cs_tradeoff_lm.pdf` |

本文と字体をそろえるなら、`ch3_heine.tex` の `\includegraphics` を `_lm` の方に変える。
