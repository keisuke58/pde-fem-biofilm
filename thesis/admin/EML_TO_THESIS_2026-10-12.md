# EML 論文で直したこと → 修論に入れる（2026-10-12、クラウドから）

修論の第 3 章（`chapters/ch3_heine.tex`、Heine の 5 種推定）は、**10 月はじめの BMB 原稿の時点の内容**のまま。
その後、EML 論文（Extreme Mechanics Letters 特集号、Tmcmc202601 `docs/revision/EML_submission/`）で推定をやり直し、書き方も直した。
修論もこれに合わせる。**数値は EML 原稿が正**（DS・DH は最終段 (iv) がまだ計算中で暫定）。

- EML 原稿: `keisuke58/Tmcmc202601`、ブランチ `claude/gate-off-map-check`、`docs/revision/EML_submission/EML_manuscript_nishioka.tex`（本文）・`EML_supplementary.tex`（補足）
- EML のタイトル（確定）: *Identifiability of species interactions in a living material: Bayesian calibration of an extended Hamilton principle multiphase model*
- 共著者レビュー（9/25、Szymon・Rumjhum）への回答書: 同じフォルダの `response_to_coauthors.tex`。**修論でも同じ指摘を受けないように**、下の対応はそこからも来ている。

---

## 1. 枠組み・言い方（第 3 章の冒頭、第 1 章 §1.x「Identifying the interactions」、第 6 章）

| 今の修論 | EML での形 | 理由 |
|---|---|---|
| "Bayesian identification of the species interaction matrix" / "full-discovery strategy"（ch3 l.4–5, l.90） | 相互作用は**モデルの実効パラメータ**。「データで何が決まり、何が決まらないか」を成分ごとに言う | 共著者レビュー（一般指摘・C1–C2）: 実験で相互作用を決めたように読める書き方を避ける |
| 生物学的な解釈を「示した」 | 文献と**矛盾しない**、と書く。機構の証明とは書かない | C1–C2・C8 |
| 相互作用の向き（Vei→Pg など、ch3 l.101–103 の表） | **対称行列なので向きはない**。Vei–Pg のように無向で書く。文献の向きは「解釈の参考」として別に示す | A4・C14 |
| "physics-constrained" | 使わない（PINN と誤解される）。"Bayesian calibration of a model derived from the extended Hamilton principle" | タイトル決定時のユーザー判断 |

## 2. 推定の段（ch3 §「Two-phase fix-ψ strategy」l.134–150、表 l.224–246、収束 l.263）

**Phase 1 / Phase 2 の 2 段の書き方をやめて、4 段にする。**

| EML | 中身 |
|---|---|
| (i) pilot | 組成だけ、広い箱 |
| (ii) p1 | (i) の事後平均 ±4σ に狭めた箱 |
| (iii) p2 | 生存率（viability）の通路を足し、**箱は (i) に戻す**。事後の質量が箱の端（外側 5% に 20% 以上、片側）にたまった成分は、その端を広げてやり直す |
| (iv) ult | (iii) の事後 ±4σ |

- **どの段も 3 seed で回し、判定 (a)–(e) を通ったものだけ使う**: (a) logL の再計算一致、(b) β=1 まで到達、(c) 粒子あたり・次元あたりの移動回数、(d) seed 間の max logL の幅 ≤ 1、(e) seed 間の中央値の差 ≤ 0.5 事後 sd（同定されない成分は外す・多峰は山ごと）。
- 「fix-ψ は upper bound」「Phase 1 vs Phase 2 の MAP 比較で検証」の論法は消す。
- 粒子数は条件・段で違う（N_p = 1000〜20000）。EML 補足の表 S10（段ごとのコスト）を使う。
- **ゲート（Hill のゲート、K_hill）は使っていない**（K_hill = 0）。修論にゲートの記述があれば外す（修論の中での説明は可。共著者向けメールでは書かない約束）。

## 3. 同定性（新しい節を足す。ch3 §Results の「Inferred interaction matrices」の前か中）

- **分類の規則**: 事後 sd が「箱の一様分布の sd」の 0.8 倍未満 → **同定された**。そうでなく、90% 区間が 0 を含まない → **符号だけ決まる（s.o.）**。どちらでもない → **同定されない（n.i.）**。
- **結果**: 15 成分のうち同定されたのは、CS 4・CH 6・DS 8・DH 9（DS・DH は暫定）。共生（commensal）の 2 条件は初期定着菌の間だけで、重なるところは符号が一致。
- **独立な確認（Fisher 情報量）**: MAP での尤度の曲率で決まる方向の数は 4・4・8・5（同じ順）。事後の分類とおおむね合う。
- EML の図 3（同定性マップ）・図 4（箱の中での位置で描いた周辺事後）・補足表 S6 をそのまま使える。
- **合成データでの回収テスト**（補足に入れる予定、GPU で計算待ち）: 結果が出たら一文足す。

## 4. 時間の対応・初期値（ch3 l.47 付近「Day 1 is consumed as the initial condition」）

EML 本文の書き方に合わせる:
- モデルの時刻は**接種からの時間**、t = s·d、s = 0.2375/21。
- 1 日目の組成（最初の測定）を t = 0 の初期状態に使う。
- HOBIC（流れの反応器）では、24 h 後に接種した培地を無菌の培地に替える（Heine et al. 2025）ので、成長の始まりを 1 日目と見ることもできる。
- 1 日目から始めた場合（t = s(d−1)）との違いは、MAP の組成誤差で最大 0.008（再推定なし）。

## 5. pH（ch3 l.49「held-out validation targets」、§heine_validation、l.304 の表「R² = 0.778 (held-out)」、ch6 l.86「R² = 0.78」）

**「検証（validation）」ではなく「確認（check）」と書く。数値も差し替え。**
- pH は推定に使っていない（λ_pH = 0）。
- 実測の組成 → 実測の pH の線形回帰を**データだけで 1 回**決めて固定（切片＋So・An・Fn・Pg、Vei 基準）: pH = 6.367 − 0.095 So + 0.208 An + 0.695 Fn + 2.484 Pg、**R² = 0.93**（CH・DH の 90 点）。
- それを事後予測の組成に当てる: **DH は R² = 0.63・RMSE 0.12（暫定）**。**CH は約 +0.13 ずれる**（モデルが Fn・Pg を 4–5% に保つのに、データは 1% 未満）。
- 旧版の R² = 0.778 は使わない。

## 6. 当てはまり（ch3 §Quantitative fit、表 tab:heine_rmse・tab:heine_r2）

- RMSE は **最終段 (iv) の MAP の組成誤差**: CS 0.179、CH 0.076、DS・DH は暫定（EML 本文と補足の値）。
- **CS の 0.179 は限界として書く**: 誤差はほぼ So と Vei（0.31・0.24、ほかは 0.05 未満）。So が 1 日目 0.18 → 3 日目 0.70 と急に増える。組成と生存率を同時には合わせられず、事後は妥協点に落ちる（チャネルの重みを条件ごとに変えて合わせることはしていない）。
- MAE は載せない（A2）。R² の説明が逆になっていた箇所に注意（A3）。
- 「ln Ẑ ≤ max ln L」を全 run で確認済み（A1）と一言。

## 7. 消す・小さくするもの（EML で外したもの）

- UMAP・相関の表・条件間の転移（cross-condition transfer、l.284–285、l.434–435）: EML では外した。**修論に残すなら「補助的」と明記**し、検証のようには書かない（C5・C9–C10）。
- 「HMC/NUTS は GPU で失敗する」（l.164）: 主張しない。速度比は自分の実装の実測（約 150×）だけ（C11）。
- gLV との比較表（l.277–313）の「pH prediction R² = 0.778 (held-out)」の行は、上の 5 に合わせて直す。

## 8. 限界の節（ch3 の Discussion または ch6）

EML の Discussion の限界をそのまま使える:
- 独立なガウス誤差（組成データなので区間がやや狭いかもしれない、B4）。
- 相互作用係数は実効量（代謝物の交換や pH などをまとめたもの）。
- 共生・ディスバイオシスの比較は、菌株（Veillonella の種、Pg の株）と培養の違いも混ざる（C15）。静置 vs 流れの比較は同じ群集の中なので混ざらない。
- 尤度は割合だけなので、バイオマスの絶対量はモデルの予測で、FEM との結合ではその事後のばらつきを伝える必要がある（第 4 章との接続）。

## 9. 数値・図の入手先

| もの | 場所（Tmcmc202601） |
|---|---|
| 最終 MAP（CS・CH 確定、DS・DH 待ち） | `docs/revision/generated/final_theta_MAP/` |
| 全 run の粒子 | `docs/revision/generated/final_runs/` |
| 図（暫定版） | `docs/revision/EML_submission/figures/provisional_2026-10-11/` |
| 本文の数値 | `docs/revision/generated/paper_numbers.json` |

**DS・DH は最終段が通ったら数値・図が変わる。** 修論に入れるのは、それを待ってからでもよい（章の構成と書き方だけ先に直す）。
