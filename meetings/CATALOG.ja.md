# スライドと図のカタログ（25 Nov 2026）

これまでに作ったスライドと図を，再利用しやすいように並べたものです．
ファイルの中身は変えていません．調べた日は 25 Nov 2026 です．

目次：
1. スライド一式（どれが何で，どう作るか）
2. 再利用しやすいスライド（1 枚ずつ）
3. 図（`assets/`）と作るスクリプト
4. 再利用するときの決まり

**状態の書き方**

- 送付済み：相手に送った記録がある．
- 発表済み：打ち合わせで見せた記録がある．「（推定）」は，その日の打ち合わせ用に作ったが，見せた記録までは見つからなかったもの．
- 未送付：まだ渡していない．
- 温存：今は出さない．理由を書いてある．
- 不明：記録から分からない．

**温存するもの（オリバーにはまだ渡さない；ユーザー，4〜9 Nov と 25 Nov）**

- 歯頸部の窓（実際の歯の形）
- Klempt 2024 式 30 の力学の項 $-K\phi_{3D}(\mathbf I:\mathbf C_e-3)$ と K の掃引
- 非対称の初期配置（修論と課題研究報告には入っている．オリバーへのスライドには入れない）
- 四種の例（Klempt 2026 case 4；2027 年 6 月ごろに出す．修論にも入れない）
- 時間で振動する栄養の供給（25 Nov）
- 25 Nov に足した「ほかに試した設定」の図（3.9）：まだ誰にも見せていない．出す前にユーザーに確かめる

---

## 1. スライド一式

### 1.1 pptxgenjs のスライド（`meetings/`）

どれも `node` と pptxgenjs で作る．pptxgenjs は `meetings/2026-11_coupled_examples/build/node_modules` にある．
出力先の `build/` は git に入っていない．

| スライド | 言語・枚数・ノート | 相手と日 | 状態 |
|---|---|---|---|
| `meetings/2026-11_coupled_examples/deck_en.js` | 英語，20 枚（`WITH_SINE=1` で 21 枚），ノートあり（EN と JA） | オリバー，Meisam，Felix．11 月 | 6 Nov に 18 枚の版を送る予定だった（実際に送ったかは不明）．25 Nov に 2 枚とノートを追加した今の版は未送付．付録 C（振動する栄養）は**温存** |
| `meetings/2026-11_coupled_examples/deck_ja.js`（`separate`） | 日本語，12 枚，ノートあり | 村松先生，2 Nov | 発表済み（README：分離の方は 2 Nov に発表） |
| `meetings/2026-11_coupled_examples/deck_ja.js`（`mixed`） | 日本語，14 枚，ノートあり | 村松先生，30 Nov | 未送付（30 Nov に発表する予定） |
| `meetings/2026-11_coupled_examples/deck_geom_en.js` | 英語，10 枚，ノートなし | オリバー，Meisam，Felix | 6 Nov に送る予定だった（不明）．README では 12 月初めに会って渡す予定 |
| `meetings/2026-11_coupled_examples/deck_k26c4_en.js` | 英語，10 枚，ノートなし | 未定（2027 年 6 月ごろ） | **温存**（四種の例） |
| `meetings/2026-11-25_oliver/deck_oliver_en.js` | 英語，7 枚，ノートあり（EN と JA） | オリバー，25 Nov | 25 Nov の打ち合わせ用．見せたかは不明 |
| `meetings/2027-03_kadaikenkyu/deck_kadai_ja.js` | 日本語，9 枚，ノートに原稿 | 慶應の課題研究（5 分の録画），締切 12 Mar 2027 | 未送付．7 枚目に非対称の初期配置の図（オリバーには**温存**） |
| `meetings/keio_deformation/deck_deform_ja.js` | 日本語，11 枚，ノートあり | 慶應の続き（村松先生） | **温存**（式 30 の力学の項，K の掃引） |
| `meetings/keio_living_growth/deck_living_ja.js` | 日本語，7 枚，ノートあり | 慶應の続き（村松先生） | 未送付．修論には入れない内容（生きている菌からの成長）．6 枚目の「次にやること」に歯頸部の窓の話が出る |

作り方と出力：

| スライド | コマンド | 出力 |
|---|---|---|
| `deck_en.js` | `bash meetings/2026-11_coupled_examples/build.sh en` | `meetings/2026-11_coupled_examples/build/coupled_examples_en.pptx` |
| `deck_en.js`（付録 C つき） | `WITH_SINE=1 bash meetings/2026-11_coupled_examples/build.sh en` | 同じファイル（21 枚） |
| `deck_ja.js` | `bash meetings/2026-11_coupled_examples/build.sh ja` | `build/report_separate_ab_academic.pptx`，`build/report_mixed_ab_academic.pptx` |
| `deck_geom_en.js` | `bash meetings/2026-11_coupled_examples/build.sh geom` | `build/geometries_en.pptx` |
| `deck_k26c4_en.js` | `bash meetings/2026-11_coupled_examples/build.sh k26c4` | `build/four_species_en.pptx` |
| 上の 4 つ全部 | `bash meetings/2026-11_coupled_examples/build.sh all` | 上と同じ |
| `deck_oliver_en.js` | `cd meetings/2026-11-25_oliver && NODE_PATH=../2026-11_coupled_examples/build/node_modules node deck_oliver_en.js` | `meetings/2026-11-25_oliver/build/oliver_1125_en.pptx` |
| `deck_kadai_ja.js` | `cd meetings/2027-03_kadaikenkyu && NODE_PATH=../2026-11_coupled_examples/build/node_modules node deck_kadai_ja.js` | `meetings/2027-03_kadaikenkyu/build/kadaikenkyu_ja.pptx` |
| `deck_deform_ja.js` | `bash meetings/keio_deformation/build.sh` | `meetings/keio_deformation/build/deform_ja.pptx` |
| `deck_living_ja.js` | `cd meetings/keio_living_growth && NODE_PATH=../2026-11_coupled_examples/build/node_modules node deck_living_ja.js` | `meetings/keio_living_growth/build/living_growth_ja.pptx` |

メモ：

- リポジトリの一番上から 1 行で作るときは `NODE_PATH=meetings/2026-11_coupled_examples/build/node_modules node <deck>.js`．
- `build.sh` は式の画像（`eq/*.tex` → `build/eq/*.png`）と図の PNG も作る．pdflatex，pdftoppm，node が要る．
  `node_modules` がなければ `build.sh` が `npm install pptxgenjs` をする．
- `deck_kadai_ja.js` は `build/eml_fit.png`，`build/eml_ident.png`（`keio_report/` の `eml_*.pdf` を `pdftoppm -r 220 -png -singlefile` で PNG にしたもの）と，
  `meetings/2026-11_coupled_examples/build/eq/eqs.png` を使う．先に `build.sh` を回しておく．
- `deck_kadai_ja.js` の冒頭のコメントは `NODE_PATH` なしの `node deck_kadai_ja.js` になっている．pptxgenjs が見つからないときは上の形で．
- `build/` に古い PDF（`coupled_examples_en.pdf` は 18 ページ）が残っている．今の PPTX とは枚数が違う．

### 1.2 LaTeX（beamer）のスライド

英語のものは `pdflatex <file>.tex` を 2 回，日本語のものは `lualatex <file>.tex` を 2 回．
IKMHIWI03 では `.\build_slides.ps1 -Tex <file>.tex [-MaxPages n]`（エラー，はみ出し，ページ数を調べる）．
出力は同じ名前の `.pdf`．PDF のページ数は題のページを含む（枚数 = フレーム数 + 1）．

| スライド | 言語・フレーム数・ノート | 相手と日 | 状態 |
|---|---|---|---|
| `slides_defense.tex` | 英語，38（本編と付録 A〜D で 20，Backup 1〜18） | コロキウム，16 Dec 2026 13:00，IKM 320 | 未発表．原稿は `speaker_notes_defense.tex`（EN と JA，`lualatex`） |
| `slides_qa_ja.tex` | 日本語，9 | 自分用（口頭試問の想定問答） | 自分用 |
| `slides_progress_1024.tex` | 英語，15 | オリバー，Meisam．10 月下旬 | 不明（4 Nov のメールの下書きに添付の候補として出てくる） |
| `slides_progress_2026-10.pptx` | 英語，15 枚，ノートあり | 同じ内容の PowerPoint 版 | 不明．作るスクリプトは見つからない |
| `slides_coupling_1020.tex` | 英語，9 | オリバー，Meisam．20 Oct の要点 | 発表済み（推定） |
| `slides_1020.tex` | 英語，30 | オリバー，Meisam．20 Oct（19 Oct から変更） | 発表済み（推定）．原稿は `speaker_notes_1019.tex` |
| `slides_3d_demo.tex` | 英語，15 | 20 Oct．3 次元の例 | 不明．原稿は `speaker_notes_3d_demo.md` |
| `slides_3d_demo.pptx` | 英語，11 枚，ノートあり（動画入り） | 同じ内容の PowerPoint 版 | 不明．`python make_slides_3d_demo_pptx.py` |
| `slides_1005.tex` | 英語，19 | 村松先生，Meisam．5 Oct | 発表済み（推定）．原稿は `speaker_notes_1005.tex`．付録 A〜C が付録の手本 |
| `slides_1001_oliver.tex` | 英語，11 | オリバー，1 Oct | 発表済み（推定）．古い内容（4 条件の話）が多い |
| `slides_1102_nem_ja.tex` | 日本語，19 | 村松先生，2 Nov（NEM の基礎） | 不明（2 Nov に使ったかは記録なし） |
| `meetings/2026-11-02_muramatsu/slides_1102_muramatsu.tex` | 日本語，11 | 村松先生，2 Nov（画面共有用，送らない） | 発表済み（推定） |
| `slides_muramatsu_ja.tex` | 日本語，9 | 村松先生，10 月（基本） | 不明．古い内容（4 条件，接着の修正）：そのままは使わない |
| `slides_kiso_ja.tex` | 日本語，5 | 自分用，村松研での説明 | 自分用 |
| `slides_usermat_ja.tex` | 日本語，6 | 自分用（USERMAT の中身） | 自分用 |
| `slides_volume_ja.tex` | 日本語，7 | 自分用（Heine らの体積データ，14 Oct） | 自分用．図は `python ansys_usermat/apdl/volume_figs.py` |

メモ：

- `.\build_slides.ps1 -All` は `slides_1001_oliver.tex`，`slides_1005.tex`，`slides_muramatsu_ja.tex` の 3 つを作る．
- `slides_1102_muramatsu.tex` は `meetings/2026-11-02_muramatsu/` の中で `lualatex` を回す（図のパスが `../../assets/`）．
- `slides_progress_1024.tex` の動画は `\movie` のリンク．送るときは `assets/anim_coupled_tvd_*_ab_inherit.mp4` を PDF の隣に置く．
- Abaqus は修論の内容ではない．`slides_muramatsu_ja.tex` と `slides_1102_nem_ja.tex` の Abaqus の話は「慶應での続き」としてだけ使う．

---

## 2. 再利用しやすいスライド

書き方：スライド（番号か題）｜何を見せるか｜使う図．

pptxgenjs のスライドの番号は，今の版（`WITH_SINE` なし）での順番．題で探す方が確実．
beamer は題で探す（`\begin{frame}{題}`）．

### pptxgenjs のスライドを新しいスライドに写す方法

- 1 枚は `{ const s = frame('題'); ... }` のブロック 1 つ．このブロックを新しい `.js` に写す．
- 一緒に写すもの：冒頭の `require`，`pres` と色の定数（`NAVY`，`INK`，`GREY`，`RULE`，`ACC`，`TM`；日本語版は `JP`），
  ヘルパーの `frame`，`bullets`，`cap`，`img`（`pngSize` も），`fit`（`deck_kadai_ja.js`，`deck_deform_ja.js`，`deck_living_ja.js` にある），
  表を使うなら `H`，`C`，`L`，動画を使うなら `vidSlide`（`deck_geom_en.js`）．最後の `pres.writeFile(...)` も．
- パスの定数（`R` = リポジトリの一番上，`B` = `build/`，`EQ` = `build/eq/`，`A` や `AS` = `assets/`）は新しい場所に合わせて直す．
- ノート：日本語版とオリバー用は各ブロックの中の `s.addNotes('...')`．
  `deck_en.js` は最後の `NOTES` 配列（1 枚につき `[英語, 日本語]`）から `addNotes` で入れる．枚数が合わないとエラーになるので，
  スライドを足したら `NOTES` にも 1 行足す．
- 式の画像は `meetings/2026-11_coupled_examples/eq/*.tex`（`eqs`，`point`，`vm`，`gci`，`pe`）と
  `meetings/keio_deformation/eq/*.tex`（`eq30`，`flow`）．それぞれの `build.sh` が PNG にする．

### 2.1 モデルと式

| スライド | 何を見せるか | 図 |
|---|---|---|
| `deck_en.js` 3「2. Governing equations」／`deck_ja.js`「2. 支配方程式」 | 式 34〜36 と連成（前線の項 $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$，$\dot c$ つきの式 35） | `eq/eqs.tex` → `build/eq/eqs.png` |
| `deck_en.js` 4「2. … the point model and the antibiotic」／`deck_ja.js`「2. 支配方程式（続き）：点モデルと抗菌薬の効き方」 | 点モデルの式と，抗菌薬がどこで効くか | `eq/point.tex` → `build/eq/point.png` |
| `deck_en.js` 6「3. Constants and numerics」／`deck_ja.js`「3. 材料定数と数値条件」 | 定数の表（論文の値と，この研究で決めた値の区別） | 表だけ |
| `deck_en.js`「Appendix A. Notation」 | 記号の一覧 | なし |
| `deck_en.js`「Appendix B. Definitions of the reported quantities」／`deck_ja.js`「付録 A. 評価量の定義」 | von Mises，GCI，セルのペクレ数の定義 | `eq/vm.tex`，`eq/gci.tex`，`eq/pe.tex` |
| `slides_defense.tex`「Approach」 | 連成の全体図と 3 つの結論 | `assets/fig_coupling_schematic_agreed.png` |
| `slides_defense.tex`「Both models from the extended Hamilton principle」 | 2 つのモデルが同じ原理から出ること | なし |
| `slides_defense.tex`「Appendix A: notation」「Appendix B: growth field and mechanics」「Appendix C: point model and one coupled substep」 | 付録の手本（記号，場と力学，点モデルと 1 ステップ） | なし |
| `slides_defense.tex`「Backup 18: what is solved where, in one step」 | 1 ステップで何をどこで解くか | なし |
| `slides_coupling_1020.tex`「The loop at one Gauss point」「Equations」「Time: one clock for both models」 | 20 Oct に合意した連成の要点 | なし |
| `slides_progress_1024.tex`「The coupling agreed on 20 October」 | 合意した連成の図 | `assets/fig_coupling_schematic_agreed.png` |
| `slides_1102_muramatsu.tex`「実装：外側の 1 ステップの流れ」 | 実装の流れ（日本語，1 枚） | `assets/fig_coupling_impl_flow_slide_ja.pdf` |
| `deck_kadai_ja.js` 5「4. 連成：点モデルを有限要素の各積分点で動かす」 | 連成の図（日本語） | `keio_report/fig_coupling_schematic_agreed.png` |
| `deck_oliver_en.js`「Appendix: notation and equations」 | 記号と式を 1 枚に | なし |

### 2.2 点モデルだけ

| スライド | 何を見せるか | 図 |
|---|---|---|
| `deck_en.js` 5「2. … the point model alone, Fritsch et al. 2025, case I」 | 場なしの点モデル：抗菌薬なし，$T^*=1/3$ から，最初から，の 3 通り | `assets/fig_fritsch_case1.png` |
| `slides_defense.tex`「Part 1: Bayesian calibration of the point model」「Part 1: what the data determine」 | 5 菌種の推定と，データが決める成分 | `thesis/figures/eml_fit.pdf`，`eml_ident.pdf`，`eml_checks.pdf` |
| `slides_defense.tex`「Backup 1」〜「Backup 4」 | 停留条件，対称な相互作用，4 段と 3 つの種，事後分布 | `thesis/figures/eml_seed_maps.pdf`，`eml_posteriors.pdf` |
| `deck_kadai_ja.js` 3「2. 相互作用の推定：TMCMC と GPU 計算」，4「3. 推定結果：データが決める相互作用」 | 推定の方法と結果（日本語） | `build/eml_ident.png` |
| `deck_k26c4_en.js`「6. Against the point model alone」（**温存**） | 四種：連成と点モデルだけの比較 | `assets/fig_four_species_living_vs_0d.png` |

### 2.3 2 つの連成の例（分けた配置と混ぜた配置，抗菌薬あり）

| スライド | 何を見せるか | 図 |
|---|---|---|
| `deck_en.js` 2「1. Objective and set-up」／`deck_ja.js`「1. 目的と問題設定」 | 目的と設定 | なし |
| `deck_en.js` 7 と 10「…: initial and boundary conditions」／`deck_ja.js`「4. 初期条件と境界条件」 | 初期条件と境界条件の図 | `assets/fig_setup_separate_ab.png`，`assets/fig_setup_mixed_ab.png` |
| `deck_en.js` 8 と 11「…: results with antibiotic」／`deck_ja.js`「5. 結果：組成・場・応力の分布」 | 3 次元の 5 つのパネル（(d) は $\alpha-1$） | `assets/fig_coupled_3d_panels_{separate,mixed}_ab_alpha.pdf`（`build.sh` が PNG にする） |
| `deck_en.js` 9 と 12「…: time evolution (T* = 0 to 1)」／`deck_ja.js`「6. 時間発展」 | 動画（5 つのパネル） | `assets/anim_coupled_tvd_{separate,mixed}_ab_inhliv_transc_alpha_5panel.mp4` と `_end.png` |
| `deck_en.js` 13「10. Seed means over time」 | シードの平均の時間変化 | `assets/fig_coupled_time_history.png` |
| `deck_en.js` 14「11. Seed means at T* = 1」／`deck_ja.js`「定量評価：シード領域の平均」 | 数値の表 | 表だけ |
| `deck_en.js` 17「13. Discussion」／`deck_ja.js`「考察」 | なぜ応力がほとんど変わらないか | なし |
| `deck_ja.js`（`mixed`）「前回（11/2）の補足 1：分布は x，z について対称」 | 切り取らない図と断面で対称性を見せる | `assets/fig_uncut_views.png`，`assets/fig_symmetry_sections.png` |
| `deck_ja.js`（`mixed`）「前回（11/2）の補足 2：応力は生存量ではなく成長 α で決まる」 | 抗菌薬の後も応力が増え続けること | `assets/fig_stress_living_over_time.png` |
| `slides_defense.tex`「Part 3: two species with nutrient and antibiotic」 | コロキウム用の 1 枚（(d) は抗菌薬） | `assets/fig_coupled_3d_panels_separate_ab.pdf` |
| `deck_kadai_ja.js` 6「5. 抗菌薬を含む 2 菌種の例」 | 同じ例（日本語） | `keio_report/fig_coupled_3d_panels_separate_ab.png` |
| `deck_en.js`「Appendix C. A nutrient supply that oscillates in time」（**温存**，`WITH_SINE=1` のときだけ） | 振動する栄養：点モデルだけと連成，表つき | `assets/fig_nutrient_sine.png` |

### 2.4 格子収束と制限関数

| スライド | 何を見せるか | 図 |
|---|---|---|
| `deck_en.js` 15「12. Grid convergence of growth and stress」 | 64〜192³ の表と GCI，superbee で 4.6 % | 表だけ |
| `slides_defense.tex`「Part 3: the front term sets the accuracy of the stress」 | 前線の項 $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$ と応力の収束 | `papers/coupling_paper/fig_grid_examples.pdf` |
| `slides_defense.tex`「Backup 17: slope limiters for the front term」 | 制限関数の教科書的な図（Sweby 図と 1 次元の移流） | `assets/fig_limiters_general.png` |
| `slides_defense.tex`「Part 2: convergence and stability」「Appendix D: stability of the explicit field update」 | 検証問題の収束と，陽解法の安定 | `assets/fig_wp3_convergence.png`，`assets/fig_wp3_stability.png` |
| `slides_progress_1024.tex`「The front term … is an advection」「Grid convergence: the stress needs a second-order front term …」 | 前線の項が移流であること，2 次が要ること | なし |
| `slides_1102_muramatsu.tex`「検証：格子収束と前線の項 …」 | 同じ内容（日本語） | なし |
| `deck_kadai_ja.js` 7「6. 数値の信頼性と形状」（非対称の図つき：オリバーには**温存**） | 制限関数の表，2 段の格子，柱，配置の影響 | `thesis/figures/fig_asymmetric_starts.png` |

### 2.5 要素（ANSYS）と Python の比較

| スライド | 何を見せるか | 図 |
|---|---|---|
| `deck_en.js` 16「12b. Your element and the Python code」／`deck_oliver_en.js`「2. Your element and the Python code」 | 同じ例での表：要素は 8〜24³，Python は 64〜192³，外挿値と GCI | 表だけ |
| `deck_oliver_en.js`「3. I checked the Python code」「4. Three questions about your element」 | Python の点検と，要素についての 3 つの質問 | なし |
| `slides_progress_1024.tex`「What the comparison found in the element and in Python」 | 比較で見つかったこと | なし |
| `slides_defense.tex`「Part 2: verification against an exact solution」 | 厳密解との比較 | `assets/fig1005_exact.png` |
| `slides_defense.tex`「Part 2: growth-induced stresses with the published stiffness」 | 論文の剛性での応力 | `assets/fig_element_beta002.png` |
| `slides_defense.tex`「Backup 8」「Backup 10」「Backup 12」 | 式 36 と要素の成長変数，表面の点の修正，要素の解法の流れ | `assets/fig1005_whole_model.png`，`assets/flow_oliver_solution_loop.png` |
| `slides_coupling_1020.tex`「First ANSYS run: what is checked」 | 点モデルの追跡の確認 | `assets/fig_living_trace_check.png` |

### 2.6 抗菌薬と，生きている菌からの成長

| スライド | 何を見せるか | 図 |
|---|---|---|
| `deck_living_ja.js` 2「1. 問題：今の式では死んだバイオフィルムも成長する」 | 式 36 の限界 | なし |
| `deck_living_ja.js` 3「2. 式：生きている量で成長を決める」 | $\dot\alpha=k_\alpha\phi_{3D}\sum_i\bar\phi_i$ | なし |
| `deck_living_ja.js` 4「3. 結果：抗菌薬が成長と応力に届く（128³）」 | 表（$\alpha-1$ と von Mises の変化） | 表だけ |
| `deck_living_ja.js` 5「4. 時間変化（64³）」 | 式 36 との比較，時間変化 | `assets/fig_growth_living_over_time.png` |
| `slides_defense.tex`「Backup 5: sensitivity of the coupled examples」 | 抗菌薬の入れ方などの感度 | なし |
| `slides_1020.tex`「Added since 5 October: growth from the living biomass」「Growth from the living biomass with synthetic data」 | 10 月の最初の結果（古い世代） | `assets/fig_coupling_schematic_live.png`，`assets/fig_scheme_c_synthetic.png` |

生きている菌からの成長は修論の本文に入れない（Outlook の 1 段落だけ）．慶應の続きの材料．

### 2.7 形状：柱，迷路，歯

| スライド | 何を見せるか | 図 |
|---|---|---|
| `deck_geom_en.js`「2. Rigid obstacles: set-up and growth」 | 剛な 4 本の柱（Klempt 2024 の 4.3 節） | `assets/fig_obstacle_columns_3d.pdf` |
| `deck_geom_en.js`「2. Rigid obstacles (cont.): initial and boundary conditions」「4. Maze: initial and boundary conditions」 | 初期条件と境界条件の表 | 表だけ |
| `deck_geom_en.js`「3. Rigid obstacles: time evolution and stress」 | 動画と数値 | `assets/anim_coupled_tvd_separate_ab_inhliv_obst_transc_alpha_5panel.mp4` |
| `deck_geom_en.js`「5. Maze: growth through the corridors」 | 迷路（4.4 節） | `assets/fig_maze_k24_h4.pdf`（`build.sh` が PNG にする） |
| `deck_geom_en.js`「6. Maze with antibiotic: the rule for newly reached cells decides」 | 新しく届いたセルの規則の比較（動画） | `assets/anim_maze_rules.mp4`，`assets/anim_maze_rules_end.png` |
| `slides_defense.tex`「Part 3: rigid bars and a maze」 | コロキウム用の 1 枚 | `assets/fig_obstacle_columns_3d.pdf`，`papers/coupling_paper/fig_maze_rules.pdf` |
| `slides_defense.tex`「Backup 15: outlook, a biofilm on a real molar」 | 実際の臼歯（10 月の例，連成なし） | `assets/anim_real_tooth_resolved_base_3d.png` |
| `slides_3d_demo.tex`「3. A biofilm on a real molar」〜「5. Nutrient and the growth limit on the molar」 | 臼歯の上のバイオフィルム（10 月の例） | `assets/fig_real_tooth_submodel.pdf`，`assets/fig_real_tooth_resolved.pdf` |
| `slides_1102_muramatsu.tex`「慶應での研究計画：本命の論文（2027年9月投稿）」 | 歯頸線の計画 | `assets/fig_cervical.png` |

歯頸部の窓（**温存**）を使ったスライドはまだない．図は `assets/fig_tooth_window_*.png`，ノートは
`ansys_usermat/apdl/results/2026-10-3d/TOOTH_WINDOW.md`．

### 2.8 力学：式 30 の項（全部 **温存**）

| スライド | 何を見せるか | 図 |
|---|---|---|
| `deck_deform_ja.js` 2「1. 一方向の連成と双方向の連成」 | 連成の向きの図 | `meetings/keio_deformation/eq/flow.tex` |
| `deck_deform_ja.js` 3「2. 式 30 の力学の項」 | $-K\phi_{3D}(\mathbf I:\mathbf C_e-3)$ | `meetings/keio_deformation/eq/eq30.tex` |
| `deck_deform_ja.js` 5「4. 結果：自由に育つ例（32³，K = 10⁵）」 | 断面での変化 | `assets/fig_mech_feedback_section_K100000.png` |
| `deck_deform_ja.js` 6「4b. 結果：K の掃引（64³〜192³，自由に育つ例）」 | K の掃引 | `assets/fig_mech_K_sweep.png` |
| `deck_deform_ja.js` 7「4c. 確認：時間刻みと発散の原因（32³，分けた例）」 | 発散が数値によること | 表だけ |
| `deck_deform_ja.js` 8「5. 結果：剛な柱で閉じ込めた例（K = 10⁵）」 | 閉じ込めたときの正のフィードバック | `assets/fig_mech_feedback_section_bars_K100000_clip.png` |
| `deck_deform_ja.js` 9「6. 先行研究との関係」 | 位置づけ | なし |

この項の形は Felix と話してから先に進める（CLAUDE.md）．

### 2.9 四種（全部 **温存**，2027 年 6 月ごろ）

| スライド | 何を見せるか | 図 |
|---|---|---|
| `deck_k26c4_en.js`「2. From two to four species: what changes」 | 2 種から 4 種で変わること | 表だけ |
| `deck_k26c4_en.js`「3. …」「4. Time evolution, …」 | 動画（分けた配置，混ぜた配置） | `assets/anim_coupled_tvd_{separate,mixed}_ab_inhliv_k26c4_transc_alpha_5panel.mp4` |
| `deck_k26c4_en.js`「5. Share of each species」 | 種ごとの割合 | `assets/fig_four_species_shares_{separate,mixed}_ab.png` |
| `deck_k26c4_en.js`「7. Seed means over time」 | 時間変化 | `assets/fig_four_species_time_history.png` |

### 2.10 コロキウムの Q&A 用の予備スライド（`slides_defense.tex`）

| スライド | 何を見せるか | 図 |
|---|---|---|
| Backup 1〜4 | 較正：停留条件，対称性，4 段と 3 つの種，事後分布 | `thesis/figures/eml_seed_maps.pdf`，`eml_posteriors.pdf` |
| Backup 5，6 | 連成の例の感度，$\beta$ の感度 | なし |
| Backup 7 | Klempt 2024 の 2 つの成長則 | `assets/fig_klempt2024_alpha_law.png` |
| Backup 9 | メッシュと非圧縮性 | なし |
| Backup 11 | Klempt 2024 の 4.1 節の再現 | `assets/fig_klempt2024_fig3_blend.png` |
| Backup 13 | 検証問題の 3 次元の図 | `assets/anim_wp3_seed_3d.png` |
| Backup 14 | 壁での応力と恒常圧の成長則（今後） | `assets/anim_law_3d.png` |
| Backup 16 | 較正した 5 菌種を要素で動かした例 | `assets/anim_cal5_ch.png` |
| Backup 17，18 | 制限関数，1 ステップで何をどこで解くか | `assets/fig_limiters_general.png` |

質問と答えの文は `slides_qa_ja.tex`（日本語）と `QUIZ_COLLOQUIUM.ja.md`．

---

## 3. 図（`assets/`）

書き方：図｜作るスクリプト｜使っている場所．拡張子が 2 つあるものは `pdf/png`．
スクリプトの場所を書いていないものは `ansys_usermat/apdl/`．
「修論」は `thesis/chapters/*.tex`（`main.tex` が読む章），「論文」は `papers/coupling_paper/manuscript/coupling_paper.tex`，
「報告」は `keio_report/kadaikenkyu_nishioka.tex`（自分の写しを `keio_report/` に持つ）．

多くのスクリプトは `results/` や `figdata/` の計算結果を読む．結果がないと作り直せない．

### 3.1 連成の全体図と設定

| 図 | スクリプト | 使用先 |
|---|---|---|
| `fig_coupling_schematic_agreed.png` | `ansys_usermat/coupling_schematic.py --agreed` | `slides_defense`，`slides_progress_1024`，`slides_1102_nem_ja`，`slides_1102_muramatsu`，`deck_kadai_ja`，論文，報告 |
| `fig_coupling_schematic_agreed_thesis.png` | `ansys_usermat/coupling_schematic.py --agreed --thesis` | 修論 1 章 |
| `fig_coupling_schematic.png`，`fig_coupling_schematic_live.png` | `ansys_usermat/coupling_schematic.py` | `slides_1005`，`slides_1020`（古い連成） |
| `fig_coupling_impl_flow_slide_ja.pdf/png` | `ansys_usermat/coupling_impl_flow_slide_ja.py` | `slides_1102_muramatsu` |
| `fig_setup_separate_ab.pdf/png` | `fig_setup_separate_ab.py` | `deck_en`，`deck_ja` |
| `fig_setup_mixed_ab.pdf/png` | `fig_setup_separate_ab.py mixed` | `deck_en`，`deck_ja`（mixed） |
| `fig_model_setup.png` | `ansys_usermat/model_setup_fig.py` | `slides_1005`，`slides_1020`，`slides_defense` |
| `fig_implementation_tools.pdf/png` | `fig_implementation_tools.py` | `slides_1020` |

### 3.2 点モデルと栄養

| 図 | スクリプト | 使用先 |
|---|---|---|
| `fig_fritsch_case1.pdf/png` | `fig_fritsch_case1.py` | `deck_en`（25 Nov） |
| `fig_nutrient_sine.pdf/png` | `fig_nutrient_sine.py` | `deck_en` 付録 C．**温存** |

### 3.3 2 つの連成の例

| 図 | スクリプト | 使用先 |
|---|---|---|
| `fig_coupled_3d_panels_separate_ab.pdf/png` | `fig_coupled_3d_panels.py separate_ab` | 修論 4 章，論文，報告，`slides_defense`，`slides_progress_1024`，`slides_1102_muramatsu`，`deck_kadai_ja` |
| `fig_coupled_3d_panels_mixed_ab.pdf/png` | `fig_coupled_3d_panels.py mixed_ab` | 修論 4 章，論文 |
| `fig_coupled_3d_panels_{separate,mixed}_ab_alpha.pdf/png` | `fig_coupled_3d_panels.py`（(d) を $\alpha-1$ にした版；引数は冒頭の説明に書かれていない） | `deck_en`，`deck_ja`（`build.sh` 経由） |
| `anim_coupled_tvd_{separate,mixed}_ab_inhliv_transc_alpha_5panel.mp4`，`_end.png` | `anim_coupled_tvd_3d.py <case> inhliv transc` | `deck_en`，`deck_ja` |
| `fig_coupled_time_history.pdf/png` | `fig_coupled_time_history.py` | `deck_en` |
| `fig_stress_living_over_time.png` | `fig_stress_living_over_time.py` | `deck_ja`（mixed） |
| `fig_symmetry_sections.pdf/png` | `fig_symmetry_sections.py` | `deck_ja`（mixed） |
| `fig_uncut_views.png` | `fig_uncut_views.py` | `deck_ja`（mixed） |
| `anim_coupled_tvd_mixed_ab_inherit.mp4`，`_inherit_end.png` | `anim_coupled_tvd_3d.py` の古い版（今の出力名は `_5panel` つき） | `slides_progress_1024` |

### 3.4 格子収束，制限関数，検証

| 図 | スクリプト | 使用先 |
|---|---|---|
| `fig_limiters_general.pdf/png` | `fig_limiters_general.py` | `slides_defense` Backup 17 |
| `fig_limiter_grid.pdf/png` | `fig_reference_style.py` | 修論 4 章（`thesis/figures/` の写し），論文 |
| `fig_hydrostatic.pdf/png`，`fig_spacetime_y.pdf/png` | `fig_reference_style.py` | 論文 |
| `fig_edge_vs_stress.pdf/png` | `fig_edge_vs_stress.py` | 論文 |
| `fig_wp3_convergence.pdf/png`，`fig_wp3_stability.pdf/png` | `wp3_figs.py` | 修論 4 章，`slides_defense` |
| `fig_nem_stability.pdf/png` | `nem_stability.py`（修論の付録のコメントによる） | 修論の付録 |
| `fig1005_exact.pdf/png`，`fig1005_whole_model.pdf/png` | `figs_1005.py` | 修論 4 章，`slides_1005`，`slides_1020`，`slides_defense` |
| `fig1005_element220.png` | `figs_1005.py` | `slides_1005`，`slides_1020` |
| `fig_element_beta002.pdf/png` | `fig_element_beta002.py` | 修論 4 章，`slides_defense` |
| `fig_mesh_convergence.pdf/png` | `fig_mesh_convergence.py` | `slides_1020` |
| `fig_living_trace_check.pdf/png` | `fig_living_trace_check.py` | `slides_coupling_1020` |
| `v222_closed_form.pdf/png` | `plot_closed_form.py` | 修論 4 章，`slides_defense` |
| `v222_convergence.pdf/png` | `plot_convergence.py` | 修論 4 章 |
| `growth_cylinder_alpha_sweep.pdf/png` | `plot_cylinder_alpha_sweep.py` | 修論 4 章 |
| `fig_klempt2024_alpha_law.pdf/png` | `JAXFEM/klempt2024_alpha_law.py` | 修論 4 章，`slides_defense` |
| `fig_klempt2024_fig3.png` | `JAXFEM/klempt2024_fig3.py` | `slides_1005`，`slides_1020` |
| `fig_klempt2024_fig3_blend.png` | `JAXFEM/klempt2024_fig3.py`（出力名を引数で渡す；推定） | 修論の付録，`slides_defense` |
| `fig_klempt2024_table4_blend.pdf/png` | `JAXFEM/klempt2024_table4_fig.py`（同上；推定） | 修論の付録 |

`papers/coupling_paper/` にあって `assets/` にない図：`fig_grid_examples.pdf`（`papers/coupling_paper/fig_grid_examples.py`；`slides_defense`，論文），
`fig_convergence.pdf`，`fig_validation_volume.pdf`（論文），`fig_maze_rules.pdf`（`papers/coupling_paper/fig_maze_rules.py`；`slides_defense`）．

### 3.5 形状：柱，迷路，歯

| 図 | スクリプト | 使用先 |
|---|---|---|
| `fig_obstacle_columns_3d.pdf/png` | `fig_obstacle_columns_3d.py` | 修論 4 章，論文，`slides_defense`，`slides_progress_1024`，`deck_geom_en` |
| `anim_coupled_tvd_separate_ab_inhliv_obst_transc_alpha_5panel.mp4`，`_end.png` | `anim_coupled_tvd_3d.py separate_ab inhliv obst transc` | `deck_geom_en` |
| `fig_maze_k24_h4.pdf` | `papers/coupling_paper/fig_maze_k24.py --h 4`（出力は `papers/coupling_paper/`；`assets/` のものは写し） | 修論 4 章，論文，`slides_progress_1024`，`deck_geom_en` |
| `anim_maze_rules.mp4`，`anim_maze_rules_end.png` | `papers/coupling_paper/anim_maze_rules.py` | `deck_geom_en` |
| `anim_maze_rules_end_crop.png` | 不明（`_end.png` を切り取ったもの） | `slides_progress_1024` |
| `fig_tooth_window_setup.pdf/png` | `fig_tooth_window_setup.py` | 論文（`papers/coupling_paper/` の写し）．**温存** |
| `fig_tooth_window_section_*.png`（10 枚） | `fig_tooth_window_section.py --nf 64 --window mixed …` | 論文は `papers/coupling_paper/fig_tooth_window_section.png` を使う．スライドでは未使用．**温存** |
| `fig_tooth_window_thickness.pdf/png` | `fig_tooth_window_thickness.py` | 未使用．**温存** |
| `anim_real_tooth_resolved_base_3d.mp4/png` | `anim_real_tooth_resolved.py` | `slides_defense`，`slides_3d_demo` |
| `anim_real_tooth_resolved_gate_3d`，`anim_real_tooth_resolved_nutrient_gate_3d`（`mp4/png`） | `anim_real_tooth_resolved.py` | `slides_3d_demo` |
| `fig_real_tooth_resolved.pdf/png` | `fig_real_tooth_resolved.py` | `slides_3d_demo` |
| `fig_real_tooth_submodel.pdf/png` | `real_tooth_submodel.py` | `slides_3d_demo` |
| `fig_wall_convergence.pdf/png` | `wall_stress_convergence.py` | `slides_3d_demo` |
| `anim_wp3_seed_3d.mp4/png` | `anim_wp3_seed_3d.py` | `slides_defense`，`slides_3d_demo` |
| `anim_seed_corner_3d`，`anim_seed_edge_3d`（`mp4/png`） | `anim_seed_offcentre_3d.py corner|edge` | `slides_3d_demo` |
| `anim_law_3d`，`anim_two_3d`（`mp4/png`） | `anim_scenarios_3d.py`（`law`，`two`） | `slides_defense`（law），`slides_3d_demo` |
| `anim_cervical_f4_3d.mp4/png` | `p3/uel/anim_cervical_3d.py` | `slides_3d_demo` |
| `fig_cervical_f4.pdf/png` | `p3/uel/fig_cervical_f4.py` | `slides_3d_demo` |
| `fig_cervical.png` | `p3/uel/fig_cervical.py` | `slides_1102_muramatsu` |
| `fig_tooth_homeostatic.png` | `keio_wp2/fig_tooth.py` | `slides_1102_nem_ja` |
| `fig_nutrient_homeostatic.png` | `keio_wp2/fig_nutrient.py` | `slides_1102_nem_ja` |

非対称の初期配置の図は `assets/` ではなく `thesis/figures/fig_asymmetric_starts.pdf/png`（`fig_asymmetric_starts.py`，出力先は引数で渡す）：
修論 4 章，報告，`deck_kadai_ja`．オリバーには**温存**．

### 3.6 力学の項，生きている菌からの成長，四種

| 図 | スクリプト | 使用先 |
|---|---|---|
| `fig_mech_K_sweep.pdf/png` | `fig_mech_K_sweep.py` | `deck_deform_ja`．**温存** |
| `fig_mech_feedback_section_K100000.pdf/png` | `fig_mech_feedback_section.py --nf 32 --K 100000` | `deck_deform_ja`．**温存** |
| `fig_mech_feedback_section_bars_K100000_clip.pdf/png` | `fig_mech_feedback_section.py`（柱と切り詰めの引数つき；名前は冒頭の説明に書かれていない） | `deck_deform_ja`．**温存** |
| `fig_growth_living_over_time.pdf/png` | `fig_growth_living_over_time.py` | `deck_living_ja` |
| `fig_scheme_c_synthetic.pdf/png` | `scheme_c_synthetic.py` | `slides_1020`（古い世代） |
| `fig_four_species_shares_{separate,mixed}_ab.png` | `fig_four_species_shares.py separate_ab|mixed_ab` | `deck_k26c4_en`．**温存** |
| `fig_four_species_living_vs_0d.png` | `fig_four_species_living_vs_0d.py` | `deck_k26c4_en`．**温存** |
| `fig_four_species_time_history.png` | `fig_four_species_time_history.py` | `deck_k26c4_en`．**温存** |
| `anim_coupled_tvd_{separate,mixed}_ab_inhliv_k26c4_transc_alpha_5panel.mp4`，`_end.png` | `anim_coupled_tvd_3d.py <case> inhliv k26c4 transc` | `deck_k26c4_en`．**温存** |

### 3.7 較正した 5 菌種と，古いスライドの図

| 図 | スクリプト | 使用先 |
|---|---|---|
| `fig_thesis_cal5.pdf/png` | `fig_thesis_cal5.py` | 修論 4 章，報告 |
| `fig_cal5_appendix.pdf/png` | `fig_cal5_appendix.py` | `slides_1020` |
| `anim_cal5_ch`（`gif/mp4/pdf/png`） | `anim_cal5_ch.py` | `slides_defense` Backup 16 |
| `fig_composition_clock.pdf/png` | `ansys_usermat/composition_clock_and_map.py` | 修論 4 章 |
| `fig_two_way_step2_thesis.pdf/png` | `ansys_usermat/two_way_step2.py` | 修論 4 章 |
| `fig_composition_transport_slide.png` | `ansys_usermat/composition_transport_check.py` | `slides_1005`，`slides_1020` |
| `fig_fields_time.pdf/png` | `fig_fields_time.py` | `slides_1020` |
| `fig_nem_demo.png` | `ansys_usermat/nem_demo.py` | `slides_1102_nem_ja` |
| `flow_oliver_solution_loop.png` | 不明（`ch5_flow/` の TikZ の図から；修論は `ch5_flow/flow_oliver_solution_loop_print.pdf` を使う） | `slides_defense` Backup 12 |
| `condition_total_phi.png` | `plot_condition_total_phi.py` | `slides_1001_oliver`，`slides_muramatsu_ja` |
| `heine_species_composition.png` | `plot_heine_composition.py`（リポジトリの一番上） | `slides_kiso_ja` |
| `fig_vol_volume.pdf`，`fig_vol_live.pdf`，`fig_vol_species.pdf`，`fig_vol_fit.pdf`，`fig_vol_stress.pdf` | `volume_figs.py` | `slides_volume_ja` |

修論の 3 章の図（`eml_*.pdf`，`dieckow_loo_*.pdf`）は `thesis/figures/` にあり，`assets/` にはない．

### 3.8 使われていない図

`assets/` には，どのスライドにも修論にも論文にも出てこない図が 150 個ほどある．主なもの：

- 新しいもの：`fig_tooth_window_thickness`，`fig_tooth_window_section_*`（**温存**），
  `anim_coupled_tvd_separate_ab_inhliv_b1edge_seedx-1_transc_alpha_5panel` と `…_bcorner_…`（非対称の配置の動画，**温存**），
  `anim_coupled_tvd_mixed_ab_inhliv_obst_transc_alpha_5panel`（混ぜた配置と柱），
  `anim_coupled_tvd_*_inhliv_transc_alpha_hist_5panel`（時間変化のグラフつきの動画），
  `fig_profiles_y`，`fig_cost_breakdown`（`fig_reference_style.py`），`fig_coupled_3d_panels_separate_ab_startpos`，
  `fig_coupled_examples`．
- 連成の説明用：`fig_coupling_equations`，`fig_coupling_equations_ja`，`fig_coupling_glossary_en`，`fig_coupling_glossary_ja`，
  `fig_coupling_impl_flow_ja`，`fig_coupling_impl_flow_slide_en`，`fig_coupling_loop_simple_ja`，`fig_coupling_schematic_keio`，
  `fig_coupling_schematic_thesis`，`fig_point_model_klempt2026`．
- 古い世代（使う前に中身を確かめる）：`fig_meeting1019_*`，`fig_composition_*`，`fig_fuel_front_*`，`fig_obstacle_columns`，
  `fig_obstacle_columns_nf32`，`fig_maze`，`anim_real_tooth_*`（`resolved` でないもの），`anim_nutrient_*`，`anim_u1_16_*`，
  `v222_conditions*`，`flow_*`，`umat_flow_*`，`algo_flow_*`，`fig_surrogate_*`，`fig_abaqus_klempt2024_table3*`（Abaqus：修論の内容ではない）．

### 3.9 ほかに試した設定（25 Nov に図にしたもの，全部**温存**）

計算は前からあったもの（4 菌種の栄養の 3 本だけ 25 Nov の新しい計算）．どれもスライドにはまだ入れていない．

| 図 | スクリプト | 中身 | 状態 |
|---|---|---|---|
| `fig_antibiotic_schedules.pdf/png` | `fig_antibiotic_schedules.py` | 抗菌薬の与え方：一定，パルス（周期 2/9 の前半に 2），正弦 $1+\sin(2\pi s)$，$T^*=1/3$ から．64³，分けた配置．(a) 与える値，(b) 種子の $\Sigma_i\bar\phi_i$，(c) 種子の抗菌薬．$\alpha-1$ と vM の終わりの値は標準出力と `papers/coupling_paper/SINE_SCHEDULES.ja.md` | **温存**（生きている菌からの成長の版の数字も出すので慶應の続き） |
| `fig_stiff_species.pdf/png` | `fig_stiff_species.py` | 菌種で硬さを変える $E=E(\phi_{3D})\sum_i w_i m_i$（`--stiff-species`）：$m_A,m_B$ = (1,1), (2,1), (1,2), (0.5,1), (1,0.5)，分けた配置と混ぜた配置の種子 vM．$\alpha-1$ は変わらない．表は `papers/coupling_paper/EXPERIMENT_COMPARISON.ja.md` §4 | **温存**（式を変えるので慶應の続き．種ごとの硬さのデータはない） |
| `fig_nutrient_sine_4sp.pdf/png` | `fig_nutrient_sine_4sp.py` | `fig_nutrient_sine` と同じ並びを 4 菌種（Klempt 2026 case 3，A の 1/2 なし）で：上は点モデルだけ（$c^*$ = 100，$50+50\sin(500t)$，50），下は連成 64³ の種子 | **温存**（四種の例と振動する栄養の両方） |

---

## 4. 再利用するときの決まり（CLAUDE.md から）

- **一人称は I**．we，our は使わない．日本語は「私」．
- **言い回し**：強い決まり文句（strictly，exactly as published など）とダッシュの多用を避ける．太字は大事な数字だけ．
- **図の文字は Times**：図のスクリプトは `ansys_usermat/figstyle.py` の `figstyle.apply()` を呼ぶ．図の中に題は入れない．パネルの (a)，(b) は説明文と合わせる．
- **記号は論文どおり**：
  - 場の式では $\phi_{3D}$，点モデルでは $\phi_i$，$\psi_i$，$\bar\phi_i=\phi_i\psi_i$．`\phi` を使う（`\varphi` は使わない）．
  - 生きている体積分率は $\sum_i\bar\phi_i$．「L」とは書かない．
  - 成長は $\alpha-1$，$\mathbf F_g=\alpha\mathbf I$，$\alpha(0)=1$．成長の速さは $k_\alpha$．
  - 種の割合は $\phi_1/(\phi_1+\phi_2)$ と書く．新しい記号を作らない．
- **項は名前だけで書かない**：式 34 の前線の項は必ず $-r\,\frac{c}{k+c}\,\mathbf n_c\cdot\nabla\phi_{3D}$（$\mathbf n_c=\nabla c/|\nabla c|$）と式で書く．
- **式 35 は $\dot c$ つき**で書く．消費は $g\,\phi_{3D}\,(\sum_i\bar\phi_i)\,c$．
- **点モデルは「論文どおりに使う」**（used as published）．「black box」とは書かない．
- **内部の名前を出さない**：段の番号，`prop(28)` のモード，計算の名前（`u8_` など）はスライドに書かない．何をした計算かを書く．
- **論文にない値は印をつける**（$f$，$\phi_{\min}$，配置，換算した $\beta$ など）．仮定は仮定と書き，感度の計算を添える．
- **単位を書く**．Klempt 2024 の値は $E=10$ Pa，$\nu=0.49$．
- **Abaqus は修論の内容ではない**．「慶應での続き」としてだけ触れる．
- **日本語では「バイオフィルム」**と書く（「膜」だけにしない）．オリバーは「オリバー」．Meisam は Dr. Soleimani（教授ではない）．
- **付録に記号と式**を入れる（手本は `slides_1005.tex` の付録 A〜C）．
- **3 次元の図と動画**：向きは 1 つ（y が上，カメラ (3.6, 2.6, 4.4)，x, y, z > 0 の 1/8 を切り取る），座標軸を必ず入れる．
  動画の 5 つのパネルは (a) 菌種 A の割合，(b) 生きている体積分率，(c) 栄養 $c$，(d) 成長 $\alpha-1$，(e) von Mises 応力．
  抗菌薬を入れる時刻をフレームに書く．
- **時間変化のグラフを添える**：3 次元の図や動画の隣に，シードの平均の時間変化（`fig_coupled_time_history.py`，
  `fig_stress_living_over_time.py`，`fig_four_species_shares.py`）．
- **指導教員向けは 20 ページ以内**．LaTeX のエラーなし，はみ出しなし．
- **送る前に作って，出来上がったページを目で見る**（`build_slides.ps1`，または `pdflatex`／`lualatex`，pptx は開いて確かめる）．
  数字は修論 4 章の表と同じか確かめる（`scripts/audit_thesis_numbers.py`）．メールは送る前にユーザーに確かめる．
- **温存のものは送らない**：歯頸部の窓，式 30 の力学の項と K の掃引，非対称の初期配置，四種の例，振動する栄養．
  `deck_en.js` は `WITH_SINE` をつけずに作る．
- **Felix の実装と理論は非公開**：スライドに中身を書かない．書いてよいのは，受け取ったこと，照合していること，公開された論文の内容だけ．
- **結果の世代を混ぜない**：10 月の古いスライド（`slides_1001_oliver.tex`，`slides_1005.tex`，`slides_1020.tex`，`slides_muramatsu_ja.tex`）の
  数字は今の結果と違う．構成や付録の手本として使い，数字は今のスライドから取る．
