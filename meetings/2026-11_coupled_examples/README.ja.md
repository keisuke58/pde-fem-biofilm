# 2 つの連成の例（分離・混合、抗菌薬あり）のスライド（1 Nov 2026）

- `deck_en.js`：英語版、分離と混合を 1 つにまとめたもの（16 枚）。**オリバー、Meisam、Felix に再来週ごろ（11 月中旬）送る**（ユーザー、1 Nov）。
- `deck_ja.js`：日本語版、例ごとに 1 つ（11 枚、発表者ノートつき）。分離の方は 2 Nov に村松先生へ発表。
- 作り方：`bash meetings/2026-11_coupled_examples/build.sh [en|ja|all]` → `build/`（git に入れない）。
  式の画像は `eq/*.tex`、図は `assets/fig_setup_*_ab.png`（`ansys_usermat/apdl/fig_setup_separate_ab.py [mixed]`）と
  `assets/fig_coupled_3d_panels_*_ab.pdf`、動画は `ansys_usermat/apdl/anim_coupled_tvd_3d.py <case> inhliv transc`。

## 送る前に確かめること
1. ~~fifa の混合 256³~~ 反映済み（1 Nov、表 2 の下の文と日本語版の考察。論文の付録 B.4 も）。
2. 数字は修論 4 章の `tab:ch5-two-species` と `tab:ch5-gci` と同じか（`scripts/audit_thesis_numbers.py`）。
3. CLAUDE.md の「指導教員向けの文書」の決まり：一人称 I、AI っぽい言い回しなし、内部のラベルなし、論文にない値は印つき（f、φ_min、配置）。
4. メールは送る前にユーザーに確認する（このファイルでは送らない）。PPTX は 12 MB 前後なので、添付か Drive のリンクか決める。

## 障害物と迷路のスライド（2 Nov 追加）
- `deck_geom_en.js` → `build/geometries_en.pptx`（英語，10 枚，`bash build.sh geom`）．Klempt 2024 の 4.3 節（剛体の棒）と 4.4 節（迷路），境界条件の表つき．
- **渡す時期**：12 月の初め（ハノーファーにいる 11/28〜12/2），会って渡す（Junker 先生との面談があればそこで）．11 月中旬のメールは 2 つの例のスライドだけ．12/10 以降は避ける（コロキウムの準備）．
- 渡す前に：障害物の動画を α 入りの新しい版に差し替える（障害物ありの 64³ の計算待ち），「棒の位置と迷路は図から読み取った」と一言添える．
