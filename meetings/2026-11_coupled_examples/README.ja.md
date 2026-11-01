# 2 つの連成の例（分離・混合、抗菌薬あり）のスライド（1 Nov 2026）

- `deck_en.js`：英語版、分離と混合を 1 つにまとめたもの（16 枚）。**オリバー、Meisam、Felix に再来週ごろ（11 月中旬）送る**（ユーザー、1 Nov）。
- `deck_ja.js`：日本語版、例ごとに 1 つ（11 枚、発表者ノートつき）。分離の方は 2 Nov に村松先生へ発表。
- 作り方：`bash meetings/2026-11_coupled_examples/build.sh [en|ja|all]` → `build/`（git に入れない）。
  式の画像は `eq/*.tex`、図は `assets/fig_setup_*_ab.png`（`ansys_usermat/apdl/fig_setup_separate_ab.py [mixed]`）と
  `assets/fig_coupled_3d_panels_*_ab.pdf`、動画は `ansys_usermat/apdl/anim_coupled_tvd_3d.py <case> inhliv transc`。

## 送る前に確かめること
1. fifa の混合 256³ の結果が来ていれば、表 2（格子収束）と考察の 1 文（「256³ で確認中」）を更新する。GCI は 4 格子で出し直す。
2. 数字は修論 4 章の `tab:ch5-two-species` と `tab:ch5-gci` と同じか（`scripts/audit_thesis_numbers.py`）。
3. CLAUDE.md の「指導教員向けの文書」の決まり：一人称 I、AI っぽい言い回しなし、内部のラベルなし、論文にない値は印つき（f、φ_min、配置）。
4. メールは送る前にユーザーに確認する（このファイルでは送らない）。PPTX は 12 MB 前後なので、添付か Drive のリンクか決める。
