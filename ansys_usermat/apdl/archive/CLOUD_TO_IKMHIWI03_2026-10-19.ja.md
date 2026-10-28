# クラウド → IKMHIWI03（2026-10-19）：第 4 章の 5 菌種 8³ を最終の CH MAP で 1 本

## 背景
修論 4.7.1（`thesis/chapters/ch4_ansys.tex` の「The calibrated five-species model in the element」）は
`w8_CH_g1_s015_raw_d1`（CH MAP 6e18e0b = seed 123 の最良の粒子、−26.52）を使っている。
第 3 章の最終の CH MAP は 10/16 に差し替えた seed 7 の磨いた点（−25.87）。差は CH で同定されない成分だけで、
点モデルの軌道の差は組成で最大 0.031、全体積で 0.006、生きた体積で 0.037（EML 側の計算）。
成長と応力は一方向連成なので変わらないはず。組成の数字（So 0.74 → 0.54 など）を最終の MAP にそろえるための 1 本。

## やること（8³、1〜2 分）
```powershell
git pull origin master          # Portable Git、CLAUDE.md の注意どおり
.\ansys_usermat\apdl\run_calibrated.ps1 -Conditions CH -Day1 -WorkDir F:\biofilm_upf_cal5d1 `
    -MapSha latest -Tag final_d1 -Push
```
- `-MapSha latest`：`map_pin.txt` の固定を無視して、Tmcmc202601 の `final_theta_MAP/CH.json` の最新の commit を使う。
  ログ（`F:\biofilm_upf_cal5d1\_cal5.log`）で、使った sha が 10/16 の差し替えの commit（6e18e0b ではない）であることを確認。
- ジョブ名は `w8_CH_g1_s015_final_d1`。出力 `ansys_usermat/apdl/results/2026-10-cal5/w8_CH_g1_s015_final_d1.json` と `_pm.json` が push される。
- 終わったら `.lock` が残っていないか確認（再起動のあとは CLAUDE.md の注意どおり）。

## 確認と報告
1. シードの α − 1 平均と von Mises 平均が `w8_CH_g1_s015_raw_d1` と全桁一致すること（6.752e−4、4.68e−4 Pa）。違ったら報告。
2. シードの Gauss 点での割合（So、An、Vd、Fn、Pg）の開始と T* = 1、生きた体積 Σφψ の開始と終わり。
   今の本文：So 0.74 → 0.54、Vd 0.15、An 0.13、Pg 0.08、Fn 約 0.09、Σφψ 0.81 → 0.24。差が 0.03 を超えたら報告に書く。
3. Day 1 の組成（phi_init.json）が 6e18e0b のものと同じかどうか。

push されたら、クラウドで第 4 章の数字・図（`fig_thesis_cal5.py`）・監査スクリプトを新しい JSON に切り替える。
