# クラウド → IKMHIWI03（2026-10-19、2 本目）：円柱の図のデータをクラウドに渡す

## 背景
修論の図 4.4（`assets/growth_cylinder_alpha_sweep.png`、2 層の円筒に α − 1 = 0.01〜0.2 の成長）は、
印刷すると文字が約 6.5pt で小さい。作図はクラウドで本文の幅（16 cm）に描き直す。
元の結果（ANSYS の出力）は IKMHIWI03 の `F:\biofilm_upf_kusepy\` にしかないので、
**縮めたデータ（α − 1 ごとの、弧に沿った u_r）だけを JSON にして push してほしい**。数 KB。

## やること（数分、ANSYS は回さない）
```powershell
git pull origin master
python ansys_usermat\apdl\plot_cylinder_alpha_sweep.py `
    0.01=F:/biofilm_upf_kusepy/gcr_g0.txt 0.02=F:/biofilm_upf_kusepy/asw_0p02.txt ... `
    -o assets\growth_cylinder_alpha_sweep.png
```
- 引数は、10/8〜11 に図 4.4 を描き直したとき（コミット 52efd6c）と同じ出力ファイル。α − 1 = 0.01、0.02、0.05、0.1、0.2 の 5 本
  （今の図の右の凡例と同じ）。ファイル名が違ったら、`F:\biofilm_upf_kusepy\` の中で探して合わせる。
- スクリプトは今回、引数があると縮めたデータを `ansys_usermat/apdl/results/cylinder_alpha_sweep.json` に書き出す
  ように直してある（図も新しい形で描く）。
- **push するのは `ansys_usermat/apdl/results/cylinder_alpha_sweep.json` と `assets/growth_cylinder_alpha_sweep.png`**。
  ファイルを名前で指定してコミット（`git add -A` はしない）。

## 確認
- 出力に α − 1 ごとの「max u_r」が 5 行出ること。今の図の左のパネルの値（0.2 で最大）と同じ桁であること。
- 図の文字が読める大きさで、凡例が線に重なっていないこと。重なっていたら、そのまま push してよい（クラウドで直す）。

push されたら、クラウドで JSON から図を描き直し、キャプションと合わせて修論に入れる。
