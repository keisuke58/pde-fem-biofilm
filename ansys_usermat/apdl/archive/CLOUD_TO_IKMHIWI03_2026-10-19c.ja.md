# クラウド → IKMHIWI03（2026-10-19 c）：A-18 5 菌種 CH の場を 40 時刻で（コロキウムのアニメーション用）

## 背景
コロキウム用に、第 4 章の検証問題のアニメーション（`assets/anim_wp3_seed.gif`、`anim_wp3_seed.py`）を作った。
本命の結果（較正した 5 菌種 CH、最終の MAP、Day 1 の組成）でも同じものを作りたい。今ある場は
`fields_w16_CH_g1_s015_raw.json` の 4 時刻だけで、しかも古い MAP（6e18e0b）なので足りない。

## やること（16³ を 1 本、10〜25 分。後処理 数分）
```powershell
git pull origin master          # Portable Git、CLAUDE.md の注意どおり
.\ansys_usermat\apdl\run_calibrated.ps1 -Conditions CH -Day1 -Mesh 16 -WorkDir F:\biofilm_upf_cal5d1 `
    -MapSha latest -Tag final_d1 -Push
```
- 8³ の `w8_CH_g1_s015_final_d1`（10/19 に回したもの）と同じ設定で、メッシュだけ 16³。ジョブ名は `w16_CH_g1_s015_final_d1`。
- ログで、使った MAP の sha が 8³ の final_d1 と同じであることを確認。
- 基準デッキは `OUTRES,ALL,ALL` なので、`.rst` に全結果セットが残るはず（A-15 と同じ。16³ で 1 GB 前後、F: の空きを確認）。

### 後処理（解き直しなし、標準の ANSYS222 でよい）
新しいマクロ `ansys_usermat/apdl/post_fields_dense.mac`（`post_fields_time.mac` と同じ中身で、時刻を T* = 0.025, 0.05, …, 1.0 の
40 個にしたもの）。A-15 と同じやり方で：
```
/POST1
FILE,w16_CH_g1_s015_final_d1,rst
/INPUT,post_fields_dense,mac
```
→ `fields_1.csv` … `fields_40.csv`。続けて
```powershell
python ansys_usermat\apdl\fields_time.py w16_CH_g1_s015_final_d1 <CSV のあるフォルダ> 40
```
→ `ansys_usermat/apdl/results/2026-10-cal5/fields_w16_CH_g1_s015_final_d1_dense.json`（`fields_time.py` に引数 N を足した。
N を省くと今までどおり 4 時刻）。この JSON と、`run_calibrated.ps1 -Push` が出す `w16_CH_g1_s015_final_d1.json`、`_pm.json` を push。
CSV は push しない。

## 確認と報告
1. JSON の `t` が 0.025 刻みで 40 個あるか（その時刻のセットがなければ、最も近いセットの時刻が入る。重複があれば報告）。
2. T* = 1 の alpha と seqv が `w16_CH_g1_s015_final_d1.json` の all_stress と全要素で一致すること（A-15 と同じ確認）。
3. シード（`seed_BIOFILM1`）の平均 φ、α − 1、vM [Pa] を T* = 0.25、0.5、0.75、1.0 で表に（A-15 の表と同じ形）。
   古い MAP の 16³（A-15 の表）と比べて、α − 1 と vM は一方向連成なので全桁一致するはず。違ったら報告。
4. 点モデルの割合（So、An、Vd、Fn、Pg）の始めと T* = 1（8³ の final_d1：So 0.74 → 0.55、Σφψ 0.81 → 0.23）。16³ で 0.01 以上違えば報告。

push されたら、クラウドで `anim_wp3_seed.py` と同じ形のアニメーションを作る（φ、α − 1、vM に加えて、
シードの点モデルの割合の時間変化を横に並べる予定）。
