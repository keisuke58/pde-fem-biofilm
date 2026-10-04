# 10/5（月）IKMHIWI03 手順書（ミーティング後、今週ANSYSはこの日だけ）

上から順に。時間切れならそこで止めて、やった所まで下の「記録」に書く。

**注意: ここで使うスクリプト（refine_deck.py、compare_mesh.py、neighbour_sign.py ほか）は
まだ master になく、作業ブランチ `claude/plan-next-hxjjve` にだけある。** いつもの作業ツリー
（改行コードの差分が480件ある）は触らず、別フォルダにブランチを取り出して使う:

```powershell
$env:Path = "C:\Users\nishioka\git\cmd;C:\Users\nishioka\git\mingw64\bin;" + $env:Path
git fetch origin claude/plan-next-hxjjve
git worktree add F:\pfb_1005 FETCH_HEAD      # 初回だけ。2回目以降は cd F:\pfb_1005; git pull origin claude/plan-next-hxjjve
cd F:\pfb_1005
. .\dev-env.ps1
$R = "<いつもの作業ツリー>\ansys_usermat\apdl\results\2026-10-01_paper_values"   # 8³ の結果CSVのある場所
```

`git fetch` が `couldn't set 'refs/remotes/origin/...'` と出ても FETCH_HEAD は正しい
（CLAUDE.md の rename-lock）。ANSYS の作業フォルダはいつもどおり `F:\biofilm_upf_wired`。

---

## 1. B-bar の確認（5分、ANSYS実行なし）

8³ の論文値の実行（ds_pv_eq36）の出力 `ds.dat` か `out.txt` を開いて、要素タイプ1の KEYOPT を見る。

```powershell
Select-String -Path F:\biofilm_upf_wired\out*.txt -Pattern "KEYOPT|SOLID185" | Select-Object -First 20
```

- **KEYOPT(2) = 0（または記載なし＝既定値0）** → B-bar。修論の「ロッキングなし」はそのままでよい。
- **KEYOPT(2) = 1 などそれ以外** → 完全積分または縮退積分。その場で記録し、クラウドで修論を直す（8³ の平均応力は B-bar 比で 3.7 倍になり得る）。

出力に出ていなければ、デッキに `ETLIST` を1行足して1回だけ流す（数秒）。

## 2. 16³ の実行（本命、計算30分〜1時間の見込み）

`<BASE>` は 8³ の論文値の実行に使った相方の stage-1 デッキ（figs_1005.py の冒頭と同じもの）。

```powershell
python ansys_usermat\apdl\refine_deck.py F:\biofilm_upf_wired\<BASE>.dat F:\biofilm_upf_wired\ds16_base.dat --n 16
#   表示: old element 220 -> new elements [2151, 2152, 2167, 2168, 2407, 2408, 2423, 2424]
python ansys_usermat\apdl\make_wired_deck.py F:\biofilm_upf_wired\ds16_base.dat F:\biofilm_upf_wired\ds16_pv_eq36.dat `
    --set K_LOCAL1=1e-3 --set K_LOCAL2=0 --set MY_BIOSTART2=0.0 `
    --set YOUNG_BIO=1e-5 --set POISSON_BIO=0.49 --set YOUNG_VOID=-1e-3 `
    --props 7=1e-3,28=1 --post both --post-elem 2151
.\ansys_usermat\apdl\run_apdl.ps1 -Deck ds16_pv_eq36.dat -WorkDir F:\biofilm_upf_wired
```

計算を流している間に 3〜5 をやる。終わったら:

```powershell
python ansys_usermat\apdl\compare_mesh.py "$R\all_stress_ds_pv_eq36.csv" `
    F:\biofilm_upf_wired\all_stress_ds16_pv_eq36.csv --grid8 8 --track 220
```

- シードが **32 要素と 256 要素**、平均 α が同じであること（違えばデッキの対応が間違い）。
- 16³/8³ の比が Python と CalculiX の予測 **0.63（von Mises）、0.67（平均応力）** に近いか。
- 表示をコピーしておく（下の記録に貼る）。

## 3. 隣接要素の応力の符号（10分、ANSYS実行なし）

```powershell
python ansys_usermat\apdl\plot_3d.py "$R\all_stress_ds_pv_eq36.csv" --grid 8 --size 2.0 --out assets\stress3d_eq36.png
python ansys_usermat\apdl\neighbour_sign.py --csv "$R\all_stress_ds_pv_eq36.csv" --grid 8
```

- `plot_3d` が **OK** と表示したときだけ `neighbour_sign` の結果を信じる（要素番号の仮定の確認）。
- Python の予測（8³）: シード −3.8e−5 Pa、第1層 +5.8e−6 Pa（56 % が引張）、第2層 +1.6e−6 Pa（79 %）。
  引張なら Klempt の「ring of tension」と一致。

## 4. Nut1 が材料ルーチンで読めるか（15分、見るだけ）

相方の材料ルーチン（`F:\biofilm_upf_wired\Usermat_P21-V21_v222.F` ほか）で、Gauss 点の栄養 `Nut1`
（`USSFin` で解いている）が材料呼び出しの中から読めるかを探す。

```powershell
Select-String -Path F:\biofilm_upf_wired\*.F -Pattern "Nut1|GetVals|USSFin" | Select-Object -First 40
```

- 読める → 次はブリッジの4つ目の入力（`c_rel`、クラウド側は実装済み）につなぐだけ。双方向連成ステップ1をANSYSで示せる。
- 読めない → 「読めない、理由」を記録。修論では「要素側の対応が必要」と書く。

読めた場合、差し込みコードは用意済み（4 Oct、テスト4件）。変数名を指定して貼り直すだけで、
`prop(33) = c_ref > 0` の点で点モデルの栄養が c* = c*_0 · min(Nut1/c_ref, 1) になる:

```powershell
python ansys_usermat\apdl\paste_fragments.py F:\biofilm_upf_wired\Usermat_P21-V21_v222.F --nut-var <変数名>
#   例: --nut-var Sdp_nut1_n  （配列なら vGdp_Nut1_n(ID) の形も可、72桁以内）
```

再ビルド（`build_wired.ps1`）と case 6 の実行は今回は時間がないので、変数名だけ記録すれば十分。
`prop(33)` を付けなければ（または `--nut-var` なしで貼れば）、これまでと同じ動作。

## 5. PR #56 の自動フッターを消す（1分）

`.env`（トークン）は**いつもの作業ツリー**にあるので、そちらで実行する（スクリプトは master にもある）。

```powershell
python scripts\strip_pr_ai_footer.py --pr 56          # 確認だけ
python scripts\strip_pr_ai_footer.py --pr 56 --apply  # 書き換え
```

## 6. 時間があれば: bio1 + bio1 の修正

`PAPER_CHECK_KLEMPT2024.md` の「Added 3 Oct」の項。1種の実行がその行を通るか確かめ、通るなら
修正・再ビルド（`link_v222.ps1`、`usermat_py_hook.f` を先にコンパイル）・論文値デッキの再実行。
半日かかり得るので、他が全部終わったときだけ。

---

## 記録（その場で書く）

| 項目 | 結果 | メモ |
|---|---|---|
| 1. KEYOPT(2) | | |
| 2. 16³ シード平均 vM / p、8³ との比 | | |
| 3. 第1層・第2層の p と引張の割合 | | |
| 4. Nut1 | 読める / 読めない | |
| 5. フッター削除 | | |
| 6. bio1+bio1 | | |

記録した表と、2・3 の画面表示をそのままクラウドのセッションに貼れば、こちらで修論4章に反映して
push する（IKMHIWI03 からブランチへ push する必要はない）。
