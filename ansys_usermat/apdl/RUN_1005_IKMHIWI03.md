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
| 1. KEYOPT(2) | 0（既定、B-bar） | デッキは `et,1,185` のみ、KEYOPT 行なし。出力に「elastoplastic なので KEYOPT(2)=2 推奨」の NOTE が出るが設定は 0 |
| 2. 16³ シード平均 vM / p、8³ との比 | vM 8.444e−5 Pa（比 1.12）、p −2.808e−5 Pa（比 0.81） | シード 32 / 256、平均 α 1.096e−3 / 1.093e−3（一致）。予測 0.63 / 0.67 とは合わない。シード外の最大 vM は 1.714e−5 → 3.109e−5（比 1.81）。要素 220 の 8 子要素 [2151…2424]：vM 3.99e−6 → 4.15e−5（10.4）、p −8.57e−5 → −7.35e−5（0.86）。167 s、エラー 0 |
| 3. 第1層・第2層の p と引張の割合 | シード −3.48e−5、第1層 +5.36e−6（56 %）、第2層 +1.53e−6（79 %）、それより外 +1.79e−7（70 %） | `plot_3d` は CHECK FAILED。原因はシードが 32 要素あり、最大 α 要素の面隣接がシード内部になること（番号の誤りではない）。下の番号の注記を参照。層の値は有効 |
| 4. Nut1 | 読める | `Usermat_P21-V21_v222.F` 403–404 行：`Sdp_nut1_n = GetVals(ofs_dp_nut1_n+(ID-1))`。USSFin が 3475 行で `SetVals(ofs_dp_nut1_n,…)`。exec の差し込み（743 行〜）はその後ろ → `--nut-var Sdp_nut1_n` |
| 5. フッター削除 | 失敗 | `--apply` が 403 "Resource not accessible by personal access token"。PAT に Pull requests の書き込み権限がない。PR は変更なし |
| 6. bio1+bio1 | 修正済み（10/1）、作業不要 | `Usermat_P21-V21_v222.F` 616 行 `Sdp_sumBio = Sdp_bio1_n + Sdp_bio2_n`（コメント「typo fixed … 2026-10-01」）、剛性の呼び出し（1220 行）で使用。ONE_SPECIES_COUPLING.md「Results 2026-10-01」のとおり 10/1 以降の exe はすべて修正入り（現行 exe 10/2 14:19）。論文値の 8³ 実行（10/1 19:53、10/2 14:11）も修正後。数値でも確認: φ² 剛性の Python が ANSYS のシード vM と 0.1 % 以内で一致（2φ なら剛性 4 倍）。COUPLING_STATUS.md の行 6 と PAPER_CHECK の「Added 3 Oct」は done にしてよい |

**16³ が予測と合わない原因（10/5、IKMHIWI03 で調査）。** 力学の解き方ではなく、成長場の違いが原因。
- `mesh_study_seed.py` と同じ B-bar ソルバに **ANSYS の要素ごとの α（SVAR 84）** を入れると、シード vM が
  ANSYS と一致する（8³ 7.505e−5 / 7.508e−5、16³ 8.449e−5 / 8.444e−5、要素ごとの相関 0.997–1.000）。
  φ（剛性 φ²）を ANSYS の値にしても変化は 2 % 以下。シード外の vM は理想場でも ANSYS と一致する
  （4.96 / 4.98e−6、4.36 / 4.73e−6）。
- 予測（0.63 / 0.67）は、シード内で α が一様（1.1e−3）という理想化した場を前提にしている。ANSYS の α は
  シード表面の 1 要素層で低い（φ がシード表面から void 側へ漏れる）: 8³ は内部 1.0999e−3 / 表面 1.0951e−3
  （幅 0.43 %）、16³ は内部 1.1002e−3 / 表面 1.0925e−3、角で 1.0815e−3（幅 1.70 %）。シードは剛性
  1e−3 の void に囲まれてほぼ自由に膨張するので、この内部と表面の α の段差 Δα による非適合ひずみ
  （E·Δα ≈ 10 Pa × 7.5e−6）がシード vM を決める。16³ では段差が大きくなるので vM は増える（比 1.12）。
  シード外の隣の要素の α も 2.3e−6 → 8.4e−6 に増える（漏れが大きい）。
- 平均応力 p は段差の影響が小さい（比 0.81。理想場のときは 0.67）。修論で引用するなら p の方が安定。
  Python に ANSYS の α を入れたときの p は ANSYS と 8–11 % 違う（Python は要素平均の α、ANSYS は積分点ごとの値）。
- 未確認: 表面の α の段差がメッシュで収束するか（φ の拡散が物理的な長さで決まるのか、NEM の近傍幅で
  決まるのか）。
- **32³ は IKMHIWI03 では解けない。** 32768 要素、積分点 262144。最初のサブステップの USSFin で
  `PARDISO phase 22 nut2 error: -2`（メモリ不足、RAM 31.5 GB、ANSYS 21 GB 使用）→
  `*** FATAL *** ERROR IKM: Not converged solution after refinments.`（Ussfin 3304 行）、848 s。
  `run_wired.ps1` は FATAL では集計行が出ないので exit 0 と報告していた → FATAL を検出して exit 3 にする修正を入れた。
- 代わりに 24³（13824 要素、シード 864、最大 7 GB）を実行中。3 点（8³ / 16³ / 24³）の表:

| mesh | seed vM | seed p | 外側 vM | α 内部 | α 表面 | α 最小 | 段差 % | 隣の void の α | Python（ANSYS の α）vM / p | Python（理想場）vM / p |
|---|---|---|---|---|---|---|---|---|---|---|
| 8³ | 7.508e−5 | −3.478e−5 | 4.98e−6 | 1.0999e−3 | 1.0951e−3 | 1.0951e−3 | 0.43 | 3.08e−6 | 7.505e−5 / −3.744e−5 | 5.703e−5 / −3.768e−5 |
| 16³ | 8.444e−5 | −2.808e−5 | 4.73e−6 | 1.0996e−3 | 1.0894e−3 | 1.0815e−3 | 0.93 | 7.90e−6 | 8.449e−5 / −2.491e−5 | 3.589e−5 / −2.519e−5 |
| 24³ | 実行中 | | | | | | | | | |

（Pa。「表面」はシード要素のうち面のどれかが void に接するもの。「内部」は 6 面ともシード。）

**要素番号の注記（10/5）。** 相方のデッキの要素番号は x→y→z の順ではない。デッキの節点から
求めた重心で確かめると、8³ は「x 逆順が最速、次に z、最後に y」、refine_deck の 16³ は
「z 最速、次に y、最後に x」。どちらも x→y→z を前提にした配置の鏡映と軸の入れ替え（等長）なので、
距離で決める層（neighbour_sign）とシード平均（compare_mesh）は正しい。正しくないのは次の二つ:
`compare_mesh --grid8 8 --track 220` は 220 を x = −0.125 に置くため 16³ の別の要素
[1639, 1640, …] を拾う（上の表は真の位置 x = +0.125 で取り直した値）。また `plot_3d --grid 8` の
図は軸が入れ替わっている。重心の列がない 8³ の CSV は、デッキから重心を付けて読むのが確実。

記録した表と、2・3 の画面表示をそのままクラウドのセッションに貼れば、こちらで修論4章に反映して
push する（IKMHIWI03 からブランチへ push する必要はない）。
