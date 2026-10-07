# 慶應のサーバーへの引継ぎ（2026年10月7日）

慶應で続ける計算（`KEIO_PLAN.ja.md` の WP1〜WP6）を、慶應の Linux サーバーの Abaqus で回すための
メモ。IKMHIWI03（Windows）でしか確かめていないことは、そう書いた。

---

## 1. 最初にすること（clone）

- 古い clone は使わない。2026年8月20日に履歴を書き換えた（著者の修正）ので、pull ではなく
  clone し直す。
- clone したら、git の名前と hook を入れる（`CLAUDE.md`「Git」の節）：
  ```
  git config --local user.name  "keisuke nishioka"
  git config --local user.email "128669518+keisuke58@users.noreply.github.com"
  cp scripts/pre-commit-no-ai-identity.sh .git/hooks/pre-commit && chmod +x .git/hooks/pre-commit
  cp scripts/commit-msg-no-ai-trailer.sh  .git/hooks/commit-msg  && chmod +x .git/hooks/commit-msg
  ```
- 作業ブランチ：`claude/plan-next-hxjjve`（10月7日の時点。master にはまだ入れていない）。
- Python 3 と numpy、scipy（入力ファイルを作るスクリプトと比べるスクリプトは、これだけで動く）。

## 2. git にあるもの・ないもの

| もの | 場所 | 慶應で |
|---|---|---|
| Abaqus の入力ファイルを作るスクリプト | `abaqus_composition/make_cube_inp.py`、`make_klempt_inp.py`、`make_implant_inp.py` | そのまま使える |
| UMAT / UMATHT / UEL（1ファイルに生成） | `abaqus_composition/make_umat.py --native`（fragment、`growth_from_phi.f`、`ansys_usermat/coupling/ecology_native.f` を取り込む） | そのまま使える |
| 結果の要約（Abaqus） | `abaqus_composition/results_1006/`、`results_1007/`（各ファイルの先頭に、作った時のコマンド） | 比較の基準 |
| 結果の JSON（ANSYS） | `ansys_usermat/apdl/results/` | 比較の基準 |
| WP2 の試作（軸対称、Python） | `keio_wp2/` | そのまま使える |
| `.odb`、`.dat` などの Abaqus の出力 | IKMHIWI03 の `F:\abaqus_work\`（git の外） | 持っていかない。要約から作り直せる |
| ANSYS の exe と、オリバーのソース | IKMHIWI03 の `F:\biofilm_upf_*`（git の外） | **持っていかない。** オリバーのソースはコミットしない |
| Felix の実装と理論の章 | IKMHIWI03 の `F:\felix_private\`（git の外） | **Felix に聞くまで持っていかない**（下の §5） |
| 一回限りの実行スクリプト | IKMHIWI03 の `F:\abaqus_work\_*.ps1`（git の外、Windows 用） | 移さない。手順は §3 の形で同じ |

## 3. Linux で Abaqus を回す

`abaqus_composition/run_comp.sh` は `run_comp.ps1 -Native` と同じことをする（Fortran 版の点モデル、
マテリアルサーバーなし）。**Linux ではまだ一度も動かしていない**（文法の確認だけ）。

```
export WORKROOT=$HOME/abaqus_work          # 既定。ホームが小さければ作業用のディスクに
python3 abaqus_composition/make_implant_inp.py abaqus_composition/wp2_imp_free_ph01.inp \
    --case 2sp_case6 --cons 1 --first-order --front 10 --blend 0.5 --dt 0.002 --T 1.0 \
    --base free --ph 0.1 --ph-local
abaqus_composition/run_comp.sh abaqus_composition/wp2_imp_free_ph01.inp 2sp_case6 4
python3 abaqus_composition/summarize_implant.py $WORKROOT/comp_wp2_imp_free_ph01/wp2_imp_free_ph01.dat \
    abaqus_composition/wp2_imp_free_ph01.json
```

1本の流れは、入力ファイルを作る → `run_comp.sh` → 比べる、の3段。比べるスクリプトは形で選ぶ：

| 形 | 作る | 比べる |
|---|---|---|
| 立方体（ANSYS と同じ問題） | `make_cube_inp.py --json <ANSYS の JSON>` | `compare_ansys.py <dat> <JSON>` |
| Klempt 2024 の 4.1 / 4.2 | `make_klempt_inp.py --case fig4_corner\|fig7_high\|fig7_low` | `compare_paper.py <dat> --case ...` |
| インプラント・歯 | `make_implant_inp.py` | `summarize_implant.py`、WP2 は `compare_wp2.py` |

Windows との違いで気をつけること：
- **`mp_mode=threads`。** UMAT・UMATHT・UEL は栄養の場を Fortran の module で共有している。
  CPU を MPI のプロセスに分けると共有できない。`run_comp.sh` はこれを付けている。
- Windows で必要だった `/libs:dll /threads`（`abaqus_v6.env`）は Linux では付けていない。
  コンパイラは Abaqus の設定どおり（`abaqus information=environment` で確かめる）。
- 点モデルの定数は環境変数 `BIOFILM_ECO_CASE`（`write_eco_cfg.py <case>` が書くファイル）。
  `run_comp.sh` が設定する。
- 長い実行は `nohup` か、サーバーのジョブ管理（ある場合）で回す。

## 4. 最初の確かめ（Linux で同じ結果が出るか）

1. `abaqus verify -user_std` で、ユーザーサブルーチンのコンパイルが通るか。
2. 短い1本（上の `wp2_imp_free_ph01`、IKMHIWI03 で 4 CPU・7 分）を回し、要約を
   `abaqus_composition/results_1007/wp2_imp_free_ph01.txt` と比べる。IKMHIWI03 では 1 CPU と
   4 CPU で全桁一致したので、コンパイラが違っても差は丸めの程度のはず。大きく違えば原因を調べてから
   先に進む。
3. 次に Klempt 2024 の 4.1（`make_klempt_inp.py --case fig4_corner`、20³）と `results_1006/` の
   同じ名前のファイル。

## 5. 守ること（引継ぎで特に気をつけるもの）

- **Felix のコードと理論の章**は、本人に聞かずに誰とも共有しない約束（`CLAUDE.md`）。慶應の
  サーバーに置くのも、先に Felix に聞く。リポジトリには照合の結果だけを書く。
- オリバーのソースはコミットしない。ANSYS の作業は IKMHIWI03 に残る。
- コミットにも PR にも、AI の署名を付けない。`.env` と PAT は共用のサーバーに置かない（push は SSH の鍵で）。
- ファイルは名前を指定してステージする（`git add -A` は使わない）。

## 6. IKMHIWI03 に残るもの・帰国（12月31日）までに

- ANSYS（オリバーの要素）の計算は IKMHIWI03 だけで回せる。帰国後も VS Code のトンネルで
  届く間は使える。
- WP3 の B（前線の項）は判定に落ちて始めていない（`ansys_usermat/apdl/FRONT_TERM_FIX.md`）。
  原因は前線の項の離散化で、慶應で続ける。
- WP3 の A（シードの成長、8³〜24³）は IKMHIWI03 で実行中。JSON は自動で push される。
- 帰国前に、`F:\abaqus_work\` のうち要約がまだ git に入っていない結果がないか確かめる。
