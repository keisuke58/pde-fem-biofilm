# 慶應のサーバーへの引継ぎ（2026年10月7日）

慶應で続ける計算（`KEIO_PLAN.ja.md` の WP1〜WP6）を、慶應の Linux サーバーの Abaqus で回すための
メモ。IKMHIWI03（Windows）でしか確かめていないことは、そう書いた。

---

## 1. 最初にすること（clone）

- 古い clone は使わない。2026年8月20日に履歴を書き換えた（著者の修正）ので、pull ではなく
  clone し直す。
- **クラウドからの指示は `CLOUD_TO_FIFA.ja.md`**（pull のたびに読む）。
- サーバーは **fifa**。clone したら、まず `bash scripts/setup_keio_server.sh` を1回走らせる（git の名前と
  hook を入れ、Python・numpy・scipy・Abaqus のコマンドを確かめる）。中身は次と同じ：
  ```
  git config --local user.name  "keisuke nishioka"
  git config --local user.email "128669518+keisuke58@users.noreply.github.com"
  cp scripts/pre-commit-no-ai-identity.sh .git/hooks/pre-commit && chmod +x .git/hooks/pre-commit
  cp scripts/commit-msg-no-ai-trailer.sh  .git/hooks/commit-msg  && chmod +x .git/hooks/commit-msg
  ```
- 作業ブランチ：`master`（10月7日に PR #56 をマージした）。fifa では `git checkout master && git pull` で最新になる。
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

## 6. IKMHIWI03 は12月中旬まで（10月7日に分かった）

IKMHIWI03 が使えるのは**12月中旬まで**。そのあとは慶應のサーバーだけで、**慶應には少なくとも
半年ほど（2027年6月ごろまで）ANSYS がない**。そのあいだ、オリバーの要素（NEM）の計算はどこでも
できない。慶應に ANSYS が入ったあとにオリバーの要素を使えるかは、オリバーに聞く（ソースは
持っていかない約束のまま）。

12月中旬までに IKMHIWI03 で済ませること：
- **ANSYS の計算を全部終える。** 慶應の論文で使う NEM 側の数字（WP3：メッシュ 8³〜24³、
  時間刻み、シード表面の α の段差）は、ここで出して JSON で git に入れる。あとから足せない。
  - WP3 の A（シードの成長）：実行中。JSON は自動で push される。
  - WP3 の B（前線の項）：判定に落ちて始めていない（`ansys_usermat/apdl/FRONT_TERM_FIX.md`、
    原因は前線の項の離散化）。12月中旬までに直して回すか、ANSYS 側は諦めて Abaqus だけにするかを
    決める。決めるのは11月中が目安（回す時間が要る）。
  - 足りない実行がないか、WP3 の論文の図と表の案から逆算して確かめる（11月中）。
- **結果を git に入れる。** ANSYS は `export_runs_json.py` で JSON に、Abaqus は要約の `.txt` に。
  `F:\` にしかない結果が残っていないか、最後に確かめる。
- 慶應のサーバーで §4 の確かめを、IKMHIWI03 が使えるうちに済ませる（食い違ったときに
  IKMHIWI03 で回し直して比べられる）。

進め方（10月7日の案）：

| 時期 | IKMHIWI03 でやること |
|---|---|
| 〜10月12日 | 修論だけ。表面の点の修正の影響を見て、修論での扱いを決める |
| 10月13〜20日（旅行の前） | 前線の項を使わない WP3 の実行を登録する：シード表面の α の段差 × β × メッシュ（8³〜24³）× 時間刻み。旅行中に回って自動で push される |
| 旅行中 | 機械は回り続ける。止まったら旅行のあとで回し直す |
| 11月2日〜11月末 | 村松先生と WP3 の範囲を決め、図と表の案から足りない実行を出して登録する。前線の項を ANSYS で直すのは**1週間だけ**試し、だめなら Abaqus だけにする |
| 12月1日〜中旬 | 最後の実行、JSON の書き出し、`F:\` にしか残っていない結果がないか確かめる。慶應のサーバーで §4 の確かめ |

WP3 の NEM と Galerkin の比較は、前線の項がなくても書ける（β の拡散とシードの成長だけで、
Rudolf et al. 2025 が比べたのも拡散）。前線の項は WP3 の必須ではない。

12月中旬のあと：
- 計算は Abaqus（慶應のサーバー）だけ。KEIO_PLAN の WP1、WP2、WP4〜WP6 はもともと Abaqus。
- WP3 で NEM の計算を足したくなったら、NEM の演算子（Rudolf et al. 2025 の Eq. 4–13）を Python で
  書くしかない。オリバーのソースは持っていかない。

## 7. P3 の計算を Abaqus で（2026年10月15日）

方針（`P3_OUTLINE.ja.md`）：P3 の計算の主役は**方式 (c)**（生きている量で成長、κ）。その裏付けと次の段階が
**一つの原理のモデル**（`p3/unified_model.pdf`）。どちらも慶應の Abaqus で回せる形にした。

**方式 (c)（UMAT、キューブ）：**
- `python abaqus_composition/make_cube_inp.py OUT.inp --json <ANSYS の JSON> --theta-json CH.json --growth --kappa 8.135097e-3`
  （mode 9。`--growth` なしは mode 7、(a) の比較用）。点モデルは `BIOFILM_ECO_CASE`（`write_eco_cfg.py --theta-json ... --phi-init ...`）。
- θ は Tmcmc202601 の `final_theta_MAP/CH.json`（6e18e0b 以降、生データ・2 通路）。初期値は同じ commit の `phi_init.json`
  （合計 1.001 のまま渡す。推定と同じ）。κ = 8.135097e−3 は新しい θ で (c) のシードの成長が (a) と同じになる値。
- ANSYS（IKMHIWI03）での (c) と (a) の vM の差：8³ で −7.0 %（16³、24³ は B4）。Abaqus では同じプログラムの中の比で比べる（A-14）。

**一つの原理のモデル（UEL、`p3/uel/`、README）：** Abaqus 2024 で動き、模擬実行と 4〜5 桁一致した（IKMHIWI03、15日）。
- ソース：`python p3/uel/make_abaqus_source.py` → `unified_uel_abaqus.for`（全部固定形式、ifort）。
- 入力：`make_column_inp.py`（柱の検査）、`make_seed_inp.py n outdir [--small] [--nu 0.49 --md] [--antibiotic] [--cur]`。
- 環境変数：`BIOFILM_ECO_CASE`（点モデル）、`BIOFILM_UEL_GP`（積分点の結果のファイル、絶対パス）。`cpus=1`。
- Abaqus で分かった注意（入力の生成にもう入っている）：
  - 自由度 11、12（φ、c）は `*STATIC` では使えない → 過渡の `*COUPLED TEMPERATURE-DISPLACEMENT`（`DELTMX` なし = 固定の刻み）、
    `*USER ELEMENT, UNSYMM`；
  - 過渡の連成には比熱が要る → 同じ節点にダミーの C3D8T（材料 1e−9）；
  - 補給の面（c = 1）は INIT のステップから与える（成長のステップの 1 増分目で切り替えると、8³ で成長が 12 % 下がる）；
  - Windows では Fortran の OPEN/WRITE のため `/libs:dll`（`run_abaqus.ps1 -LibsDll`）が要った。Linux の ifort では未確認。
- 期待値（模擬実行、`p3/uel/emulate_abaqus.py`）：柱 NT11 0.776、NT12 0.936。シード 8³ の vM 平均（φ > 0.5）0.01300（有限ひずみ）。
- **慶應で最初にすること：** fifa で `uel_column` と `uel_seed8` を回し、上の値と比べる（Linux の ifort で同じになるか）。

一つの原理のモデルを論文にするかは、Felix と Meisam に相談してから決める（Felix の未発表の 2 菌種の仕事と重なるため）。

## 8. 結合モデルの Abaqus 要素（abaqus_living、2026年10月24日）

慶應のサーバーでも Claude を SSH ＋リモートコントロールで使える（ユーザー、10月24日）。IKMHIWI03 の続きの Abaqus は
ここで回す。要素は `abaqus_living/living_uel.f90`（ガラーキン、前線の項は SUPG、新しい点のルール layout 2/3）。
検証はすべて IKMHIWI03 で済んでいる（製造解、Python との GCI、迷路、`ansys_usermat/apdl/IKMHIWI03_TO_CLOUD.ja.md`
A-46〜A-78）。

1. **ソース**：`python3 abaqus_living/make_abaqus_source.py` が固定形式の `abaqus_living/living_uel_abaqus.for` を書く
   （git にもある。`living_uel.f90` を変えたら作り直す）。点モデルは `ansys_usermat/coupling/ecology_native.f` が中に入る。
2. **環境変数**：`BIOFILM_ECO_CASE=<repo>/abaqus_living/eco_u3_a100.txt`（2 種、Fritsch case I、α\*_0 = 100）、
   `BIOFILM_UEL_GP=<作業ディレクトリ>/living_gp.txt`（ガウス点の記録）、`BIOFILM_UEL_EVERY=400`（記録の間隔、増分）。
   **CPU は 1 つ**（layout 3 と記録は module の配列を共有する。threads でも安全とは確かめていない）。
3. **最初の確かめ**（コンパイラの違い）：
   ```
   python3 abaqus_living/make_living_inp.py 16 sep lv16qs2_sep_ab.inp --quarter --supg --dt 0.00125
   abaqus job=lv16qs2_sep_ab user=abaqus_living/living_uel_abaqus.for cpus=1 interactive
   python3 abaqus_living/summary_living.py $BIOFILM_UEL_GP 16 --json lv16qs2_sep_ab_summary_fifa.json
   ```
   IKMHIWI03 の値（`abaqus_living/lv16qs2_sep_ab_summary.json`、T\* = 1）：種の Σφ̄ 0.1494、c 0.9137、a 0.9851、
   α−1 8.554e−4、von Mises 6.168e−4 Pa。丸めの程度（1e−6 相対）を超えて違えば先に進まない。
   （`--quarter --supg --dt 0.00125` が IKMHIWI03 のデッキと同じかは、`diff` で `abaqus_living/lv16qs2_sep_ab.inp` と比べる。）
4. **次の計算**：4 種（Klempt 2026 Table 2 case 4、抗菌薬で耐性の種 4 が残る例）。点モデルのファイルは
   `abaqus_composition/eco_klempt2026_t2c4_a100.txt`（行列 A は**式 19 の 1/2 なし**：論文の図はこれでしか再現しない、A-77）。
   Python の参照値（64³、TVD）：分けた例 Σφ̄ 0.762、種 1・2 の割合 0.193、α−1 9.66e−4、vM 6.22e−4 Pa（A-78 に GCI）。
   要素は eco_nact = 4 をそのまま扱う（Python と同じ Fortran の点モデル）。layout 1（分けた例）は種 2 を角に置くので、
   4 種では「種 1・2 が中央、種 3・4 が角」にする変更が要る（`start_state`、Python の `GROUP_A`）。
5. **迷路**：`make_maze_inp.py --tend 60 --supg [--ab] [--inherit]`、比べるのは `maze_compare_a67.py`
   （Python の結果は `ansys_usermat/apdl/results/2026-10-maze/`）。layout 3（受け継ぐ）は 1 増分ごとに全ガウス点を
   見るので、layout 2 の約 5 倍遅い。
6. **結果の置き場**：要約の JSON だけ git に入れ（`living_gp.txt` は大きい）、`keio/*` のブランチに push する。
