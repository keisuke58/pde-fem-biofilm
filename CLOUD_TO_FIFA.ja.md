# クラウドから fifa（慶應のサーバー）への指示

fifa の Claude は、チェーンを始める前に `git pull origin master` してこのファイルを読む。新しい項目が上に足される。

---


## 2026年10月16日（インプラントの扇形でメッシュ収束。急がない）

歯のデータの出典の報告と、扇形での Abaqus と Python の一致（4e−7）をありがとう。歯冠の STL は残し、出典と
CC BY-NC-SA 4.0 を `THIRD_PARTY.md` に書いた（クラウド）。OpenJaw の元の STL とそこから作ったメッシュは、
今まで通りリポジトリに入れない。

1. **メッシュ収束（Python の模擬実行で）：** `make_implant_uel_inp.py` で 2×6×8、4×12×16、6×18×24（nr×nt×nz）を
   `--bulge 0` と `--bulge 1` で。時間がかかりすぎるもの（目安 3 時間以上）は止めて、そこまでを報告。
2. **比べる量：** 内側の層（シード）の α − 1 と vM/E の平均、**角を除いた**内部の平均（上端 z = 2.0 と下端 z = 0 の列、
   扇形の切り口の列を除く）、最大 vM/E とその位置。
3. **最大 vM/E の位置を確かめる：** 今は上端の列（z = 1.875）に出ている。固定面（r = 2.05 で u = 0）と自由な上面の角は
   特異点になりうるので、メッシュを細かくすると最大値が上がり続けるかを見る。上がり続けるなら、立方体の角と同じく
   比べる量から外す。
4. 結果（`results/implant_<nr>x<nt>x<nz>[_b1].json`）と表を `keio/uel-implant` に push し、
   `KEIO_CHAIN_LOG.md` に報告。最後に細かいメッシュで 1 本だけ Abaqus でも回して、模擬実行と合うか確かめる。

## 2026年10月15日 深夜（歯の形状の準備：データの出典と、UEL をインプラントの形へ。急がない、24³ のあとで）

P3 で歯やインプラントの形に進む準備。IKMHIWI03 で A-18（硬い板と固定境界の比較）を回している。歯は固定境界
（バイオフィルムが乗る面で u = 0）で置き換えられる見込み。

**1. 歯のデータの出典を調べる（最優先、計算なし）**
- fifa の上で、歯の形状データを探す：`OpenJaw_Dataset`（`biofilm_3tooth_assembly.py` が `external_tooth_models/OpenJaw_Dataset/Patient_1`
  の `P1_Tooth_23/30/31` を読む）、`tier2b_real/_crown_real_watertight.stl` の元になったファイル、ほかの STL や歯のメッシュ。
  例：`find ~ / -iname '*openjaw*' -o -iname '*tooth*.stl' -o -iname '*crown*.stl' 2>/dev/null`。
- 見つけたものについて報告に：置き場所、ファイル数と大きさ、**出典**（論文、URL、DOI）、**ライセンス**（README、LICENSE、
  データの利用条件）、患者データなら匿名化と再配布の可否。分からなければ「分からない」と書く。
- **データそのものはコミットしない。** このリポジトリは公開。ライセンスが再配布を認めていると確かめられるまで、
  リポジトリに入れるのは出典の記録だけ（`THIRD_PARTY.md` に追記する案を報告に書く。追記はクラウドでやる）。
- `tier2b_real/_crown_real_watertight.stl` はすでに公開リポジトリに入っている。出典が分かったら、ライセンス上問題ないかも書く。

**2. UEL をインプラントの扇形に載せる（1 のあと）**
- 形は `abaqus_composition/make_implant_inp.py` と同じ（半径 2.05 mm の首、厚さ 0.25、高さ 2.0、90° の扇形、六面体）。
  内側の面（チタン）は u = 0、扇形の切り口は対称面、栄養は上面（z = 2.0）で c = 1、最初のシードは内側の 1 層で φ = 0.9。
- `p3/uel/make_implant_uel_inp.py` を作る：`make_implant_inp.py` のメッシュ（節点と六面体）を使い、要素と手順は
  `p3/uel/make_seed_inp.py` と同じ（UNSYMM、VARIABLES = 128、INIT と GROW の 2 段、`--md`、ν = 0.49）。
- まず `p3/uel/emulate_abaqus.py` と同じやり方（ctypes）で Python から回し、そのあと Abaqus で回して一致を確かめる
  （立方体のときは 4〜5 桁一致）。小さいメッシュ（例：nr 2、nt 6、nz 8）から。
- 報告：シードの α − 1 と vM/E の平均、最大値とその位置、Python と Abaqus の差。`keio/uel-implant` ブランチに push。
- `--bulge 1`（歯冠に近い形）は、扇形で一致を確かめてから。

## 2026年10月15日 夜（UEL の 24³ と、Abaqus での確認。急がない）

ANSYS（IKMHIWI03）と UEL の独立な比較（`p3/uel/README.md` の「Cross-code check with ANSYS」）で、シード内部の vM の差が
8³ −2.7 %、16³ −1.2 % と縮んでいる。ANSYS は 24³・32³ で 6.5e−3 に収束した。UEL の 24³ がほしい（UEL は Linux 用なので
IKMHIWI03 では回せない）。

1. `git pull origin master`。UEL の共有ライブラリを作る（`p3/uel/driver.py` の docstring の 3 行。gfortran）。
2. `python3 p3/uel/make_seed_inp.py 8 /tmp/sf`、
   `BIOFILM_ECO_CASE=/tmp/sf/uel_eco.txt python3 p3/uel/seed_frozen_reference.py 24`（数時間かかるかもしれない。
   nohup で。メモリが足りなければ止めて報告）。
3. 出力 `p3/uel/results/seed_frozen_24.json` を `keio/uel-check` ブランチで push。報告に、シードの α − 1 平均と vM 平均、
   半径 0.2 以内の要素の vM 平均（`elements` の重心から、中心 (0.5, 0.5, 0) からの距離 ≤ 0.2）を書く。
   期待：半径 0.2 以内の vM は ANSYS の 24³ の 6.536e−3 と 1 % 程度で合う。
4. 時間があれば、下の 10 月 15 日の項（Abaqus で柱とシード）も。

## 2026年10月15日（P3 の計算を Abaqus で：次に fifa を使うときに）

`KEIO_SERVER_HANDOFF.ja.md` §7 を読む。急がない（11月以降でよい）。やること：
1. `git pull`（master）。`python3 p3/uel/make_abaqus_source.py`、`python3 p3/uel/make_column_inp.py 20 p3/uel/abaqus_in`、
   `python3 p3/uel/make_seed_inp.py 8 p3/uel/abaqus_in`。
2. `export BIOFILM_ECO_CASE=$PWD/p3/uel/abaqus_in/uel_eco.txt BIOFILM_UEL_GP=$PWD/p3/uel/abaqus_in/<job>_gp.txt`、
   `abaqus job=uel_column user=p3/uel/unified_uel_abaqus.for cpus=1 interactive`（`p3/uel/abaqus_in` で）。
3. 期待：柱は x = 0 で NT11 0.7753、NT12 0.9357（IKMHIWI03）。`uel_seed8` の φ > 0.5 の要素の vM 平均 0.01300、最大 0.01637。
4. 結果（`.sta` と `_gp.txt` の gzip）を `keio/uel-check` ブランチで push し、合ったかを書く。リンクで止まったら、
   最初のエラーを書いて止める（Windows では `/libs:dll` が要った）。

## 2026年10月8日 夜（lab の設定が済んだ）

- ユーザーが fifa で `git remote add lab ...`、`git ls-remote lab`、`git push lab master:main` を実行し、通った
  （lab の `main` ができた）。SSH の鍵と組織の書き込み権限は足りている。
- fifa の Claude は `mmc-research-group` への操作が権限で止められているので、**lab への push はユーザーが手で行う**。
  fifa の Claude は、master を pull したあとで lab が遅れていれば `KEIO_CHAIN_LOG.md` に「`git push lab master:main` が必要」と
  書くだけでよい（自分では実行しない）。
- 10月8日、ユーザーが `git push lab origin/master:main` で #64 まで lab に送った（同期済み）。
- **lab に送るときは、ローカルの master を通さない**（fifa のローカル master は clone 以来古いまま）：
  `git fetch origin master && git push lab origin/master:main`（fast-forward、force は不要）。

---

## 2026年10月8日（村松研の git ができた）

ユーザーが村松研の組織に **`mmc-research-group/nishioka-biofilm-fem`**（private、空、既定ブランチ `main`）を作った。
慶應の作業はここにも置く。

1. **remote は `lab` という名前で足す。`origin`（公開の pde-fem-biofilm）は変えない。** GitHub の画面が出す
   `git remote add origin ...` はそのまま打たないこと（origin が上書きされ、チェーンの push 先が変わる）。
   ```
   git remote add lab git@github.com:mmc-research-group/nishioka-biofilm-fem.git
   git push lab master:main
   ```
2. **SSH の鍵で push する**（HTTPS とトークンは共用サーバーに置かない）。fifa の公開鍵（`~/.ssh/id_ed25519.pub`
   など）が、ユーザーの GitHub アカウントに登録されていて、そのアカウントが組織 `mmc-research-group` で書き込み
   できるかを `ssh -T git@github.com` と `git ls-remote lab` で確かめる。通らなければ、何が足りないかを
   `KEIO_CHAIN_LOG.md` に書く（鍵の登録はユーザーがする）。
3. **当面の運用：** チェーンの結果は今までどおり `origin` の `keio/*` に push し、master へのマージはクラウドが
   する。fifa は master を pull したら、`git push lab master:main` で lab の `main` を master に合わせる。
   lab だけに置くもの（未発表の原稿など）の分け方は、1月の村松先生との面談で決める。
4. **どちらの remote にも入れないもの：** Felix のコードと理論、オリバーのソース（`CLAUDE.md`）。

## 2026年10月8日 朝（ac94d34 への返事）

- 1・4・5 の対応を確認した。master に入れた（`run_comp.sh` は master の scratch の扱いと fifa の seed の
  コピーの両方が入った形で自動マージされた。チェーン 20261008b が終わってから pull すればよい）。
- hand-off §1 の作業ブランチの記述は、master ではすでに直してある（「作業ブランチ：master」）。
- `~/.config/gh/hosts.yml` の gh の認証は、ユーザーに `gh auth logout` を頼んだ。チェーンは使わないので急ぎではない。

## 2026年10月8日（PR #58 の未解決 2 件と、報告先・追従の判断）

1. **報告先：PR コメントをやめる。** チェーンの結果は `KEIO_CHAIN_LOG.md` と要約ファイル（`results_keio_*`）だけに書く。
   `gh` の認証は fifa（共用サーバー）に置かない。push は SSH の鍵で `keio/*` のブランチへ。
   **master へのマージはクラウドが行う**（定期確認で `keio/*` を拾う）。fifa から master には直接入れない。
2. **master に追従する。** チェーンを始める前に毎回 `git pull origin master` し、master を土台にして作業ブランチを作る。
3. **`keio_dispatch.py` の転送は作らない。** 次に回すものは、クラウドが `scripts/keio_runs/` に manifest を置いて
   master に入れる。fifa は pull したときに、まだ回していない manifest があれば起動する（`KEIO_CHAIN_LOG.md` に
   起動と結果を記録するので、回したかどうかはそこで分かる）。git だけで回るので、新しい認証は要らない。
4. **requirements.txt は fifa の版に緩めてよい**（numpy 2.4.2、matplotlib 3.10.8、pytest 9.0.2）。上げるのは、
   fig4_corner 一式を回し直して全桁一致を確かめてから。
5. **hand-off §1 の古い記述：** どの文かを `KEIO_CHAIN_LOG.md` に書いておくこと。クラウドが直す。
6. 実行中の **20261008b**（fig4_corner の残り 5 本）はそのまま続けてよい。終わったら結果を `keio/*` に push。
