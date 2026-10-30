# クラウドから fifa（慶應のサーバー）への指示

fifa の Claude は、チェーンを始める前に `git pull origin master` してこのファイルを読む。新しい項目が上に足される。

---

## 2026年10月29日 夜（追加：今の列が終わってから 2 本、論文の表 2 用）

成長則なし（今の式）で、64³・外側の刻み 0.0025 の 2 本。論文の表 2（今は刻み 0.025）を差し替え、成長則の列の 64³ の比較相手にもなる
（この 2 本だけクラウドにも master にもない。IKMHIWI03 にはあるが JSON がない）。
```
P="--fast --tvd --wall-grad mirror --gc 1 --gab 1 --alpha0 100 --front 1 --tab 0.333333 --phi-min 0.001 --v2 --fortran --ansys-bc --inherit-living --dtc 0.0025 --nf 64"
python3 coupled_living_consumption.py $P --only separate_ab
python3 coupled_living_consumption.py $P --mix-equal --only mixed_ab
```
確かめ：種の平均の α−1 と vM が修論の表 `tab:ch5-gci` の 64³（8.931 / 6.000、8.04 / 2.82 ×1e−4）と 0.1 % 以内で合うこと。
JSON を `results/2026-10-3d/` に追加して push、`KEIO_CHAIN_LOG.md` に Σφ̄ᵢ、A の割合、c、α−1、vM を。

## 2026年10月29日 夕方（新しい依頼：生きている量で成長する則の格子列。慶應での続き、修論には入れない）

`ansys_usermat/apdl/coupled_living_consumption.py` に `--growth-living`（既定は切）を足した：式 36 を dα/dt = k_α φ_3D Σφ̄ᵢ にする
（生きている菌だけが成長する、Wanner & Gujer 1986 と同じ考え）。今の成長則との比較の格子列を fifa で回す。

**準備**：`git pull origin master`。`pip install pyamg`（無ければ）。`--fortran` は `ecology_native.f` を ctypes で使う
（`ansys_usermat/coupling/ecology_native_ctypes.py` が自分でビルドする。できなければ `--fortran` を外して、外したことを記録）。

**共通のオプション**（修論の本番と同じ）：
```
P="--fast --tvd --wall-grad mirror --gc 1 --gab 1 --alpha0 100 --front 1 --tab 0.333333 --phi-min 0.001 --v2 --fortran --ansys-bc --inherit-living --dtc 0.0025"
```
**回すもの**（上から 1 本ずつ、nohup か tmux で切り離して）：`--growth-living` を付けて、4 つの例 × 3 格子
- 例：`--only separate`、`--only separate_ab`、`--mix-equal --only mixed`、`--mix-equal --only mixed_ab`
- 格子：`--nf 32`（指定なし）→ `--nf 64` → `--nf 128`（128³ は長いので最後。1 本の時間を記録）
- 比べる相手（成長則なし）は既にある：32³〜128³ は `results/2026-10-3d/` の `..._front1_dtc0.0025_..._inhliv_...json` と
  修論の表 `tab:ch5-gci`（抗菌薬ありの分離・混合（等量）の 32³〜128³ は A-89 の値）。

**追加（同じ日の夜）：案 A も同じ格子列で**。`--growth-living-energy`（エネルギーの結合項を −k_α α φ Σφ̄ᵢ にした形：式 36 と
式 34 の生成項 k_α α の両方に Σφ̄ᵢ がかかる）。上の `--growth-living`（案 B、式 36 だけ）と同じ 4 例 × 32³/64³/128³ を、
案 B の各格子が終わったら同じ格子で回す（順番：32³ の B → 32³ の A → 64³ の B → 64³ の A → 128³ の B → 128³ の A）。
JSON は `_growlivE` 付き。

**報告**（`KEIO_CHAIN_LOG.md`、結果の JSON は `results/2026-10-3d/` に `_growliv` 付きで、`keio/uel-implant` に push）：
例ごと・格子ごとに、種の平均の Σφ̄ᵢ、c、α−1、vM と、成長則なしとの比（α−1 と vM）。3 格子がそろったら GCI（Celik 2008、Fs 1.25）。
UEL・Abaqus の作業は止めてよい（終わっている）。

## 2026年10月29日（28 日の結果への返事。論文用、急がない）

29 日の報告（NaN は bulge 1 の最初の増分の反復 4 で残差に入る、T/240 なら通る、bulge 0 は PASS、圧力で止まる成長則 4/4 PASS）を受けて。
UEL のガード（残差と状態が非有限なら PNEWDT = 0.25）は妥当なので、そのまま使う（master に取り込んだ）。

- **直したこと（クラウド）**：`make_implant_uel_inp.py --ninc n` が Abaqus のデッキの GROW のデータ行に届いていなかった。今は届く
  （`--ninc 240` で `0.0009895833333, 0.2375, 9.6639e-07, 0.0009895833333`、手で書いた行と同じ）。`git pull origin master` のあとで使う。

上から 1 本ずつ：
1. **同じ刻みどうしで比べる**：エミュレーション `make_implant_uel_inp.py --run --nr 6 --nt 18 --nz 24 --seed-layers 3 --bulge 1 --ninc 240`
   と、ガード付きの UEL の Abaqus（同じオプションでデッキを作り直す）。`compare_implant_gp.py` の差が 1e−5 以下になるかを記録。
   なれば、28 日の FAIL は刻みの違いだけ。
2. 1. が通らなければ、両方の増分ごとの記録（エミュレーションの log と .sta）の最初に違う増分を探して記録。原因の推測は書くが、直さない。
3. 余裕があれば：bulge 1、T/60 の最初の試行で Newton の反復 3 の修正量が 0.3 に跳ねる理由の手がかり（その要素の J、det J の最小、
   bulge で歪んだ要素か）。記録だけ。

copaam は今は使わない（ユーザー）。

## 2026年10月28日（6×18×24、bulge 1 の Abaqus が NaN になる件の切り分け。急がない、論文用）

28 日の報告（`KEIO_CHAIN_LOG.md`、2×6×8 は PASS、6×18×24・bulge 1・`--seed-layers 3` は全点 NaN）への返事。修論には入らない（慶應での続き、統合論文用）。
**原因が分かっても UEL のソースは直さず、結果をそのまま記録する**（直すのはクラウド）。上から順に、1 本ずつ。

1. **どの増分で NaN が入るか**：同じデッキで、Gauss 点の出力を最初の 5 増分は毎増分にして回す。NaN が最初に出る増分、
   要素番号、その増分の .sta の行（反復数、カットバックの有無）を記録。最初の増分で入るなら、その要素の φ、c、α−1 の直前の値も。
2. **刻みを小さくする**：`make_implant_uel_inp.py abaqus_in --nr 6 --nt 18 --nz 24 --seed-layers 3 --bulge 1 --ninc 240`
   （dt = T/240 = 9.9e−4、エミュレーションのカットバック後と同じ）。通れば `compare_implant_gp.py` でエミュレーションと比べ、
   種の vM/E と最大値の差を記録。通らなければ 1. と同じ情報を。
3. **bulge の切り分け**：`--bulge 0` で同じ 6×18×24（`--seed-layers 3`、`--ninc` なし）。2×6×8 が bulge 0 で通っているので、
   格子と bulge のどちらが効くかを分ける。
4. 8×24×32 のエミュレーション（17:44 開始）は、そのまま続けて結果を記録（止まったら増分番号と stderr の最後の 30 行）。
5. 恒常圧の Abaqus の確認と Lean / Mathlib の確認は、1.〜4. の後で。

`.msg` の「negative eigenvalues」の警告は 2×6×8 でも出るので、それだけでは原因と見なさない。

## 2026年10月17日（重要：UEL の形状関数の勾配の誤りを直した。扇形の計算は直した版で）

`p3/uel/unified_uel.f90` は形状関数の勾配を ∂N/∂x = Ji^T ∂N/∂ξ で計算していた。正しくは Ji ∂N/∂ξ（Ji = Jm^{-1}、
Jm(a,i) = ∂x_i/∂ξ_a）。軸に沿った直方体（立方体の計算すべて）では同じなので、立方体の結果と ANSYS との比較（A-16）は
変わらない。**扇形（インプラント・歯）は要素が回転しているので、10月16日までの扇形の結果はすべて誤り**で、クラウドで
直した版でやり直した（`p3/uel/README.md` の「Shape gradient of rotated and skewed elements」）。

- 下の 10月16日の依頼（8×24×32、Abaqus の確認）は、**`git pull origin master` のあとライブラリと Abaqus 用のソースを
  作り直してから**回す（`p3/uel/driver.py` の 3 行、`make_abaqus_source.py`）。直す前の版で回したものがあれば、捨てて
  回し直す。比べる相手の json はクラウドで作り直してある。
- 以前の「Abaqus と Python が扇形で 4e−7 で一致」は、両方が同じ誤りを持っていたので、コードの一致の確認としてだけ有効。
- 収束しなかった 3 本（成長則なし、栄養を外側から）は、この誤りが原因だった。直した版ではすべて収束する。

## 2026年10月16日 夜（Lean と Mathlib：成長則の証明の確認。計算の合間に、ただし手順は省かない）

クラウドで、恒常圧の成長則（陰的）が 1 つの材料点では p_h を超えないことを Lean で証明した（`p3/uel/lean/HomeostaticGate.lean`、
Lean 本体だけ、クラウドで確認済み）。実数そのもので証明する Mathlib 版（`p3/uel/lean/mathlib/`）は、クラウドではディスクが
足りず**未確認**。論文の付録に使う前に、fifa で次の手順どおりに確認してほしい。**どの手順も、結果をそのまま記録する。証明や
文を直さない**（直すのはクラウド）。

**0. 準備**（空き 8 GB 以上、ネットワーク。`df -h` で確認して記録）
```
git pull origin master
curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
source ~/.elan/env
cd p3/uel/lean/mathlib
cat lean-toolchain                 # leanprover/lean4:v4.15.0 のはず
```

**1. Mathlib を取って、証明をビルドする**
```
lake exe cache get                 # ビルド済みの Mathlib（約 5 GB）。失敗したら 2 回まで再試行、だめならログを残して止める
time lake build                    # エラー 0 で終わること。警告は全文を記録
```
**合格：** `lake build` の終了コード 0、`error` の行が 0。

**2. 公理の確認（`sorry` やずるがないこと）**
```
lake env lean CheckAxioms.lean
```
**合格：** 4 つの定理それぞれが `propext`、`Classical.choice`、`Quot.sound` だけに依存し、`sorryAx` が出ない。

**3. 否定の試験（本当に検査しているか）**
```
lake env lean ShouldFail.lean ; echo "exit=$?"
```
**合格：** **エラーになること**（偽の主張「陽的でも超えない」を、陰的の証明の写しで示そうとしている）。エラーにならなければ、
検査が働いていないので、そこで止めて報告。

**4. 本体だけの版も同じ環境で**
```
lean ../HomeostaticGate.lean ; echo "exit=$?"     # exit=0、エラー 0
```

**5. 報告**（`KEIO_CHAIN_LOG.md` に、`keio/lean-check` に push）
- Lean と Mathlib の版（`lean --version`、`lake-manifest.json` の mathlib の rev）
- 手順 1〜4 の合否と、出力の要点（警告、エラーの全文、`#print axioms` の出力の全文）
- かかった時間とディスク使用量
- `.lake/` と `lake-manifest.json` 以外の生成物はコミットしない（`.lake/` は `.gitignore` 済み。`lake-manifest.json` は
  コミットしてよい：版を固定するため）

エラーが出た場合（Mathlib の名前の変更、`norm_num` が閉じない、など）は、エラーの全文と行番号を報告に貼るだけでよい。

## 2026年10月16日 深夜・追記（恒常圧の成長則に陰的な帰還を足した。Abaqus の確認は陰的でも回せる）

下の項の帰還（増分の始めの圧力、陽的）に加えて、増分の終わりの圧力で帰還を評価する陰的な扱いを入れた
（`p3/uel/unified_uel.f90` の `implicit_gate`：Gauss 点ごとの局所 Newton と、陰関数の微分を含む一貫した接線）。
切り替えは **PROPS(20) = 形 + 10**（10/11/12 が 0/1/2 の陰的版）、デッキは `make_implant_uel_inp.py --p-ref 0.05 --gate-implicit`
（タグ `_ph0.05_gi`、模擬実行の結果は `results/implant_2x6x8_ph0.05_gi.json`）。既定（PROPS(20) < 10）は変えていないので、
下の項の陽的な確認はそのまま。Abaqus の確認を陰的でも回すなら、下の手順 2 の `--p-ref 0.05` を `--p-ref 0.05 --gate-implicit`
に、ジョブ名と比較先の `.json` を `_ph0.05_gi` に読み替えるだけ（1 増分の Newton は 3 回、2 次収束。帰還が 0 と 1 の間の
Gauss 点では p が p_h を下から追うので、刻みが違っても差は小さいはず）。局所 Newton が収束しない要素は PNEWDT = 0.5 で
カットバックを頼む（模擬実行では一度も起きていない）。詳しくは `p3/uel/README.md` の「Implicit gate」。

## 2026年10月16日 深夜（恒常圧の成長則を UEL に入れた：Abaqus で 1 本確かめる。急がない、8×24×32 のあとで）

WP2（`keio_wp2/`）の恒常圧の成長則 α̇ = k_α φ min(1, max(0, 1 − p/p_h))（p は圧縮を正とする圧力）を
`p3/uel/unified_uel.f90` に入れた。PROPS(19) = p_h/E、PROPS(20) = 形（0：p_h 一定、1：p_h (φ² + f)、2：1/(1 + p/p_h)）。
**PROPS が 16 個のままなら off で、今までの結果はビット単位で同じ**（2×6×8 の bulge 0/1、凍結シード 8³ で差 0）。
増分の始めの圧力（SVARS(113:128) の p、VARIABLES = 128 が要る）を使う陽的な扱いなので、接線はそのまま（成長の項に
g を掛けるだけ）。1 要素の検証は `p3/uel/homeostatic_test.py`、模擬実行の結果は `p3/uel/README.md` の
「Homeostatic-pressure growth law」、デッキは `make_implant_uel_inp.py --p-ref 0.05`（タグ `_ph0.05`）。

1. `git pull origin master`、ライブラリと `unified_uel_abaqus.for` を作り直す（`python3 make_abaqus_source.py`。
   固定形式への変換は gfortran で確かめてある）。
2. Abaqus で 2×6×8、bulge 0、p_h/E = 0.05、外面は自由（速い。模擬実行は 8 秒）：
   ```
   cd p3/uel && python3 make_implant_uel_inp.py abaqus_in --p-ref 0.05
   export BIOFILM_ECO_CASE=$PWD/abaqus_in/uel_eco.txt BIOFILM_UEL_GP=$PWD/abaqus_in/uel_implant2x6x8_ph0.05_gp.txt
   (cd abaqus_in && abaqus job=uel_implant2x6x8_ph0.05 user=../unified_uel_abaqus.for cpus=1 interactive)
   python3 compare_implant_gp.py abaqus_in/uel_implant2x6x8_ph0.05_gp.txt results/implant_2x6x8_ph0.05.json
   ```
   期待：増分の刻みが同じ（60 増分、カットバックなし）なら、Gauss 点の α − 1 と vM は 1e−6 以下で合う（PASS）。
   帰還は増分の始めの圧力で決まるので、Abaqus が別の刻みで切ると差は時間離散化の分（1〜2 %）になる。その場合は
   `.sta` の増分列を報告に貼る。
3. 時間があれば外面を固定した方（`--p-ref 0.05 --outer fixed`、タグ `_ph0.05_of`、`results/implant_2x6x8_ph0.05_of.json`）
   も同じ手順で。こちらは閉じ込めが強く、帰還の効きが大きい（上面の角の vM/E が 0.36 → 0.08）。
4. 6×18×24 の固定シードに帰還を入れた模擬実行（`--nr 6 --nt 18 --nz 24 --seed-layers 3 --p-ref 0.05`、bulge 0 と 1。
   1 本 35 分程度）は、クラウドでは回していない。8×24×32 のあとに時間があれば。
5. 結果（`_gp.txt` の gzip、`.sta`、PASS/FAIL）を `keio/uel-implant` に push し、`KEIO_CHAIN_LOG.md` に報告。

## 2026年10月16日 夜（bulge 1 の 6×18×24 が収束した。自動の増分制御を入れた。急がない）

夕方の項の未解決（固定シード、bulge 1 の 6×18×24 が最初の増分で止まる）が片付いた。詳しくは `p3/uel/README.md` の
「Implant sector: mesh convergence」の後半。

- **原因：** 要素のつぶれでも歯冠の形でもない。最初の成長増分で、出発状態（α = 1、u = 0、空の層は床の剛性 f = 1e−3）の
  連成接線が、増分 dt = T/60 のちょうど近くでほぼ特異になる（1/dt と反応の項の固有値がそこで 0 を横切る）。Newton の
  最初の更新が dt/4 で 1e−3 mm、dt/2 で 2.4e−3 mm なのに、dt では 0.97 mm（0.04 mm の要素に）で 94 要素が一度に裏返る。
  dt/4 では 3 回で二次収束する。bulge 0 の同じ段では同じ現象が弱く出ていた（0.027 mm、line search 5 回）。
- **対処（Abaqus と同じ自動増分）：** `emulate_abaqus.march`。dt は T/60 から始め、それが上限。増分の Newton が失敗
  （残差が非有限、30 回、要素の PNEWDT < 1）したら同じ出発状態から dt/4 でやり直す（1 増分に 5 回まで）。2 増分続けて
  4 回以内で収束したら dt を 1.5 倍（上限 T/60）。最後の増分は T で終わる。失敗がなければ従来の固定 60 増分と同じで、
  収束済みの結果は変わらない（2×6×8 の bulge 0 と 1 で要素ごとの値が同一）。json に `incrementation`（増分数、
  カットバック数、反復数、最小 dt、増分の記録）が入る。
- **デッキも同じ設定になった**（`make_implant_uel_inp.py` が書く）：GROW 段が
  `*COUPLED TEMPERATURE-DISPLACEMENT, DELTMX=10`（φ と c は [0, 1] なので届かない。自動増分に切り替えるためだけ）、
  データ行 `dt0 = T/60, T, DTMIN = dt0/4^5, DTMAX = dt0`、`*CONTROLS, PARAMETERS=TIME INCREMENTATION` に
  `20, 30, , 30, , , , 5`（反復 30 回、カットバック 5 回）。縮小率 0.25 と拡大率 1.5 は Abaqus の既定値。
  夕方の項の `NINC` を 120 にする案と `*CONTROLS, PARAMETERS=LINE SEARCH` を足す案は**不要になった**（デッキを作り直せば
  自動増分が入る。古い `.inp` は使わない）。
- **結果（固定シード、3 段）：** コア（0.25 ≤ z ≤ 1.75、15° ≤ θ ≤ 75°）の α − 1 は bulge 0 で 0.01795 / 0.01467 /
  0.01432（次数 2.9）、bulge 1 で 0.01490 / 0.01206 / 0.01166（次数 2.5）。コアの vM/E は bulge 0 で 0.0354 / 0.0364 /
  0.0361、bulge 1 で 0.0325 / 0.0321 / 0.0311（段ごとに 1〜3 %、次数は読めない）。最大 vM/E は上面の角（固定面の隣、
  シードと空の層の境目）で段ごとに上がり続ける（0.27 → 0.78 → 1.22、bulge 1 は 0.29 → 0.98 → 1.58）。角の特異点なので
  報告はするが比べる量から外す。bulge 1 はコアの α − 1 を 17〜19 %、vM/E を 8〜14 % 下げる。
  6×18×24 の bulge 1 は 64 増分、カットバック 1 回、36 分（クラウド、4 コアに 4 本）。

**fifa でやること（夕方の項の 2〜4 を置き換える）：**
1. `git pull origin master`、ライブラリを作り直す。`results/` に json がある段は飛ばす（2×6×8、4×12×16、6×18×24 の
   固定シードは bulge 0 と 1 とも済み）。
2. 残りは 8×24×32 の固定シード（`--seed-layers 4`）の bulge 0 と 1。1 本ずつ nohup で（6144 要素。クラウドの目安では
   2〜3 時間、カットバックが入ればもう少し）：
   ```
   cd p3/uel && python3 make_implant_uel_inp.py /tmp/imp && export BIOFILM_ECO_CASE=/tmp/imp/uel_eco.txt
   python3 make_implant_uel_inp.py --run --nr 8 --nt 24 --nz 32 --seed-layers 4
   python3 make_implant_uel_inp.py --run --nr 8 --nt 24 --nz 32 --seed-layers 4 --bulge 1
   python3 implant_convergence.py --fixed-seed; python3 implant_convergence.py --fixed-seed --bulge 1
   ```
   `5 cutbacks did not converge` で止まったら、その増分番号と stderr の最後の 30 行を報告に書いて止める（無理に通さない）。
3. Abaqus で 1 本（6×18×24、固定シード、**bulge 1**。カットバックが Abaqus でも同じに働くかを見る）：
   ```
   python3 make_implant_uel_inp.py abaqus_in --nr 6 --nt 18 --nz 24 --seed-layers 3 --bulge 1
   python3 make_abaqus_source.py
   export BIOFILM_ECO_CASE=$PWD/abaqus_in/uel_eco.txt BIOFILM_UEL_GP=$PWD/abaqus_in/uel_implant6x18x24_b1_sl3_gp.txt
   (cd abaqus_in && abaqus job=uel_implant6x18x24_b1_sl3 user=../unified_uel_abaqus.for cpus=1 interactive)
   python3 compare_implant_gp.py abaqus_in/uel_implant6x18x24_b1_sl3_gp.txt results/implant_6x18x24_b1_sl3.json
   ```
   期待：`.sta` の最初の増分がカットバック（attempt 2、dt 9.9e−4）のあと収束し、以後は 1.5 倍ずつ戻って 3.96e−3 に
   戻る。増分の刻みが模擬実行と同じなら Gauss 点の値は 1e−6 以下で合うはず（PASS）。刻みが違っていれば（Abaqus が
   別の増分で切った、など）差は時間離散化の分（1e−3 程度）になるので、`.sta` の増分列を報告に貼る。
4. json（5 MB 以下）、`.sta`、表を `keio/uel-implant` に push し、`KEIO_CHAIN_LOG.md` に報告（4 段のコアの値と次数、
   最大値が角に留まるか、Abaqus の増分列）。

## 2026年10月16日 夕方（扇形のメッシュ収束：クラウドで出た分と、fifa に残す分。急がない）

クラウド（gfortran + Python の模擬実行）で 2×6×8、4×12×16、6×18×24 を回した（`p3/uel/README.md` の
「Implant sector: mesh convergence」、`p3/uel/results/implant_*.json`。表は `python3 p3/uel/implant_convergence.py
[--bulge 1] [--fixed-seed]`）。分かったこと：

- **内側 1 層のシードはメッシュと一緒に薄くなる**（0.125 → 0.0625 → 0.042 mm）ので問題そのものが変わり、収束しない
  （シードの vM/E が −51 %、−37 %）。収束を見るには `--seed-layers f`（f = nr/2。シードの厚さを 0.125 mm で固定、
  ファイル名に `_sl<f>`）を使う。朝の依頼の 1〜2 はこの固定シードで読み替える。
- 固定シードでは、最大 vM/E が上面と固定面の角（シードと空の層の境目）で、細かくするほど上がる（0.27 → 0.78）。角の
  特異点なので、比べる量は**シードのコア**（0.25 ≤ z ≤ 1.75、15° ≤ θ ≤ 75° のシード要素。`implant_convergence.py` が
  json の要素データから計算する）にする。コアは 2 段で α − 1 が −18 %、vM/E が +2.6 %。
- 6×18×24 の固定シードでは模擬実行の Newton が最初の増分で発散した（力学だけ。要素は PNEWDT = 1 のまま）。
  `emulate_abaqus.py` の `increment` に、残差が増えたら刻みを半分にする line search を入れた（Abaqus なら
  `*CONTROLS, PARAMETERS=LINE SEARCH`）。収束した解は同じ（`tests/test_uel_implant.py`）。bulge 1 の 6×18×24 は
  line search を入れても最初の増分で 30 回で収束しない（未解決）。

**fifa でやること：**
1. `git pull origin master`、ライブラリを作り直す（`p3/uel/driver.py` の 3 行）。
2. クラウドで終わらなかった段を回す（1 本ずつ nohup で。クラウドの目安：6×18×24 は 25〜40 分、8×24×32 は 2〜3 時間。
   すでに `results/` にある json は飛ばす）：
   ```
   cd p3/uel && python3 make_implant_uel_inp.py /tmp/imp && export BIOFILM_ECO_CASE=/tmp/imp/uel_eco.txt
   for f in 3 4; do python3 make_implant_uel_inp.py --run --nr $((2*f)) --nt $((6*f)) --nz $((8*f)) --seed-layers $f; done
   python3 implant_convergence.py --fixed-seed
   ```
   bulge 1 の固定シード（`--bulge 1 --seed-layers 3`）は、最初の増分で止まったら `make_implant_uel_inp.py` の
   `NINC` を 120 にして（増分を半分に）試し、それでも止まれば報告だけ。
3. 細かい段を 1 本、Abaqus でも（6×18×24、固定シード、bulge 0）：
   ```
   python3 make_implant_uel_inp.py abaqus_in --nr 6 --nt 18 --nz 24 --seed-layers 3    # uel_implant6x18x24_sl3.inp
   python3 make_abaqus_source.py
   export BIOFILM_ECO_CASE=$PWD/abaqus_in/uel_eco.txt BIOFILM_UEL_GP=$PWD/abaqus_in/uel_implant6x18x24_sl3_gp.txt
   (cd abaqus_in && abaqus job=uel_implant6x18x24_sl3 user=../unified_uel_abaqus.for cpus=1 interactive)
   python3 compare_implant_gp.py abaqus_in/uel_implant6x18x24_sl3_gp.txt results/implant_6x18x24_sl3.json
   ```
   期待：2×6×8 と同じく α − 1、vM/E とも 1e−6 以下（PASS）。Abaqus の Newton が最初の増分で止まったら、GROW 段に
   `*CONTROLS, PARAMETERS=LINE SEARCH` を足して再実行し、そのことを報告に書く。
4. json と表を `keio/uel-implant` に push し、`KEIO_CHAIN_LOG.md` に報告（コアの α − 1 と vM/E の段ごとの値と次数、
   最大値の位置が角に留まるか）。

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

## 29 Oct 深夜：式 35 は ċ 付きが基本

- **これからの Python の本は全部 `--transient-c 1e8` を付ける**（後退オイラー、計算時間は同じ）。
- **成長則 A/B の系列（4 例 × 32³/64³/128³、`--growth-living` と `--growth-living-energy`）と 64³ の基準 2 本も `--transient-c 1e8` を付けて回し直す。**
  今走っている極限の系列は止めずに最後まで（ċ との照合に使う）。空いている CPU とメモリがあればすぐ並べて、なければ続けて。
  論文の表は ċ 付きの方を使う。JSON のタグは `_transc1e+08` が付く。報告では同じ格子の ċ と極限の差（生きている量、c、α−1、von Mises）を並べる。
- Abaqus：`abaqus_living/make_living_inp.py` は既定で `--transient-c 1e8`（PROPS(25)）。`--transient-c 0` で前と同じ極限。
  `living_uel_abaqus.for` を作り直した（接線テスト 1e−10）。Abaqus で走らせるなら 16³ で ċ あり・なしの 1 組を先に。

## 30 Oct 朝：混合の例の応力を 256³ で（論文用、急がないが成長則の系列の次）

背景（クラウド、41d3b65 の後の診断）：混合の例・抗菌薬ありの von Mises は 32/64/128/192³ で 2.00/2.82/3.30/3.52e−4 Pa（次数 0.43、GCI 41 %）。
前線の項を切ると 32³→64³ で 6.07→6.06e−4（収束）、ε² = 1e−8 は 1e−6 と同じ。原因は前線の項の縁の解像度：前線の速さ v = r c/(k+c) ≈ 4.9 mm/T*、
β = 0.02 mm²/T* で縁の厚さ β/v ≈ 0.004 mm、セルのペクレ数 v h/(2β) は 192³ で 1.28、**256³ で 0.96（初めて 1 を下回る）**。

**回すもの**（成長則 A/B の ċ 付きの系列が終わってから、1 本ずつ）：
1. `--nf 256 --front 1 --dtc 0.0025 --phi-min 0.001 --v2 --fortran --ansys-bc --inherit-living --alpha0 100 --gab 1 --gc 1 --tab 0.333333
   --tvd --wall-grad mirror --fast --transient-c 1e8 --mix-equal --only mixed_ab`
2. 同じで `--mix-equal` を外す（Fritsch の割合）
- 先に 1 本目の最初の 10 ステップでメモリと 1 ステップの時間を測り、`KEIO_CHAIN_LOG.md` に。192³ の約 3.2 倍の時間（IKMHIWI03 で 192³ は 7〜9 時間）と見込む。
- 報告：種の平均の α−1、vM、Σφ̄、c（`KEIO_CHAIN_LOG.md`、JSON は `keio/*` に push）。

## 30 Oct：成長則の系列を受けて（クラウド）

- 取り込んで master に merge した（c84e9be）。ありがとう。
- **足りない基準は回さなくてよい**：
  - 極限の 128³ の基準は修論の表 `tab:ch5-gci` の値（分離・抗菌薬あり α−1 9.112e−4、vM 6.284e−4；混合（等量）・抗菌薬あり vM 3.30e−4）を使う。
    128³ の vM の比は分離・抗菌薬あり 0.83、混合・抗菌薬あり 0.70（成長則 A、B とも）。
  - ċ 付きの抗菌薬なしの基準は、極限の基準との差が 6e−4 以下なので極限の基準で代える。
- vM の GCI が大きいのは、32〜128³ ではセルのペクレ数が 1 より大きいため（上の「30 Oct 朝」）。**次は 256³ の混合の 2 本**（上の項）。
