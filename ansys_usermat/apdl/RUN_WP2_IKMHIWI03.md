# WP2：恒常圧の成長則（IKMHIWI03、2026年10月7日）

依頼：`IKMHIWI03_TO_SUBMISSION.ja.md` §5。試作と結果：`keio_wp2/README.ja.md` 第4〜7段。

    α̇ = k_α φ max(0, 1 − p/p_h)，p = −tr σ / 3（圧縮が正）
    p_h = P_h E k_α T*（E 基準）または P_h E(φ² + f) k_α T*（局所剛性），T* = 1

この則は仮定で、論文から取ったものではない。

## 実装（既定は off、今の結果は変わらない）

- 入れた場所：共有の fragment `callsite/phi_mode_exec.inc` と `growth_from_phi.f`
  （`BIOFILM_HOMEOSTATIC_FACTOR`）。この fragment は、オリバーの要素の usermat
  （`paste_fragments.py`）にも Abaqus の UMAT（`make_umat.py`）にも同じものが入るので、
  ANSYS と Abaqus は同じ形になる。
- 依頼には `usermat_biofilm.f` と書かれていたが、そこで α が更新されるのは生態モード
  （kUseEcology）だけである。修論と慶應の実行で Eq. 36 を積分しているのは上の
  fragment なので、そちらに入れた。
- 定数（ANSYS の prop、Abaqus の constants、番号は同じ）：
  - 49 = p_ref = P_h E k_α T*（応力の単位）。0 なら off。
  - 50 = 0：p_h = p_ref（E 基準）、1：p_h = p_ref (φ² + 51)（局所剛性）。
  - 51 = f。
- p は、ホストが渡す増分の始めの応力から計算する（ANSYS の `stress`、Abaqus の
  `STRESS`。どちらも前のステップで収束した値）。α の増分に g = max(0, 1 − p/p_h) を掛ける。
  状態変数は 95 = p、96 = g。モード 3（オリバー自身の α）には掛けない。
- Abaqus の入力：`make_implant_inp.py --ph P_h [--ph-local]`。p_ref は P_h × 定数 43（E）×
  k_α で計算する。ほかに `--first-order`（消費 gφc、試作と同じ）と
  `--base free|uz|clamped`（下面）を足した。既定のままなら前の入力と同じになる。

## 確認

| 確認 | 状態 |
|---|---|
| gfortran で Abaqus の UMAT、拘束した1点（`tests/test_homeostatic_growth.py`）：off なら Eq. 36 のまま。on なら α が漸化式 α_{n+1} = α_n + kφΔt max(0, 1 − p(α_n)/p_h) と 1e−9 で一致し、閉じた形の p(α*) = p_h で止まる。p_h の2つの形のどちらでも | 済、pass |
| 全体のテスト（`run_tests.ps1`） | 済、all green |
| ANSYS、オリバーの要素：8³ の立方体、全体が φ = 1、全節点を固定、T* = 5（`hp_clamped_check.py`、exe `F:\biofilm_upf_hp`）。off / E 基準 / 局所剛性の3本 | WP3 の判定の直後に自動で回る（`_chain_hp_check.log`）。報告は `results/2026-10-wp2/hp_clamped.txt` に push |

## 実行（Abaqus、インプラントと歯）

ANSYS にはこの形（カラー部・歯冠）のデッキがないので、Abaqus で回す。比べる相手は
`keio_wp2/results_nutrient_bc.json`（試作の bc_sweep）。

共通の設定は試作に合わせた：2菌種の case 6（成長は Eq. 36、全体の φ）、消費 1（gφc）、
前線の項 r = 10、w = 0.5、β = 0.02、dt = 0.002（dt 0.01 では前線が1ステップで1要素以上進み、
発散した）、T* = 1、Fortran の点モデル、4 CPU。

| 形 | 下面 | 帰還 | ジョブ名 |
|---|---|---|---|
| インプラント（ri 2.05、栄養は上面） | free / uz / clamped | なし / P_h = 0.1（局所剛性） | `wp2_imp_<下面>_{none,ph01}` |
| 歯（ri 2.05、bulge 1、栄養は外面） | free / uz / clamped | なし / P_h = 0.1（局所剛性） | `wp2_tooth_<下面>_{none,ph01}` |

12本で、1本 7〜8 分、合わせて約 1.5 時間。スクリプトは `F:\abaqus_work\_wp2_1007.ps1`、
ログは `_wp2_1007.log`。結果は1本ごとに `abaqus_composition/results_1007/<job>.txt`
（`summarize_implant.py`）に書いて push する。ANSYS とは別のライセンスで、ANSYS の
チェーンと同時に回しても止めない。

比べる量：周方向応力の最大とその高さ、α − 1 の平均。試作の単位は Pa（E = 10 Pa）、
Abaqus の出力は MPa（E = 1e−5 MPa）なので、比べるときに 1e6 倍する。
