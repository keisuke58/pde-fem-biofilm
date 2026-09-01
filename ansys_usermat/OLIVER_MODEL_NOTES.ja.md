# Oliver から受け取った ANSYS モデル — 中身の解析

[English](OLIVER_MODEL_NOTES.md) | **日本語**

2026-09-01 に Oliver Höchel から受け取った2つの配布物を解析した記録:
Workbench プロジェクト `BiofilmImplementation.wbpz` と、UPF ソース一式
`Nishioka_Hoechel.zip`。**どちらもコミットしていません** — 合わせて約26 MB の
バイナリと、他グループのソースであり、再配布する立場にないためです。ファイルを
開き直さずに統合方針を議論できるよう、内容をここに記録します。

記載はすべてファイルから直接読み取ったものです（デッキ `ds.dat` / `solve.out`
と Fortran ソース）。コードが明言していない解釈には、その旨を明記しています。

再現用のノートブック → [`oliver_model_analysis.ipynb`](oliver_model_analysis.ipynb)
（下記の数値・表をファイルから再導出します）。

---

## まとめ（先に結論）

| 層 | Oliver 側 | 状態 |
|---|---|---|
| パラメータ受け渡し | `USolBeg` の `parevl` → 共通ブロック | ✅ |
| 微分演算子 | NEM（重み付き最小二乗） | ✅ 検証機構つき |
| 場の求解 | `USSFin` で PARDISO ×3 | ✅ |
| 生態モデル | **2菌種 × 2栄養**、Monod + 相互作用 | ✅ 論文のスキーム |
| 材料（弾性） | `AceGenNeoHookV04`（バイオフィルム/空隙ブレンド） | ✅ |
| 材料（粘性） | ガラス用のみ（コメントアウト） | ❌ バイオフィルム用は無し |
| **成長 `Fg=(1+α)I`** | — | ❌ **こちらの貢献部分** |

**補完関係です。** 彼らは場の求解機構と n=2 のバイオフィルム弾性則を持ち、
こちらは較正済み n菌種生態モデル・検証済み粘性則・成長運動学を持っています。

---

## 1. 何が届いたか

### `.wbpz`（Workbench プロジェクト）

| 項目 | 値 |
|---|---|
| ANSYS | **2024 R2**（build 24.2） |
| 要素 | `SOLID185`、8ガウス点、`NLGEOM,ON` |
| メッシュ | 18,750 要素 / 23,556 節点 |
| 形状 | `PRJ11_TestCube` — **テスト立方体**（歯・インプラントではない） |
| 状態変数 | `TB,STATE,1,,100` |
| プロパティ | `TB,USER,1,1,1` — 1個、かつ `TBDATA` 無し |
| 求解 | 11 サブステップ、各2反復で収束、**エラー0** |

材料パラメータは `TBDATA` ではなく **APDL パラメータ（`*SET`）** で渡されます。
だから `TB,USER` のプロパティが1個で足りています。

### `Nishioka_Hoechel.zip`（UPF ソース一式）

`ANSYS-Pool/` に Fortran ソース、オブジェクト、ビルド済み `libansuser.so`、
ビルドスクリプトが入っています。構成則ルーチンは **AceGen 生成**
（Mathematica の記号処理 → Fortran、`sms.h` が AceGen ランタイムヘッダ）。

---

## 2. ⚠️ `usermat` の引数リストがリリース依存で、変わっている

両ソースから直接数えた結果:

| | 引数数 | `cutFactor` 以降 |
|---|---|---|
| **2024 R2**（Oliver） | **41** | `pVolDer, hrmflg, var3, var4, var5, var6, var7` |
| **v222**（このリポジトリ） | **42** | `var1, var2, var3, var4, var5, var6, var7, var8` |

2024 R2 では予約枠 `var1`/`var2` が名前付き引数 `pVolDer(3)`（体積ポテンシャルの
J による1〜3階微分）と `hrmflg`（調和解析フラグ）になり、**`var8` は削除**。

したがって `usermat_biofilm.f` は**そのままでは 2024 R2 でビルドできません** —
引数を1個多く宣言しており、`var8` が実引数リストの外を読みます。`README.md` が
「引数リストはリリース依存」と警告していた通りのことが、具体的に起きています。

**ビルド方式も違います**: `ANSUSERSHARED`（Linux・共有ライブラリ・クラスタ、
`ifx`/`icc`）であり、こちらの `ANSCUST.BAT`（Windows・カスタム `ANSYS.exe`）
ではありません。

→ 実機で試す手順は
[`apdl/V222_PORT_INSTRUCTIONS.md`](apdl/V222_PORT_INSTRUCTIONS.md)

---

## 3. 全体の処理の流れ

```
ANSYS の求解ループ（SOLID185, NLGEOM,ON）
 │
 ├─ USolBeg ......... 解析開始時に1回
 │     · parevl 約150回 = APDL パラメータを /usercm/ 共通ブロックへ
 │       （↑ これが prop() を使わない理由）
 │     · InitVals            — データアリーナ確保
 │     · NEM_CreateData_Init — メッシュフリー近傍演算子を構築
 │
 ├─ usermat ......... ガウス点ごと・平衡反復ごと
 │     · GetVals / GetTMP  — プールから当該点の状態を読む
 │     · CALL AceGenNeoHookV04 → 応力, dsdePl   ← 我々の則が入る場所
 │     · SetVals          — 書き戻す
 │
 └─ USSFin .......... サブステップごと
       · 組み立て + PARDISO → 温度場
       · 組み立て + PARDISO → Bio/Nut 場（2系統）
       · CalcLaserIntegralOMP / CALCPYRO（レーザー・パイロメータ＝ガラス用）
```

**交互解法（staggered / operator-split）**です。ANSYS が力学を、`USSFin` が
輸送場を PARDISO（Intel MKL 疎行列直接解法）で解き、サブステップごとに交互に
進みます。**輸送場は ANSYS の自由度ではなく**、UPF 自前のデータプールにのみ
存在します（`TB,STATE` が100個必要な理由）。

`Usermat_*.F` は `mpif.h` を include し、`GetVals`/`SetVals`/`SetNEM` という
**MPI 共有データプール**の API を宣言しています（フォルダ名 "ANSYS-Pool" の由来）。

---

## 4. NEM の正体 — 重み付き最小二乗の微分演算子

名前に反して Natural Element Method ではなく、散在点上の **移動最小二乗（MLS）
微分演算子**、すなわち一般化有限差分です。

`CalcDMat` が各点で、近傍約30点（`NEIGHBOR_CNT`）にわたる `mD(9, 近傍数)` を作ります:

```fortran
! ガウス核（バンド幅 BETA_STAR）、WMAT_THRESHOLD で疎化
mW(ii,ii) = EXP(-0.5*((4.0*SQRT(dx**2+dy**2+dz**2))/BETA_STAR)**2)
if (mW(ii,ii) <= WMAT_THRESHOLD) mW(ii,ii) = 0.0
...
CALL InversGauss(mResult3, 9)      ! 9×9 モーメント行列を反転
```

9行の意味（`AssembleSparse` の使い方から確定）:

| 行 | 意味 | 根拠 |
|---|---|---|
| 1–3 | ∂xx, ∂yy, ∂zz | `mD(1)+mD(2)+mD(3)` が Laplacian として組まれる |
| 4–6 | 混合2階微分 | `ASSEMBLE_KEY=1` では未使用 |
| 7–9 | ∂x, ∂y, ∂z | `Vdp_Dx_T` 等へ |

組み立て時に**自己係数 = 近傍係数の総和の符号反転**としており、演算子が定数を
消す（0次整合性）ことが保証されます。`DEBUG_KEY` の解析的テスト関数群
（`φ=(x²+y²+z²)/6`、cos、sin、exp、境界勾配チェック）が検証手段です。

**JAXFEM との対比**: こちらは構造格子の有限差分で微分を取ります。あちらは散在点への
最小二乗フィット。同じ微分作用素を別経路で得ており、彼らの方式だと場を FE 積分点に
直接載せられます。

---

## 5. 生態モデルは論文のもの、ただし n = 2

`usercm.inc` と `Ussfin` の場の更新を合わせると、モデルが読み取れます。
**2菌種 × 2栄養**です（当初「1バイオフィルム + 2栄養」と誤読していました）:

```
sGdp_Bio1start,  sGdp_Bio2start        ! 2菌種
sGdp_Nut1start,  sGdp_Nut2start        ! 2栄養
sGdp_MaxGrowth11 .. 22                 ! 2×2 最大増殖速度（菌種 × 栄養）
sGdp_HalfVelo11  .. 22                 ! 2×2 半飽和定数
sGdp_Interaction12, sGdp_Interaction21 ! 菌種間相互作用
```

増殖項（`Ussfin` 約1900行目）:

```fortran
GrowthBio1 = SQRT(Sdp_LapBio1**2) *
     &(  ( (sGdp_MaxGrowth11 + sGdp_Interaction12 * vGdp_Bio2_n(ID))
     &      * vGdp_Nut1_n(ID) ) / (sGdp_HalfVelo11 + vGdp_Nut1_n(ID))
     &  + ( (sGdp_MaxGrowth21 + sGdp_Interaction12 * vGdp_Bio2_n(ID))
     &      * vGdp_Nut2_n(ID) ) / (sGdp_HalfVelo21 + vGdp_Nut2_n(ID)) )
```

項ごとに読むと:

- **Monod 動力学**（栄養ごと）— `μ_max·S/(K_s + S)`、`HalfVelo` が `K_s`。
- **相互作用が増殖速度そのものをずらす** — `MaxGrowth + Interaction12·Bio2`。
  別項を足すのではなく「どれだけ速く増えるか」を相手菌種が変える形。これが
  published paper の *"novel interaction scheme"* です
  （Klempt, Geisler, Soleimani et al., *Archive of Applied Mechanics* **96**,
  164 (2026), doi:10.1007/s00419-026-03160-y。`biofilm_3tooth_refs.bib` に
  `Klempt2026ContinuumBacterialGrowth` として登録済み）。
- **`|∇²Bio|` を前係数に** — 場が曲がっている場所＝フロントで増殖が起きる。
  `JAXFEM/` の Allen–Cahn 界面項と同じ役割。
- **走化性的な配向** — `OriBio = Σ_j OriWeight_j · NormDot(∇Nut_j, ∇Bio)`。
- **ペナルティによる [0,1] 拘束** —
  `-Penalty·( max(0, Bio-1) + min(0, Bio) )`。3つのコードで3通り:
  ここはペナルティ、`hamilton_ode_jax.py` は対数バリア、ガラスモデルは
  ロジスティックシグモイド。

バイオフィルムの更新には `!Biofilm / lokales Biofilm Update (explizit, TEST)`
とあり、陽的時間積分・暫定扱いです。

> **統合にとって最も重要な発見。** 彼らの枠組みは「生物学を教え込む必要のある
> 汎用 PDF ソルバ」ではなく、**既に論文の生態モデルを実装済み**です。ただし
> **n = 2**。こちらは **n = 5**（`hamilton_ode_jax.py`）と**一般 n**
> （`hamilton_ode_jax_nsp.py`）。`Interaction12`/`Interaction21` は、こちらが
> TMCMC で較正した行列 `A` の非対角成分に対応します。
>
> つまり埋めるべき差は「モデルを移植する」より狭く、**2 → 5 菌種への一般化**、
> **較正済み相互作用行列の供給**、**成長運動学の追加**の3点です。

---

## 6. 材料 — 動いているのはバイオフィルム側

有効な構成則の呼び出しは1つだけです:

```fortran
CALL AceGenNeoHookV04(Vdp_AceGen, defGrad, stress, dsdePl,
     &   sGdp_YoungBio,   sGdp_YoungVoid,
     &   sGdp_PoissonBio, sGdp_PoissonVoid,
     &   Sdp_sumBio, Sdp_sumLocal, sedEl, ID)
```

バイオフィルム剛性と空隙剛性を局所バイオフィルム量でブレンドする Neo-Hookean
（`L` 接尾辞は *leer* ＝ 空）。**仮置きではなくバイオフィルム専用材料**で、
2026-08-05 生成とプール内で最新です（V02: 1/8 → V03: 7/29 → V04: 8/5、
署名は3つとも同一なので、インターフェイスは1月に確定して中身だけ調整中）。

コメントアウトされている方が**ガラス**です:

```fortran
!------Matmodell Tobi Start
!      IF(Sdp_Phi .EQ. 0.0D0)THEN
!      !Air Phase
!      CALL AGStressP21V07(... Sdp_T_n, vGdp_Th_Expans ...)
!------Matmodell Tobi End
```

### ⚠️ この `sAlpha` は成長 α ではありません

名前が完全に誤解を招くので明記します:

```fortran
Sdp_sumBio   = Sdp_bio1_n + Sdp_bio1_n              -> sBiofilm
Sdp_sumLocal = (Sdp_locbio1_n + Sdp_locbio2_n)/2    -> sAlpha
```

宣言のコメントは `!Summe biofilm/local Biofilm`。`sAlpha` は**局所バイオフィルム
平均**であって成長変数ではありません。**プール内のどこにも成長運動学は存在せず**、
`Fg = (1+α)I` は依然としてこちらが持ち込む部分です。

### バグの可能性

`Sdp_sumBio = Sdp_bio1_n + Sdp_bio1_n` は `bio1` を自分自身に足しています。
次の行が `locbio1` と `locbio2` を正しく平均しているので、`bio1 + bio2` の
書き間違いに見えます。こちらではビルドできないため断定はせず、Oliver への
質問として扱っています。

---

## 7. Mathematica ノートブックの位置づけ

`BiofilmTSMMathematica20250117_4species__changing.nb`（Mathematica 12.2、
2025-01-17）は、Fortran ルーチンの **AceGen ソースではありません**（`SMS*`
呼び出しが皆無）。同じ増殖モデルの **Mathematica 直接実装**で、構成は:

```
TSM-Growth
  ├── Aufstellen der Gleichungen   （方程式の構築）
  ├── Speichern der Variablen
  ├── Newton-Raphson               （求解）
  ├── Ausgabe
  └── Speichern der Ergebnisse
```

パラメータ語彙がこのリポジトリと共通です — `Kp1`（26回）と `Eta`（96回）は
どちらも `JAXFEM/hamilton_ode_jax.py` のパラメータ名で、`Kp1` は対数バリア定数。
**4菌種**版です。

系統としては: **論文 → この Mathematica ノートブック（n=4）→ Oliver の Fortran
（n=2、ANSYS 内）**、そしてこのリポジトリの `hamilton_ode_jax.py`（n=5）と
`hamilton_ode_jax_nsp.py`（一般 n）が同じ木の JAX 側の枝。したがって
**「方程式が何か」を知るにはノートブック、「ANSYS でどう解くか」を知るには
Oliver の Fortran** が良い参照先です。

---

## 8. `usermat_biofilm.f` との関係

| | Oliver 側 | このリポジトリ |
|---|---|---|
| スコープ | **非局所** — 場を自分で解く | **局所**構成則のみ |
| 成長ドライバ | 内部で解く | `ustatev(10) = α` として**外から与える** |
| 状態変数 | 100 | 10（`kStateMat=1` で14） |
| プロパティ | APDL パラメータ経由（`TB,USER` は1個） | `TB,USER` に6〜7個 |
| 微分 | メッシュフリー近傍演算子 | 不要（点局所） |
| ANSYS | 2024 R2 / Linux | v222 / Windows で検証 |

**統合の2案**（指導教員案件であって、コーディングの判断ではありません）:

- **(A) 構成則を彼らの枠組みへ移植する。** 場は彼らの NEM が供給し、こちらは
  検証済みの `Fg=(1+α)I` 成長 + Mooney-Rivlin/`D1` + 粘性応答を提供。
  `crosscheck/` の 0 ULP 一致がそのまま付いてきます。移植先は
  `CALL AceGenNeoHookV04(...)` の位置で、`Vdp_Cv_n` が粘性状態の枠として
  既に確保済みです。この場合 `usermat_biofilm.f` は成果物ではなく、移植先を
  照合する参照実装になります。
- **(B) 分離したままにする。** こちらは `JAXFEM/` から α を受け取る局所則、
  彼らは自前で場を解く。答える問いが別になるので、直接比較はできません。

---

## 9. Oliver への質問（優先順）

1. ~~**USERMAT の Fortran ソース**~~ — **受領済み**。
2. **どのリリース／プラットフォームに合わせるか。** 彼らは Linux / 2024 R2
   （クラスタ、`ANSUSERSHARED`）、こちらは Windows / v222（`ANSCUST.BAT`）。
   署名が違う（41 vs 42）ので、1つのソースでは両立しません。クラスタの
   アカウントをもらえるかが現実的な分岐点です。
3. **場を計算するのはどちらか。** 彼らの NEM が解くなら、こちらの α 場写像
   （`ustatev(10)`）は冗長になり、(A) が実際の統合経路になります。ただし
   これは彼らの判断であって、こちらが勝手に前提にすべきではありません。
4. **バイオフィルム材料はどこまで行く想定か。** `AceGenNeoHookV04` は
   バイオフィルム専用だが**純弾性**。粘性則は予定にあるか（それがこちらの
   持ち込める部分）。関連して `Sdp_bio1_n + Sdp_bio1_n` は誤記か。
5. **AceGen の Mathematica ノートブックはもらえるか。** 構成則は機械生成なので、
   生成後の `.f` を手で直しても次の生成で消えます。`Fg` を入れるならノートブックが
   編集場所になります。

---

*両配布物とも 2026-09-01 に解析。元ファイル内の内部パスやクラスタのユーザ名は
意図的に再掲していません — このリポジトリは公開であり、ソースは他グループの
コードです。*
