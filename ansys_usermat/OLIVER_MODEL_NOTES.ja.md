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
| 生態モデル | **2菌種 × 2栄養**、Monod + 相互作用 | ⚠️ 論文の式ではない(§5の訂正参照) |
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

## 5. ⚠️ 訂正(2026-09-08): 生態モデルは「論文のもの」ではなかった — 構造が似ているだけ

以下(2026-09-01 時点の記述)は「`sGdp_Interaction12`/`21` が Klempt/Geisler/
Soleimani/Junker の *A continuum multi-species biofilm model with a novel
interaction scheme*(arXiv:2509.01274、出版版 AAM 96, 164 (2026)、
doi:10.1007/s00419-026-03160-y)の"novel interaction scheme"そのものだ」と
主張していますが、**2026-09-08 に arXiv PDF を直接読んで確認したところ、この
主張は成立しません**。論文の実際の式(Eq.10 エネルギー密度、Eq.16–18 強形式):

```
Ψ  = -½ c*  φ̄·A·φ̄  +  ½ α* ψ·B·ψ                     (Eq. 10)
Ia_i := a_ii·φ̄_i + Σ_j a_ij·φ̄_j   = (A·φ̄)_i             (Eq. 16/17 由来)
Eq. 16 (φ_i):  0 = -c*·ψ_i·Ia_i + η_i(φ̇_iψ_i²+φ̄_iψ̇_i+φ̇_i) + γ
Eq. 17 (ψ_i):  0 = -c*·φ_i·Ia_i + α*·ψ_i·b_i + η_i(ψ̇_iφ_i²+φ̄_iφ̇_i) + γ
```

**A は明示的に対称行列**("the symmetric growth coefficient matrix A")で
**対角成分(自己相互作用)も含み**、c\* は栄養を表す**単一スカラー**で
**線形**に効きます。つまり論文の相互作用は**加法的・双線形**(`Ia = A·φ̄` が
散逸由来の残差に足し込まれる形)。

対して Oliver の `Interaction12·Bio2` は、Monod飽和曲線の最大増殖率を
**乗法的**にずらす項(`(MaxGrowth + Interaction12·Bio2)·Nut/(HalfVelo+Nut)`、
栄養場2本、それぞれ独立にMonod飽和)で、**対角自己項がなく**、
`Interaction12`と`Interaction21`が別々の定数として立てられていて対称性
(`A_ij=A_ji`)も要求していません。これは記法違いの同じ式ではなく、**構造的に
別の式**です。Oliver の増殖項の構造(Monod動力学 × `|∇²Bio|`界面局在化 ×
走化性的配向)はむしろ、もう一つ別のKlempt論文 — Klempt, Soleimani,
Wriggers, Junker, *A Hamilton principle-based model for diffusion-driven
biofilm growth*, Biomech Model Mechanobiol 23, 2091–2113 (2024)
(`Klempt2024DiffusionDrivenGrowth`、`CITATION_AUDIT.md` §F1b 参照) —
の記述(このリポジトリの`JAXFEM/felix_complete_reproduction.py`が
"Eq.34–36"として再現しているAllen-Cahn+logistic-Monod+走化性)に近いです。
`Interaction12/21`自体がどちらの論文でも式番号付きで確認できたわけではなく、
Oliver/Felixが2024論文の増殖構造に独自に足した項である可能性があります。

**実務上の意味:**

- **このリポジトリ自身のn=5生態モデル**(`hamilton_ode_jax.py` /
  `jax_hamilton_0d_5species_demo.py` / `ecology_jax.py` — TMCMCが較正し、
  `usermat_biofilm.f`の`kUseEcology`ブリッジが実際に走らせているもの)は、
  **Eq.10/16–18と項レベルでほぼ完全一致**することを確認しました。これが
  学術的に正しい経路であり、既に実ANSYSで動作確認済み
  (`coupling/README.md`の2026-09-07状態記録)— Oliverの`Ussfin`ではありません。
- 2026-09-07にOliverの`InteractionIJ`をn=5へ一般化して見つけたしきい値
  (0.0005で安定・0.0007で発散、2026-09-08に実ANSYSで再現確認済み)は
  **Oliver独自の(論文の式ではない)増殖則拡張の、本物で再現性のある性質**
  ですが、これは「論文モデルの臨界結合強度」ではなく「このANSYS/NEM
  パイプライン独自の拡張の安定限界」と呼ぶのが正確です。詳細は
  `apdl/N3_GROWTH_TRIAL.md`の対応する訂正を参照。
- TMCMCで較正した行列`A`をそのまま`InteractionIJ`に流用するのは
  **カテゴリーエラー**です(較正対象の式が違うため)。

---

### `usercm.inc` と `Ussfin` の場の更新(以下は元の記述、上記訂正を踏まえて読んでください)

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
  別項を足すのではなく「どれだけ速く増えるか」を相手菌種が変える形。
  **これは published paper の *"novel interaction scheme"* ではありません**
  （Klempt, Geisler, Soleimani et al., *Archive of Applied Mechanics* **96**,
  164 (2026), doi:10.1007/s00419-026-03160-y — 上記§5の訂正を参照）。
  論文の実際の式は加法的・双線形(`Ia = A·φ̄`、対称行列、対角自己項あり)で、
  乗法的なMonodシフトではありません。一致するのはこちらの
  `hamilton_ode_jax.py`/`ecology_jax.py`の方です。
- **`|∇²Bio|` を前係数に** — 場が曲がっている場所＝フロントで増殖が起きる。
  `JAXFEM/` の Allen–Cahn 界面項と同じ役割。
- **走化性的な配向** — `OriBio = Σ_j OriWeight_j · NormDot(∇Nut_j, ∇Bio)`。
- **ペナルティによる [0,1] 拘束** —
  `-Penalty·( max(0, Bio-1) + min(0, Bio) )`。3つのコードで3通り:
  ここはペナルティ、`hamilton_ode_jax.py` は対数バリア、ガラスモデルは
  ロジスティックシグモイド。

バイオフィルムの更新には `!Biofilm / lokales Biofilm Update (explizit, TEST)`
とあり、陽的時間積分・暫定扱いです。

> **2026-09-08訂正により失効。** 元々この段落は「既に論文の生態モデルを
> 実装済み」「`Interaction12`/`Interaction21`はTMCMC較正済み行列`A`の非対角
> 成分に対応する」と主張していましたが、どちらも成立しません。Oliverの増殖
> 構造はMonod/栄養駆動で対角自己項も対称性の要求もなく、こちらの`A`が較正
> されている論文の`Ia = A·φ̄`式とは構造的に別物です。**したがってTMCMC較正済み
> の`A`をそのまま`InteractionIJ`へ流用するのはカテゴリーエラーであり、単純な
> 差し替えではありません。** 論文に忠実な生態モデルをANSYS内で走らせる必要が
> あるなら、既に検証済みの経路は`ecology_jax.py`自身の`Ia = A @ (phi*psi)`を
> material callへ配線すること(`coupling/README.md`2026-09-07の状態記録で
> 完了済み)であり、Oliverの別系統のMonod相互作用項を一般化することでは
> ありません。Oliver自身の項を「論文モデルのn=2版」としてではなく、
> 「意図的に別の増殖率ヒューリスティック」として維持する価値があるかどうかは
> 彼に確認すべき論点で、このリポジトリ側で決められることではありません。

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
