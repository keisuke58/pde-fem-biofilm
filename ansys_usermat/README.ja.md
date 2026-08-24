# UMAT／USERMAT 解説 ― バイオフィルム成長・粘弾性構成則（ANSYS 移植版）

[English](README.md) | **日本語**

`ansys_usermat/usermat_biofilm.f` は、検証済みの **Abaqus UMAT**
（`umat_biofilm_visco.f` / `umat_biofilm_visco_phase2.f`）を **ANSYS Mechanical APDL
の `USERMAT`** へ移植したものです。Felix／IKM が持つ既存の ANSYS FE モデルの
現象論的な材料則を置き換え、**同一の成長・粘弾性構成則** を各ガウス点で
呼び出せるようにします（提案する修論テーマの出発点）。

- 構成則：`F = Fe·Fv·Fg`、成長は等方 `Fg = (1+α)I`
- **構成則の代数（式）は Abaqus 版の一行一行の写し**：Neo-Hookean 偏差応力 +
  `D1` 圧力項、後退オイラーの粘性更新、F 摂動による整合接線。
- **ANSYS 固有なのは「インターフェイス」だけ**。物理・アルゴリズムは同一。

## Abaqus UMAT ↔ ANSYS USERMAT 対応表

| Abaqus | ANSYS `usermat` | 意味 |
|---|---|---|
| `DFGRD1` / `DFGRD0` | `defGrad` / `defGrad_t` | 3×3 変形勾配テンソル |
| `STRESS(NTENS)` | `stress(ncomp)` | Cauchy 応力 |
| `DDSDDE` | `dsdePl(ncomp,ncomp)` | 材料接線（ヤコビアン） |
| `STATEV` | `ustatev(nStatev)` | 状態変数 |
| `PROPS` | `prop(nProp)` | 材料プロパティ |
| `DTIME` | `dTime` | 時間増分 |
| `PNEWDT < 1` | `keycut = 1`（+ `cutFactor`） | 増分カットバック要求 |
| `SSE` / `SPD` | `sedEl` / `sedPl` | ひずみエネルギー／散逸 |
| 並び `11,22,33,12,13,23` | 並び **`11,22,33,12,23,13`** | ⚠️ せん断の 5↔6 が入れ替わり |

**最大の落とし穴は応力成分の並び順**。ここでは `VI/VJ` の Voigt マップ
（`data VI /1,2,3,1,2,1/`, `VJ /1,2,3,2,3,3/`）で吸収しています。

## プロパティと状態変数

```
prop(1)=C10  prop(2)=C01  prop(3)=D1  prop(4)=eta  prop(5)=mtype  prop(6)=kUsePy
prop(7)=kStateMat
ustatev(1:9)=Fv（行優先の 3×3）   ustatev(10)=alpha（成長ドライバ）
ustatev(11:14)=C10,C01,D1,eta     （積分点ごと、kStateMat=1 のときのみ使用）
```

- **成長ドライバ `alpha`** は JAXFEM の α 場を各積分点へ写像したもの
  （`TB,STATE` やユーザ場で初期化、あるいは時間発展）。
- `kUsePy=1` にすると、インライン Fortran 則の代わりに
  **Python マテリアルフック**（後述）を選択します。
- `kStateMat=1` にすると、材料定数を `prop(1:4)` ではなく
  **積分点ごとの `ustatev(11:14)`** から読みます（後述）。

## 組成依存剛性 E(φ)（`kStateMat=1`）

モデルの第2の脚（`RESEARCH_MODEL.md` §3）。剛性は成長場 α を*経由せず並列に*
効きます。定数を `prop(1:4)` に固定していると全ガウス点が同じ剛性になるため、
4条件を区別するのが α だけになり、**本研究で最大の力学的差 ―
E が約 995 Pa（健全）〜 32 Pa（病的）の約31倍 ―** が丸ごと抜け落ちていました。

組成は CLSM の**測定値**であって解析中に発展する量ではないので、これらの定数は
計算開始前に確定しています。[`coupling/composition_to_material.py`](coupling/composition_to_material.py)
が一度だけ計算し（φ → `material_models.py` の `E(φ)`/`DI` → `C10, C01, D1, eta`）、
初期状態として配る `TB,USER`／`TB,STATE` ブロックを出力します。
**増分ごとの Python 呼び出しは一切発生しません** ― この経路は高速なインライン
Fortran コアの内側で完結します。（ソケットブリッジ `kUsePy=1` は、係数ではなく
構成則**そのもの**を差し替えるという別の問題を解くものです。）

```bash
python ansys_usermat/coupling/composition_to_material.py --phi 0.2,0.2,0.2,0.2,0.2
python ansys_usermat/coupling/composition_to_material.py --E 32 --di 0.85 --apdl
```

`ustatev(11) <= 0` は「未初期化」と解釈して `prop(1:4)` にフォールバックします
（`Fv` に対して `INIT_FV_IF_ZERO` が既に使っている「ゼロ＝未設定」イディオムと同じ）。
設定ミス時に剛性ゼロで黙って走るのではなく prop の材料に縮退します。
なお ANSYS の `TB,STATE` は**材料ごと**に効くので、空間的に組成が変わる場合は
組成ビンごとに材料を分ける必要があります。

[`tests/test_composition_material.py`](../tests/test_composition_material.py)
で end-to-end 検証済み: 材料 A の定数を `prop` 経由で流した結果と、`prop` に
材料 B を置いたまま状態変数経由で A を流した結果が一致すること（＝状態側が
確実に上書きしていること）、および固定変形で4条件が
**CH / DH / CS / DS = 566 / 177 / 555 / 20 Pa**（max |σ|）になること ―
prop 定数モデルでは表現できない差です。

## Python マテリアルフック（ガウス点ごと）

修論の核心的な成果 ― 論文で較正した **Python** 材料モデルを各ガウス点で
呼び出す ― は、ソース中の `PYTHON MATERIAL HOOK` として**実装・end-to-end
検証済み**です。仕組みは `ISO_C_BINDING` ／ローカルソケット橋渡しで、
`(defGrad, Fv_old, alpha, dTime, prop)` を Python へ送り、
`(stress, Fv_new, dsdePl)` を受け取り、Abaqus→ANSYS の Voigt 並び替え
（`MAP6`）を経て `usermat()` 自身の `stress`/`ustatev`/`dsdePl` へ書き戻し
ます。`kUsePy=1` でこの経路が有効になり、Python サーバへ接続できない・
応答が不正な場合は求解を失敗させず検証済みのインラインコアへフォール
バックします（アーキ図：`ch5_flow/flow_python_material_hook`）。

ブリッジ本体は [`coupling/`](coupling/README.md) にあります：Python 側
（`material_server.py` ― NumPy コア＋F摂動接線＋ソケットサーバ）、通信
プロトコル、Fortran 側フック（`usermat_py_hook.f`）。バイパス用ドライバ
ではなく**実際の `usermat()` エントリポイント**を通して
`tests/test_usermat_kusepy_e2e.py` が検証しており、`usermat_biofilm.f` +
`usermat_py_hook.f` + `biofilm_py_eval.c` をスタンドアロンドライバへ
ビルドし、弾性／粘性／Mooney-Rivlin の各ケースで `kUsePy=1` と
`kUsePy=0` を突き合わせます ― 応力と更新後の粘性状態は数値精度で一致、
整合接線（`dsdePl`）は浮動小数点誤差レベルで一致します（両側とも同じ
F摂動方式、`PERT=1e-7`、を使うため）。サーバ未接続時のフォールバック
経路（`PYOK=.false.` → インラインコア、クラッシュしない）も同テストで
確認しています。

この `dsdePl` 突き合わせを追加する過程で実際のバグを1件発見・修正しました：
Python 側の 6×6 ヤコビアンは行優先（NumPy／C 順）でワイヤに乗りますが、
Fortran の `RESHAPE` は列優先で詰めるため、`usermat_py_hook.f` の単純な
`reshape(d36,[6,6])` は**転置行列**を黙って返していました。対称に近い
弾性ケースでは見えず、粘性・Mooney-Rivlin ケースで初めて（符号反転を
伴う）大きな食い違いとして露見 ― 明示的な `transpose(...)` で修正済みです。

## ANSYS でのビルド・使用（概略）

```
! ANSYS のユーザプログラマブル機能ツールチェーン（ANSUSERSHARED / usermat ビルド）
! でコンパイル・リンクし、モデル内で：
TB, USER, 1, 1, 6         ! プロパティ 6 個
TBDATA, 1, C10, C01, D1, eta, mtype, kUsePy
TB, STATE, 1, , 10        ! 状態変数 10 個（Fv 1:9, alpha 10）
```

ANSYS なしでの構文チェック：

```
gfortran -c -fsyntax-only -ffixed-line-length-132 usermat_biofilm.f
```

## 検証状況

- ✅ `gfortran`（`-fsyntax-only`）で警告なくコンパイル。
- ✅ **検証済み Abaqus UMAT とビット単位で一致**（20 変形状態で
  `|Δσ|=|ΔFv|=|ΔJe|=0`）。`crosscheck/` が両方の実 Fortran コアをコンパイルして比較。
  さらに敵対的スイープ 8017 ケースでも 0 ULP、フレーム不変性残差 4.9e-17。
  等方成長パッチ、整合接線 vs 中心差分 **2.97e-8**、ANSYS のせん断並び
  （`s12,s23,s13`）も確認済み。
- ✅ **ANSYS MAPDL 2022 R2 (v222) にて動作・収束確認済み。**
  インターフェイス引数（`var0`, `var1..var8`, `tsstif`, `epsZZ`）、自動時間刻み
  制御（`keycut`/`cutFactor`）、および接線剛性（`dsdePl`）の整合性を、
  `SOLID185` + `NLGEOM,ON` の一軸引張ベンチマークで実機検証。
  **引数リストはリリース依存**なので、別バージョンへ移る際はここを最初に
  確認すること（`var0` は `coords` と `defGrad_t` の間、`var1..var8` は
  `cutFactor` の後ろ）。
- 元となる Abaqus コアは検証済み（接線 vs FD ~2.4e-8、パッチ試験 13/13）。

## 注意点／次のステップ

1. ~~対象 ANSYS リリースの正確な `usermat` 引数リストを確認する。~~
   **2022 R2 (v222) で完了。** 他リリースへ移る際は再確認。
2. ~~ANSYS が F 摂動による `∂σ/∂ε` と同じ規約を期待するか確認。~~
   **収束したことで確認済み** ― ヤコビアンの規約が違えば Newton 収束が
   悪化または失敗するが、実際に収束している。
3. ~~成長項のソルバ内検証。~~ **2026-08-19/20 に完了。** 完全拘束単一要素
   （`F = I`、FE 求解なしに答えが予測できる）の4ケース（弾性/粘性 ×
   α=0.05/0.20）が閉形式と一致、加えて `KEYOPT(1,2)` ∈ {0,1,3}（B-bar／
   enhanced／simplified enhanced strain）のスイープでも体積ロッキングなし。
   詳細・実測値は [`apdl/README.md`](apdl/README.md) と
   [`apdl/RUNBOOK.md`](apdl/RUNBOOK.md)。
   - 副産物として、円筒シェル（歯表面を模した二層構造）での拘束成長の
     予備検討で **α=0.01 は収束、α=0.015 から歪みエラー**、収束する
     α=0.01 でも外周変位が一様でなく**2山パターン**を示すことを確認
     （[`assets/growth_cylinder_bulge.png`](../assets/growth_cylinder_bulge.png)、
     `apdl/extract_cylinder_bulge.py`）。座屈的挙動の可能性があるが、
     真の固有値解析では未確認 ― 興味深い所見として記録、結論はまだ。
4. ~~`PYTHON MATERIAL HOOK`（ISO_C_BINDING／ソケット）を実装する。~~
   **完了** ― 上記のとおり、実際の `usermat()` エントリポイントを通した
   `tests/test_usermat_kusepy_e2e.py` で end-to-end 検証済み。インライン
   コアは検証基準・自動フォールバックとして残す。プロトコルとガウス点
   1回あたりの通信内容は
   [`ch5_flow/flow_python_material_hook`](../ch5_flow/README.md) に図解。
5. Python 側の材料モデルは、いまも検証済み Fortran 則の NumPy 写し
   （`material_server.py` の `stress_core`）であり、較正済み JAX モデル
   （`JAXFEM/` ／ `material_models.py`）ではまだない ― 同じインターフェイス
   の裏側を差し替えるだけなので、配線の話ではなくローカルな置き換え。
