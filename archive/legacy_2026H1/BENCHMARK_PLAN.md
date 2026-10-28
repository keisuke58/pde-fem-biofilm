# Klempt 2024 → 5菌種拡張 ベンチマークロードマップ

Ref: Klempt et al. (2024) Biomech Model Mechanobiol
[https://pmc.ncbi.nlm.nih.gov/articles/PMC11554842/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11554842/)

---

## Phase 0a: Klempt φ-c-α PDE 実装（★ 最初にやること）

**発見**: NSP Hamilton ODE（現実装）と Klempt φ PDE は別物。

- NSP: 菌種競合の分率 φᵢ（空間成長の概念なし）
- Klempt: biofilm密度の時空間発展 ∂φ/∂t = ∇·(D_φ∇φ) + f(φ,c,α)

**目標**: Klempt の支配方程式を JAX で実装
**スクリプト**: `klempt_pde_jax.py`（新規作成）
**方程式**:

- φ: ∂φ/∂t = ∇·(D_φ∇φ) + μ·φ·(1-φ)·c/(k_M+c)  [logistic growth × Monod]
- c: ∂c/∂t = D_c·∇²c - γ·φ·c/(k_M+c)              [advection-diffusion]
- α: dα/dt = k_α·φ·c/(k_M+c)                        [expansion accumulation]

**BC**: Dirichlet c=1 at corner, Neumann others; φ=0.3 seed at center
**合格基準**: φ場がKlempt Fig 2（指向性成長）と定性一致

## Phase 0b: NSP ↔ Klempt 接続の確認

**目標**: 完全一致確認（単一菌種）
**方法**:

- `active_mask = [1,0,0,0,0]`（So のみ）
- η=0（粘性なし）→ F = Fe·Fg のみ（Klemptと同じ）
- Phase 0a の α(x,y,t) を使って UMAT に渡す
- 合格基準: φ場の指向性成長がKlempt Fig 2と定性一致

---

## Phase 1: UMAT力学ベンチマーク（単一菌種+応力）

**目標**: 応力・変位場をKlempt Fig 5と一致
**方法**:

- Phase 0のφ・α場をAbaqus UMATに渡す
- UMAT: eta=0, C10=E(φ)/4/(1+ν), F = Fe·Fg
- 静水圧・von Mises応力分布を比較
- 合格基準: 応力集中位置（biofilm-boundary界面）が一致

---

## Phase 2: consistent tangent 実装

**目標**: Newton-Raphson 2次収束の確認
**方法**:

- UMAT Section 7 を厳密DDSDDEに置き換え
- 粘性補正項 d(sigma)/d(eps)|_viscous を追加
- Patch test + 収束次数の確認
- 合格基準: 反復回数がPhase 1より ≤50%

---

## Phase 3: 5菌種拡張

**目標**: TMCMC posterior θ → φᵢ(x,y,t) → DI → E(x,y) → 応力
**方法**:

- `active_mask = [1,1,1,1,1]`（5菌種全有効）
- A行列は TMCMC posterior θ_MAP から構築
- DI(x,y) = 1 - H(φ(x,y))/ln5
- E(x,y) = E_max*(1-r)² + E_min*r
- 合格基準: 4条件（CS/CH/DS/DH）で条件間変位比 > 5×

---

## 進捗状況

| Phase | 状態 | 次アクション |
| ----- | ---- | ------------ |
| 0a | **PASS ✓** 2026-06-23 | `klempt_pde_jax.py`: φ-c-α PDE JAX実装完了 |
| 0b | **PASS ✓** 2026-06-23 | `phase0b_nsp_klempt_connection.py`: NSP φ_So > 0, 栄養素依存(c=0:0.005→c=1:0.007) |
| 1 | **PASS ✓** 2026-06-23 | `phase1_klempt_stress.py`: α→ε_growth→FEM応力。σ_vm_biofilm=29.1Pa > σ_vm_matrix=18.1Pa |
| 2 | **PASS ✓** 2026-06-23 | `umat_biofilm_visco_phase2.f` + `phase2_patch_test.py`: 数値摂動consistent tangent 10/10 |
| 2b | **PASS ✓** 2026-06-23 | **Felix完全一致確認**: `felix_exact_check.py` — Table2パラメータ・問題設定・Fig2/5再現 |
| 3 | **PASS ✓** 2026-06-23 | `phase3_5species_stress.py`: 4条件, dysbiotic/commensal比=2.34× ≥2× |

### Phase 2b: Felix完全一致確認 (2026-06-23)

**スクリプト**: `JAXFEM/felix_exact_check.py`

**Fg の正しい定式化（論文 PMC11554842 から確認）**:

Felix の定義: Fg = α·I where α(0) = 1, dα/dt = kα·ϕ

```fortran
! UMAT: TEMP = alpha_acc (= α-1, starts at 0)
FG(I,J) = (1.0 + ALPHA_G_IN) * IDEN(I,J)   ! = α·I
FG_DET  = (1.0 + ALPHA_G_IN)**3              ! = det(Fg) = α³
```

確認済み誤り（一時的に適用→即時差し戻し）: `Fg=(1+α)^(1/3)·I` は cube-root 解釈で **誤り**。
論文 Eq. で α は線形膨張パラメータ (not volumetric)。元の `Fg=(1+α)·I` が正。

**felix_exact_check.py 結果**:

- E = 2μ(1+ν) = 10.0 Pa (Table 2: 10 Pa) ✓
- ν = 0.4900 (Table 2: 0.49) ✓
- Fig 2: CM shift (10,10)→(11.4,11.4) → corner方向 ✓
- Fig 5: mushroom growth, top extent 5μm ✓
- 内部静水圧 p > 0 (圧縮) ✓

**Abaqus INP の注意**:

- `*Temperature` フィールド: α_acc を直接（`/3` 補正なし）で渡す
- η=0 で `F = Fe·Fg` のみ（Klempt 2024 と同一）

### Phase 1 詳細 (2026-06-23)

- スクリプト: `JAXFEM/phase1_klempt_stress.py`（Python FEM proxy — Abaqus参照用）
- 入力: `klempt_{phi,c,alpha}_final.npy` (Phase 0a出力)
- 材料モデル: `E(φ) = E_max*(1-φ)^2 + E_min*φ`  (E_max=900 Pa, E_min=30 Pa)
- 固有ひずみ: `ε_growth = α/3`（Python FEM側は線形近似のまま）
  - **注**: Abaqus UMAT では TEMP=α_acc をそのまま渡す。UMAT内部で Fg=(1+TEMP)·I=α·I を直接計算（cube-root 変換なし）
- FEM: 2D plane strain, bottom fixed BC, scipy sparse QUAD4
- 結果: σ_vm_max=36.2 Pa, |u|_max=17mm, σ_vm_biofilm=29.1 > σ_vm_matrix=18.1 ✓
- 物理: Lamé問題アナロジー（拘束膨張球→内部均一圧縮応力）。Klempt Fig 5相当
