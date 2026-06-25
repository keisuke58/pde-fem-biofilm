# JAXFEM — Klempt 2024 実装進捗・結果ログ

最終更新: 2026-06-23

---

## フェーズ一覧

| Phase | 状態 | スクリプト | 目的 |
|-------|------|-----------|------|
| 0a | **PASS ✓** | `klempt_pde_jax.py` | φ-c-α PDE JAX実装（Klempt 2024 支配方程式） |
| 1  | **PASS ✓** | `phase1_klempt_stress.py` | α→eigenstrain→FEM応力（力学連成確認） |
| 0b | **実装済 ✓** (2026-06-26 確認) | `phase0b_nsp_klempt_connection.py` | NSP Hamilton ODE（5菌種版）φ_i(x,y)→α接続 |
| 2  | **実装済 ✓** (2026-06-26, commit 83660be) | `umat_biofilm_visco.f` (BIOFILM_STRESS_CORE + F摂動), `umat_tangent_test/`, `phase2_patch_test.py` | consistent tangent DDSDDE（exact, Sun et al. 2008 摂動法。Fortran/Python 検証 13/13: 弾性/粘弾性/MR とも vs FD ~2.4-2.9e-8。実 Abaqus 1要素ジョブ COMPLETED: maxIter2/cutback0/増分自動拡大。PNEWDT cutback・MR 対応済） |
| 3  | **実装済 ✓** (2026-06-26 確認) | `phase3_5species_stress.py`, `phase3b_voigt_stress.py` | TMCMC posterior θ_MAP → DI(x,y) → E(x,y) → 5菌種版応力 |

---

## Phase 0a — φ-c-α PDE JAX実装

**スクリプト**: `klempt_pde_jax.py`  
**状態**: PASS ✓

### 実装内容

Klempt et al. (2024) Biomech Model Mechanobiol の支配方程式をJAXで実装。

| 変数 | 支配方程式 |
|------|-----------|
| `φ(x,y,t)` | `dφ/dt = D_φ·Δφ + μ·φ·(1−φ)·c/(k_M+c)` |
| `c(x,y,t)` | `dc/dt = D_c·Δc − γ·φ·c/(k_M+c)` |
| `α(x,y,t)` | `dα/dt = k_α·φ·c/(k_M+c)` |

境界条件（Klempt Case 1相当）：
- `φ`：全壁 Neumann（ゼロフラックス）
- `c`：右上コーナー Dirichlet `c=1`、他は Neumann
- `α`：全壁 Neumann

### 保存ファイル

| ファイル | shape | 値域 |
|---------|-------|-----|
| `klempt_phi_final.npy` | `(50,50)` | 0.00044 – 0.733 |
| `klempt_c_final.npy` | `(50,50)` | 0.0069 – 1.0 |
| `klempt_alpha_final.npy` | `(50,50)` | 0.00022 – 0.759 |

---

## Phase 1 — α→eigenstrain→FEM応力

**スクリプト**: `phase1_klempt_stress.py`  
**状態**: PASS ✓

### パイプライン

```
klempt_phi_final.npy  \
klempt_c_final.npy     → E(φ) + ε_growth(α) → solve_2d_fem() → σ_vm
klempt_alpha_final.npy /
```

- `ε_growth = α / 3`（等方膨張：体積膨張=3×線形）
- `E(φ) = E_max·(1−φ)² + E_min·φ`（DI類比）
- `η = 0`：弾性のみ（Klempt 準静的線形弾性と同一条件）

### 物理確認結果

| 領域 | σ_vm [Pa] |
|------|---------|
| バイオフィルム内（α > 5%最大値） | **29.1 Pa** |
| 基質（マトリクス、α ≤ 5%最大値） | **18.1 Pa** |

**合格基準 PASS**：`σ_vm_biofilm (29.1) > σ_vm_matrix (18.1)`

**物理解釈**（Lamé問題アナロジー）：  
均一膨張球が周囲に拘束されると内部に均一圧縮応力が集中する。  
バイオフィルムが基質に拘束されて成長する場合も同じ機構——Klempt Fig 5相当の応力分布を再現。

---

## 次ステップの候補

### Phase 0b（推奨：次の一手）
NSP Hamilton ODE（5菌種版）の実装。`φ_i(x,y)` → `α`の計算 → Klemptパイプラインと接続確認。  
修論の多菌種拡張（Nishioka thesis Ch.3）との架け橋。

### Phase 2（Soleimani会議後）
UMAT `DDSDDE`（consistent tangent）の厳密化。  
現状は弾性接線のみの近似版。Newton–Raphson収束の改善が目的。  
→ `soleimani_questions.md` Q2を会議で確認してから実装判断。

### Phase 3
TMCMC posterior θ_MAP → `DI(x,y)` → `E(x,y)` → 5菌種版応力計算。  
修論 Ch.3 Bayes推定結果を力学場に接続するエンドツーエンドパイプライン。  
Phase 0b と Phase 1 がどちらも PASS してから着手。
