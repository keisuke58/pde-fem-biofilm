# JAXFEM — Klempt 2024 バイオフィルム PDE スイート

Klempt et al. (2024) *Biomech Model Mechanobiol* の支配方程式を
**JAX + NumPy + SciPy** で再実装した純数値シミュレーション群。
Abaqus/ANSYS のメッシュは使わず、有限差分（FD）で φ-c-α 場を解く。

## 位置づけ

```
このフォルダ (JAXFEM)          →  ../../../nife/masterarbeit_ansys_fem/
数式の正しさを確かめる実験台       修論本番: 同じ構成則を ANSYS Gauss 点へ結合
φ-c-α フォワード計算のみ         coupling_prototype/ の UMAT・USERMAT
```

## スクリプト一覧

### 核心 — Klempt 2024 再現

| スクリプト | 内容 | 状態 |
|-----------|------|------|
| `felix_complete_reproduction.py` | Eq.34–36 完全再現（Allen-Cahn + logistic-Monod + 走化性） | **PASS** |
| `felix_exact_check.py` | 厳密ベンチマーク（パラメータ一致確認） | **PASS** |
| `klempt_pde_jax.py` | φ-c-α PDE JAX実装（Phase 0a） | **PASS** |
| `benchmark_klempt2024.py` | Klempt Table 2 全パラメータ自動検証 | PASS |

### Hamilton 原理 PDE コア

| スクリプト | 内容 |
|-----------|------|
| `core_hamilton_1d.py` / `core_hamilton_2d.py` | Hamilton 原理ベース PDE（1D/2D） |
| `core_hamilton_1d_nutrient.py` / `core_hamilton_2d_nutrient.py` | 栄養場連成版 |
| `core_rd.py` | 反応拡散（参照実装） |

### 応力・逆問題

| スクリプト | 内容 |
|-----------|------|
| `phase1_klempt_stress.py` | α → eigenstrain → FEM応力（Phase 1 PASS） |
| `phase3_5species_stress.py` | 5菌種版応力計算 |
| `solve_stress_2d.py` / `solve_stress_2d_vem.py` | 2D応力ソルバー（FEM/VEM） |
| `param_estimation_1d.py` / `param_estimation_hamilton_1d.py` | 1D パラメータ推定 |

### デモ・可視化

| スクリプト | 内容 |
|-----------|------|
| `demo_hamilton_1d.py` / `demo_hamilton_2d.py` | Hamilton PDE デモ |
| `pinn_1d_demo.py` | PINN デモ（参考） |
| `thesis_style.py` | 論文スタイルのプロット設定（共通） |

## 保存済み中間結果

```
klempt_phi_final.npy    (50,50) — φ 最終フィールド
klempt_c_final.npy      (50,50) — 栄養場 c
klempt_alpha_final.npy  (50,50) — 固化変数 α
```

## フェーズ状況

| Phase | スクリプト | 状態 |
|-------|-----------|------|
| 0a | `klempt_pde_jax.py` | ✅ PASS |
| 1 | `phase1_klempt_stress.py` | ✅ PASS |
| 0b | `phase0b_nsp_klempt_connection.py` | ✅ 実装済 (2026-06-26 確認) |
| 2 | `umat_biofilm_visco.f` + `umat_tangent_test/` + `phase2_patch_test.py` | ✅ 実装済 (2026-06-26, exact consistent tangent, Fortran検証 vs FD ~2.9e-8) |
| 3 | `phase3_5species_stress.py`, `phase3b_voigt_stress.py` | ✅ 実装済 (2026-06-26 確認) |

詳細な結果ログ → [RESULTS.md](RESULTS.md)

## 実行

```bash
cd /home/nishioka/IKM_Hiwi/FEM/JAXFEM

# Klempt 2024 Eq.34-36 完全再現（12000ステップ・約2分）
python felix_complete_reproduction.py

# 厳密ベンチマーク
python felix_exact_check.py

# Phase 1: 応力フィールド確認
python phase1_klempt_stress.py
```

## 関連

- 修論本番フォルダ: [`../../../nife/masterarbeit_ansys_fem/`](../../../nife/masterarbeit_ansys_fem/)
- Klempt 2024 DB: [`../../notes/klempt2024_hamilton_biofilm_JP.pdf`](../../notes/klempt2024_hamilton_biofilm_JP.pdf)
- 上位 FEM スクリプト群: [`../`](../) （Abaqus `.inp` / `.py` 一式）
