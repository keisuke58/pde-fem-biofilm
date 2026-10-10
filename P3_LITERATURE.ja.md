# P3 の文献調査：生菌量で成長を駆動すること、死菌と EPS、種ごとの力学（2026年10月10日）

要素のモード 9（`ansys_usermat/apdl/callsite/phi_mode_exec.inc`、tmcmc202601 の
`docs/fem_coupling/P3_point_model_alpha_design.md`）の式
ln(1 + α) の増分 = κ χ ∫ Σ φ_i ψ_i dt を、論文でどう支えるかの調査。
修論（12 日凍結）には入れない。P3 の本文の材料。

**確かめ方の印：** 「本文」＝式や文を全文で読んで確かめた。「要旨」＝要旨や二次資料だけで、式や数値は未確認。

## 1. 体積の成長を菌の生産で決める先例（あり、直接使える）

| 文献 | 何を言っているか | 確かめ方 |
|---|---|---|
| Feng et al. 2021, Bull. Math. Biol.（`references/`） | 2 種の口腔バイオフィルム。体積分率の和が 1（Eq. 14）として、成長速度の発散が種ごとの生産の和：∇·u = g₁/ρ_s + g₂/ρ_v（Eq. 15a）。成長は「バイオフィルム内のバイオマスの生産と減少に駆動される移流」（2.3 節） | 本文 |
| Alpkvist and Klapper 2007, Bull. Math. Biol. 69:765 | 多次元・多種の連続体モデル。成分は菌の種、EPS、不活性バイオマス（2 節）。各成分は非圧縮（Eq. 4）、体積分率の和が 1（Eq. 6）、全成分が同じ速度 u で移流（Eq. 7）、Darcy 則 u = −λ∇p（Eq. 8）。質量保存 ∂ϑ_i/∂t + ∇·(uϑ_i) = g_i/ρ_i*（Eq. 9）を和にとって **∇·u = Σ g_i/ρ_i***、圧力は −λ∇²p = Σ g_i/ρ_i*（Eq. 10） | 本文（Springer の全文） |
| Wanner and Gujer 1986, Biotechnol. Bioeng. 28; Wanner and Reichert 1996, 同 49:172 | 1 次元多種モデルの原典。移流速度が新しいバイオマスの生産と直接結びつく。崩壊は不活性物質を生み、それが膜の中にたまる | 要旨・二次資料 |

**モード 9 との対応。** 物質点で見ると d ln J/dt = ∇·u なので、Wanner–Gujer 型の式は「物質点の体積の伸び率 = その点の生産率の和」と読める。
モード 9 は、生産率を g_i/ρ_i = κ φ_i ψ_i（生菌量に比例する 1 次の成長、全種で同じ κ）と置いたものに当たる。これで本文に一文書ける：
体積の成長を菌の生産の和で決めるのは Wanner–Gujer 以来の標準の仮定で、口腔の 2 種でも Feng 2021 が使っている。

**Klempt 2024 との違い（書くときに要る）。** Klempt 2024 の局所式（Table 1）は ln α の増分が場の φ の増分に比例する形。
Wanner–Gujer 型では、物質点の φ は和が 1 のまま動かず、体積は生産率で伸びる。モード 9 は後者の形で、駆動を点モデルの生菌量にした。
Klempt 2024 の形から離れることは、仮定として明記する（修論では離れないと決めた。P3 は別の判断）。

**Alpkvist and Klapper 2007 の例 1（3.1 節）がいちばん近い。** 活性バイオマス ϑ₁ と不活性 ϑ₂（休止細胞と EPS）で、同じ密度 ρ*。
g₁ = ϑ₁ρ*(μ c/(k_O2 + c) − k_d − k_i)、g₂ = ϑ₁ρ* k_i（k_i は活性から不活性への移り、k_d は CO₂ への分解）。
和をとると ∇·u = ϑ₁(μ c/(k_O2 + c) − k_d)：**体積の伸びは活性バイオマスだけに比例し、不活性への移り k_i は体積を変えない**（死んだ分は不活性として残る）。
モード 9 の「生菌量 Σφ_iψ_i に比例して伸び、死菌は体積を保つ」はこの形そのもの。

**違いとして残る点。**
- Alpkvist 2007 の伸び率には栄養の Monod 因子 c/(k + c) が掛かり、分解 k_d で縮みうる。モード 9 の α の式には栄養の因子がなく（栄養は点モデルの c* を通してだけ効く）、縮まない。
  栄養の因子を α の式にも入れるかは、P3 の感度の一つにする。
- Wanner–Gujer と Feng は、生産率が種ごとに違う（Monod の μ_i、栄養や乳酸に依存）。モード 9 は κ が 1 つで、種の違いは点モデルの φ_i ψ_i の動きから入る。
  種ごとの κ_i にすると、それは種重み（10 日に外した仮定）と同じ自由度になる。
- Wanner–Gujer と Feng は Euler 系の移流で境界が動く。モード 9 は Lagrange 系の成長テンソル（F_g = α I、Rodriguez et al. 1994 の乗法分解）で、境界の形は場の φ が決める。

## 2. 死菌と EPS は体積を占めるか

| 文献 | 何を言っているか | 確かめ方 |
|---|---|---|
| Klempt et al. 2026, Arch. Appl. Mech. 96:164（arXiv 2509.01274） | φ_i は種 i が占める体積（Eq. 6）、ψ_i はそのうち生きている割合（Eq. 7）、φ̄_i = φ_i ψ_i は生きた菌の体積（Eq. 8）、死菌は τ_i = φ_i(1 − ψ_i)。Σφ = 1（Eq. 9）なので、ψ が下がっても死菌は φ_i の中で体積を占め続ける。力学は使わず「単一の物質点」のモデル | 本文（Europe PMC の全文、前半） |
| Wanner and Reichert 1996; Reichert and Wanner 1997 | 崩壊したバイオマスは不活性物質になり膜の中にたまる（体積を占める）。成長による体積の膨張がバイオマスを表面へ押し出す | 要旨 |
| Laspidou and Rittmann 2004/2005（UMCCA） | 活性菌・不活性バイオマス・EPS を別成分として扱い、下層が上層の 5〜10 倍密になる。活性汚泥（SRT 20 日）で活性 約 28 %、不活性 約 40 %、EPS 約 31 % | 要旨（口腔ではない） |
| Feng 2021（2.1 節） | EPS と孔は陽には扱わない（Alpkvist and Klapper 2007 と同じ）。EPS を別成分にするなら Xavier et al. 2005 | 本文 |
| Asally et al. 2012, PNAS（B. subtilis） | 局所的な細胞死が圧縮応力の逃げ道になり、しわを作る | 要旨・二次資料 |

**モード 9 との対応。**
- 「死菌は体積を占めるが、生産はしない」は Klempt 2026 の変数の定義そのもので、Wanner–Gujer の不活性バイオマスとも同じ。
  だから駆動を Σφ_i ψ_i（生産する量）にし、Σφ_i（占める量）にしないのは筋が通る。10 日に決めた生菌形はこれで支えられる。
- **EPS は別成分にしていない**（点モデルも場も）。Feng 2021 と Alpkvist and Klapper 2007 と同じ単純化として限界に書く。
  EPS の体積は暗に φ に含まれる（生菌の生産に比例してできる、と読める）。
- **死菌の力学は生菌と同じとしている。** Asally 2012 は死菌の場所が力学的に違う（応力を逃がす）ことを示すので、これも限界として書く。

## 3. 種ごと・条件ごとの力学の実測（少ない。仮定のまま）

| 文献 | 何を言っているか | 確かめ方 |
|---|---|---|
| Pattem et al. 2018, Sci. Rep. 8:5691 | 唾液由来の多種バイオフィルム（ハイドロキシアパタイト上、3・5 日）を OCT と AFM で。スクロース 5 % は 0.1 % よりヤング率が低く、付着力が高い（p < 0.0001）。差は EPS の量による、と著者（EPS は測っていない、と本人たちも書いている） | 本文（Europe PMC の全文） |
| Vinogradov et al. 2004, Biofilms（S. mutans） | 粘弾性流体のようにふるまい、Burger モデルで合う | 要旨 |
| Waters et al. 2014（NIST、S. mutans） | せん断レオメータでのその場測定。基質は多糖・DNA・細胞の残骸 | 要旨 |
| Gloag et al. 2021, Front. Cell. Infect. Microbiol.（S. gordonii） | 回転円板レオメトリー。アルギニンの効果を見る | 要旨 |
| Peterson et al. 2015, FEMS Microbiol. Rev. 39:234 | 総説。粘弾性が抗菌薬の浸透と除去を決める | 要旨 |
| 大腸菌の基質変異（PMC12264770、2025） | 基質の組成（curli、セルロース）が剛性を強く決める | 要旨（口腔ではない） |

**Pattem 2018 のヤング率（Table 1、分布の 2 つの山、kPa）。** AFM、10 µm のホウケイ酸の球、Hertz、ν = 0.5、力 約 6 nN、1 条件 768 本の力-距離曲線。± が SD か SEM かは書かれていない。

| 条件 | 1 つ目の山 | 2 つ目の山 |
|---|---|---|
| スクロース 0.1 %、3 日 | 14.35 ± 1.75 | 41.05 ± 0.97 |
| スクロース 5 %、3 日 | 0.55 ± 0.02 | 2.57 ± 0.17 |
| スクロース 0.1 %、5 日 | 1.17 ± 0.08 | 46.07 ± 1.87 |
| スクロース 5 %、5 日 | 0.56 ± 0.06 | 6.66 ± 0.46 |

- 口腔の多種バイオフィルムの局所の剛性は **0.5〜46 kPa**。Klempt 2024 の E = 10 Pa はこれより 2〜4 桁小さく、Felix が 10 月 5 日に勧めた 10 kPa はこの範囲に入る。
- AFM は µm スケールの局所の値（菌体に当たれば硬い山、EPS なら柔らかい山）。FEM の連続体の剛性とは尺度が違うことは書き添える。

**分かったこと。**
- Fn や Pg を含む口腔バイオフィルムで、**種ごとの剛性や、Commensal と Dysbiotic の力学の違いを測った文献は見つからなかった**（3 回の検索、うち 2 回は深い検索）。
- 一番使えるのは「EPS の量が剛性を変える」（Pattem 2018、大腸菌の研究）こと。種ごとの差は、種ごとの EPS の作り方を通して入ると読むのが自然だが、その量的なデータもない。
- **結論：種ごとの剛性は P3 でも仮定のまま、感度として出す。** 条件の差は、まず成長（モード 9）から出す。剛性は 4 条件で共通にする。

## 4. P3 の本文で使う一文（下書き、英語）

- 成長の式：Volumetric growth driven by the production of biomass, ∇·u = Σ g_i/ρ_i (Alpkvist and Klapper 2007,
  Eq. 10), is the standard assumption of continuum biofilm models since Wanner and Gujer (1986) and is used for two
  oral species by Feng et al. (2021). In Alpkvist and Klapper's active–inert example the expansion is proportional to
  the active biomass alone, inactivation leaving the volume unchanged. At a
  material point it reads d ln J/dt = Σ g_i/ρ_i; I take the production proportional to the living volume of the
  calibrated point model, g_i/ρ_i = κ φ_i ψ_i, with one κ for all species and conditions.
- 死菌：Dead bacteria keep their volume, τ_i = φ_i(1 − ψ_i) (Klempt et al. 2026), but do not produce; growth is therefore
  driven by Σ φ_i ψ_i, not by Σ φ_i.
- 限界：EPS is not a separate constituent (as in Alpkvist and Klapper 2007 and Feng et al. 2021), dead cells are
  mechanically the same as living ones, and the stiffness is the same in all conditions: species-resolved mechanical
  data for these oral species were not found.

## 5. まだやること

- ~~Alpkvist and Klapper 2007 の式~~：済み（Springer の全文、Eq. 4、6〜11、例 1）。
- ~~Pattem 2018 のヤング率~~：済み（Europe PMC の全文、Table 1）。
- Wanner and Gujer 1986 の原典は図書館経由（Alpkvist 2007 が同じ仮定を式で書いているので、引用は Alpkvist 2007 で足りる）。
- 剛性を Pattem 2018 の範囲（kPa）にする感度：4 条件で共通の E を 10 Pa と 10 kPa の 2 通り。

## 出典

- Feng et al. 2021：`references/Feng2021_symbiotic_biofilm_BMB.txt`（2.1、2.3、2.5 節、Eq. 1、14、15）
- Klempt et al. 2026：<https://europepmc.org/article/PMC/PMC13391730>、<https://arxiv.org/abs/2509.01274>
- Alpkvist and Klapper 2007：<https://doi.org/10.1007/s11538-006-9168-7>（全文は Springer、大学のネットワーク経由）、<https://scholarworks.montana.edu/handle/1/13263>
- Pattem et al. 2018 の全文：<https://europepmc.org/article/PMC/PMC5890245>
- Wanner–Gujer 系の説明：<https://ejde.math.txstate.edu/conf-proc/10/e1/eberl-tex>、Reichert and Wanner 1997 の要旨 <https://acnpsearch.tweb-dev.unibo.it/singlejournalindex/10151725>
- Laspidou and Rittmann：<https://www.cambridge.org/core/journals/biofilms/article/modeling-biofilm-complexity-by-including-active-and-inert-biomass-and-extracellular-polymeric-substances/9D8EFB750C37EA851E275E3BA9879610>、<https://link.springer.com/article/10.1007/s11783-009-0008-5>
- Asally et al. 2012：<https://pmc.ncbi.nlm.nih.gov/articles/PMC3503208>
- Pattem et al. 2018：<https://pmc.ncbi.nlm.nih.gov/articles/PMC5890245/>
- Vinogradov et al. 2004：<https://scholarworks.montana.edu/handle/1/13389>
- Gloag et al. 2021：<https://www.frontiersin.org/articles/10.3389/fcimb.2021.784388/pdf>
- Peterson et al. 2015：<https://pmc.ncbi.nlm.nih.gov/articles/PMC4398279>
- 基質の組成と剛性（大腸菌）：<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12264770/>
- 口腔の多種バイオフィルムの構造：<https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.01716/pdf>
