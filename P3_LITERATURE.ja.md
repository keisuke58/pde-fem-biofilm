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

## 5. 方法の注意（10 月 10 日の実行から）

1. **デッキの値と設定ファイルの値をビット一致で比べない。** APDL が 16 桁の小数を読んだ値は、Fortran が同じ文字列を読んだ値と最後のビットで違うことがある。
   点モデルの θ の照合が全点で落ち、最初の組み立てで keycut になった。照合は相対許容差つき（設定ファイルの `theta_tol 1e-12`）にした。
2. **一方向の連成では、組成は応力に出ない。** 場が量と成長を決め、点モデルが組成だけを決める形（組成のモード）では、較正済みの 5 種（CH）と
   2 種の case 6 で応力が全桁同じだった。組成を応力に出すには、点モデルから成長へ戻す式が要る（モード 9、生菌量の形）。
   成長の総量が同じでも分布が違えば応力は変わる（CH、α − 1 のシード平均を揃えて von Mises −8 %）。成長の駆動の形は結果に効く仮定として書く。
3. **応力が E に比例するのは条件つき。** 剛性 (φ² + f)·E_bio、ν 固定、外力なし、変位の拘束だけ、成長が応力によらない、のとき大変形でもちょうど比例する
   （CH で 10 Pa と 10 kPa、全要素で 1000 倍）。帰還、決まった大きさの外力、剛性の決まった基材との接触のどれかがあると崩れる。
   立方体の 4 条件の比は E によらないが、歯の形状では E を 2 通り回す。

## 6. まだやること

- ~~Alpkvist and Klapper 2007 の式~~：済み（Springer の全文、Eq. 4、6〜11、例 1）。
- ~~Pattem 2018 のヤング率~~：済み（Europe PMC の全文、Table 1）。
- Wanner and Gujer 1986 の原典は図書館経由（Alpkvist 2007 が同じ仮定を式で書いているので、引用は Alpkvist 2007 で足りる）。
- ~~剛性を Pattem 2018 の範囲（kPa）にする感度~~：CH で済み（10 月 10 日、`results/2026-10-cal5/w8_CH_cal_k0882{,_E10k}.json`）。
  E_bio を 10 Pa から 10 kPa にすると、全要素の von Mises がちょうど 1000 倍（要素ごとの比 999.9994〜1000.0004）、α は全桁同じ。
  シードの von Mises 平均は 4.32e−4 Pa から 0.432 Pa。
  **理由（小ひずみだからではない）：** 剛性は (φ² + f)·E_bio で ν は固定なので、応力は E_bio × (F だけの関数)。外力がなく境界は変位の拘束だけ、
  荷重は成長だけなので、釣り合いの式は E_bio で割っても同じになり、F は E_bio によらない。だから大変形の超弾性でも応力はちょうど E_bio に比例する
  （比のずれ 4e−7 は Newton の収束誤差。粘性が効いていないことの確認にもなる）。
  **この立方体・一方向の連成では、4 条件の応力の比は E によらない。** 崩れるのは次のとき：
  応力から成長への帰還（WP2 の恒常圧）、決まった大きさの外力（咬合力、流体圧）、剛性の決まった相手との接触（歯・インプラントの基材）、
  ν や f を E と一緒に変えないとき。歯の形状に移すときは基材との剛性の比が効くので、E を 2 通り回す。

## 7. 応力が成長を抑える則（恒常圧）の先例（10月11日、要旨レベル）

WP2 の帰還の形は α̇ = k_α φ max(0, 1 − p/p_h)（`RUN_WP2_IKMHIWI03.md`、`tests/test_homeostatic_growth.py`）。

| 文献 | 内容 | 確かめ方 |
|---|---|---|
| Basan, Risler, Joanny, Sastre-Garau, Prost 2009, HFSP J. 3(4):265–272 | 組織ごとに「恒常圧」p_h があり、そこで分裂と死が釣り合う。後の総説では成長率を (p_h − p) に線形とする。**WP2 の形はこれの片側版**（p > p_h で縮まない） | 要旨と二次資料 |
| Helmlinger et al. 1997, Nat. Biotechnol. | アガロースに閉じ込めた腫瘍スフェロイドは 45〜120 mmHg（6〜16 kPa）で成長が止まる。実験の根拠 | 要旨 |
| Ambrosi and Mollica 2002, Int. J. Eng. Sci. | F = F_e F_g の成長の連続体の枠組みで、応力が成長に効く腫瘍モデル。F_g = αI の形も同じ | 二次資料だけ。本文を読む |
| Delarue et al. 2016, Nat. Phys. 12:762 | 閉じ込めた酵母の集団で、増殖が接触圧（最大 約 0.65 MPa）を作り、圧で成長が遅くなる。**微生物で圧が成長を抑える直接の実験** | 要旨 |
| Weady, Shelley ほか 2024, arXiv:2405.10158 | 増殖する細胞集団で、成長が局所の応力で抑えられる連続体理論。応力感度で模様の間隔が決まる | 要旨 |
| Willmann and Wall 2022（arXiv:2110.07200） | バイオフィルムの多物理モデル。表面成長が栄養の流束で進み、表面の法線・接線応力で減る（侵食） | 要旨。成長則の出典は本文で確かめる |

**読み：** 「圧が成長を抑える」は腫瘍（Basan、Helmlinger）と微生物（Delarue）で実験・理論ともに先例がある。
WP2 の則は Basan 型を片側にしたもので、**形は先例から引ける。p_h の値は口腔バイオフィルムでの実測がなく仮定**。
感度（p_h を 2〜3 通り）が必須。Meisam に聞くこと：p_h を腫瘍や酵母の値から取るか、無次元で振るか。

## 8. 歯・インプラントの上のバイオフィルムの応力を計算した先行研究（10月11日）

**見つからなかった。** 歯科の有限要素は歯・歯根膜・骨・インプラントの応力が中心で、バイオフィルムは入っていない
（2023 年の Dental Materials の総説、147 本でもバイオフィルムの層はない）。近いもの：

| 文献 | 内容 | 違い |
|---|---|---|
| 2015、J. Biomech.（Biofilm and saliva affect the biomechanical behavior of dental implants） | バイオフィルムと唾液でインプラント部品の摩擦係数が変わり、FEA の応力分布が変わる | バイオフィルムは接触の性質としてだけ入る |
| S. mutans バイオフィルムの水滴除去の研究（Rmaile ほか、J. Dent. Res. 2014 と思われる） | ヤング率 約 380 Pa、粘弾性。CFD（有限要素）で壁のせん断応力 | 流体による除去で、成長の応力ではない。著者・数値は原典で確かめる |
| Mandal, Mottram, McGinty 2025, arXiv:2509.18065 | 医療インプラントの上のバイオフィルム成長（増殖菌・persister・死菌・EPS、栄養依存） | 力学なし、単一種 |
| Laspidou, Rittmann, Karamanos 2005 | UMCCA の密度から場所ごとのヤング率を作り FE に渡す | 一方向、口腔ではない |
| Klempt et al. 2024 | Hamilton 原理の成長と応力（この研究の土台） | 帰還なし、立方体と簡単な形状 |

**読み：** 「多種の組成から成長を作り、歯・インプラントの形状で応力を出し、応力を成長に戻す」組み合わせは、
見た範囲では先例がない。P3 の新しさはここに置ける。ただし検索は要旨レベルで、Scopus か Google Scholar で
"biofilm" AND ("finite element" OR "continuum") AND ("tooth" OR "implant" OR "dental") をもう一度引いて確かめる。

## 9. Meisam に聞くこと（10月11日のメモ、未送信）

修論の提出（10月20日）が終わってから聞く。修論のコメントへの返信には混ぜない。
**一度に全部は聞かない（ユーザー、10月11日）。** 会ったときに口頭で 1 問目（論文の形）だけ聞き、残りは話の流れで聞く。下の一覧は自分の整理用。

1. **論文の形：** 応力の値を予測する論文ではなく、「どの仮定が応力の分布を決めるか」を調べる感度の論文にしてよいか。
   検証できる実測（歯の上のバイオフィルムの応力）はない。
2. **恒常圧 p_h の値：** 口腔バイオフィルムの実測はない。腫瘍（Helmlinger 1997、6〜16 kPa）や酵母（Delarue 2016）の値を借りるか、
   無次元にして 2〜3 通り振るか。
3. **種ごとの硬さ：** 種ごとの実測はない（Pattem 2018 は条件ごとの 0.6〜41 kPa だけ）。全種同じ E にして、
   組成は成長だけを通して効く、とするので足りるか。
4. **E の値：** Klempt 2024 の 10 Pa と、Felix が勧める 10 kPa のどちらを基準にするか。立方体では応力は E に比例するが、
   歯の形状で基材との接触や帰還があると比例しなくなるので、2 通り回す予定。
5. **成長の駆動：** 点モデルの生きている量 Σφ_iψ_i で成長させる形（κ は Eq. 36 の成長量に合わせた、論文の値ではない）でよいか。
   先例は Alpkvist and Klapper 2007。
6. **投稿先と共著：** BMMB か J. Mech. Behav. Biomed. Mater. でよいか。Meisam（と Felix）に共著に入ってもらえるか。

英語の短い版（メール用の下書き）：

> For a paper after the thesis, I would like to study which modelling assumptions set the stress
> distribution in an oral biofilm on a tooth or implant, rather than predict stress values, since no
> measured stresses exist for comparison. Open points: (1) is a sensitivity study of this kind a
> suitable form; (2) the homeostatic pressure p_h, borrowed from tumour or yeast data or varied in
> dimensionless form; (3) one Young's modulus for all species; (4) E = 10 Pa or 10 kPa as the base
> case; (5) growth driven by the living biomass of the point model; (6) the journal and co-authors.

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
