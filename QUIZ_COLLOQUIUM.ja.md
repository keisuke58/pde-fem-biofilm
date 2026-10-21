# コロキウム用のクイズ（2026年10月11日〜）

ユーザーが「クイズ」「テスト」と言ったら、ここから 5 問出す。前に出していない問題を優先し、
答えは聞かれてから出す。答えは一言で言える形を基本にし、数字は補足にとどめる。
出した日と、答えられなかった問題を下の記録に足していく（答えられなかった問題は次にもう一度出す）。

## A. 全体の流れ（やさしい）

| # | 問題 | 答え |
|---|---|---|
| A1 | 扱ったバイオフィルムの菌は何種類か | 5 種（So, An, Vei, Fn, Pg） |
| A2 | 応力の計算に使ったソフトは | ANSYS（2022 R2、Oliver の要素を USERMAT で） |
| A3 | 組成のパラメータを実験に合わせた方法は | TMCMC（ベイズ推定） |
| A4 | 比べた応力の種類は | von Mises 応力 |
| A5 | 土台の論文は | Klempt et al. 2024（Hamilton 原理の成長モデル） |
| A6 | 実験データは誰の何のデータか | Heine らの 4 条件（CS, CH, DS, DH）の組成 |
| A7 | 4 条件の C/D、S/H は何の違いか | C = 共生的、D = 病的（dysbiotic）、S = 静置、H = 流れあり（HOBIC） |
| A8 | Abaqus の計算は修論に入っているか | 入っていない。慶應での続き |

## B. 成長モデル（基本）

| # | 問題 | 答え |
|---|---|---|
| B1 | 成長すると体はどの方向に膨らむか | 全方向に同じだけ（等方、F_g = αI） |
| B2 | α の最初の値は | 1（成長量は α − 1） |
| B3 | 成長の速さを決める記号は | k_α |
| B4 | β は何の係数か | φ の拡散の係数（Eq. 34） |
| B5 | β = 0.02 mm²/T* はどこから | 論文の β = 2 を 2 mm の立方体に直した（仮定） |
| B6 | Oliver の例の β（1e−4）を使わない理由 | 拡散の長さがメッシュより小さく、応力が収束しない |
| B7 | バイオフィルムのヤング率（Klempt 2024） | E = 10 Pa（ν = 0.49） |

## C. 数値の検証（修論の芯）

| # | 問題 | 答え |
|---|---|---|
| C1 | 誤差がいちばん小さいメッシュは | 24³（シード平均で 1 % 弱） |
| C2 | 誤差は何と比べたか | 別に書いた正確な解（有限体積法、細かいメッシュから外挿） |
| C3 | 応力は何で確かめたか | 別に書いた hex8 の有限要素コード |
| C4 | 時間刻みを大きくしすぎると | 隣の点で符号が交互になって壊れる（計算が止まる） |
| C5 | 短い計算（T* = 1.1）で見えない理由 | 壊れるモードが最初は小さすぎて、育つ前に終わる |
| C6 | 安定の限界がメッシュでばらつく理由 | 30 番目の近傍の選び方（同じ距離の 7 点から 5 点） |
| C7 | 近傍を 32 点にすると | 対称になり限界が上がる。精度はほぼ同じ |
| C8 | 近傍は多いほど精度がよいか | よくならない。精度を決めるのはメッシュ |

## D. 5 種モデルと結果

| # | 問題 | 答え |
|---|---|---|
| D1 | 5 種モデルは要素のどこで計算されるか | Gauss 点ごと |
| D2 | 場が成長量を決める形では、組成は応力に出るか | 出ない（全桁同じ） |
| D3 | 生きている量で成長させると応力は | 8 % 下がった（成長の総量はそろえてある） |
| D4 | 総量が同じなのに応力が変わる理由 | 成長の場所の分け方（ばらつき）が変わったから |
| D5 | κ は論文の値か | 違う。成長量を Eq. 36 にそろえて決めた仮の値 |
| D6 | 応力は E に比例するか | この立方体ではちょうど比例。帰還や外力があると崩れる |

## E. 事務

| # | 問題 | 答え |
|---|---|---|
| E1 | コロキウムは成績に入るか | 入らない。合否のみ（Studienleistung） |
| E2 | 成績は何で決まるか | 修論本体。審査員 2 人の平均 |
| E3 | 審査員は | Junker 先生（主）、Wangenheim 先生（副） |
| E4 | コロキウムの日時と形式 | 2026年12月16日（水）13:00、IKM 320 号室、発表 20 分＋質疑 10〜20 分 |

## F. 英語で答える基本（10月13日、日本語と英語）

「Q&A」と言われたら、ここの問題を答えつき・英語と日本語で見せる（基本から順に）。新しい問答もここに足す。

答えは英語で言えるように練習する。日本語は意味の確認用。

| # | 問題 | 答え（EN ／ JA） |
|---|---|---|
| F1 | What is the main contribution? | I brought a calibrated multispecies point model into a finite element growth model, and verified each step. ／ 較正した多菌種の点モデルを有限要素の成長モデルに組み込み、各段階を検証した |
| F2 | Why a Bayesian method? | The data are sparse, so I need to know which interactions the data determine and which they leave open. ／ データが少ないので、何が決まり何が決まらないかを事後分布で示す |
| F3 | How many interactions are determined? | Five, six, nine and eight of the fifteen in CS, CH, DS and DH; in the commensal conditions only the early colonisers. ／ 15 個のうち CS 5、CH 6、DS 9、DH 8 個。常在では初期定着菌の間だけ |
| F4 | Why is A symmetric? | Only the symmetric part enters the free energy, so directed effects cannot be represented; a limitation. ／ 自由エネルギーには対称部分しか入らない。一方向の作用は表せない（限界） |
| F5 | How do you know the ANSYS implementation is correct? | It reproduces the closed-form solution, and the growth field agrees with a finite-volume reference within about one percent. ／ 閉じた形の解を再現し、有限体積法の基準解と 1 % 程度で一致 |
| F6 | Why did you run to T* = 5? | Near the limit the unstable mode starts very small; short runs look smooth even when the step is too large. ／ 限界近くでは不安定モードが小さく、短い計算では見えない |
| F7 | Is the pH a validation? | It was not used in the inference, so it is a check of the predicted composition, not a full validation. ／ 推定に使っていないので、組成の予測の確認。完全な検証ではない |
| F8 | What is the biggest limitation? | No measured stress of a growing oral biofilm; the stresses are verified, not validated. ／ 応力の実測がないので、検証はできても妥当性確認はできていない |
| F9 | What comes next? | At Keio: stress feeding back into growth, on tooth and implant geometries. ／ 慶應で、応力から成長への帰還を入れ、歯とインプラントで計算 |
| F10 | Why does the stress in a biofilm matter? | Stress decides where the biofilm deforms or detaches, and around implants it acts on the tissue. ／ 応力は、バイオフィルムがどこで変形したり剥がれたりするかを決め、インプラントの周りでは組織に作用するから |
| F11 | What is the extended Hamilton principle here? | A variational principle: the growth and the composition follow from making the sum of energy rate, dissipation and constraints stationary. ／ 変分原理。エネルギーの変化率、散逸、拘束の和を停留させると、成長と組成の式が出る |
| F12 | What is TMCMC? | A sampler that moves particles from the prior to the posterior through a series of tempered distributions. ／ 粒子を、事前分布から段階的に温度を下げた分布を通して事後分布へ移すサンプラー |
| F13 | What does alpha mean? | The growth stretch: the growth part of the deformation is alpha times the identity, and alpha starts at one. ／ 成長による伸び。成長の変形が alpha 倍の単位テンソルで、alpha は 1 から始まる |
| F14 | What is the IKM element? | An ANSYS element of the institute that solves the biofilm field with a nodal integration method; I coupled the growth law to it. ／ 研究所の ANSYS の要素で、バイオフィルムの場を節点積分法で解く。そこに成長則をつないだ |
| F15 | Where does the point model enter the element? | At every Gauss point, next to the field, so the composition is known at each integration point. ／ 各 Gauss 点で場と並べて解くので、各点で組成が分かる |
| F16 | Why beta = 0.02 square millimetres per unit time? | It is the value of Klempt et al. 2024 converted to the two-millimetre cube; the conversion is my assumption, so I show the sensitivity. ／ Klempt 2024 の値を 2 mm の立方体に換算したもの。換算は仮定なので感度も示した |
| F17 | Why is the stiffness so low, ten pascals? | It is the value of Klempt et al. 2024; in this setting the stress is proportional to E, so ratios between conditions do not depend on it. ／ Klempt 2024 の値。この設定では応力は E に比例するので、条件間の比は E によらない |
| F18 | Why ANSYS and not your own code? | The element of the institute runs in ANSYS, and the user material lets me add the growth law without changing the solver. ／ 研究所の要素が ANSYS で動いていて、ユーザー材料ならソルバーを変えずに成長則を足せるから |
| F19 | What does "identified" mean? | The posterior standard deviation is below 0.8 times that of a uniform distribution on the prior box, so the data narrow it. ／ 事後の標準偏差が、事前の箱の一様分布の 0.8 倍より小さい。データで狭まっている |
| F20 | Why are so few entries identified in the commensal conditions? | Fn and Pg stay below one percent there, so the data carry little information on their interactions. ／ 常在では Fn と Pg が 1 % 未満なので、その相互作用の情報がデータにほとんどない |
| F21 | Why four stages and three seeds? | Stages make the tempering cheaper and more stable; three seeds show whether a run has really converged. ／ 段に分けると計算が安く安定する。3 seed で本当に収束したかが分かる |
| F22 | What does the Fn knockout show? | In the model, removing Fn removes the late growth of Pg in the dysbiotic flow condition; a co-culture without Fn could test it. ／ モデルで Fn を除くと、dysbiotic の流れ条件で Pg の終盤の増加が消える。Fn 抜きの共培養で確かめられる |
| F23 | Can you compare commensal and dysbiotic directly? | Not fully: the two communities also differ in the Veillonella species and the Pg strain. ／ 完全にはできない。Veillonella の種と Pg の株も違うから |
| F24 | What is lambda in the stability analysis? | beta times the time step over the mesh size squared; above a limit the explicit update becomes unstable. ／ β × 時間刻み ÷ メッシュ幅の 2 乗。限界を超えると陽的な更新が不安定になる |
| F25 | Why does the limit depend on the neighbours? | The thirtieth neighbour lies in a shell of seven equal points, and which five are taken makes the stencil slightly asymmetric. ／ 30 番目の近傍は同じ距離の 7 点の殻にあり、どの 5 点を取るかで少し非対称になる |
| F26 | Is the coupling one-way or two-way? | Two-way through the nutrient: the fields give the point model its nutrient and antibiotic, and only its living biomass consumes them; the nutrient then steers the spreading. The stress does not act back. ／ 栄養を通して双方向。場が点モデルに栄養と抗菌薬を渡し、生きている菌だけがそれを消費する。栄養が広がりを決める。応力は戻らない（20 Oct 更新） |
| F27 | What was the hardest part? | Getting the convergence of the inference and of the element right; several runs that looked fine failed the checks. ／ 推定と要素の収束をきちんと確かめること。よく見えた計算がいくつも判定に落ちた |
| F28 | Why do you not scale the point model to the field? | The point model is formulated on a material point, and its authors say its void is an auxiliary variable, not a physical void. So I use it as published, with its own constraint, and give it only the nutrient and the antibiotic. ／ 点モデルは物質点のモデルで、空隙は補助変数（物理的な空隙ではない）と著者が書いている。だから発表された形のまま使い、栄養と抗菌薬だけを渡す |
| F29 | How do you link the time of the two models? | I choose T* = 1 as the duration of the point model's data. Then the point-model time per increment follows; it is not fitted. With synthetic data the coupled point model then reproduces the point model alone to 7e-7. ／ T* = 1 を点モデルのデータの期間にする。すると 1 増分あたりの点モデルの時間が決まり、推定しない。合成データで単体と 7e−7 で一致 |
| F30 | How does the point model act back on the field? | Only its living biomass, the sum of phi_i psi_i, consumes the nutrient. The nutrient then steers the spreading through the front term, and so the growth and the stress. ／ 生きている菌（Σφ_iψ_i）だけが栄養を消費する。栄養が前線の項で広がりを決め、成長と応力に届く |
| F31 | How do the nutrient values reach the point model? | In the neighboured element method the fields live at the integration points, so the material routine reads the value at its own Gauss point; there is no interpolation from nodes. Inside an increment I interpolate linearly in time between the previous and the current value. ／ NEM では場が積分点の値なので、材料ルーチンは自分の Gauss 点の値をそのまま読む（節点からの補間はない）。増分の中では前と今の値のあいだを時間で線形補間 |
| F32 | Where does the antibiotic come from? | The element solves a second nutrient for a second species, which I do not use. I use that field for the antibiotic, with the same diffusion and uptake by the living biomass. The calibration data contain no antibiotic, so I show it with synthetic data. ／ 要素の、使っていない 2 つめの栄養の場を抗菌薬に使う。拡散と生きた菌による取り込みは同じ形。較正データに抗菌薬はないので合成データで示す |
| F33 | How did you check the coupled point model in ANSYS? | I recomputed every step in Python from the same start state and the same inputs; they agree to 2e-12 in every step, once the Newton iteration is converged (up to 50 iterations at the first step). ／ 同じ初期状態と入力で全ステップを Python で計算し直し、Newton を収束させれば（最初のステップは最大 50 回）全ステップ 2e−12 で一致 |

## G. 覚えておく基本の数値（10月13日）

| # | 何の値 | 値（EN ／ JA） |
|---|---|---|
| G1 | Species and conditions | Five species, four conditions (CS, CH, DS, DH), days 1 to 21 ／ 5 菌種、4 条件、1〜21 日 |
| G2 | Interaction unknowns | Fifteen entries of a symmetric five-by-five matrix ／ 対称 5×5 行列の 15 成分 |
| G3 | Identified entries | 5, 6, 9, 8 of 15 (CS, CH, DS, DH; DH from stage iii) ／ 15 個中 5、6、9、8 個（DH は段階 (iii)） |
| G4 | Nutrient scale | c* = 25, fixed because only c* times A enters ／ c* = 25。c*×A の積でしか入らないので固定 |
| G5 | GPU speedup | 170 times per step, 147 times end to end ／ 1 段あたり 170 倍、全体で 147 倍 |
| G6 | Stiffness (Klempt 2024) | E = 10 Pa, Poisson ratio 0.49, shear modulus 3.3557 Pa ／ E = 10 Pa、ν = 0.49、μ = 3.3557 Pa |
| G7 | Diffusion of phi | beta = 0.02 mm² per unit time (converted, an assumption) ／ β = 0.02 mm²/T*（換算した仮定） |
| G8 | Test geometry | Two-millimetre cube, seed of 32 elements on the 8³ mesh ／ 2 mm の立方体、8³ でシード 32 要素 |
| G9 | Meshes | 8³, 16³, 24³ elements ／ 8³、16³、24³ |
| G10 | Growth error of the element | Seed mean +4.1, +1.3, +0.7 % on 8³, 16³, 24³ ／ シード平均 +4.1、+1.3、+0.7 % |
| G11 | Closed-form check | Constrained cube, alpha − 1 = 0.05: stress −1.0193e−4 MPa in all three directions ／ 拘束した立方体、α − 1 = 0.05 で応力 −1.0193e−4 MPa |
| G12 | Stability limit (lambda) | 0.062 on 8³; 0.12 to 0.14 on 16³; 0.12 to 0.16 on 24³ (30 neighbours) ／ 8³ で 0.062、16³ で 0.12〜0.14、24³ で 0.12〜0.16 |
| G13 | Runs to T* = 5 | 24³: smooth at 0.130, stops at 0.158 ／ 24³：0.130 で滑らか、0.158 で停止 |
| G14 | Five-species growth in the element | Same mean growth as Eq. 36 by choice of the rate; seed von Mises stress 8 % lower ／ 成長の平均は Eq. 36 にそろえた。シードの von Mises 応力は 8 % 低い |
| G15 | Colloquium | 16 Dec 2026, 13:00, IKM Room 320; 20 min talk, 10 to 20 min questions ／ 12 月 16 日 13:00、320 号室、発表 20 分、質疑 10〜20 分 |

## H. 数値ごとに聞かれそうな追加の質問（10月13日）

| # | 数値 | 予想される質問 | 答え（EN ／ JA） |
|---|---|---|---|
| H1 | 5 species | Why these five? | They are the defined consortium of Heine et al., from early colonisers to the pathogen Pg. ／ Heine らが決めた菌群で、初期定着菌から病原菌 Pg まで並ぶ |
| H2 | 15 unknowns | With 25 data points, is that not too many? | Yes, which is why I report which entries the data determine instead of a single best fit. ／ 多い。だから最適値一つではなく、どの成分が決まるかを報告した |
| H3 | 5, 6, 9, 8 | What about the entries that are not identified? | I do not interpret them; they stay at the prior. ／ 解釈しない。事前分布のまま |
| H4 | c* = 25 | Does the choice of c* change the results? | No, only c* times A enters, so another c* rescales A. ／ 変わらない。c*×A でしか入らないので、A の大きさが変わるだけ |
| H5 | 147 times | Why not use a gradient-based sampler? | I did not compare them; I used TMCMC because it runs all particles in parallel on the GPU. ／ 比較していない。粒子を GPU で並列に回せる TMCMC を使った |
| H6 | E = 10 Pa | Real biofilms are stiffer, up to kPa. | Yes; in this setting the stress scales with E, so the pattern holds and only the magnitude changes. ／ その通り。この設定では応力は E に比例するので、分布は同じで大きさだけ変わる |
| H7 | beta = 0.02 | What happens with a much smaller beta? | The diffusion length falls below the mesh size and the seed stress does not converge. ／ 拡散の長さがメッシュより小さくなり、応力が収束しない |
| H8 | 2 mm cube | Is a cube realistic? | No, it is a test geometry for verification; teeth and implants come at Keio. ／ 現実的ではない。検証用の形で、歯とインプラントは慶應で |
| H9 | 8³ to 24³ | Why not finer meshes? | The NEM element is expensive; for the reference I used a separate code up to 192³. ／ 要素の計算が重い。基準解は別のコードで 192³ まで解いた |
| H10 | +0.7 % | Is that converged? | It falls with refinement, about first order for the stress; 24³ is within a few percent of the limit. ／ 細かくすると減る（応力はほぼ 1 次）。24³ で極限から数 % 以内 |
| H11 | −1.0193e−4 MPa | Why is the stress negative? | The cube is constrained, so the growth is pushed back and the stress is compressive. ／ 拘束されているので、成長が押し戻されて圧縮になる |
| H12 | lambda limits | Which time step did you use? | Below the limit with a margin, lambda at most 0.128 on 16³ and 0.072 on 24³. ／ 余裕を持って限界より下（16³ で 0.128 以下、24³ で 0.072 以下） |
| H13 | 0.130 / 0.158 | Why not an implicit update? | It is the element's existing explicit scheme; an implicit version is future work. ／ 要素にもともとある陽的な方法だから。陰的な版は今後の課題 |
| H14 | −8 % | Is that difference significant? | It is one mesh and one rate chosen by me, so it shows the effect, not a precise value. ／ メッシュ 1 つ、速さは私が決めた値なので、効果の有無を示すもので、精密な値ではない |
| H15 | 20 min | (準備) | Main slides 16, appendix for questions. ／ 本編 16 枚、質問には付録で答える |

## I. DH が段階 (iii) の理由（10月19日、日本語と英語）

要点：最後の段階で、3 つの seed が違うモードに分かれた。当てはまりの良さは変わらない。データが決めきれないことを、制約として書いた。

| # | 予想される質問 | 答え（EN ／ JA） |
|---|---|---|
| I1 | Why is DH reported from stage (iii), not the final stage? | At stage (iv) the three seeds settled in different posterior modes, so criterion (e) failed: the medians of a22 differed by 0.79 posterior standard deviations, and the limit is 0.5. Stage (iii) passes all criteria, so I report that one. ／ 段階 (iv) では 3 つの seed が別々のモードに落ち着き、基準 (e) を満たさなかった（a22 の中央値が事後標準偏差の 0.79 倍ずれた。上限は 0.5 倍）。段階 (iii) は全部の基準を満たすので、そちらを使った |
| I2 | Does that make the DH results unreliable? | The fit is good: the viable composition has an RMSE of 0.073. What the split shows is that the data do not fix the size of the late Pg increase, and I do not interpret that part. ／ 当てはまりは良い（生菌の組成の RMSE 0.073）。分かれたことから分かるのは、後半の Pg の増え方の大きさをデータが決められないということ。そこは解釈していない |
| I3 | Why did the final stage split? | The final stage uses the most data. Two sets of coefficients explain the DH data about equally well, and the seeds found different ones. ／ 最後の段階は使うデータが一番多い。DH のデータをほぼ同じくらい説明できる係数の組が 2 つあり、seed ごとに違う方を見つけた |
| I4 | Would more particles or a longer run fix it? | More particles would sample both modes better, but they cannot decide between them. That needs more measurements in the late phase. ／ 粒子を増やせば両方のモードをよく拾えるが、どちらかには決められない。決めるには後半の測定点がもっと必要 |
| I5 | Does this affect the stress results? | No. The finite element runs of Chapters 4 and 5 use the CH posterior, which passes at the final stage. ／ 影響しない。第 4・5 章の有限要素の計算は CH の結果を使っていて、CH は最後の段階まで基準を満たしている |
| I6 | Which DH results can I trust? | The entries reported as identified: eight in DH, for example the self-limitation of Vei and Fn. ／ 同定できたと書いた成分（DH では 8 つ。例えば Vei と Fn の自己抑制） |

## 記録

| 日付 | 出した問題 | 答えられなかった問題 |
|---|---|---|
| 10月11日 | B1, C1, C4, C5, D3（言い換え）／ A1〜A5 | 第 1 回（数字の多い版）は難しすぎた。一言の答えに切り替え |
| 10月13日 | A6, B5, C2, D2, E1 | （答えなし、次に再出題） |
| 10月13日 | F1, F3, F5, F6, F8（英語で） | （答え待ち） |
| 10月13日 | Q&A：F10〜F18 を表示 | （読む用） |
| 10月13日 | Q&A：F19〜F27 を表示 | （読む用） |
| 10月13日 | 基本の数値 G1〜G15 を表示 | （読む用） |
| 10月13日 | 数値ごとの予想質問 H1〜H15 を表示 | （読む用） |
| 10月19日 | Q&A：DH が段階 (iii) の理由 I1〜I6 を表示 | （読む用） |
| 10月19日 | A6, C2, D2（再出題）, B2, I1 | 答えを表示（次回また出す） |
| 10月19日 | B5, E1（再出題）, C3, D4, I3 | 答えを表示（次回また出す） |
| 10月19日 | A7, B6, C6, D5, E2 | （答え待ち） |
| 10月21日 | A8, B3, C8, D6, A7（再出題） | （答え待ち） |
