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

答えは英語で言えるように練習する。日本語は意味の確認用。

| # | 問題 | 答え（EN ／ JA） |
|---|---|---|
| F1 | What is the main contribution? | I brought a calibrated multispecies point model into a finite element growth model, and verified each step. ／ 較正した多菌種の点モデルを有限要素の成長モデルに組み込み、各段階を検証した |
| F2 | Why a Bayesian method? | The data are sparse, so I need to know which interactions the data determine and which they leave open. ／ データが少ないので、何が決まり何が決まらないかを事後分布で示す |
| F3 | How many interactions are determined? | Four, six, eight and nine of the fifteen; in the commensal conditions only the early colonisers. ／ 15 個のうち 4、6、8、9 個。常在では初期定着菌の間だけ |
| F4 | Why is A symmetric? | Only the symmetric part enters the free energy, so directed effects cannot be represented; a limitation. ／ 自由エネルギーには対称部分しか入らない。一方向の作用は表せない（限界） |
| F5 | How do you know the ANSYS implementation is correct? | It reproduces the closed-form solution, and the growth field agrees with a finite-volume reference within about one percent. ／ 閉じた形の解を再現し、有限体積法の基準解と 1 % 程度で一致 |
| F6 | Why did you run to T* = 5? | Near the limit the unstable mode starts very small; short runs look smooth even when the step is too large. ／ 限界近くでは不安定モードが小さく、短い計算では見えない |
| F7 | Is the pH a validation? | It was not used in the inference, so it is a check of the predicted composition, not a full validation. ／ 推定に使っていないので、組成の予測の確認。完全な検証ではない |
| F8 | What is the biggest limitation? | No measured stress of a growing oral biofilm; the stresses are verified, not validated. ／ 応力の実測がないので、検証はできても妥当性確認はできていない |
| F9 | What comes next? | At Keio: stress feeding back into growth, on tooth and implant geometries. ／ 慶應で、応力から成長への帰還を入れ、歯とインプラントで計算 |

## 記録

| 日付 | 出した問題 | 答えられなかった問題 |
|---|---|---|
| 10月11日 | B1, C1, C4, C5, D3（言い換え）／ A1〜A5 | 第 1 回（数字の多い版）は難しすぎた。一言の答えに切り替え |
| 10月13日 | A6, B5, C2, D2, E1 | （答え待ち） |
