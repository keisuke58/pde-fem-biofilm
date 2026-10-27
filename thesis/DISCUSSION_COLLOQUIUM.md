# Discussion of the main results for the colloquium (16 Dec 2026)

Purpose: the "why" behind the main numbers of the thesis, for the 20-minute talk and the questions
(IKM Room 320, Prof. Junker and the IKM committee). Written 27 Oct 2026, after the decision to add no
new results. Every number is taken from the thesis (`thesis/chapters/ch4_ansys.tex`, `ch3_heine.tex`,
`ch6_conclusion.tex`); numbers that are only in the paper draft (`papers/coupling_paper/manuscript/
coupling_paper.tex`, Section 2.4 and 4) are marked **[paper draft]**. Statements derived here by algebra
are marked **[quick check]**, readings that go beyond what is computed **[interpretation]**.

Each item: (1) the observation with its number, (2) the mechanism from the equations, (3) a one-minute
spoken explanation (EN, then JA), (4) the likely follow-up question with a short answer. The follow-up
questions are also in `QUIZ_COLLOQUIUM.ja.md`, section K.

## 0. The equations the items refer to

Field (Klempt et al. 2024, Eqs. 34 to 36, consumption by the living biomass as agreed, 20 Oct 2026):

$$\dot\phi_{3D}=\beta\,\nabla^2\phi_{3D}+k_\alpha\,\alpha-r\,\frac{c}{k+c}\,\mathbf n_c\cdot\nabla\phi_{3D},\qquad
\mathbf n_c=\frac{\nabla c}{|\nabla c|}\quad\text{(Eq. 34)}$$

$$\dot c=d\,\nabla^2c-g\,\phi_{3D}\Bigl(\sum_i\bar\phi_i\Bigr)c,\qquad \bar\phi_i=\phi_i\psi_i\quad\text{(Eq. 35)}$$

$$\dot\alpha=k_\alpha\,\phi_{3D}\quad\text{(Eq. 36)},\qquad \mathbf F_g=\alpha\mathbf I,\quad E(\phi_{3D})=E\,(\phi_{3D}^2+f).$$

Point model (Klempt et al. 2026, thesis Eq. ch5-point), second line, for the living fraction:

$$\eta_i\,(\phi_i\psi_i\,\dot\phi_i+\phi_i^2\,\dot\psi_i)=c^*\phi_i\sum_jA_{ij}\phi_j\psi_j-\alpha^*b_i\,\psi_i-(\text{barrier}).$$

Constants of the two-species examples (thesis §4.7.1, Klempt 2024 Table 2 converted to the 2 mm cube,
Fritsch et al. 2025 case I): $k_\alpha=10^{-3}\,T^{*-1}$, $\beta=0.02$ mm$^2/T^*$, $r=10$ mm$/T^*$, $k=1$,
$g\ell^2/d=1$ ($\ell=1$ mm), $E=10$ Pa, $\nu=0.49$, $f=10^{-3}$; $\mathbf A=\begin{pmatrix}1&0.1\\0.1&1\end{pmatrix}$,
$b=(1,2)$, $\eta=(1,2)$, $c^*_0=\alpha^*_0=100$, antibiotic from $T^*=1/3$, $s=0.1$.

The path from the composition to the stress, which most items use (paper draft Eq. path; the same
statement is in thesis §4.6 "From the point model to the field"):

$$\sum_i\bar\phi_i\ \xrightarrow{\ \text{Eq. 35}\ }\ c\ \xrightarrow{\ -r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}\ \text{of Eq. 34}\ }\ \phi_{3D}\ \xrightarrow{\ \text{Eq. 36}\ }\ \alpha\ \longrightarrow\ \boldsymbol\sigma .$$

---

## 1. Why the composition hardly changes the growth

**(1) Observation.** Table 4 (grid convergence, `tab:ch5-gci`): seed mean of $\alpha-1$ at $T^*=1$,
extrapolated, $9.16\times10^{-4}$ without and $9.22\times10^{-4}$ with the antibiotic in the separate
example (0.7 %), $8.56$ against $8.54\times10^{-4}$ in the mixed one (0.2 %), while the living volume
fraction of the seed falls from $0.965$ to $0.152$ and from $0.961$ to $0.038$ (`tab:ch5-two-species`).
Conclusion: "the growth changes by less than 2 %". In the earlier one-way five-species run (§4.7.3) the
seed mean of $\alpha-1$ is $6.752\times10^{-4}$ to all printed digits for the CH MAP, for equal shares
and for the two-species runs.

**(2) Mechanism.**
- Eq. 36 contains $\phi_{3D}$ only. Neither $\phi_i$, $\psi_i$ nor $c$ appears in it. The point-model
  state enters the field equations at one place only, the factor $\sum_i\bar\phi_i$ in Eq. 35, and the
  nutrient enters Eq. 34 only in the term $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$. So the
  composition reaches $\alpha$ only through the position of the edge of the biofilm.
- Bound **[quick check]**: with $0\le\phi_{3D}\le1$, integrating Eq. 36 gives
  $\alpha(T)-1=k_\alpha\int_0^T\phi_{3D}\,dt\le k_\alpha T$, so at $T^*=1$ no composition can raise
  $\alpha-1$ above $10^{-3}$. The seed holds $\phi_{3D}=1$ from the start and reaches $9.16\times10^{-4}$,
  92 % of the bound; the remaining 8 % is the part of the seed mean where the edge has $\phi_{3D}<1$.
  A change of the composition can therefore raise the seed growth by at most about 9 % and lower it
  only through the edge.
- The antibiotic lowers $\psi_i$ (the living fraction), not $\phi_{3D}$: dead biomass stays in the
  field and keeps growing at the same rate (a limitation, §6.2). Hence "about 1 %".
- In the one-way five-species run the path is cut at both ends (no front term, consumption not taken
  from the point model), so the digits are identical; that run demonstrates that the calibrated model
  runs in the element, not an effect on the stress.

**(3) Spoken, EN.** "The growth law is the third equation of Klempt et al. 2024, d alpha by dt equals
k alpha times phi 3D. It counts all biomass, living or dead, and the composition does not appear in it.
The point model reaches the field only through the consumption: only living biomass consumes nutrient.
The nutrient then moves the edge of the biofilm through the front term, minus r c over k plus c times
n c dot grad phi 3D. That is the only path from the composition to the growth. Inside the seed phi 3D is
one from the start, so alpha minus one is bounded by k alpha times T star, ten to the minus three at
T star one, and the seed is already at 92 percent of that bound. The antibiotic removes most of the
living biomass, from 0.97 to 0.15, but it lowers psi, not phi 3D, so the growth changes by about one
percent. This is a property of the published equations, not of my numerics, and it tells us what a model
needs if the composition is to change the mechanics: a growth rate that counts the living biomass."

**(3) Spoken, JA.** 「成長則は Klempt 2024 の第 3 式、$\dot\alpha=k_\alpha\phi_{3D}$ です。生きていても死んでいても
全部のバイオマスを数え、組成は式に出てきません。点モデルが場に届くのは消費の項だけで、生きている菌だけが
栄養を消費します。栄養は式 34 の項 $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$ を通してバイオフィルムの縁を動かします。
組成から成長への道はこれ一本です。シードの中では $\phi_{3D}$ が最初から 1 なので、$\alpha-1$ は $k_\alpha T^*$、
$T^*=1$ で $10^{-3}$ を超えられず、シードはすでにその 92 % にあります。抗菌薬は生きている量を 0.97 から 0.15 に
下げますが、下げるのは $\psi_i$ で $\phi_{3D}$ ではないので、成長は約 1 % しか変わりません。これは私の数値計算の
性質ではなく、発表された式の性質です。組成が力学を変えるモデルには、生きている量で成長する成長則が要る、
ということを示しています。」

**(4) Follow-up.** *"Then the point model does not matter for the stress at all?"* With the published
equations it matters only through the front speed $r\,c/(k+c)$, which lies between 0 and $r/2$, so the
effect is bounded. That is the result: the coupling is the base on which a second path (growth or front
speed from $\sum_i\bar\phi_i$, species-dependent stiffness) can be added and measured. An estimate from the
time series of the seed **[paper draft, estimate]**: with $\dot\alpha=k_\alpha\phi_{3D}\sum_i\bar\phi_i$
the seed growth would fall by about 40 % (separate) and 65 % (mixed) under the antibiotic instead of 1 %.
The thesis keeps that only in the outlook (user's decision, 15 Oct).

---

## 2. Why the antibiotic's effect on the converged seed stress is within the grid uncertainty

**(1) Observation.** Separate example, seed von Mises stress extrapolated to the grid limit
(`tab:ch5-gci`): $6.55\times10^{-4}$ Pa without and $6.70\times10^{-4}$ Pa with the antibiotic, GCI 0.7 %
and 8.3 %. The 8.3 % band of the antibiotic value ($\pm0.56\times10^{-4}$ Pa) contains the value without.
On one grid the difference is about 10 %: $6.73$ against $6.00\times10^{-4}$ Pa on $64^3$ at the step
$0.0025$, $6.56$ against $5.68\times10^{-4}$ Pa at the step $0.025$ (`tab:ch5-two-species`). The two
series run in opposite directions: $7.44,\ 6.73,\ 6.585$ (falling, $p=2.3$) without, $5.52,\ 6.00,\ 6.28$
(rising, $p=0.75$) with the antibiotic.

**(2) Mechanism.**
- Physical side: the antibiotic lowers $\sum_i\bar\phi_i$ from 0.97 to 0.15, so less nutrient is consumed
  and $c$ in the seed rises from $0.54$ to $0.91$ (Table `tab:ch5-two-species`). The front speed
  $r\,c/(k+c)$ with $k=1$ rises from $0.35\,r$ to $0.48\,r$, by 37 % **[quick check]**. That moves the
  edge of the biofilm faster, but the stress of the seed is set by the jump of $\alpha$ across the edge of
  the seed, and $\alpha$ inside the seed is $k_\alpha\int\phi_{3D}dt$ with $\phi_{3D}=1$ in both cases.
  The height of the jump hardly changes; the edge moves. **[interpretation]** The first layer around the
  seed is crossed by the front within about $h/(r\,c/(k+c))\approx0.125/3.5\approx0.04\,T^*$ on $16^3$
  **[quick check]**, long before the antibiotic starts at $T^*=1/3$, so the faster front after the onset
  acts far from the seed.
- Numerical side: the two cases have different discretisation errors. Without antibiotic the coarse grids
  overestimate the stress (series falling), with antibiotic they underestimate it (series rising); the
  signs of the series are computed facts, the reason is **[interpretation]**: after the onset the edge
  moves faster and the second-order upwind scheme smears a faster edge more on a coarse grid, which lowers
  the jump and the stress, so the stress rises with refinement. A comparison on one grid therefore
  compares two different numerical errors, and it can show 10 % with either sign. The yardstick is the
  grid convergence index of Celik et al. 2008 (three grids, observed order, extrapolated value, safety
  factor 1.25): the difference of 2 % lies inside the 8.3 % band.
- The same holds for the time step: at $\Delta t=0.025$ the separate stress is about 6 % below the limit
  and the mixed one 15 % above (§4.7.1 "Time step"), while $\alpha-1$ changes by less than 0.3 %.

**(3) Spoken, EN.** "The antibiotic removes most of the living biomass, so less nutrient is consumed and
the nutrient in the seed rises from 0.54 to 0.91. Through the front term, minus r c over k plus c times
n c dot grad phi 3D, the edge of the biofilm moves faster, by 37 percent. But the stress of the seed comes
from the jump of alpha across the edge of the seed, and alpha inside the seed is the same with and without
antibiotic, because phi 3D is one there. So the edge moves, but the jump keeps its height. Numerically,
on a single grid I see 10 percent, but the two cases converge from opposite sides: without the antibiotic
the stress falls with refinement, with it the stress rises. A comparison on one grid compares two
different discretisation errors. After extrapolation with the grid convergence index the two values are
6.55 and 6.70 times ten to the minus four pascal, with an uncertainty of 8 percent on the second, so the
difference is within the uncertainty. I report the band, not the single-grid number."

**(3) Spoken, JA.** 「抗菌薬は生きている菌の大半を取り除くので、消費される栄養が減り、シードの栄養は 0.54 から 0.91
に上がります。式 34 の項 $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$ を通して、縁の速さは 37 % 上がります。しかし
シードの応力はシードの縁での $\alpha$ の段差で決まり、シードの中の $\alpha$ は抗菌薬のあるなしで同じです（そこでは
$\phi_{3D}=1$ だから）。縁は動きますが、段差の高さは変わりません。数値的には、一つの格子で見ると 10 % の差が出ますが、
二つの場合は反対側から収束します。抗菌薬なしは細かくすると下がり、ありは上がる。一つの格子での比較は、二つの
違う離散化誤差を比べていることになります。格子収束指標で外挿すると $6.55$ と $6.70\times10^{-4}$ Pa、後者の不確かさが
8 % なので、差は不確かさの中です。私は一格子の数字ではなく幅を報告します。」

**(4) Follow-up.** *"Why trust an extrapolation with an observed order of 0.75?"* The three grids are
in the asymptotic range (the ratio of the differences is consistent with a single order), the index is
the standard of Celik et al. 2008 for exactly this situation, and I quote the band, 8.3 %, as the
uncertainty, not the extrapolated value alone. Where the grids are not in the asymptotic range (mixed with
the Fritsch start fractions) I give no extrapolation. The sign of a one-grid difference is not robust;
the band is.

---

## 3. Why the mixed seed carries about half the stress of the separate one

**(1) Observation.** Seed von Mises stress at $T^*=1$: $3.50$ against $6.56\times10^{-4}$ Pa ($64^3$,
step 0.025, `tab:ch5-two-species`), ratio 0.53; extrapolated $4.50$ against $6.55\times10^{-4}$ Pa,
ratio 0.69 (`tab:ch5-gci`). The mixed growth is 7 % lower ($8.56$ against $9.16\times10^{-4}$). The
largest $\phi_{3D}$ in the seed at the end rises with the grid, $0.60,\ 0.69,\ 0.73$ on $32^3$, $64^3$,
$128^3$; the thesis: "the front carries intermediate values of $\phi_{3D}$ at a nearly discontinuous
edge", and the stress "is set by the jump of $\alpha$ across the edge of the seed, which the front
sharpens".

**(2) Mechanism.**
- The stress is a jump quantity: inside the seed $\alpha-1=k_\alpha\int\phi_{3D}dt$ is large, outside
  it is small, and the soft surroundings ($E\phi_{3D}^2$) resist the swelling of the seed only at that
  jump. A lower or wider jump gives a lower von Mises stress.
- In both examples the field seed starts at $\phi_{3D}=1$ (the 0.55 of the mixed example is the point
  model's $\phi_1+\phi_2$; see "Points to check" below). In the separate example the maximum of
  $\phi_{3D}$ stays at 1 to the end on every grid; in the mixed example it falls below 1 from about
  $T^*=0.4$ and ends at $0.60/0.69/0.73$ (recorded maxima of the result files, `results/2026-10-3d`,
  the numbers quoted in the thesis). So in the mixed example the block of $\phi_{3D}=1$ is diluted: the
  term $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$ transports $\phi_{3D}$ along $\mathbf n_c$, and
  where the direction field converges the block is compressed below the resolution of the grid and the
  scheme smears it; the maximum then rises with the grid **[interpretation]**. With a lower $\phi_{3D}$ in
  the seed, (a) $\alpha-1$ is lower (−7 %, so the jump is lower), (b) the edge is wide, and (c) the
  stiffness $E\phi_{3D}^2$ of the seed is lower ($0.73^2=0.53$ **[quick check]**). All three lower the
  stress; the thesis names (b).
- Why the block is kept in one example and diluted in the other is not analysed in the thesis. The only
  difference between the two examples in the field equations is where species B consumes the nutrient
  (the four edge columns in the separate example), which changes $\nabla c$ and so $\mathbf n_c$ around
  the seed **[interpretation]**.

**(3) Spoken, EN.** "The seed stress comes from the jump of alpha across the edge of the seed. In the
separate example the seed keeps phi 3D equal to one to the end, so the jump is sharp. In the mixed example
the front term, minus r c over k plus c times n c dot grad phi 3D, carries the biofilm away from the seed
in a way that leaves intermediate values: the largest phi 3D in the seed is only 0.6 to 0.73 at the end,
and it still rises with the grid. A lower phi 3D means less growth, a wider edge and, through the
stiffness E times phi squared, a softer seed. All three lower the stress, so the mixed seed carries about
half, 3.5 against 6.6 times ten to the minus four pascal on the same grid. The same nearly discontinuous
edge is why the mixed stress converges slowly, with an observed order of 0.75 and a grid index of 26 to
29 percent, whereas the growth converges within 3 percent in every example."

**(3) Spoken, JA.** 「シードの応力はシードの縁での $\alpha$ の段差から来ます。分けて置く例ではシードの $\phi_{3D}$ が
最後まで 1 なので段差が鋭い。混ぜて置く例では、項 $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$ がバイオフィルムを
シードから運び出し、中間の値を残します。シード内の最大の $\phi_{3D}$ は終わりに 0.6 から 0.73 で、格子とともにまだ
上がります。$\phi_{3D}$ が低いと成長が少なく、縁が広く、剛性 $E\phi_{3D}^2$ でシードが柔らかい。三つとも応力を下げるので、
混ぜた例は同じ格子で $3.5$ 対 $6.6\times10^{-4}$ Pa、ほぼ半分です。同じ不連続に近い縁が、混ぜた例の応力の収束が遅い
理由でもあります（次数 0.75、GCI 26〜29 %）。成長量はどの例でも 3 % 以内で収束します。」

**(4) Follow-up.** *"Is half the stress a physical result or a numerical one?"* The direction is
physical: a seed whose edge carries intermediate $\phi_{3D}$ loads itself less. The size is uncertain:
the extrapolated ratio is 0.69 with a 29 % band on the mixed value, so I say "about half on the same grid,
between a half and two thirds after extrapolation".

---

## 4. Why the stress converges slowly but the growth fast

**(1) Observation.** Growth: GCI 0.7 to 3.3 % in all examples, observed orders 1.2 to 1.6. Stress: 0.7 %
(separate, no antibiotic, $p=2.3$), 8.3 % (with antibiotic, $p=0.75$), 26 to 29 % (mixed, $p\approx0.75$).
Cell Péclet number $r h/(2\beta)$ about 31 on $16^3$ (§4.7.1). With a constant added diffusion $r h/4$ the
element's seed stress on $16^3$ is about a quarter of the converged value; with upwind gradients it is
$2.14,\ 4.56,\ 5.15\times10^{-4}$ Pa on $8^3$, $16^3$, $24^3$, extrapolated $5.8\times10^{-4}$ Pa (GCI
16.5 %), inside the band of the second-order schemes. In the one-species seed problem the growth converges
with an observed order of 1.7 in $h$ and the stress "only to about first order, because the stiffness jumps
at the surface of the seed" (§4.7).

**(2) Mechanism.**
- Growth is an integral quantity: the seed mean of $\alpha-1$ is $k_\alpha$ times the time integral of
  the seed mean of $\phi_{3D}$. A scheme that smears the edge over a few cells moves a little $\phi_{3D}$
  across the edge of the seed but changes the integral by a small amount (the explicit update conserves the
  mass of $\phi$ up to the source; the domain integral of $\alpha-1$ is within 0.3 % on every mesh in the
  one-species check).
- Stress is a derivative quantity: it is set by the jump of $\alpha$, and of the stiffness
  $E\phi_{3D}^2$, across one layer of elements at the edge of the seed. A smeared edge lowers the jump
  directly. A derivative loses one order of accuracy against the integral, which is what the observed
  orders show (0.75 against 1.2 to 1.6; about 1 against 1.7 in the one-species case).
- $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$ is an advection with the speed $r\,c/(k+c)$ against
  the diffusion $\beta$. **[quick check]** $Pe_h=r h/(2\beta)=10\times0.125/0.04=31$ on $16^3$: a central
  scheme oscillates for $Pe_h>2$. The constant added diffusion $r h/4$ is the smallest that brings the
  effective cell Péclet number to about 2 ($r h/(2(\beta+rh/4))\to2$), but it equals $0.31$ mm$^2/T^*$ on
  $16^3$, sixteen times $\beta$, with a diffusion length $\sqrt{0.31\,T^*}\approx0.56$ mm over the run,
  4.5 elements: the edge of the seed is smeared, hence a quarter of the stress, and the smearing shrinks only
  linearly with $h$ (still rising on $96^3$ in Python). A second-order upwind scheme (MUSCL, minmod)
  adds diffusion only where the solution is not smooth, so the edge stays within one or two cells.
- The time step acts in the same way: the stress is sensitive to the direction $\mathbf n_c$ after the
  onset, when $\sum_i\bar\phi_i$ is small and $\nabla c$ weak, while the growth is not (6 % and 15 %
  against 0.3 %).

**(3) Spoken, EN.** "The growth and the stress converge differently because they are different kinds of
quantity. The growth of the seed is an integral, k alpha times the time integral of phi 3D, and an
integral forgives a smeared edge. The stress is a jump, the jump of alpha and of the stiffness across the
edge of the seed, and a smeared edge lowers that jump directly, so the stress loses one order of accuracy.
The edge is moved by the front term, minus r c over k plus c times n c dot grad phi 3D, which is an
advection: its cell Péclet number is about 31 on the 16 cubed grid, so a central scheme oscillates and
needs stabilisation. A constant added diffusion of r h over four stabilises it, but on 16 cubed that is
sixteen times beta, with a diffusion length of about half a millimetre, and the seed stress comes out at a
quarter of the converged value. With a second-order upwind scheme the diffusion is added only at the edge,
and the stress converges: 0.7 percent without antibiotic, 8 percent with it, and 26 to 29 percent for the
mixed seed, whose edge is nearly discontinuous."

**(3) Spoken, JA.** 「成長量と応力の収束が違うのは、量の種類が違うからです。シードの成長量は積分、$k_\alpha$ と
$\phi_{3D}$ の時間積分で、積分はならされた縁を許します。応力は段差、シードの縁での $\alpha$ と剛性の段差で、ならされた
縁は段差を直接下げるので、応力は精度の次数を一つ失います。縁を動かすのは式 34 の項 $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$
で、これは移流です。セルのペクレ数は $16^3$ で約 31、中心差分は振動するので安定化が要ります。一定の人工拡散 $rh/4$ で
安定にはなりますが、$16^3$ では $\beta$ の 16 倍、拡散の長さが約 0.5 mm で、シードの応力は収束値の 4 分の 1 になります。
2 次の風上法では縁でだけ拡散を足すので、応力が収束します。抗菌薬なしで 0.7 %、ありで 8 %、縁がほぼ不連続な混ぜた
シードで 26〜29 % です。」

**(4) Follow-up.** *"Why not an implicit or a Galerkin scheme in the element?"* The explicit NEM update
is the element as delivered; its stability limit is $\beta\Delta t/h^2\le0.12$ to $0.16$ (Appendix, Bloch
analysis), and I corrected three numerical details of it without changing its equations. A second-order
stabilisation that stays stable on fine meshes inside the element is the first outlook item; the values
of the thesis come from the Python code, which reaches $128^3$.

---

## 5. Where the antibiotic acts: the balance in the point model

**(1) Observation.** Living volume fraction of the seed $0.965\to0.152$ (separate) and $0.961\to0.038$
(mixed); species B is removed first; "where the nutrient is low, A dies as well, because the point model
compares $c^*=c^*_0c$ with $\alpha^*b_i$"; the nutrient in the seed rises from $0.536$ to $0.911$
(separate) and from $0.608$ to $0.973$ (mixed) (`tab:ch5-two-species`). In the pauses of a pulsed dose
the living fraction recovers from $0.17$ to $0.30$ while growth and stress change by less than 1 %.

**(2) Mechanism.**
- Second line of Eq. ch5-point at fixed $\phi_i$, away from the barriers, for a species that dominates
  its point ($A_{ii}\phi_i^2\psi_i\gg\sum_{j\ne i}A_{ij}\phi_j\psi_j$) **[quick check]**:
  $$\eta_i\phi_i^2\dot\psi_i=c^*A_{ii}\phi_i^2\psi_i-\alpha^*b_i\psi_i
  \quad\Longrightarrow\quad
  \frac{\dot\psi_i}{\psi_i}=\frac{c^*A_{ii}}{\eta_i}-\frac{\alpha^*b_i}{\eta_i\,\phi_i^2}.$$
  The nutrient term does not depend on how much of the species is present; the antibiotic term grows as
  $1/\phi_i^2$. $\psi_i$ rises where $c^*A_{ii}\phi_i^2>\alpha^*b_i$ and decays exponentially where it does
  not. With $c^*=c^*_0c$, $\alpha^*=\alpha^*_0a$ and $c^*_0=\alpha^*_0$, species A ($b_1=1$) declines at the
  rate $(a/\phi_1^2-c)\,c^*_0s/\eta_1$ per unit of $T^*$.
- *Where the nutrient is low A dies too*: the sign of $a/\phi_1^2-c$. Next to the supply face $c\approx a$
  and a dense species survives; inside, $c$ is lower and the balance tips. This is the spatial pattern of
  Figure `fig:ch5-3d`, which a model without space cannot show.
- *Why the mixed seed dies faster*: $1/\phi_i^2$. **[paper draft, derived]** In the separate seed A holds
  $\phi_1\approx0.91$ with $c=0.91$, $a=0.98$ at $T^*=1$: rate $(0.98/0.91^2-0.91)\times10=2.7$ per $T^*$;
  over the $2/3\,T^*$ after the onset, $e^{-2.7\cdot2/3}=0.17$, consistent with $\psi_1$ from 0.98 to 0.16.
  In the mixed seed A holds only $\phi_1\approx0.57$: rate $(0.98/0.57^2-0.97)\times10\approx21$ per $T^*$,
  gone within $0.2\,T^*$. Ten times faster, at the same nutrient and antibiotic, because the species shares
  its point.
- *Why B goes first* **[quick check, nuance]**: with $b=(1,2)$ and $\eta=(1,2)$ the antibiotic term
  $\alpha^*b_i/(\eta_i\phi_i^2)$ is the same for both species at equal $\phi_i$, but B's nutrient term
  $c^*A_{22}/\eta_2$ is half of A's, so B's balance tips earlier; and once A takes the share of the mixed
  seed (0.86 without antibiotic), $\phi_2$ is the smaller fraction and $1/\phi_2^2$ the larger. The thesis
  says "the more sensitive one"; in the talk say "B is twice as sensitive and, with its larger
  viscosity, recovers at half the rate, so its balance tips first".
- *Why the nutrient rises*: Eq. 35 consumes with $\kappa=\phi_{3D}\sum_i\bar\phi_i$. With $d=10^8$
  mm$^2/T^*$ the field is quasi-steady, $\nabla^2c=(\Lambda^2/\ell^2)\kappa\,c$ with the Thiele modulus
  $\Lambda^2=g\ell^2/d=1$, depletion length $\ell/(\Lambda\sqrt\kappa)$: 1 mm for living biofilm,
  $1/\sqrt{0.15}=2.6$ mm at $\sum_i\bar\phi_i=0.15$ **[quick check]**. A slab of living biofilm between
  the supply face and the far wall has $c=\cosh((\ell-y)/\lambda)/\cosh(2\ell/\lambda)$, at the centre
  $0.41$ for $\kappa=1$ and $0.82$ for $\kappa=0.15$ **[paper draft; checked]**; the seed means 0.54 and
  0.91 lie above these because the biofilm outside the seed does not fill the cube at $\phi_{3D}=1$. Dead
  biofilm consumes nothing, so the nutrient (and the antibiotic, same equation) rise together.

**(3) Spoken, EN.** "In the point model the living fraction of a species follows a balance: a nutrient
term, c star times A i i, against an antibiotic term, alpha star times b i divided by phi i squared. The
nutrient term does not care how much of the species is there; the antibiotic term grows as one over phi i
squared. So a species that fills its point is killed slowly, and the same species thinly present is killed
at once. This explains three results. First, where the antibiotic acts: next to the supply face the
nutrient is high and species A survives, inside it is low and A dies as well. Second, why the mixed seed
dies ten times faster: each species holds about half the volume, so the antibiotic term is about three
times larger, and the decline rate goes from about 3 to about 20 per unit time. Third, why B goes first:
it is twice as sensitive, and with its larger viscosity it recovers at half the rate. And when the biomass
dies, nothing consumes the nutrient any more, so the nutrient in the seed rises from 0.54 to 0.91; the
depletion length grows from one to 2.6 millimetres."

**(3) Spoken, JA.** 「点モデルでは、各菌種の生きている割合は釣り合いで決まります。栄養の項 $c^*A_{ii}$ と、抗菌薬の項
$\alpha^*b_i/\phi_i^2$ です。栄養の項はその菌がどれだけいるかによらず、抗菌薬の項は $1/\phi_i^2$ で大きくなります。
だから点を満たしている菌はゆっくり殺され、同じ菌が薄くいると一気に殺されます。これで三つの結果が説明できます。
一つ目、抗菌薬が効く場所。供給面の近くでは栄養が高く A が生き残り、内側では栄養が低く A も死にます。二つ目、
混ぜたシードが 10 倍速く死ぬ理由。各菌が体積の半分ほどしか持たないので抗菌薬の項が 2〜3 倍になり、減る速さが
単位時間あたり約 3 から約 20 になります。三つ目、B が先に死ぬ理由。感受性が 2 倍で、粘性が大きいぶん回復が半分の
速さだからです。そしてバイオマスが死ぬと栄養を消費するものがなくなり、シードの栄養は 0.54 から 0.91 に上がります。
枯渇の長さは 1 mm から 2.6 mm に伸びます。」

**(4) Follow-up.** *"Is the balance a threshold? Do you get living and dead regions?"* No. $\psi_i$
decays exponentially rather than reaching zero, and $\phi_i$ evolves at the same time by the first line,
so the outcome is graded: a living fraction between 0 and 1 at every point, as in panel (b) of the 3D
figures. *"Does dead biofilm still spread?"* Yes; the speed $r\,c/(k+c)$ depends on the nutrient only. I
kept the published equation and state it as a limitation.

---

## 6. Why the rule for new points decides the recovery in the maze but not in the cube

**(1) Observation.** Cube, separate example with antibiotic: living volume fraction of the seed $0.147$
with the start values of the nearer seed against $0.152$ with the neighbours' living composition (a few
per cent). Maze (Fig. 11 of Klempt 2024, 2D, $4\,\mu$m cells, antibiotic at the entrance from $T^*=20$):
the living fraction of the colonised region falls to $0.12$ at $T^*=24$ and recovers to $0.35$ at $T^*=60$
with the neighbours' living composition; with the start state it stays at $0.04$. Without antibiotic 0.92.

**(2) Mechanism.**
- The point model has no transport. Eq. 34 moves $\phi_{3D}$ without carrying the point state, so the
  composition travels only through the state a newly reached point is given. The rule acts only on points
  reached *after* the antibiotic onset.
- Cube: with $r=10$ mm$/T^*$ the biofilm fills the cube between the bars by $T^*=0.2$ (§4.7.2) and
  reaches the walls before the antibiotic starts at $T^*=1/3$ **[paper draft: 92 to 95 % of the cells hold
  $\phi_{3D}\ge0.5$ at $T^*=1$]**. Few points are reached after the onset, and the seed itself is never a
  new point, so the seed values differ by a few per cent only.
- Maze: $r=200\,\mu$m$/T^*$ in a $560\,\mu$m maze over $T^*=60$; the biofilm keeps spreading after the
  onset and most of the colonised area at the end was reached after it **[paper draft]**. The fate of a new
  point is then set by the balance of item 5: a point that inherits the neighbours' amounts starts dense
  ($\phi_i$ near the established value, $a/\phi_i^2\approx1.2$) and is killed slowly; a point started from
  the colony's start state ($\phi_i=0.25$ to $0.30$) has $a/\phi_i^2\approx16$ and dies at once
  **[paper draft, quick check: $1/0.25^2=16$, $1/0.9^2=1.2$]**. Hence 0.35 against 0.04.
- The rule is an assumption of this work; I chose inheritance because the biofilm spreads through the
  growth of its living cells, and I report the alternative.

**(3) Spoken, EN.** "The point model has no space, so the composition travels with the biofilm only
through the state that a newly reached point is given. The papers do not fix that state. I give a new point
the living composition of its neighbours; the alternative is the start state of the colony. In the cube
this hardly matters, because with the front speed of the paper the biofilm has reached the walls before
the antibiotic starts, so almost no points are reached after the onset, and the seed values change by a
few percent. In the maze the biofilm keeps spreading after the antibiotic, and most of the colonised area
at the end was reached after the onset. There the rule decides: a point that inherits its neighbours
starts dense, and by the one over phi squared balance it is killed slowly, so the living fraction recovers
to 0.35; a point started from the thin start state dies at once, and it stays at 0.04. So the rule is an
open modelling choice, and I show where it matters."

**(3) Spoken, JA.** 「点モデルには空間がないので、組成は新しく到達した点に与える状態を通してしかバイオフィルムと一緒に
動きません。論文はその状態を決めていません。私は新しい点に隣の生きている組成を与え、別案はコロニーの初期状態です。
立方体ではほとんど効きません。論文の前線の速さでは抗菌薬が始まる前に壁まで届いていて、投与後に到達する点が
ほとんどなく、シードの値は数 % しか変わりません。迷路では投与後もバイオフィルムが広がり続け、最後に占めている面積の
大半は投与後に到達した点です。そこでは規則が決めます。隣から受け継ぐ点は密に始まり、$1/\phi_i^2$ の釣り合いで
ゆっくり殺されるので、生きている割合は 0.35 まで回復します。薄い初期状態から始める点は一気に死に、0.04 のままです。
規則はモデル化の未決の選択で、どこで効くかを示しました。」

**(4) Follow-up.** *"Which rule is right?"* Neither paper says; a species-resolved transport would
replace the rule. Inheritance follows from the biofilm spreading by the growth of living cells; the start
state is what the ANSYS runs used by position. I report both and where they differ.

---

## 7. Why rigid bars raise the seed stress about tenfold

**(1) Observation.** Four bars after Fig. 9 of Klempt 2024 (0.25 mm square, 0.5 to 0.75 mm from the
centre, no biofilm, no nutrient, frictionless sliding), $32^3$, no antibiotic: seed von Mises stress
$7.1\times10^{-3}$ against $7.4\times10^{-4}$ Pa without bars in the separate example (about 10 times),
$9.5\times10^{-4}$ against $2.1\times10^{-4}$ Pa in the mixed one; the growth $\alpha-1$ of the seed is 9
and 21 % lower; the biofilm fills the cube between the bars by $T^*=0.2$.

**(2) Mechanism.**
- Growth is $\mathbf F_g=\alpha\mathbf I$, stress comes only from $\mathbf F_e=\mathbf F\mathbf F_g^{-1}$.
  A body that can grow freely takes $\mathbf F=\alpha\mathbf I$, $\mathbf F_e=\mathbf I$, zero stress
  (verified: $1.9\times10^{-11}$ MPa, seven decades below the constrained case, §4.4.1). Fully
  constrained isotropic growth gives $J_e=\alpha^{-3}$ and the pure pressure
  $\boldsymbol\sigma=\frac{2}{D_1}(J_e-1)\mathbf I$ (Eq. `ch5-constrained`), with zero deviatoric and
  zero von Mises stress.
- Without bars the seed is held only by its soft surroundings, $E\phi_{3D}^2$, so it swells almost freely
  and the von Mises stress comes from the jump at its edge, $7.4\times10^{-4}$ Pa. The bars block the
  swelling in some directions but not in others (the biofilm slides along them), so $\mathbf F_e$ becomes
  anisotropic and a deviatoric stress appears. **[quick check, order of magnitude]** Growth $\alpha-1$
  constrained in one direction with the others free gives $\sigma=-E(\alpha-1)$ and a von Mises stress of
  $E(\alpha-1)=10\ \text{Pa}\times8.3\times10^{-4}\approx8\times10^{-3}$ Pa, the order of the
  $7.1\times10^{-3}$ Pa observed; the fully constrained pressure would be
  $K\cdot3(\alpha-1)=166.7\times2.5\times10^{-3}\approx0.4$ Pa, but it carries no von Mises stress.
  So the tenfold rise is the step from "constrained by soft surroundings" to "constrained by rigid walls in
  some directions".
- The mixed seed rises less (4.5 times) because its $\alpha$ and stiffness are lower (item 3).
- The 9 and 21 % lower growth: the thesis gives no mechanism. **[interpretation]** The bars hold no
  nutrient, so $c$ and $\mathbf n_c$ around the seed change and the term
  $-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$ carries more $\phi_{3D}$ out of the seed; by Eq. 36 the
  seed mean of $\alpha-1$ follows.
- Single-grid values: on $32^3$ the stress without bars is 14 % above its extrapolated value, small against
  the factor 10.

**(3) Spoken, EN.** "Growth itself causes no stress. A body that can grow freely simply becomes larger,
with F equal to alpha times the identity and zero stress; I verified that in ANSYS, the free case gives
stress at round-off. Stress appears where growth is constrained. In the cube the seed is held only by the
soft surroundings, whose stiffness scales with phi 3D squared, so it swells almost freely and the stress
comes from the jump at its edge, seven times ten to the minus four pascal. The rigid bars block the
swelling in some directions, and then the elastic deformation is no longer isotropic and a deviatoric
stress appears. Its order is the modulus times the growth, 10 pascal times 8 times ten to the minus four,
about 8 times ten to the minus three pascal, which is what I see: 7 times ten to the minus three, about
tenfold. For a fully constrained seed the closed form gives a pure pressure with zero von Mises stress, so
the bars are the case in between: constrained in some directions, free in others."

**(3) Spoken, JA.** 「成長そのものは応力を生みません。自由に成長できる物体は $\mathbf F=\alpha\mathbf I$ で大きくなるだけで
応力はゼロです。ANSYS で確かめ、自由な場合の応力は丸め誤差でした。応力は成長が拘束される所に出ます。立方体では
シードは剛性が $\phi_{3D}^2$ に比例する柔らかい周囲にしか支えられず、ほぼ自由に膨らむので、応力は縁の段差から来る
$7\times10^{-4}$ Pa です。剛な棒は膨らみをいくつかの方向で止めるので、弾性変形が等方でなくなり偏差応力が出ます。その
大きさは弾性率と成長量の積、10 Pa と $8\times10^{-4}$ で約 $8\times10^{-3}$ Pa、見ている $7\times10^{-3}$ Pa、約 10 倍と
合います。完全に拘束したシードは閉じた形では純粋な圧力で von Mises 応力はゼロなので、棒はその中間、ある方向は拘束、
ある方向は自由、という場合です。」

**(4) Follow-up.** *"Why is the growth lower with the bars?"* The thesis reports the numbers (9 and
21 %) on one grid and gives no mechanism; the bars hold no nutrient, so the nutrient field and the
direction $\mathbf n_c$ around the seed change, and with them how much $\phi_{3D}$ the term
$-r\frac{c}{k+c}\mathbf n_c\cdot\nabla\phi_{3D}$ carries out of the seed. I would say that as a reading,
not as a result.

---

## 8. Chapter 3: identifiability, DH at stage (iii), the pH check

**(1) Observation.** 5, 6, 9 and 8 of the 15 entries of $\mathbf A$ identified in CS, CH, DS, DH; in the
commensal conditions only entries among *So*, *An*, *Vei*. DH reported from stage (iii) because at stage
(iv) criterion (e) fails (medians of $a_{22}$ differ by 0.79 posterior standard deviations across seeds,
limit 0.5). pH: regression fitted once on measured pH and composition ($R^2=0.88$), applied to the
predicted composition; DH follows the measured rise from 6.4 to 7.0 ($R^2=0.71$, RMSE 0.097, 85 % inside
the 90 % band); CH matches the level (RMSE 0.082) but not the variation ($R^2<0$).

**(2) Mechanism.**
- *Identified* means: pooled posterior standard deviation below 0.8 times that of a uniform distribution
  on the final box; the data narrow the entry. An entry can be narrowed only if the data respond to it.
  In CS and CH *Fn* and *Pg* stay below 1 % from day 3 on, so the terms $c^*A_{ij}\phi_j\psi_j$ with
  $j\in\{Fn,Pg\}$ are small in every equation and the 9 entries involving *Fn* or *Pg* cannot be narrowed
  **[quick check: 9 of 15 involve species 4 or 5, leaving 6 among So, An, Vei; CH identifies all 6, CS 5
  plus one sign-only]**. In DS and DH *Fn* and *Pg* reach 10 to 27 %, and 9 and 8 entries are identified,
  now including theirs. With 25 to 30 data points per condition and 15 unknowns, part of the matrix stays
  at the prior; stating which part is the result.
- *DH at stage (iii)*: the final stage uses the full likelihood on the narrowed box; two sets of
  coefficients explain the DH data about equally well (the size of the late *Pg* increase: seed MAPs give
  ratios 9.6 to 18 against the measured 14.9, 90 % interval 0.6 to 32), the three seeds settled in
  different modes, and criterion (e) is exactly the test for that. Stage (iii) passes all criteria, so I
  report it and state the split as a limitation. The fit is not worse (RMSE 0.073 of the viable
  composition).
- *pH check*: the pH did not enter the likelihood. The regression is a map from composition to pH fitted
  on the data; applying it to the posterior composition tests whether the predicted composition carries
  the information that the measured pH carries. In DH it does (the rise follows the composition change);
  in CH the composition hardly changes, so a regression on it cannot reproduce the variation of the pH,
  only its level. It is a check of the predicted composition, not a validation of a pH model.

**(3) Spoken, EN.** "The result of the calibration is a statement, entry by entry, of what the data
determine. An entry is identified when its posterior is narrower than 0.8 of the uniform prior on the
box. In the commensal conditions only five and six of the fifteen entries are identified, all among the
early colonisers, because Fusobacterium and Porphyromonas stay below one percent there and the data carry
no information on their interactions. In the dysbiotic conditions, where they reach 10 to 27 percent,
nine and eight are identified. For DH I report the third stage, because at the fourth the three seeds
settled in different posterior modes: the data do not fix the size of the late Pg increase, and my
criterion (e) is exactly the test that catches it. The pH was not used in the inference. A regression from
measured composition to measured pH, applied to the predicted composition, follows the measured rise in DH
with an R squared of 0.71, so the predicted composition carries that information. It is a check, not a
validation."

**(3) Spoken, JA.** 「較正の結果は、成分ごとにデータが何を決めるかの記述です。事後分布の幅が箱の一様分布の 0.8 倍より
狭ければ同定とします。常在の条件では 15 個中 5 個と 6 個しか同定できず、全部が初期定着菌の間です。Fn と Pg が 1 % 未満で、
その相互作用の情報がデータにないからです。病的な条件では 10〜27 % に達するので 9 個と 8 個です。DH は段階 (iii) を
報告します。段階 (iv) では 3 つの seed が別のモードに落ち着いたからで、後半の Pg の増え方の大きさをデータが決められず、
基準 (e) はまさにそれを捉える判定です。pH は推定に使っていません。測定した組成から pH への回帰を予測した組成に当てると、
DH では測定した上昇を $R^2=0.71$ で追うので、予測した組成がその情報を持っています。確認であって検証ではありません。」

**(4) Follow-up.** *"Why 0.8 and not another threshold?"* It is a stated convention: below 0.8 the
posterior is clearly narrower than the prior; the recovery test on synthetic data reproduced the
classification for 14 to 15 of 15 entries, which supports the threshold, while the intervals themselves
are read as a lower bound on the uncertainty (true value inside the 90 % interval for only 2 of the 6
identified CH entries).

---

## Numbers to keep at hand

| what | value | source |
|---|---|---|
| bound of the seed growth | $k_\alpha T^*=10^{-3}$ at $T^*=1$; seed at $9.16\times10^{-4}$ (92 %) | Eq. 36, `tab:ch5-gci` |
| growth with/without antibiotic | $9.22/9.16$ (separate), $8.54/8.56\times10^{-4}$ (mixed) | `tab:ch5-gci` |
| living fraction of the seed | $0.965\to0.152$, $0.961\to0.038$ | `tab:ch5-two-species` |
| nutrient in the seed | $0.536\to0.911$, $0.608\to0.973$ | `tab:ch5-two-species` |
| seed stress, extrapolated | $6.55/6.70\times10^{-4}$ Pa (separate), $4.50/3.99$ (mixed) | `tab:ch5-gci` |
| GCI growth / stress | 0.7 to 3.3 % / 0.7, 8.3, 26 to 29 % | `tab:ch5-gci` |
| one grid, antibiotic effect | $6.73\to6.00$ ($64^3$, 0.0025); $6.56\to5.68$ (0.025) | tables |
| cell Péclet number | $rh/(2\beta)\approx31$ on $16^3$ | §4.7.1 |
| constant added diffusion | $rh/4$: a quarter of the stress on $16^3$ | §4.7.1 |
| largest $\phi_{3D}$ in the mixed seed | 0.60, 0.69, 0.73 on $32^3$, $64^3$, $128^3$ | §4.7.1 |
| front speed factor | $c/(1+c)$: $0.35\to0.48$ (37 %) | quick check |
| depletion length | $\ell/\sqrt{\kappa}$: 1 mm $\to$ 2.6 mm | paper draft, checked |
| decline rate of A | 2.7 (separate) and 21 (mixed) per $T^*$ | paper draft, checked |
| maze | living 0.92; with antibiotic 0.12 at $T^*=24$, 0.35 at 60 (inherit) / 0.04 (start state) | §4.7.2 |
| cube, rule for new points | 0.147 against 0.152 | §4.7.1 |
| bars | $7.1\times10^{-3}$ vs $7.4\times10^{-4}$ Pa; mixed $9.5$ vs $2.1\times10^{-4}$ Pa | §4.7.2 |
| Ch. 3 | 5/6/9/8 of 15; DH stage (iii), $a_{22}$ 0.79 sd; pH DH $R^2=0.71$ | Ch. 3 |

## Points to check in the thesis text (found while writing this; not changed here)

1. Caption of Figure `fig:ch5-bars` (ch4_ansys.tex, "mixed example, which starts at $\phi_{3D}=0.55$")
   and the same wording in the paper draft (caption of the geometry figure and §4.5 "when they start mixed
   at $\phi_{3D}=0.55$"): in the Python code the field seed starts at $\phi_{3D}=1$ in both examples
   (`phi = (seed | edges).astype(float)`), and the result files record a domain maximum of $\phi_{3D}$ of
   1.000 at $T^*=0$ for the mixed runs. The 0.55 is the point model's $\phi_1+\phi_2$ (the recorded
   living fraction of the seed starts at 0.549). The maximum falls below 1 only from $T^*\approx0.4$ on.
   The isosurface level 0.2 of that figure may still be the right choice for the mixed example at
   $T^*=0.2$ (on $32^3$ with bars), but "starts at 0.55" is not what the field does.
2. "Species B, the more sensitive one ($b_2=2$), is removed first": with $\eta_2=2$ the antibiotic term
   $\alpha^*b_i/(\eta_i\phi_i^2)$ is the same for both species; B goes first because its nutrient term
   $c^*A_{22}/\eta_2$ is half of A's and because its volume fraction in the mixed seed is the smaller one.
   The observation is right; "because it is more sensitive" is half the reason. Wording only.
3. "The largest $\phi_{3D}$ in the seed still rises with the grid (0.60, 0.69, 0.73)": the recorded
   quantity in the result files is the maximum over the whole domain at the end of the run. Fine if the
   maximum lies in the seed, which the figures suggest; the sentence could say "in the domain".
