# 引き継ぎ：IKMHIWI03 → 別の PC の Claude（2026-10-20 21 時）

この文書だけ読めば続けられるように書いた。詳しい記録は `ansys_usermat/apdl/IKMHIWI03_TO_CLOUD.ja.md` の **A-19**。
ブランチは `claude/plan-next-hxjjve`（IKMHIWI03 の結果はここに来る。master に入れるのはクラウドだけ、CLAUDE.md のとおり）。

## 今日（20 日）の IKMHIWI03 の作業

1. **10/19 の依頼 3 本は完了・push 済み**
   - 8³ CH を最終 MAP（Tmcmc202601 `baa0a37`、10/16）で：`results/2026-10-cal5/w8_CH_g1_s015_final_d1*.json`。
     成長と応力は旧 MAP（6e18e0b）と全桁一致。T\* = 1 の割合 So 0.553、An 0.132、Vd 0.147、Fn 0.088、Pg 0.080、
     Σφψ 0.811 → 0.231。Day 1 の組成は旧 MAP と同じ。
   - 円柱の図 4.4 のデータ：`results/cylinder_alpha_sweep.json`（max u_r は α − 1 = 0.2 で 3.56e−2）。
   - 16³ CH、40 時刻の場（アニメーション用）：`results/2026-10-cal5/fields_w16_CH_g1_s015_final_d1_dense.json`。
     旧 MAP の 16³ と α − 1・vM・φ・c が全桁一致。
2. **オリバーと決めた結合の ANSYS 実装**（`CLOUD_TO_IKMHIWI03_2026-10-20.ja.md`）：A-19 を参照。
   - 点モデルはブラックボックス（usermat モード 11）、c\* と α\* は場から時間で線形補間、
     消費 g φ L c（L = Σφ_i ψ_i）、前線項 |∇φ|、2 つめの栄養の場を抗菌薬に。
   - 8³ の確認 6 本は完了・push 済み（`results/2026-10-living/`）。
   - 16³ の本番 8 本は IKMHIWI03 で連鎖中（切り離し、終わると自動で push）。
     ログ `F:\biofilm_upf_living\_living_prod16.log`。

## 21 日夜までの追記（A-20〜A-22、詳しくは `IKMHIWI03_TO_CLOUD.ja.md`）

- **ビルドは `F:\biofilm_upf_living4` が最終**：living の A-19 の変更に加えて
  - usermat（repo の `callsite/phi_mode_exec.inc`）：最初の荷重ステップは場だけを解く（prop(60) = 1、`make_living_deck.py --start-step 0.0025`）。
  - USSFin（オリバーの要素、repo には入れない）：**L を毎回読む 1 行**（今までは初期化の中で 1 回だけ読んでいて L が固定だった、A-22）。
  - 点モデルのファイル：`newton 50 1e-20`（最初のステップの Newton を収束させる、A-20）。
  - 差分：`F:\oliver_element_patches_2026-10-21\`（**ユーザーが OneDrive などに写す**）。
- **使える結果は `v2f_*` と `fsens_*` だけ**（`results/2026-10-living/`）。lv8_*、v2_8_*、v2s_*、sens_*、i8* は L が固定の版で使わない。
  - 8³：v2f_8_sep_ab、v2f_8_mix_ab、v2f_8_sep、v2f_8_mix。感度 fsens_phimin001、fsens_dt00125、fsens_gab4、fsens_a10、fsens_r10（r = 10 は
    要素が潰れて止まる）、閾値 fsens_a40、fsens_a60。
  - 16³：v2f_16_mix_ab → v2f_16_sep_ab → v2f_16_mix → v2f_16_sep（1 本 約 70 分）。40 時刻の場 fields_v2f_16_{mix,sep}_ab_dense.json。
  - 24³：v2f_24_sep_ab（メッシュ収束、終わるかどうかは時間次第）。
- 照合：`check_living_trace.py <job> --eco abaqus_composition/eco_fritsch_2sp_a50.txt`（点モデルだけを見る。L の受け渡しは見えない）。

## 別の PC で次にやること（優先順）

1. **check_living_trace.py を回す**（クラウド、Python と jax が要る）：
   `python ansys_usermat/apdl/check_living_trace.py lv8_sep`（ほかに lv8_sep_l1、lv8_sep_ab、lv8_mix）。
   PASS なら、点モデルを ANSYS の中で正しく呼べている。
2. **L = 1 の 1 本が L ≈ 0.97 の 1 本と大きく違う件**を Python 版（`coupled_living_consumption.py`、同じ設定）と比べる。
   シードの φ 1.0 → 0.5、c 0.76 → 0.36（A-19 の 2.）。前線項の cos の向きの切り替えが効いているように見える。
3. **r（10 mm/T\*）は ANSYS では解けない**：8³ で 0.3 分で要素が潰れた。r/10 だけを修論に使うか、
   USSFin の前線項を風上差分にするか決める（決めるのはユーザー）。
4. 16³ の本番の JSON が push されたら、Python 版と表で比べる（シードと四隅の L・A の割合・c・a、α − 1、vM）。
5. **確認 2（L = 1 で今の消費と全桁一致）は成り立たない**：オリバーの要素の今の消費は g φ（0 次）で、g φ c ではなかった。
   修論の書き方（「要素の消費を g φ から g φ L c に変えた」）を Meisam・オリバーに伝えるか検討。

## repo に入っていないもの（公開 repo なので入れない）

- **オリバーの要素の改変版**（USSFin・USolBeg・usercm.inc・Usermat の変更）は `F:\biofilm_upf_living\` にだけある。
  変えた場所の説明は A-19 にある。ANSYS で作り直すには、オリバーの元のソースに同じ変更を手で入れる必要がある。
- **退避フォルダ** `ansys_usermat/apdl/f_backup_2026-10-20/`（9 MB、デッキの zip を含む）：IKMHIWI03 のローカルだけ。
  自動モードの判定で git に追加できなかった。
- オリバーの要素とフェリックスのファイル（`F:\felix_private`、4 GB）：**ユーザーが OneDrive などに手でコピーする**
  （Claude からのコピーは判定で止められた）。
- IKMHIWI03 は 21 日まで。ANSYS の作業はこれで終わり。別の PC には ANSYS が無い前提。

## 作業の決まり（CLAUDE.md の要点）

- 返答は日本語で短く。文書の言語はそのまま。
- commit に AI の痕跡を残さない（Co-Authored-By などを付けない）。ファイルは名前で指定して stage する（`git add -A` 禁止）。
- フェリックスのコードと理論は repo に入れない。オリバーの要素のコードも同じ。
- 修論・スライドの表記：φ、α（F_g = αI、成長は α − 1）、k_α、点モデルは φ_i、ψ_i、φ̄_i = φ_i ψ_i。一人称は I。
