# 修論の成績 1.0 に効く要因 — 規程ベース

> 作成 2026-08-13。推測ではなく `_thesis_forms/` の一次資料から。
> 出典: MPO §17 / §20、`04_Standards_Abschlussarbeiten.pdf`、`05,06_Verfahrensrichtlinien`

---

## 0. 前提 — 修論は「授業6個分」

DD は **90 LP で修了**。修論は **30 LP**。

**MPO §20(3) Satz 2:**

> soweit nicht in den Anlagen … besondere Gewichtungen ausgewiesen sind,
> **gewichtend die darin aufgeführten Leistungspunkte verwendet**

Anlage 1.4 に特別な重みの指定はない → **LP がそのまま重み**になる。

```
授業1つ  =  5 LP
修論     = 30 LP  ＝  授業 6個分
```

成績のつく授業が12個（60 LP）なら:

```
総合成績 = 授業12個 ＋ 修論（6個分）= 18個分の加重平均
         → 修論だけで 1/3
```

### ⚠️ 分母に入るのは「成績のついた科目」だけ

**§20(3) Satz 1:** `das arithmetische Mittel der Noten aller **benoteten**
Prüfungsleistungen`

**Studienleistung（合否のみ）は分母に入らない。**
履修リストに `SL` 表記の科目があるため、無評価の分だけ**修論の比率は上がる**。

| 成績付き科目 | 修論の比率 |
|---|---|
| 60 LP | 33% |
| 45 LP | **40%** |

→ **どの科目が benotet か QIS で確認する**（`python3 qis_fetch.py --grades`）。

### 帰結：下がったときに取り返せない

修論が **1.3 → 1.7**（0.4 悪化）＝ **6個分が一斉に 0.4 悪化**。
取り返すには残り12個を**全部 0.2 ずつ良くする**必要がある。
**試験は全部終わっているので不可能。**

> **修論は「稼ぐ場所」ではなく「絶対に崩してはいけない場所」。**

---

## 1. 🔴 最重要の構造：1.0 は「2人とも 1.0」でしか出ない

**MPO §17(2):**

> Wird eine Prüfungsleistung durch zwei Prüfende bewertet, … **Die Note errechnet
> sich in diesem Fall aus dem Durchschnitt der von den Prüfenden festgesetzten
> Einzelnoten.**

**MPO §20(3) Satz 4（§17(2) に準用）:**

> **nur die erste Dezimalstelle** hinter dem Komma … alle weiteren Stellen werden
> **ohne Rundung gestrichen**.

| Erstprüfer | Zweitprüfer | 平均 | 記載される成績 |
|---|---|---|---|
| 1.0 | 1.0 | 1.00 | **1.0** |
| 1.0 | 1.3 | 1.15 | **1.1** |
| 1.3 | 1.3 | 1.30 | **1.3** |
| 1.0 | 1.7 | 1.35 | **1.3** |
| 1.3 | 1.7 | 1.50 | 1.5 |

> ### ⚠️ 帰結：**Zweitprüfer が成績の半分を持つ。**
> しかも Zweitprüfer は**別研究所**（様式の要求）＝ **経過を知らずに完成品だけを読む。**
> Erstprüfer は文脈を補完して読んでくれるが、**Zweitprüfer は補完してくれない。**
>
> → **「事情を知らない読者が、本文だけで納得できるか」が 1.0 の分岐点。**

---

## 2. 狙うべき閾値は 1.0 ではなく 1.3

**MPO §20(4):**

> Lautet die Gesamtnote **mindestens 1,3 oder besser** und ist die Masterarbeit mit
> der Note **mindestens 1,3** bewertet, so wird das Prädikat **„mit Auszeichnung"**
> verliehen.

**証書に載る称号の閾値は 1.3。** 1.0 と 1.3 の差は数字だけで、称号は変わらない。
逆に **1.5 になると称号を失う**（＝ 1.3 と 1.7 の組み合わせで即アウト）。

> **実効目標：両審査者を 1.3 以内に収める。**
> 「1.0 を狙う」より「1.7 を1つも出さない」ほうが設計として正しい。

### 称号は「鍵が2つ」— 修論は独立した拒否権を持つ

§20(4) の条件は **AND**。片方でも欠けると称号は付かない。

| | 条件 |
|---|---|
| **鍵1** | **総合成績 ≤ 1.3** |
| **鍵2** | **修論の成績 ≤ 1.3**（単独で） |

**§20(3) Satz 4** は小数第2位以下を**切り捨て**（四捨五入しない）
→ 鍵1 の実条件は加重平均 **1.40 未満**（1.39 でも 1.3 と表示され通る）。

**鍵2 には平均が効かない。修論が 1.5 以上になった時点で、
授業が全部 1.0 でも称号は付かない。**

| 修論 | 鍵2 | 鍵1（授業12個の平均） | 称号 |
|---|---|---|---|
| **1.0** | ✅ | **1.60 未満なら可**（余裕あり） | 条件次第 |
| **1.3** | ✅ | **1.45 未満なら可** | 条件次第 |
| 1.5 | ❌ | — | **不可** |
| 1.7 | ❌ | — | **不可** |

> ⚠️ **訂正（初版の誤り）:** 初版は「修論 1.7 なら授業平均 1.25 未満で可」と
> 書いていたが**誤り**。それは鍵1 だけの計算で、**鍵2 が独立に落ちる**。
> **修論が 1.3 を超えた瞬間、授業側で挽回する道は存在しない。**

---

## 3. §17(1) が定める 1.0 の水準

> 1,0; 1,3 = „sehr gut" = **eine besonders hervorragende Leistung**

LUH は „hervorragend" ではなく **„besonders hervorragend"** と書いている。
2.0（„erheblich über den durchschnittlichen Anforderungen"）との距離は大きい。

---

## 4. 学部が明文化している「中身」の基準（04_Standards）

### ✅ 必須条件

> Ausgangspunkt ist eine **wissenschaftliche Fragestellung** oder ein
> **ungelöstes ingenieurwissenschaftliches Problem**.
> **Das Ausführen einer „Abarbeitungsliste" reicht nicht aus.**

**「やったことリスト」は不可**と名指しで書かれている。
パイプラインを何本組んだかではなく、**どの問いに答えたか**で構成する。

### 🔴 修士に固有の要求 — ここが 1.0 の分水嶺

> **Besonders bei Masterarbeiten ist eine kritische Reflexion der Ergebnisse
> erforderlich.**

この文書の中で **„Besonders bei Masterarbeiten"** と限定されているのはこの1点だけ。
＝ **学部との差はここで測られる。**

計算系の修論が落とすのは典型的にここ:
- 何が**再現できなかった**か
- どの仮定が**効いていた**か（感度・反証）
- 手法が**通用しない範囲**はどこか
- 数字が良かったのは**手法のおかげか、データのおかげか**

> RL のレポートで RND を落とした判断を書いたのと同じ筋。**あれを本文の各章でやる。**

### その他

> Dies gilt hinsichtlich der Vorgehensweise und eingesetzten Methoden, der zu
> berücksichtigenden Theorien sowie der **Darstellungsform** in der schriftlichen
> Arbeit und der **formalen Kriterien (Zitation, Sprache)**.

体裁・引用・言語も明示的な評価対象。

---

## 5. 🎯 いま取れる最大のレバー：Bewertungsbogen を要求する

**04_Standards:**

> In dem Gespräch sollten die Studierenden auch über die **vom Institut
> festgelegten Bewertungskriterien** informiert werden. **Ein vom Institut
> bereitgestellter Bewertungsbogen wird empfohlen. Er sollte den Studierenden
> zugänglich gemacht werden.**

**評価シートを学生に見せることが学部の指針として明記されている。**
＝ **要求してよい。** これがあれば採点基準が推測からチェックリストに変わる。

**Ausgabe（9月中旬）のタイミングで聞くのが自然。**

---

## 6. 規程が用意している「無料の保険」

04_Standards が推奨している手続き。使わないと損。

| いつ | 何 | 効果 |
|---|---|---|
| 着手時 | **Exposé**（問い・手法・作業計画・日程・目標成果） | 問いの形を先に審査者と合意できる |
| 着手時 | 作業計画・指導計画を**書面化して双方署名** | 期待値のズレを潰す |
| **最初の1/3が過ぎた頃** | **Feedbackgespräch**＋Sachbestandsbericht | 構成・文献・章ドラフトへの是正が早期に入る |
| 提出前 | **Rohfassung を Betreuer と議論** | 提出前に指摘を回収できる |

> 3ヶ月で走る場合、「最初の1/3」は **Ausgabe から約4週間 ＝ 10月中旬**。
> ここで一度フィードバックを取ることを、最初から計画に入れる。

---

## 7. まとめ

0. **修論 = 授業6個分。総合成績の 1/3（無評価科目が多ければそれ以上）。**
   **下がっても後から埋め戻せない。**
1. **Zweitprüfer が半分。文脈を知らない読者に本文だけで通じること。**
2. **目標は「両者 1.3 以内」。1.7 を1つも出さない設計。**
   称号は鍵2つの AND。**修論が 1.3 を超えたら授業で挽回できない**（拒否権）。
   総合側の実条件は **1.40 未満**（切り捨てのため 1.39 でも可）。
3. **「やったことリスト」ではなく「問い → 答え」の構成。**
4. **限界・反証・効いた仮定の批判的考察を各章に置く**（修士固有の要求）。
5. **Bewertungsbogen を要求する**（学部指針が学生に見せろと書いている）。
6. **10月中旬の Feedbackgespräch と提出前の Rohfassung レビューを予定に入れる。**

---

## 関連

- `_thesis_forms/README.md` — Zulassung・提出・採点の手続き
- `MPO_summary.md` — MPO 全体のまとめ
- `THESIS_TIMELINE_RULES.md`
