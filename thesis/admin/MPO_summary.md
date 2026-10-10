# LUH Maschinenbau M.Sc. — Master-Prüfungsordnung (MPO) 修論関連まとめ

出典：[LUH Fakultät für Maschinenbau MPO](https://www.uni-hannover.de/fileadmin/luh/content/pruefungsamt/formulare/maschb/po/f_maschbau_mpo.pdf)
（原本 `notes/LUH_Maschinenbau_MPO.pdf`／復号版 `notes/LUH_Maschinenbau_MPO_decrypted.pdf`、20頁、底本 Verkündungsblatt 11/2015＝2015-08-25）
作成 2026-06-04。§番号・独語原語つき。

---

## A. 修論モジュールの基本 §7(1), §2, Anlage 1.4
- **構成**：修論本体（Masterarbeit）＋**場合により口頭試問**（mündliche Prüfungsleistung, Anlage 1.4 による）。
- **30 LP**（Leistungspunkte）。全課程 120 LP・4学期・標準修業 **2年**（§2）。
- **目的**：所定期間内に専門分野の問題を **独力で科学的方法** により扱えることを示す。

## B. 着手・テーマ §7(2)(3), §12(3)
- **Zulassung（受験許可）が必要**。許可要件＝修士課程で **最低 60 LP 取得** ＋ Anlage 1.4 の追加要件。例外は Prüfungsausschuss が判断。
- **テーマ交付は Zulassung 後のみ**。
- テーマは **着手期間の最初の 1/3 以内に 1 回だけ返却可**（Rückgabe）。返却後 6 か月以内に再登録しないと、第一審査者が指定テーマを送付。

## C. 期限 §7(4), §15(2) ★重要
- **テーマ交付（Ausgabe）から 6 か月以内** に提出（**書面＋電子** 両方）。
- 審査：両審査者が **6 週（遅くとも 10 週）以内** に採点。
- **延長**：重大事由（wichtige Gründe）を **遅滞なく書面で** 申請 → Prüfungsausschuss が判断。
  - 📗 **手続き指針（06_Verfahrensrichtlinien_EN）で場合分けが確定（2026-08-18 追記）**：
    **病気以外**（近親者の死亡・家族の介護・競技スポーツ・犯罪被害）＝ **一度きり・最大4週間**／
    **病気**（診断書必須）＝ **最大 1/3**（修論なら2ヶ月）。
  - **納品遅延・指導教員の休暇・機器故障は正当事由に当たらない**（研究所側の責任）。
  - 🔴 **延長申請は学生と指導教員の連名でしか出せない。**
- **提出徒過・着手後の撤退 = 不合格（nicht bestanden）**。ただし重大事由が認められれば「未実施」扱い。

## D. 提出時の宣誓 §7(5) — Selbstständigkeitserklärung
独力執筆／記載以外の出典・手段不使用／引用箇所明示／同一・類似形で他試験機関に未提出。
→ **修論末尾の Selbstständigkeitserklärung が既に対応済**（`main.tex`）。

## E. 言語 §7(7)
独語が原則、**審査者（Junker）の合意で英語可**、例外で他言語。**独語 Zusammenfassung の義務規定なし** → **英語のみで規定準拠**（現状の英語 Abstract のみで可、独語 Zusammenfassung は `\iffalse` で無効化済）。

## F. 採点 §17, §7(6)
- 修論は **2 名の審査者**（zwei Prüfende）。**両者が ausreichend(4.0) 以上で合格、点数は両者の平均**。
  > ⚠️ **訂正（2026-08-13）：「主査 Junker・副査 Soleimani」は成立しない。**
  > Zulassung 様式が **Erst- と Zweitprüfer は異なる研究所**から来ることを要求しており、
  > 両先生とも IKM（`@ikm.uni-hannover.de`）。**IKM 外の副査が必要。**
  > → `_thesis_forms/README.md`
- 成績段階：**1.0/1.3 sehr gut｜1.7–2.3 gut｜2.7–3.3 befriedigend｜3.7–4.0 ausreichend｜5.0 nicht ausreichend**。
- モジュールに口頭試問が含まれる場合は Anlage 1.4 の重みで合成（口頭の形式・時間は Absprache で決定）。

## G. 総合成績 §20
- Gesamtnote ＝ 採点済 Prüfungsleistung の **LP 加重平均**（**小数第1位まで、以降切り捨て・四捨五入なし**）。
- ≤1.5 sehr gut／≤2.5 gut／≤3.5 befriedigend／≤4.0 ausreichend／>4.0 不合格。
- **mit Auszeichnung（優等）**：総合 ≤1.3 **かつ** 修論 ≤1.3。

## H. 再試・最終不合格 §14(2)
- 不合格の修論は **1 回のみ再試可**。再度不合格＝ **endgültig nicht bestanden**（全体不合格）。

## I. 証書 §21
- 証書類は **独語で発行＋Prüfungsamt が英訳** を発行。修論タイトルと成績を記載。

---

## J. Anlage 1.4 — 修論モジュールの中身（★口頭試問あり）
| 項目 | 内容 |
|---|---|
| Semester | 4 |
| Zulassung 要件 | **60 LP ＋ Studienarbeit ＋ 20週 Praktikum** |
| **Prüfungsleistung**（採点対象） | **Masterarbeit のみ** |
| **Studienleistung**（合否のみ） | **Präsentation** |
| LP | 30 |

> ⚠️ **訂正（2026-08-13）：発表は成績に入らない。**
> 以前ここに「Prüfungsleistung = Masterarbeit ＋ Präsentation（採点される発表）」と
> 書いていたが誤り。Anlage 1.4 の表で **Präsentation は Studienleistung の列**にある
> ＝ **bestanden / nicht bestanden の合否判定のみ**。
> §7(6) は「モジュールが**複数の Prüfungsleistung** を含む場合に Anlage 1.4 の重みで合成」
> と定めるが、Anlage 1.4 の Prüfungsleistung は **Masterarbeit 1つだけ**。
> → **モジュール成績＝修論本体の成績**（両審査者の平均、§17）。

- 発表自体は**必須**（落とすとモジュールが完成しない）。Anlage 2.1 定義：
  *„Das Kolloquium umfasst eine Präsentation und eine anschließende Diskussion über
  Fragestellung, Methoden und Ergebnisse der Arbeit.“*（発表＋問い・手法・結果の議論）。
- **形式・時間は審査者と Absprache（合意）で決定**。
- 発表は **Bearbeitungszeit の期間内**に実施（Verfahrensrichtlinien）。**採点の後ではない。**
- → **defense/Kolloquium の準備は必須ワークストリーム**（overview デッキが土台）。

## 🔴 このまとめの底本は古い（2026-08-18 追記）

**適用されるのは 2024-10-01 施行版**（§24(2)：在籍者全員に適用）。
このまとめは **2015 年版**を底本にしているので、以下の点で実際と食い違う：

| 箇所 | 2015 版（このまとめ） | **2024 版（正）** |
|---|---|---|
| 修論モジュールの附則 | Anlage **1.4** | **Anlage 1.7** |
| 実務経験 §4(3) | 合計 **20 週** | **mindestens 12 Wochen** |
| §14(3) | 毎学期15 LP **＋累積 15×学期数** | **毎学期15 LP のみ**（累積は削除） |
| 副査の所属 | §5 のみ（学部の縛りが読めた） | **§7(8) 新設：主査だけが機械工学部必須。副査は制限なし** |
| §15 の延長 | §15(2) 最大1/3 | **§15(6)** 最大1/3（条番号が移動） |
| §7(6) | Modul**handbuch** | Modul**katalog** |

→ **差分の全体は `notes/regelwerk_gaps.md`。**
→ 現物は `_thesis_forms/03c_Pruefungsordnung_MPO2024_Aenderung.pdf`。

## ⚠️ 確認推奨（未確定）
1. 🟡 QIS の `Exam Regulations Version: 2017` 表記。2024 版は形式上「2015 版の改正」なので
   システム上のキーが古いままと思われるが、Müller 氏に1行で確認しておく。
2. これは Maschinenbau M.Sc. の規定。**国際/ダブルディグリー特別プログラム**だと別 PO の可能性。
3. 口頭試問の有無・形式は **Anlage 1.4（モジュール記述）＋ Modulhandbuch** で確定。
4. 英語執筆への **Junker の合意**（§7(7)、通常は黙示的）。

## この修論への当てはめ（現状）
- 英語のみ＝規定準拠 ✅／Selbstständigkeitserklärung 記載済 ✅／タイトル確定 ✅
- 2審査者＝**Junker（主査・IKM）・Wangenheim（副査・IDS）で確定**（2026-08-18 Soleimani メール）。
  日常指導は Soleimani＋Geisler（Betreuer 欄）。Wangenheim は Engineering Dynamics and Vibrations の試験官でもある
- 6か月期限・60 LP 要件・優等条件（総合&修論 ≤1.3）を意識して計画。
