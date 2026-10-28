# IKMHIWI03 の git 未管理分のバックアップ計画（12 月末まで）

IKMHIWI03 は 2026 年 12 月末まで（CLAUDE.md）。git に入っていない大事なものを zip にして、ユーザーが別の場所に移す
（ユーザー、10 月 22 日）。F: の git 未管理分は全部で約 1 TB、大半は作り直せる（.rst / .odb 約 400 GB、点モデルの記録 csv、
行列の書き出し、ソルバーの一時ファイル）。

## やり方

```
.\scripts\backup_untracked_ikmhiwi03.ps1            # 下見：zip ごとのファイル数と大きさ
.\scripts\backup_untracked_ikmhiwi03.ps1 -Go -Rst   # F:\backup_2026-12 に zip を書く（-Rst：本番 u8_ の .rst も）
```

最後の計算が終わってから（12 月中旬）、`-Go` で作り、ユーザーが移す。下見の大きさ（22 日、圧縮前）：

| zip | 中身 | 大きさ（圧縮前） | 置き場所 |
|---|---|---|---|
| `private_felix.zip` | `F:\felix_private` そのまま（フェリックスのコードと博士論文の章） | 約 4.0 GB | **非公開のみ**（フェリックスの依頼、CLAUDE.md） |
| `private_element.zip` | 各ビルドのオリバーの要素のソース（.f .F .inc .h）、`F:\oliver_element_patches*` | 約 32 MB | **非公開のみ** |
| `work_text.zip` | 各作業フォルダのデッキ（.dat）、スクリプト、マクロ、点モデルの設定、保存形ラプラシアン（作った所のものだけ）、実行と連鎖の記録、自動処理（auto_*）、小さい JSON；Abaqus の living_gp.txt、.sta、.msg、.dat、.inp；10 月 20 日のバックアップ | 約 4 GB（テキストなので圧縮で 1 GB 前後の見込み） | どこでも（秘密なし。ただし公開しない方がよい） |
| `rst_main.zip`（任意） | 本番 u8_ の .rst | 約 4 GB（24³ が終われば増える） | どこでも |

入れないもの（作り直せる、または一時ファイル）：本番以外の .rst / .odb、comp_trace.csv、phi_trace.csv、front_dump、
nem_dump、20 MB を超える out_*.txt、ソルバーの一時ファイル（.esav .full .emat .r00* .nr* .page .osav .stt .mdl ...）。
`.env`（PAT）は入れない（秘密）。

## そのほか

- 結果の数字（fields_*.json、*_pm.json、Python の json）と、文書・スクリプト・パッチは git にある。
- ANSYS のライセンス、Intel Fortran、Abaqus はこの PC と大学のもの。慶應では Abaqus と Python だけ（CLAUDE.md）。
- 慶應で ANSYS を使えないので、オリバーの要素の再ビルドは IKMHIWI03 でしかできない。12 月までに要るビルドは作っておく。
