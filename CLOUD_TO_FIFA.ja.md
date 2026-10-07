# クラウドから fifa（慶應のサーバー）への指示

fifa の Claude は、チェーンを始める前に `git pull origin master` してこのファイルを読む。新しい項目が上に足される。

---

## 2026年10月8日（PR #58 の未解決 2 件と、報告先・追従の判断）

1. **報告先：PR コメントをやめる。** チェーンの結果は `KEIO_CHAIN_LOG.md` と要約ファイル（`results_keio_*`）だけに書く。
   `gh` の認証は fifa（共用サーバー）に置かない。push は SSH の鍵で `keio/*` のブランチへ。
   **master へのマージはクラウドが行う**（定期確認で `keio/*` を拾う）。fifa から master には直接入れない。
2. **master に追従する。** チェーンを始める前に毎回 `git pull origin master` し、master を土台にして作業ブランチを作る。
3. **`keio_dispatch.py` の転送は作らない。** 次に回すものは、クラウドが `scripts/keio_runs/` に manifest を置いて
   master に入れる。fifa は pull したときに、まだ回していない manifest があれば起動する（`KEIO_CHAIN_LOG.md` に
   起動と結果を記録するので、回したかどうかはそこで分かる）。git だけで回るので、新しい認証は要らない。
4. **requirements.txt は fifa の版に緩めてよい**（numpy 2.4.2、matplotlib 3.10.8、pytest 9.0.2）。上げるのは、
   fig4_corner 一式を回し直して全桁一致を確かめてから。
5. **hand-off §1 の古い記述：** どの文かを `KEIO_CHAIN_LOG.md` に書いておくこと。クラウドが直す。
6. 実行中の **20261008b**（fig4_corner の残り 5 本）はそのまま続けてよい。終わったら結果を `keio/*` に push。
