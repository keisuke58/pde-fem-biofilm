# クラウドから fifa（慶應のサーバー）への指示

fifa の Claude は、チェーンを始める前に `git pull origin master` してこのファイルを読む。新しい項目が上に足される。

---

## 2026年10月8日 夜（lab の設定が済んだ）

- ユーザーが fifa で `git remote add lab ...`、`git ls-remote lab`、`git push lab master:main` を実行し、通った
  （lab の `main` ができた）。SSH の鍵と組織の書き込み権限は足りている。
- fifa の Claude は `mmc-research-group` への操作が権限で止められているので、**lab への push はユーザーが手で行う**。
  fifa の Claude は、master を pull したあとで lab が遅れていれば `KEIO_CHAIN_LOG.md` に「`git push lab master:main` が必要」と
  書くだけでよい（自分では実行しない）。

---

## 2026年10月8日（村松研の git ができた）

ユーザーが村松研の組織に **`mmc-research-group/nishioka-biofilm-fem`**（private、空、既定ブランチ `main`）を作った。
慶應の作業はここにも置く。

1. **remote は `lab` という名前で足す。`origin`（公開の pde-fem-biofilm）は変えない。** GitHub の画面が出す
   `git remote add origin ...` はそのまま打たないこと（origin が上書きされ、チェーンの push 先が変わる）。
   ```
   git remote add lab git@github.com:mmc-research-group/nishioka-biofilm-fem.git
   git push lab master:main
   ```
2. **SSH の鍵で push する**（HTTPS とトークンは共用サーバーに置かない）。fifa の公開鍵（`~/.ssh/id_ed25519.pub`
   など）が、ユーザーの GitHub アカウントに登録されていて、そのアカウントが組織 `mmc-research-group` で書き込み
   できるかを `ssh -T git@github.com` と `git ls-remote lab` で確かめる。通らなければ、何が足りないかを
   `KEIO_CHAIN_LOG.md` に書く（鍵の登録はユーザーがする）。
3. **当面の運用：** チェーンの結果は今までどおり `origin` の `keio/*` に push し、master へのマージはクラウドが
   する。fifa は master を pull したら、`git push lab master:main` で lab の `main` を master に合わせる。
   lab だけに置くもの（未発表の原稿など）の分け方は、1月の村松先生との面談で決める。
4. **どちらの remote にも入れないもの：** Felix のコードと理論、オリバーのソース（`CLAUDE.md`）。

## 2026年10月8日 朝（ac94d34 への返事）

- 1・4・5 の対応を確認した。master に入れた（`run_comp.sh` は master の scratch の扱いと fifa の seed の
  コピーの両方が入った形で自動マージされた。チェーン 20261008b が終わってから pull すればよい）。
- hand-off §1 の作業ブランチの記述は、master ではすでに直してある（「作業ブランチ：master」）。
- `~/.config/gh/hosts.yml` の gh の認証は、ユーザーに `gh auth logout` を頼んだ。チェーンは使わないので急ぎではない。

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
