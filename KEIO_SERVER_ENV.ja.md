# 慶應のサーバー（fifa）の環境（2026年10月8日に確かめた）

`KEIO_SERVER_HANDOFF.ja.md` の §3 を実際に動かして分かったことの記録。clone した直後に
`scripts/setup_keio_server.sh` を走らせたあと、この文書を読む。

ホスト名は `fifa`（`fifawc`、172.17.36.98、ユーザー `nishioka`）。**共用のサーバー**で、
ほかに machida・kikkawa・blanloeil・murata が VS Code Server を動かしている。重い実行を
始める前に `uptime` と `nproc` を見る。

---

## 1. 手元の環境

| もの | 中身 |
|---|---|
| CPU | 12 コア |
| メモリ | 62 GiB |
| Abaqus | **2024**（`~/SIMULIA/EstProducts/2024`、Site ID 200000000010029、RELr426） |
| Fortran | `ifort` 2021.9.0（oneAPI 2023.1.0）。`ifx` と `gfortran` もある |
| Python | `~/IKM_Hiwi/.venv_jax/bin/python3`（3.12.0）が `python3`。`.bashrc` の PATH で入る |
| JAX | 0.9.0.1（CPU）。`requirements.txt` の注記どおり別に入っている |
| ANSYS | **ない**（`nvidia-smi` もない。GPU なし） |

`abaqus verify -user_std` は PASS（2026年10月8日）。ユーザーサブルーチンは Linux でも通る。

## 2. ディスク（**重要**）

```
/home   148T  使用 65T  空き 83T   44%   <- ここで作業する
/       160G  使用 150G 空き 11G   94%   <- /tmp がここ。危ない
```

Abaqus の scratch は既定で `/tmp`（`-tmpdir /tmp/nishioka_<job>_<pid>`）だった。`/` の空きは
11 GB しかなく、**約15人が使う共用のサーバー**なので埋めると全員が止まる。20³ の 1 本で 132 MB
だったので、`KEIO_PLAN` の 40³・48³ だと 1 本で 2 GB 近くになる。

**2026年10月8日に `~/abaqus_v6.env` を作って `/home` に移した**（`scratch =
"/home/nishioka/abaqus_work/scratch"`）。`abaqus information=environment` で確認できる。
これは私の Abaqus の実行すべてに効く。

`WORKROOT` は既定の `$HOME/abaqus_work` でよい（`/home` に 83T ある）。

`abaqus_composition/run_comp.sh` も、念のため scratch と `TMPDIR` を `$WORKROOT/scratch` に明示している
（`~/abaqus_v6.env` がない別のアカウントや別のマシンでも `/tmp` を使わないため。場所は `ABQ_SCRATCH` で変えられる）。
どちらも同じ `/home/nishioka/abaqus_work/scratch` を指すので食い違いはない。

## 3. 並列（**間違えると静かに壊れる**）

`abaqus information=environment` の既定は **`mp_mode=MPI`**。UMAT・UMATHT・UEL は栄養の場を
Fortran の module で共有しているので、MPI のプロセスに分けると共有できない。
`run_comp.sh` は `mp_mode=threads` を明示している。**自分で `abaqus` を直接呼ぶときは必ず付ける。**

## 4. Windows と同じ数字が出る理由

`wp2_imp_free_ph01` を Linux で回し、IKMHIWI03（Windows）の
`abaqus_composition/results_1007/wp2_imp_free_ph01.txt` と**全桁一致**した
（`results_keio_1008/wp2_imp_free_ph01.txt`）。hand-off §4 は「丸めの程度の差」を見込んでいたが
差はゼロ。理由は Abaqus が `compile_fortran` を再現性のある設定で固定しているため：

```
ifort -fp-model precise -no-fma -fimf-arch-consistency=true -prec-div -prec-sqrt
      -fp-speculation=safe -fprotect-parens
```

コンパイラが変わっても結果が動かないので、今後 Windows 側と食い違いが出たら
**丸めではなく実装の違い**を疑う。

## 5. パッケージの版（`requirements.txt` と違う）

| もの | requirements.txt | fifa に入っている |
|---|---|---|
| numpy | 2.4.6 | **2.4.2** |
| scipy | 1.17.1 | 1.17.1 |
| matplotlib | 3.11.0 | **3.10.8** |
| openpyxl | 3.1.5 | 3.1.5 |
| pytest | 9.1.1 | **9.0.2** |

`setup_keio_server.sh` は import できるかだけを見るので、この食い違いは warn にならない。
上の全桁一致は numpy 2.4.2 で出ているので、**版を上げると一致が崩れるかもしれない**。
揃えるかどうかは未決（§6）。

## 6. 走らせ方（§3 の実例）

```
export WORKROOT=$HOME/abaqus_work
cd ~/IKM_Hiwi/pde-fem-biofilm

# 1本 = 入力を作る -> run_comp.sh -> 比べる
python3 abaqus_composition/make_implant_inp.py abaqus_composition/<job>.inp <...>
abaqus_composition/run_comp.sh abaqus_composition/<job>.inp 2sp_case6 4
python3 abaqus_composition/summarize_implant.py \
    $WORKROOT/comp_<job>/<job>.dat abaqus_composition/<job>.json > <job>.txt
```

気をつけること：

- `summarize_implant.py` の**第2引数は入力**（`make_*_inp.py` が書くメッシュ情報の JSON）。
  出力先ではない。表は stdout に出るので、リダイレクトして `results_*/` に置く。
  hand-off §3 のコマンド列はここが読み取りにくい。
- `.inp` は `.gitignore` の `*.inp` で除外。メッシュ情報の `.json` も結果ではないので
  コミットしない。git に入れるのは要約の `.txt` だけ（`results_1006/`・`results_1007/` と同じ）。
- 長い実行は `setsid nohup ... &`。ジョブ管理システムはない。
- `.bashrc` が `abaqus` を shell 関数で包んでいて、`$HOME` 直下からの実行を拒否する。
  必ずプロジェクトのディレクトリに `cd` してから呼ぶ。
- **ほかの実行と重ねると壁時計時間は比べられない。** 数字は決まるが時間は伸びる。

## 7. 消費の項（gφ か gφc か）

`make_klempt_inp.py` のフラグの意味を間違えやすい：

| フラグ | 意味 |
|---|---|
| `--first-order` | **消費 gφc（1次）。** これがないと論文の gφ（0次） |
| `--felix` | 消費ではなく**前線項**の実装差（第一著者の現行実装：正のところだけ加算、φ<1 で (1−φ) 倍） |
| `--pen` | φ を [0, 1] に抑える剛性ペナルティ。`--pen 0` で上限なし |

`KEIO_PLAN.ja.md` §3 の「4.1 は gφc でしか論文に合わない」（RMS φ 0.066 / c 0.025）は
`results_1006/fx41_fo_p0.txt` = `--first-order --pen 0`、`--felix` なし。
`--felix` を付けた `ff41_fo` は RMS φ 0.252 で、同じ gφc でも一致しない。

## 8. テスト（2026年10月8日）

```
python3 -m pytest tests/ -q
473 件すべて通る（うち 3 件 xfail）、14 分ほど
```

55 個のテストファイルが fifa でそのまま通る。JAX を使うものも含めて、Python 側は
手を入れずに動く。

## 9. まだ確かめていないこと

- `run_comp.sh` 以外の経路（`make_cube_inp.py` + `compare_ansys.py`、`compare_paper.py`）は
  2026年10月8日の時点で Linux で未実行。
- 40³・48³ の大きいメッシュは未実行（§2 の scratch を先に片付けること）。
- `keio_wp2/`（軸対称の Python 試作）は未実行。

## 10. 自動で回す（`scripts/run_chain_keio.py`）

`ansys_usermat/apdl/run_chain.ps1`（IKMHIWI03 の ANSYS 用）の Linux/Abaqus 版。
1本ずつ回して要約を作り、基準と比べ、要約だけを commit して push し、PR にコメントを書く。
**次に何を回すかの判断は PR で行う**（この文書を書いた時点では、クラウドの Claude セッションが
定期的に PR を読んで決める構成にしている）。

```
python3 scripts/run_chain_keio.py --name 1008b \
    --runs scripts/keio_runs/1008b_fig4_corner.json --pr 58 --push
```

呼ぶと即座に戻り、PID とログの場所（`$WORKROOT/_chain_<name>.log`）を出す。
`setsid` で切り離すので、端末を閉じてもセッションが落ちても止まらない。

**共用サーバー向けのガード**（約15人が使うので意図的に遠慮している）：

| 既定 | 中身 |
|---|---|
| `--max-cpus 4` | 12 コアの 1/3。`--cpus` はこれで頭打ちになる |
| `--nice 10` | Abaqus を `nice` 経由で起動 |
| `--max-load 6.0` | 1分 load がこれ以上なら2分ごとに待つ。`--load-timeout 120` 分で諦める |
| `--min-free-root 8` | `/` の空きが 8 GiB を切ったら**走らせない**（`/` が埋まると全員が止まる） |
| `--min-free-work 100` | 作業ディスクの空きが 100 GiB 未満なら走らせない |
| 逐次のみ | 並列実行はしない |
| 後片付け | 要約を書いたら `.odb`・`.sim`・scratch・`/tmp` の残りを消す |

マニフェスト（`scripts/keio_runs/*.json`）の `{inp}`・`{dat}`・`{json}` は、入力ファイル、
終わったジョブの `.dat`、`make_*_inp.py` が書くメッシュ情報 JSON に置き換わる。
`reference` を省くと比較せずに記録だけする。

`--dry-run` で入力生成だけを試せる（マニフェストの検算に使う）。
