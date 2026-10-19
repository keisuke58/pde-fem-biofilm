# 修論用に描き直した図 3.2・3.6（2026-10-19、EML 側セッション）

どちらも幅 16 cm（本文幅）で作ってある。`\includegraphics[width=\textwidth]` なら拡大縮小なし（文字は印刷時に目盛り 8 pt・軸ラベル 9 pt・列見出し 10 pt、Times = Nimbus Roman を埋め込み）。
上のタイトルとパネル内の凡例は外し、凡例は図の下に 1 つ。縦軸は 0〜1、目盛り 0・0.5・1。

| ファイル | 大きさ | 元（Tmcmc202601） |
|---|---|---|
| `figures/eml_seed_maps.pdf`（図 3.6） | 160 × 180 mm | `tools/plot_seed_maps.py --thesis` |
| `figures/eml_fit.pdf`（図 3.2） | 160 × 185 mm | `tools/fig2_v2.py` の `make(..., thesis=True)` |

## 図 3.6 のキャプション用: seed ごとの max logL と組成 RMSE

- max logL は最終段の粒子の最大値（polish 前）。RMSE は seed の MAP の生きた菌の組成（φψ/Σφψ）と、データの中央値の差。
- 段は CS・CH・DS が段 (iv)、DH は段 (iii)（段 (iv) は seed の一致の基準 (e) を満たさないため）。

| 条件 | seed 42 | seed 123 | seed 7 | 事後の中央値の RMSE |
|---|---|---|---|---|
| CS (stage iv) | −32.91 / 0.092 | −32.86 / 0.090 | −33.08 / 0.093 | 0.094 |
| CH (stage iv) | −26.84 / 0.060 | −26.52 / 0.060 | −27.11 / 0.066 | 0.086 |
| DS (stage iv) | −96.80 / 0.057 | −96.84 / 0.057 | −96.94 / 0.057 | 0.057 |
| DH (stage iii) | −51.41 / 0.073 | −51.54 / 0.076 | −51.77 / 0.073 | 0.083 |

（各セル = max logL / RMSE）

キャプション例（英語）:
> MAP estimate of each seed (coloured lines; 42, 123, 7), 10–90% band (grey) and median (dashed) of the posterior prediction pooled over the three seeds, and measured viable composition (dots). Stage (iv) for CS, CH and DS; stage (iii) for DH. The seeds' maximum log-likelihoods and composition RMSEs are listed in Table X / given in the text.

## 図 3.2 のキャプション用

列見出しの RMSE は外したので、キャプションか本文に書く。最終の MAP（3 seed の polish 後の最大）の RMSE は生きた菌 / 全菌 / 生きた割合で
CS 0.087 / 0.039 / 0.23、CH 0.059 / 0.078 / 0.12、DS 0.058 / 0.080 / 0.094、DH 0.073 / 0.159 / 0.07〜0.08（`EML_numbers_2026-10-19.md` と同じ）。
凡例: 実線 = MAP の生きた菌の組成 φψ/Σφψ、破線 = MAP の全菌の組成 φ/Σφ、帯 = 事後予測の 50%・80%、
塗りの点 = 生きた菌のデータ（中央値）、白抜き = 全菌のデータ（中央値）、箱 = 生きた菌の反復（四分位）。
