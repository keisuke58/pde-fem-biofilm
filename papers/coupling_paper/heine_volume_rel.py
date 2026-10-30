#!/usr/bin/env python3
"""heine_volume_rel.py -- biofilm volume of Heine et al. 2025 (Front. Oral Health 6:1649419, Fig. 2A) relative to day 1:
median per condition and day divided by the day-1 median, with a 90 % bootstrap band. The raw workbook
("Zusammenfassung HOBIC+Static.xlsx", from N. Heine) is not redistributed: keep it outside git and pass its path.

    python papers/coupling_paper/heine_volume_rel.py /path/to/Zusammenfassung_HOBIC+Static.xlsx [--out rel.json]

Sheets: HOBIC22 = CH, Static22 = CS, HOBIC24 = DH, Static24 = DS (checked against volume_long.csv of Tmcmc202601).
Day 1 of HOBIC22 is listed twice (40x objective, then "Umgerechnet auf 63er Objektiv"); the 63x-converted rows are used,
the objective of the other days.
"""
import json
import random
import statistics as st
import sys

import openpyxl

NAME = {"HOBIC22": "CH", "Static22": "CS", "HOBIC24": "DH", "Static24": "DS"}
DAYS = [1, 3, 6, 10, 15, 21]


def main():
    wb = openpyxl.load_workbook(sys.argv[1], data_only=True)
    out = {}
    for ws in wb:
        if ws.title not in NAME:
            continue
        rows = list(ws.iter_rows(min_row=3, values_only=True))
        cut = [i for i, r in enumerate(rows) if isinstance(r[0], str)]
        cols = []
        for j in range(len(DAYS)):
            v = [r[j] for i, r in enumerate(rows) if isinstance(r[j], (int, float)) and (j > 0 or not cut or i > cut[0])]
            cols.append(v)
        m0 = st.median(cols[0])
        random.seed(1)
        rec = {}
        for d, c in zip(DAYS, cols):
            bs = sorted(st.median(random.choices(c, k=len(c))) / st.median(random.choices(cols[0], k=len(cols[0])))
                        for _ in range(2000))
            rec[d] = dict(rel=st.median(c) / m0, lo=bs[100], hi=bs[1900], n=len(c))
        out[NAME[ws.title]] = rec
        print(NAME[ws.title], " ".join(f"d{d}:{r['rel']:.2f}[{r['lo']:.2f},{r['hi']:.2f}]" for d, r in rec.items()))
    if "--out" in sys.argv:
        json.dump(out, open(sys.argv[sys.argv.index("--out") + 1], "w"), indent=1)


if __name__ == "__main__":
    main()
