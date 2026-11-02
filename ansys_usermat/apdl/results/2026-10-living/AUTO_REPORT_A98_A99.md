# Automatic report from IKMHIWI03 (A-98, A-99), written by F:\auto_night\day20.ps1

Written 2026-11-02 15:12. Not checked by hand; numbers copied from the run logs.

## A-98: 192^3 without antibiotic, dc/dt against the limit (seed averages, T* = 1)

| example | quantity | limit | with dc/dt | difference |
|---|---|---|---|---|
| separate | alpha-1 | 9.135281e-04 | 9.135430e-04 | +0.002 % |
| separate | vM [Pa] | 6.571409e-04 | 6.571118e-04 | -0.004 % |
| separate | living_seed | 9.647470e-01 | 9.647469e-01 | -0.000 % |
| separate | share_a_seed | 9.780549e-01 | 9.780549e-01 | -0.000 % |
| separate | c_seed | 5.405786e-01 | 5.405761e-01 | -0.000 % |
| separate | phi_3D max (domain, all t) | 1.000199e+00 | 1.000199e+00 | +0.000 % |
| mixed (equal) | alpha-1 | 8.452865e-04 | 8.453000e-04 | +0.002 % |
| mixed (equal) | vM [Pa] | 3.894018e-04 | 3.894270e-04 | +0.006 % |
| mixed (equal) | living_seed | 9.608212e-01 | 9.608247e-01 | +0.000 % |
| mixed (equal) | share_a_seed | 8.548060e-01 | 8.549427e-01 | +0.016 % |
| mixed (equal) | c_seed | 6.063786e-01 | 6.063755e-01 | -0.001 % |
| mixed (equal) | phi_3D max (domain, all t) | 1.000155e+00 | 1.000155e+00 | +0.000 % |

day18 log:

```
2026-10-30 11:31:45 start (C: free 6.8 GB)
2026-10-30 11:34:27 start (C: free 6.8 GB)
2026-10-30 11:36:44 start (C: free 6.8 GB)
2026-10-30 14:53:08 nf 192 --only separate: [separate, antibiotic False] 11757 s; colonised 1769472 cells; seed: living 0.965, share A 0.978, c 0.541 (exit 0) (196 min; C: free 6.0 GB)
2026-11-02 15:09:09 nf 192 --mix-equal --only mixed: [mixed, antibiotic False] 260134 s; colonised 1769472 cells; seed: living 0.961, share A 0.855, c 0.606 (exit 0) (4,336 min; C: free 5.9 GB)
2026-11-02 15:09:09 end
```

## A-99: mirrored walls, wider weights at edges (x1.5) and corners (x2)

Manufactured solution: L2 errors at T* = 0.05 (references: mt8_8 phi 9.47e-4, c 8.95e-4; mb16_16 3.37e-4, 2.15e-4).
mb24_24 = the A-70 operator (factors 1) at 24^3; me = widened. mf16_16_sep_ab: the coupled 16^3 example (mu16 settings).

day19 log:

```
2026-10-30 12:42:07 waiting for F:\auto_night\day18.log
2026-11-02 15:09:31 start
2026-11-02 15:09:54 me8_8 (0 min): 2026-11-02 15:09:54 DONE me8_8 exit -1, errors , 0.4 min / 2026-11-02 15:09:54   me8_8 FAILED, not exported
2026-11-02 15:10:16 me16_16 (0 min): 2026-11-02 15:10:16 DONE me16_16 exit -1, errors , 0.3 min / 2026-11-02 15:10:16   me16_16 FAILED, not exported
2026-11-02 15:10:38 mf16_16_sep_ab (0 min): 2026-11-02 15:10:38 DONE mf16_16_sep_ab exit -1, errors , 0.3 min / 2026-11-02 15:10:38   mf16_16_sep_ab FAILED, not exported
2026-11-02 15:11:00 me24_24 (0 min): 2026-11-02 15:11:00 DONE me24_24 exit -1, errors , 0.3 min / 2026-11-02 15:11:00   me24_24 FAILED, not exported
2026-11-02 15:11:22 mb24_24 (0 min): 2026-11-02 15:11:22 DONE mb24_24 exit -1, errors , 0.3 min / 2026-11-02 15:11:22   mb24_24 FAILED, not exported
2026-11-02 15:11:22 end
```

run_living log of me8_8:

```
2026-11-02 15:09:32 start: me8_8:8:--layout mix --front 1 --v2 --start-step 0.00025 --dt-ussfin 0.1 --cons 1 --surf-nearest --lap-file --mms --deltim 0.00025 --tend 0.05 --eco eco_u3_a100.txt --mirror --artdiff 8
2026-11-02 15:09:33   deck: artdiff: MY_BETA1 0.02 -> 0.645 (h = 0.25 mm, r = 10)
2026-11-02 15:09:33   deck: ! make_living_deck.py: v2, layout mix, antibiotic False, front r = 10 mm/T*, L = 1 False, old consumption False, g = 1, g_a 0.25, phi_min 2, start step 0.00025, outer dt 0.00025, USSFin dt 0.1, end T* 0.05, quasi-static c; BIOFILM1 32 elements, NUTRIENT2 (all outer faces) 296 elements
2026-11-02 15:09:33 START me8_8 (8^3: --layout mix --front 1 --v2 --start-step 0.00025 --dt-ussfin 0.1 --cons 1 --surf-nearest --lap-file --mms --deltim 0.00025 --tend 0.05  --mirror --artdiff 8; point model eco_u3_a100.txt)
2026-11-02 15:09:54 DONE me8_8 exit -1, errors , 0.4 min
2026-11-02 15:09:54   me8_8 FAILED, not exported
2026-11-02 15:09:54 end
```

run_living log of me16_16:

```
2026-11-02 15:09:55 start: me16_16:16:--layout mix --front 1 --v2 --start-step 0.00025 --dt-ussfin 0.1 --cons 1 --surf-nearest --lap-file --mms --deltim 0.00025 --tend 0.05 --eco eco_u3_a100.txt --mirror --artdiff 16 --trace-elems 1912,8
2026-11-02 15:09:55   deck: artdiff: MY_BETA1 0.02 -> 0.3325 (h = 0.125 mm, r = 10)
2026-11-02 15:09:55   deck: ! make_living_deck.py: v2, layout mix, antibiotic False, front r = 10 mm/T*, L = 1 False, old consumption False, g = 1, g_a 0.25, phi_min 2, start step 0.00025, outer dt 0.00025, USSFin dt 0.1, end T* 0.05, quasi-static c; BIOFILM1 256 elements, NUTRIENT2 (all outer faces) 1352 elements
2026-11-02 15:09:55 START me16_16 (16^3: --layout mix --front 1 --v2 --start-step 0.00025 --dt-ussfin 0.1 --cons 1 --surf-nearest --lap-file --mms --deltim 0.00025 --tend 0.05  --mirror --artdiff 16; point model eco_u3_a100.txt)
2026-11-02 15:10:16 DONE me16_16 exit -1, errors , 0.3 min
2026-11-02 15:10:16   me16_16 FAILED, not exported
2026-11-02 15:10:16 end
```

run_living log of mf16_16_sep_ab:

```
2026-11-02 15:10:16 start: mf16_16_sep_ab:16:--layout sep --ab --t-ab 0.333333 --front 1 --v2 --start-step 0.0025 --dt-ussfin 0.1 --cons 1 --cons-ab 1 --phi-min 0.001 --transient-c 1e8 --nut1-init-all --surf-nearest --lap-file --eco eco_u3_a100.txt --trace-elems 1,1912 --deltim 0.00125 --front-upwind --mirror
2026-11-02 15:10:17   deck: ! make_living_deck.py: v2, layout sep, antibiotic True from T* 0.333333, front r = 10 mm/T*, L = 1 False, old consumption False, g = 1, g_a 100000000.0, phi_min 0.001, start step 0.0025, outer dt 0.00125, USSFin dt 0.1, end T* base, dc/dt with d = 1e+08; BIOFILM1 384 elements, NUTRIENT2 (all outer faces) 1352 elements
2026-11-02 15:10:17 START mf16_16_sep_ab (16^3: --layout sep --ab --t-ab 0.333333 --front 1 --v2 --start-step 0.0025 --dt-ussfin 0.1 --cons 1 --cons-ab 1 --phi-min 0.001 --transient-c 1e8 --nut1-init-all --surf-nearest --lap-file   --deltim 0.00125 --front-upwind --mirror; point model eco_u3_a100.txt)
2026-11-02 15:10:38 DONE mf16_16_sep_ab exit -1, errors , 0.3 min
2026-11-02 15:10:38   mf16_16_sep_ab FAILED, not exported
2026-11-02 15:10:38 end
```

run_living log of me24_24:

```
2026-11-02 15:10:38 start: me24_24:24:--layout mix --front 1 --v2 --start-step 0.00025 --dt-ussfin 0.1 --cons 1 --surf-nearest --lap-file --mms --deltim 0.00025 --tend 0.05 --eco eco_u3_a100.txt --mirror --artdiff 24
2026-11-02 15:10:39   deck: artdiff: MY_BETA1 0.02 -> 0.228333 (h = 0.0833333 mm, r = 10)
2026-11-02 15:10:39   deck: ! make_living_deck.py: v2, layout mix, antibiotic False, front r = 10 mm/T*, L = 1 False, old consumption False, g = 1, g_a 0.25, phi_min 2, start step 0.00025, outer dt 0.00025, USSFin dt 0.1, end T* 0.05, quasi-static c; BIOFILM1 864 elements, NUTRIENT2 (all outer faces) 3176 elements
2026-11-02 15:10:39 START me24_24 (24^3: --layout mix --front 1 --v2 --start-step 0.00025 --dt-ussfin 0.1 --cons 1 --surf-nearest --lap-file --mms --deltim 0.00025 --tend 0.05  --mirror --artdiff 24; point model eco_u3_a100.txt)
2026-11-02 15:11:00 DONE me24_24 exit -1, errors , 0.3 min
2026-11-02 15:11:00   me24_24 FAILED, not exported
2026-11-02 15:11:00 end
```

run_living log of mb24_24:

```
2026-11-02 15:11:00 start: mb24_24:24:--layout mix --front 1 --v2 --start-step 0.00025 --dt-ussfin 0.1 --cons 1 --surf-nearest --lap-file --mms --deltim 0.00025 --tend 0.05 --eco eco_u3_a100.txt --mirror --artdiff 24
2026-11-02 15:11:02   deck: artdiff: MY_BETA1 0.02 -> 0.228333 (h = 0.0833333 mm, r = 10)
2026-11-02 15:11:02   deck: ! make_living_deck.py: v2, layout mix, antibiotic False, front r = 10 mm/T*, L = 1 False, old consumption False, g = 1, g_a 0.25, phi_min 2, start step 0.00025, outer dt 0.00025, USSFin dt 0.1, end T* 0.05, quasi-static c; BIOFILM1 864 elements, NUTRIENT2 (all outer faces) 3176 elements
2026-11-02 15:11:02 START mb24_24 (24^3: --layout mix --front 1 --v2 --start-step 0.00025 --dt-ussfin 0.1 --cons 1 --surf-nearest --lap-file --mms --deltim 0.00025 --tend 0.05  --mirror --artdiff 24; point model eco_u3_a100.txt)
2026-11-02 15:11:22 DONE mb24_24 exit -1, errors , 0.3 min
2026-11-02 15:11:22   mb24_24 FAILED, not exported
2026-11-02 15:11:22 end
```

