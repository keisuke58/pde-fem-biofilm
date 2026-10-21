# Draft (not sent): two observations on the element, to Oliver (cc Felix?)

Subject: Two observations on the element: neighbour selection and the phi time step

Hi Oliver,

while checking the coupled runs I found two things in the element that I wanted to ask you about. I have not
changed the code for the thesis runs.

1. Neighbour selection and symmetry

Results of mirror-symmetric set-ups are not symmetric. With only the diffusion of phi switched on (8^3, seed in the
centre), phi differs by up to 0.0085 between mirrored elements at T* = 0.25; in the Chapter 4 runs (16^3) by up to
0.9 % in x. As far as I can see, GetNhgb starts from a radius estimated from NEIGHBOR_CNT, enlarges it by 10 % until
there are enough points, and then takes exactly NEIGHBOR_CNT points in element order. The last distance shell is
therefore cut by the numbering at every Gauss point, and the neighbour sets of mirrored points are not mirrored.

| neighbours | max difference in phi, x | z | run time |
|---|---|---|---|
| first 30 in element order (deck value) | 8.5e-3 | 6.4e-3 | 3.1 min |
| first 120 | 4.6e-4 | 2.0e-3 | 7.1 min |
| first 400 | 2.8e-4 | 1.1e-3 | 13.1 min |
| nearest 30, a shell that straddles 30 left out (test build) | 8.6e-4 | 7.0e-4 | 2.7 min |

The mean phi in the seed changes by only 0.5 % between these, so for the thesis I keep 30 and state the asymmetry
as a limitation. Even the distance-sorted version keeps about 8e-4, so there seems to be a second, smaller source.

2. phi and the inner USSFin step

In USSFin the nutrient is solved with the inner step, but phi (Bio1..5, Bioloc) is advanced with the outer
increment at every inner step. With MY_USSFIN_DT = 0.1 and outer increments up to 0.1, as in all my Chapter 4 runs,
there is one inner step and this is exact. With a smaller inner step (I had set 0.0025 for the coupled runs) phi
advances (outer / inner) times too fast: in a pure diffusion test the seed mean at T* = 0.25 is 0.352 with ten inner
steps and 0.764 with one. I am rerunning the coupled runs with one inner step. Was the outer increment intended
there, or should it be the inner step?

Two questions in short:
1. Is the cut of the last shell in element order intended, and would a distance-sorted selection be in line with
   the operator?
2. Should phi be advanced with the inner step in USSFin?

(Felix: does your later version select the neighbours and advance phi in the same way? A yes or no is enough.)

Best regards,
Keisuke
