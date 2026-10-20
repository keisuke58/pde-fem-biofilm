# Draft (not sent): the neighbour cap in the NEM operator, to Oliver (cc Felix?)

Subject: Neighbour selection in the NEM element and a small asymmetry

Hi Oliver,

while checking the coupled runs I found that results of symmetric set-ups are not symmetric. The seed, the nutrient
layer and the corner blocks are mirror-symmetric in x and z, but on the 8^3 mesh phi differs by up to 0.0085 between
mirrored elements at T* = 0.25 with only the diffusion of phi switched on (no consumption, no front term). The front
term then amplifies it to about 0.09 by T* = 0.75. In the Chapter 4 runs (16^3) the difference is up to 0.9 % in x.

As far as I can see, the cause is the neighbour search (GetNhgb in NEM_UserData_P21_V05.F): it goes through the
points in element order and keeps the first NEIGHBOR_CNT points found within the search radius. Near the surface
there are more than 30 points in the radius, so which ones are dropped depends on the numbering, and the neighbour
sets of mirrored points are not mirrored. Changing only the deck value:

| NEIGHBOR_CNT | max difference in phi, x | z | run time |
|---|---|---|---|
| 30  | 8.5e-3 | 6.4e-3 | 3.1 min |
| 60  | 2.5e-3 | 4.2e-3 | 3.4 min |
| 120 | 4.6e-4 | 2.0e-3 | 7.1 min |

A larger cap removes most of the x part; some z part remains, so there may be a second ordering effect (perhaps the
surface points or the extended list).

For the thesis I keep NEIGHBOR_CNT = 30, as in all earlier runs, and add one 16^3 run with 120 to show how much the
main numbers change. Two questions:
1. Is the cap meant to be reached near the surface, or should it be larger than the number of points in the radius?
2. Would sorting the points in the radius by distance before taking them be in line with how the operator is meant
   to work? I have not changed the code.

(Felix: does your later version select the neighbours in the same way? No need to send code, a yes or no is enough.)

Best regards,
Keisuke
