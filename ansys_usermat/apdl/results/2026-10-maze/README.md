# Maze runs (Klempt et al. 2024, Sect. 4.4) with the agreed coupling

`maze_coupled_2d.py`, figure `fig_maze.py` -> `assets/fig_maze.pdf`. 2D, cell 2.5 um, maze 560 um across (own design:
perfect polar maze, seed 1, walls 10 um; the paper's maze is not given as data). Constants: Table 6 of Klempt 2024
(d = 1e13 um^2/T*, beta = 10 um^2/T*, k_alpha = 1e-6, r = 200 um/T*), the rest Table 2 (g = 1e8, k = 1); point model
Fritsch case I (c*_0 = 100, s = 0.1), alpha*_0 = 100, phi_min = 0.001; colony at the entrance (both species, 0.25/0.30),
nutrient held in the centre (r < 20 um).

- `maze_h2.5_tend8.json`: no antibiotic, T* = 8 (the front needs about T* = 7 for the ~500 um path at the paper's r;
  the paper gives its maze times only in per cent). The point model therefore runs over eight data windows (s = 0.1).
- `maze_h2.5_tend8_ab.json`: antibiotic held at the entrance from T* = 8/3, same diffusivity and uptake as the nutrient.
