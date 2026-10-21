#!/usr/bin/env python3
"""Build the 20 um version of the thesis (the user, 22 Oct 2026): main_um.tex + chapters_um/*.tex.

The computations use the element's example cube of 2 mm (deck in /UNITS,MPA); all constants were converted from
Klempt et al. 2024, Table 2 (20 um cube) by keeping beta/L^2, d/L^2, r/L and g, so the equations are the same. This
script writes the same text with every length on the paper's 20 um cube (x 1/100: mm -> 10 um, mm^2 -> 100 um^2);
the numbers of beta, d and r then are the paper's own. The mm version (main.tex) stays as it is.

    python thesis/chapters_um/make_um.py   (then build main_um.tex)
"""
from pathlib import Path

T = Path(__file__).resolve().parents[1]
U = r"$~$\mu$m"
SUBS = {
    "ch4_ansys.tex": [
        ("because the width $\\beta^*=0.2$~mm of the Gaussian weight is fixed in\nmillimetres.",
         "because the width $\\beta^*=2" + U + " of the Gaussian weight is a fixed\nlength."),
        ("$d=10^{10}~\\mu$m$^2/T^*$\non their $20~\\mu$m cube, converted like $\\beta$ to $d=10^8$~mm$^2/T^*$ on the $2$~mm cube, the nutrient adjusts within",
         "$d=10^{10}~\\mu$m$^2/T^*$\non the $20~\\mu$m cube used here, the nutrient adjusts within"),
        ("$\\ell=1$~mm. The antibiotic", "$\\ell=10" + U + ". The antibiotic"),
        ("with $d=1$~mm$^2/T^*$ the living", "with $d=100" + U + "$^2/T^*$ the living"),
        ("The model is a 2~mm cube instead of the paper's 20~$\\mu$m one. For the\ndiffusion coefficient of~\\eqref{eq:ch5-field} I use\n$\\beta=0.02$~mm$^2/T^*$, the value $\\beta=2$ of Klempt et al.\\ 2024, Table~2,\nconverted to the 2~mm cube. The conversion is an assumption: the paper gives\nno length and time scale, and its first author describes $\\beta$ as chosen\nfor sensible results and dependent on the mesh. The framework's example\ninput, $\\beta=10^{-4}$~mm$^2/T^*$, is reported as a sensitivity case.",
         "The model is the paper's cube of 20~$\\mu$m edge length. The framework's example model is the same cube\n"
         "written in millimetres (\\texttt{/UNITS,MPA}, edge 2~mm); its runs are reported here with every length divided\n"
         "by 100 and the constants scaled with it ($\\beta$ and $d$ by $10^{-4}$, $r$ by $10^{-2}$), which leaves the\n"
         "equations and the stresses unchanged. Axes of figures made from these runs are still in the deck's millimetres\n"
         "(1~mm there is 10~$\\mu$m here). For the\ndiffusion coefficient of~\\eqref{eq:ch5-field} I use\n"
         "$\\beta=2" + U + "$^2/T^*$, the value of Klempt et al.\\ 2024, Table~2.\n"
         "The paper gives\nno length and time scale, and its first author describes $\\beta$ as chosen\nfor sensible results and dependent on the mesh. The framework's example\n"
         "input, $\\beta=10^{-2}" + U + "$^2/T^*$, is reported as a sensitivity case."),
        ("(\\texttt{/UNITS,MPA}, lengths in mm), not values from the paper: $r=100$~mm$/T^*$,\n$k=0.1$, $d=1$~mm$^2/T^*$",
         "(scaled to the 20~$\\mu$m cube), not values from the paper: $r=1000" + U + "$/T^*$,\n$k=0.1$, $d=100" + U + "$^2/T^*$"),
        ("$0.4$~mm from the centre on $8^3$", "$4" + U + " from the centre on $8^3$"),
        ("and $0.4$~mm from it;", "and $4" + U + " from it;"),
        ("From $\\beta=0.005$ to $0.1$~mm$^2/T^*$ the mean", "From $\\beta=0.5$ to $10" + U + "$^2/T^*$ the mean"),
        ("For $\\beta\\ge0.02$ the diffusion\nlength $\\sqrt{\\beta T^*}\\ge0.14$~mm is resolved",
         "For $\\beta\\ge2" + U + "$^2/T^*$ the diffusion\nlength $\\sqrt{\\beta T^*}\\ge1.4" + U + " is resolved"),
        ("at $\\beta=0.02$ and $0.05$ and within $16\\,\\%$ at\n$\\beta=0.1$.", "at $\\beta=2$ and $5$ and within $16\\,\\%$ at\n$\\beta=10$."),
        ("$\\beta=10^{-4}$~mm$^2/T^*$, the layer is $0.01$~mm thick", "$\\beta=10^{-2}" + U + "$^2/T^*$, the layer is $0.1" + U + " thick"),
        ("for different $\\beta$ [Pa]; $\\beta=0.02$ is the value", "for different $\\beta$ [Pa]; $\\beta=2" + U + "$^2/T^*$ is the value"),
        ("$16^3$ at $\\beta=0.1$ with", "$16^3$ at $\\beta=10$ with"),
        ("$\\beta$ [mm$^2/T^*$] & $0.005$ & $0.01$ & $0.02$ & $0.05$ & $0.1$ \\\\",
         "$\\beta$ [$\\mu$m$^2/T^*$] & $0.5$ & $1$ & $2$ & $5$ & $10$ \\\\"),
        ("Table~\\ref{tab:ch5-beta} for $\\beta=0.05$ and $0.1$ are", "Table~\\ref{tab:ch5-beta} for $\\beta=5$ and $10$ are"),
        ("example input $\\beta=10^{-4}$, case~6", "example input $\\beta=10^{-2}" + U + "$^2/T^*$, case~6"),
        ("$y=-1$~mm and consumed in the seed; its strength", "$y=-10" + U + " and consumed in the seed; its strength"),
        ("$\\Lambda^2=g\\ell^2/d$ with $\\ell=1$~mm,", "$\\Lambda^2=g\\ell^2/d$ with $\\ell=10" + U + ","),
        ("Table~2, converted to the cube like $\\beta$ (an assumption), is weaker", "Table~2, is weaker"),
        ("held at $y=-1$~mm\n", "held at $y=-10" + U + "\n"),
        ("$\\beta=10^{-4}$~mm$^2/T^*$ they do not converge", "$\\beta=10^{-2}" + U + "$^2/T^*$ they do not converge"),
        ("~mm$^2/T^*$", "GENERIC"),   # the remaining beta = 0.02 mm^2/T* (below)
    ],
    "ch6_conclusion.tex": [("$\\beta=0.02$~mm$^2/T^*$", "$\\beta=2" + U + "$^2/T^*$")],
    "appendix_klempt2024.tex": [("Chapter~\\ref{ch:ansys}, Table~2 converted to the 2~mm cube, is\n        therefore an assumption.",
                                 "Chapter~\\ref{ch:ansys}, the value of Table~2, is\n        therefore an assumption.")],
    "appendix_stability.tex": [("$\\beta^*=0.2$~mm. For", "$\\beta^*=2" + U + ". For"),
                               ("$\\beta^*$ is a length in millimetres,", "$\\beta^*$ is a fixed length,"),
                               ("    $\\beta^*=0.2$~mm; shaded", "    $\\beta^*=2" + U + "; shaded"),
                               ("$\\beta=0.1$~mm$^2/T^*$ to $T^*=1.1$", "$\\beta=10" + U + "$^2/T^*$ to $T^*=1.1$"),
                               ("$\\beta^*=0.2$~mm as in the element", "$\\beta^*=2" + U + " as in the element")],
}
for name, subs in SUBS.items():
    s = (T / "chapters" / name).read_text()
    for old, new in subs:
        if new == "GENERIC":
            s = s.replace("$\\beta=0.02$~mm$^2/T^*$", "$\\beta=2" + U + "$^2/T^*$")
            continue
        n = s.count(old)
        if n == 0:
            raise SystemExit(f"{name}: not found: {old[:70]!r}")
        s = s.replace(old, new)
    (T / "chapters_um" / name).write_text(s)
m = (T / "main.tex").read_text()
for name in SUBS:
    m = m.replace("\\input{chapters/" + name[:-4] + "}", "\\input{chapters_um/" + name[:-4] + "}")
m = m.replace("width [mm] of the weight", "width [$\\mu$m] of the weight")
(T / "main_um.tex").write_text(m)
print("wrote main_um.tex and", ", ".join(SUBS))
