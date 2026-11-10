#!/usr/bin/env python3
"""Renders the display equations of README.md as SVG images: dark text on a white background.

GitHub draws README math as native MathML, whose look depends on the math fonts of each browser;
SVG images look the same everywhere. The white background keeps them readable whatever the theme:
a <picture> with a prefers-color-scheme source follows the operating system, not the GitHub theme,
and showed light text on a white page when the two differed.
Needs Matplotlib and a LaTeX installation (amsmath, amssymb).

  python3 docs/eq/make_equations.py            # all, from the repository root
  python3 docs/eq/make_equations.py width rstar   # only these
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"text.usetex": True, "text.latex.preamble": r"\usepackage{amsmath,amssymb}",
                            "svg.fonttype": "path", "font.size": 15})

EQUATIONS = {
    "width": r"w(x) = \min_{y \in O} d(x, y), \qquad x \in F",
    "free": r"V_r = \{\, x \in F : w(x) \ge r \,\}, \qquad r = w_{\max},\, \dots,\, 1",
    "signed": r"s(x) = \begin{cases} \phantom{-}w(x), & x \in F \\ -\min_{y \in F} d(x, y), & x \in O \end{cases}"
              r"\qquad V_r = \{\, x : s(x) \ge r \,\}",
    "betti": r"\beta_k(r) = \#\{\, \mathrm{bars\ of\ } H_k : r_d < r \le r_b \,\}",
    "rstar": r"r^{*} = \max\bigl\{\, r : \operatorname{rank}\bigl(H_1(V_r) \to H_1(\mathbb{T}^2)\bigr) = 2 \,\bigr\}",
    "entropy": r"E_k = -\sum_i p_i \ln p_i, \qquad p_i = \frac{\ell_i}{\sum_j \ell_j}, \quad \ell_i = r_{b,i} - r_{d,i}",
}
COLOR = "#1f2328"  # GitHub's light-theme text colour

out = Path(__file__).resolve().parent
selected = sys.argv[1:] or list(EQUATIONS)
for name in selected:
    tex = EQUATIONS[name]
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, rf"$\displaystyle {tex}$", color=COLOR)
    fig.savefig(out / f"{name}.svg", facecolor="white", bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    print(out / f"{name}.svg")
