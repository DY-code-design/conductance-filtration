# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig5.py
#   sha256(src) : 90d2ed2c358490fa
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 5. Computed on the fly; no result file is read.
import pathlib
import numpy as np
import matplotlib
import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle
from scipy.stats import kendalltau
import _data

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()
C_H, C_D, C_G, C_X = "#7b8794", "#3b6ea5", "#c0392b", "#e07b39"
NS = 14
rng = np.random.default_rng(3)
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6})


def panel_a(ax):
    P = {"A": (0.22, 0.74, 0.155), "B": (0.52, 0.82, 0.155),
         "C": (0.46, 0.44, 0.195), "D": (0.83, 0.58, 0.128),
         "E": (0.17, 0.27, 0.118), "F": (0.75, 0.19, 0.145)}
    for k, (x, y, r) in P.items():
        ax.add_patch(Circle((x, y), r, fc="#dfe6ee", ec="#4a5568", lw=0.9, zorder=2))
        ax.text(x, y, k, ha="center", va="center", fontsize=9.2,
                color="#2d3748", zorder=6,
                bbox=dict(fc="white", ec="none", pad=0.9, alpha=0.92))
    edges = [("A", "B", 1, "$e_1$", +1), ("B", "C", 1, "$e_2$", -1),
             ("C", "D", 0, "$e_3$", +1), ("C", "F", 0, "$e_4$", +1),
             ("A", "E", 0, "$e_5$", -1)]
    D_LAB = 0.090
    for a, b, contact, lab, s in edges:
        xa, ya, _ = P[a]; xb, yb, _ = P[b]
        kw = dict(color=C_G, lw=1.8, zorder=4)
        if not contact:
            kw["dashes"] = (3, 2)
        ax.plot([xa, xb], [ya, yb], **kw)
        dx, dy = xb - xa, yb - ya
        n = np.hypot(dx, dy)
        ux, uy = -dy / n, dx / n
        ax.text((xa + xb) / 2 + s * D_LAB * ux, (ya + yb) / 2 + s * D_LAB * uy, lab,
                ha="center", va="center", fontsize=9.0, color=C_G, zorder=5,
                bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.9))
    ax.plot([], [], "-", color=C_G, lw=1.6, label="contact branch")
    ax.plot([], [], "--", color=C_G, lw=1.6, label="near-contact branch")
    ax.legend(loc="lower left", frameon=False, fontsize=8.4, handlelength=1.6,
              borderpad=0.1, labelspacing=0.25, bbox_to_anchor=(-0.02, -0.02))
    ax.set_xlim(0, 1); ax.set_ylim(0.02, 1.02); ax.set_aspect("equal"); ax.axis("off")
    pass
    ax.text(0.5, -0.02, r"every edge carries $h$, $D$ and $g$ at the same time",
            transform=ax.transAxes, ha="center", va="top", fontsize=8.4, color="#4a5568")


def main():
    fig, ax = plt.subplots(figsize=(3.5, 3.4))
    panel_a(ax)
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig5_which_axis.{ext}", dpi=220, bbox_inches="tight")
    print("wrote", OUT / "fig5_which_axis.pdf")


if __name__ == "__main__":
    main()
