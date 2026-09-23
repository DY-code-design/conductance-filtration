# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/_branchfig.py
#   sha256(src) : 1b31b399f6b76e75
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Shared drawing parts for Fig. 2-4. Reads pair_conductance; does not change it.
import pathlib
import numpy as np
import matplotlib
import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch
import _data
import L0_engine.thermal_sheaf_filtration as tsf

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()
C_F, C_M, C_I, C_T = "#dfe6ee", "#f6e7c9", "#c0392b", "#4a5568"
C_B = "#7e57c2"
CS = ["#08306b", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef"]
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})


def save(fig, stem):
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"{stem}.{ext}", dpi=220, bbox_inches="tight")
    print("wrote", OUT / f"{stem}.pdf")


def sphere(ax, x, y, r, lab=None):
    ax.add_patch(Circle((x, y), r, fc=C_F, ec="#4a5568", lw=1.0, zorder=3))
    if lab:
        ax.text(x, y, lab, ha="center", va="center", fontsize=7, zorder=4)

def wire(ax, x0, x1, y, c="#4a5568", lw=0.9):
    ax.plot([x0, x1], [y, y], "-", color=c, lw=lw, zorder=2)


def vwire(ax, x, y0, y1, c="#4a5568", lw=0.9):
    ax.plot([x, x], [y0, y1], "-", color=c, lw=lw, zorder=2)


def res(ax, xc, y, w, lab, fc="#ffffff", ec="#4a5568", fs=6.0, hh=0.028, lw=0.9):
    ax.add_patch(Rectangle((xc - w / 2, y - hh), w, 2 * hh, fc=fc, ec=ec,
                           lw=lw, zorder=4))
    ax.text(xc, y, lab, ha="center", va="center", fontsize=fs, zorder=5)


def node(ax, x, y, lab, side="left", fs=6.8):
    ax.plot([x], [y], "o", ms=5.2, mfc="#1f2937", mec="none", zorder=6)
    ax.text(x + (-0.020 if side == "left" else 0.020), y, lab,
            ha="right" if side == "left" else "left", va="center",
            fontsize=fs, zorder=6)

def curve(tab, cfg, r, hs):
    out = []
    for h in hs:
        g, _ = tsf.pair_conductance(float(h), r, r, 1, 1, tab, cfg)
        out.append(np.nan if g is None else g)
    return np.array(out)
