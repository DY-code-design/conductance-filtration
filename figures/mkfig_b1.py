# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig_b1.py
#   sha256(src) : 5171efb3f3b070dd
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. B.1. Concept figure; no data is read.
import math
import pathlib
import matplotlib
import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()
C_S, C_V = "#2c6fbb", "#4a5568"
C_PATH, C_ELD = "#c0392b", "#94a3b8"
NOTE_Y = -0.50
plt.rcParams.update({"font.size": 7.2, "axes.linewidth": 0.6})

P = {"a": (0.0, 1.55), "b": (0.95, 0.95), "c": (2.05, 0.95), "v": (3.0, 1.55),
     "s1": (0.6, 0.0), "s2": (2.4, 0.0)}
E = [("a", "b", 1.0), ("c", "v", 2.0), ("b", "c", 2.5),
     ("b", "s1", 3.0), ("c", "s2", 4.0)]


def draw_graph(ax, contract=False, path=None, bottleneck=None):
    pos = dict(P)
    pos["q"] = (1.5, 0.0)

    def pt(n):
        if contract and n in ("s1", "s2"):
            return pos["q"]
        return pos[n]

    for a, b, t in E:
        x1, y1 = pt(a)
        x2, y2 = pt(b)
        on_path = path and ((a, b) in path or (b, a) in path)
        ax.plot([x1, x2], [y1, y2], "-", lw=2.6 if on_path else 1.1,
                color=C_PATH if on_path else "#9aa5b1", zorder=1)
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        is_bn = bool(bottleneck and (a, b) == bottleneck)
        L = math.hypot(x2 - x1, y2 - y1)
        nx, ny = -(y2 - y1) / L, (x2 - x1) / L
        sgn = 1.0 if ((mx - 1.5) * nx + (my - 0.9) * ny) > 0 else -1.0
        off = 0.24 if is_bn else 0.155
        ax.text(mx + sgn * nx * off, my + sgn * ny * off, f"{t:g}",
                fontsize=6.6, color="#333", ha="center", va="center")
        if is_bn:
            ax.plot([mx], [my], "o", ms=13, mfc="none", mec=C_PATH, mew=1.4,
                    zorder=3)
    for n, (x, y) in P.items():
        if contract and n in ("s1", "s2"):
            continue
        if n.startswith("s"):
            ax.add_patch(Rectangle((x - 0.13, y - 0.13), 0.26, 0.26,
                                   fc=C_S, ec="none", zorder=2))
            ax.text(x, y - 0.34, n, fontsize=6.6, ha="center", color=C_S)
        else:
            ax.plot([x], [y], "o", ms=9, color=C_V, zorder=2)
            ax.text(x, y + 0.17, n, fontsize=7, ha="center")
    if contract:
        x, y = pos["q"]
        ax.add_patch(Rectangle((x - 0.15, y - 0.15), 0.30, 0.30,
                               fc=C_S, ec="none", zorder=2))
        ax.text(x, y - 0.36, "q (S contracted)", fontsize=6.4, ha="center",
                color=C_S)
    ax.set_xlim(-0.5, 3.7)
    ax.set_ylim(-0.85, 2.05)
    ax.axis("off")


def main():
    fig, axs = plt.subplots(1, 3, figsize=(6.25, 2.5))

    ax = axs[0]
    draw_graph(ax)
    ax.set_title("(a)  relative persistence,  S = {s1, s2}", loc="left",
                 fontsize=7.4)
    ax.text(1.6, NOTE_Y, "[1, 3]  dies by reaching S", fontsize=6.2,
            ha="center", va="top", color=C_S)
    ax.text(1.6, NOTE_Y - 0.167, "[2, 2.5]  dies by the elder rule",
            fontsize=6.2, ha="center", va="top", color="#6b7683")

    ax = axs[1]
    draw_graph(ax, contract=True)
    ax.set_title("(b)  quotient graph: S contracted to q", loc="left",
                 fontsize=7.4)
    ax.text(1.6, NOTE_Y, "same edges, hence the same\nunion–find as (a)",
            fontsize=6.2, ha="center", va="top", color="#555", linespacing=1.55)

    ax = axs[2]
    draw_graph(ax, path={("v", "c"), ("b", "c"), ("b", "s1")},
               bottleneck=("b", "s1"))
    ax.set_title("(c)  the widest path from v to S", loc="left", fontsize=7.4)
    ax.text(1.6, NOTE_Y,
            "$\\theta_v=\\max(2,\\ 2.5,\\ 3)=3$\ncircle: the bottleneck edge",
            fontsize=6.2, ha="center", va="top", color=C_PATH, linespacing=1.55)

    fig.tight_layout(pad=0.6, w_pad=0.7)
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig_b1_equivalence.{ext}", dpi=220,
                    bbox_inches="tight", pad_inches=0.03, facecolor="white")
    print("wrote", OUT / "fig_b1_equivalence.pdf")


if __name__ == "__main__":
    main()
