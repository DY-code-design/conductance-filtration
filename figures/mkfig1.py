# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig1.py
#   sha256(src) : 0fbfb8ee72988cf2
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 1. Concept figure; no data is read.
import pathlib
import matplotlib
import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()
C_STR, C_MOD, C_FIL, C_SOL, C_ACT = "#4a5568", "#a07b32", "#2c6fbb", "#c0392b", "#2f855a"
C_LINE = "#64748b"
plt.rcParams.update({"font.size": 8.5})


PAD = 0.010
BOXES = []
ENDS = []


def check_arrows():
    bad = []
    for ex, ey, what in ENDS:
        for bx, by, bw, bh in BOXES:
            if (abs(ex - bx) < bw / 2 + PAD) and (abs(ey - by) < bh / 2 + PAD):
                bad.append(f"{what}: endpoint ({ex}, {ey}) falls inside box ({bx}, {by})")
    assert not bad, "arrow head hidden:\n  " + "\n  ".join(bad)
    print(f"{len(ENDS)} arrow heads / {len(BOXES)} boxes -- none hidden inside a box")


def box(ax, x, y, w, h, title, body, colour, sec=None, fs=8.6):
    BOXES.append((x, y, w, h))
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.010,rounding_size=0.018",
                                fc="white", ec=colour, lw=1.3, zorder=3))
    pad_top = (h - (0.082 + (body.count("\n") + 1) * 0.049)) / 2
    ax.text(x, y + h / 2 - pad_top, title, ha="center", va="top", fontsize=fs,
            color=colour, zorder=4)
    ax.text(x, y + h / 2 - pad_top - 0.052, body, ha="center", va="top", fontsize=7.7,
            color="#2d3748", zorder=4, linespacing=1.5)
    if sec:
        ax.text(x + w / 2 - 0.012, y - h / 2 + 0.020, sec, ha="right", va="bottom",
                fontsize=7.4, color=colour, style="italic", zorder=4)


def arrow(ax, p0, p1, colour="#94a3b8", rad=0.0, lw=1.2, ls="-", shrinkA=2,
          what="arrow"):
    ENDS.append((p1[0], p1[1], what + " head"))
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=11,
                                 color=colour, lw=lw, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}", zorder=2,
                                 shrinkA=shrinkA, shrinkB=2))


def main():
    TOP, BOT = 0.75, 0.10
    fig, ax = plt.subplots(figsize=(6.12, 2.454))

    box(ax, 0.125, TOP - 0.1375, 0.22, 0.275, "packing",
        "$\\phi$, size distribution,\nfiller, surface treatment", C_STR)
    box(ax, 0.435, TOP - 0.1375, 0.29, 0.275, "branch-conductance model",
        "contact / near-contact / bond\n$\\Rightarrow$ one $g_e$ per edge", C_MOD, "§2", fs=8.2)
    box(ax, 0.790, TOP - 0.1375, 0.30, 0.275, "conductance filtration",
        "critical conductance level $\\theta^*$,\nsink-relative barcode\n"
        r"no solve — $O(|E|\log|E|)$", C_FIL, "§3", fs=8.4)

    box(ax, 0.790, BOT + 0.1375, 0.30, 0.275, "sparse solve",
        "$P^{\\rm share}=\\partial\\log G_{\\rm eff}/\\partial\\log g_e$,\n"
        "exact finite change\none solve + one factorisation", C_SOL, "§4", fs=8.4)
    box(ax, 0.435, BOT + 0.1375, 0.29, 0.275, "what to change",
        "which edge, which mechanism,\nhow much it returns", C_ACT, "§4, §6")
    box(ax, 0.125, BOT + 0.1375, 0.22, 0.275, "re-evaluate", "closes the loop", C_ACT)

    arrow(ax, (0.230, 0.6125), (0.277, 0.6125), C_LINE)
    arrow(ax, (0.560, 0.6125), (0.628, 0.6125), C_LINE)
    arrow(ax, (0.790, 0.500), (0.790, 0.388), C_LINE)
    arrow(ax, (0.680, 0.2375), (0.594, 0.2375), C_LINE)
    arrow(ax, (0.300, 0.2375), (0.250, 0.2375), C_LINE)
    arrow(ax, (0.125, 0.330), (0.125, 0.460), C_LINE, ls=(0, (4, 2)))

    check_arrows()
    ax.set_xlim(0, 1); ax.set_ylim(0.06, 0.79); ax.axis("off")
    fig.subplots_adjust(left=0.005, right=0.995, top=0.99, bottom=0.02)
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig1_design_loop.{ext}", dpi=220, bbox_inches="tight")
    print("wrote", OUT / "fig1_design_loop.pdf")


if __name__ == "__main__":
    main()
