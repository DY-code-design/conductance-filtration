# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig6.py
#   sha256(src) : a8ee2b0fba0c1102
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 6. Computed on the fly; no result file is read.
# The axes are ranks; tau is the Kendall rank correlation between the order by
# centre distance D and the order by theta = -log10 g.
import pathlib

import numpy as np
import matplotlib

import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kendalltau
import _data

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()
C_OK, C_NG, C_T = "#2c6fbb", "#c0392b", "#4a5568"
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})


def ranks(v):
    v = np.asarray(v, float)
    o = np.argsort(v, kind="stable")
    r = np.empty(len(v)); r[o] = np.arange(len(v))
    for x in np.unique(v):
        m = v == x
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def panel_rank(ax, rows, tag, title, colour, note):
    D = np.array([r["D"] for r in rows])
    th = np.array([-np.log10(r["g"]) for r in rows])
    t = kendalltau(D, th).statistic
    rd, rt = ranks(D) / len(D), ranks(th) / len(th)
    ax.plot([0, 1], [0, 1], "-", color="#94a3b8", lw=1.0, zorder=1)
    ax.scatter(rd, rt, s=9, color=colour, alpha=0.55, lw=0, zorder=3)
    ax.set_xlabel(r"rank by $D$   (centre-to-centre distance)", fontsize=8)
    ax.set_ylabel(r"rank by $\theta=-\log_{10}g$", fontsize=8)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.set_title(f"{tag}  {title}", loc="left", fontsize=8.6)
    ax.text(0.04, 0.96, rf"$\tau(D,\theta)={t:+.3f}$",
            transform=ax.transAxes, va="top", fontsize=8.0, color=C_T)
    ax.grid(alpha=0.22, lw=0.5)
    return t


def main():
    r_mono = _data.edge_table(_data.build(family="phi_060", n=400, seed=0))
    r_bi = _data.edge_table(_data.build(family="bimodal_5to1", n=400, seed=0))

    fig, axs = plt.subplots(1, 2, figsize=(6.25, 3.05))
    tm = panel_rank(axs[0], r_mono, "(a)", "monodisperse, one material", C_OK,
                    "every edge on the diagonal:\n"
                    r"the two orders coincide")
    tb = panel_rank(axs[1], r_bi, "(b)", "bimodal 5:1", C_NG,
                    "same material; only the size\n"
                    r"distribution is changed")
    fig.tight_layout(pad=0.9, w_pad=2.0)
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig6_size_breaks_order.{ext}", dpi=220,
                    bbox_inches="tight")

    print(f"tau(D,theta): monodisperse {tm:+.4f} / bimodal {tb:+.4f}")
    print("wrote", OUT / "fig6_size_breaks_order.pdf")


if __name__ == "__main__":
    main()
