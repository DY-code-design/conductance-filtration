# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig7.py
#   sha256(src) : d3e835b37b52aeee
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 7. Reads phi = 0.600, seed 0 from data/results_t1.csv
# (D* = tstar_G2, theta* = tstar_C).
import pathlib

import numpy as np
import pandas as pd
import matplotlib

import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()

def load_t1():
    d = pd.read_csv(_paths.poc_csv("T1_filtration_compare", "results_t1.csv"))
    assert set(d.engine_variant) == {"holm"}, set(d.engine_variant)
    assert set(d.baseline_fp) == {"014d59a0d7b143ba"}, set(d.baseline_fp)
    rows = []
    for filler, matrix, treat in (("AlN", "epoxy", "untreated"),
                                  ("SiO2", "epoxy", "untreated"),
                                  ("AlN", "high_k_resin", "treated")):
        r = d[(d.family == "phi_060") & (d.seed == 0) & (d.filler == filler)
              & (d.matrix == matrix) & (d.treat == treat)]
        assert len(r) == 1, (filler, matrix, treat, len(r))
        rows.append(r.iloc[0])
    Dstar = np.array([r.tstar_G2 for r in rows])
    tstar = np.array([r.tstar_C for r in rows])
    assert len({round(v, 6) for v in Dstar}) == 1, Dstar
    return Dstar, tstar
C_D, C_G = "#3b6ea5", "#c0392b"
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})


def panel_d(ax_top, ax_bot):
    mats = ["AlN /\nepoxy", "SiO$_2$ /\nepoxy", "AlN / high-$\\kappa$\nresin, treated"]
    Dstar, tstar = load_t1()
    x = np.arange(3)

    ax_top.plot(x, Dstar, "o-", color=C_D, lw=1.5, ms=6)
    ax_top.set_ylim(19.0, 19.5)
    ax_top.set_yticks([19.0, 19.25, 19.5])
    ax_top.tick_params(axis="y", labelsize=7.8)
    ax_top.set_ylabel(r"$D^*$ [$\mu$m]", color=C_D, fontsize=8.7, labelpad=1)
    ax_top.tick_params(axis="y", colors=C_D)
    ax_top.set_xticks(x); ax_top.set_xticklabels([])
    ax_top.set_xlim(-0.45, 2.45)
    ax_top.text(1.0, 19.32, "identical to 6 decimals", ha="center", va="bottom",
                fontsize=8.3, color=C_D, style="italic")

    ax_bot.plot(x, tstar, "s-", color=C_G, lw=1.5, ms=6)
    ax_bot.set_ylim(2.95, 5.35)
    ax_bot.set_ylabel(r"$\theta^*=-\log_{10}g^*$", color=C_G, fontsize=8.7, labelpad=1)
    ax_bot.tick_params(axis="y", colors=C_G, labelsize=7.8)
    ax_bot.set_xticks(x); ax_bot.set_xticklabels(mats, fontsize=7.9)
    ax_bot.set_xlim(-0.45, 2.45)
    for xi, v in zip(x, tstar):
        off = {0: (-8, 9), 1: (13, 0), 2: (-14, 2)}[int(xi)]
        ha = {0: "center", 1: "left", 2: "right"}[int(xi)]
        ax_bot.annotate(f"{v:.3f}", (xi, v), textcoords="offset points",
                        xytext=off, ha=ha, va="center",
                        fontsize=8, color=C_G)
    ax_bot.annotate("", xy=(2.28, float(tstar.min())), xytext=(2.28, float(tstar.max())),
                    arrowprops=dict(arrowstyle="<->", color=C_G, lw=1.0))
    ax_bot.text(2.22, 4.02, f"{tstar.max() - tstar.min():.3f}\ndecades", fontsize=8, color=C_G,
                va="center", ha="right", linespacing=1.3)
    for a in (ax_top, ax_bot):
        a.spines[["top", "right"]].set_visible(False)


def main():
    fig = plt.figure(figsize=(4.1, 4.0))
    gs = gridspec.GridSpec(2, 1, figure=fig, hspace=0.12, height_ratios=[1, 1.35])
    panel_d(fig.add_subplot(gs[0]), fig.add_subplot(gs[1]))
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig7_material_only.{ext}", dpi=220, bbox_inches="tight")
    print("wrote", OUT / "fig7_material_only.pdf")


if __name__ == "__main__":
    main()
