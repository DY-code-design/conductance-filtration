# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig14.py
#   sha256(src) : f9764cad6cb92290
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source. In addition, 1 top-level definitions that the figures and
# verification scripts of the paper do not use were omitted (they are listed
# by name in rb_manifest); the remaining definitions are unchanged.
# ---------------------------------------------------------------------------
# Fig. 14. Reads data/results_d1_demo.csv and data/results_d1_edges.csv.
# The four material conditions sharing a diameter and a seed use the same
# structure. The treatment ratio is a median over seeds, not a ratio of medians.
import pathlib

import numpy as np
import pandas as pd
import matplotlib

import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()
CSV = _paths.poc_csv("D1_public_demo", "results_d1_demo.csv")
EDG = _paths.poc_csv("D1_public_demo", "results_d1_edges.csv")

COL = {"AlN": "#2c6fbb", "SiO2": "#c0392b",
       "contact": "#2c6fbb", "gap": "#e8b923", "top": "#c0392b", "bg": "#c8d0da",
       "interface": "#2c6fbb", "matrix": "#e8b923", "filler": "#7d8794"}
plt.rcParams.update({"font.size": 7.6, "axes.linewidth": 0.55})

ROWS = [(dd, f, t) for f in ("AlN", "SiO2") for dd in (5.0, 25.0)
        for t in ("untreated", "treated")]
COMBOS = [(dd, f) for f in ("AlN", "SiO2") for dd in (5.0, 25.0)]


def panel_a(ax, d, med):
    pair = d.pivot_table(index=["d_um", "filler", "seed"], columns="treat",
                         values="k_eff")
    pair["ratio"] = pair["treated"] / pair["untreated"]
    rat = pair.groupby(["d_um", "filler"])["ratio"].median()
    for i, (dd, f) in enumerate(COMBOS):
        ku = med.loc[(dd, f, "untreated"), "k_eff"]
        kt = med.loc[(dd, f, "treated"), "k_eff"]
        for s in sorted(d.seed.unique()):
            ax.plot([i, i], [pair.loc[(dd, f, s), "untreated"],
                             pair.loc[(dd, f, s), "treated"]],
                    "-", color=COL[f], lw=0.6, alpha=0.30, zorder=1)
        ax.plot([i, i], [ku, kt], "-", color=COL[f], lw=1.0, alpha=0.55)
        ax.plot(i, ku, "o", color="white", mec=COL[f], mew=1.1, ms=5, zorder=3)
        ax.plot(i, kt, "o", color=COL[f], ms=5, zorder=3)
        ax.annotate(f"$\\times${rat.loc[(dd, f)]:.3f}", (i, max(ku, kt)),
                    textcoords="offset points", xytext=(0, 8),
                    ha="center", fontsize=6.8, color=COL[f])
    ax.set_yscale("log"); ax.set_ylim(1.0, 70.0)
    ax.yaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 20, 50]))
    fmt = ScalarFormatter(); fmt.set_scientific(False)
    ax.yaxis.set_major_formatter(fmt)
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks(range(len(COMBOS)))
    ax.set_xticklabels([f"{f}\n{dd:.0f} " + r"$\mu$m" for dd, f in COMBOS],
                       fontsize=6.8, linespacing=1.4)
    ax.set_xlim(-0.6, len(COMBOS) - 0.4)
    ax.set_ylabel(r"$k_{\rm eff}$  [W/(m$\cdot$K)]", fontsize=8)
    ax.set_title(r"(a)  $k_{\rm eff}$ and the treatment ratio", loc="left",
                 fontsize=8.4)
    def pair(fill):
        return tuple(Line2D([], [], marker="o", ls="none", ms=5, mew=1.1,
                            mec=COL[f], mfc=(COL[f] if fill else "white"))
                     for f in ("AlN", "SiO2"))

    ax.legend([pair(False), pair(True)], ["untreated", "treated"],
              handler_map={tuple: HandlerTuple(ndivide=None, pad=0.55)},
              frameon=False, fontsize=6.8, loc="center left",
              handlelength=2.0, handletextpad=0.6)
    ax.tick_params(labelsize=7.0)


def panel_b(ax, med):
    y = np.arange(len(ROWS))[::-1]
    sc = np.array([med.loc[r, "share_contact"] for r in ROWS])
    sg = np.array([med.loc[r, "share_gap"] for r in ROWS])
    ax.barh(y, sc, height=0.62, color=COL["contact"], label="contact")
    ax.barh(y, sg, height=0.62, left=sc, color=COL["gap"],
            label="near-contact")
    for yy, a in zip(y, sc):
        ax.text(a - 0.02, yy, f"{a:.3f}", ha="right", va="center",
                fontsize=6.6, color="white")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{f} {dd:.0f}" + r"$\mu$m " + ("tr" if t == "treated" else "un")
                        for dd, f, t in ROWS], fontsize=6.6)
    ax.set_xlim(0, 1.06); ax.set_ylim(-0.7, 9.1)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("share of the total dissipation", fontsize=8)
    ax.legend(frameon=False, fontsize=7.4, loc="upper center", ncol=2,
              handlelength=1.4, handletextpad=0.6, columnspacing=2.0)
    ax.set_title("(b)  dissipation share by branch type", loc="left", fontsize=8.4)
    ax.tick_params(labelsize=7.0)


def panel_c(ax, e):
    P = e["P_share"].values
    k = max(1, int(round(0.01 * len(P))))
    cut = np.sort(P)[-k]
    hot = P >= cut
    ax.scatter(e["x_mid"][~hot], e["y_mid"][~hot], s=1.1, c="black",
               lw=0, alpha=0.45, rasterized=True)
    ax.scatter(e["x_mid"][hot], e["y_mid"][hot], s=13, c=COL["top"],
               lw=0, alpha=0.9, rasterized=True)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$  [$\mu$m]", fontsize=8)
    ax.set_ylabel(r"$y$  [$\mu$m]", fontsize=8)
    ax.set_title("(c)  spatial distribution of the dissipation", loc="left",
                 fontsize=8.4)
    print(f"(c) top 1% = {k} of {len(P):,} edges, {100*P[hot].sum():.0f}% of the dissipation")
    ax.tick_params(labelsize=7.0)


def panel_d(ax, med):
    x = np.arange(len(ROWS))
    th = np.array([med.loc[r, "theta_star"] for r in ROWS])
    ax.bar(x, th, width=0.62, color=[COL[f] for _, f, _ in ROWS], alpha=0.85)
    ax.set_ylabel(r"$\theta^*$", fontsize=8)
    ax.set_ylim(0, 6.2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{f} {dd:.0f} " + ("tr" if t == "treated" else "un")
                        for dd, f, t in ROWS], fontsize=6.4, rotation=90)
    ax.set_xlim(-0.6, len(ROWS) - 0.4)
    ax.tick_params(labelsize=7.0)
    ax.set_title(r"(d)  critical conductance level $\theta^*$  ($d$ in $\mu$m)",
                 loc="left", fontsize=8.4)


def panel_e(ax, med):
    x = np.arange(len(ROWS))
    marks = [("gain_interface", "interface", r"interface  (resistances $\to 0$)"),
             ("gain_matrix", "matrix", r"matrix  ($k_{\rm m}\to\infty$)"),
             ("gain_filler", "filler", r"filler  ($\kappa_{\rm p}\to\infty$)")]
    w = 0.26
    for k, (col, lab, leg) in enumerate(marks):
        v = np.array([med.loc[r, col] for r in ROWS])
        ax.bar(x + (k - 1) * w, v, width=w * 0.92, color=COL[lab], label=leg,
               alpha=0.9)
    ax.axhline(1.0, color="#94a3b8", lw=0.9, ls="--",
               label=r"$\Gamma_m=1$  (no change)")
    ax.set_yscale("log"); ax.set_ylim(0.8, 1500)
    ax.yaxis.set_major_locator(FixedLocator([1, 10, 100, 1000]))
    fmt = ScalarFormatter(); fmt.set_scientific(False)
    ax.yaxis.set_major_formatter(fmt)
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks(x)
    ax.set_xticklabels([f"{f} {dd:.0f}" + r"$\mu$m" + "\n" + t
                        for dd, f, t in ROWS], fontsize=6.4, linespacing=1.35)
    ax.set_xlim(-0.6, len(ROWS) - 0.4)
    ax.set_ylabel(r"$\Gamma_m$", fontsize=8, labelpad=1.0)
    ax.tick_params(labelsize=7.0)
    ax.set_title(r"(e)  idealisation ratio $\Gamma_m$ per mechanism", loc="left",
                 fontsize=8.4)
    ax.legend(frameon=False, fontsize=6.6, ncol=4, loc="upper left",
              handlelength=1.4, columnspacing=1.1)


def main():
    d = pd.read_csv(CSV)
    e = pd.read_csv(EDG)
    for df in (d, e):
        assert "engine_version" in df.columns, "refusing a CSV with no engine stamp"
    med = d.groupby(["d_um", "filler", "treat"]).median(numeric_only=True)

    fig = plt.figure(figsize=(5.8, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.72, 1.02, 0.66],
                          wspace=0.42, hspace=0.42)
    panel_a(fig.add_subplot(gs[0, 0]), d, med)
    panel_b(fig.add_subplot(gs[0, 1]), med)
    panel_c(fig.add_subplot(gs[1, 0]), e)
    panel_d(fig.add_subplot(gs[1, 1]), med)
    panel_e(fig.add_subplot(gs[2, :]), med)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.965, bottom=0.075)
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig14_demo.{ext}", dpi=220, bbox_inches="tight")
    print("wrote", OUT / "fig14_demo.pdf")


if __name__ == "__main__":
    main()
