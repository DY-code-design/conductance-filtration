# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig12.py
#   sha256(src) : ba8ccd1529053958
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 12. Reads data/results_poc5_grid.csv. The population follows the
# emt_scope column of the data, not a filter written here.
import pathlib
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.transforms import blended_transform_factory
import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()
CSV = _paths.poc_csv("PoC5_lattice_emt", "results_poc5_grid.csv")
COL = {"SC": "#c0392b", "BCC": "#2c6fbb", "FCC": "#2f855a"}
Z = {"SC": 6, "BCC": 8, "FCC": 12}
plt.rcParams.update({"font.size": 7.8, "axes.linewidth": 0.6})


C_REF = "#48566a"


def emt_label(ax):
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(0.992, 1.0, "EMT", transform=tr, ha="right", va="center",
            fontsize=7.6, color=C_REF, zorder=6,
            bbox=dict(fc="white", ec="none", pad=1.2))


def main():
    d = pd.read_csv(CSV)
    assert set(d.engine_variant) == {"holm"}, d.engine_variant.unique()
    assert set(d.baseline_fp) == {"014d59a0d7b143ba"}, d.baseline_fp.unique()
    assert "emt_scope" in d.columns, "no emt_scope column: the data predates the population definition"
    s = d[d.emt_scope].copy()
    assert len(s) == 3078, len(s)
    s["ratio"] = s.k_net / s.k_nan

    fig, axs = plt.subplots(1, 2, figsize=(6.25, 2.6))

    ax = axs[0]
    for lat in ("SC", "BCC", "FCC"):
        gp = s[s.lattice == lat]
        ax.semilogy(gp.phi_rel, gp.ratio, "o", ms=2.6, alpha=0.25,
                    color=COL[lat], mec="none")
        m = gp.groupby("phi_rel")["ratio"].median()
        ax.semilogy(m.index, m.values, "o-", ms=3.6, lw=1.1, color=COL[lat],
                    label=f"{lat} ($z$={Z[lat]})")
    ax.axhline(1, color=C_REF, lw=1.0, ls=(0, (4, 2.5)), zorder=1)
    emt_label(ax)
    ax.set_xlabel(r"$\phi/\phi_{\max}$", fontsize=8)
    ax.set_ylabel(r"$k_{\rm net}\,/\,k_{\rm EMT}$", fontsize=8)
    ax.set_title("(a)  deviation vs. packing, by lattice", loc="left", fontsize=8.6)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.grid(alpha=0.25, lw=0.5)

    ax = axs[1]
    for lat in ("SC", "BCC", "FCC"):
        gp = s[s.lattice == lat]
        m = gp.groupby("kappa")["ratio"].median()
        q1 = gp.groupby("kappa")["ratio"].quantile(0.25)
        q3 = gp.groupby("kappa")["ratio"].quantile(0.75)
        ax.fill_between(m.index, q1.values, q3.values, color=COL[lat],
                        alpha=0.12, lw=0)
        ax.semilogx(m.index, m.values, "o-", ms=3.6, lw=1.1, color=COL[lat],
                    label=lat)
    ax.axhline(1, color=C_REF, lw=1.0, ls=(0, (4, 2.5)), zorder=1)
    emt_label(ax)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\kappa_f/k_m$", fontsize=8)
    ax.set_ylabel(r"$k_{\rm net}\,/\,k_{\rm EMT}$   (median, IQR band)", fontsize=8)
    ax.set_title("(b)  ...and vs. the material contrast", loc="left", fontsize=8.6)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.grid(alpha=0.25, lw=0.5)

    for a in axs:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    fig.tight_layout(pad=0.9, w_pad=2.0)
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig12_emt_deviation.{ext}", dpi=220,
                    bbox_inches="tight", pad_inches=0.05, facecolor="white")
    print("wrote", OUT / "fig12_emt_deviation.pdf", f"({len(s)} conditions)")


if __name__ == "__main__":
    main()
