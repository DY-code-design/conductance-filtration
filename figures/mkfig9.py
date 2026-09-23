# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig9.py
#   sha256(src) : 43b7f191c9ca47fb
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 9. Reads data/results_f4_meta.csv and data/results_f4_theta_sweep.csv.
import pathlib
import numpy as np
import pandas as pd
import matplotlib
import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()
C_B, C_G, C_S = "#4a5568", "#2c6fbb", "#e07b39"
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})


def main():
    d = pd.read_csv(_paths.poc_csv("F4_theta_sweep", "results_f4_theta_sweep.csv")).sort_values("theta")
    m = pd.read_csv(_paths.poc_csv("F4_theta_sweep", "results_f4_meta.csv")).iloc[0]
    assert "engine_version" in d.columns and "engine_version" in m.index
    ts, gs = float(m.theta_star_percolation), float(m.g_star_percolation)
    assert abs(gs - 10 ** (-ts)) / gs < 1e-3, (gs, 10 ** (-ts))
    assert f"{gs:.2e}" == "4.58e-04", f"g* moved: {gs:.3e} (the paper, Sec. 3.6, has 4.58e-4)"
    print(f"{m.family} {m.filler}/{m.matrix}/{m.treat}  N={m.n_particle} E={m.n_edge} "
          f"phi={m.phi}  engine={m.engine_version}")
    print(f"theta* = {ts:.4f}  g* = {gs:.3e} W/K  G_eff(full) = {m.G_eff_full:.4e} W/K")

    fig, ax = plt.subplots(figsize=(5.8, 3.6))

    ax.plot(d.theta, d.beta0, "-", color=C_B, lw=1.6, label=r"$\beta_0(\theta)$")
    ax.set_xlabel(r"$\theta=-\log_{10}g$   (larger $\theta$ = weaker edges included)", fontsize=8.6)
    ax.set_ylabel(r"$\beta_0(\theta)$   number of components", color=C_B, fontsize=8.6)
    ax.tick_params(axis="y", colors=C_B)
    ax.set_ylim(0, d.beta0.max() * 1.06)

    ax2 = ax.twinx()
    pos = d[d.G_eff > 0]
    ax2.semilogy(pos.theta, pos.G_eff, "-", color=C_G, lw=1.8, label=r"$G_{\rm eff}(\theta)$")
    ax2.set_ylabel(r"$G_{\rm eff}(\theta)$   [W/K]", color=C_G, fontsize=8.6)
    ax2.tick_params(axis="y", colors=C_G, labelsize=8)
    ax2.set_ylim(pos.G_eff[pos.G_eff > 0].min() * 0.4, float(m.G_eff_full) * 3.2)

    HEAD = 1.34
    ax.set_ylim(0, ax.get_ylim()[1] * HEAD)
    b2, t2 = ax2.get_ylim()
    ax2.set_ylim(b2, t2 * (t2 / b2) ** (HEAD - 1.0))

    ax.axvline(ts, color=C_S, lw=1.2, ls="--", zorder=1)
    ytop = ax.get_ylim()[1]
    ax.text(ts + 0.05, ytop * 0.855,
            rf"$\theta^*={ts:.3f}$   $g^*={gs/1e-4:.2f}\times10^{{-4}}$ W/K",
            fontsize=8.6, color=C_S, ha="left", va="top")
    C_REG = "#3f4d60"
    x0 = ax.get_xlim()[0]
    ax.axvspan(x0, ts, color="#e7edf5", zorder=0)
    ax.set_xlim(left=x0)
    ax.text(ts - 0.02, ytop * 0.955, "no path\nto sink",
            ha="right", va="top", fontsize=8.8, color=C_REG, linespacing=1.35)
    ax.text(ts + 0.05, ytop * 0.955, "percolating",
            ha="left", va="top", fontsize=8.8, color=C_REG)


    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8.0, loc="upper right",
              handlelength=1.8)
    ax.set_title(f"{m.filler}/{m.matrix}, $\\phi$={float(m.phi):.2f}, N={int(m.n_particle)}, "
                 f"$|E|$={int(m.n_edge)}", loc="left", fontsize=8.5, color="#4a5568")
    for a in (ax, ax2):
        a.spines["top"].set_visible(False)
    ax.grid(alpha=0.22, lw=0.5)
    fig.tight_layout(pad=0.8)
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig9_theta_sweep.{ext}", dpi=220, bbox_inches="tight")
    print("wrote", OUT / "fig9_theta_sweep.pdf")


if __name__ == "__main__":
    main()
