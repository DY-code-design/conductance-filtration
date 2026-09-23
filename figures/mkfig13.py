# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig13.py
#   sha256(src) : 7de8bdec5caab40f
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 13. Reads data/results_ensemble_theta.csv and data/results_y7_holdout.csv.
# The fit and the band are computed from the data and asserted against the
# rounded values printed in the paper, so the figure fails if the data move.
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
CSV = _paths.poc_csv("E_ensemble", "results_ensemble_theta.csv")
HOLDOUT = _paths.poc_csv("E_ensemble", "results_y7_holdout.csv")
A_PUB, B_PUB, RES_PUB = -1.31, 1.96, 0.250
C_ALL, C_PHI06, C_FIT = "#94a3b8", "#c0392b", "#2c6fbb"
LAB_ALL = "the other 42 structures"
LAB_06 = r"$\phi=0.600$ (27 structures)"
plt.rcParams.update({"font.size": 7.2, "axes.linewidth": 0.6})


def main():
    d = pd.read_csv(CSV)
    assert set(d.engine_variant) == {"holm"}, d.engine_variant.unique()
    assert set(d.baseline_fp) == {"014d59a0d7b143ba"}, d.baseline_fp.unique()
    assert len(d) == 69, len(d)
    lg = np.log10(d.G_eff)
    is06 = np.isclose(d.phi, 0.6, atol=0.005)
    assert int(is06.sum()) == 27, int(is06.sum())

    A, B = (float(x) for x in np.polyfit(d.theta_star, lg, 1))
    h = pd.read_csv(HOLDOUT)
    assert set(h.baseline_fp) == {"014d59a0d7b143ba"}, h.baseline_fp.unique()
    assert int(h.n.sum()) == 69, int(h.n.sum())
    RES = float(np.sqrt((h.n * h.res_rms ** 2).sum() / h.n.sum()))
    for got, pub, name in ((A, A_PUB, "a"), (B, B_PUB, "b"), (RES, RES_PUB, "RES")):
        assert round(got, 3 if name == "RES" else 2) == pub,\
            f"{name}: {got!r} from the data but {pub!r} in the paper (one of them is stale)"

    fig, axs = plt.subplots(1, 2, figsize=(6.25, 2.7), sharey=True)

    ax = axs[0]
    ax.plot(d.phi[~is06], lg[~is06], "o", ms=4, color=C_ALL, mec="none",
            alpha=0.85, label=LAB_ALL)
    ax.plot(d.phi[is06], lg[is06], "o", ms=4.5, color=C_PHI06, mec="none",
            label=LAB_06)
    ax.plot([0.6, 0.6], [lg[is06].min(), lg[is06].max()], "-", color=C_PHI06, lw=0.8)
    ax.legend(frameon=False, fontsize=6.8, loc="upper left", handletextpad=0.4,
              borderpad=0.2)
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$\log_{10}G_{\rm eff}$")
    ax.set_title(r"(a)  $\log_{10}G_{\rm eff}$ vs. $\phi$", loc="left", fontsize=8)
    ax.grid(alpha=0.25, lw=0.5)

    ax = axs[1]
    xs = np.linspace(d.theta_star.min() - 0.05, d.theta_star.max() + 0.05, 10)
    ax.fill_between(xs, A * xs + B - RES, A * xs + B + RES, color=C_FIT,
                    alpha=0.12, lw=0)
    ax.plot(xs, A * xs + B, "-", color=C_FIT, lw=1.2)
    ax.plot(d.theta_star[~is06], lg[~is06], "o", ms=4, color=C_ALL, mec="none",
            alpha=0.85, label=LAB_ALL)
    ax.plot(d.theta_star[is06], lg[is06], "o", ms=4.5, color=C_PHI06, mec="none",
            label=LAB_06)
    ax.legend(frameon=False, fontsize=6.8, loc="upper right", handletextpad=0.4,
              borderpad=0.2)
    ax.annotate(f"$\\log_{{10}}G_{{\\rm eff}}={A:.2f}\\,\\theta^*+{B:.2f}$\n"
                f"band: $\\pm${RES:.3f} decades", xy=(2.95, -4.1),
                fontsize=6.8, color=C_FIT)
    ax.set_xlabel(r"$\theta^*$")
    ax.set_title(r"(b)  $\log_{10}G_{\rm eff}$ vs. $\theta^*$", loc="left",
                 fontsize=8)
    ax.grid(alpha=0.25, lw=0.5)

    for a_ in axs:
        a_.spines["top"].set_visible(False)
        a_.spines["right"].set_visible(False)
    fig.tight_layout(pad=0.9, w_pad=1.4)
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig13_theta_estimate.{ext}", dpi=220,
                    bbox_inches="tight", pad_inches=0.05, facecolor="white")
    print("wrote", OUT / "fig13_theta_estimate.pdf")


if __name__ == "__main__":
    main()
