# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig11.py
#   sha256(src) : 18322413cb0687e8
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 11. Reads data/results_s2_scaling.csv. The critical level (b) is not
# timed: the measurement is not in the data.
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
CSV = _paths.poc_csv("S2_scaling", "results_s2_scaling.csv")
COL = {"a": "#4a5568", "b": "#2c6fbb", "c": "#2f855a", "dd": "#c0392b", "ds": "#e07b39"}
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})


def pick(df, *names):
    for n in names:
        if n in df.columns:
            return n
    raise KeyError(names)


def main():
    d = pd.read_csv(CSV)
    assert "engine_version" in d.columns, "refusing a CSV with no engine stamp"
    print("N:", sorted(d.N.unique()), " status:", sorted(d.d_dense_status.unique()))
    cN = pick(d, "N")
    cE = pick(d, "n_edge")
    ca = pick(d, "t_a_build_s")
    cb = pick(d, "t_b_barcode_s")
    cc = pick(d, "t_c_pshare_s")
    cds = pick(d, "t_d_sparse_topk_s")
    cdd = pick(d, "t_d_dense_s")
    d.loc[d["d_dense_status"] != "ok", cdd] = float("nan")
    g = d.groupby(cN).median(numeric_only=True).reset_index()
    E = d.groupby(cN)[cE].median()

    fig, axs = plt.subplots(1, 2, figsize=(6.25, 2.62))

    ax = axs[0]
    N = g[cN].values
    tot = g['t_abc_s'].values
    for c, lab, col, mk in ((ca, "preprocessing", COL["a"], "o"),
                            (cb, "(a) barcode", COL["b"], "s"),
                            (cc, "(c) $P^{\\rm share}$", COL["c"], "^"),
                            (cds, "(d) finite change, sparse", COL["ds"], "d")):
        ax.loglog(N, g[c], mk + "-", color=col, lw=1.3, ms=5, label=lab)
    dm = (g[cdd].notna() & (g[cdd] > 0)).values
    ax.loglog(N[dm], g[cdd].values[dm], "v--", color=COL["dd"], lw=1.3, ms=6,
              label="(d) finite change, dense")
    nwall = N[~dm].min() if (~dm).any() else None
    if nwall is not None:
        ax.axvline(nwall * 0.72, color=COL["dd"], lw=1.0, ls=":")
        ax.text(nwall * 0.80, 2.6e-3,
                "dense: N/A beyond\n$\\tt{max\\_n}$", ha="left", va="bottom",
                fontsize=7, color=COL["dd"], linespacing=1.4)
    ax.loglog(N, tot, "-", color="#94a3b8", lw=2.4, alpha=0.55, zorder=1,
              label="pre + (a) + (c)")
    ref = tot[0] * (N / N[0])
    ax.loglog(N, ref, ":", color="#cbd5e0", lw=1.0)
    ax.text(N[-2], ref[-2] * 0.45, "slope 1", fontsize=6.8, color="#94a3b8")
    ax.set_xlabel("$N$   particles", fontsize=8)
    ax.set_ylabel("time [s]  (median of 3)", fontsize=8)
    ax.set_title("(a)  time per procedure", loc="left", fontsize=8.6)
    ax.grid(alpha=0.25, lw=0.5)

    ax = axs[1]
    frac = (100 * g[ca].values / tot)
    ax.plot(N, frac, "o-", color=COL["a"], lw=1.5, ms=5,
            label="preprocessing share  (left)")
    ax.set_xscale("log")
    ax.set_ylim(0, 100)
    ax.set_xlabel("$N$   particles", fontsize=8)
    ax.set_ylabel("preprocessing share of pre + (a) + (c)  [%]", color=COL["a"], fontsize=8)
    ax.tick_params(axis="y", colors=COL["a"])
    ax.axhline(50, color="#cbd5e0", lw=0.8, ls="--")
    biggest_pre = g[ca].values > np.maximum(g[cb].values, g[cc].values)
    assert biggest_pre[0], ("preprocessing is not the largest stage even at the "
                            "smallest N; the premise of this figure no longer holds")
    n_last = N[:1 + (len(N) - 1 if biggest_pre.all()
                     else int(np.argmin(biggest_pre)) - 1)][-1]

    ax2 = ax.twinx()
    norm = (g["b_over_ElogE"] / g["b_over_ElogE"].iloc[0]).values
    raw_t = g[cb].values
    ax2.plot(N, norm, "s--", color=COL["b"], lw=1.4, ms=5,
             label=r"$t_b/(|E|\log|E|)$  (right)")
    ax2.set_ylabel(r"$t_b\,/\,(|E|\log|E|)$  (rel. $N$=200)",
                   color=COL["b"], fontsize=7.2)
    ax2.tick_params(axis="y", colors=COL["b"])
    ax2.set_ylim(0, 2.2)
    ax2.axhline(1.0, color=COL["b"], lw=0.7, ls=":")
    ax.set_title("(b)  share of preprocessing in the total time", loc="left",
                 fontsize=8.6)
    for a in (ax, ax2):
        a.spines["top"].set_visible(True)
    ax.grid(alpha=0.25, lw=0.5)

    fig.tight_layout(pad=0.9, w_pad=2.4, rect=[0, 0.185, 1, 1])
    twin = {1: ax2}
    for k, (a, nc) in enumerate(zip(axs, (2, 1))):
        h, l = a.get_legend_handles_labels()
        if k in twin:
            h2, l2 = twin[k].get_legend_handles_labels()
            h, l = h + h2, l + l2
        pos = a.get_position()
        fig.legend(h, l, loc="upper center", ncol=nc, frameon=False, fontsize=7.6,
                   columnspacing=1.4, handlelength=2.0, labelspacing=0.35,
                   bbox_to_anchor=((pos.x0 + pos.x1) / 2, 0.178))
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig11_scaling.{ext}", dpi=220, bbox_inches="tight")
    print(f"(a) share {frac.min():.1f}-{frac.max():.1f}%   "
          f"t_b/(ElogE) end {norm[-1]:.2f}x spread {norm.max()/norm.min():.1f}x   "
          f"total at N={int(N[-1])}: {tot[-1]:.2f}s")
    print("wrote", OUT / "fig11_scaling.pdf")


if __name__ == "__main__":
    main()
