# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig4.py
#   sha256(src) : ca859e85580c9cbe
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 4. Computed on the fly; no result file is read.
from _branchfig import plt, np, CS, C_I, curve, save

C_REF = "#48566a"
import _data


def panel_a(ax):
    hs2 = np.linspace(-1.5, 2.5, 600)
    for tr, c, ls, lab in (("untreated", "#2c6fbb", "-", "untreated"),
                           ("treated", "#4a5568", "--", "treated")):
        d2 = _data.build(n=40, seed=1, treat=tr)
        g = curve(d2["table"], d2["cfg"], 10.0, hs2)
        m = ~np.isnan(g)
        ax.semilogy(hs2[m], g[m], ls, color=c, lw=1.6, label=lab)
        up = int((np.diff(g[m]) > 0).sum())
        print(f"(a) {tr:10s}: increasing steps {up}/{len(g[m])-1}")
    ax.axvline(0, color=C_REF, lw=1.2, ls=(0, (2.5, 2)), zorder=1)
    ax.text(0.20, 0.025, "contact branch", transform=ax.transAxes,
            ha="center", fontsize=7.9, color=C_REF)
    ax.text(0.72, 0.025, "near-contact branch", transform=ax.transAxes,
            ha="center", fontsize=7.9, color=C_REF)
    ax.set_xlabel(r"$h$  [$\mu$m]   ($h<0$: overlap)", fontsize=8.6)
    ax.set_ylabel(r"$g_e$  [W/K]", fontsize=8.6)
    ax.legend(frameon=False, fontsize=7.8, loc="upper right")
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_title("(a)", loc="left", fontsize=9.2)


def panel_b(ax):
    base = _data.build(n=40, seed=1)
    tab, cfg = base["table"], base["cfg"]
    hs = np.logspace(-3, np.log10(3.0), 260)
    ds = [2, 5, 10, 20, 50, 100]
    gh = {}
    for c, dd in zip(CS, ds):
        g = curve(tab, cfg, dd / 2, hs)
        ax.loglog(hs, g, "-", color=c, lw=1.4, label=f"{dd}")
        gh[dd] = g
    h0 = 0.1
    i0 = int(np.argmin(abs(hs - h0)))
    col = np.array([gh[dd][i0] for dd in ds])
    dec = np.log10(np.nanmax(col) / np.nanmin(col))
    ax.axvline(h0, color=C_REF, lw=1.0, ls=(0, (4, 2.5)), zorder=1)
    ax.annotate("", xy=(h0, np.nanmin(col)), xytext=(h0, np.nanmax(col)),
                arrowprops=dict(arrowstyle="<->", color=C_I, lw=1.2))
    ax.text(h0 * 1.4, 5.8e-5, f"{dec:.2f} decades\nat the same $h$",
            fontsize=7.8, color=C_I, ha="left", va="center", linespacing=1.4)
    ax.set_xlabel(r"surface separation  $h$  [$\mu$m]", fontsize=8.6)
    ax.set_ylabel(r"$g_e$  [W/K]", fontsize=8.6)
    ax.set_ylim(bottom=1.6e-7)
    leg = ax.legend(title=r"$d$ [$\mu$m]", frameon=False, fontsize=7.6,
                    title_fontsize=7.6, loc="lower left", ncol=3,
                    columnspacing=0.8, handlelength=1.2, labelspacing=0.25)
    leg._legend_box.align = "left"
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_title("(b)", loc="left", fontsize=9.2)
    print(f"(b) at h={h0} um, sweeping d from 2 to 100 um spans {dec:.2f} decades")


def main():
    fig, axs = plt.subplots(1, 2, figsize=(6.12, 2.85))
    panel_a(axs[0])
    panel_b(axs[1])
    fig.subplots_adjust(left=0.085, right=0.985, top=0.92, bottom=0.15,
                        wspace=0.30)
    save(fig, "fig4_conductance_axis")


if __name__ == "__main__":
    main()
