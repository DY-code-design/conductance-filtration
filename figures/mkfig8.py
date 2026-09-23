# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig8.py
#   sha256(src) : 2829e720686bfd00
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 8. Concept figure on a schematic network of 12 nodes. All values are
# chosen for the figure, not measured.
import pathlib

import matplotlib

import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()
C_SNK, C_LIVE, C_DEAD, C_OFF = "#2c6fbb", "#c0392b", "#9aa5b1", "#dde3ea"
C_LBL = "#55606d"
C_TXT, C_ELDER = "#2d3748", "#e07b39"
plt.rcParams.update({"font.size": 8.4, "axes.linewidth": 0.6})

POS = {
    "a": (0.10, 0.86), "b": (0.42, 0.92), "c": (0.74, 0.84),
    "d": (0.06, 0.55), "e": (0.38, 0.60), "f": (0.72, 0.52),
    "g": (0.20, 0.28), "h": (0.55, 0.30), "i": (0.86, 0.26),
    "j": (0.12, 0.05), "k": (0.50, 0.04), "l": (0.88, 0.05),
}
SINK = ["j", "k", "l"]
EDGES = [("a", "b", 1.0), ("b", "c", 1.4), ("a", "d", 1.2), ("b", "e", 2.6),
         ("c", "f", 1.1), ("d", "e", 3.4), ("e", "f", 2.2), ("d", "g", 2.8),
         ("e", "h", 3.8), ("f", "i", 1.6), ("g", "h", 2.0), ("h", "i", 2.4),
         ("g", "j", 1.8), ("h", "k", 3.0), ("i", "l", 1.3)]
SNAP = [1.5, 2.5, 3.2]


def components(theta):
    parent = {n: n for n in POS}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    for u, v, t in sorted(EDGES, key=lambda e: e[2]):
        if t <= theta:
            parent[find(u)] = find(v)
    comp = {}
    for n in POS:
        comp.setdefault(find(n), []).append(n)
    return {r: (m, any(x in SINK for x in m)) for r, m in comp.items()}


def draw_net(ax, theta, title, show_elder=False):
    for u, v, t in EDGES:
        on = t <= theta
        ax.plot(*zip(POS[u], POS[v]), "-", lw=1.9 if on else 0.8,
                color=C_TXT if on else C_OFF, zorder=1)
    cmp_ = components(theta)
    for r, (members, has_sink) in cmp_.items():
        col = C_SNK if has_sink else C_LIVE
        for n in members:
            ax.plot(*POS[n], "o", ms=6.2, color=col, mec="white", mew=0.8, zorder=3)
    for n in SINK:
        ax.plot(*POS[n], "o", ms=14, mfc="none", mec=C_SNK, mew=1.5, zorder=2)
    if show_elder:
        mid = ((POS["b"][0] + POS["c"][0]) / 2, (POS["b"][1] + POS["c"][1]) / 2)
        ax.annotate("elder rule\ncloses one bar",
                    xy=mid, xytext=(0.34, 0.24), fontsize=9.6,
                    color=C_ELDER, linespacing=1.35, ha="center",
                    bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0),
                    arrowprops=dict(arrowstyle="->", color=C_ELDER, lw=0.9))
    ax.set_xlim(-0.06, 1.02); ax.set_ylim(-0.10, 1.02)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, loc="left", fontsize=11.0)


def main():
    fig = plt.figure(figsize=(5.55, 7.15))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.18, 1.18, 0.98],
                          hspace=0.25, wspace=0.02,
                          left=0.012, right=0.995, top=0.985, bottom=0.050)

    ax = fig.add_subplot(gs[0, 0])
    for u, v, t in EDGES:
        ax.plot(*zip(POS[u], POS[v]), "-", lw=1.4, color=C_TXT, zorder=1)
        mx = (POS[u][0] + POS[v][0]) / 2; my = (POS[u][1] + POS[v][1]) / 2
        ax.text(mx, my, f"{t:.1f}", fontsize=9.0, color=C_TXT, ha="center",
                va="center", bbox=dict(fc="white", ec="none", pad=0.5))
    for n in POS:
        ax.plot(*POS[n], "o", ms=6.2, color="#8c98a8", mec="white", mew=0.8, zorder=3)
    for n in SINK:
        ax.plot(*POS[n], "o", ms=14, mfc="none", mec=C_SNK, mew=1.5, zorder=2)
    ax.axhline(-0.03, color=C_SNK, lw=2.0)
    ax.text(0.5, -0.135, r"$S$  (sink side)", ha="center", fontsize=10.2, color=C_SNK)
    ax.set_xlim(-0.06, 1.02); ax.set_ylim(-0.14, 1.02)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(r"(a)  edge weights $\theta$", loc="left", fontsize=11.0)

    for k, th in enumerate(SNAP):
        ax = fig.add_subplot(gs[(k + 1) // 2, (k + 1) % 2])
        draw_net(ax, th, rf"(b{k+1})  $\theta={th}$", show_elder=(k == 0))
    nets_bottom = min(a.get_position().y0 for a in fig.axes)
    fig.text(0.5, nets_bottom - 0.004,
             "red: not yet connected to $S$   /   blue: connected to $S$",
             fontsize=10.2, color=C_TXT, ha="center", va="top")

    ax = fig.add_subplot(gs[2, :])
    births, deaths = [], []
    seen = {}
    order = sorted(EDGES, key=lambda e: e[2])
    parent = {n: n for n in POS}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    birth, has_s, bars = {}, {n: n in SINK for n in POS}, []
    for u, v, t in order:
        for n in (u, v):
            birth.setdefault(n, t)
        ru, rv = find(u), find(v)
        if ru == rv:
            continue
        su, sv = has_s[ru], has_s[rv]
        if su != sv:
            loser = rv if su else ru
            bars.append((birth[loser], t, "sink"))
        elif not su and not sv:
            young = ru if birth[ru] > birth[rv] else rv
            bars.append((birth[young], t, "elder"))
        b = min(birth[ru], birth[rv])
        parent[ru] = rv
        has_s[rv] = su or sv
        birth[rv] = b
    for y, (b, d, kind) in enumerate(sorted(bars, key=lambda x: (x[1], x[0]))):
        col = C_ELDER if kind == "elder" else C_SNK
        if d - b < 1e-9:
            ax.plot([d], [y], "|", ms=9, mew=2.0, color=col)
        else:
            ax.plot([b, d], [y, y], "-", lw=3.0, color=col, solid_capstyle="butt")
    for th in SNAP:
        ax.axvline(th, color="#b4bcc6", lw=0.9, ls=":")
        ax.text(th, len(bars) + 0.4, f"$\\theta$={th}", fontsize=9.2,
                color=C_LBL, ha="center")
    ax.set_xlim(-0.1, 4.0); ax.set_ylim(-1, len(bars) + 1.6)
    ax.set_yticks([]); ax.set_xlabel(r"$\theta=-\log_{10}g$", fontsize=10.6)
    ax.tick_params(labelsize=8.2)
    ax.set_title("(c)  the sink-relative barcode", loc="left", fontsize=11.0)
    nz = sum(1 for b, d, _ in bars if d - b < 1e-9)
    ax.legend(handles=[Line2D([], [], color=C_SNK, lw=3,
                              label="closed on reaching $S$"),
                       Line2D([], [], color=C_ELDER, lw=3,
                              label="closed by the elder rule"),
                       Line2D([], [], color=C_LBL, marker="|", ls="none", ms=9,
                              mew=2, label=f"zero length ({nz} of {len(bars)})")],
              frameon=False, fontsize=9.2, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)

    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig8_relative_persistence.{ext}", dpi=220,
                    bbox_inches="tight")
    print(f"bars: {len(bars)}  (sink {sum(1 for b in bars if b[2]=='sink')}, "
          f"elder {sum(1 for b in bars if b[2]=='elder')})")
    print("wrote", OUT / "fig8_relative_persistence.pdf")


if __name__ == "__main__":
    main()
