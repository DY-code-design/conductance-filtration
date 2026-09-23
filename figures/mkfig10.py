# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig10.py
#   sha256(src) : c4fa0fa804bdfef8
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 10. The pseudocode is transcribed from the engine (sweep,
# percolation_threshold, solve_dirichlet, shorted_edge_resistance);
# this script only draws it.
import pathlib

import matplotlib

import _pdffont
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = pathlib.Path(__file__).parent
import _paths
OUT = _paths.out_dir()
PREVIEW = _paths.preview_dir()

C_EDGE, C_FACE, C_TITLE, C_CODE = "#c8d0da", "#f7f9fb", "#2d3748", "#1f2933"

BARCODE = r"""input: E with θ = -log10 g;  sink set S
sort E by g descending   # θ ascending
for (u, v, θ) in E:
    ru, rv <- find(u), find(v)
    for r in {ru, rv} with no birth yet:
        birth[r]   <- θ
        at_sink[r] <- (r in S)
    if ru == rv: continue
    if at_sink[ru] != at_sink[rv]:
        d <- the root not at the sink
    else if neither is at the sink:
        d <- the later-born root   # elder
    else: d <- none
    if d: emit bar (birth[d], θ)
    r <- union(ru, rv)
    at_sink[r] <- at_sink[ru] or at_sink[rv]
    birth[r]   <- min(birth[ru], birth[rv])"""

THETA = r"""input: E with θ;  source, sink
sort E by g descending   # θ ascending
for (u, v, θ) in E:
    union(u, v)
    if find(source) == find(sink):
        θ* <- θ        # g* = 10^-θ*
        return θ*"""

SHARE = r"""input: E with g;  source and sink
   I = interior nodes, B = boundary
L <- laplacian(g)
u[B] <- 1 at source, 0 at sink
solve L[I,I] u[I] = -L[I,B] u[B]
for e = (i, j) in E:
    P[e] <- g[e] * (u[i] - u[j])^2
P_share <- P / sum(P)"""

FINITE = r"""input: E with g;  k;  change Δ
merge source and sink into one node
L_s <- laplacian(g) of that graph
factorise L_s once
for e = (i, j) in top-k by P_share:
    b_e <- unit(i) - unit(j)
    x   <- backsolve(L_s, b_e)
    Rt[e] <- b_e^T x      # R tilde
    du_e  <- u[i] - u[j]  # from (c)
    G_new[e] <- G_eff(g)
        + Δ du_e^2 / (1 + Δ Rt[e])"""

PANELS = [
    ("(a)  sink-relative barcode", BARCODE),
    ("(b)  critical conductance level", THETA),
    ("(c)  dissipation share", SHARE),
    ("(d)  finite change for the top $k$", FINITE),
]


def box(fig, x, y, w, h, title, code, fs, tfs):
    fig.patches.append(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.010",
        transform=fig.transFigure, facecolor=C_FACE, edgecolor=C_EDGE,
        linewidth=0.55, zorder=0))
    fig.text(x + 0.012, y + h - 0.022, title, fontsize=tfs, color=C_TITLE,
             va="top", ha="left")
    fig.text(x + 0.016, y + h - 0.062, code, fontsize=fs, color=C_CODE,
             va="top", ha="left", family=MONO, linespacing=LINESP)


MONO = "Noto Sans Mono"
FS, TFS = 7.4, 9.0
LINESP = 1.62
W_IN = 6.25
PAD_TOP_IN = 0.30
PAD_BOT_IN = 0.13
GAP_IN = 0.10


def _width_in(code):
    f = plt.figure(figsize=(8, 8))
    tx = f.text(0, 0, code, fontsize=FS, family=MONO, linespacing=LINESP)
    bb = tx.get_window_extent(f.canvas.get_renderer())
    plt.close(f)
    return bb.width / f.dpi


def _need_in(code):
    f = plt.figure(figsize=(8, 8))
    t = f.text(0, 0, code, fontsize=FS, family=MONO, linespacing=LINESP)
    bb = t.get_window_extent(f.canvas.get_renderer())
    plt.close(f)
    return PAD_TOP_IN + bb.height / f.dpi + PAD_BOT_IN


def main():
    COLS = [[0, 1], [2, 3]]
    need = [_need_in(c) for _, c in PANELS]
    tot = [sum(need[i] for i in col) + GAP_IN * (len(col) - 1) for col in COLS]
    H_IN = max(tot)
    fig = plt.figure(figsize=(W_IN, H_IN))
    for col, x, w in zip(COLS, (0.005, 0.545), (0.520, 0.450)):
        rest = H_IN - tot[COLS.index(col)]
        gap_extra = 0.0
        base = sum(need[i] for i in col)
        y = 1.0
        for i in col:
            h = (need[i] + rest * need[i] / base) / H_IN
            y -= h
            inner = w * W_IN - 0.028 * W_IN - 0.06
            got = _width_in(PANELS[i][1])
            assert got <= inner, (
                f"{PANELS[i][0]}: the code is {got:.2f} in wide, past the {inner:.2f} in "
                f"inside the box. Wrap the line or widen the column; do not delete it")
            box(fig, x, y, w, h, *PANELS[i], fs=FS, tfs=TFS)
            y -= (GAP_IN + gap_extra) / H_IN
    for ext, dst in (("pdf", OUT), ("png", PREVIEW)):
        fig.savefig(dst / f"fig10_pseudocode.{ext}", dpi=220,
                    bbox_inches="tight", pad_inches=0.02, facecolor="white")
    print(f"wrote {OUT / 'fig10_pseudocode.pdf'}  ({W_IN:.2f} x {H_IN:.2f} in)")


if __name__ == "__main__":
    main()
