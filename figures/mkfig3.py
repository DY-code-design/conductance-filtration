# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig3.py
#   sha256(src) : ea34c6fed71321be
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 3. Concept figure; no data is read.
from _branchfig import (plt, np, Circle, Rectangle, FancyBboxPatch,
                        C_I, C_T, C_B, wire, vwire, res, node, save)

FS_TAG = 9.2
FS_RES = 8.6
FS_NOTE = 8.4


def draw(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0.075, 1.0); ax.axis("off")
    xL, xR = 0.070, 0.935
    W_THR, X_THR = 0.104, 0.0875
    X_BUS = 0.175

    y = 0.545
    ax.text(0.012, y + 0.185, r"(b)  contact:  $g_{\rm contact}=g_{\rm spot}+g_{\rm ann}$", fontsize=FS_TAG, color=C_T)
    node(ax, xL, y, "$i$"); node(ax, xR, y, "$j$", side="right")
    wire(ax, xL, xL + X_BUS, y)
    res(ax, xL + X_THR, y, W_THR, r"$R_{{\rm thr},i}$", fc="#efe8f8", ec=C_B,
        fs=7.4, hh=0.042)
    wire(ax, xR - X_BUS, xR, y)
    res(ax, xR - X_THR, y, W_THR, r"$R_{{\rm thr},j}$", fc="#efe8f8", ec=C_B,
        fs=7.4, hh=0.042)
    xa, xb = xL + X_BUS, xR - X_BUS
    vwire(ax, xa, y - 0.112, y + 0.095); vwire(ax, xb, y - 0.112, y + 0.095)
    wire(ax, xa, xb, y + 0.095)
    for k, L in enumerate([r"$R_{{\rm con},i}$", r"$R_{\rm int}$", r"$R_{{\rm con},j}$"]):
        res(ax, xa + 0.094 + k * 0.164, y + 0.095, 0.140, L, fs=FS_RES, hh=0.040)
    ax.text(xb + 0.012, y + 0.095, "spot\n$1/g_{\\rm spot}$", fontsize=7.8,
            color="#7d8794", va="center", linespacing=1.25)
    wire(ax, xa, xb, y - 0.112)
    res(ax, (xa + xb) / 2, y - 0.112, 0.250, r"$R_{\rm ann}$",
        fc="#fdf3e2", ec="#c9a227", fs=FS_RES, hh=0.040)
    ax.text(xb + 0.012, y - 0.112, "annulus\n$1/g_{\\rm ann}$", fontsize=7.8,
            color="#a07b32", va="center", linespacing=1.25)

    y = 0.885
    ax.text(0.012, y + 0.090, r"(a)  near-contact:  $g_{\rm prox}$", fontsize=FS_TAG, color=C_T)
    node(ax, xL, y, "$i$"); node(ax, xR, y, "$j$", side="right")
    wire(ax, xL, xR, y)
    res(ax, xL + X_THR, y, W_THR, r"$R_{{\rm thr},i}$", fc="#efe8f8", ec=C_B,
        fs=7.4, hh=0.042)
    res(ax, xR - X_THR, y, W_THR, r"$R_{{\rm thr},j}$", fc="#efe8f8", ec=C_B,
        fs=7.4, hh=0.042)
    res(ax, (xL + xR) / 2, y, 0.540,
        r"$R_{\rm gap}\,(\,h+a_{K,i}\!+\!a_{K,j}+h_{\rm BOB}\,)$",
        fc="#eaf1fa", ec="#2c6fbb", fs=FS_RES, hh=0.046)

    y = 0.205
    ax.text(0.012, y + 0.090, r"(c)  bond:  $g_{\rm bond}$", fontsize=FS_TAG, color=C_T)
    node(ax, xL, y, "$i$"); node(ax, xR, y, "$j$", side="right")
    wire(ax, xL, xR, y)
    res(ax, (xL + xR) / 2, y, 0.520, r"$R_{\rm axial}=1/g_{\rm bond}=L/(\kappa_h\pi R_{\min}^2)$",
        fc="#eef2f6", ec="#4a5568", fs=FS_RES, hh=0.046)
    ax.text(0.50, y - 0.080,
            r"no $R_{\rm thr}$, no interface",
            ha="center", va="top", fontsize=FS_NOTE, color="#4a5568")


def main():
    fig, ax = plt.subplots(figsize=(5.0, 2.958))
    draw(ax)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.99, bottom=0.02)
    save(fig, "fig3_branch_circuits")


if __name__ == "__main__":
    main()
