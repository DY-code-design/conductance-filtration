# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/mkfig2.py
#   sha256(src) : 6186ae6fd3e24804
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Fig. 2. Concept figure; no data is read.
from _branchfig import (plt, np, Circle, Rectangle, FancyBboxPatch,
                        C_F, C_M, C_I, C_T, sphere, save)


def arc(ax, cx, cy, r, t0, t1, ls, lw=1.0, color="#4a5568", z=5):
    th = np.radians(np.linspace(t0, t1, 240))
    ax.plot(cx + r * np.cos(th), cy + r * np.sin(th), ls, color=color, lw=lw,
            zorder=z)


def gap_fill(ax, xc, d, rr, y, y0, y1):
    Y = y + np.linspace(y0, y1, 160)
    dx = np.sqrt(np.maximum(rr ** 2 - (Y - y) ** 2, 0.0))
    ax.fill(np.concatenate([xc - d + dx, (xc + d - dx)[::-1]]),
            np.concatenate([Y, Y[::-1]]), fc="#ecd9ad", ec="none", zorder=4)


def draw(ax):
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, fc=C_M,
                           ec="none", alpha=0.40, zorder=0))
    ax.text(0.015, 0.938, "matrix", fontsize=8.6, color="#a07b32", va="top")
    y, rr = 0.80, 0.070
    D_NEAR, D_CON, D_BOND = 0.092, 0.0587, 0.0587
    R_OUT = rr * 0.90
    DXO = np.sqrt(rr ** 2 - R_OUT ** 2)
    cols = [
        (0.175, "near", "one annular\nintegral"),
        (0.500, "contact", "solid spot in parallel\nwith the matrix annulus"),
        (0.825, "bond", "axial bulk\nconduction"),
    ]
    for k, (xc, name, sub) in enumerate(cols):
        if k == 0:
            d = D_NEAR
            sphere(ax, xc - d, y, rr); sphere(ax, xc + d, y, rr)
            gap_fill(ax, xc, d, rr, y, -R_OUT, R_OUT)
            for s in (+1, -1):
                ax.plot([xc - d + DXO, xc + d - DXO], [y + s * R_OUT] * 2, ":",
                        color="#a07b32", lw=0.8, zorder=5)
            ax.annotate("", xy=(xc - (d - rr), y), xytext=(xc + (d - rr), y),
                        arrowprops=dict(arrowstyle="<->", color=C_I, lw=1.0),
                        zorder=7)
            ax.text(xc, y + 0.092, "$h$", ha="center", fontsize=9.0, color=C_I)

        elif k == 1:
            d = D_CON
            a = np.sqrt(rr ** 2 - d ** 2)
            delta = 2 * (rr - d)
            sphere(ax, xc - d, y, rr); sphere(ax, xc + d, y, rr)
            for s in (+1, -1):
                gap_fill(ax, xc, d, rr, y, s * a, s * R_OUT)
                ax.plot([xc - d + DXO, xc + d - DXO], [y + s * R_OUT] * 2, ":",
                        color="#a07b32", lw=0.8, zorder=5)
            Y = y + np.linspace(-a, a, 160)
            dx = np.sqrt(np.maximum(rr ** 2 - (Y - y) ** 2, 0.0))
            ax.fill(np.concatenate([xc + d - dx, (xc - d + dx)[::-1]]),
                    np.concatenate([Y, Y[::-1]]),
                    fc="#e8a49b", ec=C_I, lw=1.1, zorder=6)
            yd = y - rr - 0.026
            for sx in (-1, +1):
                ax.plot([xc + sx * delta / 2] * 2, [y, yd - 0.004], ":",
                        color=C_I, lw=0.7, zorder=7)
                ax.annotate("", xy=(xc + sx * (delta / 2 + 0.030), yd),
                            xytext=(xc + sx * delta / 2, yd),
                            arrowprops=dict(arrowstyle="<-", color=C_I, lw=0.9))
            ax.text(xc, yd - 0.026, "$\\delta$", ha="center", va="top",
                    fontsize=9.0, color=C_I)

        else:
            d = D_BOND
            a = np.sqrt(rr ** 2 - d ** 2)
            for sx in (-1, +1):
                ax.add_patch(Circle((xc + sx * d, y), rr, fc=C_F, ec="none",
                                    zorder=3))
            Y = y + np.linspace(-a, a, 160)
            dxl = np.sqrt(np.maximum(rr ** 2 - (Y - y) ** 2, 0.0))
            ax.fill(np.concatenate([xc + d - dxl, (xc - d + dxl)[::-1]]),
                    np.concatenate([Y, Y[::-1]]),
                    fc="#c3d0e2", ec="none", zorder=4)
            th = np.degrees(np.arctan2(a, d))
            arc(ax, xc - d, y, rr, th, 360 - th, "-")
            arc(ax, xc + d, y, rr, 180 + th, 540 - th, "-")
            arc(ax, xc - d, y, rr, -th, th, ":")
            arc(ax, xc + d, y, rr, 180 - th, 180 + th, ":")

        ax.text(xc, 0.640, name, ha="center", va="top", fontsize=9.4)
        ax.text(xc, 0.563, sub, ha="center", va="top", fontsize=8.4,
                linespacing=1.5)
    ax.set_xlim(0, 1); ax.set_ylim(0.460, 0.955); ax.set_aspect("equal")
    ax.axis("off")


def main():
    fig, ax = plt.subplots(figsize=(4.9, 2.423))
    draw(ax)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    save(fig, "fig2_three_branches")


if __name__ == "__main__":
    main()
