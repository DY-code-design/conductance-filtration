# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/_data.py
#   sha256(src) : 6fbbd723c7014820
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Packing generation and edge extraction for Fig. 4, Fig. 5 and Fig. 6.
import numpy as np
import _paths
_paths.add_lib()
import _bootstrap
import ensemble as EN
import L0_engine.thermal_sheaf_filtration as tsf


def build(family="phi_060", filler="AlN", matrix="epoxy", treat="untreated",
          n=400, seed=0, d_um=20.0):
    d = EN.build(family=family, d_um=d_um, seed=seed, treatment=treat,
                 matrix=matrix, filler=filler, n=n)
    return d


def edge_table(d):
    Gt = d["Gt"]; pos = d["pos"]; dia = d["dia"]; ids = d["ids"]
    idx = {int(a): k for k, a in enumerate(ids)}
    rows = []
    for u, v, dd in Gt.edges(data=True):
        if u < 0 or v < 0:
            continue
        i, j = idx[int(u)], idx[int(v)]
        D = float(np.linalg.norm(pos[i] - pos[j]))
        h = D - 0.5 * (dia[i] + dia[j])
        rows.append(dict(u=int(u), v=int(v), h=h, D=D,
                         g=float(dd["g"]), branch=dd.get("branch", "")))
    return rows


if __name__ == "__main__":
    d = build()
    r = edge_table(d)
    print(len(r), "edges;  overlap frac =",
          round(sum(1 for x in r if x["h"] <= 0) / len(r), 3))
    print(sorted(set(x["branch"] for x in r)))
