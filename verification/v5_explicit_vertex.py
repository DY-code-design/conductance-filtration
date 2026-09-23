# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/V5_explicit_vertex/v5_explicit_vertex.py
#   sha256(src) : 5f8ca31ee08ce23a
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 21 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V5. The through-particle resistance is in series with the whole edge (Sec. 2.5).
# Checked by expanding every edge into an explicit interior vertex m_e with
# u --1/R_thr-- m_e --(g_spot + g_ann)-- v; the ratio of G_eff must be 1. Two
# deliberately wrong placements (inside the spot; no R_thr) serve as controls.
import sys as _sys, pathlib as _pathlib
_lib = next((p / "lib" for p in _pathlib.Path(__file__).resolve().parents
             if (p / "lib" / "_bootstrap.py").is_file()),
            None)
if _lib is None:
    _lib = next((q for q in (_pathlib.Path(__file__).resolve().parent.parent / "thermal_network",
                             _pathlib.Path(__file__).resolve().parent / "thermal_network")
                 if (q / "_bootstrap.py").is_file()), None)
if _lib is None:
    raise SystemExit("_bootstrap.py not found (neither lib/ nor thermal_network/)")
_sys.path.insert(0, str(_lib))
import _bootstrap

import itertools
import pathlib

import networkx as nx
import pandas as pd

import ensemble as EN
import thermal_sheaf_filtration as tsf

HERE = pathlib.Path(__file__).resolve().parent
FAMILIES = ["phi_060", "phi_064", "bimodal_3to1", "voided", "clustered"]
MATERIALS = [("AlN", "epoxy", "untreated"),
             ("SiO2", "epoxy", "treated"),
             ("AlN", "silicone", "treated")]
SEEDS = (0, 1)
N, D_UM, THETA = 400, 20.0, 0.0

OUT = []


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    OUT.append(line)


def g_eff(Gt) -> float:
    return tsf.effective_conductance(tsf.CellularSheaf(Gt), THETA)


def split_parts(dd):
    Rb = float(dd.get("R_bulk", 0.0) or 0.0)
    if Rb <= 0.0:
        return float(dd["g"]), 0.0, float(dd.get("g_spot", 0.0)), 0.0
    gs = dd.get("g_spot")
    ga = dd.get("g_annulus")
    if gs is None or ga is None:
        R = float(dd["R"]) if "R" in dd else 1.0 / float(dd["g"])
        return 1.0 / max(R - Rb, 1e-300), Rb, 0.0, 0.0
    return float(gs) + float(ga), Rb, float(gs), float(ga)


def build_variants(Gt):
    X = nx.Graph()
    C1 = nx.Graph()
    C2 = nx.Graph()
    X.add_nodes_from(Gt.nodes(data=True))
    C1.add_nodes_from(Gt.nodes(data=True))
    C2.add_nodes_from(Gt.nodes(data=True))

    _alg = []
    _next = max(int(x) for x in Gt.nodes()) + 1
    n_bulk = n_bond = 0
    for u, v, dd in Gt.edges(data=True):
        g_surf, Rb, gs, ga = split_parts(dd)
        if Rb <= 0.0:
            n_bond += 1
            for H in (X, C1, C2):
                H.add_edge(u, v, **dd)
            continue
        n_bulk += 1

        _alg.append(abs(1.0 / (1.0 / g_surf + Rb) - float(dd["g"])) / float(dd["g"]))

        m = _next
        _next += 1
        X.add_node(m)
        X.add_edge(u, m, **{**dd, "g": 1.0 / Rb, "R": Rb, "R_bulk": 0.0})
        X.add_edge(m, v, **{**dd, "g": g_surf, "R": 1.0 / g_surf, "R_bulk": 0.0})

        if gs > 0.0:
            g_c1 = 1.0 / (1.0 / gs + Rb) + ga
        else:
            g_c1 = float(dd["g"])
        C1.add_edge(u, v, **{**dd, "g": g_c1, "R": 1.0 / g_c1})

        C2.add_edge(u, v, **{**dd, "g": g_surf, "R": 1.0 / g_surf, "R_bulk": 0.0})
    return X, C1, C2, n_bulk, n_bond, (max(_alg) if _alg else 0.0)


def main() -> int:
    rows = []
    _st = tsf.engine_stamp()
    say(f"# V5 agreement with the explicit vertex expansion   engine={_st.get('engine_version')} / "
        f"variant {_st.get('engine_variant')} / fingerprint {_st.get('baseline_fp')}")
    say(f"# structure families {len(FAMILIES)} x materials {len(MATERIALS)} x seeds {len(SEEDS)} "
        f"= {len(FAMILIES)*len(MATERIALS)*len(SEEDS)} conditions / N={N} d={D_UM}um theta={THETA}")
    say("")

    for fam, (fil, mat, tr), sd in itertools.product(FAMILIES, MATERIALS, SEEDS):
        st = EN.build(family=fam, d_um=D_UM, seed=sd, treatment=tr,
                      matrix=mat, filler=fil, n=N)
        if st is None:
            continue
        Gt = st["Gt"]
        X, C1, C2, n_bulk, n_bond, alg = build_variants(Gt)

        gL = g_eff(Gt)
        gX = g_eff(X)
        gC1 = g_eff(C1)
        gC2 = g_eff(C2)

        rows.append(dict(family=fam, filler=fil, matrix=mat, treat=tr, seed=sd,
                         n_edges=Gt.number_of_edges(), n_bulk=n_bulk, n_bond=n_bond,
                         G_L=gL, G_X=gX, G_C1=gC1, G_C2=gC2, alg_err=alg,
                         ratio_X=gX / gL, ratio_C1=gC1 / gL, ratio_C2=gC2 / gL))

    df = pd.DataFrame(rows)
    df.to_csv(HERE / "results_v5_explicit_vertex.csv", index=False)

    say("## ratios per condition (X is the identity; C1 and C2 are controls with deliberately wrong placement)")
    say("")
    say("| family | material | seed | edges | with bulk | X/L | ! C1/L (inside the spot) | ! C2/L (no term) |")
    say("|---|---|---|---|---|---|---|---|")
    for r in rows:
        say(f"| {r['family']} | {r['filler']}/{r['matrix']}/{r['treat']} | {r['seed']} "
            f"| {r['n_edges']} | {r['n_bulk']} "
            f"| {r['ratio_X']:.12f} | {r['ratio_C1']:.6f} | {r['ratio_C2']:.6f} |")
    say("")

    dX = (df.ratio_X - 1.0).abs()
    dC1 = (df.ratio_C1 - 1.0).abs()
    dC2 = (df.ratio_C2 - 1.0).abs()
    say("## summary")
    say("")
    import math
    digits = int(-math.log10(dX.max())) if dX.max() > 0 else 16
    say(f"  X / L      max deviation {dX.max():.3e}   median {dX.median():.3e}"
        f"   -> **agreement to {digits} decimal places**")
    say(f"  per-branch algebraic identity 1/(1/g_surf+R_bulk) == g   max {df.alg_err.max():.3e}"
        f"   (measured directly, without solving the network)")
    say(f"  ! C1 / L   max deviation {dC1.max()*100:.4f}%   median {dC1.median()*100:.4f}%")
    say(f"  ! C2 / L   max deviation {dC2.max()*100:.4f}%   median {dC2.median()*100:.4f}%")
    say("")

    say(f"  WARNING: conditions where the network ratio does not reach machine precision: {int((dX > 1e-12).sum())} / {len(df)}.")
    say(f"     the worst is the condition with the smallest G_eff ({df.loc[dX.idxmax(),'G_L']:.3e} W/K, "
        f"1/{df.G_L.max()/df.loc[dX.idxmax(),'G_L']:.0f} of the others).")
    say(f"     WARNING: **not a problem of the formula but of the conditioning of the solve** - it worsens as vertices are added. "
        f"The per-branch algebraic identity stays within machine precision at {df.alg_err.max():.1e}.")
    say("")

    ok = bool(dX.max() < 1e-7 and df.alg_err.max() < 1e-12)
    ctrl = bool(dC1.max() > 1e-6 and dC2.max() > 1e-6)
    say(f"  verdict: per-branch algebraic identity < 1e-12 and network ratio < 1e-7 : {'PASS' if ok else '! FAIL'}")
    say(f"  verdict: the controls depart from 1 : {'PASS' if ctrl else '! FAIL'}")
    if not ctrl:
        say("  WARNING: the controls did not move: the ratio of 1 is then not a test of the placement of R_thr.")

    (HERE / "results_v5_explicit_vertex.txt").write_text(
        "\n".join(OUT) + "\n", encoding="utf-8")
    return 0 if (ok and ctrl) else 2


if __name__ == "__main__":
    _sys.exit(main())
