# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/PoC4a_identities/poc4a_identities.py
#   sha256(src) : 1cdc6dfa2834f270
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 22 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V7. Identities checked on the quantities the implementation itself produces:
# sum_e P_e = G_eff, 1/R_eff(src, snk) = G_eff and Foster's theorem
# sum_e g_e R_eff(e) = n - 1, each with and without the terminal edges.
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
_bootstrap.tee_stdout("results_poc4a.txt")

import numpy as np
import networkx as nx

import synth_packing as SP
import thermal_sheaf_filtration as tsf
from thermal_sheaf_filtration import (Material, MaterialTable, AnalysisConfig,
                                      SOURCE_NODE, SINK_NODE)

SEED = 0
SEED_STRUCT = 42
D_UM = SP.D_REPRESENTATIVE_UM


def _laplacian_pinv(H, nodes):
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    L = np.zeros((n, n))
    for u, v, d in H.edges(data=True):
        i, j, g = idx[u], idx[v], d["g"]
        L[i, i] += g; L[j, j] += g; L[i, j] -= g; L[j, i] -= g
    return np.linalg.pinv(L), idx


def _reff(P, idx, a, b):
    i, j = idx[a], idx[b]
    return float(P[i, i] + P[j, j] - 2 * P[i, j])


def main():
    np.random.seed(SEED)

    st = SP.build(d_um=D_UM, seed=SEED_STRUCT)
    G, Gt, sheaf, cfg = st["G"], st["Gt"], st["sheaf"], st["cfg"]
    gv = st["g_virtual"]

    print(f"structure: N={G.number_of_nodes()} E={G.number_of_edges()} "
          f"(+2 virtual terminals, terminal edges {Gt.number_of_edges()-G.number_of_edges()})  "
          f"g_virtual={gv:.3e}")

    theta = 0.0
    T_map, Q = tsf.solve_dirichlet(sheaf, theta)
    G_eff = tsf.effective_conductance(sheaf, theta)
    print(f"G_eff = {G_eff:.12e} W/K   (Q={Q:.12e})\n")

    print("=== I1: sum_e P_e = G_eff (Sec. 4.3, item 2) ===")
    P_all = 0.0
    P_virtual = 0.0
    for u, v, d in Gt.edges(data=True):
        if d["g"] < theta:
            continue
        if u not in T_map or v not in T_map:
            continue
        p = d["g"] * (T_map[u] - T_map[v]) ** 2
        P_all += p
        if u in (SOURCE_NODE, SINK_NODE) or v in (SOURCE_NODE, SINK_NODE):
            P_virtual += p
    err_all = abs(P_all - G_eff) / G_eff
    print(f"  I1a all edges (terminal edges included)  sum P = {P_all:.12e}   rel.err = {err_all:.3e}"
          f"   {'PASS' if err_all < 1e-10 else 'FAIL'}")

    flux = tsf.edge_flux_table(sheaf, theta)
    P_real = float(flux["P"].sum())
    print(f"  I1b edge_flux_table (physical edges only) sum P = {P_real:.12e}   "
          f"sum P/G_eff = {P_real/G_eff:.9f}")
    print(f"      dissipation in the terminal edges = {P_virtual:.6e} "
          f"({P_virtual/G_eff*100:.4f}% of G_eff)")
    print(f"      -> P_share is normalised over the physical edges only (sum P_share = "
          f"{flux['P_share'].sum():.12f}).")
    print(f"      sum_e P_e = G_eff holds exactly with the terminal edges included; over the physical "
          f"edges alone, {P_virtual/G_eff*100:.3f}% of G_eff is dissipated in the terminal edges.\n")

    red, summ = tsf.edge_redundancy(sheaf, theta, max_n=cfg.redundancy_max_n)
    print("=== I2: 1/R_eff^{source-sink} = G_eff ===")
    Rss = summ["R_eff_source_sink"]
    err_i2 = abs(1.0 / Rss - G_eff) / G_eff
    print(f"  R_eff = {Rss:.12e}   1/R_eff = {1/Rss:.12e}   rel.err = {err_i2:.3e}"
          f"   {'PASS' if err_i2 < 1e-8 else 'FAIL'}\n")

    H = nx.Graph()
    for u, v, d in Gt.edges(data=True):
        if d["g"] >= theta:
            H.add_edge(u, v, g=d["g"])
    comp = nx.node_connected_component(H, SOURCE_NODE)
    Hc = H.subgraph(comp)
    nodes = list(Hc.nodes())
    Pinv, idx = _laplacian_pinv(Hc, nodes)
    n = len(nodes)

    print("=== I4: sum_e g_e R_eff(e) = n - 1 (Foster's theorem) ===")
    foster_all = sum(d["g"] * _reff(Pinv, idx, u, v) for u, v, d in Hc.edges(data=True))
    err_i4 = abs(foster_all - (n - 1)) / (n - 1)
    print(f"  I4a all edges  sum g R = {foster_all:.10f}   n-1 = {n-1}   rel.err = {err_i4:.3e}"
          f"   {'PASS' if err_i4 < 1e-8 else 'FAIL'}")
    crit_sum = float(red["crit"].sum())
    print(f"  I4b sum of the crit column (physical edges only) = {crit_sum:.10f}   "
          f"deficit = {n-1-crit_sum:.4f} = the {Hc.number_of_edges()-len(red)} terminal edges")
    print(f"      -> the implemented crit equals the Spielman-Srivastava leverage score. "
          f"The terminal edges must be added to check Foster.\n")

    print("=== X (pre-check for PoC-1): is edge_redundancy.R_eff_uv the value on the shorted circuit? ===")
    Hs = nx.Graph()
    MERGED = -99
    for u, v, d in Hc.edges(data=True):
        a = MERGED if u in (SOURCE_NODE, SINK_NODE) else u
        b = MERGED if v in (SOURCE_NODE, SINK_NODE) else v
        if a == b:
            continue
        if Hs.has_edge(a, b):
            Hs[a][b]["g"] += d["g"]
        else:
            Hs.add_edge(a, b, g=d["g"])
    nodes_s = list(Hs.nodes())
    Pinv_s, idx_s = _laplacian_pinv(Hs, nodes_s)

    sub = red.nlargest(5, "crit")
    print(f"  {'edge':>16} {'R_eff_uv(impl.)':>18} {'R~_e(shorted)':>18} {'ratio':>10}")
    for _, r in sub.iterrows():
        u, v = int(r["u"]), int(r["v"])
        Rt = _reff(Pinv_s, idx_s, u, v)
        print(f"  {f'({u},{v})':>16} {r['R_eff_uv']:18.6e} {Rt:18.6e} "
              f"{r['R_eff_uv']/Rt:10.4f}")
    ratios = []
    for _, r in red.iterrows():
        u, v = int(r["u"]), int(r["v"])
        ratios.append(r["R_eff_uv"] / _reff(Pinv_s, idx_s, u, v))
    ratios = np.array(ratios)
    same = np.allclose(ratios, 1.0, rtol=1e-8)
    print(f"  ratio over all {len(ratios)} edges: median {np.median(ratios):.6f} / "
          f"min {ratios.min():.6f} / max {ratios.max():.6f}")
    print(f"  -> R_eff_uv {'equals' if same else '**does not equal**'} the shorted version."
          f"{'' if same else ' The Sec. 4.4 closed form needs the shorted R~_e computed separately (implemented in PoC-1).'}")


if __name__ == "__main__":
    main()
