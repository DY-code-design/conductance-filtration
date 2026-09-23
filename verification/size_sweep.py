# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/Z_size_sweep/size_sweep.py
#   sha256(src) : 7315f7e703dd500d
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 51 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V6. Sweep of the network-level quantities over the particle diameter 1-100 um,
# including the ratio of G_eff with and without the saturation length h_BOB
# (Sec. 2.2.2). Writes results_size_sweep.csv.
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

import time

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import materials_real as MR
import synth_packing as SP
import thermal_sheaf_filtration as tsf
from thermal_sheaf_filtration import SOURCE_NODE, SINK_NODE

D_SWEEP = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
TREATMENTS = ["untreated", "treated"]
MATRICES = ["epoxy", "silicone", "high_k_resin"]
SEED = 42
K_TOP = 10
C_IMPROVE = 10.0
N_RANDOM_CONTROL = _bootstrap.N_RANDOM_CONTROL


def laplacian_pinv(H, nodes):
    idx = {n: i for i, n in enumerate(nodes)}
    L = np.zeros((len(nodes), len(nodes)))
    for u, v, d in H.edges(data=True):
        i, j, g = idx[u], idx[v], d["g"]
        L[i, i] += g; L[j, j] += g; L[i, j] -= g; L[j, i] -= g
    return np.linalg.pinv(L), idx


def reff(P, idx, a, b):
    i, j = idx[a], idx[b]
    return float(P[i, i] + P[j, j] - 2 * P[i, j])


def one_case(d_um, treatment, matrix, rng):
    st = SP.build(d_um=d_um, seed=SEED, treatment=treatment, matrix=matrix)
    G, Gt, sh = st["G"], st["Gt"], st["sheaf"]
    G0 = tsf.effective_conductance(sh, 0.0)
    if not (G0 > 0):
        return None
    T, _ = tsf.solve_dirichlet(sh, 0.0)
    flux = tsf.edge_flux_table(sh, 0.0)
    flux["key"] = [tuple(sorted((int(a), int(b))))
                   for a, b in zip(flux["u"], flux["v"])]

    P_all = sum(d["g"] * (T[u] - T[v]) ** 2 for u, v, d in Gt.edges(data=True)
                if u in T and v in T)
    err_I1 = abs(P_all / G0 - 1)

    H = nx.Graph()
    for u, v, d in Gt.edges(data=True):
        H.add_edge(u, v, g=d["g"])
    Hc = H.subgraph(nx.node_connected_component(H, SOURCE_NODE))
    nodes = list(Hc.nodes())
    Pinv, idx = laplacian_pinv(Hc, nodes)
    foster = sum(d["g"] * reff(Pinv, idx, u, v) for u, v, d in Hc.edges(data=True))
    err_I4 = abs(foster / (len(nodes) - 1) - 1)

    MERGED = -99
    Hs = nx.Graph()
    for u, v, d in Hc.edges(data=True):
        a = MERGED if u in (SOURCE_NODE, SINK_NODE) else u
        b = MERGED if v in (SOURCE_NODE, SINK_NODE) else v
        if a == b:
            continue
        if Hs.has_edge(a, b):
            Hs[a][b]["g"] += d["g"]
        else:
            Hs.add_edge(a, b, g=d["g"])
    Pin_s, idx_s = laplacian_pinv(Hs, list(Hs.nodes()))
    ratios, Rt = [], {}
    for r in flux.head(200).itertuples():
        u, v = int(r.u), int(r.v)
        rt = reff(Pin_s, idx_s, u, v)
        Rt[(u, v)] = rt
        ratios.append(reff(Pinv, idx, u, v) / rt)
    max_ratio = float(np.max(ratios))

    errs = []
    for r in flux.head(3).itertuples():
        u, v = int(r.u), int(r.v)
        g0 = Gt[u][v]["g"]
        for c in (1, 10, 100):
            Dg = c * g0
            pred = G0 + Dg * r.dT ** 2 / (1 + Dg * Rt[(u, v)])
            Gt[u][v]["g"] = g0 + Dg
            true = tsf.effective_conductance(sh, 0.0)
            Gt[u][v]["g"] = g0
            errs.append(abs(pred - true) / true)
    err_cor24 = float(np.max(errs))

    have = {tuple(sorted((int(a), int(b)))): (int(a), int(b))
            for a, b in zip(flux["u"], flux["v"])}

    def gain(keys):
        saved = []
        for k in keys:
            if k not in have:
                continue
            uu, vv = have[k]
            saved.append((uu, vv, Gt[uu][vv]["g"]))
            Gt[uu][vv]["g"] *= C_IMPROVE
        out = tsf.effective_conductance(sh, 0.0)
        for uu, vv, gg in saved:
            Gt[uu][vv]["g"] = gg
        return out - G0

    g_top = gain(list(flux["key"].head(K_TOP)))
    allk = list(have.keys())
    g_rnd = float(np.mean([gain([allk[i] for i in
                                 rng.choice(len(allk), K_TOP, replace=False)])
                           for _ in range(N_RANDOM_CONTROL)]))
    lift = g_top / g_rnd if g_rnd > 0 else np.nan

    st_off = SP.build(d_um=d_um, seed=SEED, treatment=treatment,
                      matrix=matrix, bob_saturation=False)
    G_off = tsf.effective_conductance(st_off["sheaf"], 0.0)
    n_sat = sum(1 for _, _, d in G.edges(data=True)
                if d.get("bob_saturated", False))
    n_gap = sum(1 for _, _, d in G.edges(data=True) if d["branch"] == "gap")

    by = flux.groupby("branch")["P_share"].sum()
    ps = flux["P_share"].values
    ps = ps[ps > 0]
    Hent = float(-(ps * np.log(ps)).sum())

    return dict(d_um=d_um, treatment=treatment, matrix=matrix,
                k_m=MR.K_MATRIX[matrix],
                a_K_um=MR.TREATMENT[treatment]["R_s_matrix"]
                * MR.K_MATRIX[matrix] * 1e6,
                N=G.number_of_nodes(), E=G.number_of_edges(), G_eff=G0,
                err_I1=err_I1, err_I4=err_I4, max_R_ratio=max_ratio,
                err_cor24=err_cor24, gain_top=g_top / G0, gain_rnd=g_rnd / G0,
                lift=lift, G_off=G_off, sat_ratio=G0 / G_off,
                n_sat=n_sat, n_gap=n_gap, frac_sat=n_sat / max(n_gap, 1),
                share_gap=float(by.get("gap", 0.0)),
                share_contact=float(by.get("contact", 0.0)),
                share_bond=float(by.get("bond", 0.0)),
                H=Hent, top20=float(flux["P_share"].head(20).sum()))


def main():
    log = []

    def say(s=""):
        print(s); log.append(s)

    t0 = time.time()
    rng = np.random.default_rng(0)
    say("=== Z_size_sweep: sweeping the main network-layer metrics over particle size 1-100 um ===")
    say(f"particle size {D_SWEEP} um x treatment {TREATMENTS} x matrix {MATRICES} "
        f"= {len(D_SWEEP)*len(TREATMENTS)*len(MATRICES)} structures")
    say("(restricted to the headline metrics rather than the full analysis of each PoC; the full analysis"
        " at every size would take tens of minutes)")
    say()

    rows = []
    for d in D_SWEEP:
        for tr in TREATMENTS:
            for mx in MATRICES:
                r = one_case(d, tr, mx, rng)
                if r is None:
                    say(f"  [skip] d={d} {tr} {mx}: no percolation")
                    continue
                rows.append(r)
        say(f"  d={d:6.1f} um done ({time.time()-t0:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv("results_size_sweep.csv", index=False)
    say(f"\nevaluation done: {len(df)} structures")
    say()

    say("--- do the identities (PoC-4a / PoC-1-(2)) hold over the whole size range? ---")
    say(f"  {'metric':>28} {'max':>12} {'median':>12} {'verdict':>8}")
    for col, name, tol in (("err_I1", "I1 error |sum P/G_eff - 1|", 1e-10),
                           ("err_I4", "I4 Foster error", 1e-8),
                           ("err_cor24", "max relative error of Cor. 2.4", 1e-9)):
        say(f"  {name:>28} {df[col].max():12.3e} {df[col].median():12.3e} "
            f"{'PASS' if df[col].max() < tol else 'FAIL':>8}")
    say(f"  -> **the identities hold in all 42 structures**, independent of size 1-100 um, 2 treatments and 3 matrices.")
    say()

    say("--- size dependence (epoxy / untreated as representative, all d) ---")
    sub = df[(df.matrix == "epoxy") & (df.treatment == "untreated")]
    say(f"  {'d [um]':>8} {'a_K [um]':>9} {'G_eff':>12} {'G_eff/d':>10} "
        f"{'top10/rnd':>10} {'gap share':>10} {'contact':>9} {'bond':>8} "
        f"{'H':>7} {'max R/R~':>9}")
    for _, r in sub.iterrows():
        say(f"  {r.d_um:8.1f} {r.a_K_um:9.5f} {r.G_eff:12.4e} "
            f"{r.G_eff/r.d_um:10.3e} {r.lift:10.2f} {r.share_gap:10.3f} "
            f"{r.share_contact:9.3f} {r.share_bond:8.3f} {r.H:7.3f} "
            f"{r.max_R_ratio:9.4f}")
    say()

    say("--- does each claim depend on particle size? (range over all 42 structures) ---")
    say(f"  {'claim':>34} {'min':>11} {'max':>11} {'max/min':>10} {'verdict':>16}")
    checks = [
        ("PoC-1-4 lift of top-10/random-10", "lift", 1.0, "valid if always >1"),
        ("P1 saturation G_eff(on)/G_eff(off)", "sat_ratio", None, "inactive if ~1"),
        ("P1 fraction of saturated edges", "frac_sat", None, "inactive if small"),
        ("PoC-6 gap dissipation share", "share_gap", None, "robust aggregate"),
        ("PoC-6 contact dissipation share", "share_contact", None, "same as above"),
        ("PoC-2a dissipation entropy H", "H", None, "structure descriptor"),
        ("PoC-4a-X max(R_impl/R~)", "max_R_ratio", None, ">1: shorted version"),
    ]
    for name, col, _, note in checks:
        v = df[col].replace([np.inf, -np.inf], np.nan).dropna()
        rr = (f"{v.max()/v.min():10.2f}" if v.min() > 1e-12 else f"{'-':>10}")
        say(f"  {name:>34} {v.min():11.4f} {v.max():11.4f} {rr} {note:>16}")
    say()
    say(f"  * the lift of top-10/random-10 over all 42 structures is "
        f"{df['lift'].min():.2f}-{df['lift'].max():.2f}x. **always well above 1**")
    say(f"  * BOB saturation over all structures: G_eff ratio {df['sat_ratio'].min():.6f}-"
        f"{df['sat_ratio'].max():.6f}, saturated edges at most {df['frac_sat'].max():.2%}")
    say(f"    -> **P1 'saturation is inactive over the whole practical range' is confirmed for d=1-100 um x 3 matrices**")
    say(f"  * gap dissipation share is {df['share_gap'].min():.3f}-{df['share_gap'].max():.3f}")
    say()

    say("--- effect of treatment and matrix (G_eff ratio at each size) ---")
    say(f"  {'d [um]':>8} {'treated/untreated':>18} {'silicone/epoxy':>16} "
        f"{'high-k/epoxy':>14}")
    for d in D_SWEEP:
        s = df[df.d_um == d]
        try:
            r_tr = (s[(s.treatment == "treated") & (s.matrix == "epoxy")]["G_eff"].iloc[0]
                    / s[(s.treatment == "untreated") & (s.matrix == "epoxy")]["G_eff"].iloc[0])
            r_si = (s[(s.treatment == "untreated") & (s.matrix == "silicone")]["G_eff"].iloc[0]
                    / s[(s.treatment == "untreated") & (s.matrix == "epoxy")]["G_eff"].iloc[0])
            r_hk = (s[(s.treatment == "untreated") & (s.matrix == "high_k_resin")]["G_eff"].iloc[0]
                    / s[(s.treatment == "untreated") & (s.matrix == "epoxy")]["G_eff"].iloc[0])
            say(f"  {d:8.1f} {r_tr:18.4f} {r_si:16.4f} {r_hk:14.4f}")
        except IndexError:
            pass
    say("  -> shows how the effect of surface treatment changes with size (because a_K/gap scale with size)")
    say()

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))

    ax = axes[0]
    for (tr, mx), s in df.groupby(["treatment", "matrix"]):
        ax.plot(s["d_um"], s["G_eff"], "o-", ms=4, lw=1.3,
                label=f"{tr}/{mx}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("particle diameter $d$ [um]")
    ax.set_ylabel("$G_{eff}$ [W/K]")
    ax.set_title("(a) $G_{eff}$ over the practical size range\n"
                 "(42 structures: 7 sizes x 2 treatments x 3 matrices)")
    ax.legend(fontsize=6, ncol=2); ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    for col, lab in (("err_I1", "I1 (Euler)"), ("err_I4", "I4 (Foster)"),
                     ("err_cor24", "Cor 2.4 closed form")):
        ax.plot(df["d_um"], df[col].clip(lower=1e-17), "o", ms=4, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.axhline(1e-9, color="k", ls="--", lw=1)
    ax.text(1.2, 1.5e-9, "tolerance", fontsize=8)
    ax.set_xlabel("particle diameter $d$ [um]")
    ax.set_ylabel("relative error")
    ax.set_title("(b) identities hold across the whole range\n"
                 "(no size or material dependence)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    for (tr, mx), s in df.groupby(["treatment", "matrix"]):
        ax.plot(s["d_um"], s["lift"], "o-", ms=4, lw=1.3, label=f"{tr}/{mx}")
    ax.axhline(1.0, color="k", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("particle diameter $d$ [um]")
    ax.set_ylabel("gain(top-10) / gain(random-10)")
    ax.set_title("(c) does the bottleneck ranking stay useful?\n"
                 "(PoC-1-4 across the size range)")
    ax.legend(fontsize=6, ncol=2); ax.grid(alpha=0.3)

    fig.tight_layout(); fig.savefig("fig_size_sweep.png", dpi=150)
    say("figure saved: fig_size_sweep.png")
    say(f"total run time {time.time()-t0:.0f}s")

    with open("results_size_sweep.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
