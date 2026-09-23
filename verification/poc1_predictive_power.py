# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/PoC1_predictive_power/poc1_predictive_power.py
#   sha256(src) : 442309ef603cfc47
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 39 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V10, V12. On one full-size network (N ~ 1000, |E| ~ 2500): first-order convergence of
# the dissipation share to d log G_eff / d log g_e (V10, Sec. 4.3), the finite-change
# closed form (Sec. 4.4), and the departure of the first-order prediction from a
# re-solve when the top k edges are improved by a factor C (V12, section 3 of the output).
# Writes results_poc1.txt.
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
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import synth_packing as SP
import thermal_sheaf_filtration as tsf
from thermal_sheaf_filtration import (Material, MaterialTable, AnalysisConfig,
                                      SOURCE_NODE, SINK_NODE)

SEED = 0
D_UM = SP.D_REPRESENTATIVE_UM

THETA = 0.0
OUTFIG = "fig_poc1.png"


def build(d_um=None):
    d_um = D_UM if d_um is None else d_um
    st = SP.build(d_um=d_um, seed=42)
    return st["G"], st["Gt"], st["sheaf"], st["cfg"]


def shorted_reff(Gt, theta=THETA):
    MERGED = -99
    H = nx.Graph()
    for u, v, d in Gt.edges(data=True):
        if d["g"] < theta:
            continue
        H.add_edge(u, v, g=d["g"])
    comp = nx.node_connected_component(H, SOURCE_NODE)
    Hc = H.subgraph(comp)
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
    nodes = list(Hs.nodes()); idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    L = np.zeros((n, n))
    for u, v, d in Hs.edges(data=True):
        i, j, g = idx[u], idx[v], d["g"]
        L[i, i] += g; L[j, j] += g; L[i, j] -= g; L[j, i] -= g
    P = np.linalg.pinv(L)
    dP = np.diag(P)

    def R(a, b):
        i, j = idx[a], idx[b]
        return float(dP[i] + dP[j] - 2 * P[i, j])
    return R, comp


def main():
    log = []

    def say(s=""):
        print(s); log.append(s)

    np.random.seed(SEED)
    t_start = time.time()

    G, Gt, sheaf, cfg = build()
    G_eff0 = tsf.effective_conductance(sheaf, THETA)
    flux = tsf.edge_flux_table(sheaf, THETA)
    say(f"structure: N={G.number_of_nodes()} E={G.number_of_edges()} "
        f"(terminal edges {Gt.number_of_edges()-G.number_of_edges()})")
    say(f"G_eff0 = {G_eff0:.12e} W/K,  physical edges {len(flux)}")
    say()

    flux["P_share_theory"] = flux["P"] / G_eff0
    say(f"sum P_share(CSV, normalised over physical edges) = {flux['P_share'].sum():.12f}")
    say(f"sum P_share(theory, /G_eff)    = {flux['P_share_theory'].sum():.12f}  "
        f"<- short of 1 by the terminal-edge part {1-flux['P_share_theory'].sum():.2e}")
    say()

    say("=== 0. difference between the shorted R~_e and the implemented R_eff_uv ===")
    Rfun, comp = shorted_reff(Gt)
    red, summ = tsf.edge_redundancy(sheaf, THETA, max_n=cfg.redundancy_max_n)
    red = red.set_index(["u", "v"])
    Rt, Rimpl = [], []
    for _, r in flux.iterrows():
        u, v = int(r["u"]), int(r["v"])
        Rt.append(Rfun(u, v))
        key = (u, v) if (u, v) in red.index else (v, u)
        Rimpl.append(red.loc[key, "R_eff_uv"] if key in red.index else np.nan)
    flux["R_tilde"] = Rt
    flux["R_impl"] = Rimpl
    ratio = (flux["R_impl"] / flux["R_tilde"]).dropna()
    say(f"  R_impl/R~: median {ratio.median():.6f} / max {ratio.max():.6f} "
        f"/ edges with >1.001: {(ratio>1.001).sum()} ({(ratio>1.001).mean()*100:.1f}%)")
    say(f"  -> the Sec. 4.4 closed form needs R~ (shorted). R~ is used from here on.")
    say()

    say("=== 1. first-order sensitivity  P_share_e = dlog G_eff/dlog g_e (Sec. 4.3, item 3) ===")
    say("  WARNING: on the metric: P_share spans 20 decades (min ~1e-20).")
    say("    for edges with P_share ~ 1e-9 or below, Delta log G_eff is buried in the sparse-solver rounding (~1e-12),")
    say("    so the relative error cannot be evaluated in principle. The verdict is therefore")
    say("    restricted to the edge set carrying 99.9% of the dissipation; the all-edge distribution is also shown.")
    cum = flux["P_share_theory"].cumsum() / flux["P_share_theory"].sum()
    n_eval = int((cum < 0.999).sum()) + 1
    say(f"    -> edges in the verdict set {n_eval} / {len(flux)} "
        f"(P_share >= {flux['P_share_theory'].iloc[n_eval-1]:.2e})")
    say()
    say(f"  {'delta':>8} | {'verdict set':>28} | {'all edges':>28}")
    say(f"  {'':>8} | {'max':>13} {'median':>13} | {'max':>13} {'median':>13}")
    sens_rows = []
    for delta in (1e-1, 1e-2, 1e-3):
        errs = np.empty(len(flux))
        t0 = time.time()
        for k, r in enumerate(flux.itertuples()):
            u, v = int(r.u), int(r.v)
            g0 = Gt[u][v]["g"]
            Gt[u][v]["g"] = g0 * (1 + delta)
            Gnew = tsf.effective_conductance(sheaf, THETA)
            Gt[u][v]["g"] = g0
            num = np.log(Gnew / G_eff0) / np.log(1 + delta)
            errs[k] = abs(num - r.P_share_theory) / max(r.P_share_theory, 1e-300)
        sens_rows.append((delta, errs))
        e = errs[:n_eval]
        say(f"  {delta:8.0e} | {np.max(e):13.3e} {np.median(e):13.3e} | "
            f"{np.max(errs):13.3e} {np.median(errs):13.3e}"
            f"   ({time.time()-t0:.1f}s)")
    d_arr = np.array([d for d, _ in sens_rows])
    e_med = np.array([np.median(e) for _, e in sens_rows])
    slope = np.polyfit(np.log(d_arr), np.log(e_med), 1)[0]
    e_max = np.array([np.max(e[:n_eval]) for _, e in sens_rows])
    slope_max = np.polyfit(np.log(d_arr), np.log(e_max), 1)[0]
    say(f"  log-log slope of the error in delta: median {slope:.3f} / max over the verdict set "
        f"{slope_max:.3f}   (1.0 means first-order accuracy as expected, i.e. a second-order residual)")
    ok1 = e_max[-1] < 1e-3
    say(f"  verdict (delta=1e-3, max rel.err over the verdict set < 1e-3): "
        f"{'PASS' if ok1 else 'FAIL'}  (= {e_max[-1]:.3e})")
    say()

    say("=== 2. closed form for a finite change (Sec. 4.4)  G' = G + Delta*dT^2/(1+Delta*R~_e) ===")
    n = len(flux)
    picks = [0, 1, 2, n // 4, n // 2, 3 * n // 4, n - 1]
    say(f"  {'edge':>14} {'rank':>6} {'P_share':>10} {'D/g':>6} "
        f"{'G_pred':>16} {'G_true':>16} {'rel.err':>10}")
    max_err_exact = 0.0
    max_err_impl = 0.0
    for rank in picks:
        r = flux.iloc[rank]
        u, v = int(r["u"]), int(r["v"])
        g0 = Gt[u][v]["g"]
        for c in (1, 10, 100):
            D = c * g0
            Gpred = G_eff0 + D * r["dT"] ** 2 / (1 + D * r["R_tilde"])
            Gt[u][v]["g"] = g0 + D
            Gtrue = tsf.effective_conductance(sheaf, THETA)
            Gt[u][v]["g"] = g0
            err = abs(Gpred - Gtrue) / Gtrue
            max_err_exact = max(max_err_exact, err)
            if np.isfinite(r["R_impl"]):
                Gp2 = G_eff0 + D * r["dT"] ** 2 / (1 + D * r["R_impl"])
                max_err_impl = max(max_err_impl, abs(Gp2 - Gtrue) / Gtrue)
            say(f"  {f'({u},{v})':>14} {rank:6d} {r['P_share']:10.3e} {c:6d} "
                f"{Gpred:16.9e} {Gtrue:16.9e} {err:10.2e}")
    say(f"  --> max relative error with R~ (shorted) = {max_err_exact:.3e}"
        f"   {'PASS (machine precision)' if max_err_exact < 1e-10 else 'FAIL'}")
    say(f"  --> with R_eff_uv (implemented, not shorted)   = {max_err_impl:.3e}"
        f"   <- direct evidence that the shorted version is required")
    say()

    say("=== 3. design rule: improving the top k edges at once (xC) ===")
    say(f"  {'k':>4} {'C':>4} {'top-k gain':>12} {'random-k':>12} {'bottom-k':>12} "
        f"{'top/rand':>10} {'1st-order':>12} {'pred/meas':>10}")
    design_rows = []
    rng = np.random.default_rng(SEED)
    for k in (1, 4, 10, 50):
        for C in (2, 10):
            def apply(idxs):
                saved = []
                for i in idxs:
                    u, v = int(flux.iloc[i]["u"]), int(flux.iloc[i]["v"])
                    saved.append((u, v, Gt[u][v]["g"]))
                    Gt[u][v]["g"] *= C
                out = tsf.effective_conductance(sheaf, THETA)
                for u, v, g in saved:
                    Gt[u][v]["g"] = g
                return out
            top = apply(range(k))
            bot = apply(range(n - k, n))
            rnds = [apply(rng.choice(n, k, replace=False))
                    for _ in range(_bootstrap.N_RANDOM_CONTROL)]
            rnd = float(np.mean(rnds))
            pred = G_eff0 * np.exp(flux["P_share_theory"].iloc[:k].sum() * np.log(C))
            gain_top = top / G_eff0 - 1
            gain_rnd = rnd / G_eff0 - 1
            gain_bot = bot / G_eff0 - 1
            say(f"  {k:4d} {C:4d} {gain_top:12.4%} {gain_rnd:12.4%} {gain_bot:12.4%} "
                f"{(gain_top/max(gain_rnd,1e-300)):10.1f} {pred/G_eff0-1:12.4%} "
                f"{(pred-G_eff0)/(top-G_eff0):10.3f}")
            design_rows.append((k, C, gain_top, gain_rnd, gain_bot,
                                pred / G_eff0 - 1))
    say("  -> pred/meas > 1 is the over-estimate of first-order theory (concavity and the sum rule)."
        "The larger k and C, the larger the departure: this is the range of validity of first-order theory.")
    say()

    say("=== 4. cross sensitivity of two edges d^2log G/dlog g_e dlog g_f (how first-order theory breaks) ===")
    say("  a single edge is concave (over-predicted), but improving top edges together can be under-predicted.")
    say("  edges in series raise each other's sensitivity (cross term > 0), which is physically natural.")
    say(f"  {'pair':>22} {'joint Dlog G':>13} {'sum (1st)':>13} {'cross term':>13} {'sign':>6}")
    for a in range(3):
        for b in range(a + 1, 4):
            ua, va = int(flux.iloc[a]["u"]), int(flux.iloc[a]["v"])
            ub, vb = int(flux.iloc[b]["u"]), int(flux.iloc[b]["v"])
            C = 2.0
            g_a, g_b = Gt[ua][va]["g"], Gt[ub][vb]["g"]
            Gt[ua][va]["g"] = g_a * C
            G_a = tsf.effective_conductance(sheaf, THETA)
            Gt[ub][vb]["g"] = g_b * C
            G_ab = tsf.effective_conductance(sheaf, THETA)
            Gt[ua][va]["g"] = g_a
            G_b = tsf.effective_conductance(sheaf, THETA)
            Gt[ub][vb]["g"] = g_b
            d_ab = np.log(G_ab / G_eff0)
            d_sum = np.log(G_a / G_eff0) + np.log(G_b / G_eff0)
            cross = d_ab - d_sum
            say(f"  {f'({ua},{va})x({ub},{vb})':>22} {d_ab:13.6f} {d_sum:13.6f} "
                f"{cross:13.3e} {'+' if cross>0 else '-':>6}")
    say()

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    ax = axes[0]
    for delta, errs in sens_rows:
        ax.hist(np.log10(np.maximum(errs, 1e-16)), bins=50, alpha=0.55,
                label=f"$\\delta$={delta:g}")
    ax.set_xlabel("log$_{10}$ relative error of $P^{share}$ vs finite difference")
    ax.set_ylabel("edges")
    ax.set_title(f"(a) 1st-order sensitivity\nslope of median error = {slope:.2f}")
    ax.legend()

    ax = axes[1]
    pr, tr = [], []
    sub = list(range(0, n, max(1, n // 300)))
    for i in sub:
        r = flux.iloc[i]
        u, v = int(r["u"]), int(r["v"])
        g0 = Gt[u][v]["g"]; D = 10 * g0
        pr.append(G_eff0 + D * r["dT"] ** 2 / (1 + D * r["R_tilde"]))
        Gt[u][v]["g"] = g0 + D
        tr.append(tsf.effective_conductance(sheaf, THETA))
        Gt[u][v]["g"] = g0
    pr = np.array(pr); tr = np.array(tr)
    ax.plot(tr / G_eff0, pr / G_eff0, "o", ms=3.5, alpha=0.6)
    lim = [min(tr.min(), pr.min()) / G_eff0, max(tr.max(), pr.max()) / G_eff0]
    ax.plot(lim, lim, "k--", lw=1)
    ax.set_xlabel("$G_{eff}$ recomputed / $G_{eff}^0$")
    ax.set_ylabel("$G_{eff}$ predicted (Cor 2.4) / $G_{eff}^0$")
    ax.set_title(f"(b) closed form, $\\Delta=10g_e$ ({len(sub)} edges)\n"
                 f"max rel.err = {np.max(np.abs(pr-tr)/tr):.1e}")

    ax = axes[2]
    ks = [1, 4, 10, 50]
    for C, mk in ((2, "o"), (10, "s")):
        gt = [d[2] for d in design_rows if d[1] == C]
        gr = [d[3] for d in design_rows if d[1] == C]
        gb = [d[4] for d in design_rows if d[1] == C]
        ax.plot(ks, gt, mk + "-", label=f"top-k, C={C}")
        ax.plot(ks, gr, mk + "--", label=f"random-k, C={C}", alpha=0.6)
        ax.plot(ks, gb, mk + ":", label=f"bottom-k, C={C}", alpha=0.6)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("k (number of improved edges)")
    ax.set_ylabel("$\\Delta G_{eff}/G_{eff}^0$")
    ax.set_title("(c) design rule: bottleneck ranking")
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(OUTFIG, dpi=150)
    say(f"figure saved: {OUTFIG}")
    say(f"total run time {time.time()-t_start:.1f}s")

    import datetime as _dt
    _st = tsf.engine_stamp()
    stamp = (f"\n=== run identification ===\n"
             f"  engine_version : {_st['engine_version']}\n"
             f"  run time        : {_dt.datetime.now().isoformat(timespec='seconds')}\n"
             f"  SEED / structure seed : {SEED} / 42 (deterministic; the same engine gives the same numbers)\n")

    with open("results_poc1.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(log) + "\n" + stamp)


if __name__ == "__main__":
    main()
