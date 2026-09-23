# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/T1_filtration_compare/t1_filtration_compare.py
#   sha256(src) : 8382fb6da59de0ec
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 64 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# Comparison of the conductance filtration (-log10 g) with two geometric
# filtrations on the same edge set: surface separation h and centre distance D
# (Sec. 3.1, Fig. 6, Fig. 7). Writes results_t1.csv; imported by
# y4_theta_separation.py (V15, V16).
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

import pathlib
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

import ensemble as EN
import materials_real as MR
import thermal_sheaf_filtration as tsf

D_UM = 20.0
SEEDS = (0, 1)
CASES = [
    ("phi_060",      "AlN", "epoxy", "untreated", "monodisperse, single material"),
    ("phi_064",      "AlN", "epoxy", "untreated", "monodisperse, single material (high loading)"),
    ("bimodal_3to1", "AlN", "epoxy", "untreated", "bimodal 3:1 (CV 0.47)"),
    ("bimodal_5to1", "AlN", "epoxy", "untreated", "bimodal 5:1 (CV 0.55)"),
    ("poly_wide",    "AlN", "epoxy", "untreated", "wide continuous distribution (CV 0.51)"),
    ("coat_mix",     "AlN", "epoxy", "untreated", "30% low-kappa coated particles"),
    ("clustered",    "AlN", "epoxy", "untreated", "agglomerated"),
    ("voided",       "AlN", "epoxy", "untreated", "spherical void"),
    ("stratified",   "AlN", "epoxy", "untreated", "stratified"),
    ("fiber_z",      "AlN", "epoxy", "untreated", "fibres along the heat flow"),
    ("phi_060",      "SiO2", "epoxy", "untreated", "monodisperse, low-kappa filler"),
    ("phi_060",      "AlN", "high_k_resin", "treated",  "monodisperse, high-kappa matrix with surface treatment"),
]
OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); OUT.append(s)


class _UF:
    def __init__(self, nodes):
        self.p = {n: n for n in nodes}
    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]; x = self.p[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return ra
        self.p[rb] = ra
        return ra


def filtration(nodes, edges, sink_node, source_node, absolute=False):
    es = sorted(edges, key=lambda e: e[2])
    uf = _UF(nodes)
    birth, has_sink = {}, {}
    birth_edge = {}
    bars, pairing = [], []
    theta_star, crit = None, None
    for k_e, (u, v, w) in enumerate(es):
        ru, rv = uf.find(u), uf.find(v)
        for r, nd in ((ru, u), (rv, v)):
            if r not in birth:
                birth[r] = w
                birth_edge[r] = k_e
                has_sink[r] = (nd == sink_node) and not absolute
        if ru == rv:
            continue
        su, sv = has_sink.get(ru, False), has_sink.get(rv, False)
        bu, bv = birth[ru], birth[rv]
        rn = uf.union(ru, rv)
        if su != sv:
            loser = rv if su else ru
            bars.append((birth[loser], w))
            pairing.append((birth_edge[loser], k_e))
        elif not su and not sv:
            young = ru if bu > bv else rv
            bars.append((birth[young], w))
            pairing.append((birth_edge[young], k_e))
        has_sink[rn] = su or sv
        if birth[ru] <= birth[rv]:
            birth[rn], birth_edge[rn] = birth[ru], birth_edge[ru]
        else:
            birth[rn], birth_edge[rn] = birth[rv], birth_edge[rv]
        if theta_star is None and uf.find(source_node) == uf.find(sink_node):
            theta_star, crit = w, (u, v)
    pers = [d - b for b, d in bars]
    roots = {uf.find(n) for n in nodes}
    return dict(n_bars=len(bars),
                max_pers=(max(pers) if pers else 0.0),
                theta_star=theta_star, crit_edge=crit,
                n_essential=len(roots),
                pairing=frozenset(pairing),
                bars=sorted(bars),
                edge_order=[(u, v) for u, v, _ in es])


def run_case(family, filler, matrix, treat, label, seed):
    st = EN.build(family, d_um=D_UM, seed=seed, treatment=treat,
                  matrix=matrix, filler=filler)
    if st is None:
        return None
    G, Gt, sh = st["G"], st["Gt"], st["sheaf"]
    pos, dia, info = st["pos"], st["dia"], st["info"]
    idx = {int(v): k for k, v in enumerate(st["ids"])}
    G0 = tsf.effective_conductance(sh, 0.0)
    flux = tsf.edge_flux_table(sh, 0.0)
    rank = {}
    for k, r in enumerate(flux.itertuples()):
        rank[frozenset((int(r.u), int(r.v)))] = k + 1

    real, w_h, w_D, w_C = [], [], [], []
    term = []
    for u, v, d in Gt.edges(data=True):
        if u < 0 or v < 0:
            term.append((u, v))
            continue
        i, j = idx[int(u)], idx[int(v)]
        Dij = float(np.linalg.norm(pos[i] - pos[j]))
        h = Dij - 0.5 * (dia[i] + dia[j])
        real.append((int(u), int(v)))
        w_h.append(max(h, 0.0))
        w_D.append(Dij)
        w_C.append(-float(np.log10(d["g"])))
    w_h, w_D, w_C = np.array(w_h), np.array(w_D), np.array(w_C)

    tau_D = float(kendalltau(w_D, w_C).statistic)
    tau_h = float(kendalltau(w_h, w_C).statistic)
    rho_D = float(spearmanr(w_D, w_C).statistic)
    frac_ov = float((w_h <= 0).mean())

    nodes = list(Gt.nodes())
    SRC, SNK = tsf.SOURCE_NODE, tsf.SINK_NODE
    res = {}
    for tag, w in (("G1", w_h), ("G2", w_D), ("C", w_C)):
        wmin = float(min(w.min(), 0.0)) - 1.0
        edges = [(u, v, float(x)) for (u, v), x in zip(real, w)]
        edges += [(u, v, wmin) for u, v in term]
        res[tag] = filtration(nodes, edges, SNK, SRC)

    out = dict(family=family, label=label, filler=filler, matrix=matrix,
               treat=treat, seed=seed, phi=float(info.get("phi", np.nan)),
               z_contact=float(info.get("z_contact", np.nan)),
               n_edges=len(real), frac_overlap=frac_ov,
               n_contact=int(info.get("n_contact", 0)),
               n_gap=int(info.get("n_gap", 0)),
               n_bond=int(info.get("n_bond", 0)),
               G_eff=G0, tau_G2_C=tau_D, tau_G1_C=tau_h, rho_G2_C=rho_D)
    for tag in ("G1", "G2", "C"):
        r = res[tag]
        out[f"nbars_{tag}"] = r["n_bars"]
        out[f"maxpers_{tag}"] = r["max_pers"]
        out[f"tstar_{tag}"] = r["theta_star"]
        ce = r["crit_edge"]
        out[f"crit_{tag}"] = str(ce)
        out[f"critrank_{tag}"] = (rank.get(frozenset(map(int, ce)), np.nan)
                                  if ce else np.nan)
    out["crit_same_G2_C"] = int(out["crit_G2"] == out["crit_C"])
    out["crit_same_G1_C"] = int(out["crit_G1"] == out["crit_C"])

    def _named(r):
        eo = r["edge_order"]
        return frozenset((eo[b], eo[d]) for b, d in r["pairing"])
    n2, nc = _named(res["G2"]), _named(res["C"])
    out["pair_jaccard_G2_C"] = len(n2 & nc) / max(len(n2 | nc), 1)
    out["pair_same_G2_C"] = int(n2 == nc)

    wmin = float(min(w_C.min(), 0.0)) - 1.0
    edges_abs = [(u, v, float(x)) for (u, v), x in zip(real, w_C)]
    edges_abs += [(u, v, wmin) for u, v in term]
    ab = filtration(nodes, edges_abs, SNK, SRC, absolute=True)
    out["nbars_abs"] = ab["n_bars"]
    out["maxpers_abs"] = ab["max_pers"]
    out["ncomp_end"] = ab["n_essential"]
    pa, pb = ab["pairing"], res["C"]["pairing"]
    out["pair_jac_abs_rel"] = len(pa & pb) / max(len(pa | pb), 1)
    ba = np.array(ab["bars"]); bb = np.array(res["C"]["bars"])
    out["diag_same_abs_rel"] = int(ba.shape == bb.shape and
                                   np.allclose(ba, bb, atol=1e-12))
    out["diag_maxdiff_abs_rel"] = (float(np.abs(ba - bb).max())
                                   if ba.shape == bb.shape else np.nan)

    for tag, pr in (("rel", "relative_h0"), ("leg", "legacy")):
        h = tsf.ConductanceFiltration(Gt, st["cfg"]).sweep(pairing=pr)
        pv = h["persistence_to_sink"].dropna().values
        sv = h["sink_theta"].dropna().values
        out[f"eng_nbars_{tag}"] = int(len(pv))
        out[f"eng_maxpers_{tag}"] = float(pv.max()) if len(pv) else 0.0
        out[f"eng_medpers_{tag}"] = float(np.median(pv)) if len(pv) else 0.0
        out[f"eng_tstar_{tag}"] = (float(-np.log10(sv.min()))
                                   if len(sv) else np.nan)
    return out


def main():
    rows = []
    for (fam, fil, mat, tr, lab) in CASES:
        for sd in SEEDS:
            try:
                r = run_case(fam, fil, mat, tr, lab, sd)
            except Exception as e:
                say(f"  !! {fam}/{fil}/{mat}/{tr}/seed{sd}: {type(e).__name__}: {e}")
                continue
            if r:
                rows.append(r)
                say(f"  {lab:28s} seed{sd}  n_e={r['n_edges']:5d} "
                    f"tau(G2,C)={r['tau_G2_C']:+.4f} tau(G1,C)={r['tau_G1_C']:+.4f} "
                    f"ov={r['frac_overlap']:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(pathlib.Path(__file__).with_name("results_t1.csv"), index=False)

    say("")
    say("=" * 100)
    say("T1: ordinary connected-component filtration (geometric) vs conductance filtration")
    say("=" * 100)
    say(f"structure: equivalent diameter d={D_UM:.0f} um, family catalogue, seed {SEEDS}, "
        f"edges {df['n_edges'].min()}-{df['n_edges'].max()}")
    say("")
    say("--- Table 1: agreement of edge order (Kendall tau) and degeneracy of G1 ---")
    say(f"{'case':30s} {'phi':>6s} {'edges':>6s} {'overlap':>10s} "
        f"{'tau(G2,C)':>10s} {'tau(G1,C)':>10s}")
    g = df.groupby(["label", "family", "filler", "matrix", "treat"], sort=False)
    for k, sub in g:
        say(f"{k[0]:30s} {sub['phi'].mean():6.3f} {sub['n_edges'].mean():6.0f} "
            f"{sub['frac_overlap'].mean():10.3f} {sub['tau_G2_C'].mean():+10.4f} "
            f"{sub['tau_G1_C'].mean():+10.4f}")
    say("")
    say("--- Table 2: percolation thresholds and critical edges ---")
    say("WARNING: the table below is a seed average. **Do not use it for representative values** (for same-configuration claims use seed 0 from the CSV)")
    say(f"{'case':30s} {'eps*(G1)':>10s} {'D*(G2)':>10s} {'theta*(C)':>12s} "
        f"{'crit same':>10s} {'crit-edge P rank (C)':>20s}")
    for k, sub in g:
        say(f"{k[0]:30s} {sub['tstar_G1'].mean():10.4f} "
            f"{sub['tstar_G2'].mean():10.3f} {sub['tstar_C'].mean():12.4f} "
            f"{sub['crit_same_G2_C'].mean():10.2f} "
            f"{sub['critrank_C'].mean():20.1f}")
    say("")
    say("--- Table 3: barcode (number of bars, longest bar) ---")
    say(f"{'case':30s} {'nbars G1':>9s} {'G2':>7s} {'C':>7s} "
        f"{'max G1':>9s} {'max G2':>9s} {'max C':>9s}")
    for k, sub in g:
        say(f"{k[0]:30s} {sub['nbars_G1'].mean():9.1f} {sub['nbars_G2'].mean():7.1f} "
            f"{sub['nbars_C'].mean():7.1f} {sub['maxpers_G1'].mean():9.4f} "
            f"{sub['maxpers_G2'].mean():9.3f} {sub['maxpers_C'].mean():9.3f}")
    say("")
    say("--- comparison with predictions ---")
    mono = df[(df.family.isin(["phi_060", "phi_064"]))]
    m1 = mono[(mono.filler == "AlN") & (mono.matrix == "epoxy")]
    say(f"T1-1 monodisperse, single material: tau(G2,C)=1.000: measured "
        f"{m1['tau_G2_C'].min():.6f}-{m1['tau_G2_C'].max():.6f} "
        f"-> {'a-priori prediction: met' if m1['tau_G2_C'].min() > 0.9999 else 'a-priori prediction: not met'}")
    hi = df[df.phi >= 0.55]
    say(f"T1-2 eps*(G1)=0 for phi>=0.55: measured max {hi['tstar_G1'].max():.4f} "
        f"-> {'a-priori prediction: met' if hi['tstar_G1'].max() <= 1e-12 else 'a-priori prediction: not met'}")
    bi = df[df.family == "bimodal_5to1"]
    say(f"T1-3 tau<0.95 for bimodal 5:1: measured {bi['tau_G2_C'].mean():.4f} "
        f"-> {'a-priori prediction: met' if bi['tau_G2_C'].mean() < 0.95 else 'a-priori prediction: not met'}")
    cm = df[df.family == "coat_mix"]
    say(f"T1-4 coat_mix: tau<1 and above bimodal: measured {cm['tau_G2_C'].mean():.4f} "
        f"-> {'a-priori prediction: met' if (cm['tau_G2_C'].mean() < 0.9999 and cm['tau_G2_C'].mean() > bi['tau_G2_C'].mean()) else 'a-priori prediction: not met'}")
    nb = df[["nbars_G1", "nbars_G2", "nbars_C"]].values
    say(f"T1-5 number of bars identical on the 3 axes: max per-row difference {int(np.abs(nb - nb[:, :1]).max())} "
        f"-> {'a-priori prediction: met' if np.abs(nb - nb[:, :1]).max() == 0 else 'a-priori prediction: not met'}")
    say(f"T1-6 critical-edge agreement rate (G2,C): overall {df['crit_same_G2_C'].mean():.3f} "
        f"/ monodisperse, single material {m1['crit_same_G2_C'].mean():.3f}")
    say(f"T1-7 spread of theta* across families: C {np.ptp(df['tstar_C']):.3f} decades / "
        f"G2 {np.ptp(df['tstar_G2']):.3f} um")
    say("")
    say("--- Table 4: agreement of the pairing (T1-8) and difference from ordinary H_0 (T1-9/10) ---")
    say(f"{'case':30s} {'pair jaccard(G2,C)':>18s} {'exact':>8s} "
        f"{'rel H0 == ordinary H0':>20s} {'pair jaccard':>12s} {'ncomp end':>9s}")
    for k, sub in g:
        say(f"{k[0]:30s} {sub['pair_jaccard_G2_C'].mean():18.4f} "
            f"{sub['pair_same_G2_C'].mean():8.2f} "
            f"{sub['diag_same_abs_rel'].mean():20.2f} "
            f"{sub['pair_jac_abs_rel'].mean():12.4f} "
            f"{sub['ncomp_end'].mean():9.1f}")
    say("")
    say("--- Table 5: comparison with the legacy setting that does not close bars (engine switch) ---")
    say(f"{'case':30s} {'nbars(rel)':>11s} {'(legacy)':>9s} "
        f"{'max(rel)':>10s} {'(legacy)':>9s} {'med(rel)':>10s} {'(legacy)':>9s}")
    for k, sub in g:
        say(f"{k[0]:30s} {sub['eng_nbars_rel'].mean():11.1f} "
            f"{sub['eng_nbars_leg'].mean():9.1f} "
            f"{sub['eng_maxpers_rel'].mean():10.3f} "
            f"{sub['eng_maxpers_leg'].mean():9.3f} "
            f"{sub['eng_medpers_rel'].mean():10.3f} "
            f"{sub['eng_medpers_leg'].mean():9.3f}")
    say("")
    say(f"T1-8 monodisperse, single material: pairing agrees exactly: measured "
        f"{m1['pair_same_G2_C'].mean():.2f} (agreement rate {m1['pair_jaccard_G2_C'].mean():.4f})"
        f" -> {'a-priori prediction: met' if m1['pair_same_G2_C'].mean() == 1.0 else 'a-priori prediction: not met'}")
    say(f"    pairing agreement rate for bimodal / coated mix: "
        f"bimodal_3to1 {df[df.family=='bimodal_3to1']['pair_jaccard_G2_C'].mean():.4f} / "
        f"bimodal_5to1 {bi['pair_jaccard_G2_C'].mean():.4f} / "
        f"poly_wide {df[df.family=='poly_wide']['pair_jaccard_G2_C'].mean():.4f} / "
        f"coat_mix {cm['pair_jaccard_G2_C'].mean():.4f}")
    say(f"T1-9/10 (after the replacement) relative H_0(X,S) and ordinary H_0 (elder rule, sink not treated"
        f" specially) agree: **persistence diagrams fully identical {df['diag_same_abs_rel'].mean():.3f}**"
        f" (max difference in bar values {df['diag_maxdiff_abs_rel'].max():.2e}), pairing agreement rate "
        f"{df['pair_jac_abs_rel'].mean():.4f}, "
        f"number of bars {df['nbars_C'].mean():.1f} = {df['nbars_abs'].mean():.1f}, "
        f"longest {df['maxpers_C'].mean():.4f} = {df['maxpers_abs'].mean():.4f}")
    say(f"    => because the virtual terminals are inserted first, the sink component is always the eldest;"
        f" hence relative H_0 coincides with standard H_0 (= the standard toolkit can be used as is)")
    say(f"T1-11 difference from the legacy setting that does not close bars: number of bars "
        f"{df['eng_nbars_leg'].mean():.1f} -> {df['eng_nbars_rel'].mean():.1f}"
        f" (x{df['eng_nbars_rel'].mean()/df['eng_nbars_leg'].mean():.2f}), "
        f"longest bar {df['eng_maxpers_leg'].mean():.3f} -> "
        f"{df['eng_maxpers_rel'].mean():.3f} decades")

    pathlib.Path(__file__).with_name("results_t1.txt").write_text(
        "\n".join(OUT) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
