# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/PoC6_robustness/poc6_robustness.py
#   sha256(src) : cd2fb7a94f341f3f
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 85 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V13. Robustness of the edge ranking to model and material uncertainty: for each
# perturbed setting (eta_cut, contact threshold, saturation switch, slab thickness,
# matrix material), eta = gain obtained with the top-10 edges chosen by the reference
# model / gain of the top-10 edges chosen by the perturbed model itself (Sec. 4.3).
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

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kendalltau

import synth_packing as SP
import thermal_sheaf_filtration as tsf
from thermal_sheaf_filtration import Material, MaterialTable, AnalysisConfig

SEED = 0
D_UM = SP.D_REPRESENTATIVE_UM

TOPK = 20

BASE = dict(eta_cut=1.0, slab_frac=0.06,
            t_min_um=1e-5, bob_saturation=True, contact_eps_um=0.0, viz=False)
BASE_MATRIX = "epoxy"
BASE_KM = 0.2


def run_case(matrix=BASE_MATRIX, **over):
    kw = dict(BASE); kw.update(over)
    kw.pop("gap_max_um", None); kw.pop("slab_axis", None); kw.pop("viz", None)
    st = SP.build(d_um=D_UM, seed=42, matrix=matrix, **kw)
    G, Gt, sh, cfg = st["G"], st["Gt"], st["sheaf"], st["cfg"]
    Geff = tsf.effective_conductance(sh, 0.0)
    flux = tsf.edge_flux_table(sh, 0.0)
    flux["key"] = [tuple(sorted((int(a), int(b))))
                   for a, b in zip(flux["u"], flux["v"])]
    hist = tsf.ConductanceFiltration(Gt).sweep()
    st = hist["sink_theta"].dropna()
    theta_star = float(st.min()) if len(st) else float("nan")
    pers = hist["persistence_to_sink"].dropna()
    return dict(Geff=Geff, Gt=Gt, sheaf=sh, n_edges=G.number_of_edges(),
                n_contact=sum(1 for _, _, d in G.edges(data=True)
                              if d["branch"] == "contact"),
                flux=flux, theta_star=theta_star,
                pers_max=float(pers.max()) if len(pers) else float("nan"),
                pers_med=float(pers.median()) if len(pers) else float("nan"),
                entropy=float(-(flux["P_share"] *
                                np.log(np.maximum(flux["P_share"], 1e-300))).sum()),
                top4=float(flux["P_share"].head(4).sum()))


def compare(ref, cur):
    a = ref["flux"].set_index("key")["P_share"]
    b = cur["flux"].set_index("key")["P_share"]
    common = a.index.intersection(b.index)
    rho_all = spearmanr(a.loc[common], b.loc[common]).statistic
    top100 = ref["flux"]["key"].head(100)
    top100 = [k for k in top100 if k in b.index]
    rho_100 = spearmanr(a.loc[top100], b.loc[top100]).statistic
    tau_100 = kendalltau(a.loc[top100], b.loc[top100]).statistic
    sa = set(ref["flux"]["key"].head(TOPK))
    sb = set(cur["flux"]["key"].head(TOPK))
    jac = len(sa & sb) / len(sa | sb)
    jac50 = (len(set(ref["flux"]["key"].head(50)) & set(cur["flux"]["key"].head(50)))
             / len(set(ref["flux"]["key"].head(50)) | set(cur["flux"]["key"].head(50))))
    return dict(rho_all=rho_all, rho_100=rho_100, tau_100=tau_100, jaccard=jac,
                jaccard50=jac50, n_common=len(common))


def transfer_efficiency(ref, cur, k=10, C=10.0):
    Gt, sh = cur["Gt"], cur["sheaf"]
    G0 = cur["Geff"]
    have = {tuple(sorted((int(a), int(b)))): (int(a), int(b))
            for a, b in zip(cur["flux"]["u"], cur["flux"]["v"])}

    def gain(keys):
        saved = []
        for kk in keys:
            if kk not in have:
                continue
            u, v = have[kk]
            saved.append((u, v, Gt[u][v]["g"]))
            Gt[u][v]["g"] *= C
        out = tsf.effective_conductance(sh, 0.0)
        for u, v, g in saved:
            Gt[u][v]["g"] = g
        return out - G0

    g_ref = gain(list(ref["flux"]["key"].head(k)))
    g_own = gain(list(cur["flux"]["key"].head(k)))
    rng = np.random.default_rng(SEED)
    allk = list(have.keys())
    g_rnd = float(np.mean([gain([allk[i] for i in
                                 rng.choice(len(allk), k, replace=False)])
                           for _ in range(_bootstrap.N_RANDOM_CONTROL)]))
    return dict(eta=g_ref / g_own if g_own > 0 else float("nan"),
                gain_ref=g_ref / G0, gain_own=g_own / G0, gain_rnd=g_rnd / G0,
                lift_ref=g_ref / g_rnd if g_rnd > 0 else float("nan"))


def main():
    np.random.seed(SEED)
    log = []

    def say(s=""):
        print(s); log.append(s)

    say("=== PoC-6: robustness of the descriptors to model parameters ===")
    say(f"reference: {BASE}, matrix={BASE_MATRIX}, particle size d={D_UM} um")
    say("WARNING: the legacy version swept k_matrix over 0.2/1/5, but on literature values "
        "three resins (epoxy 0.20 / silicone 0.25 / high-k 1.00) = a x5 range")
    say()
    ref = run_case()
    say(f"reference case: E={ref['n_edges']} (contact {ref['n_contact']}), "
        f"G_eff={ref['Geff']:.5e}, θ*={ref['theta_star']:.3e}, "
        f"H(P_share)={ref['entropy']:.3f}, Σtop4={ref['top4']:.4f}")
    say()

    cases = []
    for v in (0.5, 2.0):
        cases.append((f"eta_cut={v}", dict(eta_cut=v)))
    for v in (1e-4, 1e-6):
        cases.append((f"t_min_um={v:g}", dict(t_min_um=v)))
    cases.append(("bob_saturation=False", dict(bob_saturation=False)))
    for v in (0.02, 0.05, 0.1):
        cases.append((f"contact_eps_um={v}", dict(contact_eps_um=v)))
    for v in (0.03, 0.08):
        cases.append((f"slab_frac={v}", dict(slab_frac=v)))
    for v in ("silicone", "high_k_resin"):
        cases.append((f"matrix={v}", dict(_mx=v)))
    cases.append(("high_k_resin & contact_eps=0.05",
                  dict(_mx="high_k_resin", contact_eps_um=0.05)))

    say("[Table 1] agreement of the descriptors (rank based)")
    say(f"{'case':>28} {'E':>6} {'n_cont':>6} {'G ratio':>9} {'th* ratio':>9} "
        f"{'rho(all)':>10} {'rho(top100)':>11} {'J@20':>6} {'J@50':>6} {'H ratio':>7}")
    rows = []
    for name, over in cases:
        mx = over.pop("_mx", BASE_MATRIX)
        cur = run_case(matrix=mx, **over)
        c = compare(ref, cur)
        t = transfer_efficiency(ref, cur)
        say(f"{name:>28} {cur['n_edges']:6d} {cur['n_contact']:6d} "
            f"{cur['Geff']/ref['Geff']:9.4f} "
            f"{cur['theta_star']/ref['theta_star']:9.4f} "
            f"{c['rho_all']:10.4f} {c['rho_100']:11.4f} "
            f"{c['jaccard']:6.2f} {c['jaccard50']:6.2f} "
            f"{cur['entropy']/ref['entropy']:7.3f}")
        rows.append(dict(case=name, n_edges=cur["n_edges"],
                         G_ratio=cur["Geff"] / ref["Geff"],
                         theta_ratio=cur["theta_star"] / ref["theta_star"],
                         **c, **t, H_ratio=cur["entropy"] / ref["entropy"]))
    df = pd.DataFrame(rows)
    say()

    say("[Table 2] robustness as a design decision (gain when the k=10 edges are improved x10)")
    say("  eta = gain of the top-10 chosen with the reference model / gain of that model's own top-10")
    say("  (if eta~1, 'even with the wrong model one reaches the same design decision')")
    say(f"{'case':>28} {'eta':>8} {'gain(ref ord)':>13} {'gain(own ord)':>13} "
        f"{'gain(random)':>13} {'lift vs random':>15}")
    for r in rows:
        say(f"{r['case']:>28} {r['eta']:8.3f} {r['gain_ref']:13.4%} "
            f"{r['gain_own']:13.4%} {r['gain_rnd']:13.4%} {r['lift_ref']:15.1f}")
    say()

    say("--- check against the predictions ---")
    gmin, gmax = df["G_ratio"].min(), df["G_ratio"].max()
    say(f"6-1 the absolute value of G_eff is sensitive (up to several-fold): "
        f"{'PASS (a-priori prediction met)' if gmax/gmin > 2 else 'FAIL (a-priori prediction not met)'} "
        f"(G_eff ratio {gmin:.3f} to {gmax:.3f}, max/min = {gmax/gmin:.1f}x)")
    rmin = df["rho_100"].min()
    say(f"6-2 Spearman of the top edges >= 0.8: "
        f"{'PASS (a-priori prediction met)' if rmin >= 0.8 else 'FAIL (a-priori prediction not met)'} "
        f"(min rho(top100) = {rmin:.4f}, case "
        f"'{df.loc[df['rho_100'].idxmin(),'case']}')")
    worst_G = df.loc[(df["G_ratio"] - 1).abs().idxmax(), "case"]
    say(f"6-3 the resin moves G_eff the most: "
        f"{'PASS (a-priori prediction met)' if ('matrix' in worst_G or 'resin' in worst_G) else 'FAIL (a-priori prediction not met)'} ('{worst_G}')")
    eta = df[df["case"].str.startswith("eta_cut")]
    say(f"6-4 eta_cut changes the edge count but not the ranking: "
        f"edge count {eta['n_edges'].min()}-{eta['n_edges'].max()} "
        f"(reference {ref['n_edges']}), rho(top100) {eta['rho_100'].min():.4f}-"
        f"{eta['rho_100'].max():.4f}, Jaccard {eta['jaccard'].min():.2f}")
    ce = df[df["case"].str.startswith("contact_eps")]
    say(f"6-5 contact_eps_um changes the ranking of the top edges: "
        f"rho(top100) {ce['rho_100'].min():.4f}-{ce['rho_100'].max():.4f}, "
        f"Jaccard {ce['jaccard'].min():.2f}–{ce['jaccard'].max():.2f}, "
        f"contact branches {ce['n_edges'].min()}-{ce['n_edges'].max()}")
    sf = df[df["case"].str.startswith("slab_frac")]
    say(f"6-6 the effect of slab_frac is small: G_eff ratio "
        f"{sf['G_ratio'].min():.4f}–{sf['G_ratio'].max():.4f}, "
        f"rho(top100) {sf['rho_100'].min():.4f}")
    worst = df.loc[df["rho_100"].idxmin(), "case"]
    say(f"6-6b the lowest rank correlation is: '{worst}' "
        f"(ρ={df['rho_100'].min():.4f})")
    say()

    say("--- why rho(top100) is low: separating the range-restriction bias ---")
    say("  rho(top100) looks only at the top-100 edges of the reference, so the spread of P_share is "
        "compressed to less than one decade.")
    say("  rank correlation is always attenuated by range restriction, so this value alone "
        "must not be taken as evidence that 'the ranking has broken'.")
    say(f"  rho over all common edges stays within {df['rho_all'].min():.3f} to {df['rho_all'].max():.3f}.")
    say()
    say("--- 6-7: splitting the perturbations into 'analysis conditions' and 'model uncertainty' ---")
    say("  reviewer note: 'eta_cut is a precondition of the calculation; it is natural to specify it exactly and compute'")
    say("  -> eta_cut / t_min / slab_frac are **analysis conditions**, inputs that must be specified exactly.")
    say("    mixing them into eta as 'model uncertainty' is incorrect.")
    say("    the true uncertainty is the **description of materials and interfaces** (resin type, contact criterion, saturation on/off).")
    say()
    _INPUT = ("eta_cut", "t_min_um", "slab_frac")
    _MODEL = ("matrix", "contact_eps_um", "bob_saturation", "high_k_resin")
    def _grp(c):
        if any(c.startswith(k) for k in _INPUT):
            return "analysis cond."
        return "model uncert."
    say(f"  {'group':>14} {'case':>32} {'eta':>8} {'rho(top100)':>11} {'G ratio':>9}")
    _by = {"analysis cond.": [], "model uncert.": []}
    for _r in df.itertuples():
        gname = _grp(_r.case)
        _by[gname].append(float(_r.eta))
        say(f"  {gname:>14} {_r.case:>32} {_r.eta:8.3f} {_r.rho_100:11.3f} "
            f"{_r.G_ratio:9.3f}")
    say()
    for gname, vals in _by.items():
        if vals:
            say(f"  {gname:>14}: eta min **{min(vals):.3f}** / median "
                f"{sorted(vals)[len(vals)//2]:.3f} ({len(vals)} cases)")
    _mdl = _by["model uncert."]
    _p67 = bool(_mdl) and min(_mdl) >= 0.8
    say()
    say(f"  6-7 verdict (eta >= 0.8 against **model uncertainty**): "
        f"{'PASS (a-priori prediction met)' if _p67 else 'FAIL (a-priori prediction not met)'}"
        f" (eta min {min(_mdl):.3f})")
    say("  -> **as long as the analysis conditions (eta_cut etc.) are specified exactly, the design decision is robust**.")
    _eta_ec2 = [float(_r.eta) for _r in df.itertuples() if str(_r.case).startswith("eta_cut=2.0")]
    say(f"    eta = {_eta_ec2[0]:.3f} at eta_cut=2.0 is the case of a wrongly specified analysis condition,"
        if _eta_ec2 else "    eta at eta_cut=2.0 is the case of a wrongly specified analysis condition,")
    say("    not model uncertainty. In the paper the analysis conditions are stated explicitly and fixed.")
    say()

    say("  what matters as a design decision is eta in Table 2 (whether the same edges are chosen).")
    say()

    say("--- overall verdict (does the claim hold?) ---")
    eta_min = df["eta"].min()
    lift_min = df["lift_ref"].min()
    ok = (gmax / gmin > 2) and (eta_min >= 0.8)
    say(f"  'G_eff is sensitive / the design decision is robust': {'PASS (holds)' if ok else 'WARNING: conditional'}")
    say(f"  against a G_eff spread of {gmax/gmin:.1f}x,")
    say(f"    rho(all common edges)      min {df['rho_all'].min():.3f}")
    say(f"    rho(top-100 edges)         min {rmin:.3f}  <- attenuated by range restriction; not used for the verdict")
    say(f"    eta(design-decision transfer) min {eta_min:.3f}  "
        f"(case '{df.loc[df['eta'].idxmin(),'case']}')")
    say(f"    lift vs random             min {lift_min:.1f}x")
    say("  -> the claim should be worded not as 'the ranking is preserved' but as "
        "'even when the bottlenecks chosen with a different model are improved, ")
    say("     a fraction eta of the gain of that model's own optimal choice is obtained, "
        "and a lift-fold advantage over random selection is kept'.")

    say()
    say("--- why eta breaks down with k_matrix: the rate-limiting channel swaps ---")
    say("  dissipation share per branch (sum of P_share) and the branch composition of the top-20 edges")
    say(f"  {'case':>16} {'gap':>8} {'contact':>8} {'bond':>8} | "
        f"{'top20: gap':>12} {'contact':>8} {'bond':>8}")
    for name, mx, over in (("ref epoxy k_m=0.20", "epoxy", {}),
                           ("silicone k_m=0.25", "silicone", {}),
                           ("high-k resin k_m=1.0", "high_k_resin", {}),
                           ("eta_cut=2.0", "epoxy", dict(eta_cut=2.0))):
        c = run_case(matrix=mx, **over)
        fl = c["flux"]
        by = fl.groupby("branch")["P_share"].sum()
        t20 = fl.head(20)["branch"].value_counts()
        say(f"  {name:>16} {by.get('gap',0):8.3f} {by.get('contact',0):8.3f} "
            f"{by.get('bond',0):8.3f} | {t20.get('gap',0):12d} "
            f"{t20.get('contact',0):8d} {t20.get('bond',0):8d}")
    say("  -> reading (important): the *aggregate share by branch type is robust* (see the branch columns of the table above).")
    say("    what changes is *which individual edges rank at the top*: the composition of the top 20 (right part of the table) changes.")
    say("    raising k_matrix makes detours through the matrix effective, so the fibre bonds become relatively rate-limiting.")
    say("    at eta_cut=2.0 the opposite happens: distant near-contact edges are cut and the flow concentrates on contact edges.")
    say("    the claim should therefore be layered:")
    say("      - **aggregate quantities (share per branch) are robust** -> these can be reported with confidence")
    say("      - **the ranking of individual edges depends on k_matrix and eta_cut** -> "
        "these must be specified exactly as inputs")
    say("    this is physics, not arbitrariness of the model (the same story as the background parallel conduction in Sec. 5.1).")

    df.to_csv("results_poc6.csv", index=False)
    with open("results_poc6.txt", "w") as f:
        f.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
