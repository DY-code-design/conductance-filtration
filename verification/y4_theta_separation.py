# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/Y4_theta_separation/y4_theta_separation.py
#   sha256(src) : 36770bd7ea246df3
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 58 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V15, V16. Pairs of structures from the 23 families with nearly equal G_eff:
# whether theta_star separates them (V15) and the spread of the idealised gains
# of Sec. 4.5 across the same pairs (V16). Writes results_y4_pairs.csv.
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

import sys
import time

import numpy as np
import pandas as pd

import ensemble as EN
import materials_real as MR
import quantities24 as Q
import thermal_sheaf_filtration as tsf

_bootstrap.add_case_path("T1_filtration_compare")
from t1_filtration_compare import filtration

D_UM = 20.0
FILLER, MATRIX = "AlN", "epoxy"
SEEDS = (0, 1)
K_TOP, BOOST = 20, 10.0
N_RANDOM_CONTROL = _bootstrap.N_RANDOM_CONTROL
TOL = 0.01

PAIRS = [
    (1, "uniform vs agglomerated", "phi_060", "clustered"),
    (1, "uniform vs spherical void", "phi_060", "voided"),
    (1, "monodisperse vs bimodal 5:1", "phi_064", "bimodal_5to1"),
    (1, "bimodal vs wide continuous", "bimodal_3to1", "poly_wide"),
    (1, "uniform vs stratified (strong)", "phi_060", "strat_180"),
    (1, "fibres z vs fibres xy", "fiber_z", "fiber_xy"),
    (2, "stratified weak vs medium", "strat_025", "strat_055"),
    (2, "stratified medium vs strong", "strat_055", "strat_100"),
    (2, "stratified strong vs very strong", "strat_100", "strat_180"),
    (2, "segregation none vs weak", "seg_000", "seg_033"),
    (2, "segregation weak vs medium", "seg_033", "seg_067"),
    (2, "segregation medium vs full", "seg_067", "seg_100"),
    (3, "phi0.50 vs phi0.60", "phi_050", "phi_060"),
    (3, "phi0.55 vs phi0.64", "phi_055", "phi_064"),
    (3, "phi0.50 vs phi0.64", "phi_050", "phi_064"),
    (3, "phi0.64 vs bimodal0.70", "phi_064", "bimodal_5to1"),
]

CSV = _pathlib.Path("results_y4_pairs.csv")
TXT = _pathlib.Path("results_y4_pairs.txt")

DUP_LABEL = "monodisperse vs bimodal 5:1"


def build(fam, seed, knob):
    t = MR.TREATMENT["untreated"]
    tab = tsf.MaterialTable(
        materials={k: tsf.Material(f"m{k}", kappa=MR.KAPPA[FILLER],
                                   R_s_contact=t["R_s_contact"] * knob,
                                   R_s_matrix=t["R_s_matrix"] * knob)
                   for k in (1, 2)},
        k_matrix=MR.K_MATRIX[MATRIX])
    cfg = Q.Q24Config(eta_cut=1.0, slab_axis=2, slab_frac=0.06)
    ids, typ, pos, dia, mol, info = EN.make(fam, D_UM, seed)
    G, Gt, sh, acfg = Q._net(ids, typ, pos, dia, mol, tab, cfg, D_UM)
    return dict(G=G, Gt=Gt, sh=sh, ids=np.asarray(ids), pos=np.asarray(pos),
                dia=np.asarray(dia), info=info)


def geff(fam, seed, knob):
    return float(tsf.effective_conductance(build(fam, seed, knob)["sh"], 0.0))


KNOB0, KNOB_MAX = 1.0, 1.0e4


def match(fam_a, fam_b, seed, rs0=None):
    rs0 = KNOB0 if rs0 is None else rs0
    ga, gb = geff(fam_a, seed, rs0), geff(fam_b, seed, rs0)
    if ga <= 0 or gb <= 0:
        return None
    hi_fam, target = (fam_a, gb) if ga >= gb else (fam_b, ga)
    lo, hi = rs0, KNOB_MAX
    f_lo, f_hi = geff(hi_fam, seed, lo), geff(hi_fam, seed, hi)
    if not (min(f_lo, f_hi) <= target <= max(f_lo, f_hi)):
        return None
    mid = lo
    for _ in range(45):
        mid = np.sqrt(lo * hi)
        f = geff(hi_fam, seed, mid)
        if abs(f / target - 1) < TOL:
            break
        if f > target:
            lo = mid
        else:
            hi = mid
    else:
        return None
    return ((mid, rs0) if hi_fam == fam_a else (rs0, mid)), hi_fam


def axes(st):
    G, Gt, sh = st["G"], st["Gt"], st["sh"]
    pos, dia, ids = st["pos"], st["dia"], st["ids"]
    idx = {int(v): k for k, v in enumerate(ids)}
    G0 = float(tsf.effective_conductance(sh, 0.0))
    flux = tsf.edge_flux_table(sh, 0.0)

    real, w_h, w_D, w_C, term = [], [], [], [], []
    for u, v, d in Gt.edges(data=True):
        if u < 0 or v < 0:
            term.append((u, v)); continue
        i, j = idx[int(u)], idx[int(v)]
        D = float(np.linalg.norm(pos[i] - pos[j]))
        real.append((int(u), int(v)))
        w_h.append(max(D - 0.5 * (dia[i] + dia[j]), 0.0))
        w_D.append(D)
        w_C.append(-float(np.log10(d["g"])))
    out = dict(G_eff=G0, n_edge=len(real))
    nodes = list(Gt.nodes())
    for tag, w in (("G1", np.array(w_h)), ("G2", np.array(w_D)),
                   ("C", np.array(w_C))):
        wmin = float(min(w.min(), 0.0)) - 1.0
        edges = [(u, v, float(x)) for (u, v), x in zip(real, w)]
        edges += [(u, v, wmin) for u, v in term]
        r = filtration(nodes, edges, tsf.SINK_NODE, tsf.SOURCE_NODE)
        out[f"tstar_{tag}"] = r["theta_star"]

    rng = np.random.default_rng(20260813)

    def boost(pairs):
        saved = []
        for u, v in pairs:
            saved.append((u, v, Gt[u][v]["g"])); Gt[u][v]["g"] *= BOOST
        g = tsf.effective_conductance(sh, 0.0)
        for u, v, x in saved:
            Gt[u][v]["g"] = x
        return g / G0 - 1.0

    k = min(K_TOP, len(real))
    gt = boost([(int(x.u), int(x.v)) for x in flux.head(k).itertuples()])
    gr = float(np.mean([boost([real[i] for i in
                               rng.choice(len(real), k, replace=False)])
                        for _ in range(N_RANDOM_CONTROL)]))
    out["gain_ratio"] = gt / gr if gr > 0 else np.nan
    return out


def main() -> None:
    want = {int(a) for a in sys.argv[1:]} or {1, 2, 3}
    rows, done = [], set()
    if CSV.is_file():
        old = pd.read_csv(CSV)
        rows = old.to_dict("records")
        done = set(zip(old["label"], old["seed"].astype(int)))
        print(f"(resume) {len(done)} conditions already done")
    t0 = time.time()
    for kind, label, fa, fb in PAIRS:
        if kind not in want:
            continue
        for seed in SEEDS:
            if (label, seed) in done:
                continue
            m = match(fa, fb, seed)
            if m is None:
                rows.append(dict(kind=kind, label=label, fam_a=fa, fam_b=fb,
                                 seed=seed, matched=False))
                print(f"  [skip] {label} seed{seed}: G_eff cannot be matched")
                pd.DataFrame(rows).assign(**tsf.engine_stamp()).to_csv(CSV, index=False)
                continue
            (rs_a, rs_b), tuned = m
            sa, sb = build(fa, seed, rs_a), build(fb, seed, rs_b)
            A, B = axes(sa), axes(sb)
            r = dict(kind=kind, label=label, fam_a=fa, fam_b=fb, seed=seed,
                     matched=True, tuned=tuned,
                     knob_a=float(rs_a), knob_b=float(rs_b),
                     knob_ratio=float(max(rs_a, rs_b) / min(rs_a, rs_b)),
                     phi_a=float(sa["info"]["phi"]), phi_b=float(sb["info"]["phi"]),
                     G_eff_a=A["G_eff"], G_eff_b=B["G_eff"],
                     geff_ratio=B["G_eff"] / A["G_eff"])
            for k_ in ("tstar_C", "tstar_G2", "tstar_G1", "gain_ratio"):
                r[f"{k_}_a"], r[f"{k_}_b"] = A[k_], B[k_]
                r[f"{k_}_ratio"] = (B[k_] / A[k_]) if A[k_] else np.nan
            rows.append(r)
            pd.DataFrame(rows).assign(**tsf.engine_stamp()).to_csv(CSV, index=False)
            print(f"  kind{kind} {label:<20} seed{seed}  G_eff ratio={r['geff_ratio']:.4f}  "
                  f"theta*={r['tstar_C_ratio']:.3f}  D*={r['tstar_G2_ratio']:.3f}  "
                  f"eps*={r['tstar_G1_ratio']:.3f}  lift={r['gain_ratio_ratio']:.3f}  "
                  f"({time.time()-t0:.0f}s)")
    print("== requested kinds done ==")


def report() -> None:
    df = pd.read_csv(CSV)
    ok = df[df.matched == True]
    log = []

    def say(s=""):
        print(s); log.append(s)

    def sep(x, lo=0.8, hi=1.25):
        return "separated" if (x < lo or x > hi) else ("not separated" if 0.9 <= x <= 1.1 else "undecided")

    say("=" * 104)
    say("Y4: separability of theta* (12 or more pairs)")
    say("=" * 104)
    say(f"pairs {df.label.nunique()} / conditions matched to equal G_eff {len(ok)} / {len(df)}"
        f" (yield {len(ok)/len(df):.0%}) / engine={df['engine_version'].iloc[0]}")
    say("verdict: ratio below 0.8 or above 1.25 = **separated** / 0.9-1.1 = **not separated** / in between = undecided")
    say()
    say(f"{'k':>2} {'pair':<22} {'n':>2} {'Geff rat':>8} {'theta*(C)':>8} {'':<5} "
        f"{'D*(G2)':>8} {'':<5} {'eps*(G1)':>8} {'':<5} {'lift':>8} {'':<5}")
    counts = {c: 0 for c in ("tstar_C", "tstar_G2", "tstar_G1", "gain_ratio")}
    for (kind, label), g in ok.groupby(["kind", "label"], sort=False):
        vals = {c: g[f"{c}_ratio"].median() for c in counts}
        for c, v in vals.items():
            if sep(v) == "separated":
                counts[c] += 1
        say(f"{kind:>2} {label:<22} {len(g):>2} {g.geff_ratio.median():8.4f} "
            + " ".join(f"{vals[c]:8.3f} {sep(vals[c]):<5}" for c in
                       ("tstar_C", "tstar_G2", "tstar_G1", "gain_ratio")))
    n_pair = ok.label.nunique()
    say()
    say(f"--- number of separated pairs (out of all {n_pair}) ---")
    for c, name in (("tstar_C", "theta*(perc.)"), ("tstar_G2", "D*(centre dist.)"),
                    ("tstar_G1", "eps*(surf. gap)"), ("gain_ratio", "lift")):
        say(f"  {name:<16} **{counts[c]} / {n_pair}**")
    say()
    for kind in sorted(ok.kind.unique()):
        sub = ok[ok.kind == kind]
        n = sub.label.nunique()
        cs = sum(1 for _, g in sub.groupby("label")
                 if sep(g["tstar_C_ratio"].median()) == "separated")
        say(f"  kind{kind}: pairs separated by theta* {cs} / {n}")

    say()
    say("=" * 104)
    say(f"* {ok[ok.label != DUP_LABEL].label.nunique()} independent pairs after removing the duplicate"
        f" (kind 1 '{DUP_LABEL}' is the same family pair as kind 3 'phi0.64 vs bimodal0.70')")
    say("=" * 104)
    ok15 = ok[ok.label != DUP_LABEL]
    c15 = {c: 0 for c in counts}
    tc = []
    for (kind, label), g in ok15.groupby(["kind", "label"], sort=False):
        vals = {c: g[f"{c}_ratio"].median() for c in c15}
        tc.append(vals["tstar_C"])
        for c, v in vals.items():
            if sep(v) == "separated":
                c15[c] += 1
    n15 = ok15.label.nunique()
    say(f"--- number of separated pairs (out of {n15} independent) ---")
    for c, name in (("tstar_C", "theta*(perc.)"), ("tstar_G2", "D*(centre dist.)"),
                    ("tstar_G1", "eps*(surf. gap)"), ("gain_ratio", "lift")):
        say(f"  {name:<16} **{c15[c]} / {n15}**")
    say()
    say("--- range of the theta*(perc.) ratio (* the value depends on the aggregation unit; the text must state the unit) ---")
    say(f"  per **pair** (seed median, {n15} pairs)     : "
        f"**{min(tc):.4f} - {max(tc):.4f}**   <- this is what the text quotes")
    say(f"  per condition (per seed, {len(ok15)} rows)   : "
        f"{ok15.tstar_C_ratio.min():.4f} - {ok15.tstar_C_ratio.max():.4f}")
    say()
    kmin = ok15.loc[(ok15.knob_ratio - 1).abs().idxmin()]
    g_kmax = ok15.groupby("label")["knob_ratio"].median().idxmax()
    kmax = ok15[ok15.label == g_kmax]
    say("--- interfacial thermal resistance multiplier (* not an input: result of the bisection for equal G_eff) ---")
    say(f"  pair moved least: {kmin.label:<22} x{kmin.knob_ratio:6.3f} -> "
        f"theta* ratio {ok15[ok15.label==kmin.label].tstar_C_ratio.median():.4f}")
    say(f"  pair moved most:  {g_kmax:<22} x{kmax.knob_ratio.median():6.1f} -> "
        f"theta* ratio {kmax.tstar_C_ratio.median():.4f}")
    TXT.write_text("\n".join(log) + "\n", encoding="utf-8")


if __name__ == "__main__":
    if "--report" in sys.argv:
        sys.argv.remove("--report"); report()
    else:
        main()
        report()
