# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/PoC3_stability/poc3_stability.py
#   sha256(src) : cbfe2e22a25eb6df
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 32 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# Perturbation stability of the barcode: d_B between the diagrams before and after
# a perturbation of the edge conductances, against the bound |Delta theta|_inf
# (Sec. 3.5): V19 (edge-weight perturbation). Also imported by poc3c_fixed_complex.py.
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
import time

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ensemble as EN
import materials_real as MR
import thermal_sheaf_filtration as tsf

FAMILIES = ["phi_060", "phi_064", "voided", "clustered", "stratified"]
SEEDS = (0, 1)
MATS = [("AlN", "epoxy", "untreated"), ("SiO2", "silicone", "untreated")]
JIT = [1e-4, 1e-3, 1e-2]
EPS = [0.01, 0.1, 1.0]
D_UM = 20.0

OUT = []
def say(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    OUT.append(line)


from persistence_distance import (
    TAU, K_CAP, _trim, _matchable, bottleneck)

def diagram(ids, typ, pos, dia, mol, tab, cfg, g_scale=None, rng=None):
    G = tsf.build_thermal_network(ids, typ, pos, dia, cfg, table=tab, mols=mol)
    src, snk = tsf.identify_slabs(ids, pos, cfg)
    if not src or not snk:
        return None, None
    Gt, _ = tsf.attach_virtual_terminals(G, src, snk, cfg)
    logs = None
    if g_scale is not None:
        logs = []
        for u, v, d in Gt.edges(data=True):
            f = 10.0 ** rng.uniform(-g_scale, g_scale)
            logs.append(abs(np.log10(f)))
            d["g"] = d["g"] * f
    hist = tsf.ConductanceFiltration(Gt).sweep()
    h = hist.dropna(subset=["sink_theta"])
    D = [(np.log10(max(b, 1e-300)), np.log10(max(s, 1e-300)))
         for b, s in zip(h["birth_theta"], h["sink_theta"])]
    return D, (max(logs) if logs else 0.0)


def main():
    t0 = time.time()
    cfg = tsf.AnalysisConfig(gap_max_um=0.5 * D_UM, eta_cut=1.0, slab_axis=2,
                             slab_frac=0.06, viz=False)
    rows = []
    for fam, seed, (f, m, t) in itertools.product(FAMILIES, SEEDS, MATS):
        ids, typ, pos, dia, mol, info = EN.make(fam, D_UM, seed)
        tab = MR.pair(f, t, "AlN", t, m)
        D0, _ = diagram(ids, typ, pos.copy(), dia, mol, tab, cfg)
        if D0 is None:
            continue
        base = dict(family=fam, seed=seed, filler=f, matrix=m, n_bars=len(D0))
        for e in EPS:
            rng = np.random.default_rng(100 + seed)
            D1, dmax = diagram(ids, typ, pos.copy(), dia, mol, tab, cfg,
                               g_scale=e, rng=rng)
            if D1 is None:
                continue
            db = bottleneck(D0, D1)
            rows.append(dict(kind="P-edge", level=e, d_B=db, bound=dmax,
                             ratio=(db / dmax if dmax > 0 else np.nan), **base))
        for j in JIT:
            rng = np.random.default_rng(200 + seed)
            p2 = pos + rng.normal(0.0, j * D_UM, size=pos.shape)
            D1, _ = diagram(ids, typ, p2, dia, mol, tab, cfg)
            if D1 is None:
                continue
            rows.append(dict(kind="P-geo", level=j, d_B=bottleneck(D0, D1),
                             bound=np.nan, ratio=np.nan, **base))
        for lv in ("low", "high"):
            tab2 = MR.pair(f, t, "AlN", t, m)
            for k in tab2.materials:
                mm = tab2.materials[k]
                sc = 10.0 ** (0.65 if lv == "high" else -0.65)
                mm.R_s_contact = mm.R_s_contact * sc
                mm.R_s_matrix = mm.R_s_matrix * sc
            D1, _ = diagram(ids, typ, pos.copy(), dia, mol, tab2, cfg)
            if D1 is None:
                continue
            rows.append(dict(kind="P-mat", level=(0.65 if lv == "high"
                                                  else -0.65),
                             d_B=bottleneck(D0, D1), bound=np.nan,
                             ratio=np.nan, **base))
    df = pd.DataFrame(rows)
    df.to_csv("results_poc3.csv", index=False)

    say("=" * 92)
    say("PoC-3: perturbation stability (bottleneck distance, own implementation)")
    say("=" * 92)
    say(f"conditions **{len(df)}** (families {df['family'].nunique()} x seed "
        f"{df['seed'].nunique()} x materials {df['filler'].nunique()} x 3 perturbation types)")
    say(f"bars per persistence diagram {df['n_bars'].min()}-{df['n_bars'].max()}"
        f" (bottleneck computed only over bars with persistence >= {TAU} decades)")
    say(f"run time {time.time()-t0:.0f}s")
    say()
    say("--- S1: does 1-Lipschitz (d_B <= ||Delta log g||_inf) hold under P-edge? ---")
    pe = df[df.kind == "P-edge"]
    ok = (pe["d_B"] <= pe["bound"] + 2 * TAU + 1e-12)
    ok_strict = (pe["d_B"] <= pe["bound"] * (1 + 1e-9))
    say(f"  threshold TAU={TAU} decades (truncation-error tolerance 2*TAU={2*TAU})")
    say(f"  d_B <= ||Delta log g||_inf + 2*TAU : holds {int(ok.sum())} / {len(pe)}")
    say(f"  d_B <= ||Delta log g||_inf (strict)   : holds {int(ok_strict.sum())} / {len(pe)}")
    say(f"  excess max(d_B - bound) = "
        f"{(pe['d_B'] - pe['bound']).max():+.5f} decades")
    say(f"  verdict: {'PASS holds within truncation error' if ok.all() else 'FAIL violations present'}")
    say()
    say("--- S2: ratio to the bound d_B/||Delta log g||_inf ---")
    say(f"  {'eps [dec]':>12} {'n':>4} {'d_B med':>10} {'bound med':>10} "
        f"{'ratio med':>9} {'ratio max':>9}")
    for lv, g in pe.groupby("level"):
        say(f"  {lv:12.2f} {len(g):4d} {g['d_B'].median():10.4f} "
            f"{g['bound'].median():10.4f} {g['ratio'].median():9.4f} "
            f"{g['ratio'].max():9.4f}")
    say(f"  overall median ratio **{pe['ratio'].median():.4f}**"
        f" (within the 0.1-1.0 band: "
        f"{'PASS' if 0.1 <= pe['ratio'].median() <= 1.0 else 'FAIL'})")
    say()
    say("--- S3: response to P-geo (position jitter delta/d) ---")
    pg = df[df.kind == "P-geo"]
    say(f"  {'delta/d':>10} {'n':>4} {'d_B med':>10} {'d_B max':>10}")
    for lv, g in pg.groupby("level"):
        say(f"  {lv:10.0e} {len(g):4d} {g['d_B'].median():10.4f} "
            f"{g['d_B'].max():10.4f}")
    mm = pg.groupby("level")["d_B"].median()
    if len(mm) >= 2 and (mm > 0).all():
        sl = np.polyfit(np.log10(mm.index.values), np.log10(mm.values), 1)[0]
        say(f"  log-log slope **{sl:.3f}** (within the 0.5-1.5 band: "
            f"{'PASS' if 0.5 <= sl <= 1.5 else 'FAIL'})")
    say()
    say("--- S4: P-mat (R_c varied over the literature range +-0.65 decades) ---")
    pm = df[df.kind == "P-mat"]
    say(f"  d_B median **{pm['d_B'].median():.4f} decades** / max "
        f"{pm['d_B'].max():.4f} decades")
    say(f"  below 0.3 decades: {'PASS' if pm['d_B'].max() < 0.3 else 'FAIL'}")
    say()
    say("--- S5: differences between structure families ---")
    say(f"  {'family':>12} {'P-edge(eps=.1)':>14} {'P-geo(1e-3)':>12} "
        f"{'P-mat':>9}")
    for fam, g in df.groupby("family"):
        a = g[(g.kind == "P-edge") & (g.level == 0.1)]["d_B"].median()
        b = g[(g.kind == "P-geo") & (g.level == 1e-3)]["d_B"].median()
        c = g[g.kind == "P-mat"]["d_B"].median()
        say(f"  {fam:>12} {a:14.4f} {b:12.4f} {c:9.4f}")
    sp = df[df.kind == "P-mat"].groupby("family")["d_B"].median()
    say(f"  spread of P-mat across families x{sp.max()/max(sp.min(),1e-12):.2f}"
        f" (within 3x: {'PASS' if sp.max()/max(sp.min(),1e-12) <= 3 else 'FAIL'})")
    say()

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))
    a = ax[0]
    a.scatter(pe["bound"], pe["d_B"], s=18, alpha=.7)
    lim = [1e-3, max(pe["bound"].max(), 1) * 1.3]
    a.plot(lim, lim, "r--", lw=1, label=r"$d_B=\|\Delta\log g\|_\infty$")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel(r"$\|\Delta\log_{10} g\|_\infty$ [decades]")
    a.set_ylabel(r"$d_B$ [decades]")
    a.set_title("(a) S1: 1-Lipschitz stability")
    a.legend(fontsize=7); a.grid(alpha=.3, which="both")
    a = ax[1]
    for lv, g in pg.groupby("level"):
        a.scatter([lv] * len(g), g["d_B"], s=16, alpha=.6)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel(r"jitter $\delta/d$"); a.set_ylabel(r"$d_B$ [decades]")
    a.set_title("(b) S3: geometric perturbation")
    a.grid(alpha=.3, which="both")
    a = ax[2]
    dat = [df[(df.kind == k)]["d_B"].values
           for k in ("P-edge", "P-geo", "P-mat")]
    a.boxplot(dat, tick_labels=["P-edge", "P-geo", "P-mat"])
    a.set_yscale("log"); a.set_ylabel(r"$d_B$ [decades]")
    a.set_title("(c) perturbation type comparison")
    a.grid(alpha=.3, axis="y", which="both")
    fig.tight_layout(pad=1.3); fig.savefig("fig_poc3.png", dpi=140)
    say("figure saved: fig_poc3.png")
    pathlib.Path("results_poc3.txt").write_text("\n".join(OUT) + "\n",
                                                encoding="utf-8")


if __name__ == "__main__":
    main()
