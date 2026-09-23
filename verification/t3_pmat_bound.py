# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/PoC3_stability/t3_pmat_bound.py
#   sha256(src) : ef2e3b8922eb4ee2
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 29 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V20. Propagation of the uncertainty in the filler-filler interfacial resistance R_c
# (literature range +-0.65 decades, i.e. 1.30 decades in total) to the edge conductances
# and to the barcode (Sec. 7.3, Table 9): the bound |Delta log g|_inf, the bottleneck
# distance d_B between the barcodes, and their ratio, for 5 families x 2 seeds x 6 fillers
# x 2 directions = 120 conditions. Writes results_t3.csv / results_t3.txt.
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
import numpy as np
import pandas as pd

import ensemble as EN
import materials_real as MR
import thermal_sheaf_filtration as tsf
from poc3_stability import bottleneck, TAU

D_UM = 20.0
FAMILIES = ("phi_060", "phi_064", "clustered", "stratified", "voided")
SEEDS = (0, 1)
FILLERS = ("AlN", "Al2O3", "SiO2", "diamond", "hBN", "ZnO")
MATRIX = "epoxy"
TREAT = "untreated"
OUT = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)


def edges_g(ids, typ, pos, dia, mol, tab, cfg):
    G = tsf.build_thermal_network(ids, typ, pos, dia, cfg, table=tab, mols=mol)
    src, snk = tsf.identify_slabs(ids, pos, cfg)
    if not src or not snk:
        return None, None
    Gt, _ = tsf.attach_virtual_terminals(G, src, snk, cfg)
    gs = {frozenset((u, v)): d["g"] for u, v, d in Gt.edges(data=True)
          if u >= 0 and v >= 0}
    hist = tsf.ConductanceFiltration(Gt, cfg).sweep()
    h = hist.dropna(subset=["sink_theta"])
    D = [(np.log10(max(b, 1e-300)), np.log10(max(s, 1e-300)))
         for b, s in zip(h["birth_theta"], h["sink_theta"])]
    return gs, D


def main():
    cfg = tsf.AnalysisConfig(gap_max_um=0.5 * D_UM, eta_cut=1.0, slab_axis=2,
                             slab_frac=0.06, viz=False)
    rows = []
    for fam, seed, fil in itertools.product(FAMILIES, SEEDS, FILLERS):
        ids, typ, pos, dia, mol, info = EN.make(fam, D_UM, seed)
        tab0 = MR.pair(fil, TREAT, "AlN", TREAT, MATRIX)
        g0, D0 = edges_g(ids, typ, pos.copy(), dia, mol, tab0, cfg)
        if g0 is None:
            continue
        for lv, sc in (("low", 10.0 ** -0.65), ("high", 10.0 ** 0.65)):
            tab1 = MR.pair(fil, TREAT, "AlN", TREAT, MATRIX)
            for k in tab1.materials:
                mm = tab1.materials[k]
                mm.R_s_contact = mm.R_s_contact * sc
                mm.R_s_matrix = mm.R_s_matrix * sc
            g1, D1 = edges_g(ids, typ, pos.copy(), dia, mol, tab1, cfg)
            if g1 is None:
                continue
            common = set(g0) & set(g1)
            dlog = np.array([abs(np.log10(g1[k] / g0[k])) for k in common])
            bound = float(dlog.max())
            db = bottleneck(D0, D1)
            rows.append(dict(family=fam, seed=seed, filler=fil, level=lv,
                             n_edges=len(common),
                             bound=bound, dlog_med=float(np.median(dlog)),
                             d_B=db, ratio=(db / bound if bound > 0 else np.nan),
                             ok=int(db <= bound + 2 * TAU + 1e-12),
                             ok_strict=int(db <= bound * (1 + 1e-9))))
        say(f"  {fam:>12} seed={seed} {fil:>8}  bound={rows[-1]['bound']:.4f} "
            f"d_B={rows[-1]['d_B']:.4f} ratio={rows[-1]['ratio']:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(pathlib.Path(__file__).with_name("results_t3.csv"), index=False)

    say("")
    say("=" * 92)
    say("T3: upper bound |Delta log g|_inf under a perturbation of R_c over the literature range (1.30 decades)")
    say("=" * 92)
    say(f"conditions {len(df)} (families {len(FAMILIES)} x seeds {len(SEEDS)} x "
        f"fillers {len(FILLERS)} x 2 directions) / edges {df.n_edges.min()}-{df.n_edges.max()}")
    say("")
    say("--- all conditions ---")
    say(f"  |Delta log g|_inf (bound)   median {df.bound.median():.4f} decades / "
        f"5–95% {df.bound.quantile(.05):.4f}–{df.bound.quantile(.95):.4f} / "
        f"max {df.bound.max():.4f}")
    say(f"  median |Delta log g|         median {df.dlog_med.median():.4f} decades"
        f" (= how much a typical edge moves)")
    say(f"  d_B                          median {df.d_B.median():.4f} decades / max {df.d_B.max():.4f}")
    say(f"  ratio d_B / bound            median {df.ratio.median():.4f} / "
        f"5–95% {df.ratio.quantile(.05):.4f}–{df.ratio.quantile(.95):.4f}")
    say(f"  bound holds                  strictly {df.ok_strict.sum()}/{len(df)} / "
        f"with tolerance {df.ok.sum()}/{len(df)}")
    say(f"  compression 1.30 decades / bound  median **{1.30/df.bound.median():.1f} x**")
    say("")
    say("--- by filler (in order of alpha = kappa/k_m) ---")
    say(f"{'filler':>10s} {'alpha':>8s} {'bound med':>10s} {'d_B med':>10s} "
        f"{'ratio':>7s} {'compr.':>8s}")
    for fil in sorted(FILLERS, key=lambda f: MR.KAPPA[f]):
        sub = df[df.filler == fil]
        b = sub.bound.median()
        say(f"{fil:>10s} {MR.KAPPA[fil]/MR.K_MATRIX[MATRIX]:8.0f} {b:10.4f} "
            f"{sub.d_B.median():10.4f} {sub.ratio.median():7.3f} "
            f"{1.30/b if b>0 else np.inf:8.1f}")
    say("")
    say("--- comparison with the predictions ---")
    b = df.bound.median()
    say(f"T3-1 |Delta log g|_inf in 0.05-0.12 decades: measured median {b:.4f} -> "
        f"{'a-priori prediction: met' if 0.05 <= b <= 0.12 else 'a-priori prediction: not met'}")
    _v = ("a-priori prediction: met" if df.ok_strict.all()
          else f"WARNING: exceeded in part strictly (with tolerance {df.ok.sum()}/{len(df)})")
    say(f"T3-2 bound holds in all conditions: strictly {df.ok_strict.sum()}/{len(df)} -> {_v}")
    r = df.ratio.median()
    say(f"T3-3 median ratio in 0.7-1.0: measured {r:.4f} -> "
        f"{'a-priori prediction: met' if 0.7 <= r <= 1.0 else 'a-priori prediction: not met'}")
    bb = df.groupby("filler").bound.median().sort_values()
    say(f"T3-4 compression depends on the filler: bound min {bb.index[0]} {bb.iloc[0]:.4f} -> "
        f"max {bb.index[-1]} {bb.iloc[-1]:.4f} (x{bb.iloc[-1]/bb.iloc[0]:.2f})"
        f" -> {'a-priori prediction: met' if bb.iloc[-1]/bb.iloc[0] > 1.5 else 'a-priori prediction: not met'}")
    say(f"    bound for diamond {df[df.filler=='diamond'].bound.median():.4f} vs "
        f"SiO₂ {df[df.filler=='SiO2'].bound.median():.4f}")
    pathlib.Path(__file__).with_name("results_t3.txt").write_text(
        "\n".join(OUT) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
