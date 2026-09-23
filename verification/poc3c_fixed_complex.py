# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/PoC3_stability/poc3c_fixed_complex.py
#   sha256(src) : aed526ea0396e53d
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 19 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V19. Stability bound under coordinate perturbation: with the edge set frozen
# (fixed_pairs) the particles are jittered and d_B <= |Delta theta|_inf,
# theta = -log10 g, is checked (Sec. 3.5).
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

import numpy as np
import pandas as pd

import ensemble as EN
import materials_real as MR
import thermal_sheaf_filtration as tsf
from poc3_stability import bottleneck, TAU

FAMILIES = ["phi_060", "phi_064", "voided", "clustered", "stratified"]
SEEDS = (0, 1)
MATS = [("AlN", "epoxy", "untreated"), ("SiO2", "silicone", "untreated")]
JIT = [1e-4, 1e-3, 1e-2]
D_UM = 20.0

OUT = []
def say(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    OUT.append(line)


def theta_map(Gt):
    return {tuple(sorted((int(u), int(v)))): -np.log10(max(d["g"], 1e-300))
            for u, v, d in Gt.edges(data=True)}


def diagram(Gt, pairing):
    h = tsf.ConductanceFiltration(Gt).sweep(pairing=pairing)
    h = h.dropna(subset=["sink_theta"])
    return [(-np.log10(max(b, 1e-300)), -np.log10(max(s, 1e-300)))
            for b, s in zip(h["birth_theta"], h["sink_theta"])]


def main():
    t0 = time.time()
    rows = []
    for fam, seed, (f, m, t) in itertools.product(FAMILIES, SEEDS, MATS):
        ids, typ, pos, dia, mol, info = EN.make(fam, D_UM, seed)
        tab = MR.pair(f, t, "AlN", t, m)
        for fixed in (False, True):
            cfg = tsf.AnalysisConfig(gap_max_um=0.5 * D_UM, eta_cut=1.0,
                                     slab_axis=2, slab_frac=0.06, viz=False,
                                     fixed_complex=fixed,
                                     fixed_complex_gap_um=(0.7 * D_UM
                                                           if fixed else 0.0))
            fp = None
            if fixed:
                cfg_enum = tsf.AnalysisConfig(
                    gap_max_um=0.5 * D_UM, eta_cut=1.0,
                    fixed_complex=True, fixed_complex_gap_um=0.55 * D_UM)
                fp = tsf.pairs_from_positions(pos, dia, cfg_enum)

            def build(pp):
                G = tsf.build_thermal_network(ids, typ, pp, dia, cfg,
                                              table=tab, mols=mol,
                                              fixed_pairs=fp)
                src, snk = tsf.identify_slabs(ids, pp, cfg)
                if not src or not snk:
                    return None
                Gt, _ = tsf.attach_virtual_terminals(G, src, snk, cfg)
                return Gt

            Gt0 = build(pos.copy())
            if Gt0 is None:
                continue
            th0 = theta_map(Gt0)
            D0 = {pg: diagram(Gt0, pg) for pg in ("relative_h0", "legacy")}
            for j in JIT:
                rng = np.random.default_rng(300 + seed)
                Gt1 = build(pos + rng.normal(0.0, j * D_UM, size=pos.shape))
                if Gt1 is None:
                    continue
                th1 = theta_map(Gt1)
                common = set(th0) & set(th1)
                dsym = len(set(th0) ^ set(th1))
                bound = max((abs(th0[k] - th1[k]) for k in common), default=0.0)
                n_oor = sum(1 for _u, _v, _d in Gt1.edges(data=True)
                            if _d.get("out_of_range"))
                for pg in ("relative_h0", "legacy"):
                    D1 = diagram(Gt1, pg)
                    rows.append(dict(family=fam, seed=seed, filler=f,
                                     matrix=m, jitter=j, fixed=fixed,
                                     pairing=pg, n_edge=len(th0),
                                     d_sym=dsym, bound=bound, n_oor=n_oor,
                                     d_B=bottleneck(D0[pg], D1)))
    df = pd.DataFrame(rows)
    df.to_csv("results_poc3c.csv", index=False)
    df["ratio"] = df["d_B"] / df["bound"].replace(0, np.nan)
    df["ok"] = df["d_B"] <= df["bound"] + 2 * TAU + 1e-12
    df["ok_strict"] = df["d_B"] <= df["bound"] * (1 + 1e-9)

    say("=" * 96)
    say("PoC-3c: do premises (A) fixed complex + (C) elder rule make P-geo Lipschitz as well?")
    say("=" * 96)
    say(f"conditions **{len(df)}** / run {time.time()-t0:.0f}s")
    say()
    say("--- is the complex really fixed? (symmetric difference of the edge sets before/after perturbation) ---")
    for fx, g in df.groupby("fixed"):
        say(f"  fixed_complex={str(fx):>5}: edges {int(g['n_edge'].median())}, "
            f"symmetric difference median **{int(g['d_sym'].median())}** / max "
            f"{int(g['d_sym'].max())}, out-of-range edges (g_floor) max "
            f"**{int(g['n_oor'].max())}**")
    say()
    say("--- 1-Lipschitz verdict (d_B <= ||Delta theta||_inf) ---")
    say(f"  {'cplx':>6} {'pairing':>13} {'strict':>8} {'+2*TAU':>8} "
        f"{'excess max':>10} {'ratio med':>9} {'ratio max':>9}")
    for (fx, pg), g in df.groupby(["fixed", "pairing"]):
        say(f"  {str(fx):>6} {pg:>13} "
            f"{f'{int(g.ok_strict.sum())}/{len(g)}':>8} "
            f"{f'{int(g.ok.sum())}/{len(g)}':>8} "
            f"{(g['d_B']-g['bound']).max():+10.5f} {g['ratio'].median():9.4f} "
            f"{g['ratio'].max():9.4f}")
    say()
    say("--- by jitter (fixed complex, relative H_0) ---")
    sub = df[(df.fixed) & (df.pairing == "relative_h0")]
    say(f"  {'delta/d':>8} {'n':>4} {'|Dtheta| med':>12} {'d_B med':>10} {'ratio med':>9}")
    for j, g in sub.groupby("jitter"):
        say(f"  {j:8.0e} {len(g):4d} {g['bound'].median():12.5f} "
            f"{g['d_B'].median():10.5f} {g['ratio'].median():9.4f}")
    say()
    a = df[(df.fixed) & (df.pairing == "relative_h0")]
    b = df[(~df.fixed) & (df.pairing == "relative_h0")]
    say(f"**Conclusion**")
    say(f"  fixed complex + relative H_0 : 1-Lipschitz "
        f"{'PASS holds' if bool(a['ok'].all()) else 'FAIL violations present'}"
        f" (strict {int(a.ok_strict.sum())}/{len(a)})")
    say(f"  without fixing the complex: "
        f"{'PASS' if bool(b['ok'].all()) else 'FAIL violations present'}"
        f" (strict {int(b.ok_strict.sum())}/{len(b)}, symmetric difference median "
        f"{int(b['d_sym'].median())} edges enter or leave)")
    pathlib.Path("results_poc3c.txt").write_text("\n".join(OUT) + "\n",
                                                 encoding="utf-8")


if __name__ == "__main__":
    main()
