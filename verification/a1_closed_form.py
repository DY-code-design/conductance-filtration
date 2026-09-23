# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/A1_closed_form_ensemble/a1_closed_form.py
#   sha256(src) : 6197f60ef6aac309
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 33 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V11. Over the 69 structures of the ensemble: maximal relative error of the
# finite-change closed form (Sec. 4.4) against a full re-solve, the error made when
# the ordinary effective resistance is substituted for the shorted R~_e, and the lift
# of the top-k edges over k random edges.
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

import ensemble as EN
import thermal_sheaf_filtration as tsf

SEEDS = (0, 1, 2)
K_VERIFY = 20
DELTA_OVER_G = (1.0, 10.0, 100.0)
K_LIST = (1, 5, 10, 20, 50)
BOOST = 10.0
N_RANDOM_CONTROL = _bootstrap.N_RANDOM_CONTROL

CSV = _pathlib.Path("results_a1_closed_form.csv")
CSV_K = _pathlib.Path("results_a1_gain_k.csv")
TXT = _pathlib.Path("results_a1_closed_form.txt")


def one(st, rng):
    Gt, sh, G = st["Gt"], st["sheaf"], st["G"]
    G0 = float(tsf.effective_conductance(sh, 0.0))
    flux = tsf.edge_flux_table(sh, 0.0)
    real = [(int(r.u), int(r.v)) for r in flux.itertuples()]
    dT = {frozenset(k): float(r.dT) for k, r in zip(real, flux.itertuples())}

    top = real[:K_VERIFY]
    rnd = [real[i] for i in rng.choice(len(real), min(K_VERIFY, len(real)),
                                       replace=False)]
    tgt = top + rnd
    rt = tsf.shorted_edge_resistance(sh, 0.0, edges=tgt, method="sparse")
    RT = {frozenset((int(r.u), int(r.v))): float(r.R_tilde) for r in rt.itertuples()}
    red, _ = tsf.edge_redundancy(sh, 0.0, max_n=10 ** 9)
    RE = {frozenset((int(r.u), int(r.v))): float(r.R_eff_uv) for r in red.itertuples()}

    errs, errs_w, ratios = [], [], []
    for u, v in tgt:
        k = frozenset((u, v))
        if k not in RT or k not in dT:
            continue
        for c in DELTA_OVER_G:
            D = c * Gt[u][v]["g"]
            pred = tsf.predict_finite_change(G0, dT[k], RT[k], D)
            sv = Gt[u][v]["g"]; Gt[u][v]["g"] = sv + D
            true = float(tsf.effective_conductance(sh, 0.0)); Gt[u][v]["g"] = sv
            errs.append(abs(pred / true - 1))
            if k in RE:
                pw = tsf.predict_finite_change(G0, dT[k], RE[k], D)
                errs_w.append(abs(pw / true - 1))
        if k in RE:
            ratios.append(RE[k] / RT[k])

    def boost(pairs):
        saved = []
        for u, v in pairs:
            saved.append((u, v, Gt[u][v]["g"])); Gt[u][v]["g"] *= BOOST
        g = float(tsf.effective_conductance(sh, 0.0))
        for u, v, x in saved:
            Gt[u][v]["g"] = x
        return g / G0 - 1.0

    gains = {}
    for k in K_LIST:
        kk = min(k, len(real))
        gt = boost(real[:kk])
        gr = float(np.mean([boost([real[i] for i in
                                   rng.choice(len(real), kk, replace=False)])
                            for _ in range(N_RANDOM_CONTROL)]))
        gains[k] = dict(gain_top_pct=gt * 100, gain_rand_pct=gr * 100,
                        gain_ratio=(gt / gr if gr > 0 else np.nan))

    info = st["info"]
    base = dict(family=info.get("family"), seed=info.get("seed"),
                phi=float(info.get("phi", np.nan)),
                N=G.number_of_nodes(), E=G.number_of_edges(), G_eff=G0)
    row = dict(base,
               n_edge_tested=len(errs) // len(DELTA_OVER_G),
               err_closed_max=max(errs), err_closed_med=float(np.median(errs)),
               err_wrong_max=(max(errs_w) if errs_w else np.nan),
               err_wrong_med=(float(np.median(errs_w)) if errs_w else np.nan),
               ratio_med=(float(np.median(ratios)) if ratios else np.nan),
               ratio_max=(max(ratios) if ratios else np.nan),
               ratio_min=(min(ratios) if ratios else np.nan),
               frac_ratio_lt_1p001=(float(np.mean([r < 1.001 for r in ratios]))
                                    if ratios else np.nan))
    krows = [dict(base, k=k, **v) for k, v in gains.items()]
    return row, krows


def main() -> None:
    rows, krows, done = [], [], set()
    if CSV.is_file():
        old = pd.read_csv(CSV); rows = old.to_dict("records")
        done = set(zip(old.family, old.seed.astype(int)))
        krows = pd.read_csv(CSV_K).to_dict("records") if CSV_K.is_file() else []
        print(f"(resume) done {len(done)} structures")
    t0 = time.time()
    remaining = [f for f in EN.FAMILY_ORDER
                 if any((f, s) not in done for s in SEEDS)]
    if len(remaining) < len(EN.FAMILY_ORDER):
        print(f"(resume) families to generate {len(remaining)}/{len(EN.FAMILY_ORDER)}: "
              f"{', '.join(remaining) if remaining else 'none'}")
    for st in EN.ensemble(families=remaining, seeds=SEEDS):
        key = (st["info"].get("family"), int(st["info"].get("seed")))
        if key in done:
            continue
        r, kr = one(st, np.random.default_rng(1000 + key[1]))
        rows.append(r); krows += kr
        pd.DataFrame(rows).assign(**tsf.engine_stamp()).to_csv(CSV, index=False)
        pd.DataFrame(krows).assign(**tsf.engine_stamp()).to_csv(CSV_K, index=False)
        print(f"  {r['family']:<14} seed{r['seed']}  closed form {r['err_closed_max']:.2e}  "
              f"wrong R {r['err_wrong_max']:.2e}  ratio {r['ratio_med']:.6f}  "
              f"({time.time()-t0:.0f}s)")
    print("== done ==")


def report() -> None:
    d = pd.read_csv(CSV); k = pd.read_csv(CSV_K)
    log = []

    def say(s=""):
        print(s); log.append(s)

    say("=" * 96)
    say("A1: population check of the Sec. 4.4 closed form (finite change) + k dependence of the lift")
    say("=" * 96)
    say(f"population {len(d)} structures (phi {d.phi.min():.3f}-{d.phi.max():.3f}) / "
        f"per structure: top {K_VERIFY} + random {K_VERIFY} edges x Delta/g={DELTA_OVER_G} / "
        f"engine={d['engine_version'].iloc[0]}")
    say()
    say("--- A-1: relative error of the closed form (shorted R~) ---")
    say(f"  per-structure maximum: min {d.err_closed_max.min():.3e} / "
        f"median {d.err_closed_max.median():.3e} / **max {d.err_closed_max.max():.3e}**")
    say(f"  per-structure median: median {d.err_closed_med.median():.3e}")
    for th in (1e-10, 1e-11, 1e-12):
        say(f"    structures with max error < {th:.0e}: {int((d.err_closed_max < th).sum())} / {len(d)}")
    say()
    say("--- A-2: * when the ordinary R_eff is substituted by mistake ---")
    say(f"  per-structure maximum: min {d.err_wrong_max.min():.3e} / "
        f"median {d.err_wrong_max.median():.3e} / max {d.err_wrong_max.max():.3e}")
    say(f"  => gap to the closed form in decades (median to median): "
        f"{np.log10(d.err_wrong_max.median() / d.err_closed_max.median()):.1f} decades")
    say()
    say("  distribution of the ratio R_eff/R~ (the closer to 1, the less noticeable)")
    say(f"    median {d.ratio_med.median():.6f} / max {d.ratio_max.max():.6f} / "
        f"min {d.ratio_min.min():.6f}")
    say(f"    * fraction of edges with ratio < 1.001: median {d.frac_ratio_lt_1p001.median():.1%}"
        f"(per-structure min {d.frac_ratio_lt_1p001.min():.1%} / max {d.frac_ratio_lt_1p001.max():.1%})")
    say()
    say("--- A-3: k dependence of the lift (C=10, 69 structures) ---")
    say(f"{'k':>4} {'lift median':>12} {'min':>9} {'max':>10} "
        f"{'top-k gain[%]':>16} {'random-k gain[%]':>18}")
    for kk, g in k.groupby("k"):
        say(f"{kk:>4} {g.gain_ratio.median():12.2f} {g.gain_ratio.min():9.2f} "
            f"{g.gain_ratio.max():10.1f} {g.gain_top_pct.median():16.3f} "
            f"{g.gain_rand_pct.median():18.4f}")
    say()
    _a1 = k[k.k == 10].gain_ratio.median()
    try:
        _ens = _bootstrap.poc_csv("results_ensemble.csv", "..", "E_ensemble")
    except FileNotFoundError:
        _ens = None
    say("  WARNING: cross-check: median at k=10 against the median `lift` in `results_ensemble.csv`")
    if _ens is not None:
        _l = pd.read_csv(_ens)["lift"].median()
        say(f"     A1 k=10 = {_a1:.2f} / ensemble lift median = {_l:.2f} "
            f"/ ratio **{_a1/_l:.3f}**")
        say("     WARNING: the two are **not strictly equal** (A1 takes the top k in P^share order, "
            "ensemble picks K_TOP edges with a different random stream). **Matching order of magnitude and direction is enough.**")
    else:
        say(f"     A1 k=10 = {_a1:.2f} (no ensemble CSV, so no comparison)")
    TXT.write_text("\n".join(log) + "\n", encoding="utf-8")


if __name__ == "__main__":
    import sys
    if "--report" in sys.argv:
        report()
    else:
        main()
        report()
