# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/V8V9_log_sensitivity/v89_log_sensitivity.py
#   sha256(src) : 463cc816e0481dc3
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# V8, V9. On the N = 900 packing (family phi_060, seed 0, AlN / epoxy / untreated,
# 3,514 edges): V8, the sum of the dissipation over all edges (terminal edges
# included) against G_eff; V9, the dissipation share P_e / G_eff of every edge
# against the central difference of log G_eff with respect to log g_e (+-1e-6).
# Resumable: partial results are saved to results_v89_log_sensitivity.csv.
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

CSV = _pathlib.Path("results_v89_log_sensitivity.csv")
TXT = _pathlib.Path("results_v89_log_sensitivity.txt")
H = 1e-6
SAVE_EVERY = 200
BANDS = [(1e-2, 1.0), (1e-3, 1e-2), (1e-4, 1e-3), (1e-5, 1e-4), (0.0, 1e-5)]


def say(f, s=""):
    print(s)
    f.write(s + "\n")


def main():
    st = EN.build("phi_060", d_um=20.0, seed=0, filler="AlN", treatment="untreated", matrix="epoxy")
    Gt, sheaf = st["Gt"], st["sheaf"]
    T, Q = tsf.solve_dirichlet(sheaf, 0.0)
    n_nodes = Gt.number_of_nodes() - 2
    phys = [(u, v, d) for u, v, d in Gt.edges(data=True) if u >= 0 and v >= 0]

    P_all = sum(d["g"] * (T[u] - T[v]) ** 2 for u, v, d in Gt.edges(data=True))
    P_phys = sum(d["g"] * (T[u] - T[v]) ** 2 for u, v, d in phys)
    rel_v8 = abs(P_all - Q) / Q

    COLS = ["u", "v", "g", "share_exact", "slope_cd", "rel_err"]
    rows, done = [], set()
    if CSV.is_file():
        old = pd.read_csv(CSV)
        rows = old[COLS].to_dict("records")
        done = {(int(r["u"]), int(r["v"])) for r in rows}
    todo = [(u, v, d) for u, v, d in phys if (int(u), int(v)) not in done]
    t0 = time.time()
    for k, (u, v, d) in enumerate(todo, 1):
        g0 = d["g"]
        share = g0 * (T[u] - T[v]) ** 2 / Q
        Gt[u][v]["g"] = g0 * np.exp(H)
        Gp = tsf.effective_conductance(sheaf, 0.0)
        Gt[u][v]["g"] = g0 * np.exp(-H)
        Gm = tsf.effective_conductance(sheaf, 0.0)
        Gt[u][v]["g"] = g0
        slope = (np.log(Gp) - np.log(Gm)) / (2 * H)
        rows.append(dict(u=int(u), v=int(v), g=g0, share_exact=share, slope_cd=slope,
                         rel_err=abs(slope - share) / share if share > 0 else np.nan))
        if k % SAVE_EVERY == 0 or k == len(todo):
            pd.DataFrame(rows).assign(**tsf.engine_stamp()).to_csv(CSV, index=False)
    df = pd.DataFrame(rows)
    if len(df) < len(phys):
        print(f"resumable: {len(df)}/{len(phys)} edges done; run again to continue")
        raise SystemExit(3)

    with TXT.open("w", encoding="utf-8") as f:
        say(f, f"structure: family phi_060, seed 0, AlN/epoxy/untreated, N={n_nodes}, "
               f"E(physical)={len(phys)}, variant {tsf.active_variant()}")
        say(f, f"G_eff = {Q:.12e} W/K")
        say(f, "")
        say(f, "=== V8: sum_e P_e = G_eff ===")
        say(f, f"  sum over all edges (terminal edges included) = {P_all:.12e}   rel. diff = {rel_v8:.3e}"
               f"   {'PASS' if rel_v8 < 1e-12 else 'FAIL'}")
        say(f, f"  sum over the physical edges only               = {P_phys:.12e}   ratio to G_eff = {P_phys / Q:.9f}")
        say(f, "")
        say(f, f"=== V9: P_e/G_eff vs central difference of log G_eff w.r.t. log g_e (+-{H:g}) ===")
        say(f, f"  {'share band':>18} {'edges':>6} {'max rel err':>12} {'median':>10}")
        for lo, hi in BANDS:
            m = (df.share_exact >= lo) & (df.share_exact < hi)
            if m.any():
                say(f, f"  [{lo:8.0e}, {hi:8.0e}) {int(m.sum()):6d} {df.rel_err[m].max():12.3e} {df.rel_err[m].median():10.3e}")
        top = df.loc[df.share_exact.idxmax()]
        say(f, f"  top edge ({int(top.u)},{int(top.v)}): share {top.share_exact:.4e}, rel err {top.rel_err:.3e}")
        m3 = df.share_exact >= 1e-3
        say(f, f"  verdict (share >= 1e-3, {int(m3.sum())} edges): max rel err {df.rel_err[m3].max():.3e}"
               f"   {'PASS' if df.rel_err[m3].max() < 1e-5 else 'FAIL'}")
        say(f, "  note: for small shares the rounding of the difference (~1e-15 / h) dominates the relative error.")
        say(f, f"  ({time.time() - t0:.0f} s for {len(todo)} edges in this run)")
    print(f"wrote {CSV} / {TXT}")


if __name__ == "__main__":
    main()
