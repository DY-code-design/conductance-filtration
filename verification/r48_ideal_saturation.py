# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/R48_ideal_saturation/r48_ideal_saturation.py
#   sha256(src) : c5fd078ce884a4af
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 13 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V14. Saturation of the finite surrogates: the idealised gains of Sec. 4.5 use a
# factor of 1e6 in place of infinity (1e-30 in place of zero); the factor is swept
# to confirm that the gains have converged.
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

import copy
import numpy as np
import pandas as pd

import ensemble as EN
import materials_real as MR
import quantities24 as Q
import thermal_sheaf_filtration as tsf

FAMILY, N_PART, MATRIX, SEED = "bimodal_5to1", 2000, "epoxy", 0
CASES = [
    (25.0, "SiO2", "treated"),
    (25.0, "SiO2", "untreated"),
    (25.0, "AlN",  "treated"),
    ( 5.0, "SiO2", "treated"),
]
SCALES = [1e4, 1e6, 1e8, 1e10, 1e12]
EPS = [1e-20, 1e-30, 1e-40]

CSV = _pathlib.Path("results_r48_saturation.csv")
TXT = _pathlib.Path("results_r48_saturation.txt")
OUT = []
def say(*a):
    line = " ".join(str(x) for x in a); print(line, flush=True); OUT.append(line)


def cfg():
    return Q.Q24Config(eta_cut=1.0, slab_axis=2, slab_frac=0.06)


def ideal_table(filler, treat, kind, val):
    t = MR.pair(filler, treat, filler, treat, MATRIX)
    if kind == "interface":
        for k in t.materials:
            t.materials[k].R_s_matrix = val
            t.materials[k].R_s_contact = val
    elif kind == "matrix":
        t.k_matrix = t.k_matrix * val
    elif kind == "filler":
        for k in t.materials:
            t.materials[k].kappa = t.materials[k].kappa * val
    return t


def one(geo, d_um, filler, treat):
    ids, typ, pos, dia, mol, info = geo
    c = cfg()
    base = MR.pair(filler, treat, filler, treat, MATRIX)
    G0net, _, sh0, _ = Q._net(ids, typ, pos.copy(), dia, mol, base, c, d_um)
    G0 = float(tsf.effective_conductance(sh0, 0.0))
    rows = []
    for kind, vals in (("interface", EPS), ("matrix", SCALES), ("filler", SCALES)):
        for v in vals:
            t = ideal_table(filler, treat, kind, v)
            Gi, _, shi, _ = Q._net(ids, typ, pos.copy(), dia, mol, t, c, d_um)
            g = float(tsf.effective_conductance(shi, 0.0))
            rows.append(dict(d_um=d_um, filler=filler, treat=treat, matrix=MATRIX,
                             family=FAMILY, n_part=N_PART, seed=SEED,
                             ideal=kind, scale=v, G_base=G0, G_ideal=g,
                             gamma=g / G0,
                             n_edges_base=G0net.number_of_edges(),
                             n_edges_ideal=Gi.number_of_edges()))
    return rows


def main() -> None:
    rows = []
    for d_um, fil, tr in CASES:
        geo = EN.make(FAMILY, d_um, SEED, n=N_PART)
        rows += one(geo, d_um, fil, tr)
        say(f"  {fil}/{tr}/d={d_um:g}um done")
    df = pd.DataFrame(rows)
    df.to_csv(CSV, index=False)

    say("=" * 100)
    say("R48: does the stand-in for 'infinite / zero' in the idealised Gamma saturate?")
    say("=" * 100)
    say(f"structure {FAMILY} N={N_PART} seed={SEED} / {MATRIX} / "
        f"engine {tsf.engine_stamp()['engine_variant']} {tsf.engine_stamp()['baseline_fp']}")
    say("WARNING: the defaults in `quantities24.py` are: interface 1e-30 / matrix and particle x1e6")
    say("")
    for (d_um, fil, tr), g in df.groupby(["d_um", "filler", "treat"], sort=False):
        say(f"--- {fil} / {tr} / d={d_um:g} um ---")
        for kind, gg in g.groupby("ideal", sort=False):
            gg = gg.sort_values("scale")
            lim = gg.gamma.iloc[0] if kind == "interface" else gg.gamma.iloc[-1]
            cur = gg[gg.scale.isin([1e-30, 1e6])].gamma.iloc[0]
            cells = "  ".join(f"{s:.0e}:{v:.4f}" for s, v in zip(gg.scale, gg.gamma))
            dev = abs(cur / lim - 1)
            ed = gg.n_edges_ideal.nunique()
            say(f"  {kind:<10} {cells}")
            say(f"  {'':<10} converged value {lim:.4f} / **deviation of the default {dev:.3%}**"
                f" / distinct edge counts {ed} (1 means the edge set is unchanged)")
        say("")

    say("--- verdict ---")
    bad, worst = [], (0.0, None)
    for (d_um, fil, tr, kind), g in df.groupby(["d_um", "filler", "treat", "ideal"]):
        g = g.sort_values("scale")
        lim = g.gamma.iloc[0] if kind == "interface" else g.gamma.iloc[-1]
        cur = g[g.scale.isin([1e-30, 1e6])].gamma.iloc[0]
        dev = abs(cur / lim - 1)
        if dev > worst[0]:
            worst = (dev, f"{fil}/{tr}/d={d_um:g} {kind}")
        if dev > 0.01:
            bad.append((fil, tr, d_um, kind, dev))
    say(f"  * largest deviation between default and converged value: **{worst[1]}  {worst[0]:.3%}**")
    if bad:
        say(f"  WARNING: **{len(bad)} cases deviate from the default by more than 1%**")
        for fil, tr, d_um, kind, dev in bad:
            say(f"      {fil}/{tr}/d={d_um:g} {kind}: {dev:.3%}")
    else:
        say("  PASS **all conditions and all operations within 1% of the default** (= the stand-in is sufficient)")
    say(f"  cases in which the edge set changed with the scale factor: "
        f"{int((df.groupby(['d_um','filler','treat','ideal']).n_edges_ideal.nunique()>1).sum())} cases")
    TXT.write_text("\n".join(OUT) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
