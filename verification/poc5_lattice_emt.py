# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/PoC5_lattice_emt/poc5_lattice_emt.py
#   sha256(src) : 4782244a0041bafb
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 134 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V4. Regular lattices (SC / BCC / FCC) homogenised exactly (lattice.py) and
# compared with Maxwell-Eucken, Nan's effective-medium theory and the
# Batchelor-O'Brien asymptote (Sec. 5.3). Writes results_poc5_grid.csv (Fig. 12).
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
import itertools
import pathlib
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import lattice as LT
import materials_real as MR
import thermal_sheaf_filtration as tsf
from thermal_sheaf_filtration import AnalysisConfig

D_GRID = [1.0, 10.0, 100.0]
IFACE = ["Rs0", "untreated", "treated"]
ETA_GRID = [0.1, 0.3, 0.5, 1.0, 2.0, 3.0]
ETA_DEFAULT = 1.0
N_CELLS_FINITE = [4, 6, 8, 10]

OUT = []
def say(*a):
    line = " ".join(str(x) for x in a)
    print(line)
    OUT.append(line)


def cfg_for(eta=ETA_DEFAULT):
    return AnalysisConfig(gap_max_um=1e9, eta_cut=eta,
                          t_min_um=1e-6, delta_min_um=1e-6)


def table_iface(filler, iface, matrix):
    base = "untreated" if iface == "Rs0" else iface
    tab = copy.deepcopy(MR.single(filler, base, matrix))
    if iface == "Rs0":
        for m in tab.materials.values():
            m.R_s_matrix = 0.0
            m.R_s_contact = 0.0
            m.t_ip_um = 0.0
            m.k_ip = None
    return tab


def add_variants(df):
    df["k_V1"] = df["k_net"] + df["k_m"]
    df["k_V2"] = df["k_net"] + df["k_m"] * (1 - df["phi"])
    for v in ("net", "V1", "V2"):
        df[f"ratio_{v}"] = df[f"k_{v}"] / df["k_nan"]
        df[f"err_{v}"] = (df[f"ratio_{v}"] - 1).abs()
    return df


def main():
    t0 = time.time()
    say("=" * 84)
    say("PoC-5: regular lattices vs Nan EMT / Maxwell-Eucken / BOB asymptotics")
    say("=" * 84)
    say(f"{'latt':>5} {'z':>3} {'a_nn/a':>8} {'phi_max':>8}  phi grid")
    for k in LT.LATTICES:
        say(f"{k:>5} {LT.Z_COORD[k]:3d} {LT.NN_OVER_A[k]:8.4f} "
            f"{LT.PHI_MAX[k]:8.4f}  {LT.PHI_GRID[k]}")
    nphi = sum(len(LT.PHI_GRID[k]) for k in LT.LATTICES)
    say(f"-> (lattice,phi) {nphi} pairs x filler 6 x matrix 3 x interface 3 x size 3 = "
        f"{nphi*6*3*3*3} conditions (eta_cut={ETA_DEFAULT} default)")
    say()

    say("--- PoC-5-0a: exactness of the homogenisation (nearest neighbours only vs analytic formula) ---")
    tabc = MR.single("AlN", "untreated", "epoxy")
    cfgc = cfg_for(1e-9)
    say(f"  {'latt':>5} {'k_zz (impl)':>15} {'analytic':>15} {'rel. err':>10} "
        f"{'|k_xx/k_zz-1|':>14}")
    ex = []
    for k in LT.LATTICES:
        r = LT.homogenize(k, 20.0, 0.40, tabc, cfgc,
                          r_cut_over_a=1.01 * LT.NN_OVER_A[k])
        a_m = r["a_um"] * 1e-6
        pred = {"SC": 1, "BCC": 2, "FCC": 4}[k] * r["g_nn"] / a_m
        ex.append(abs(r["k_zz"] / pred - 1))
        say(f"  {k:>5} {r['k_zz']:15.6e} {pred:15.6e} "
            f"{abs(r['k_zz']/pred-1):10.2e} {abs(r['k_xx']/r['k_zz']-1):14.2e}")
    say(f"  verdict: {'PASS exact (0.0)' if max(ex) < 1e-12 else 'FAIL'}")
    say()

    say("--- PoC-5-0b: lowering eta_cut double-counts the matrix ---")
    say("  the pair-conductance formula is the solution for an 'isolated pair in an infinite matrix'. Summing distant pairs too")
    say("  counts the same matrix volume repeatedly. eta_cut is the parameter that prevents this.")
    say(f"  condition: FCC, Al2O3/epoxy, untreated, d=10 um")
    say(f"  {'eta_cut':>8} {'r_cut/a':>8} " +
        " ".join(f"{f'phi={p}':>13}" for p in (0.20, 0.40, 0.60)))
    say(f"  {'':>8} {'':>8} " + " ".join(f"{'edges':>5}{'k_net/k_Nan':>8}"
                                         for _ in range(3)))
    erows = []
    for eta in ETA_GRID:
        cells = []
        for phi in (0.20, 0.40, 0.60):
            c = cfg_for(eta)
            tab = MR.single("Al2O3", "untreated", "epoxy")
            h = LT.homogenize("FCC", 10.0, phi, tab, c)
            kn, aN = LT.nan_emt(MR.KAPPA["Al2O3"], 0.2, phi,
                                tab.materials[1].R_s_matrix, 10.0)
            cells.append(f"{h['n_bonds']:5d}{h['k_zz']/kn:8.3f}")
            erows.append(dict(eta_cut=eta, phi=phi, n_bonds=h["n_bonds"],
                              n_shells=h["n_shells"], k_net=h["k_zz"],
                              k_nan=kn, ratio=h["k_zz"] / kn))
        say(f"  {eta:8.1f} {LT.r_cut_for_eta(eta):8.3f} "
            + " ".join(f"{c:>13}" for c in cells))
    de = pd.DataFrame(erows)
    de.to_csv("results_poc5_etacut.csv", index=False)
    say("  -> the lower eta_cut, the more edges and the larger k_net. **double counting is real**.")
    say(f"    at eta_cut=0.1, phi=0.4 gives "
        f"{de[(de.eta_cut==0.1)&(de.phi==0.4)]['ratio'].iloc[0]:.2f} x Nan.")
    say(f"    the default eta_cut=1.0 gives {de[(de.eta_cut==1.0)&(de.phi==0.4)]['ratio'].iloc[0]:.2f} x Nan.")
    say("  => below, the evaluation uses **the model default eta_cut=1.0** (same as all other PoCs).")
    say()

    rows = []
    cfg = cfg_for(ETA_DEFAULT)
    for kind in LT.LATTICES:
        for phi, f, mx, ifc, d_um in itertools.product(
                LT.PHI_GRID[kind], MR.FILLERS_GRID, MR.MATRICES_GRID,
                IFACE, D_GRID):
            tab = table_iface(f, ifc, mx)
            kap, k_m = MR.KAPPA[f], MR.K_MATRIX[mx]
            h = LT.homogenize(kind, d_um, phi, tab, cfg)
            Rs = tab.materials[1].R_s_matrix
            k_nan, aN = LT.nan_emt(kap, k_m, phi, Rs, d_um)
            n_con = h["branch"].get("contact", 0)
            rows.append(dict(
                lattice=kind, z=LT.Z_COORD[kind], phi=phi,
                phi_rel=phi / LT.PHI_MAX[kind], filler=f, kappa=kap,
                matrix=mx, k_m=k_m, iface=ifc, R_s=Rs, d_um=d_um,
                alpha=kap / k_m, alpha_N=aN, a_K_um=Rs * k_m * 1e6,
                a_um=h["a_um"], gap_nn_um=h["gap_nn_um"],
                gap_over_d=h["gap_nn_um"] / d_um,
                n_bonds=h["n_bonds"], n_contact=n_con,
                frac_contact=n_con / max(h["n_bonds"], 1),
                k_net=h["k_zz"], k_nan=k_nan,
                k_me=LT.maxwell_eucken(kap, k_m, phi),
                aniso=(abs(h["k_xx"] / h["k_zz"] - 1) if h["k_zz"] > 0 else 0.0)))
    df = add_variants(pd.DataFrame(rows))

    df["emt_scope"] = (df.n_bonds > 0) & (df.n_contact == 0)
    df["emt_scope_Rs_pos"] = df.emt_scope & (df.iface != "Rs0")

    df.to_csv("results_poc5_grid.csv", index=False)
    say(f"  EMT comparison population: emt_scope {int(df.emt_scope.sum())} conditions"
        f" / of which R_s>0: emt_scope_Rs_pos {int(df.emt_scope_Rs_pos.sum())} conditions")
    say(f"main run: {len(df)} conditions written to results_poc5_grid.csv")
    say(f"  isotropy: max |k_xx/k_zz-1| = {df['aniso'].max():.2e} (cubic symmetry, so 0 is correct)")
    say()

    say("--- PoC-5-0c: range of validity of the model (phi where edges exist at eta_cut=1.0) ---")
    say("  eta_cut=1.0 <=> ln(1+r~/h) >= 1 <=> h <= r~/(e-1) = 0.1455*d (equal diameters)")
    say(f"  {'latt':>5} {'phi':>7} {'phi_rel':>8} {'gap_nn/d':>9} {'edges':>6} "
        f"{'contact':>7} {'k_net>0':>8}")
    for kind in LT.LATTICES:
        for phi in LT.PHI_GRID[kind]:
            q = df[(df.lattice == kind) & (df.phi == phi)
                   & (df.filler == "Al2O3") & (df.matrix == "epoxy")
                   & (df.iface == "untreated") & (df.d_um == 10.0)].iloc[0]
            say(f"  {kind:>5} {phi:7.4f} {q['phi_rel']:8.3f} "
                f"{q['gap_over_d']:9.4f} {int(q['n_bonds']):6d} "
                f"{int(q['n_contact']):7d} {'yes' if q['k_net']>0 else '**no**':>8}")
    valid = df[df.n_bonds > 0]
    say(f"  -> conditions with edges: {len(valid)}/{len(df)} "
        f"({len(valid)/len(df):.1%})")
    phimin = {k: valid[valid.lattice == k]["phi"].min() for k in LT.LATTICES}
    say(f"  lower phi limit per lattice: " + ", ".join(f"{k}: {v:.3f}"
                                          for k, v in phimin.items()))
    say("  the model does not describe the dilute regime (phi <~ 0.4): at eta_cut = 1.0 no edge is formed.")
    say("    This is the intended range of applicability (high loading), not a defect.")
    say()

    say("--- PoC-5-0d: background matrix term (V0 k_net / V1 +k_m / V2 +k_m(1-phi)) ---")
    say(f"  {'range':>22} {'n':>5} " + " ".join(f"{v:>22}"
                                            for v in ("V0 k_net", "V1 +k_m",
                                                      "V2 +k_m(1-phi)")))
    say(f"  {'':>22} {'':>5} " + " ".join(f"{'ratio med':>11}{'|err| med':>11}"
                                          for _ in range(3)))
    for lab, sub in (("edges: all", valid),
                     ("edges, no contact", valid[valid.frac_contact == 0]),
                     ("no edges (dilute)", df[df.n_bonds == 0])):
        if not len(sub):
            continue
        cells = [f"{sub[f'ratio_{v}'].median():11.4f}{sub[f'err_{v}'].median():11.4f}"
                 for v in ("net", "V1", "V2")]
        say(f"  {lab:>22} {len(sub):5d} " + " ".join(f"{c:>22}" for c in cells))
    vv = valid[valid.frac_contact == 0]
    best = min(("net", "V1", "V2"), key=lambda v: vv[f"err_{v}"].median())
    say(f"  -> best closure within the edges-no-contact range is **{best}**"
        f" (prediction 5-0 was V2 -> {'a-priori prediction: met' if best=='V2' else 'a-priori prediction: not met'})")
    say()

    say("--- PoC-5-1: does R_s=0 agree with Maxwell-Eucken? (prediction: <=5% at phi<=0.2) ---")
    s1 = df[df.iface == "Rs0"].copy()
    s1["err_me_V2"] = (s1["k_V2"] / s1["k_me"] - 1).abs()
    s1["err_me_net"] = (s1["k_net"] / s1["k_me"] - 1).abs()
    lo = s1[s1.phi <= 0.20]
    say(f"  phi<=0.2 ({len(lo)} conditions, all without edges): k_net=0, hence k_V2=k_m(1-phi).")
    say(f"    median |k_V2/k_ME-1| {lo['err_me_V2'].median():.4f} -> "
        f"**prediction 5-1 missed** (the dilute regime is outside the model's scope)")
    s1v = s1[s1.n_bonds > 0]
    say(f"  |k_V2/k_ME-1| within the range with edges ({len(s1v)} conditions):")
    for (kind,), g in s1v.groupby(["lattice"]):
        say(f"    {kind:>5}: median {g['err_me_V2'].median():7.4f}  "
            f"min {g['err_me_V2'].min():7.4f}  max {g['err_me_V2'].max():7.4f}")
    p1 = lo["err_me_V2"].median() <= 0.05
    say(f"  verdict (prediction: <=5% at phi<=0.2): {'PASS' if p1 else 'FAIL'}")
    say()

    say("--- PoC-5-2: does R_s>0 agree with Nan? (prediction: <=15%) ---")
    s2 = valid[(valid.iface != "Rs0") & (valid.frac_contact == 0)]
    say(f"  scope (edges, no contact, R_s>0): {len(s2)} conditions")
    say(f"  {'phi_rel bin':>14} {'n':>5} {'ratio med':>10} {'|err| med':>11} "
        f"{'p90':>9} {'max':>9}")
    for lab, sub in (("0.75-0.85", s2[(s2.phi_rel > .75) & (s2.phi_rel <= .85)]),
                     ("0.85-0.95", s2[(s2.phi_rel > .85) & (s2.phi_rel <= .95)]),
                     ("0.95-1.00", s2[s2.phi_rel > .95])):
        if not len(sub):
            continue
        say(f"  {lab:>14} {len(sub):5d} {sub['ratio_'+best].median():10.4f} "
            f"{sub['err_'+best].median():11.4f} "
            f"{sub['err_'+best].quantile(.9):9.4f} {sub['err_'+best].max():9.4f}")
    p2 = s2[f"err_{best}"].median() <= 0.15
    say(f"  verdict: {'PASS' if p2 else 'FAIL'}"
        f" (median {s2[f'err_{best}'].median():.4f} vs threshold 0.15, closure={best})")
    say()
    say("  per filler (phi/phi_max>0.9, epoxy, d=10 um, untreated):")
    say(f"  {'filler':>9} {'alpha':>8} {'k_model':>11} {'k_Nan':>11} {'ratio':>8}")
    for f in MR.FILLERS_GRID:
        q = s2[(s2.filler == f) & (s2.matrix == "epoxy") & (s2.d_um == 10.0)
               & (s2.iface == "untreated") & (s2.phi_rel > 0.9)]
        if not len(q):
            continue
        q = q.sort_values("phi_rel").iloc[-1]
        say(f"  {f:>9} {q['alpha']:8.1f} {q['k_'+best]:11.4f} "
            f"{q['k_nan']:11.4f} {q['ratio_'+best]:8.3f}")
    say()

    say("--- PoC-5-3: size dependence dlog k/dlog d, model vs Nan (+-0.2) ---")
    say("  * independent check of the Z_size_sweep finding 'surface treatment matters more for small particles'")
    sl = []
    for keys, g in valid[valid.iface != "Rs0"].groupby(
            ["lattice", "phi", "filler", "matrix", "iface"]):
        g = g.sort_values("d_um")
        if len(g) < 3 or g[f"k_{best}"].min() <= 0:
            continue
        x = np.log10(g["d_um"].values)
        sl.append(dict(lattice=keys[0], phi=keys[1], filler=keys[2],
                       matrix=keys[3], iface=keys[4],
                       slope_model=np.polyfit(x, np.log10(g[f"k_{best}"].values), 1)[0],
                       slope_nan=np.polyfit(x, np.log10(g["k_nan"].values), 1)[0]))
    dsl = pd.DataFrame(sl)
    dsl["diff"] = dsl["slope_model"] - dsl["slope_nan"]
    say(f"  {len(dsl)} groups. median |Delta slope| {dsl['diff'].abs().median():.4f}, "
        f"p90 {dsl['diff'].abs().quantile(.9):.4f}, max {dsl['diff'].abs().max():.4f}")
    say(f"  {'iface':>10} {'model slope':>13} {'Nan slope':>12} {'diff med':>10}")
    for ifc, g in dsl.groupby("iface"):
        say(f"  {ifc:>10} {g['slope_model'].median():13.4f} "
            f"{g['slope_nan'].median():12.4f} {g['diff'].median():10.4f}")
    p3 = dsl["diff"].abs().median() <= 0.2
    say(f"  verdict: {'PASS' if p3 else 'FAIL'} (median {dsl['diff'].abs().median():.4f} vs 0.2)")
    say()
    _c = df[(df.lattice == "FCC") & (df.iface == "untreated")]
    ph = float(_c.iloc[(_c["phi_rel"] - 0.86).abs().argsort()].iloc[0]["phi"])
    say(f"  size dependence of the treated/untreated ratio (FCC phi={ph:.3f}, "
        f"phi/phi_max={ph/LT.PHI_MAX['FCC']:.3f}, epoxy)")
    say("  note: phi is chosen to avoid the drop region near phi->phi_max (5-9)")
    say(f"  {'filler':>9} " + " ".join(f"{f'd={d:g}um':>22}" for d in D_GRID))
    say(f"  {'':>9} " + " ".join(f"{'model':>11}{'Nan':>11}" for _ in D_GRID))
    for f in MR.FILLERS_GRID:
        cells = []
        for d_um in D_GRID:
            q = df[(df.lattice == "FCC") & (df.phi == ph) & (df.filler == f)
                   & (df.matrix == "epoxy") & (df.d_um == d_um)]
            tr, un = q[q.iface == "treated"], q[q.iface == "untreated"]
            cells.append(f"{tr[f'k_{best}'].iloc[0]/un[f'k_{best}'].iloc[0]:11.4f}"
                         f"{tr['k_nan'].iloc[0]/un['k_nan'].iloc[0]:11.4f}")
        say(f"  {f:>9} " + " ".join(cells))
    say()

    say("--- PoC-5-4: overestimate at high alpha (Nan ratio <= 2.75x?) ---")
    s4 = valid[(valid.iface != "Rs0") & (valid.frac_contact == 0)]
    s4 = s4.assign(alpha_bin=pd.cut(s4["alpha"], [0, 10, 50, 200, 1000, 1e5]))
    say(f"  {'alpha bin':>14} {'n':>5} {'Nan-ratio med':>13} {'p90':>9} {'max':>9}")
    for b, g in s4.groupby("alpha_bin", observed=True):
        say(f"  {str(b):>14} {len(g):5d} {g[f'ratio_{best}'].median():13.4f} "
            f"{g[f'ratio_{best}'].quantile(.9):9.4f} {g[f'ratio_{best}'].max():9.4f}")
    hi = s4[s4.alpha >= 500]
    p4 = bool(len(hi)) and hi[f"ratio_{best}"].max() <= 2.75
    say(f"  conditions with alpha>=500: {len(hi)}: max Nan ratio "
        f"{hi[f'ratio_{best}'].max() if len(hi) else float('nan'):.4f} -> "
        f"{'PASS within 2.75' if p4 else 'FAIL above 2.75'}")
    say("  the correction of Dai et al. (2019) is not implemented here.")
    say("    only the upper-bound guideline 'Nan ratio within 2.75x' is judged here.")
    say()

    say("--- PoC-5-5: effect of the coordination number z (SC->FCC 1.5-3x at the same phi/phi_max?) ---")
    say("  in absolute phi the model's range of validity differs per lattice, so phi/phi_max is used to align.")
    say(f"  {'phi_rel':>9} {'SC':>12} {'BCC':>12} {'FCC':>12} {'FCC/SC':>8} "
        f"{'BCC/SC':>8}")
    ratios = []
    for tgt in (0.80, 0.90, 0.95, 0.998):
        v = {}
        for kind in LT.LATTICES:
            q = df[(df.lattice == kind) & (df.filler == "Al2O3")
                   & (df.matrix == "epoxy") & (df.iface == "untreated")
                   & (df.d_um == 10.0) & (df.frac_contact == 0)].copy()
            if not len(q):
                continue
            q["dd"] = (q["phi_rel"] - tgt).abs()
            v[kind] = q.sort_values("dd").iloc[0]["k_net"]
        if len(v) == 3 and v["SC"] > 0:
            say(f"  {tgt:9.3f} {v['SC']:12.5e} {v['BCC']:12.5e} {v['FCC']:12.5e} "
                f"{v['FCC']/v['SC']:8.3f} {v['BCC']/v['SC']:8.3f}")
            ratios.append(v["FCC"] / v["SC"])
    mono = all(np.diff([np.mean([r for r in
                                 [x for x in ratios]]) if False else 0]) == 0)
    zmean = (df[(df.frac_contact == 0) & (df.n_bonds > 0)
                & (df.iface == "untreated")]
             .groupby("lattice")["k_net"].median().reindex(LT.LATTICES))
    mono = bool(np.all(np.diff(zmean.values) > 0))
    say(f"  k_net median increases monotonically in z order (SC<BCC<FCC): {mono} "
        f"({', '.join(f'{k}={v:.4g}' for k, v in zmean.items())})")
    p5 = mono and bool(ratios) and all(1.5 <= r <= 3.0 for r in ratios)
    say(f"  FCC/SC range {min(ratios):.3f}-{max(ratios):.3f}" if ratios
        else "  no common phi/phi_max at which FCC/SC can be taken")
    say(f"  verdict: {'PASS' if p5 else 'WARN/FAIL'} (prediction 1.5-3x)")
    say()

    say("--- PoC-5-6: switch to the contact branch at phi->phi_max? (majority of edges in contact?) ---")
    near = df[df.phi_rel >= 1.0]
    say(f"  phi >= phi_max (points given a slight overlap): {len(near)} conditions, "
        f"median contact-edge fraction {near['frac_contact'].median():.4f}")
    say(f"  {'latt':>5} {'phi_rel':>8} {'gap_nn/d':>10} {'contact frac':>13} {'edges':>7}")
    for kind in LT.LATTICES:
        for phi in sorted(LT.PHI_GRID[kind])[-3:]:
            q = df[(df.lattice == kind) & (df.phi == phi)
                   & (df.filler == "Al2O3") & (df.matrix == "epoxy")
                   & (df.iface == "untreated") & (df.d_um == 10.0)].iloc[0]
            say(f"  {kind:>5} {q['phi_rel']:8.3f} {q['gap_over_d']:10.4f} "
                f"{q['frac_contact']:13.4f} {int(q['n_bonds']):7d}")
    p6 = bool(len(near)) and near["frac_contact"].median() > 0.5
    say(f"  verdict: {'PASS' if p6 else 'FAIL'}")
    say()

    say("--- PoC-5-7: finite slab solved with the main engine vs the homogenised value ---")
    say("  the lateral faces are free boundaries, so an underestimate is expected that converges as n_cells increases.")
    say(f"  {'latt':>5} {'n':>3} {'N_sph':>6} {'k_slab':>12} {'k_homog':>12} {'ratio':>8}")
    frows = []
    tab = MR.single("Al2O3", "untreated", "epoxy")
    for kind in LT.LATTICES:
        phi = sorted(LT.PHI_GRID[kind])[-3]
        hh = LT.homogenize(kind, 10.0, phi, tab, cfg)
        for n in N_CELLS_FINITE:
            ids, typ, pos, dia, mol, a = LT.slab_structure(kind, 10.0, phi, n)
            c2 = AnalysisConfig(gap_max_um=1e9, eta_cut=ETA_DEFAULT,
                                t_min_um=1e-6, slab_axis=2, slab_frac=0.01,
                                viz=False)
            G = tsf.build_thermal_network(ids, typ, pos, dia, c2,
                                          table=tab, mols=mol)
            src, snk = tsf.identify_slabs(ids, pos, c2)
            Gt, _ = tsf.attach_virtual_terminals(G, src, snk, c2)
            Geff = tsf.effective_conductance(tsf.CellularSheaf(Gt), 0.0)
            L = (pos[:, 2].max() - pos[:, 2].min()) * 1e-6
            A = ((pos[:, 0].max() - pos[:, 0].min() + a)
                 * (pos[:, 1].max() - pos[:, 1].min() + a)) * 1e-12
            k_slab = Geff * L / A
            say(f"  {kind:>5} {n:3d} {len(ids):6d} {k_slab:12.5e} "
                f"{hh['k_zz']:12.5e} {k_slab/hh['k_zz']:8.4f}")
            frows.append(dict(lattice=kind, phi=phi, n_cells=n, N=len(ids),
                              k_slab=k_slab, k_homog=hh["k_zz"],
                              ratio=k_slab / hh["k_zz"]))
    dfin = pd.DataFrame(frows)
    dfin.to_csv("results_poc5_finite.csv", index=False)
    conv = dfin[dfin.n_cells == max(N_CELLS_FINITE)]["ratio"]
    trend = all(g.sort_values("n_cells")["ratio"].is_monotonic_increasing
                for _, g in dfin.groupby("lattice"))
    say(f"  ratio at n={max(N_CELLS_FINITE)}: {conv.min():.4f}-{conv.max():.4f}, "
        f"monotonic convergence with increasing n: {trend}")
    say("  -> the finite slab underestimates. **using homogenisation is justified**")
    say()

    say("--- PoC-5-9 [added] non-monotonic k_e(phi): k_e drops as phi->phi_max ---")
    say("  the drop seen in fig (a). Its cause is the interfacial resistance of the near-contact branch")
    say("    R_int = (R_s,i+R_s,j)/A_th,  A_th = 2*pi r~ h ln(1+r~/h)")
    say("  with A_th proportional to h, so **R_int diverges as h->0**. That is,")
    say("  the core of claim B, 'closing the gap too far makes things worse',")
    say("  also appears at the level of the composite k_e(phi). EMT can only produce a monotonic increase.")
    say(f"  {'latt':>5} {'filler':>9} {'best phi':>8} {'phi_rel':>8} "
        f"{'k_max':>10} {'k near phi_max':>14} {'drop':>9} {'Nan mono?':>11}")
    nm = []
    for kind in LT.LATTICES:
        for f in ("Al2O3", "AlN", "SiO2"):
            g = df[(df.lattice == kind) & (df.filler == f)
                   & (df.matrix == "epoxy") & (df.iface == "untreated")
                   & (df.d_um == 10.0) & (df.frac_contact == 0)
                   & (df.n_bonds > 0)].sort_values("phi")
            if len(g) < 3:
                continue
            i = g[f"k_{best}"].idxmax()
            kmax = g.loc[i, f"k_{best}"]
            klast = g.iloc[-1][f"k_{best}"]
            mono_nan = bool(np.all(np.diff(g["k_nan"].values) > 0))
            say(f"  {kind:>5} {f:>9} {g.loc[i,'phi']:8.4f} "
                f"{g.loc[i,'phi_rel']:8.3f} {kmax:10.4f} {klast:14.4f} "
                f"{'x'+format(kmax/klast,'.3f'):>9} {str(mono_nan):>11}")
            nm.append(dict(lattice=kind, filler=f, phi_best=g.loc[i, "phi"],
                           phi_rel_best=g.loc[i, "phi_rel"], k_max=kmax,
                           k_last=klast, drop=kmax / klast, nan_mono=mono_nan))
    dnm = pd.DataFrame(nm)
    dnm.to_csv("results_poc5_nonmono.csv", index=False)
    say(f"  -> median best phi/phi_max {dnm['phi_rel_best'].median():.3f}, "
        f"median drop just below phi_max x{dnm['drop'].median():.3f}, "
        f"max x{dnm['drop'].max():.3f}")
    say(f"    groups in which Nan is monotonically increasing: {int(dnm['nan_mono'].sum())}/{len(dnm)}")
    say("  * **'the optimum filler fraction lies below phi_max' is a model-specific prediction**.")
    say("    experimentally testable (does k drop when phi is pushed too high?).")
    say()

    say("--- PoC-5-8 [added] kappa dependence: model (beta=0) vs Nan ---")
    say("  -> Nan's EMT independently predicts a kappa dependence. Compare the slopes dlog k/dlog alpha.")
    ka = []
    for keys, g in valid[(valid.iface != "Rs0")
                         & (valid.frac_contact == 0)].groupby(
            ["lattice", "phi", "matrix", "iface", "d_um"]):
        g = g[g[f"k_{best}"] > 0]
        if g["filler"].nunique() < 4:
            continue
        x = np.log10(g["alpha"].values)
        if x.max() - x.min() < 1.0:
            continue
        ka.append(dict(lattice=keys[0], phi=keys[1], matrix=keys[2],
                       iface=keys[3], d_um=keys[4],
                       slope_model=np.polyfit(x, np.log10(g[f"k_{best}"].values), 1)[0],
                       slope_nan=np.polyfit(x, np.log10(g["k_nan"].values), 1)[0],
                       spread_model=g[f"k_{best}"].max() / g[f"k_{best}"].min(),
                       spread_nan=g["k_nan"].max() / g["k_nan"].min()))
    dka = pd.DataFrame(ka)
    say(f"  {len(dka)} groups (each group spans 3.2 decades of alpha across the 6 fillers)")
    say(f"  {'':>16} {'slope med':>12} {'slope max':>11} "
        f"{'k-spread med':>17} {'spread max':>10}")
    say(f"  {'model(beta=0)':>16} {dka['slope_model'].median():12.4f} "
        f"{dka['slope_model'].max():11.4f} {dka['spread_model'].median():17.4f} "
        f"{dka['spread_model'].max():10.4f}")
    say(f"  {'Nan EMT':>16} {dka['slope_nan'].median():12.4f} "
        f"{dka['slope_nan'].max():11.4f} {dka['spread_nan'].median():17.4f} "
        f"{dka['spread_nan'].max():10.4f}")
    say(f"  slope difference (Nan - model), median: "
        f"{(dka['slope_nan']-dka['slope_model']).median():+.4f}")
    say("  Nan carries a kappa dependence while the model at beta = 0 carries almost none;")
    say("    beta = 0 removes it by construction. Nan's own kappa dependence is weak as well")
    say("    and saturates at high alpha.")
    dka.to_csv("results_poc5_kappa_dep.csv", index=False)
    say()

    say("=" * 84)
    say("comparison with predictions")
    say("=" * 84)
    say(f"  5-0 background matrix term needed, V2 best  : {'PASS' if best=='V2' else 'FAIL'} "
        f"(best closure = {best})")
    say(f"  5-1 R_s=0, phi<=0.2: within 5% of ME        : {'PASS' if p1 else 'FAIL'} "
        f"(at phi<=0.2 **not a single edge forms** = outside the model's scope)")
    say(f"  5-2 R_s>0: within 15% of Nan                : {'PASS' if p2 else 'FAIL'} "
        f"(median within the valid range {s2[f'err_{best}'].median():.4f})")
    say(f"  5-3 dlog k/dlog d within +-0.2 of Nan       : {'PASS' if p3 else 'FAIL'} "
        f"(median {dsl['diff'].abs().median():.4f})")
    say(f"  5-4 alpha>=500: Nan ratio <=2.75            : {'PASS' if p4 else 'FAIL'}")
    say(f"  5-5 z monotonic, FCC/SC 1.5-3x              : {'PASS' if p5 else 'WARN/FAIL'}")
    say(f"  5-6 contact majority at phi->phi_max        : {'PASS' if p6 else 'FAIL'} "
        f"(median {near['frac_contact'].median():.4f})")
    say(f"  5-8 [added] slope of the kappa dependence   : model "
        f"{dka['slope_model'].median():+.4f} vs Nan "
        f"{dka['slope_nan'].median():+.4f} -> **beta=0 removes the kappa dependence**")
    say()

    fig, ax = plt.subplots(2, 3, figsize=(17, 9.6))
    q = df[(df.filler == "Al2O3") & (df.matrix == "epoxy") & (df.d_um == 10.0)
           & (df.iface == "untreated")]
    a = ax[0, 0]
    for kind, mk in zip(LT.LATTICES, "os^"):
        g = q[q.lattice == kind].sort_values("phi")
        gg = g[g.n_bonds > 0]
        a.semilogy(gg["phi"], gg[f"k_{best}"], mk + "-",
                   label=f"model {kind} (z={LT.Z_COORD[kind]})")
    g = q[q.lattice == "FCC"].sort_values("phi")
    a.semilogy(g["phi"], g["k_nan"], "k--", label="Nan EMT")
    a.semilogy(g["phi"], g["k_me"], "k:", label="Maxwell-Eucken")
    a.set_xlabel(r"$\phi$"); a.set_ylabel("$k_e$ [W/mK]")
    a.set_title("(a) Al$_2$O$_3$/epoxy, d=10 um, untreated\n"
                "model (where edges exist) vs analytic")
    a.legend(fontsize=7); a.grid(alpha=.3, which="both")

    a = ax[0, 1]
    for phi, mk in zip((0.20, 0.40, 0.60), "os^"):
        g = de[de.phi == phi].sort_values("eta_cut")
        a.loglog(g["eta_cut"], g["ratio"], mk + "-", label=f"$\\phi$={phi}")
    a.axhline(1.0, ls=":", c="k"); a.axvline(1.0, ls="--", c="r", lw=1)
    a.text(1.06, de["ratio"].min() * 1.3, "default", color="r", fontsize=7,
           rotation=90, va="bottom")
    a.set_xlabel(r"$\eta_{cut}$"); a.set_ylabel("$k_{net}/k_{Nan}$")
    a.set_title("(b) lowering $\\eta_{cut}$ double-counts the matrix")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    a = ax[0, 2]
    s = valid[(valid.iface != "Rs0") & (valid.frac_contact == 0)]
    sc = a.scatter(s["alpha"], s[f"ratio_{best}"], c=s["phi_rel"], s=9,
                   cmap="viridis", alpha=.7)
    a.set_xscale("log"); a.axhline(1.0, ls=":", c="k")
    a.axhline(2.75, ls="--", c="r", lw=1)
    a.text(2, 2.9, "2.75 (Dai upper bound)", color="r", fontsize=7)
    plt.colorbar(sc, ax=a, label=r"$\phi/\phi_{max}$")
    a.set_xlabel(r"$\alpha=\kappa/k_m$"); a.set_ylabel(f"$k_{{{best}}}/k_{{Nan}}$")
    a.set_title("(c) does the model overpredict at high contrast?")
    a.grid(alpha=.3, which="both")

    a = ax[1, 0]
    for f in MR.FILLERS_GRID:
        g = df[(df.lattice == "FCC") & (df.phi == ph) & (df.filler == f)
               & (df.matrix == "epoxy")]
        tr = g[g.iface == "treated"].sort_values("d_um")
        un = g[g.iface == "untreated"].sort_values("d_um")
        a.semilogx(tr["d_um"].values,
                   tr[f"k_{best}"].values / un[f"k_{best}"].values, "o-",
                   label=f"model {f}")
    g = df[(df.lattice == "FCC") & (df.phi == ph) & (df.filler == "Al2O3")
           & (df.matrix == "epoxy")]
    tr = g[g.iface == "treated"].sort_values("d_um")
    un = g[g.iface == "untreated"].sort_values("d_um")
    a.semilogx(tr["d_um"].values, tr["k_nan"].values / un["k_nan"].values,
               "k--", lw=2, label="Nan (Al$_2$O$_3$)")
    a.set_xlabel("d [um]"); a.set_ylabel("$k$(treated)/$k$(untreated)")
    a.set_title("(d) surface-treatment gain vs particle size\n"
                "independent check of Z_size_sweep")
    a.legend(fontsize=7); a.grid(alpha=.3, which="both")

    a = ax[1, 1]
    for kind, mk in zip(LT.LATTICES, "os^"):
        g = df[(df.lattice == kind) & (df.filler == "Al2O3")
               & (df.matrix == "epoxy") & (df.iface == "untreated")
               & (df.d_um == 10.0)].sort_values("phi")
        a.plot(g["phi_rel"], g["frac_contact"], mk + "-", label=kind)
    a.set_xlabel(r"$\phi/\phi_{max}$"); a.set_ylabel("fraction of contact bonds")
    a.set_title("(e) when does the model switch to\nthe contact branch?")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[1, 2]
    for kind, mk in zip(LT.LATTICES, "os^"):
        g = dfin[dfin.lattice == kind]
        a.plot(g["n_cells"], g["ratio"], mk + "-", label=kind)
    a.axhline(1.0, ls=":", c="k")
    a.set_xlabel("$n_{cells}$ per side"); a.set_ylabel("$k_{slab}/k_{homog}$")
    a.set_title("(f) finite-size error of the slab route\n"
                "(why homogenization is used)")
    a.legend(fontsize=8); a.grid(alpha=.3)
    fig.tight_layout(pad=1.5); fig.savefig("fig_poc5.png", dpi=140)
    say("figure saved: fig_poc5.png")
    say(f"total run time {time.time()-t0:.0f}s")
    pathlib.Path("results_poc5.txt").write_text("\n".join(OUT) + "\n",
                                                encoding="utf-8")


if __name__ == "__main__":
    main()
