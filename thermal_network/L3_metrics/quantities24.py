# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L3_metrics/quantities24.py
#   sha256(src) : 1140afc0d6e868be
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 52 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical. In addition, 5 top-level definitions (and 1 module-level
# statement used only by them) that the figures and verification scripts of
# the paper do not use were omitted (they are listed by name in rb_manifest);
# the remaining definitions are unchanged.
# ---------------------------------------------------------------------------
# The 24 descriptors computed for one structure and material combination
# (schema version Q24_SCHEMA_VERSION); compute() returns them in one dict.
# theta_star = -log10 g_star is the critical conductance
# level of Sec. 3.6; the persistence quantities come from the barcode of Sec. 3.4
# and the sensitivity quantities from the dissipation share of Sec. 4.3.
from __future__ import annotations

import sys as _sys
import pathlib as _pathlib

_LIB = _pathlib.Path(__file__).resolve().parents[1]
if __package__ in (None, ""):
    for _d in ("L0_engine", "L1_materials", "L2_structures"):
        _p = str(_LIB / _d)
        if _p not in _sys.path:
            _sys.path.insert(0, _p)
if str(_LIB) not in _sys.path:
    _sys.path.insert(0, str(_LIB))

from dataclasses import dataclass, asdict

import numpy as np

from _constants import N_RANDOM_CONTROL
import materials_real as MR
import thermal_sheaf_filtration as tsf
from persistence_distance import bottleneck, TAU as _TAU

Q24_SCHEMA_VERSION = "1.1"

RC_LIT_SPAN_DECADES = 1.30


@dataclass
class Q24Config:
    gap_max_um: float | None = None
    eta_cut: float = 1.0
    slab_axis: int = 2
    slab_frac: float = 0.06
    box_L_um: float | None = None
    box_A_um2: float | None = None
    p_degrade: float = 0.10
    degrade_factor: float = 100.0
    n_mc: int = 20
    k_top: int = 20
    boost_factor: float = 10.0
    n_rand_control: int = N_RANDOM_CONTROL
    rc_halfwidth_dec: float = 0.65
    mc_seed: int = 1234
    keep_raw: bool = True
    skip_layer5: bool = False


QUANTITIES = [
    dict(no=1, layer=0, key="phi", name="packing fraction phi", unit="-", fig=False),
    dict(no=2, layer=0, key="z_contact", name="contact coordination number z", unit="edges", fig=False),
    dict(no=3, layer=0, key="edge_per_particle", name="edges per particle", unit="-", fig=False),
    dict(no=4, layer=0, key="n_isolated", name="isolated particles", unit="particles", fig=False),
    dict(no=5, layer=1, key="_hist", name="distribution of R by branch type", unit="K/W", fig=True),
    dict(no=6, layer=1, key="R_median", name="median R", unit="K/W", fig=False),
    dict(no=7, layer=1, key="R_decades", name="spread of R, log10(R95/R05)", unit="decades", fig=False),
    dict(no=8, layer=1, key="frac_contact", name="branch type counts, as a fraction", unit="-", fig=False),
    dict(no=9, layer=1, key="twoaK_over_h", name="median 2 a_K / h", unit="-", fig=False),
    dict(no=10, layer=1, key="_h", name="distribution of the gap h and the reduced radius", unit="um", fig=True),
    dict(no=11, layer=2, key="gain_interface", name="idealisation gain, interface -> 0", unit="x", fig=False),
    dict(no=12, layer=2, key="gain_matrix", name="idealisation gain, matrix -> infinity", unit="x", fig=False),
    dict(no=13, layer=2, key="gain_filler", name="idealisation gain, particle -> infinity", unit="x", fig=False),
    dict(no=14, layer=3, key="k_eff", name="G_eff and the effective thermal conductivity", unit="W/mK", fig=False),
    dict(no=15, layer=3, key="share_contact", name="dissipation share by branch type", unit="-", fig=False),
    dict(no=16, layer=3, key="cum_top5pct", name="cumulative dissipation curve", unit="-", fig=True),
    dict(no=17, layer=3, key="_flux_map", name="map of the top 5% of the dissipation", unit="-", fig=True),
    dict(no=18, layer=3, key="Reff_p95_over_p50", name="effective resistance to the sink", unit="K/W", fig=True),
    dict(no=19, layer=4, key="g_star", name="percolation threshold g*", unit="W/K", fig=False),
    dict(no=20, layer=4, key="g_allconn", name="all-connected threshold", unit="W/K", fig=False),
    dict(no=21, layer=4, key="n_bars", name="barcode / persistence diagram", unit="bars", fig=True),
    dict(no=22, layer=5, key="ret10_mean", name="retention under degradation (mean, sd, lower 5%)", unit="-", fig=True),
    dict(no=23, layer=5, key="gain_ratio", name="lift, top k edges over k random edges", unit="x", fig=True),
    dict(no=24, layer=5, key="d_B", name="d_B against the literature range of R_c", unit="decades", fig=True),
]

try:
    from q24_guidance import GUIDANCE as _GUIDANCE
    for _q in QUANTITIES:
        _q.update(_GUIDANCE.get(_q["key"], {}))
except ImportError:
    pass

_LAYER5_KEYS = ("ret10_mean", "ret10_std", "ret10_p05", "gain_top_pct",
                "gain_rand_pct", "gain_ratio", "d_B", "dlog_g_inf", "compress")


def stamp(cfg: Q24Config | None = None) -> dict:
    s = dict(tsf.engine_stamp())
    s["q24_schema_version"] = Q24_SCHEMA_VERSION
    s["persistence_tau"] = _TAU
    if cfg is not None:
        s["q24_cfg"] = asdict(cfg)
    return s


def _table(filler, matrix, treat, ideal=None, rc_scale=None):
    tab = MR.pair(filler, treat, filler, treat, matrix)
    if ideal == "interface":
        for k in tab.materials:
            tab.materials[k].R_s_matrix = 1e-30
            tab.materials[k].R_s_contact = 1e-30
    elif ideal == "matrix":
        tab.k_matrix = tab.k_matrix * 1e6
    elif ideal == "filler":
        for k in tab.materials:
            tab.materials[k].kappa = tab.materials[k].kappa * 1e6
    if rc_scale is not None:
        for k in tab.materials:
            tab.materials[k].R_s_matrix *= rc_scale
            tab.materials[k].R_s_contact *= rc_scale
    return tab


def _net(ids, typ, pos, dia, mol, tab, cfg: Q24Config, d_um):
    acfg = tsf.AnalysisConfig(
        gap_max_um=(cfg.gap_max_um if cfg.gap_max_um is not None else 0.5 * d_um),
        eta_cut=cfg.eta_cut, slab_axis=cfg.slab_axis,
        slab_frac=cfg.slab_frac, viz=False)
    G = tsf.build_thermal_network(ids, typ, pos, dia, acfg, table=tab, mols=mol)
    src, snk = tsf.identify_slabs(ids, pos, acfg)
    Gt, gv = tsf.attach_virtual_terminals(G, src, snk, acfg)
    return G, Gt, tsf.CellularSheaf(Gt), acfg


def _bars(Gt, acfg):
    h = tsf.ConductanceFiltration(Gt, acfg).sweep().dropna(subset=["sink_theta"])
    if "born_at_terminal" in h:
        h = h[~h["born_at_terminal"].astype(bool)]
    return h


def _diagram(h):
    return [(np.log10(max(b, 1e-300)), np.log10(max(s, 1e-300)))
            for b, s in zip(h["birth_theta"], h["sink_theta"])]


def compute(ids, types, positions, diameters, mols=None, *,
            filler=None, matrix=None, treat="untreated",
            table=None, cfg: Q24Config | None = None,
            d_um=None, label="", step=0, info=None) -> dict:
    cfg = cfg or Q24Config()
    pos = np.asarray(positions, float)
    dia = np.asarray(diameters, float)
    if d_um is None:
        d_um = float(np.median(dia))
    if table is None:
        if filler is None or matrix is None:
            raise ValueError("specify filler and matrix, or table")
        table = _table(filler, matrix, treat)

    km = table.k_matrix
    kappa = float(np.mean([m.kappa for m in table.materials.values()]))
    aK = float(np.mean([m.R_s_matrix for m in table.materials.values()])) * km * 1e6

    G, Gt, sh, acfg = _net(ids, types, pos, dia, mols, table, cfg, d_um)
    G0 = tsf.effective_conductance(sh, 0.0)
    flux = tsf.edge_flux_table(sh, 0.0)

    r = dict(step=step, label=label, filler=filler, matrix=matrix, treat=treat,
             alpha=kappa / km, a_K_um=aK, d_um=float(d_um))

    if info is not None and "phi" in info:
        r["phi"] = float(info["phi"]); r["phi_source"] = "generator"
    else:
        vol = float((np.pi / 6.0) * (dia ** 3).sum())
        box = pos.max(axis=0) - pos.min(axis=0) + dia.max()
        r["phi"] = float(vol / float(np.prod(box))); r["phi_source"] = "bbox_estimate"
    r["n_particle"] = int(G.number_of_nodes())
    r["n_edge"] = int(G.number_of_edges())
    r["edge_per_particle"] = r["n_edge"] / r["n_particle"]
    Rf = tsf.effective_resistance_to_sink(sh, 0.0)
    r["n_isolated"] = int(sum(1 for k, v in Rf.items()
                              if k >= 0 and not np.isfinite(v)))

    idx = {int(v): k for k, v in enumerate(np.asarray(ids))}
    Rall, brs, hs, rts = [], [], [], []
    for u, v, d in G.edges(data=True):
        Rall.append(1.0 / d["g"]); brs.append(d["branch"])
        if d["branch"] == "gap":
            hs.append(d["gap_um"])
        Ri, Rj = dia[idx[int(u)]] / 2, dia[idx[int(v)]] / 2
        rts.append(Ri * Rj / (Ri + Rj))
    Rall = np.array(Rall); brs = np.array(brs)
    hs = np.array(hs); rts = np.array(rts)
    q05, q50, q95 = np.percentile(Rall, [5, 50, 95])
    r["R_median"] = float(q50); r["R_p05"] = float(q05); r["R_p95"] = float(q95)
    r["R_decades"] = float(np.log10(q95 / q05))
    for b in ("gap", "contact", "bond"):
        r[f"n_{b}"] = int((brs == b).sum())
        r[f"frac_{b}"] = float((brs == b).mean())
    r["z_contact"] = 2.0 * float((brs == "contact").sum()) / r["n_particle"]
    r["h_median_um"] = float(np.median(hs)) if len(hs) else float("nan")
    r["h_p05_um"] = float(np.percentile(hs, 5)) if len(hs) else float("nan")
    r["h_p95_um"] = float(np.percentile(hs, 95)) if len(hs) else float("nan")
    r["rt_median_um"] = float(np.median(rts))
    r["twoaK_over_h"] = float(np.median(2 * aK / hs)) if len(hs) else float("nan")
    if cfg.keep_raw:
        r["_hist"] = {b: [float(x) for x in Rall[brs == b]]
                      for b in ("gap", "contact", "bond") if (brs == b).any()}
        r["_h"] = [float(x) for x in hs]
        r["_rt"] = [float(x) for x in rts]

    for ideal, tag in (("interface", "gain_interface"), ("matrix", "gain_matrix"),
                       ("filler", "gain_filler")):
        tab_i = (_table(filler, matrix, treat, ideal) if filler is not None
                 else _ideal_from_table(table, ideal))
        _, _, sh2, _ = _net(ids, types, pos, dia, mols, tab_i, cfg, d_um)
        r[tag] = float(tsf.effective_conductance(sh2, 0.0) / G0)

    r["G_eff"] = float(G0)
    if cfg.box_L_um is not None and cfg.box_A_um2 is not None:
        L = cfg.box_L_um * 1e-6
        A = cfg.box_A_um2 * 1e-12
        r["box_source"] = "explicit"
    else:
        ax = cfg.slab_axis
        other = [i for i in range(3) if i != ax]
        L = float(pos[:, ax].max() - pos[:, ax].min()) * 1e-6
        A = float(np.prod([pos[:, i].max() - pos[:, i].min() for i in other])) * 1e-12
        r["box_source"] = "bbox"
    r["k_eff"] = float(G0 * L / A)
    r["slab_L_um"] = L * 1e6; r["slab_A_um2"] = A * 1e12
    by = flux.groupby("branch")["P_share"].sum()
    for b in ("gap", "contact", "bond"):
        r[f"share_{b}"] = float(by.get(b, 0.0))
    ps = np.sort(flux["P_share"].values)[::-1]
    cum = np.cumsum(ps) / ps.sum()
    n = len(ps)
    r["cum_top1pct"] = float(cum[max(int(0.01 * n) - 1, 0)])
    r["cum_top5pct"] = float(cum[max(int(0.05 * n) - 1, 0)])
    r["n_for_half"] = int(np.searchsorted(cum, 0.5) + 1)
    r["frac_for_half"] = r["n_for_half"] / n
    v = np.array([x for k, x in Rf.items() if k >= 0 and np.isfinite(x)])
    r["Reff_median"] = float(np.median(v)); r["Reff_max"] = float(v.max())
    r["Reff_p95_over_p50"] = float(np.percentile(v, 95) / np.median(v))
    if cfg.keep_raw:
        r["_cum"] = [float(x) for x in cum]

    gs, crit, ga = tsf.percolation_threshold(Gt)
    r["g_star"] = float(gs)
    r["theta_star"] = float(-np.log10(gs))
    r["g_allconn"] = float(ga) if ga is not None else float("nan")
    r["theta_allconn"] = float(-np.log10(ga)) if ga is not None else float("nan")
    h = _bars(Gt, acfg)
    pv = h["persistence_to_sink"].values
    r["n_bars"] = int(len(pv))
    r["bar_max"] = float(pv.max()) if len(pv) else 0.0
    r["bar_zero_frac"] = float((pv < 1e-12).mean()) if len(pv) else float("nan")
    if cfg.keep_raw:
        r["_bars"] = [[float(-np.log10(b)), float(-np.log10(s))]
                      for b, s in zip(h["birth_theta"], h["sink_theta"])]

    if cfg.skip_layer5:
        for k in _LAYER5_KEYS:
            r[k] = float("nan")
        r.update(stamp(cfg))
        return r

    rng = np.random.default_rng(cfg.mc_seed + step)
    real = [(u, v) for u, v in Gt.edges() if u >= 0 and v >= 0]
    n_sel = max(1, int(round(cfg.p_degrade * len(real))))
    ret = np.empty(cfg.n_mc)
    for m in range(cfg.n_mc):
        pick = [real[i] for i in rng.choice(len(real), n_sel, replace=False)]
        saved = []
        for u, v in pick:
            saved.append((u, v, Gt[u][v]["g"])); Gt[u][v]["g"] /= cfg.degrade_factor
        ret[m] = tsf.effective_conductance(sh, 0.0) / G0
        for u, v, g in saved:
            Gt[u][v]["g"] = g
    r["ret10_mean"] = float(ret.mean()); r["ret10_std"] = float(ret.std(ddof=1))
    r["ret10_p05"] = float(np.percentile(ret, 5))

    def boost(pairs):
        saved = []
        for u, v in pairs:
            saved.append((u, v, Gt[u][v]["g"])); Gt[u][v]["g"] *= cfg.boost_factor
        g = tsf.effective_conductance(sh, 0.0)
        for u, v, x in saved:
            Gt[u][v]["g"] = x
        return g / G0 - 1.0

    top = [(int(x.u), int(x.v)) for x in flux.head(cfg.k_top).itertuples()]
    gt = boost(top)
    gr = float(np.mean([boost([real[i] for i in
                               rng.choice(len(real), cfg.k_top, replace=False)])
                        for _ in range(cfg.n_rand_control)]))
    r["gain_top_pct"] = gt * 100
    r["gain_rand_pct"] = gr * 100
    r["gain_ratio"] = gt / gr if gr > 0 else float("nan")

    D0 = _diagram(h)
    g_a = {frozenset((u, v)): d["g"] for u, v, d in G.edges(data=True)}
    dbs, bnds = [], []
    for sc in (10 ** -cfg.rc_halfwidth_dec, 10 ** cfg.rc_halfwidth_dec):
        tab2 = (_table(filler, matrix, treat, None, sc) if filler is not None
                else _scale_rc_from_table(table, sc))
        G2, Gt2, _, acfg2 = _net(ids, types, pos, dia, mols, tab2, cfg, d_um)
        dbs.append(bottleneck(D0, _diagram(_bars(Gt2, acfg2))))
        g_b = {frozenset((u, v)): d["g"] for u, v, d in G2.edges(data=True)}
        com = set(g_a) & set(g_b)
        bnds.append(float(np.abs(np.log10([g_b[k] / g_a[k] for k in com])).max()))
    r["d_B"] = float(max(dbs))
    r["dlog_g_inf"] = float(max(bnds))
    r["compress"] = RC_LIT_SPAN_DECADES / r["dlog_g_inf"]

    r.update(stamp(cfg))
    return r


def _ideal_from_table(table, ideal):
    import copy
    t = copy.deepcopy(table)
    if ideal == "interface":
        for k in t.materials:
            t.materials[k].R_s_matrix = 1e-30
            t.materials[k].R_s_contact = 1e-30
    elif ideal == "matrix":
        t.k_matrix = t.k_matrix * 1e6
    elif ideal == "filler":
        for k in t.materials:
            t.materials[k].kappa = t.materials[k].kappa * 1e6
    return t


def _scale_rc_from_table(table, sc):
    import copy
    t = copy.deepcopy(table)
    for k in t.materials:
        t.materials[k].R_s_matrix *= sc
        t.materials[k].R_s_contact *= sc
    return t

