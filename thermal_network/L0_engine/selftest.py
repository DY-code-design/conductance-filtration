# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L0_engine/selftest.py
#   sha256(src) : db474df566b98b29
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 76 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# Run-time self-test of the package: presence of all public symbols,
# closed-form identities of the branch model and of the finite-change formula,
# consistency of the sparse and dense solvers, and the model-variant switch.
# Run it before anything else; every check must print PASS.
import sys as _sys, pathlib as _pathlib
_LIB = _pathlib.Path(__file__).resolve().parent.parent
for _d in ("L0_engine", "L1_materials", "L2_structures"):
    _p = str(_LIB / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import sys
import time
import numpy as np

FAILED = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {name}" + (f"  ({detail})" if detail else ""))
    if not cond:
        FAILED.append(name)


try:
    import thermal_sheaf_filtration as tsf
except Exception as e:
    print(f"[FAIL] import thermal_sheaf_filtration: {e}")
    sys.exit(1)

REQUIRED = [
    "AnalysisConfig", "Material", "MaterialTable", "pair_conductance",
    "build_thermal_network", "identify_slabs", "attach_virtual_terminals",
    "CellularSheaf", "ConductanceFiltration",
    "solve_dirichlet", "effective_conductance", "edge_flux_table",
    "effective_resistance_to_sink",
    "shorted_edge_resistance",
    "geometric_contact_radius", "spot_radius_um", "predict_finite_change",
]

REQUIRED_INTERNAL = [
    "load_last_frame_with_types", "heat_proxy_map", "betti_numbers", "lambda_min_grounded",
    "grounded_slow_mode", "node_sink_theta", "risk_table",
    "BlockSheaf", "build_jet_sheaf", "kron_reduction",
    "cluster_conduction_tensor", "phase2_component_table",
    "component_sweep", "cluster_census", "plot_cluster_census",
    "z_visibility_sweep", "plot_z_visibility",
    "jet_weak_modes", "jet_weak_modes_all",
    "run_all",
    "run_all", "plot_summary", "main",
]
_HAS_PHASE2 = all(hasattr(tsf, n) for n in ("kron_reduction", "jet_weak_modes", "component_sweep",
                                             "cluster_census", "heat_proxy_map"))
missing = [s for s in REQUIRED if not hasattr(tsf, s)]
_present_internal = [s for s in REQUIRED_INTERNAL if hasattr(tsf, s)]
check("completeness: required symbols", not missing,
      f"missing: {missing}" if missing else f"{len(REQUIRED)} symbols"
      + (f" (+{len(_present_internal)} internal)" if _present_internal else ""))
if missing:
    print("\n*** the file may be truncated; obtain a fresh copy ***")
    sys.exit(1)

import networkx as nx
from thermal_sheaf_filtration import Material, MaterialTable, AnalysisConfig

cfg = AnalysisConfig(gap_max_um=10.0, slab_frac=0.08)
tab30 = MaterialTable(materials={1: Material("f", 30.0)}, k_matrix=0.2)

r_t = 5.0
a_geo_m = tsf.geometric_contact_radius(10.0, 10.0, 0.1) * 1e-6
g, at = tsf.pair_conductance(-0.1, 10.0, 10.0, 1, 1, tab30, cfg)
_a_um = a_geo_m * 1e6
_cm0 = getattr(cfg, "constriction_model", "current")
if _cm0 == "current":
    _Rcon_exp = 2.0 / (2.0 * 30.0 * a_geo_m)
elif _cm0 == "holm":
    _Rcon_exp = 2.0 / (4.0 * 30.0 * a_geo_m)
else:
    _cut0 = getattr(cfg, "constriction_cut", "R_min")
    _Rc0 = {"R_min": 10.0, "R_each": 10.0, "r_tilde": r_t}[_cut0]
    _Rcon_exp = 2.0 * np.arctan(_Rc0 / _a_um) / (2.0 * np.pi * 30.0 * a_geo_m)
check("T1 constriction-limited g_spot = 1/R_con (written out from the variant formula)",
      abs(at["g_spot"] - 1.0 / _Rcon_exp) / (1.0 / _Rcon_exp) < 1e-12,
      f'g_spot={at["g_spot"]:.6e}, expected {1.0/_Rcon_exp:.6e}, '
      f'variant {_cm0}, a_geo={_a_um:.6f} um')
_cfgL = AnalysisConfig(gap_max_um=cfg.gap_max_um, eta_cut=cfg.eta_cut,
                       contact_radius="hertz", contact_composition="exclusive",
                       prox_flux_order="series_of_parallels",
                       bulk_transport=False,
                       constriction_model="current", constriction_cut="R_min",
                       bulk_excludes_constriction=False)
_aH_m = np.sqrt(0.1 * r_t) * 1e-6
_gL, _ = tsf.pair_conductance(-0.1, 10.0, 10.0, 1, 1, tab30, _cfgL)
check("T1b legacy setting also gives g = kappa a_Hertz (backward compatibility)",
      abs(_gL - 30.0 * _aH_m) / (30.0 * _aH_m) < 1e-12)


def _R_bulk_prox_expected(kappa, R_um, cfg):
    l = R_um
    if getattr(cfg, "bulk_excludes_constriction", False):
        cut = getattr(cfg, "constriction_cut", "R_min")
        rt = R_um / 2.0
        c = {"R_min": R_um, "R_each": R_um, "r_tilde": rt}[cut]
        l = max(R_um - c, 0.0)
    return 2.0 * (l * 1e-6 / (kappa * np.pi * (R_um * 1e-6) ** 2))


g, _at3 = tsf.pair_conductance(1.0, 10.0, 10.0, 1, 1, tab30, cfg)
_hb3 = r_t / ((30.0 / 0.2) ** 2 - 1.0)
_gap_only = 2 * np.pi * 0.2 * (r_t * 1e-6) * np.log1p(r_t / (1.0 + _hb3))
_Rb3 = _R_bulk_prox_expected(30.0, 10.0, cfg)
exp = 1.0 / (1.0 / _gap_only + _Rb3)
check("T3 near-contact branch closed form (h_BOB shift + bulk term)", abs(g - exp) / exp < 1e-12,
      f"h_BOB={_hb3:.6f} um, R_bulk={_Rb3:.1f} K/W, "
      f"bulk term scales g by x{exp/_gap_only:.4f}")

from scipy.integrate import quad as _quad
_cfg_v1 = AnalysisConfig(bob_saturation=False, bulk_transport=False, gap_max_um=10.0)
_tab_v1 = MaterialTable(materials={1: Material("f", 30.0, R_s_matrix=1.11e-7)}, k_matrix=0.2)
_worst_v1 = 0.0
for _h_um in (0.001, 0.01, 0.1, 1.0, 2.0):
    _g_eng, _ = tsf.pair_conductance(_h_um, 10.0, 10.0, 1, 1, _tab_v1, _cfg_v1)
    _rt_m, _h_m, _k = 5.0e-6, _h_um * 1e-6, 0.2
    _Rs2 = 2 * 1.11e-7
    _num, _ = _quad(lambda r: 2 * np.pi * r / (_Rs2 + (_h_m + r * r / (2 * _rt_m)) / _k),
                    0.0, np.sqrt(2.0) * _rt_m, epsabs=0, epsrel=1e-13, limit=500)
    _worst_v1 = max(_worst_v1, abs(_num - _g_eng) / _g_eng)
check("V1 closed form of the near-contact branch vs numerical annular integration (5 values of h)", _worst_v1 < 1e-12,
      f"max rel {_worst_v1:.2e}")

tab_ip = MaterialTable(materials={1: Material("f", 30.0, t_ip_um=1.0,
                                              k_ip=2.0)}, k_matrix=0.2)
g_ip, _ = tsf.pair_conductance(1.0, 10.0, 10.0, 1, 1, tab_ip, cfg)
check("T4 effect of the interphase (k_eff ratio x10 dominates)",
      8.0 < g_ip / g < 10.0, f"g_ip/g = {g_ip/g:.4f}")

alpha = 30.0 / 0.2
g_sat_exp = 4 * np.pi * 0.2 * (r_t * 1e-6) * np.log(alpha)
tab30_nR = MaterialTable(materials={1: Material("f", 30.0, R_s_matrix=0.0,
                                                R_s_contact=0.0)},
                         k_matrix=0.2)
_hbT = r_t / (alpha ** 2 - 1.0)
_g_an = 2 * np.pi * 0.2 * (r_t * 1e-6) * np.log1p(r_t / _hbT)
check("T10 (V2) at h=0 and a_K=0 the shift equals the legacy BOB cap (analytic identity)",
      abs(_g_an - g_sat_exp) / g_sat_exp < 1e-12,
      f"ratio={_g_an/g_sat_exp:.12f}")
g_deep, at_deep = tsf.pair_conductance(1e-12, 10.0, 10.0, 1, 1, tab30_nR, cfg)
_g_exp_clamp = 1.0 / (1.0 / (2 * np.pi * 0.2 * (r_t * 1e-6) * np.log1p(
    r_t / (cfg.t_min_um + _hbT))) + _R_bulk_prox_expected(30.0, 10.0, cfg))
check("T10b implementation matches the closed form including the t_min clamp",
      abs(g_deep - _g_exp_clamp) / _g_exp_clamp < 1e-12,
      f"g={g_deep:.6e} / ratio to cap {g_deep/g_sat_exp:.6f}"
      f"(lowered by the t_min={cfg.t_min_um} um clamp)")

cfg_nosat = AnalysisConfig(gap_max_um=10.0, slab_frac=0.08, bob_saturation=False)
g_on, at_on = tsf.pair_conductance(1.0, 10.0, 10.0, 1, 1, tab30, cfg)
g_off, _ = tsf.pair_conductance(1.0, 10.0, 10.0, 1, 1, tab30, cfg_nosat)
tab_hi = MaterialTable(materials={1: Material("f", 20000.0)}, k_matrix=0.2)
g_on_hi, _ = tsf.pair_conductance(1.0, 10.0, 10.0, 1, 1, tab_hi, cfg)
g_off_hi, _ = tsf.pair_conductance(1.0, 10.0, 10.0, 1, 1, tab_hi, cfg_nosat)
check("T11 as alpha->inf, bob_saturation on/off agree (within 1e-9 at alpha=1e5)",
      abs(g_on_hi - g_off_hi) / g_off_hi < 1e-9,
      f"at alpha=150 the difference is {abs(g_on-g_off)/g_off:.3e}, "
      f"at alpha=1e5 it is {abs(g_on_hi-g_off_hi)/g_off_hi:.3e}")

A_th_shift = (2 * np.pi * (r_t * 1e-6) * (1.0 * 1e-6)
              * np.log1p(r_t / (1.0 + _hb3)))
check("T11b A_th is consistent with the shifted g (bulk term excluded)",
      abs(at_on["A_th_m2"] - A_th_shift) / A_th_shift < 1e-14,
      f"rel.diff={abs(at_on['A_th_m2']-A_th_shift)/A_th_shift:.2e}")

tab_lowk = MaterialTable(materials={1: Material("f", 0.5)}, k_matrix=2.0)
_, at_low = tsf.pair_conductance(1e-7, 10.0, 10.0, 1, 1, tab_lowk, cfg)
check("T12 no saturation applied for alpha<=1", not at_low["bob_saturated"],
      f"alpha={at_low['bob_alpha']:.3f}")

_G = nx.Graph()
for _p in (1, 2, 3, 4):
    _G.add_node(_p, atom_type=1)
for _u, _v, _g in ((1, 2, 2.0), (2, 3, 1.0), (3, 4, 3.0), (2, 4, 0.7)):
    _G.add_edge(_u, _v, g=_g, r_birth=0.0)
_cfg13 = AnalysisConfig(gap_max_um=10.0, slab_frac=0.08, g_virtual_factor=1e3)
_Gt, _ = tsf.attach_virtual_terminals(_G, {1}, {4}, _cfg13)
_sh = tsf.CellularSheaf(_Gt)
_G0 = tsf.effective_conductance(_sh, 0.0)
_flux = tsf.edge_flux_table(_sh, 0.0)
_rt = tsf.shorted_edge_resistance(_sh, 0.0).set_index(["u", "v"])
_errs = []
for _, _r in _flux.iterrows():
    _u, _v = int(_r["u"]), int(_r["v"])
    _key = (_u, _v) if (_u, _v) in _rt.index else (_v, _u)
    for _c in (1.0, 10.0, 100.0):
        _D = _c * _Gt[_u][_v]["g"]
        _pred = tsf.predict_finite_change(_G0, _r["dT"],
                                          _rt.loc[_key, "R_tilde"], _D)
        _g_save = _Gt[_u][_v]["g"]
        _Gt[_u][_v]["g"] = _g_save + _D
        _true = tsf.effective_conductance(_sh, 0.0)
        _Gt[_u][_v]["g"] = _g_save
        _errs.append(abs(_pred - _true) / _true)
check("T13 Corollary 2.4 closed form for finite changes (uses shorted R~)", max(_errs) < 1e-10,
      f"max rel.err={max(_errs):.2e} over {len(_errs)} cases")

_rt_d = tsf.shorted_edge_resistance(_sh, 0.0, method="dense")
_rt_s = tsf.shorted_edge_resistance(_sh, 0.0, method="sparse")
_m = _rt_d.merge(_rt_s, on=["u", "v"], suffixes=("_d", "_s"))
_rel = ((_m["R_tilde_d"] - _m["R_tilde_s"]).abs()
        / _m["R_tilde_d"].abs()).max()
check("T22 sparse-solver R~ = dense R~", len(_m) == len(_rt_d) and _rel < 1e-8,
      f"checked on {len(_m)} edges, max relative difference {_rel:.2e}")

_sel = [(int(r.u), int(r.v)) for r in _rt_d.head(5).itertuples()]
_rt_k = tsf.shorted_edge_resistance(_sh, 0.0, edges=_sel)
_mk = _rt_k.merge(_rt_d, on=["u", "v"], suffixes=("_k", "_d"))
_relk = ((_mk["R_tilde_k"] - _mk["R_tilde_d"]).abs()
         / _mk["R_tilde_d"].abs()).max()
check("T22b partial R~ computation (edges=) = all-edge computation",
      len(_mk) == len(_sel) and _relk < 1e-8,
      f"{len(_mk)} edges, max relative difference {_relk:.2e}")

_rt_dense_small = tsf.shorted_edge_resistance(_sh, 0.0, max_n=3, method="dense")
_rt_sparse_small = tsf.shorted_edge_resistance(_sh, 0.0, max_n=3, method="auto")
check("T22c legacy path returns an empty table for n>max_n; the sparse version returns values",
      len(_rt_dense_small) == 0 and len(_rt_sparse_small) == len(_rt_d),
      f"dense={len(_rt_dense_small)} rows / auto(sparse)={len(_rt_sparse_small)} rows")


import materials_real as MR

_UNI = dict(contact_radius="geometric", contact_composition="parallel",
            prox_flux_order="parallel_of_series")

_ag = [tsf.geometric_contact_radius(10.0, 10.0, _d) / np.sqrt(2 * 5.0 * _d)
       for _d in (1e-4, 1e-3, 1e-2)]
check("T14u geometric intersection circle -> sqrt(2 r~ delta) (small-delta limit)",
      max(abs(_x - 1.0) for _x in _ag) < 1e-3,
      f"a_geo/sqrt(2r~delta) = {['%.6f' % x for x in _ag]}")

_r2 = tsf.geometric_contact_radius(10.0, 10.0, 1e-3) / np.sqrt(1e-3 * 5.0)
check("T15u a_geo/a_Hertz = sqrt(2) (twice the area)", abs(_r2 / np.sqrt(2.0) - 1) < 1e-3,
      f"ratio={_r2:.6f} vs sqrt(2)={np.sqrt(2):.6f}")

_aK2 = 2 * MR.TREATMENT["untreated"]["R_s_matrix"] * 0.2 * 1e6
_rho = [tsf.spot_radius_um(_h, 5.0, _aK2) for _h in (0.0, 1e-3, 1e-2, 0.1, 1.0)]
check("T16u rho(0)=sqrt(4 r~ a_K)>0 and monotonically increasing in h",
      abs(_rho[0] / np.sqrt(4 * 5.0 * _aK2 / 2) - 1) < 1e-12
      and all(np.diff(_rho) > 0),
      f"rho(0)={_rho[0]:.5f} um, monotone={all(np.diff(_rho) > 0)}")

_tabu = MR.single("Al2O3", "untreated", "epoxy")
_cu = AnalysisConfig(gap_max_um=100.0, eta_cut=1.0, t_min_um=1e-9,
                     delta_min_um=1e-12, **_UNI)
_gp, _ = tsf.pair_conductance(1e-9, 10.0, 10.0, 1, 1, _tabu, _cu)
_gc, _ = tsf.pair_conductance(-1e-9, 10.0, 10.0, 1, 1, _tabu, _cu)
check("T17u (V3) the switch is continuous (delta->0+ and h->0+ agree)",
      abs(_gc / _gp - 1) < 1e-6,
      f"g(delta->0+)/g(h->0+) = {_gc/_gp:.9f}")

_xs = np.concatenate([np.logspace(np.log10(2.5), -9, 220),
                      -np.logspace(-9, np.log10(0.5), 220)])
_gs = []
for _x in _xs:
    _g, _ = tsf.pair_conductance(float(_x), 10.0, 10.0, 1, 1, _tabu, _cu)
    _gs.append(np.nan if _g is None else _g)
_gs = np.array(_gs)
_ok = np.isfinite(_gs)
_dg = np.diff(_gs[_ok])
_viol = int((_dg < -1e-18 * np.max(_gs[_ok])).sum())
check("T18u monotonicity (g is non-decreasing from approach to compression)", _viol == 0,
      f"violations {_viol}/{len(_dg)} points, g range "
      f"{_gs[_ok].min():.3e}-{_gs[_ok].max():.3e} W/K")

_cl = AnalysisConfig(gap_max_um=100.0, eta_cut=1.0, t_min_um=1e-9,
                     delta_min_um=1e-12, contact_radius="hertz",
                     contact_composition="exclusive",
                     prox_flux_order="series_of_parallels")
_gl = []
for _x in _xs:
    _g, _ = tsf.pair_conductance(float(_x), 10.0, 10.0, 1, 1, _tabu, _cl)
    _gl.append(np.nan if _g is None else _g)
_gl = np.array(_gl)
_okl = np.isfinite(_gl)
_violl = int((np.diff(_gl[_okl]) < -1e-18 * np.max(_gl[_okl])).sum())
check("T19u monotonicity fails under the legacy setting", _violl > 0,
      f"violations under the legacy setting {_violl}/{int(_okl.sum())-1} points")

G = nx.Graph()
for pid in (1, 2, 3, 10, 11):
    G.add_node(pid, atom_type=1)
G.add_edge(1, 2, g=2.0, r_birth=0.0)
G.add_edge(2, 3, g=1.0, r_birth=0.0)
G.add_edge(10, 11, g=0.5, r_birth=0.0)
Gt, _ = tsf.attach_virtual_terminals(G, {1}, {3}, cfg)
sheaf = tsf.CellularSheaf(Gt)
ge = tsf.effective_conductance(sheaf, 0.0)
check("P1 series G_eff=2/3 (with a floating component)", abs(ge - 2 / 3) < 1e-3,
      f"G_eff={ge:.5f}")

def _phase2_checks():
    G = nx.Graph()
    grid = {}
    for i in range(3):
        for j in range(3):
            n = 3 * i + j
            grid[(i, j)] = n
            G.add_node(n, atom_type=1, pos=np.array([float(i), 0.0, float(j)]))
    for (i, j), n in grid.items():
        for di, dj in ((1, 0), (0, 1)):
            if (i + di, j + dj) in grid:
                G.add_edge(n, grid[(i + di, j + dj)], g=1.0)
    G.add_node(100, atom_type=1, pos=np.array([-1.0, 0.0, 1.0]))
    G.add_node(101, atom_type=1, pos=np.array([3.0, 0.0, 1.0]))
    G.add_edge(100, grid[(0, 1)], g=1.0)
    G.add_edge(101, grid[(2, 1)], g=1.0)
    Lam, cut, _ = tsf.kron_reduction(G, set(grid.values()), 0.0)
    term = [c[1] for c in cut]
    i_s, i_t = term.index(100), term.index(101)
    m = Lam.shape[0]
    free = [k for k in range(m) if k not in (i_s, i_t)]
    u = np.zeros(m); u[i_s] = 1.0
    if free:
        u[free] = np.linalg.solve(Lam[np.ix_(free, free)],
                                  -Lam[np.ix_(free, [i_s])] @ np.array([1.0]))
    Q_red = -float(Lam[i_t, :] @ u)
    nodes = sorted(G.nodes()); idx = {v: k for k, v in enumerate(nodes)}
    N = len(nodes); W = np.zeros((N, N))
    for a, b, d in G.edges(data=True):
        W[idx[a], idx[b]] = W[idx[b], idx[a]] = d["g"]
    L = np.diag(W.sum(1)) - W
    bd = [idx[100], idx[101]]; Tb = np.array([1.0, 0.0])
    it = [k for k in range(N) if k not in bd]
    Ti = np.linalg.solve(L[np.ix_(it, it)], -L[np.ix_(it, bd)] @ Tb)
    T = np.zeros(N); T[it] = Ti; T[bd] = Tb
    Q_full = -float(L[idx[101], :] @ T)
    check("V3b Kron reduction preserves G_eff", abs(Q_red - Q_full) / Q_full < 1e-10,
          f"{Q_red:.6e} vs {Q_full:.6e}")

    jw = tsf.jet_weak_modes(G, nodes=set(grid.values()), eps=1e-2, n_modes=2)
    check("J1 jet kernel dimension=4 (connected cluster)", jw["kernel_dim"] == 4,
          f"kdim={jw['kernel_dim']}")
    jw_th = tsf.jet_weak_modes(G, nodes=set(grid.values()), eps=1e-2,
                               n_modes=1, theta=2.0)
    check("J-theta jet theta filter (all edges excluded -> kernel=4n)",
          jw_th["kernel_dim"] == 4 * 9, f"kdim={jw_th['kernel_dim']}")

    pos, dia, types = [], [], []
    for z in (1.0, 2.9, 4.8, 6.7, 8.6):
        pos.append([2.0, 2.0, z]); dia.append(2.0); types.append(1)
    for x, z in ((8.0, 5.0), (9.9, 5.0), (8.95, 6.6)):
        pos.append([x, 2.0, z]); dia.append(2.0); types.append(1)
    pos, dia, types = np.array(pos), np.array(dia), np.array(types)
    ids = np.arange(1, len(dia) + 1)
    Gs = tsf.build_thermal_network(ids, types, pos, dia, cfg, tab30)
    src, snk = tsf.identify_slabs(ids, pos, cfg)
    Gts, _ = tsf.attach_virtual_terminals(Gs, src, snk, cfg)
    heat = tsf.heat_proxy_map(ids, dia, cfg)
    gv = np.array([d["g"] for _, _, d in Gs.edges(data=True)])
    thetas = np.logspace(np.log10(gv.min() * 0.5), np.log10(gv.max() * 1.01), 6)[::-1]
    sw = tsf.component_sweep(Gts, thetas, heat=heat, n_min=3, jet=True)
    need_cols = {"theta", "rep", "n_nodes", "k_zz", "A", "rank_K", "S",
                 "G_boundary", "stagnation", "has_sink", "lam_weak"}
    check("SW sweep columns present", need_cols <= set(sw.columns),
          f"missing: {need_cols - set(sw.columns)}")
    cen = tsf.cluster_census(sw)
    low = cen.iloc[-1]
    check("SW census: chain=good 1 / isolated triad=z-blocked 1",
          int(low.n_good) == 1 and int(low.n_isolated_z) == 1,
          f"n_good={int(low.n_good)}, n_isolated_z={int(low.n_isolated_z)}")


if _HAS_PHASE2:
    _phase2_checks()
else:
    print("[SKIP] Phase-2 checks (Kron reduction / jet sheaf / component sweep): not part of this record")

if "--perf" in sys.argv:
    print("\n--- perf: N=2500 ---")
    rng = np.random.default_rng(0)
    nx_, ny_, nz_ = 17, 9, 17
    xs, ys, zs = (np.arange(n) * 2.05 + 1.2 for n in (nx_, ny_, nz_))
    P = np.array([[x, y, z] for x in xs for y in ys for z in zs])
    P += rng.uniform(-0.16, 0.16, P.shape)
    D = rng.uniform(1.9, 2.1, len(P))
    T_ = np.ones(len(P), int)
    I_ = np.arange(1, len(P) + 1)
    t0 = time.time()
    Gp = tsf.build_thermal_network(I_, T_, P, D, cfg, tab30)
    t1 = time.time()
    print(f"build_thermal_network: N={len(P)} E={Gp.number_of_edges()} "
          f"{t1-t0:.1f}s")
    srcp, snkp = tsf.identify_slabs(I_, P, cfg)
    Gtp, gvv = tsf.attach_virtual_terminals(Gp, srcp, snkp, cfg)
    shp = tsf.CellularSheaf(Gtp)
    t2 = time.time()
    hist = tsf.ConductanceFiltration(Gtp).sweep()
    t3 = time.time()
    print(f"UF sweep: {t3-t2:.1f}s ({len(hist)} records)")
    gvals = np.array([d['g'] for _, _, d in Gp.edges(data=True)])
    ths = np.logspace(np.log10(gvals.min() * 0.5),
                      np.log10(gvals.max()), 8)[::-1]
    t4 = time.time()
    for th in ths:
        tsf.effective_conductance(shp, th)
    t5 = time.time()
    print(f"G_eff x8 theta: {t5-t4:.1f}s")
    if _HAS_PHASE2:
        t6 = time.time()
        swp = tsf.component_sweep(Gtp, ths, heat=tsf.heat_proxy_map(I_, D, cfg),
                                  n_min=3, jet=False)
        t7 = time.time()
        print(f"component_sweep(jet=False) x8 theta: {t7-t6:.1f}s ({len(swp)} rows)")
    t8 = time.time()
    R = tsf.effective_resistance_to_sink(shp, gvals.min() * 0.5)
    t9 = time.time()
    print(f"effective_resistance_to_sink (N solves): {t9-t8:.1f}s")
    print(f"perf total: {t9-t0:.1f}s")

import ensemble as _EN
_st = _EN.build("phi_060", d_um=20.0, seed=0)
_gs, _ce, _ga = tsf.percolation_threshold(_st["Gt"])
_h = tsf.ConductanceFiltration(_st["Gt"], _st["cfg"]).sweep()
_mn = float(_h["sink_theta"].dropna().min())
check("T20 all-connected threshold = min sink_theta of the sweep",
      abs(_ga - _mn) / _mn < 1e-12, f"(relative difference {abs(_ga-_mn)/_mn:.2e})")
check("T20b percolation threshold > all-connected threshold (in g)",
      _gs > _ga, f"(percolation {_gs:.4e} > all-connected {_ga:.4e})")

_L = {}
for _fac in (10.0, 1000.0):
    _c2 = tsf.AnalysisConfig(gap_max_um=10.0, eta_cut=1.0, slab_axis=2,
                             slab_frac=0.06, viz=False, g_virtual_factor=_fac)
    _G2 = tsf.build_thermal_network(_st["ids"], _st["typ"], _st["pos"], _st["dia"],
                                   _c2, table=_st["table"], mols=_st["mol"])
    _a2, _b2 = tsf.identify_slabs(_st["ids"], _st["pos"], _c2)
    _Gt2, _ = tsf.attach_virtual_terminals(_G2, _a2, _b2, _c2)
    _h2 = tsf.ConductanceFiltration(_Gt2, _c2).sweep().dropna(subset=["sink_theta"])
    _p = _h2["persistence_to_sink"].values
    _t = _h2["born_at_terminal"].values.astype(bool)
    _L[_fac] = (float(_p.max()), float(_p[~_t].max()))
check("T21 the longest bar moves by 2 decades with the terminal factor",
      abs((_L[1000.0][0] - _L[10.0][0]) - 2.0) < 1e-9,
      f"(x10 {_L[10.0][0]:.4f} -> x1000 {_L[1000.0][0]:.4f})")
check("T21b excluding terminal-born bars, invariant to the terminal factor",
      abs(_L[1000.0][1] - _L[10.0][1]) < 1e-9,
      f"(both {_L[1000.0][1]:.4f} decades)")
_gs2, _, _ = tsf.percolation_threshold(_st["Gt"])
_h3 = tsf.ConductanceFiltration(_st["Gt"], _st["cfg"]).sweep().dropna(subset=["sink_theta"])
_gv = max(d["g"] for u, v, d in _st["Gt"].edges(data=True)
          if u < 0 or v < 0)
check("T21c pers_max = theta*(percolation) - theta_virtual",
      abs(_h3["persistence_to_sink"].max()
          - (-np.log10(_gs2) + np.log10(_gv))) < 1e-12,
      f"(pers_max {_h3['persistence_to_sink'].max():.6f})")

_kap, _a = 150.0, 2.0e-6
_R_holm_one  = 1.0 / (4.0 * _kap * _a)
_R_holm_both = 1.0 / (2.0 * _kap * _a)
check("T23 Holm one-sided x2 = two-sided total (self-consistency of the external formula)",
      abs(2 * _R_holm_one / _R_holm_both - 1) < 1e-15)

def _rcon(model, ri=10.0, rj=10.0, a_um=2.0, cut="R_min"):
    import numpy as _np
    a_m = a_um * 1e-6
    if model == "current":
        return 2.0 * (1.0 / (2.0 * _kap * a_m))
    if model == "holm":
        return 2.0 * (1.0 / (4.0 * _kap * a_m))
    Rc = {"R_min": min(ri, rj), "R_each": ri, "r_tilde": ri * rj / (ri + rj)}[cut]
    return 2.0 * (_np.arctan(Rc / a_um) / (2.0 * _np.pi * _kap * a_m))

check("T23a current is 2x the Holm two-sided total",
      abs(_rcon("current") / _R_holm_both - 2.0) < 1e-12,
      f"(ratio {_rcon('current')/_R_holm_both:.6f})")
check("T23b holm matches the Holm two-sided total",
      abs(_rcon("holm") / _R_holm_both - 1.0) < 1e-12)
check("T23c holm_finite converges to holm as R_cut/a->inf",
      abs(_rcon("holm_finite", ri=1e9, rj=1e9) / _rcon("holm") - 1.0) < 1e-8,
      f"(ratio {_rcon('holm_finite', ri=1e9, rj=1e9)/_rcon('holm'):.9f})")

_tab23 = MR.pair("AlN", "untreated", "AlN", "untreated", "epoxy")
def _gspot(model, cut="R_min"):
    c = tsf.AnalysisConfig(constriction_model=model, constriction_cut=cut,
                           contact_composition="exclusive", bulk_transport=False)
    _, at = tsf.pair_conductance(-1.0, 10.0, 10.0, 1, 1, _tab23, c)
    return at["R_con"]
check("T23d engine: current / holm = 2.000000",
      abs(_gspot("current") / _gspot("holm") - 2.0) < 1e-12,
      f"(ratio {_gspot('current')/_gspot('holm'):.9f})")
check("T23e engine: holm_finite(R_min) is smaller than holm (finite cutoff)",
      _gspot("holm_finite") < _gspot("holm"),
      f"(ratio {_gspot('holm_finite')/_gspot('holm'):.6f})")
check("T23f engine: cutoff choice r_tilde < R_min <= R_each (consistent with the ordering of the cutoff radii)",
      _gspot("holm_finite", "r_tilde") < _gspot("holm_finite", "R_min")
      <= _gspot("holm_finite", "R_each"),
      f"(r~ {_gspot('holm_finite','r_tilde'):.4g} / R_min {_gspot('holm_finite','R_min'):.4g}"
      f" / R_each {_gspot('holm_finite','R_each'):.4g})")

import os as _os
_sv = _os.environ.pop("THERMAL_VARIANT", None)
try:
    _c0 = tsf.AnalysisConfig()
    check("T24 with nothing specified, the current default applies (holm / no shortening)",
          tsf.active_variant() == "holm"
          and (_c0.constriction_model, _c0.constriction_cut,
               _c0.bulk_excludes_constriction) == ("holm", "R_min", False),
          f"({_c0.constriction_model}/{_c0.constriction_cut}/{_c0.bulk_excludes_constriction})")
    check("T24a engine_stamp contains engine_variant",
          tsf.engine_stamp().get("engine_variant") == "holm")

    _os.environ["THERMAL_VARIANT"] = "current"
    _cg = tsf.AnalysisConfig()
    check("T24g THERMAL_VARIANT=current reproduces the legacy default",
          (_cg.constriction_model, _cg.constriction_cut,
           _cg.bulk_excludes_constriction) == ("current", "R_min", False),
          f"({_cg.constriction_model}/{_cg.constriction_cut}/{_cg.bulk_excludes_constriction})")
    _os.environ.pop("THERMAL_VARIANT", None)

    _os.environ["THERMAL_VARIANT"] = "c"
    _c = tsf.AnalysisConfig()
    check("T24b THERMAL_VARIANT=c is a withdrawn variant, kept for sensitivity analysis",
          (_c.constriction_model, _c.constriction_cut, _c.bulk_excludes_constriction)
          == ("holm_finite", "r_tilde", True),
          f"({_c.constriction_model}/{_c.constriction_cut}/{_c.bulk_excludes_constriction})")
    check("T24c engine_stamp under the variant is 'c'",
          tsf.engine_stamp().get("engine_variant") == "c")
    _os.environ["THERMAL_VARIANT"] = "holm"
    _ch = tsf.AnalysisConfig()
    check("T24h THERMAL_VARIANT=holm equals the default",
          (_ch.constriction_model, _ch.bulk_excludes_constriction) == ("holm", False))
    _os.environ["THERMAL_VARIANT"] = "current"
    check("T24d an explicit caller-side setting takes precedence over the variant",
          tsf.AnalysisConfig(constriction_model="holm_finite").constriction_model
          == "holm_finite")

    _os.environ["THERMAL_VARIANT"] = "typo_xxx"
    _raised = False
    try:
        tsf.active_variant()
    except ValueError:
        _raised = True
    check("T24e an unknown variant name raises an exception (no silent fallback to current)", _raised)
finally:
    _os.environ.pop("THERMAL_VARIANT", None)
    if _sv is not None:
        _os.environ["THERMAL_VARIANT"] = _sv
check("T24f the environment is restored after the check",
      _os.environ.get("THERMAL_VARIANT", None) == _sv)

print()
if FAILED:
    print(f"*** {len(FAILED)} FAIL: {FAILED} ***")
    sys.exit(1)
print("ALL SELFTEST PASS")
