# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L0_engine/thermal_sheaf_filtration.py
#   sha256(src) : df0b5b4c58af80c5
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 66 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical. In addition, 29 top-level definitions (and 6 module-level
# statements used only by them) that the figures and verification scripts of
# the paper do not use were omitted (they are listed by name in rb_manifest);
# the remaining definitions are unchanged.
# ---------------------------------------------------------------------------
# Core of the software package accompanying the paper: the branch-conductance
# model (Sec. 2), the conductance filtration with sink-relative persistence
# (Sec. 3) and the dissipation share and finite-change formulas (Sec. 4).
# Units: lengths in um, conductivities in W/(m K), interfacial resistances per
# unit area in K m^2/W, edge conductances g_e in W/K.
# The published copy of this file contains only the routines that the figures and
# verification scripts of the paper use; the omitted definitions are listed in rb_manifest.

from __future__ import annotations

from dataclasses import dataclass, field, fields as _dc_fields

# Version tag stamped into every result CSV with the variant and the baseline fingerprint.
ENGINE_VERSION = "v4"
ENGINE_DEFAULTS_LOG = {
    "v4": "v4: through-particle resistance R_thr added in series on contact and near-contact edges (bulk_transport=True), R_thr,i = x_i/(kappa_i pi R_i^2), x_i = R_i on near-contact edges; per-edge form (a particle with z contacts carries the term z times). Sparse solver path added to shorted_edge_resistance (defaults unchanged)."
          ""
          ""
          ""
          ""
          ""
          ""
          ""
          "",
}


def engine_stamp() -> dict:
    return {"engine_version": ENGINE_VERSION,
            "engine_note": ENGINE_DEFAULTS_LOG[ENGINE_VERSION],
            "engine_variant": active_variant(),
            "baseline_fp": _baseline_fp()}


def _baseline_fp() -> str:
    global _BASELINE_FP_CACHE
    if _BASELINE_FP_CACHE is not None:
        return _BASELINE_FP_CACHE
    try:
        import hashlib
        import pathlib as _pl
        j = _pl.Path(__file__).resolve().with_name("baseline.json")
        _BASELINE_FP_CACHE = hashlib.sha256(j.read_bytes()).hexdigest()[:16]
    except Exception:
        _BASELINE_FP_CACHE = ""
    return _BASELINE_FP_CACHE


_BASELINE_FP_CACHE = None


import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial import cKDTree
from scipy import sparse
from scipy.sparse.linalg import spsolve


# All model settings. The defaults of the model switches (contact_radius, constriction_model,
# contact_composition, prox_flux_order, bob_mode, bulk_transport, barcode_pairing) are those of
# Table 1; the geometric settings (gap_max_um, slab_frac, k_filler, filename) are set by the callers.
@dataclass
class AnalysisConfig:
    filename: str = "pack_last.lammpstrj"

    k_filler: float = 25.0
    k_matrix: float = 0.2
    # absolute distance cut-off on the separation h [um] (Sec. 2.6)
    gap_max_um: float = 1.5
    # interaction cut-off: a pair is an edge only if ln(1 + r~/h) >= eta_cut (Sec. 2.6)
    eta_cut: float = 1.0
    barcode_pairing: str = "relative_h0"  # "relative_h0": persistence relative to the sink with the elder rule (Sec. 3.4)
    fixed_complex: bool = False
    fixed_complex_gap_um: float = 0.0
    # through-particle resistance R_thr in series on every edge (Sec. 2.5)
    bulk_transport: bool = True
    area_ratio_r_a: float = 1.0
    t_min_um: float = 1e-5
    delta_min_um: float = 1e-4
    bob_saturation: bool = True
    # "shift": saturation as the length h_BOB = r~/(alpha^2 - 1) added to the gap (Sec. 2.2.2)
    bob_mode: str = "shift"
    contact_eps_um: float = 0.0
    roughness_um: float = 0.0
    contact_radius: str = "geometric"  # "geometric": a is the radius of the intersection circle of the two spheres (Sec. 2.3.1)
    constriction_model: str = "holm"  # "holm": R_con,i = 1/(4 kappa_i a) (Sec. 2.3.1); other values kept for sensitivity studies
    bulk_excludes_constriction: bool = False
    constriction_cut: str = "R_min"
    contact_composition: str = "parallel"  # "parallel": g_contact = g_spot + g_ann (Sec. 2.3)
    prox_flux_order: str = "parallel_of_series"  # "parallel_of_series": g_prox = 2 pi k_gap r~ ln(1 + r~/(h + a_K,i + a_K,j + h_BOB)) (Sec. 2.2)
    prox_constriction_beta: float = 0.0

    slab_axis: int = 2
    # boundary slab thickness as a fraction of the box length; particles centred in it form src / snk (Sec. 2.1)
    slab_frac: float = 0.05
    # conductance of the edges joining the virtual terminals to the slab particles, max(g) * factor
    g_virtual_factor: float = 1e3

    n_theta_spectral: int = 25
    heat_proxy: str = "volume"
    # overlapping pairs with the same mol id (>= 1) are bond branches (Sec. 2.4)
    bond_within_mol: bool = True

    outdir: str = "results"
    viz: bool = True
    viz_proj: str = "xz"
    viz_vmax: float | None = None
    redundancy: bool = True
    redundancy_max_n: int = 1500

    def __post_init__(self):
        ov = VARIANTS[active_variant()]
        for k, v in ov.items():
            if getattr(self, k) == _CFG_FIELD_DEFAULTS[k]:
                object.__setattr__(self, k, v)


_CFG_FIELD_DEFAULTS = {f.name: f.default for f in _dc_fields(AnalysisConfig)}


# Model variants (environment variable THERMAL_VARIANT). "holm" is the default and
# the one used in the paper; the others are kept for sensitivity studies.
VARIANTS = {
    "holm":    {},
    "current": {"constriction_model": "current",
                "constriction_cut": "R_min",
                "bulk_excludes_constriction": False},
    "c":       {"constriction_model": "holm_finite",
                "constriction_cut": "r_tilde",
                "bulk_excludes_constriction": True},
}
_VARIANT_ENV = "THERMAL_VARIANT"
_DEFAULT_VARIANT = "holm"


def active_variant() -> str:
    import os
    v = os.environ.get(_VARIANT_ENV, _DEFAULT_VARIANT) or _DEFAULT_VARIANT
    if v not in VARIANTS:
        raise ValueError(
            f"{_VARIANT_ENV}={v!r} is an unknown variant. Available: {sorted(VARIANTS)}."
            "(a typo silently falling back to current would leave mixed outputs behind unnoticed, hence an exception)")
    return v


SOURCE_NODE = -1
SINK_NODE = -2


# Per-material properties: kappa [W/(m K)]; interfacial resistance per unit area
# [K m^2/W] at a solid-solid face (R_s_contact) and at the matrix (R_s_matrix).
@dataclass
class Material:
    name: str
    kappa: float
    R_s_contact: float = 0.0
    R_s_matrix: float = 0.0
    t_ip_um: float = 0.0
    k_ip: float | None = None


# Material of each particle type plus the matrix conductivity k_mat.
@dataclass
class MaterialTable:
    materials: dict
    k_matrix: float = 0.2
    contact_Rc_override: dict = field(default_factory=dict)

    def mat(self, t: int) -> Material:
        try:
            return self.materials[int(t)]
        except KeyError:
            raise KeyError(
                f"atom type {t} is not in the MaterialTable. "
                f"registered: {sorted(self.materials)}")

    def contact_Rc(self, ti: int, tj: int) -> float:
        key = tuple(sorted((int(ti), int(tj))))
        if key in self.contact_Rc_override:
            return self.contact_Rc_override[key]
        return self.mat(ti).R_s_contact + self.mat(tj).R_s_contact

    @classmethod
    def uniform(cls, types, kappa: float, k_matrix: float, **mat_kw):
        mats = {int(t): Material(f"type{int(t)}", kappa, **mat_kw)
                for t in set(int(x) for x in types)}
        return cls(materials=mats, k_matrix=k_matrix)


# Radius a of the intersection circle of two spheres overlapping by delta
# (Sec. 2.3.1): D = r_i + r_j - delta, x_i = (D^2 + r_i^2 - r_j^2)/(2D),
# a^2 = r_i^2 - x_i^2. Returns min(r_i, r_j) if one sphere lies inside the other.
def geometric_contact_radius(ri_um, rj_um, delta_um):
    D = ri_um + rj_um - delta_um
    if D <= abs(ri_um - rj_um):
        return float(min(ri_um, rj_um))
    if D <= 0:
        return float(min(ri_um, rj_um))
    x_i = (D * D + ri_um * ri_um - rj_um * rj_um) / (2.0 * D)
    a2 = ri_um * ri_um - x_i * x_i
    return float(np.sqrt(a2)) if a2 > 0.0 else 0.0


def spot_radius_um(gap_um, r_t_um, aK2_um):
    return float(np.sqrt(2.0 * r_t_um * max(gap_um + aK2_um, 0.0)))


# Conductance g_e [W/K] of one particle pair (Sec. 2). h <= 0: contact branch
# g_spot + g_ann; h > 0: near-contact branch g_prox; R_thr is then added in
# series. Returns (None, {}) beyond the interaction cut-off; the dict holds the parts of 1/g_e.
def pair_conductance(gap_um, ri_um, rj_um, ti, tj,
                     table: MaterialTable, cfg: AnalysisConfig):
    mi, mj = table.mat(ti), table.mat(tj)
    r_t_um = ri_um * rj_um / (ri_um + rj_um)
    r_t_m = r_t_um * 1e-6

    sigma = getattr(cfg, "roughness_um", 0.0)
    if sigma > 0.0:
        import copy as _copy
        cfg0 = _copy.copy(cfg)
        cfg0.roughness_um = 0.0
        parts, at = [], {}
        if gap_um > 0.0:
            g_p, a_p = pair_conductance(gap_um, ri_um, rj_um, ti, tj, table, cfg0)
            if g_p is not None:
                parts.append(g_p)
                at.update({f"gap_{k}": v for k, v in a_p.items()})
        delta_a = sigma - gap_um
        if delta_a > 0.0:
            g_c, a_c = pair_conductance(-delta_a, ri_um, rj_um, ti, tj,
                                        table, cfg0)
            if g_c is not None:
                parts.append(g_c)
                at.update({f"con_{k}": v for k, v in a_c.items()})
        if not parts:
            return None, {}
        g_tot = float(sum(parts))
        at.update({"branch": ("mixed" if len(parts) == 2
                              else ("gap" if gap_um > 0 else "contact")),
                   "R": 1.0 / g_tot, "roughness_um": float(sigma),
                   "delta_asperity_um": float(max(delta_a, 0.0)),
                   "R_con": at.get("con_R_con", 0.0) + at.get("gap_R_con", 0.0),
                   "R_int": at.get("gap_R_int", 0.0) + at.get("con_R_int", 0.0),
                   "R_gap": at.get("gap_R_gap", 0.0),
                   "gap_um_eff": float(gap_um)})
        return g_tot, at

    # ---- contact branch (h <= 0) ----
    if gap_um <= getattr(cfg, "contact_eps_um", 0.0):
        delta = max(-gap_um, cfg.delta_min_um)
        if getattr(cfg, "contact_radius", "hertz") == "geometric":
            a_um = geometric_contact_radius(ri_um, rj_um, delta)
            a_um = max(a_um, 1e-30)
        else:
            a_um = np.sqrt(delta * r_t_um)
        a_m = a_um * 1e-6
        _cm = getattr(cfg, "constriction_model", "current")
        _cut_i = _cut_j = 0.0
        if _cm == "current":
            R_con = 1.0 / (2.0 * mi.kappa * a_m) + 1.0 / (2.0 * mj.kappa * a_m)
        elif _cm == "holm":
            R_con = 1.0 / (4.0 * mi.kappa * a_m) + 1.0 / (4.0 * mj.kappa * a_m)  # constriction resistance of both particles, 1/(4 kappa_i a) + 1/(4 kappa_j a)
        elif _cm == "holm_finite":
            _cut = getattr(cfg, "constriction_cut", "R_min")
            if _cut == "R_each":
                _Ri, _Rj = ri_um, rj_um
            elif _cut == "r_tilde":
                _Ri = _Rj = r_t_um
            else:
                _Ri = _Rj = min(ri_um, rj_um)
            R_con = (np.arctan(_Ri / a_um) / (2.0 * np.pi * mi.kappa * a_m)
                     + np.arctan(_Rj / a_um) / (2.0 * np.pi * mj.kappa * a_m))
            _cut_i, _cut_j = float(_Ri), float(_Rj)
        else:
            raise ValueError(f"invalid constriction_model: {_cm!r}")
        _ra = max(getattr(cfg, "area_ratio_r_a", 1.0), 1e-12)
        R_int = table.contact_Rc(ti, tj) / (_ra * np.pi * a_m ** 2)  # interfacial resistance of the contact face, (R_c,i + R_c,j)/(pi a^2)
        R = R_con + R_int
        g_spot = 1.0 / R  # solid path: g_spot = 1/(R_con + R_int) (Sec. 2.3.1)
        attrs = {"branch": "contact", "a_um": float(a_um),
                 "a_def": getattr(cfg, "contact_radius", "hertz"),
                 "R_con": R_con, "R_int": R_int, "R_gap": 0.0,
                 "g_spot": float(g_spot), "g_annulus": 0.0}

        # matrix annulus outside the intersection circle (Sec. 2.3.2):
        #   g_ann = 2 pi k_mat r~ ln(1 + (r~ - delta)/(a_K,i + a_K,j + h_BOB)),
        # added in parallel with g_spot.
        if getattr(cfg, "contact_composition", "exclusive") == "parallel":
            k_m0 = table.k_matrix
            _ra = max(getattr(cfg, "area_ratio_r_a", 1.0), 1e-12)
            aK2_um = ((mi.R_s_matrix + mj.R_s_matrix) / _ra
                      * k_m0 * 1e6)
            aK2_um = max(aK2_um, cfg.t_min_um)
            span = r_t_um - delta
            if span > 0.0:
                h_bob0 = 0.0
                kap_h0 = 2.0 / (1.0 / mi.kappa + 1.0 / mj.kappa)
                al0 = kap_h0 / k_m0
                if getattr(cfg, "bob_saturation", True) and al0 > 1.0:
                    h_bob0 = r_t_um / (al0 ** 2 - 1.0)
                if getattr(cfg, "bob_mode", "shift") == "cap":
                    g_ann = 2.0 * np.pi * k_m0 * r_t_m * np.log1p(span / aK2_um)
                    if getattr(cfg, "bob_saturation", True) and al0 > 1.0:
                        g_cap = 4.0 * np.pi * k_m0 * r_t_m * np.log(al0)
                        g_ann = min(g_ann, g_cap)
                else:
                    g_ann = 2.0 * np.pi * k_m0 * r_t_m * np.log1p(
                        span / (aK2_um + h_bob0))
                beta0 = getattr(cfg, "prox_constriction_beta", 0.0)
                if beta0 > 0.0:
                    rho0_m = spot_radius_um(0.0, r_t_um, aK2_um) * 1e-6
                    if rho0_m > 0:
                        R_sp0 = beta0 * (1.0 / mi.kappa + 1.0 / mj.kappa)\
                            / (2.0 * rho0_m)
                        g_ann = 1.0 / (1.0 / g_ann + R_sp0)
                attrs["g_annulus"] = float(g_ann)
                R = 1.0 / (g_spot + g_ann)
                attrs["R"] = R
        attrs["rho_spot_um"] = float(a_um)
        # through-particle resistance in series with the whole edge (Sec. 2.5):
        #   R_thr = x_i/(kappa_i pi r_i^2) + x_j/(kappa_j pi r_j^2),
        # x_i = distance from the intersection circle to the centre of particle i,
        # cross-section = great circle pi r_i^2.
        if getattr(cfg, "bulk_transport", False):
            _Dc = max(ri_um + rj_um - delta, 1e-12)
            _xi = (_Dc ** 2 + ri_um ** 2 - rj_um ** 2) / (2.0 * _Dc)
            _xj = _Dc - _xi
            _li, _lj = _xi, _xj
            if getattr(cfg, "bulk_excludes_constriction", False):
                _li = max(_xi - _cut_i, 0.0)
                _lj = max(_xj - _cut_j, 0.0)
            _Rb = (_li * 1e-6 / (mi.kappa * np.pi * (ri_um * 1e-6) ** 2)
                   + _lj * 1e-6 / (mj.kappa * np.pi * (rj_um * 1e-6) ** 2))
            attrs["R_bulk"] = float(_Rb)
            attrs["bulk_len_i_um"], attrs["bulk_len_j_um"] = float(_li), float(_lj)
            attrs["x_i_um"], attrs["x_j_um"] = float(_xi), float(_xj)
            R = R + _Rb
            attrs["R"] = R
        else:
            attrs["R_bulk"] = 0.0
    # ---- near-contact branch (h > 0) ----
    else:
        gap = max(gap_um, cfg.t_min_um)
        gap_inf_um = r_t_um / (np.exp(cfg.eta_cut) - 1.0)
        # interaction cut-off (Sec. 2.6): drop the pair if ln(1 + r~/h) < eta_cut
        # or h > gap_max_um.
        if not getattr(cfg, "fixed_complex", False):
            if gap_um > min(gap_inf_um, cfg.gap_max_um):
                return None, {}
        else:
            _cap = (cfg.fixed_complex_gap_um if cfg.fixed_complex_gap_um > 0
                    else cfg.gap_max_um)
            if gap_um > _cap:
                return None, {}
        t_i, t_j = mi.t_ip_um, mj.t_ip_um
        k_i = mi.k_ip if mi.k_ip is not None else table.k_matrix
        k_j = mj.k_ip if mj.k_ip is not None else table.k_matrix
        tot_ip = t_i + t_j
        scale = min(1.0, gap / tot_ip) if tot_ip > 0 else 0.0
        L_i, L_j = t_i * scale, t_j * scale
        L_m = max(gap - L_i - L_j, 0.0)
        k_eff = gap / (L_i / k_i + L_j / k_j + L_m / table.k_matrix)  # k_gap: series-effective conductivity across the gap (equals k_mat without an interphase layer)

        ln_fac = np.log1p(r_t_um / gap)
        g_gap = 2.0 * np.pi * k_eff * r_t_m * ln_fac

        bob_sat, bob_alpha = False, float("nan")
        bob_sat_cap = None
        bob_h_shift = None
        # saturation length (Sec. 2.2.2): alpha = kappa_h/k_mat with kappa_h the
        # harmonic mean of the two particle conductivities, h_BOB = r~/(alpha^2 - 1).
        if getattr(cfg, "bob_saturation", True):
            kap_h = 2.0 / (1.0 / mi.kappa + 1.0 / mj.kappa)
            bob_alpha = kap_h / table.k_matrix
            if bob_alpha > 1.0:
                g_sat = 4.0 * np.pi * k_eff * r_t_m * np.log(bob_alpha)
                bob_sat_cap = g_sat
                bob_h_shift = r_t_um / (bob_alpha ** 2 - 1.0)
                if getattr(cfg, "bob_mode", "shift") == "cap" and g_sat < g_gap:
                    g_gap, bob_sat = g_sat, True

        R_gap = 1.0 / g_gap
        A_th = g_gap * (gap * 1e-6) / k_eff
        R_int = (mi.R_s_matrix + mj.R_s_matrix) / A_th

        beta = getattr(cfg, "prox_constriction_beta", 0.0)
        R_sp = 0.0
        if beta > 0.0 and A_th > 0.0:
            rho_m = np.sqrt(A_th / np.pi)
            R_sp = beta * (1.0 / mi.kappa + 1.0 / mj.kappa) / (2.0 * rho_m)

        # near-contact conductance (Sec. 2.2), with a_K,i = R_s,i k_gap:
        #   g_prox = 2 pi k_gap r~ ln(1 + r~/(h + a_K,i + a_K,j + h_BOB)).
        order = getattr(cfg, "prox_flux_order", "series_of_parallels")
        if order == "parallel_of_series":
            _ra = max(getattr(cfg, "area_ratio_r_a", 1.0), 1e-12)
            Rs_sum = (mi.R_s_matrix + mj.R_s_matrix) / _ra
            a_K2_um = Rs_sum * k_eff * 1e6
            h_bob = 0.0
            if bob_h_shift is not None:
                h_bob = bob_h_shift
            ln_shift = np.log1p(r_t_um / (gap + a_K2_um + h_bob))
            g_pos = 2.0 * np.pi * k_eff * r_t_m * ln_shift
            if getattr(cfg, "bob_mode", "shift") == "cap":
                if bob_sat_cap is not None and bob_sat_cap < g_pos:
                    g_pos, bob_sat = bob_sat_cap, True
                else:
                    bob_sat = False
            else:
                bob_sat = bool(h_bob > gap + a_K2_um)
            R_gap = 1.0 / g_pos
            R_int = 0.0
            A_th = g_pos * (gap * 1e-6) / k_eff
            rho_m = spot_radius_um(gap, r_t_um, a_K2_um) * 1e-6
            if beta > 0.0 and rho_m > 0.0:
                R_sp = beta * (1.0 / mi.kappa + 1.0 / mj.kappa) / (2.0 * rho_m)
            else:
                R_sp = 0.0

        # through-particle resistance with x_i = r_i on a near-contact edge (Sec. 2.5).
        R_bulk = 0.0
        if getattr(cfg, "bulk_transport", False):
            if not getattr(cfg, "bulk_excludes_constriction", False):
                R_bulk = (1.0 / (mi.kappa * np.pi * ri_um * 1e-6)
                          + 1.0 / (mj.kappa * np.pi * rj_um * 1e-6))
            else:
                _cut = getattr(cfg, "constriction_cut", "R_min")
                if _cut == "R_each":
                    _ci, _cj = ri_um, rj_um
                elif _cut == "r_tilde":
                    _ci = _cj = r_t_um
                else:
                    _ci = _cj = min(ri_um, rj_um)
                _li, _lj = max(ri_um - _ci, 0.0), max(rj_um - _cj, 0.0)
                R_bulk = (_li * 1e-6 / (mi.kappa * np.pi * (ri_um * 1e-6) ** 2)
                          + _lj * 1e-6 / (mj.kappa * np.pi * (rj_um * 1e-6) ** 2))
        R = R_gap + R_int + R_sp + R_bulk
        attrs = {"branch": "gap", "k_gap_eff": float(k_eff),
                 "R_bulk": float(R_bulk),
                 "A_th_m2": float(A_th), "flux_order": order,
                 "R_con": R_sp, "R_int": R_int, "R_gap": R_gap,
                 "rho_um": float(np.sqrt(A_th / np.pi) * 1e6) if A_th > 0 else 0.0,
                 "rho_spot_um": spot_radius_um(
                     gap, r_t_um,
                     (mi.R_s_matrix + mj.R_s_matrix) * k_eff * 1e6),
                 "bob_saturated": bool(bob_sat), "bob_alpha": float(bob_alpha),
                 "bob_h_shift_um": (float(bob_h_shift)
                                    if bob_h_shift is not None else 0.0)}
    g = 1.0 / max(R, 1e-300)
    attrs["R"] = R
    return g, attrs


def pairs_from_positions(positions, diameters, cfg: AnalysisConfig):
    radii = np.asarray(diameters) / 2.0
    r_max_center = 2.0 * radii.max() + max(cfg.gap_max_um,
                                           cfg.fixed_complex_gap_um)
    tree = cKDTree(positions)
    return sorted(tree.query_pairs(r=r_max_center))


# Builds the weighted graph G = (V, E, g): a k-d tree lists candidate pairs and
# pair_conductance classifies each pair and assigns g_e (Sec. 5.1). Edge
# attributes: gap_um, branch ("contact" / "gap" / "bond"), g and the parts of 1/g.
def build_thermal_network(ids, types, positions, diameters, cfg: AnalysisConfig,
                          table: MaterialTable | None = None, mols=None,
                          fixed_pairs=None, g_floor: float = 1e-30) -> nx.Graph:
    if table is None:
        table = MaterialTable.uniform(types, cfg.k_filler, cfg.k_matrix)

    radii = diameters / 2.0
    r_max_center = 2.0 * radii.max() + cfg.gap_max_um
    tree = cKDTree(positions)
    pairs = tree.query_pairs(r=r_max_center, output_type="ndarray")

    mol_arr = None if mols is None else np.asarray(mols).astype(int)
    G = nx.Graph()
    for k, (pid, t) in enumerate(zip(ids, types)):
        attr = {"atom_type": int(t), "pos": positions[k]}
        if mol_arr is not None:
            attr["mol"] = int(mol_arr[k])
        G.add_node(int(pid), **attr)

    _use = list(fixed_pairs) if fixed_pairs is not None else pairs
    for i_idx, j_idx in _use:
        d = float(np.linalg.norm(positions[i_idx] - positions[j_idx]))
        gap = d - (radii[i_idx] + radii[j_idx])
        g, attrs = pair_conductance(gap, radii[i_idx], radii[j_idx],
                                    types[i_idx], types[j_idx], table, cfg)
        if g is None:
            if fixed_pairs is None:
                continue
            g = g_floor
            attrs = {"branch": "gap", "R_con": 0.0, "R_int": 0.0,
                     "R_gap": 1.0 / g_floor, "R": 1.0 / g_floor,
                     "out_of_range": True}
        if (getattr(cfg, "bond_within_mol", True) and mol_arr is not None
                and mol_arr[i_idx] == mol_arr[j_idx] and mol_arr[i_idx] >= 1):
            # bond branch (Sec. 2.4): axial conduction through the cross-section,
            #   g_bond = kappa_h pi R_min^2 / L,  L = centre-to-centre distance;
            # no interface, no constriction and no additional R_thr.
            if gap < 0.0:
                _ri, _rj = float(radii[i_idx]), float(radii[j_idx])
                _mi, _mj = table.mat(types[i_idx]), table.mat(types[j_idx])
                _kap_h = 2.0 / (1.0 / _mi.kappa + 1.0 / _mj.kappa)
                _R_min_m = min(_ri, _rj) * 1e-6
                _L_m = max(d, 1e-30) * 1e-6
                _A_m2 = np.pi * _R_min_m ** 2
                _Rbulk = _L_m / (_kap_h * _A_m2)
                g = 1.0 / _Rbulk
                _a_um = geometric_contact_radius(_ri, _rj, -gap)
                attrs = {"branch": "bond", "a_um": float(_a_um),
                         "a_def": "cylinder_axial",
                         "R_con": 0.0, "R_int": 0.0, "R_gap": 0.0,
                         "R_axial": _Rbulk, "R_bulk": 0.0, "R": _Rbulk,
                         "g_spot": float(g), "g_annulus": 0.0,
                         "A_cross_m2": float(_A_m2),
                         "rho_spot_um": float(min(_ri, _rj))}
        G.add_edge(int(ids[i_idx]), int(ids[j_idx]),
                   distance_um=d, gap_um=float(gap),
                   r_birth=max(0.0, float(gap)),
                   g=float(g), **attrs)
    return G


# Source and sink: particles whose centres lie in the boundary slab at either
# face of the box along slab_axis (Sec. 2.1).
def identify_slabs(ids, positions, cfg: AnalysisConfig):
    z = positions[:, cfg.slab_axis]
    lo, hi = z.min(), z.max()
    thick = (hi - lo) * cfg.slab_frac
    source_ids = {int(i) for i, zi in zip(ids, z) if zi >= hi - thick}
    sink_ids = {int(i) for i, zi in zip(ids, z) if zi <= lo + thick}
    return source_ids, sink_ids


# Adds the virtual terminals src and snk, joined to every slab particle by an
# edge of very large conductance, so that u|src = 1 and u|snk = 0 are imposed.
def attach_virtual_terminals(G: nx.Graph, source_ids, sink_ids, cfg: AnalysisConfig):
    g_virtual = max(d["g"] for _, _, d in G.edges(data=True)) * cfg.g_virtual_factor
    H = G.copy()
    H.add_node(SOURCE_NODE, atom_type=0)
    H.add_node(SINK_NODE, atom_type=0)
    for i in source_ids:
        H.add_edge(SOURCE_NODE, i, g=g_virtual, r_birth=0.0)
    for i in sink_ids:
        H.add_edge(SINK_NODE, i, g=g_virtual, r_birth=0.0)
    return H, g_virtual


# Weighted graph Laplacian L = delta^T delta with (delta u)_e = sqrt(g_e)(u_i - u_j)
# (Sec. 3.2). laplacian(theta) uses only the edges with g_e >= theta.
class CellularSheaf:

    def __init__(self, G: nx.Graph, weight: str = "g"):
        self.G = G
        self.weight = weight
        self.nodes = sorted(G.nodes())
        self.index = {v: k for k, v in enumerate(self.nodes)}

    def coboundary(self, edge_mask=None) -> sparse.csr_matrix:
        edges = list(self.G.edges(data=True)) if edge_mask is None else edge_mask
        m, n = len(edges), len(self.nodes)
        rows, cols, vals = [], [], []
        for k, (u, v, d) in enumerate(edges):
            w = np.sqrt(d[self.weight])
            rows += [k, k]
            cols += [self.index[u], self.index[v]]
            vals += [w, -w]
        return sparse.csr_matrix((vals, (rows, cols)), shape=(m, n))

    def laplacian(self, theta: float | None = None) -> sparse.csr_matrix:
        if theta is None:
            edges = None
        else:
            edges = [(u, v, d) for u, v, d in self.G.edges(data=True)
                     if d[self.weight] >= theta]
        delta = self.coboundary(edges)
        return (delta.T @ delta).tocsr()

    def grounded_laplacian(self, L: sparse.csr_matrix, grounded_nodes) -> tuple:
        keep = [k for k, v in enumerate(self.nodes) if v not in grounded_nodes]
        return L[np.ix_(keep, keep)].tocsr(), [self.nodes[k] for k in keep]


class _UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return ra
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.rank[ra] = max(self.rank[ra], self.rank[rb] + 1)
        return ra


# Procedure (a) of Sec. 5.1: edges are inserted in decreasing order of g and
# components are merged by union-find; every record of sweep() is one bar of
# the persistence relative to the sink (Sec. 3.4). Columns birth_theta and
# sink_theta are conductances g [W/K], not -log10 g; persistence_to_sink is in decades.
class ConductanceFiltration:

    def __init__(self, G_with_terminals: nx.Graph, cfg=None):
        self.G = G_with_terminals
        self.cfg = cfg

    def thresholds(self, n: int) -> np.ndarray:
        g = np.array([d["g"] for _, _, d in self.G.edges(data=True)])
        qs = np.linspace(0.0, 1.0, n)
        return np.quantile(g, 1.0 - qs)

    def sweep(self, pairing: str | None = None) -> pd.DataFrame:
        if pairing is None:
            pairing = getattr(getattr(self, "cfg", None), "barcode_pairing",
                              None) or "relative_h0"
        edges = sorted(self.G.edges(data=True), key=lambda e: e[2]["g"],
                       reverse=True)
        # a bar born at the terminal conductance comes from the terminal edges, not the structure (born_at_terminal).
        _g_term = max((d["g"] for u, v, d in self.G.edges(data=True)
                       if u in (SOURCE_NODE, SINK_NODE)
                       or v in (SOURCE_NODE, SINK_NODE)), default=None)
        uf = _UnionFind(list(self.G.nodes()))
        birth: dict = {}
        has_sink: dict = {}
        sink_theta: dict = {}
        records = []
        for u, v, d in edges:
            theta = d["g"]
            ru, rv = uf.find(u), uf.find(v)
            for r, node in ((ru, u), (rv, v)):
                if r not in birth:
                    birth[r] = theta
                    has_sink[r] = (r == SINK_NODE) or (node == SINK_NODE)
            if ru == rv:
                continue
            su, sv = has_sink.get(ru, False), has_sink.get(rv, False)
            bu, bv = birth.get(ru, theta), birth.get(rv, theta)
            r_new = uf.union(ru, rv)
            # merging with the component that contains the sink: the other bar dies here.
            if su != sv:
                loser = rv if su else ru
                if loser not in sink_theta:
                    sink_theta[loser] = theta
                    _b = birth.get(loser, theta)
                    records.append({
                        "rep": loser,
                        "birth_theta": _b,
                        "sink_theta": theta,
                        "persistence_to_sink": float(np.log10(_b / theta)),
                        "born_at_terminal": bool(
                            _g_term is not None
                            and abs(_b - _g_term) <= 1e-12 * abs(_g_term)),
                    })
            # elder rule: when two components without the sink merge, the younger bar
            # (the one with the smaller birth g) dies (Sec. 3.4).
            elif pairing == "relative_h0" and not su and not sv:
                young_r, young_b = ((ru, bu) if bu < bv else (rv, bv))
                records.append({
                    "rep": young_r,
                    "birth_theta": young_b,
                    "sink_theta": theta,
                    "persistence_to_sink": float(np.log10(young_b / theta)),
                    "born_at_terminal": bool(
                        _g_term is not None
                        and abs(young_b - _g_term) <= 1e-12 * abs(_g_term)),
                })
            has_sink[r_new] = su or sv
            birth[r_new] = max(bu, bv) if pairing == "relative_h0"\
                else min(bu, bv)
        return pd.DataFrame(records).sort_values(
            "persistence_to_sink", ascending=False, ignore_index=True)

    def component_scores(self, theta: float, heat: dict) -> pd.DataFrame:
        eps = 1e-12
        H = nx.Graph()
        H.add_nodes_from(self.G.nodes())
        H.add_edges_from((u, v) for u, v, d in self.G.edges(data=True)
                         if d["g"] >= theta)
        rows = []
        for k, comp in enumerate(nx.connected_components(H)):
            A = set(comp)
            S = sum(heat.get(i, 0.0) for i in A)
            Gbd = sum(d["g"] for u, v, d in self.G.edges(data=True)
                      if (u in A) ^ (v in A))
            rows.append({"theta": theta, "component_id": k, "n_nodes": len(A),
                         "S": S, "G_boundary": Gbd,
                         "stagnation": S / (Gbd + eps),
                         "has_sink": SINK_NODE in A})
        return pd.DataFrame(rows)


# Procedure (b) of Sec. 5.1: g_star is the conductance at which src and snk
# first lie in one component; the critical conductance level is
# theta_star = -log10 g_star (Sec. 3.6). g_all_connected: every vertex has joined the sink.
def percolation_threshold(G_with_terminals: nx.Graph):
    es = sorted(G_with_terminals.edges(data=True), key=lambda e: -e[2]["g"])
    parent = {n: n for n in G_with_terminals.nodes()}
    size = {n: 1 for n in G_with_terminals.nodes()}
    N = G_with_terminals.number_of_nodes()

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    g_star, crit, g_all = None, None, None
    for u, v, d in es:
        a, b = find(u), find(v)
        if a != b:
            parent[b] = a
            size[a] += size[b]
        if g_star is None and find(SOURCE_NODE) == find(SINK_NODE):
            g_star, crit = d["g"], (u, v)
        if g_all is None and size[find(SINK_NODE)] == N:
            g_all = d["g"]
        if g_star is not None and g_all is not None:
            break
    return g_star, crit, g_all


# Solves L u = 0 on the interior vertices with u(src) = 1 and u(snk) = 0, on
# the component containing both terminals (one sparse solve). Returns the vertex
# temperatures and the total heat flow Q = G_eff; (None, 0.0) before percolation.
def solve_dirichlet(sheaf: CellularSheaf, theta: float):
    G = sheaf.G
    w = sheaf.weight
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    edges = [(u, v, d) for u, v, d in G.edges(data=True) if d[w] >= theta]
    H.add_edges_from((u, v) for u, v, _ in edges)
    comp = nx.node_connected_component(H, SOURCE_NODE)
    if SINK_NODE not in comp:
        return None, 0.0

    nodes = sorted(comp)
    idx = {v: k for k, v in enumerate(nodes)}
    n = len(nodes)
    rows, cols, vals = [], [], []
    for u, v, d in edges:
        if u not in comp or v not in comp:
            continue
        i, j, g = idx[u], idx[v], d[w]
        rows += [i, j, i, j]; cols += [i, j, j, i]; vals += [g, g, -g, -g]
    L = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))

    bnd = [idx[SOURCE_NODE], idx[SINK_NODE]]
    T_bnd = np.array([1.0, 0.0])
    interior = [k for k in range(n) if k not in bnd]
    T = np.zeros(n)
    T[bnd] = T_bnd
    if interior:
        A = L[np.ix_(interior, interior)]
        b = -L[np.ix_(interior, bnd)] @ T_bnd
        T[interior] = spsolve(A.tocsc(), b)
    Q = float(-(L[idx[SINK_NODE], :] @ T))
    return dict(zip(nodes, T)), max(Q, 0.0)


# G_eff(theta): effective conductance of the sub-graph with g_e >= theta.
def effective_conductance(sheaf: CellularSheaf, theta: float) -> float:
    _, Q = solve_dirichlet(sheaf, theta)
    return Q


# Effective resistance from each vertex to the sink (used by quantities24).
def _grounded_system(sheaf: CellularSheaf, theta: float):
    G = sheaf.G
    w = sheaf.weight
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    edges = [(u, v, d) for u, v, d in G.edges(data=True) if d[w] >= theta]
    H.add_edges_from((u, v) for u, v, _ in edges)
    comp = nx.node_connected_component(H, SINK_NODE)
    nodes_c = sorted(comp)
    idx_c = {v: k for k, v in enumerate(nodes_c)}
    n = len(nodes_c)
    rows, cols, vals = [], [], []
    for u, v, d in edges:
        if u not in comp or v not in comp:
            continue
        i, j, g = idx_c[u], idx_c[v], d[w]
        rows += [i, j, i, j]; cols += [i, j, j, i]; vals += [g, g, -g, -g]
    L = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
    keep = [k for k, v in enumerate(nodes_c) if v != SINK_NODE]
    L_gnd = L[np.ix_(keep, keep)].tocsc()
    solve = sparse.linalg.factorized(L_gnd)
    interior = [nodes_c[k] for k in keep]
    idx = {v: j for j, v in enumerate(interior)}
    return solve, interior, idx, comp, L_gnd


def effective_resistance_to_sink(sheaf: CellularSheaf, theta: float) -> dict:
    solve, interior, idx, comp, _ = _grounded_system(sheaf, theta)
    R = {}
    n = len(interior)
    for v in sheaf.G.nodes():
        if v in (SOURCE_NODE, SINK_NODE):
            continue
        if v not in comp:
            R[v] = np.inf
            continue
        b = np.zeros(n)
        b[idx[v]] = 1.0
        R[v] = float(solve(b)[idx[v]])
    return R


# Procedure (c) of Sec. 5.1: from one Dirichlet solve, the dissipation
# P_e = g_e (u_i - u_j)^2 of every edge and its share P_share = P_e / sum P_e.
# P_share equals the logarithmic sensitivity d log G_eff / d log g_e (Sec. 4.3) up to the
# (negligible) dissipation in the terminal edges, which is excluded from the normalisation.
def edge_flux_table(sheaf: CellularSheaf, theta: float) -> pd.DataFrame:
    T_map, Q = solve_dirichlet(sheaf, theta)
    cols = ["u", "v", "type_u", "type_v", "branch", "gap_um",
            "g", "dT", "q", "P", "P_share"]
    if T_map is None:
        return pd.DataFrame(columns=cols)
    rows = []
    for u, v, d in sheaf.G.edges(data=True):
        if d[sheaf.weight] < theta:
            continue
        if u in (SOURCE_NODE, SINK_NODE) or v in (SOURCE_NODE, SINK_NODE):
            continue
        if u not in T_map or v not in T_map:
            continue
        g = d[sheaf.weight]
        dT = T_map[u] - T_map[v]
        q = g * dT
        rows.append({"u": u, "v": v,
                     "type_u": int(sheaf.G.nodes[u].get("atom_type", 0)),
                     "type_v": int(sheaf.G.nodes[v].get("atom_type", 0)),
                     "branch": d.get("branch", ""),
                     "gap_um": d.get("gap_um", float("nan")),
                     "g": g, "dT": dT, "q": q, "P": q * dT})
    df = pd.DataFrame(rows, columns=cols[:-1])
    if len(df):
        df["P_share"] = df["P"] / max(df["P"].sum(), 1e-300)
        df = df.sort_values("P_share", ascending=False, ignore_index=True)
    return df


# Procedure (d) of Sec. 5.1: R~_e is the effective resistance between the ends
# of e in the circuit in which src and snk are identified (Sec. 4.4). One
# factorisation of the shorted Laplacian, then one back-substitution per
# requested edge (edges=...); a dense pseudo-inverse is used for small graphs.
def shorted_edge_resistance(sheaf: CellularSheaf, theta: float,
                            max_n: int = 1500, edges=None,
                            method: str = "auto") -> pd.DataFrame:
    G = sheaf.G
    w = sheaf.weight
    cols = ["u", "v", "g", "R_tilde"]
    H = nx.Graph()
    for u, v, d in G.edges(data=True):
        if d[w] >= theta:
            H.add_edge(u, v, g=d[w])
    if (SOURCE_NODE not in H) or (SINK_NODE not in H) or\
            (not nx.has_path(H, SOURCE_NODE, SINK_NODE)):
        return pd.DataFrame(columns=cols)
    Hc = H.subgraph(nx.node_connected_component(H, SOURCE_NODE))
    if method not in ("auto", "dense", "sparse"):
        raise ValueError(f"method must be one of auto/dense/sparse: {method!r}")
    use_sparse = (method == "sparse" or
                  (method == "auto" and
                   (edges is not None or Hc.number_of_nodes() > max_n)))
    if not use_sparse and Hc.number_of_nodes() > max_n:
        return pd.DataFrame(columns=cols)

    MERGED = -99
    Hs = nx.Graph()
    for u, v, d in Hc.edges(data=True):
        a = MERGED if u in (SOURCE_NODE, SINK_NODE) else u
        b = MERGED if v in (SOURCE_NODE, SINK_NODE) else v
        if a == b:
            continue
        if Hs.has_edge(a, b):
            Hs[a][b]["g"] += d["g"]
        else:
            Hs.add_edge(a, b, g=d["g"])

    nodes = list(Hs.nodes()); idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    if edges is None:
        targets = [(int(u), int(v), float(d["g"]))
                   for u, v, d in Hc.edges(data=True)
                   if u not in (SOURCE_NODE, SINK_NODE)
                   and v not in (SOURCE_NODE, SINK_NODE)]
    else:
        want = {frozenset((int(u), int(v))) for u, v in edges}
        targets = [(int(u), int(v), float(d["g"]))
                   for u, v, d in Hc.edges(data=True)
                   if frozenset((int(u), int(v))) in want
                   and u not in (SOURCE_NODE, SINK_NODE)
                   and v not in (SOURCE_NODE, SINK_NODE)]

    if not use_sparse:
        L = np.zeros((n, n))
        for u, v, d in Hs.edges(data=True):
            i, j, g = idx[u], idx[v], d["g"]
            L[i, i] += g; L[j, j] += g; L[i, j] -= g; L[j, i] -= g
        P = np.linalg.pinv(L)
        dP = np.diag(P)
        rows = [{"u": u, "v": v, "g": g,
                 "R_tilde": float(dP[idx[u]] + dP[idx[v]]
                                  - 2 * P[idx[u], idx[v]])}
                for u, v, g in targets]
        return pd.DataFrame(rows, columns=cols)

    import scipy.sparse as _sp
    import scipy.sparse.linalg as _spl

    ii, jj, vv = [], [], []
    for u, v, d in Hs.edges(data=True):
        i, j, g = idx[u], idx[v], d["g"]
        ii += [i, j, i, j]; jj += [i, j, j, i]; vv += [g, g, -g, -g]
    L = _sp.csc_matrix((vv, (ii, jj)), shape=(n, n))
    keep = np.arange(1, n)
    Lg = L[keep][:, keep].tocsc()
    solve = _spl.factorized(Lg)

    def _r(i, j):
        b = np.zeros(n - 1)
        if i > 0:
            b[i - 1] += 1.0
        if j > 0:
            b[j - 1] -= 1.0
        x = solve(b)
        xi = x[i - 1] if i > 0 else 0.0
        xj = x[j - 1] if j > 0 else 0.0
        return float(xi - xj)

    rows = [{"u": u, "v": v, "g": g, "R_tilde": _r(idx[u], idx[v])}
            for u, v, g in targets]
    return pd.DataFrame(rows, columns=cols)


# Closed form of Sec. 4.4, exact for any Delta:
#   G_eff(g + Delta e) = G_eff(g) + Delta (du_e)^2 / (1 + Delta R~_e),
# with du_e from edge_flux_table and R~_e from shorted_edge_resistance.
def predict_finite_change(G_eff: float, dT: float, R_tilde: float,
                          delta_g: float) -> float:
    return G_eff + delta_g * dT ** 2 / (1.0 + delta_g * R_tilde)


# Edge criticality g_e R_eff(u, v) in the ordinary (un-shorted) circuit.
# R_eff_uv is not R~_e and must not be used in the finite-change formula.
def edge_redundancy(sheaf: CellularSheaf, theta: float, max_n: int = 1500):
    G = sheaf.G
    cols = ["u", "v", "g", "R_eff_uv", "crit", "redundancy"]
    H = nx.Graph()
    for u, v, d in G.edges(data=True):
        if d[sheaf.weight] >= theta:
            H.add_edge(u, v, g=d[sheaf.weight])
    if (SOURCE_NODE not in H) or (SINK_NODE not in H) or (not nx.has_path(H, SOURCE_NODE, SINK_NODE)):
        return pd.DataFrame(columns=cols), {}
    comp = nx.node_connected_component(H, SOURCE_NODE)
    Hc = H.subgraph(comp)
    nodes = list(Hc.nodes()); idx = {nn: i for i, nn in enumerate(nodes)}; n = len(nodes)
    summary = {"n_nodes": int(n), "n_edges": int(Hc.number_of_edges())}
    if n > max_n:
        rows_i, cols_i, vals = [], [], []
        for u, v, d in Hc.edges(data=True):
            i, j = idx[u], idx[v]; g = d["g"]
            rows_i += [i, j, i, j]; cols_i += [i, j, j, i]; vals += [g, g, -g, -g]
        L = sparse.csr_matrix((vals, (rows_i, cols_i)), shape=(n, n))
        gnd = idx[SINK_NODE]
        keep = [i for i in range(n) if i != gnd]
        Lr = L[keep][:, keep]
        b = np.zeros(len(keep)); b[keep.index(idx[SOURCE_NODE])] = 1.0
        try:
            x = spsolve(Lr.tocsc(), b)
            Rss = float(x[keep.index(idx[SOURCE_NODE])])
        except Exception:
            Rss = float("nan")
        summary.update({"R_eff_source_sink": Rss,
                        "G_eff_check": (1.0 / Rss if Rss and Rss > 0 else float("inf")),
                        "kirchhoff_index": float("nan"),
                        "note": f"n={n}>max_n={max_n}: all-edge crit skipped (source-sink R_eff only)"})
        return pd.DataFrame(columns=cols), summary
    L = np.zeros((n, n))
    for u, v, d in Hc.edges(data=True):
        i, j = idx[u], idx[v]; g = d["g"]
        L[i, i] += g; L[j, j] += g; L[i, j] -= g; L[j, i] -= g
    P = np.linalg.pinv(L)
    dP = np.diag(P)
    rows = []
    for u, v, d in Hc.edges(data=True):
        if u in (SOURCE_NODE, SINK_NODE) or v in (SOURCE_NODE, SINK_NODE):
            continue
        i, j = idx[u], idx[v]
        Reff = float(dP[i] + dP[j] - 2 * P[i, j]); g = float(d["g"])
        crit = g * Reff
        rows.append({"u": int(u), "v": int(v), "g": g, "R_eff_uv": Reff,
                     "crit": crit, "redundancy": 1.0 - crit})
    df = (pd.DataFrame(rows).sort_values("crit", ascending=False, ignore_index=True)
          if rows else pd.DataFrame(columns=cols))
    si, sj = idx[SOURCE_NODE], idx[SINK_NODE]
    Rss = float(dP[si] + dP[sj] - 2 * P[si, sj])
    summary.update({"R_eff_source_sink": Rss,
                    "G_eff_check": (1.0 / Rss if Rss > 0 else float("inf")),
                    "kirchhoff_index": float(n * np.trace(P))})
    if len(df):
        summary["n_bridges_crit_gt_0.95"] = int((df["crit"] > 0.95).sum())
        summary["mean_crit"] = float(df["crit"].mean())
        summary["median_redundancy"] = float(df["redundancy"].median())
    return df, summary

