# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L2_structures/ensemble.py
#   sha256(src) : 5bdae1e0c6bcba9e
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 29 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical. In addition, 4 module-level statements (a command-line block,
# imports or constants) used only by definitions that are not part of this
# record were omitted (listed in rb_manifest); the remaining code is
# unchanged.
# ---------------------------------------------------------------------------
# Catalogue and generator of the structure families of Appendix C (Table C.2):
# random packings at phi = 0.50-0.64, bimodal / polydisperse mixtures, poor
# dispersion (voids, clusters), stratification / segregation, fibres, treatment
# patchiness. make(family, d_um, seed) returns one packing; build(...) also
# returns the graph with terminals and the material table.
import sys as _sys, pathlib as _pathlib
_LIB = _pathlib.Path(__file__).resolve().parent.parent
for _d in ("L0_engine", "L1_materials", "L2_structures"):
    _p = str(_LIB / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import zlib

import numpy as np
from scipy.spatial import cKDTree

import materials_real as MR
import synth_packing as SP
import thermal_sheaf_filtration as tsf

FAMILIES = {
    "phi_050": dict(group="reference series", phi=0.50),
    "phi_055": dict(group="reference series", phi=0.55),
    "phi_060": dict(group="reference series", phi=0.60),
    "phi_064": dict(group="reference series", phi=0.64),
    "bimodal_3to1": dict(group="high loading", phi=0.66, bimodal=(3.0, 0.70)),
    "bimodal_5to1": dict(group="high loading", phi=0.70, bimodal=(5.0, 0.75)),
    "poly_wide": dict(group="high loading", phi=0.64, span=0.50),
    "voided": dict(group="poor dispersion", phi=0.60, void=0.26),
    "clustered": dict(group="poor dispersion", phi=0.60, cluster=0.32),
    "stratified": dict(group="process segregation", phi=0.60, strat=0.55),
    "segregated": dict(group="process segregation", phi=0.66, bimodal=(3.0, 0.70),
                       segregate=1.0),
    "fiber_z": dict(group="anisotropy", phi=0.60, fiber=("z", 24, 6)),
    "fiber_xy": dict(group="anisotropy", phi=0.60, fiber=("xy", 24, 6)),
    "coat_mix": dict(group="non-uniform treatment", phi=0.60, coat_frac=0.30),
    "strat_025": dict(group="segregation strength", phi=0.60, strat=0.25),
    "strat_055": dict(group="segregation strength", phi=0.60, strat=0.55),
    "strat_100": dict(group="segregation strength", phi=0.60, strat=1.00),
    "strat_180": dict(group="segregation strength", phi=0.60, strat=1.80),
    "seg_000": dict(group="segregation strength", phi=0.66, bimodal=(3.0, 0.70),
                    segregate=0.0),
    "seg_033": dict(group="segregation strength", phi=0.66, bimodal=(3.0, 0.70),
                    segregate=0.33),
    "seg_067": dict(group="segregation strength", phi=0.66, bimodal=(3.0, 0.70),
                    segregate=0.67),
    "seg_100": dict(group="segregation strength", phi=0.66, bimodal=(3.0, 0.70),
                    segregate=1.0),
    "jittered_ref": dict(group="legacy baseline", kind="jittered", phi=None),
}

try:
    from family_notes import NOTES as _NOTES
    for _k, _v in _NOTES.items():
        if _k in FAMILIES:
            FAMILIES[_k].update(_v)
except ImportError:
    pass
FAMILY_ORDER = list(FAMILIES)
GROUPS = []
for _k in FAMILY_ORDER:
    _g = FAMILIES[_k]["group"]
    if _g not in GROUPS:
        GROUPS.append(_g)


N_TARGET = 900


def _diameters(n, d_mean, cfgf, rng):
    if "bimodal" in cfgf:
        ratio, vfrac_large = cfgf["bimodal"]
        wL = vfrac_large / ratio ** 3
        wS = 1.0 - vfrac_large
        nL = max(1, int(round(n * wL / (wL + wS))))
        d = np.concatenate([np.full(nL, ratio), np.full(n - nL, 1.0)])
        rng.shuffle(d)
    elif cfgf.get("span", 0.0) > 0.0:
        d = rng.lognormal(mean=0.0, sigma=cfgf["span"], size=n)
    else:
        d = np.ones(n)
    d = d * (d_mean ** 3 * n / (d ** 3).sum()) ** (1.0 / 3.0)
    return d


def _relax(pos, dia, box, n_iter, grav=0.0, push=0.5):
    r = dia / 2.0
    rmax = r.max()
    for it in range(n_iter):
        if grav > 0.0:
            pos[:, 2] -= grav * rmax * (1.0 - it / n_iter)
        tree = cKDTree(pos)
        pairs = tree.query_pairs(2.0 * rmax, output_type="ndarray")
        if len(pairs):
            i, j = pairs[:, 0], pairs[:, 1]
            dvec = pos[j] - pos[i]
            dist = np.linalg.norm(dvec, axis=1)
            ov = (r[i] + r[j]) - dist
            m = (ov > 0) & (dist > 1e-12)
            if m.any():
                u = dvec[m] / dist[m][:, None]
                shift = (push * ov[m])[:, None] * u * 0.5
                np.add.at(pos, j[m], +shift)
                np.add.at(pos, i[m], -shift)
        for k in range(3):
            pos[:, k] = np.clip(pos[:, k], r, box[k] - r)
    return pos


def _add_fibers(pos, dia, typ, mol, box, spec, d_um, rng):
    axis, n_fib, n_bead = spec
    d_bead = 1.0 * d_um
    step = 0.9 * d_bead
    new_pos, new_typ, new_mol, new_dia = [], [], [], []
    for f in range(n_fib):
        if axis == "z":
            u = np.array([0.0, 0.0, 1.0])
        else:
            th = rng.uniform(0, 2 * np.pi)
            u = np.array([np.cos(th), np.sin(th), 0.0])
        L = step * (n_bead - 1)
        lo = 0.55 * d_bead + np.abs(u) * (L / 2.0)
        hi = box - lo
        p0 = np.where(hi > lo, rng.uniform(np.minimum(lo, hi),
                                          np.maximum(lo, hi)), box / 2.0)
        for b in range(n_bead):
            new_pos.append(p0 + (b - (n_bead - 1) / 2.0) * step * u)
            new_typ.append(2)
            new_mol.append(f + 1)
            new_dia.append(d_bead)
    if not new_pos:
        return pos, dia, typ, mol
    pos = np.vstack([pos, np.array(new_pos)])
    dia = np.concatenate([dia, np.array(new_dia)])
    typ = np.concatenate([typ, np.array(new_typ, int)])
    mol = np.concatenate([mol, np.array(new_mol, int)])
    for k in range(3):
        pos[:, k] = np.clip(pos[:, k], dia / 2, box[k] - dia / 2)
    return pos, dia, typ, mol


def make(family="phi_060", d_um=20.0, seed=0, n=N_TARGET, n_iter=200):
    cfgf = FAMILIES[family]
    rng = np.random.default_rng(seed * 1009 + (zlib.crc32(family.encode()) % 997))

    if cfgf.get("kind") == "jittered":
        ids, typ, pos, dia, mol = SP.make_structure(d_um, seed=seed)
        box = pos.max(axis=0) - pos.min(axis=0) + dia.max()
    else:
        dia = _diameters(n, d_um, cfgf, rng)
        vol = (np.pi / 6.0) * (dia ** 3).sum()
        L = (vol / cfgf["phi"]) ** (1.0 / 3.0)
        box = np.array([L, L, L])
        pos = rng.uniform(0.0, 1.0, size=(n, 3)) * box

        if cfgf.get("void", 0.0) > 0.0:
            c = box * 0.5
            Rv = cfgf["void"] * L
            v = pos - c
            dv = np.linalg.norm(v, axis=1)
            m = dv < Rv
            if m.any():
                pos[m] = c + v[m] / np.maximum(dv[m], 1e-9)[:, None] * Rv * 1.05
        if cfgf.get("cluster", 0.0) > 0.0:
            nk = 6
            cen = rng.uniform(0.15, 0.85, size=(nk, 3)) * box
            idx = rng.integers(0, nk, size=n)
            pos += cfgf["cluster"] * (cen[idx] - pos)
        if cfgf.get("strat", 0.0) > 0.0:
            t = pos[:, 2] / box[2]
            p = cfgf["strat"]
            pos[:, 2] = box[2] * t ** (1.0 + p)
        seg = cfgf.get("segregate", 0.0)
        seg = 1.0 if seg is True else float(seg or 0.0)
        if seg > 0.0:
            order = np.argsort(-dia)
            zs = np.sort(pos[:, 2])
            z_full = np.empty_like(pos[:, 2])
            z_full[order] = zs
            pos[:, 2] = (1.0 - seg) * pos[:, 2] + seg * z_full

        typ = np.ones(n, int)
        mol = np.zeros(n, int)
        if cfgf.get("coat_frac", 0.0) > 0.0:
            nc = int(round(cfgf["coat_frac"] * n))
            typ[rng.choice(n, nc, replace=False)] = 2
        pos = _relax(pos, dia, box, n_iter, grav=cfgf.get("grav", 0.0))
        if "fiber" in cfgf:
            pos, dia, typ, mol = _add_fibers(pos, dia, typ, mol, box,
                                             cfgf["fiber"], d_um, rng)
        ids = np.arange(1, len(pos) + 1)

    tree = cKDTree(pos)
    pr = tree.query_pairs(dia.max(), output_type="ndarray")
    ov = np.array([])
    if len(pr):
        dist = np.linalg.norm(pos[pr[:, 1]] - pos[pr[:, 0]], axis=1)
        need = (dia[pr[:, 0]] + dia[pr[:, 1]]) / 2.0
        ov = need - dist
        ov = ov[ov > 0]
    vol = (np.pi / 6.0) * (dia ** 3).sum()
    box_v = float(np.prod(box))
    info = dict(family=family, group=cfgf["group"], purpose=cfgf.get("purpose", ""),
                claim=cfgf.get("claim", ""), d_um=d_um, seed=seed, n=len(ids),
                phi_target=cfgf.get("phi"),
                phi=vol / box_v if box_v > 0 else np.nan,
                d_cv=float(dia.std() / dia.mean()),
                n_overlap=int(len(ov)),
                ov_med_over_d=float(np.median(ov) / d_um) if len(ov) else 0.0,
                ov_max_over_d=float(ov.max() / d_um) if len(ov) else 0.0,
                box_z=float(box[2]))
    return ids, typ, pos, dia, mol, info


def build(family="phi_060", d_um=20.0, seed=0, treatment="untreated",
          matrix="epoxy", filler="Al2O3", n=N_TARGET, **cfg_kw):
    ids, typ, pos, dia, mol, info = make(family, d_um, seed, n)
    t2 = "lowk_coat" if FAMILIES[family].get("coat_frac", 0.0) > 0 else treatment
    f2 = "AlN" if "fiber" in FAMILIES[family] else filler
    tab = MR.pair(filler, treatment, f2, t2, matrix)
    d = dict(gap_max_um=0.5 * d_um, eta_cut=1.0, slab_axis=2, slab_frac=0.06,
             viz=False)
    d.update(cfg_kw)
    cfg = tsf.AnalysisConfig(**d)
    G = tsf.build_thermal_network(ids, typ, pos, dia, cfg, table=tab, mols=mol)
    src, snk = tsf.identify_slabs(ids, pos, cfg)
    if not src or not snk:
        return None
    Gt, gv = tsf.attach_virtual_terminals(G, src, snk, cfg)
    br = {}
    for _, _, dd in G.edges(data=True):
        br[dd["branch"]] = br.get(dd["branch"], 0) + 1
    info.update(N=G.number_of_nodes(), E=G.number_of_edges(),
                n_contact=br.get("contact", 0), n_gap=br.get("gap", 0),
                n_bond=br.get("bond", 0),
                z_mean=2.0 * G.number_of_edges() / max(G.number_of_nodes(), 1),
                z_contact=2.0 * br.get("contact", 0) / max(len(ids), 1),
                frac_type2=float((typ == 2).mean()))
    return dict(G=G, Gt=Gt, sheaf=tsf.CellularSheaf(Gt), cfg=cfg, table=tab,
                info=info, ids=ids, typ=typ, pos=pos, dia=dia, mol=mol)


def ensemble(families=None, seeds=(0, 1, 2), d_um=20.0, **kw):
    for f in (families or FAMILY_ORDER):
        for s in seeds:
            st = build(f, d_um=d_um, seed=s, **kw)
            if st is not None:
                yield st


if __name__ == "__main__":
    print("=" * 104)
    print("structure family catalogue (ensemble.py)")
    print("=" * 104)
    print(f"{'group':>10} {'family':>14} {'phi_tgt':>7} {'phi_act':>7} {'z_con':>6} "
          f"{'N':>5} {'E':>6} {'CV(d)':>8} {'ov/d':>7} {'type2':>6} "
          f"{'G_eff [W/K]':>12}  purpose")
    rows = []
    for st in ensemble(seeds=(0,)):
        i = st["info"]
        G0 = tsf.effective_conductance(st["sheaf"], 0.0)
        rows.append((i, G0))
        pt = i["phi_target"]
        print(f"{i['group']:>10} {i['family']:>14} "
              f"{(f'{pt:.2f}' if pt else '-'):>7} {i['phi']:7.4f} "
              f"{i['z_contact']:6.2f} {i['N']:5d} {i['E']:6d} {i['d_cv']:8.3f} "
              f"{i['ov_med_over_d']:7.4f} {i['frac_type2']:6.2f} {G0:12.4e}  "
              f"{i.get('purpose', '')}")
    print()
    hp = [(i, g) for i, g in rows if (i["phi_target"] or 0) >= 0.50]
    print(f"-> {len(rows)} families in total, of which high loading (phi>=0.50): **{len(hp)} families**")
    ph = [i["phi"] for i, _ in hp]
    zz = [i["z_contact"] for i, _ in hp]
    gg = [g for _, g in hp]
    print(f"  spread of the high-loading region: phi {min(ph):.3f}-{max(ph):.3f}, "
          f"z_contact {min(zz):.2f}-{max(zz):.2f}, G_eff x{max(gg)/min(gg):.1f}")
