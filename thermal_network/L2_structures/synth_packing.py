# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L2_structures/synth_packing.py
#   sha256(src) : 76f73e182a0e413b
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 7 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical. In addition, 1 top-level definitions (and 2 module-level
# statements used only by them) that the figures and verification scripts of
# the paper do not use were omitted (they are listed by name in rb_manifest);
# the remaining definitions are unchanged.
# ---------------------------------------------------------------------------
# Synthetic packings that are geometrically similar in the particle diameter d
# (all lengths scale with d), with materials taken from materials_real.
import sys as _sys, pathlib as _pathlib
_LIB = _pathlib.Path(__file__).resolve().parent.parent
for _d in ("L0_engine", "L1_materials", "L2_structures"):
    _p = str(_LIB / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import numpy as np

import fluxfield as ff
import materials_real as MR
import thermal_sheaf_filtration as tsf
from thermal_sheaf_filtration import AnalysisConfig

D_REF = 4.0

D_REPRESENTATIVE_UM = 20.0


def make_structure(d_um=D_REPRESENTATIVE_UM, seed=42, spacing_ratio=1.0125,
                   jitter_ratio=0.0625, add_fiber=True):
    s = d_um / D_REF
    ids, typ, pos, dia, mol = ff.make_structure(
        seed=seed,
        Lx=70.0 * s, Ly=18.0 * s, Lz=70.0 * s,
        d=d_um,
        spacing=spacing_ratio * d_um,
        void_center=(35.0 * s, 45.0 * s), void_r=13.0 * s,
        jitter=jitter_ratio * d_um,
        add_fiber=False)
    if not add_fiber:
        return ids, typ, pos, dia, mol
    p0 = np.array([12.0, 9.0, 20.0]) * s
    p1 = np.array([20.0, 9.0, 52.0]) * s
    d_bead = 1.25 * d_um
    fib = np.array([p0 + t * (p1 - p0) for t in np.linspace(0, 1, 7)])
    pos = np.vstack([pos, fib])
    typ = np.concatenate([typ, np.full(len(fib), 2, int)])
    dia = np.concatenate([dia, np.full(len(fib), d_bead)])
    mol = np.concatenate([mol, np.ones(len(fib), int)])
    ids = np.arange(1, len(pos) + 1)
    return ids, typ, pos, dia, mol


def table(filler_sphere="Al2O3", filler_fiber="AlN", treatment="untreated",
          matrix="epoxy"):
    return MR.pair(filler_sphere, treatment, filler_fiber, treatment, matrix)


def cfg(d_um=D_REPRESENTATIVE_UM, gap_max_ratio=0.375, eta_cut=1.0,
        slab_frac=0.06, **kw):
    d = dict(gap_max_um=gap_max_ratio * d_um, eta_cut=eta_cut,
             slab_axis=2, slab_frac=slab_frac, viz=False)
    d.update(kw)
    return AnalysisConfig(**d)


def build(d_um=D_REPRESENTATIVE_UM, seed=42, treatment="untreated",
          matrix="epoxy", filler_sphere="Al2O3", filler_fiber="AlN", **cfg_kw):
    tab = table(filler_sphere, filler_fiber, treatment, matrix)
    c = cfg(d_um, **cfg_kw)
    ids, typ, pos, dia, mol = make_structure(d_um, seed=seed)
    G = tsf.build_thermal_network(ids, typ, pos, dia, c, table=tab, mols=mol)
    src, snk = tsf.identify_slabs(ids, pos, c)
    Gt, gv = tsf.attach_virtual_terminals(G, src, snk, c)
    return dict(G=G, Gt=Gt, sheaf=tsf.CellularSheaf(Gt), cfg=c, table=tab,
                d_um=d_um, g_virtual=gv, ids=ids, typ=typ, pos=pos, dia=dia,
                mol=mol)

