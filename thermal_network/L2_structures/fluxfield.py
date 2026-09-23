# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L2_structures/fluxfield.py
#   sha256(src) : 57c0104a92839f24
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 1 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical. In addition, 4 top-level definitions (and 7 module-level
# statements used only by them) that the figures and verification scripts of
# the paper do not use were omitted (they are listed by name in rb_manifest);
# the remaining definitions are unchanged.
# ---------------------------------------------------------------------------
# make_structure: synthetic packing on a jittered lattice (used by the other
# generators).

import sys as _sys, pathlib as _pathlib
_LIB = _pathlib.Path(__file__).resolve().parent.parent
for _d in ("L0_engine", "L1_materials", "L2_structures"):
    _p = str(_LIB / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import numpy as np


def make_structure(seed=42, Lx=70.0, Ly=18.0, Lz=70.0, d=4.0, spacing=4.4,
                   void_center=(35.0, 45.0), void_r=13.0, jitter=0.35,
                   add_fiber=True):
    rng = np.random.default_rng(seed)
    pos, typ, dia, mol = [], [], [], []
    nx_, ny_, nz_ = int(Lx / spacing), int(Ly / spacing), int(Lz / spacing)
    for ix in range(nx_):
        for iy in range(ny_):
            for iz in range(nz_):
                p = np.array([ix, iy, iz]) * spacing + spacing / 2
                p = p + rng.normal(0, jitter, 3)
                if (p[0] - void_center[0]) ** 2 + (p[2] - void_center[1]) ** 2 < void_r ** 2:
                    continue
                pos.append(p); typ.append(1); dia.append(d); mol.append(0)

    if add_fiber:
        p0 = np.array([12.0, Ly / 2, 20.0]); p1 = np.array([20.0, Ly / 2, 52.0])
        for t in np.linspace(0, 1, 7):
            pos.append(p0 + t * (p1 - p0)); typ.append(2); dia.append(5.0); mol.append(1)

    pos = np.array(pos); typ = np.array(typ, int)
    dia = np.array(dia); mol = np.array(mol, int)
    ids = np.arange(1, len(pos) + 1)
    return ids, typ, pos, dia, mol

