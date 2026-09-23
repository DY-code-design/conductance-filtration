# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L2_structures/lattice.py
#   sha256(src) : 7c7f06719b609a16
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 7 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical. In addition, 2 top-level definitions that the figures and
# verification scripts of the paper do not use were omitted (they are listed
# by name in rb_manifest); the remaining definitions are unchanged.
# ---------------------------------------------------------------------------
# Regular lattices (SC / BCC / FCC) and their exact homogenisation (Sec. 5.3,
# V4): for a Bravais lattice with one site per primitive cell the uniform-gradient
# solution is exact, k_ab = (1/(2 V_prim)) sum_b g_b (R_b)_a (R_b)_b, so the
# comparison with effective-medium theory contains no boundary error.
# maxwell_eucken and nan_emt give the reference effective-medium values.
import sys as _sys, pathlib as _pathlib
_LIB = _pathlib.Path(__file__).resolve().parent.parent
for _d in ("L0_engine", "L1_materials", "L2_structures"):
    _p = str(_LIB / _d)
    if _p not in _sys.path:
        _sys.path.insert(0, _p)

import itertools

import numpy as np

import thermal_sheaf_filtration as tsf

BASIS = {
    "SC":  np.array([[0.0, 0.0, 0.0]]),
    "BCC": np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
    "FCC": np.array([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5],
                     [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]),
}
NN_OVER_A = {"SC": 1.0, "BCC": np.sqrt(3) / 2.0, "FCC": 1.0 / np.sqrt(2)}
Z_COORD = {"SC": 6, "BCC": 8, "FCC": 12}
PHI_MAX = {k: len(BASIS[k]) * (np.pi / 6.0) * NN_OVER_A[k] ** 3
           for k in BASIS}
LATTICES = ["SC", "BCC", "FCC"]

_BASE_PHI = {
    "SC":  [0.10, 0.20, 0.30, 0.35, 0.40, 0.45, 0.50],
    "BCC": [0.10, 0.20, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65],
    "FCC": [0.10, 0.20, 0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70],
}
PHI_GRID = {k: sorted(_BASE_PHI[k] + [round(c * PHI_MAX[k], 4)
                                      for c in (0.98, 0.998, 1.01)])
            for k in BASIS}


def lattice_constant(kind, d_um, phi):
    nb = len(BASIS[kind])
    return d_um * (nb * np.pi / (6.0 * phi)) ** (1.0 / 3.0)


def neighbors(kind, a, r_cut_over_a=1.30):
    pts = []
    m = int(np.ceil(r_cut_over_a)) + 1
    rng = range(-m, m + 1)
    for i, j, k in itertools.product(rng, rng, rng):
        for b in BASIS[kind]:
            v = (np.array([i, j, k], float) + b) * a
            if np.linalg.norm(v) > 1e-12:
                pts.append(v)
    pts = np.array(pts)
    keep = np.linalg.norm(pts, axis=1) <= r_cut_over_a * a
    return pts[keep]


def r_cut_for_eta(eta, margin=1.08):
    return margin * (1.0 + 1.0 / (4.0 * (np.exp(eta) - 1.0)))


def homogenize(kind, d_um, phi, table, cfg, r_cut_over_a=None):
    a = lattice_constant(kind, d_um, phi)
    if r_cut_over_a is None:
        r_cut_over_a = r_cut_for_eta(cfg.eta_cut)
    disp = neighbors(kind, a, r_cut_over_a)
    r = d_um / 2.0
    V_prim = (a * 1e-6) ** 3 / len(BASIS[kind])

    K = np.zeros(3)
    n_bonds = 0
    branch = {}
    g_nn, shells = None, {}
    nn_d = NN_OVER_A[kind] * a
    for v in disp:
        dist = float(np.linalg.norm(v))
        gap = dist - d_um
        g, at = tsf.pair_conductance(gap, r, r, 1, 1, table, cfg)
        if g is None:
            continue
        n_bonds += 1
        br = at.get("branch", "?")
        branch[br] = branch.get(br, 0) + 1
        shells[round(dist / a, 4)] = shells.get(round(dist / a, 4), 0) + 1
        if abs(dist - nn_d) < 1e-9 * a:
            g_nn = g
        vm = v * 1e-6
        K += g * vm ** 2
    K /= 2.0 * V_prim
    return dict(kind=kind, phi=phi, d_um=d_um, a_um=a,
                gap_nn_um=nn_d - d_um, k_zz=K[2], k_xx=K[0], k_yy=K[1],
                n_bonds=n_bonds, n_shells=len(shells),
                branch=branch, g_nn=g_nn)


def slab_structure(kind, d_um, phi, n_cells=6):
    a = lattice_constant(kind, d_um, phi)
    pos = []
    for i, j, k in itertools.product(range(n_cells), repeat=3):
        for b in BASIS[kind]:
            pos.append((np.array([i, j, k], float) + b) * a)
    pos = np.array(pos)
    n = len(pos)
    ids = np.arange(1, n + 1)
    typ = np.ones(n, int)
    dia = np.full(n, float(d_um))
    mol = np.zeros(n, int)
    return ids, typ, pos, dia, mol, a


def maxwell_eucken(kappa, k_m, phi):
    return k_m * ((kappa + 2 * k_m + 2 * phi * (kappa - k_m)) /
                  (kappa + 2 * k_m - phi * (kappa - k_m)))


def nan_emt(kappa, k_m, phi, R_s_matrix, d_um):
    aN = R_s_matrix * k_m / (d_um * 1e-6 / 2.0)
    A = kappa * (1 + 2 * aN) + 2 * k_m
    B = kappa * (1 - aN) - k_m
    return k_m * (A + 2 * phi * B) / (A - phi * B), aN


if __name__ == "__main__":
    import materials_real as MR
    from thermal_sheaf_filtration import AnalysisConfig
    print("=== lattice geometry (lattice.py) ===")
    print(f"{'latt':>5} {'z':>3} {'a_nn/a':>8} {'phi_max':>8}")
    for k in LATTICES:
        print(f"{k:>5} {Z_COORD[k]:3d} {NN_OVER_A[k]:8.4f} {PHI_MAX[k]:8.4f}")
    print()
    print("=== homogenisation check (for SC, k_zz = g_nn/a) ===")
    tab = MR.single("AlN", "untreated", "epoxy")
    cfg = AnalysisConfig(gap_max_um=1e6, eta_cut=1e-9, t_min_um=1e-9)
    for k in LATTICES:
        r = homogenize(k, 20.0, 0.40, tab, cfg, r_cut_over_a=1.01 * NN_OVER_A[k])
        a_m = r["a_um"] * 1e-6
        pred = {"SC": r["g_nn"] / a_m, "BCC": 2 * r["g_nn"] / a_m,
                "FCC": 4 * r["g_nn"] / a_m}[k]
        print(f"{k:>5} nearest neighbours only: k_zz={r['k_zz']:.6e}  "
              f"closed form={pred:.6e}  relative error={abs(r['k_zz']/pred-1):.2e}  "
              f"isotropy |k_xx/k_zz-1|={abs(r['k_xx']/r['k_zz']-1):.2e}")
    print()
    print("=== check: Nan(alpha_N=0) == Maxwell-Eucken ===")
    for phi in (0.1, 0.3, 0.5):
        a, _ = nan_emt(170.0, 0.2, phi, 0.0, 20.0)
        b = maxwell_eucken(170.0, 0.2, phi)
        print(f"  phi={phi}: Nan={a:.10f}  ME={b:.10f}  diff={abs(a/b-1):.2e}")
