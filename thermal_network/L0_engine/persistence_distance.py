# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L0_engine/persistence_distance.py
#   sha256(src) : 73a9b769a5b9daac
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Bottleneck distance d_B(D, D') = min over partial matchings of max |p - q|_inf
# between two persistence diagrams (matching to the diagonal allowed), solved by
# bisection over the matching threshold. Bars shorter than TAU are discarded on both
# sides first and at most K_CAP longest bars are kept, so the stability test reads
# d_B <= |Delta theta|_inf + 2 TAU.
import numpy as np

TAU = 0.005
K_CAP = 400

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching


def _trim(D):
    if not len(D):
        return []
    per = np.abs(np.asarray([d - b for b, d in D], float))
    keep = np.where(per >= TAU)[0]
    if len(keep) > K_CAP:
        keep = keep[np.argsort(-per[keep])[:K_CAP]]
    return [D[k] for k in keep]


def _matchable(C, pers1, pers2, eps):
    n, m = C.shape
    N = n + m
    A = np.zeros((N, N), dtype=bool)
    A[:n, :m] = (C <= eps)
    if n:
        A[np.arange(n), m + np.arange(n)] = (pers1 / 2.0 <= eps)
    if m:
        A[n + np.arange(m), np.arange(m)] = (pers2 / 2.0 <= eps)
    A[n:, m:] = True
    M = maximum_bipartite_matching(csr_matrix(A.astype(np.int8)),
                                   perm_type="column")
    return bool((M >= 0).all())


def bottleneck(D1, D2):
    D1, D2 = _trim(list(D1)), _trim(list(D2))
    if not D1 and not D2:
        return 0.0
    P1 = np.array(D1, float).reshape(-1, 2)
    P2 = np.array(D2, float).reshape(-1, 2)
    pers1 = np.abs(P1[:, 1] - P1[:, 0]) if len(P1) else np.array([])
    pers2 = np.abs(P2[:, 1] - P2[:, 0]) if len(P2) else np.array([])
    if len(P1) == 0:
        return float(pers2.max() / 2.0)
    if len(P2) == 0:
        return float(pers1.max() / 2.0)
    C = np.maximum(np.abs(P1[:, None, 0] - P2[None, :, 0]),
                   np.abs(P1[:, None, 1] - P2[None, :, 1]))
    cand = np.unique(np.concatenate([C.ravel(), pers1 / 2.0, pers2 / 2.0,
                                     [0.0]]))
    lo, hi = 0, len(cand) - 1
    if not _matchable(C, pers1, pers2, cand[hi]):
        return float(cand[hi])
    while lo < hi:
        mid = (lo + hi) // 2
        if _matchable(C, pers1, pers2, cand[mid]):
            hi = mid
        else:
            lo = mid + 1
    return float(cand[lo])
