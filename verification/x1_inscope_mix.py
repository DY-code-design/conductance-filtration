# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/B_network_layer/X1_inscope_mix/x1_inscope_mix.py
#   sha256(src) : 0dbbaeb82d7d668f
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 18 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V18. Edge order under treatment patchiness: on one monodisperse packing, 30% of
# the particles untreated and 70% treated; Kendall's tau between the centre-distance
# order and the conductance order of the edges, and the efficiency of a surrogate
# that ignores the patchiness.
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

import time

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

import ensemble as EN
import materials_real as MR
import thermal_sheaf_filtration as tsf
from thermal_sheaf_filtration import Material, MaterialTable, AnalysisConfig

FAMILY, D_UM = "phi_060", 20.0
SEEDS = (0, 1)
FILLER, MATRIX = "AlN", "epoxy"
MINOR, MAJOR = "untreated", "treated"
FRAC_MINOR = 0.30
K_TOP = 20
BOOST = 10.0

OUT_CSV = "results_x1_mix.csv"
OUT_TXT = "results_x1_mix.txt"


def table_mixed():
    def mat(name, treat):
        t = MR.TREATMENT[treat]
        return Material(name, kappa=MR.KAPPA[FILLER],
                        R_s_contact=t["R_s_contact"], R_s_matrix=t["R_s_matrix"])
    return MaterialTable(materials={1: mat("major", MAJOR), 2: mat("minor", MINOR)},
                         k_matrix=MR.K_MATRIX[MATRIX])


def table_uniform(treat):
    t = MR.TREATMENT[treat]
    m = Material("uni", kappa=MR.KAPPA[FILLER],
                 R_s_contact=t["R_s_contact"], R_s_matrix=t["R_s_matrix"])
    return MaterialTable(materials={1: m, 2: m}, k_matrix=MR.K_MATRIX[MATRIX])


def build(ids, typ, pos, dia, mol, tab):
    cfg = AnalysisConfig(gap_max_um=0.5 * D_UM, eta_cut=1.0, slab_axis=2,
                         slab_frac=0.06, viz=False)
    G = tsf.build_thermal_network(ids, typ, pos, dia, cfg, table=tab, mols=mol)
    src, snk = tsf.identify_slabs(ids, pos, cfg)
    Gt, _ = tsf.attach_virtual_terminals(G, src, snk, cfg)
    return G, Gt, tsf.CellularSheaf(Gt)


def gain(Gt, sh, pairs):
    G0 = tsf.effective_conductance(sh, 0.0)
    saved = []
    for u, v in pairs:
        saved.append((u, v, Gt[u][v]["g"])); Gt[u][v]["g"] *= BOOST
    g = tsf.effective_conductance(sh, 0.0)
    for u, v, x in saved:
        Gt[u][v]["g"] = x
    return g / G0 - 1.0


def main() -> None:
    log = []

    def say(s=""):
        print(s); log.append(s)

    t0 = time.time()
    say("=== X1: in-scope material mixing (t_ip=0, measurement-derived R_s only) ===")
    ru = MR.TREATMENT[MINOR]["R_s_matrix"]; rt = MR.TREATMENT[MAJOR]["R_s_matrix"]
    say(f"structure {FAMILY} (monodisperse) / {FILLER}/{MATRIX} / "
        f"{int(FRAC_MINOR*100)}% {MINOR} + {int((1-FRAC_MINOR)*100)}% {MAJOR}")
    say(f"R_s_matrix: {MINOR} {ru:.3e} / {MAJOR} {rt:.3e} -> **ratio {ru/rt:.3f}x**"
        f" (`lowk_coat` was {MR.TREATMENT['lowk_coat']['R_s_matrix']/rt:.1f}x)")
    say()

    rows = []
    for seed in SEEDS:
        ids, typ, pos, dia, mol, info = EN.make(FAMILY, D_UM, seed)
        ids = np.asarray(ids); pos = np.asarray(pos); dia = np.asarray(dia)
        rng = np.random.default_rng(20260813 + seed)
        n = len(ids)
        minor_idx = rng.choice(n, int(round(FRAC_MINOR * n)), replace=False)
        typ_mix = np.ones(n, dtype=int)
        typ_mix[minor_idx] = 2

        G, Gt, sh = build(ids, typ_mix, pos.copy(), dia, mol, table_mixed())
        idx = {int(v): k for k, v in enumerate(ids)}

        w_D, w_C = [], []
        for u, v, d in G.edges(data=True):
            i, j = idx[int(u)], idx[int(v)]
            w_D.append(float(np.linalg.norm(pos[i] - pos[j])))
            w_C.append(-float(np.log10(d["g"])))
        tau = float(kendalltau(w_D, w_C).statistic)
        rho = float(spearmanr(w_D, w_C).statistic)

        flux = tsf.edge_flux_table(sh, 0.0)
        real = [(int(r.u), int(r.v)) for r in flux.itertuples()]
        true_top = real[:K_TOP]
        g_true = gain(Gt, sh, true_top)
        eta = {}
        for uni in (MAJOR, MINOR):
            Gu, Gtu, shu = build(ids, np.ones(n, dtype=int), pos.copy(), dia,
                                 mol, table_uniform(uni))
            fu = tsf.edge_flux_table(shu, 0.0)
            sur_top = [(int(r.u), int(r.v)) for r in fu.itertuples()][:K_TOP]
            g_sur = gain(Gt, sh, [p for p in sur_top
                                  if Gt.has_edge(*p)])
            eta[uni] = g_sur / g_true if g_true > 0 else np.nan

        g_rnd = float(np.mean([gain(Gt, sh, [real[i] for i in
                                             rng.choice(len(real), K_TOP, replace=False)])
                               for _ in range(_bootstrap.N_RANDOM_CONTROL)])) / g_true

        rows.append(dict(seed=seed, family=FAMILY, filler=FILLER, matrix=MATRIX,
                         frac_minor=FRAC_MINOR, minor=MINOR, major=MAJOR,
                         d_cv=float(info.get("d_cv", np.nan)),
                         phi=float(info["phi"]), n_edge=G.number_of_edges(),
                         Rs_ratio=ru / rt, tau_G2_C=tau, rho_G2_C=rho,
                         eta_sur_major=eta[MAJOR], eta_sur_minor=eta[MINOR],
                         eta_random=g_rnd, gain_true=g_true))
        say(f"  seed{seed}: CV={info.get('d_cv', float('nan')):.3f} edges={G.number_of_edges()} "
            f"tau={tau:+.4f}  eta_sur({MAJOR})={eta[MAJOR]:.4f}  "
            f"eta_sur({MINOR})={eta[MINOR]:.4f}  eta_random={g_rnd:.4f}")

    df = pd.DataFrame(rows).assign(**tsf.engine_stamp())
    df.to_csv(OUT_CSV, index=False)
    say()
    say("--- summary (seed average) ---")
    say(f"  CV of particle size          : {df.d_cv.mean():.6f} (should be 0 for monodisperse)")
    say(f"  tau(D, theta)                : **{df.tau_G2_C.mean():+.4f}**"
        f" (per seed {', '.join(f'{x:+.4f}' for x in df.tau_G2_C)})")
    say(f"  eta_sur (all treated as {MAJOR}): **{df.eta_sur_major.mean():.4f}**")
    say(f"  eta_sur (all treated as {MINOR}): **{df.eta_sur_minor.mean():.4f}**")
    say(f"  random k-edge control        : {df.eta_random.mean():.4f}")
    say()
    say(f"WARNING: `coat_mix` (out of scope) gave tau=0.684 / eta_sur=0.821."
        f" The R_s ratio differs: {ru/rt:.2f}x vs 26.3x.")
    say(f"total run time {time.time()-t0:.0f}s")
    _pathlib.Path(OUT_TXT).write_text("\n".join(log) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
