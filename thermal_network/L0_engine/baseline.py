# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/L0_engine/baseline.py
#   sha256(src) : 55df50f255e10545
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 24 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical. In addition, 1 module-level statement (a command-line block,
# imports or constants) used only by definitions that are not part of this
# record were omitted (listed in rb_manifest); the remaining code is
# unchanged.
# ---------------------------------------------------------------------------
# Frozen baseline: a fixed set of pair conductances and of quantities24 descriptors,
# stored with full precision in baseline.json. `--check` recomputes
# them and lists every value that moved. The first 16 hex digits of the sha256 of
# baseline.json are the fingerprint stamped into every result CSV (baseline_fp).
from __future__ import annotations

import sys as _sys
import pathlib as _pathlib

if __package__ in (None, ""):
    _LIB = _pathlib.Path(__file__).resolve().parents[1]
    for _d in ("L0_engine", "L1_materials", "L2_structures", "L3_metrics"):
        _p = str(_LIB / _d)
        if _p not in _sys.path:
            _sys.path.insert(0, _p)

import hashlib
import json
import math
import time


import materials_real as MR
import thermal_sheaf_filtration as tsf
import quantities24 as Q

JSON = _pathlib.Path(__file__).resolve().with_name("baseline.json")

PAIR_CASES = [
    ("SiO2",    "epoxy",        "untreated", +2.0,   5.0),
    ("SiO2",    "epoxy",        "untreated", +0.05,  5.0),
    ("SiO2",    "epoxy",        "untreated", -0.02,  5.0),
    ("Al2O3",   "epoxy",        "untreated", +2.0,   5.0),
    ("Al2O3",   "epoxy",        "treated",   +0.05,  5.0),
    ("AlN",     "epoxy",        "untreated", +2.0,   5.0),
    ("AlN",     "epoxy",        "untreated", +0.05,  5.0),
    ("AlN",     "epoxy",        "untreated", -0.02,  5.0),
    ("AlN",     "epoxy",        "treated",   +0.05,  5.0),
    ("AlN",     "epoxy",        "treated",   -0.02,  5.0),
    ("AlN",     "high_k_resin", "treated",   +0.05,  5.0),
    ("AlN",     "silicone",     "untreated", +2.0,   5.0),
    ("hBN",     "epoxy",        "untreated", +0.05,  5.0),
    ("ZnO",     "epoxy",        "untreated", +0.05,  5.0),
    ("diamond", "epoxy",        "treated",   +0.05,  5.0),
    ("AlN",     "epoxy",        "untreated", +0.05,  0.5),
    ("AlN",     "epoxy",        "untreated", +0.05, 25.0),
]

Q24_CASES = [
    ("phi_060", 20.0, 0, "SiO2", "epoxy",        "untreated"),
    ("phi_060", 20.0, 0, "AlN",  "epoxy",        "treated"),
    ("phi_060", 20.0, 0, "AlN",  "high_k_resin", "treated"),
    ("voided",  20.0, 0, "AlN",  "epoxy",        "untreated"),
]
Q24_CFG = dict(skip_layer5=True, keep_raw=False)

TOL_SAME = 1e-12
TOL_NOISE = 1e-6


def _pair_row(fil, mat, tr, gap, rt):
    tab = MR.pair(fil, tr, fil, tr, mat)
    m = tab.materials[1]
    ri, rj = 2 * rt, 2 * rt
    cfg = tsf.AnalysisConfig(gap_max_um=10.0 * rt, eta_cut=1.0)
    g, info = tsf.pair_conductance(gap, ri, rj, 1, 1, table=tab, cfg=cfg)
    if g is None:
        return dict(case=f"{fil}/{mat}/{tr}/gap={gap:+g}/rt={rt:g}",
                    kappa=float(m.kappa), k_matrix=float(tab.k_matrix),
                    R_s_matrix=float(m.R_s_matrix), g=float("nan"),
                    branch="rejected", bob_saturated=0)
    row = dict(case=f"{fil}/{mat}/{tr}/gap={gap:+g}/rt={rt:g}",
               kappa=float(m.kappa), k_matrix=float(tab.k_matrix),
               R_s_matrix=float(m.R_s_matrix), g=float(g))
    for k in ("R_con", "R_int", "R_gap", "R_bulk", "R", "A_th_m2",
              "rho_um", "rho_spot_um", "k_gap_eff", "bob_h_shift_um",
              "bob_alpha", "a_um", "delta_um"):
        if k in info and isinstance(info[k], (int, float)):
            row[k] = float(info[k])
    row["branch"] = info.get("branch", "")
    row["bob_saturated"] = int(bool(info.get("bob_saturated", False)))
    return row


def measure() -> dict:
    t0 = time.time()
    out = {"stamp": Q.stamp(), "pair": [], "q24": []}
    for c in PAIR_CASES:
        out["pair"].append(_pair_row(*c))
    import ensemble as EN
    cfg = Q.Q24Config(**Q24_CFG)
    for fam, d, seed, fil, mat, tr in Q24_CASES:
        geo = EN.make(fam, d, seed)
        r = Q.compute(*geo[:5], info=geo[5], filler=fil, matrix=mat, treat=tr,
                      d_um=d, cfg=cfg, label=f"{fam}/{fil}/{mat}/{tr}")
        row = {k: v for k, v in r.items()
               if not k.startswith("_") and isinstance(v, (int, float))
               and not isinstance(v, bool)}
        row["case"] = f"{fam}/seed{seed}/{fil}/{mat}/{tr}"
        out["q24"].append(row)
    out["elapsed_s"] = round(time.time() - t0, 2)
    return out


def fingerprint() -> str:
    if not JSON.is_file():
        return "(not frozen)"
    return hashlib.sha256(JSON.read_bytes()).hexdigest()[:16]


def _index(d):
    return ({r["case"]: r for r in d["pair"]},
            {r["case"]: r for r in d["q24"]})


def check(verbose=True) -> int:
    if not JSON.is_file():
        print(f"no frozen file: {JSON}")
        print("  run `baseline.py --freeze` first")
        return 2
    old = json.loads(JSON.read_text(encoding="utf-8"))
    new = measure()
    op, oq = _index(old)
    np_, nq = _index(new)

    print("=" * 96)
    print("frozen baseline diff")
    print(f"  frozen : engine={old['stamp']['engine_version']} "
          f"schema={old['stamp'].get('q24_schema_version')}")
    print(f"  now    : engine={new['stamp']['engine_version']} "
          f"schema={new['stamp'].get('q24_schema_version')}"
          f"  (measured in {new['elapsed_s']} s)")
    print(f"  fingerprint: {fingerprint()}  <- check that it matches the recorded fixed value")
    print("=" * 96)

    moved, noise, same, missing = [], 0, 0, []
    for tag, o, n in (("pair conductance", op, np_), ("quantities", oq, nq)):
        for case, orow in o.items():
            if case not in n:
                missing.append(f"{tag} {case} (case disappeared)")
                continue
            for k, ov in orow.items():
                if k == "case" or not isinstance(ov, (int, float)):
                    continue
                if k not in n[case]:
                    missing.append(f"{tag} {case}: {k} (column disappeared)")
                    continue
                nv = n[case][k]
                if not isinstance(nv, (int, float)):
                    continue
                if (isinstance(ov, float) and math.isnan(ov)
                        and isinstance(nv, float) and math.isnan(nv)):
                    same += 1
                    continue
                rel = abs(nv - ov) / max(abs(ov), 1e-300)
                if rel <= TOL_SAME:
                    same += 1
                elif rel <= TOL_NOISE:
                    noise += 1
                else:
                    moved.append((tag, case, k, ov, nv,
                                  math.log10(max(rel, 1e-300))))
        for case in n:
            if case not in o:
                missing.append(f"{tag} {case} (case added; --freeze is required)")

    print(f"identical {same} / numerical noise (<1e-6) {noise} / "
          f"**moved {len(moved)}** / missing or added {len(missing)}")
    if missing:
        print("\n[column/case mismatches]")
        for m in missing[:20]:
            print(f"  ! {m}")
    if moved:
        moved.sort(key=lambda x: -x[5])
        print(f"\n{'kind':<10}{'quantity':<22}{'old':>14}{'new':>14}{'ratio':>10}  case")
        print("-" * 96)
        for tag, case, k, ov, nv, _ in moved:
            ratio = (nv / ov) if ov not in (0,) else float("inf")
            print(f"{tag:<10}{k:<22}{ov:>14.6g}{nv:>14.6g}{ratio:>10.4f}  {case}")
        print("-" * 96)
        keys = sorted({m[2] for m in moved})
        print(f"moved quantities: {', '.join(keys)}")
        print("\nnext: check that this list of moved quantities is what you expected;")
        print("      if it is wider than expected, stop here before re-running anything.")
    else:
        print("\nno difference. The defaults have not changed (no re-run needed).")
    return 1 if (moved or missing) else 0


def freeze():
    new = measure()
    JSON.write_text(json.dumps(new, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"frozen: {JSON}")
    print(f"  fingerprint: {fingerprint()}  <- copy this into the recorded fixed value")
    print(f"  engine={new['stamp']['engine_version']} "
          f"schema={new['stamp']['q24_schema_version']} "
          f"/ pair {len(new['pair'])} cases / quantities {len(new['q24'])} cases "
          f"/ {new['elapsed_s']} s")


def show():
    d = json.loads(JSON.read_text(encoding="utf-8"))
    print(f"engine={d['stamp']['engine_version']} "
          f"schema={d['stamp']['q24_schema_version']}")
    print("\n--- pair conductance (g [W/K]) ---")
    for r in d["pair"]:
        print(f"  {r['g']:.10e}  {r['case']}")
    print("\n--- quantities (main columns only) ---")
    ks = ["k_eff", "gain_interface", "gain_matrix", "gain_filler",
          "theta_star", "theta_allconn", "R_decades", "share_contact"]
    print("  " + "".join(f"{k:>16}" for k in ks) + "  case")
    for r in d["q24"]:
        print("  " + "".join(f"{r.get(k, float('nan')):>16.6g}" for k in ks)
              + f"  {r['case']}")


def _guard_variant() -> None:
    import os as _os
    v = tsf.active_variant()
    if _os.environ.get("THERMAL_VARIANT") in (None, "", tsf._DEFAULT_VARIANT):
        return
    if "--allow-variant" in _sys.argv:
        print(f"WARNING: running with variant '{v}' active (--allow-variant given).")
        print("  this result **cannot be compared with the frozen baseline**.")
        return
    print("=" * 96)
    print(f"x aborted: variant '{v}' is active (environment variable THERMAL_VARIANT).")
    print("  comparing the frozen baseline under a variant would report no difference")
    print("  while computing a different model. Unset the variant and run again:")
    print(f"    unset THERMAL_VARIANT && python3 {_pathlib.Path(__file__).name} --check")
    print("  if you intend to run with the variant, add --allow-variant (the comparison result is then invalid).")
    print("=" * 96)
    raise SystemExit(3)


if __name__ == "__main__":
    _guard_variant()
    if "--freeze" in _sys.argv:
        freeze()
    elif "--show" in _sys.argv:
        show()
    else:
        raise SystemExit(check())
