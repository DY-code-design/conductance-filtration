# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : figsrc/_paths.py
#   sha256(src) : 90adfb3fffc0e89f
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source.
# ---------------------------------------------------------------------------
# Path resolution: data/ holds the result files, thermal_network/ the package,
# fig/ takes the output. Missing directories are created.
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent

BUNDLED_DATA = HERE.parent / "data"

POC_LAYER = "B_network_layer"


def _find_lib():
    cands = [HERE.parent / "thermal_network"]
    cands += [p / "03_poc" / "lib" for p in HERE.parents]
    for c in cands:
        if (c / "_bootstrap.py").is_file():
            return c
    raise FileNotFoundError(
        "cannot find the package directory (the one holding _bootstrap.py); "
        "in this record it is thermal_network/, next to figures/: "
        + ", ".join(str(c) for c in cands))


LIB = _find_lib()


def add_lib():
    s = str(LIB)
    if s not in sys.path:
        sys.path.insert(0, s)
    return LIB


def poc_csv(subdir, name):
    p = BUNDLED_DATA / name
    if p.is_file():
        return p
    p = LIB.parent / POC_LAYER / subdir / name
    if p.is_file():
        return p
    raise FileNotFoundError(
        f"cannot find {name} (tried {BUNDLED_DATA / name} and"
        f" {LIB.parent / POC_LAYER / subdir / name})")


def out_dir():
    p = HERE.parent / "manuscript" / "fig"
    if p.is_dir():
        return p
    p = HERE.parent / "fig"
    p.mkdir(parents=True, exist_ok=True)
    return p


def preview_dir():
    p = HERE / "preview"
    p.mkdir(parents=True, exist_ok=True)
    return p
