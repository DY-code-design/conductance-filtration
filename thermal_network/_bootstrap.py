# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/lib/_bootstrap.py
#   sha256(src) : a7ce507bd560cf14
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 3 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical. In addition, 1 top-level definitions (and 1 module-level
# statement used only by them) that the figures and verification scripts of
# the paper do not use were omitted (they are listed by name in rb_manifest);
# the remaining definitions are unchanged.
# ---------------------------------------------------------------------------
# Import this first in every script: it puts the four layers of thermal_network/
# on sys.path and moves the working directory to the script's own folder, so
# that results are written next to the script.
import os
import pathlib
import sys

_LIB = pathlib.Path(__file__).resolve().parent
LAYERS = ("L0_engine", "L1_materials", "L2_structures", "L3_metrics")
for _d in LAYERS:
    _p = str(_LIB / _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

_caller = None
for _f in sys._current_frames().values():
    pass
import inspect
for _fi in inspect.stack():
    _fn = pathlib.Path(_fi.filename).resolve()
    if _fn != pathlib.Path(__file__).resolve() and _fn.suffix == ".py":
        _caller = _fn.parent
        break
if _caller is not None and _caller.is_dir():
    os.chdir(_caller)

ROOT = _LIB.parent


from _constants import N_RANDOM_CONTROL


def _install_stamp_hook():
    try:
        import pandas as _pd
    except Exception:
        return
    if getattr(_pd.DataFrame.to_csv, "_stamped", False):
        return
    _orig = _pd.DataFrame.to_csv

    def to_csv(self, path_or_buf=None, *a, **kw):
        try:
            name = pathlib.Path(str(path_or_buf)).name
            if (path_or_buf is not None and name.startswith("results_")
                    and name.endswith(".csv")):
                import thermal_sheaf_filtration as _tsf
                st = _tsf.engine_stamp()
                if "engine_variant" not in self.columns:
                    self = self.assign(**st)
                elif self["engine_variant"].isna().any():
                    self = self.copy()
                    for k, v in st.items():
                        if k in self.columns:
                            self[k] = self[k].fillna(v)
                        else:
                            self[k] = v
        except Exception:
            pass

        try:
            _p = pathlib.Path(str(path_or_buf))
            _atomic = (path_or_buf is not None
                       and isinstance(path_or_buf, (str, pathlib.Path))
                       and _p.suffix == ".csv"
                       and kw.get("mode", "w") == "w")
        except Exception:
            _atomic = False
        if _atomic:
            _tmp = _p.with_name(_p.name + f".part{os.getpid()}")
            try:
                r = _orig(self, _tmp, *a, **kw)
                os.replace(_tmp, _p)
                return r
            except BaseException:
                try:
                    _tmp.unlink()
                except Exception:
                    pass
                raise
        return _orig(self, path_or_buf, *a, **kw)

    to_csv._stamped = True
    _pd.DataFrame.to_csv = to_csv


_install_stamp_hook()


def tee_stdout(path):
    p = pathlib.Path(path)

    class _Tee:
        def __init__(self, *streams):
            self._streams = streams

        def write(self, s):
            for st in self._streams:
                try:
                    st.write(s)
                except ValueError:
                    pass
            return len(s)

        def flush(self):
            for st in self._streams:
                try:
                    st.flush()
                except ValueError:
                    pass

    fh = open(p, "w", encoding="utf-8")
    sys.stdout = _Tee(sys.__stdout__, fh)

    def _close():
        try:
            fh.flush(); fh.close()
        except Exception:
            pass
        sys.stdout = sys.__stdout__

    import atexit
    atexit.register(_close)
    return p


def poc_csv(name, *hints):
    cand = [ROOT / "data" / name]
    if hints:
        cand.append(pathlib.Path(*hints) / name)
    for c in cand:
        if c.is_file():
            return c.resolve()
    for c in ROOT.rglob(name):
        if "__pycache__" not in c.parts:
            return c.resolve()
    raise FileNotFoundError(f"{name} not found (searched in: {ROOT})")


def add_case_path(case):
    here = pathlib.Path.cwd()
    for c in (here.parent / case, here):
        if c.is_dir():
            p = str(c.resolve())
            if p not in sys.path:
                sys.path.insert(0, p)
            return c.resolve()
    raise FileNotFoundError(f"folder for case {case} not found")

