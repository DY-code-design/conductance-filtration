# The conductance filtration — software record

Implementation, figure scripts and result files for the paper

> **The conductance filtration: sink-relative persistence and exact edge
> sensitivity for heat transport in particle-filled composites**
> Daisuke Yasufuku (U-MAP Co., Ltd.)

This record is the software that reproduces every figure and every verification
value of that paper. It is released on its own so that it can be cited and
versioned independently of the paper. The paper's ancillary files carry only
what the text points at directly (the verification values of Table S.1 and the
definition table of the structure families); everything runnable is here.

Two limits are worth stating at the top. First, seven of the result files under
`data/` are distributed as data, and the scripts that produced them are not part
of this record: `results_f4_meta.csv`, `results_f4_theta_sweep.csv`,
`results_s2_scaling.csv`, `results_ensemble_theta.csv`, `results_y7_holdout.csv`,
`results_d1_demo.csv` and `results_d1_edges.csv`. Figures 9, 11, 13 and 14 are
redrawn from them; the numbers in them are not recomputed here. The other three
(`results_t1.csv`, `results_poc5_grid.csv`, `results_ensemble.csv`) are written
by scripts in this record. Second, definitions the paper does not use have been
left out of the distributed sources; the next section says exactly what that
means.

## Citing

```
This record        : https://doi.org/10.5281/zenodo.22898126
Paper (Japanese)   : https://doi.org/10.5281/zenodo.22760047
```

`CITATION.cff` carries the same values in machine-readable form.

The two DOIs are of different kinds. The one for the paper is a *concept* DOI:
it always resolves to the latest version of that record. The one for this
record is a *version* DOI: it points at this release and will keep pointing at
it after later releases. If you want the one that always resolves to the latest
version of this software, use the concept DOI shown on its Zenodo record page.
An English version of the paper is to appear on arXiv.

## Quick start

Requires Python 3.10 or later and the packages in `requirements.txt`
(numpy, scipy, networkx, pandas, matplotlib). Run everything from the root of
this record; nothing needs to be installed or put on `PYTHONPATH`.

```bash
pip install -r requirements.txt

# 1. the engine reproduces the frozen baseline
python3 thermal_network/L0_engine/baseline.py --check
#    -> identical 496 / numerical noise (<1e-6) 0 / **moved 0** / missing or added 0
#    -> no difference. The defaults have not changed (no re-run needed).

# 2. the run-time self-test
python3 thermal_network/L0_engine/selftest.py
#    -> ALL SELFTEST PASS

# 3. reproduce one figure
python3 figures/mkfig7.py
#    -> wrote ./fig/fig7_material_only.pdf

# 4. recompute one verification item
python3 verification/poc4a_identities.py
```

Figures are written to `fig/` and preview images to `figures/preview/`; both are
created on first use. Check 1 is the one to run first: if the fingerprint does
not reproduce, the numbers below are not the numbers in the paper.

**Tested environment.** Every step in this README was run to completion, and
every value of Table S.1 was reproduced, with Python 3.10 on Ubuntu 22.04.
Other operating systems and Python versions have not been verified.

## Layout

| Path | What it is |
|---|---|
| `thermal_network/` | The implementation. `L0_engine` (conductances, filtration, persistence, the frozen baseline), `L1_materials`, `L2_structures` (packing generation, lattices, flux field), `L3_metrics` |
| `figures/` | One script per figure of the paper, plus four shared helpers (`_paths`, `_pdffont`, `_branchfig`, `_data`) |
| `data/` | The result files the figure scripts read, flat, one directory. Seven of them are data only (see the top of this file) |
| `verification/` | The 17 scripts behind items V4–V20 of Table S.1 |
| `requirements.txt` | The packages the authors ran |
| `rb_manifest.json`, `rb_manifest.md` | Machine-generated inventory: the sha256 and the working-repository origin of every file here |
| `CITATION.cff`, `LICENSE` | Citation metadata; MIT |

Each `.py` file carries only a short explanation written for this record; the
working repository's development notes are not included. Each file's header
gives the file it was extracted from and the sha256 of that source.

Three things were done to the sources on the way out, and each file's header
says which of them applies.

1. **Comments and docstrings are removed.** With docstrings removed, the two
   abstract syntax trees are identical.
2. **Reader-facing strings were translated from Japanese.** In the files that
   print a report, the guarantee is that with every string constant masked the
   two trees are identical, so the logic is unchanged. The manifest lists, for
   each file, how many strings were translated.
3. **Top-level definitions the paper does not use were left out.** These are the
   end-to-end pipeline for LAMMPS dumps and its plotting helpers, the Phase-2
   work (jet sheaf, Kron reduction, cluster and component sweeps), the BET and
   Hertz helpers behind the disabled roughness model, and the tables and imports
   that only those use. `rb_manifest.json` lists every one of them by name, with
   its length and the reason it was left out, under `omitted`.

Taken together, the guarantee is this: **with docstrings removed, string
constants masked, and the listed definitions excluded, the syntax tree of every
file here is identical to the working copy it came from.** Nothing that any
figure or any verification item reaches was removed — the checks below were run
against these sources, not against the working copy.

One consequence is visible in the self-test: this record runs 43 checks, and
prints one `[SKIP]` line for the Phase-2 checks, which are not part of it. The
working repository runs more, because it still has the Phase-2 code.

One thing is left as it is on purpose: the `engine_note` field inside
`thermal_network/L0_engine/baseline.json` is a development record written in
Japanese. The fingerprint is the sha256 of that whole file, so editing one
character of it would change the fingerprint and break the check against the
`baseline_fp` column of every result file.

## Figures

Each script writes the PDF that appears in the paper under that number.

| Figure | Script | Data it reads |
|---|---|---|
| 1 | `figures/mkfig1.py` | — (concept figure) |
| 2 | `figures/mkfig2.py` | — |
| 3 | `figures/mkfig3.py` | — |
| 4 | `figures/mkfig4.py` | — (computed on the fly) |
| 5 | `figures/mkfig5.py` | — (computed on the fly) |
| 6 | `figures/mkfig6.py` | — (computed on the fly) |
| 7 | `figures/mkfig7.py` | `results_t1.csv` |
| 8 | `figures/mkfig8.py` | — (concept figure) |
| 9 | `figures/mkfig9.py` | `results_f4_meta.csv`, `results_f4_theta_sweep.csv` |
| 10 | `figures/mkfig10.py` | — (pseudocode) |
| 11 | `figures/mkfig11.py` | `results_s2_scaling.csv` |
| 12 | `figures/mkfig12.py` | `results_poc5_grid.csv` |
| 13 | `figures/mkfig13.py` | `results_ensemble_theta.csv`, `results_y7_holdout.csv` |
| 14 | `figures/mkfig14.py` | `results_d1_demo.csv`, `results_d1_edges.csv` |
| B.1 | `figures/mkfig_b1.py` | — |

`data/family_definitions.csv` is read by no figure: it is the source of Table
C.2 of the paper, and it is also distributed with the paper's ancillary files.

## Table S.1 — what recomputes each item

The values these commands print are listed in `verification.md` of the paper's
ancillary files. Items V1–V3 have no script of their own: they are identities
checked by the package at run time, among the 43 checks that
`thermal_network/L0_engine/selftest.py` runs. Their lines are

```
[PASS] V1 closed form of the near-contact branch vs numerical annular integration (5 values of h)
[PASS] T10 (V2) at h=0 and a_K=0 the shift equals the legacy BOB cap (analytic identity)
[PASS] T17u (V3) the switch is continuous (delta->0+ and h->0+ agree)
```

| Item | Command |
|---|---|
| V1 | `python3 thermal_network/L0_engine/selftest.py` |
| V2 | `python3 thermal_network/L0_engine/selftest.py` |
| V3 | `python3 thermal_network/L0_engine/selftest.py` |
| V4 | `python3 verification/poc5_lattice_emt.py` |
| V5 | `python3 verification/v5_explicit_vertex.py` |
| V6 | `python3 verification/size_sweep.py` |
| V7 | `python3 verification/poc4a_identities.py` |
| V8 | `python3 verification/v89_log_sensitivity.py` |
| V9 | `python3 verification/v89_log_sensitivity.py` |
| V10 | `python3 verification/poc1_predictive_power.py` |
| V11 | `python3 verification/a1_closed_form.py` |
| V12 | `python3 verification/poc1_predictive_power.py` |
| V13 | `python3 verification/poc6_robustness.py` |
| V14 | `python3 verification/r48_ideal_saturation.py` |
| V15 | `python3 verification/y4_theta_separation.py` |
| V16 | `python3 verification/y4_theta_separation.py` |
| V17 | `python3 verification/m7d_Rin_Rout.py`, `python3 verification/m7c_isothermal_test.py` |
| V18 | `python3 verification/x1_inscope_mix.py` |
| V19 | `python3 verification/poc3_stability.py`, `python3 verification/poc3c_fixed_complex.py` |
| V20 | `python3 verification/t3_pmat_bound.py` |

One further script is included because other checks call it:
`verification/t1_filtration_compare.py`.

Running the whole of `verification/` takes 15 to 25 minutes. Each script writes
its own `results_*` files next to itself, inside `verification/`. Three of them
— `a1_closed_form.py`, `v89_log_sensitivity.py` and `y4_theta_separation.py` —
keep their progress as they go, so if one is interrupted, running it again
resumes rather than starting over.

Two things in the printed output are not verification results. A line reading
`a-priori prediction: met / not met` compares the outcome with what the authors
predicted before running it; it is a record of the prediction, not a pass or
fail of the check. Labels such as `PoC-N` and `T3-1` are the authors' internal
names for the experiments, kept so that the output can be traced back to the
working notes.

## What is fixed, and what is regenerated

Every packed structure used in the paper is synthetic; no measured structure is
included. The structures themselves are not distributed: the generation rules
and the random seeds are part of the implementation, so the same population is
regenerated on the fly. Random seeds are fixed, and every result file records
the parameters, the version of the implementation that computed the
conductances, and the version of the definitions of the descriptors.

The conductances were computed with engine variant `holm`, baseline fingerprint
`014d59a0d7b143ba`. Version pinning is by `requirements.txt` together with that
fingerprint: the point is not that the code installs, but that it produces the
same numbers. `baseline.py --check` is what decides that.

The implementation carries three features the paper does not use —
intra-particle spreading resistance, surface roughness, and the specific
surface area ratio — and all three are disabled here. The reasons are given in
Section S.3 of the paper's ancillary material.

## Provenance

`rb_manifest.json` records, for every file in this record, the path it was
extracted from in the authors' working repository and its sha256.

**What you can check yourself, and what you cannot.** For `data/*.csv` and for
`thermal_network/L0_engine/baseline.json`, the recorded sha256 is the sha256 of
the file as distributed, so your own tool will reproduce it. For the `.py`
files it is not: the recorded value is the sha256 of the working-repository
source before the docstrings, strings and omitted definitions of the previous
section were applied, and that source is not distributed. Those entries record
where the file came from; they are not a checksum of what you have. The
baseline fingerprint is the check that matters for the numbers, and
`baseline.py --check` is how you run it.

Figure PDFs cannot be compared by sha256: matplotlib embeds a creation time, so
two runs of the same unchanged script differ byte for byte. Compare the
rendered pages instead.

## License

MIT (see `LICENSE`). The paper itself is CC BY 4.0.
