# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/A_model_layer/M7_bond_definition/m7d_Rin_Rout.py
#   sha256(src) : 7eb3ab095676dd51
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 12 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V17. Isothermal-node ratio R_in / R_out: the internal (through-solid) resistance of
# a body against the parallel combination of all edges attached to it, for a fibre
# modelled as a chain of six spheres (486 conditions) and for a single sphere (Sec. 7.2).
# R_in / R_out << 1 means the isothermal-node treatment is justified.
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

import itertools, numpy as np, pandas as pd, networkx as nx
import ensemble as EN, materials_real as MR, thermal_sheaf_filtration as tsf

print("=== checking, on the actual branches, the condition under which the isothermal-node approximation holds ===")
print("verdict: internal resistance of a body vs the **parallel combination of all external branches** attached to that body")
print("(comparing against a single branch is incorrect; with many branches in parallel the external resistance is correspondingly smaller)")
print()
D0=20.0
rows=[]
for f,m,t in itertools.product(MR.FILLERS_GRID,MR.MATRICES_GRID,MR.TREATMENTS_GRID):
    ids,typ,pos,dia,mol,info=EN.make('fiber_z',D0,0)
    idx={int(v):k for k,v in enumerate(ids)}
    tb=MR.pair(f,t,'AlN',t,m)
    cfg=tsf.AnalysisConfig(gap_max_um=0.5*D0,eta_cut=1.0,slab_axis=2,slab_frac=0.06,viz=False)
    G=tsf.build_thermal_network(ids,typ,pos,dia,cfg,table=tb,mols=mol)
    Gb=nx.Graph(); Gb.add_nodes_from(G.nodes())
    for u,v,dd in G.edges(data=True):
        if dd['branch']=='bond': Gb.add_edge(u,v)
    for comp in nx.connected_components(Gb):
        if len(comp)<2: continue
        Rin_model=sum(1.0/G[u][v]['g'] for u,v in Gb.subgraph(comp).edges())
        gout=sum(dd['g'] for u,v,dd in G.edges(comp,data=True)
                 if dd['branch']!='bond')
        n_out=sum(1 for u,v,dd in G.edges(comp,data=True) if dd['branch']!='bond')
        rows.append(dict(filler=f,matrix=m,treat=t,n_bead=len(comp),
                         R_in=Rin_model,R_out=1.0/gout,n_out=n_out,
                         ratio=Rin_model/(1.0/gout)))
d=pd.DataFrame(rows); d.to_csv('results_m7d.csv',index=False)
print("fiber_z, per fibre (%d conditions x %d fibres)"%(d.filler.nunique()*d.matrix.nunique()*d.treat.nunique(),len(d)//54))
print("  bead count %d / attached non-bond branches, median %d"%(d.n_bead.median(),d.n_out.median()))
print("  internal resistance R_in (model)      median %9.1f K/W"%d.R_in.median())
print("  external resistance R_out (parallel)  median %9.1f K/W"%d.R_out.median())
print("  **R_in / R_out                         median %.2f** (<<1 means the isothermal treatment is valid)"%d.ratio.median())
print()
print("  R_in/R_out by treatment:")
print(d.groupby('treat').ratio.agg(['median','min','max']).to_string())
print()
print("  R_in/R_out by kappa (by filler):")
print(d.groupby('filler').ratio.median().sort_values().to_string())
print()
print("=== comparison: a single sphere ===")
ids,typ,pos,dia,mol,info=EN.make('phi_060',D0,0)
tb=MR.pair('AlN','untreated','AlN','untreated','epoxy')
cfg=tsf.AnalysisConfig(gap_max_um=0.5*D0,eta_cut=1.0,slab_axis=2,slab_frac=0.06,viz=False)
G=tsf.build_thermal_network(ids,typ,pos,dia,cfg,table=tb,mols=mol)
rr=[]
for n in list(G.nodes())[:300]:
    gout=sum(dd['g'] for _,_,dd in G.edges(n,data=True))
    if gout<=0: continue
    Rin=1.0/(MR.KAPPA['AlN']*dia[0]*1e-6)
    rr.append(Rin*gout)
print("  sphere (AlN/epoxy, d=20um): R_in/R_out median %.3f (n=%d)"%(np.median(rr),len(rr)))
