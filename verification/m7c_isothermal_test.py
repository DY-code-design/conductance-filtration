# ---------------------------------------------------------------------------
# Public software record accompanying the paper.
# Auto-extracted from the authors' working repository; do not edit here.
#   source      : 03_poc/A_model_layer/M7_bond_definition/m7c_isothermal_test.py
#   sha256(src) : b8ea28ef618f0e3d
# Comments are a short explanation written for this record; the working
# repository's development notes and docstrings are not included. The code is
# unmodified: its abstract syntax tree, with docstrings removed, is identical
# to the source except that 15 user-facing strings (messages and labels) were
# translated from Japanese; with string constants masked the two trees are
# identical.
# ---------------------------------------------------------------------------
# V17. Isothermal-node ratio: the internal resistance of a particle (sphere or
# fibre) compared with the resistance of the edges attached to it, and the effect
# of three treatments of the bond branch on the network conductance.
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

import itertools, numpy as np, pandas as pd
import ensemble as EN, materials_real as MR, thermal_sheaf_filtration as tsf

def one_d_R(ri,rj,D,ki,kj):
    xi=(D**2+ri**2-rj**2)/(2*D); xj=D-xi
    return (np.arctanh(min(xi/ri,1-1e-15))/(ki*np.pi*ri*1e-6)
            +np.arctanh(min(xj/rj,1-1e-15))/(kj*np.pi*rj*1e-6))

print("=== (1) is the isothermal-node approximation admissible for fibres too? ===")
print("criterion: internal resistance of a body << resistance of the branches attached to that body")
R=10.0; k=170.0; L=18.0; nb=6
R_fib=(nb-1)*L*1e-6/(k*np.pi*(R*1e-6)**2)
R_sph=1.0/(k*2*R*1e-6)
tab=MR.pair('AlN','untreated','AlN','untreated','epoxy')
ri,rj=MR.radii_from_rt(5.0)
g_edge,_=tsf.pair_conductance(0.1,ri,rj,1,2,tab,tsf.AnalysisConfig(gap_max_um=10.0,eta_cut=1.0))
R_edge=1.0/g_edge
print("  internal resistance of a fibre (6 beads, length %.0fum)  %8.1f K/W"%((nb-1)*L,R_fib))
print("  internal scale 1/(kappa d) of a sphere (d=20um)         %8.1f K/W"%R_sph)
print("  resistance of a typical near-contact branch (h=0.1um)   %8.1f K/W"%R_edge)
print("  -> edge / internal = fibre %.0f x / sphere %.0f x (comparison with a single edge; the comparison with the parallel combination of all attached edges is m7d_Rin_Rout.py)"%(R_edge/R_fib,R_edge/R_sph))
print()
print("=== (2) sphere x fibre (different mol): what happens in the current model? ===")
print("  different mol, so no bond branch: ordinary contact / near-contact branches. **the fibre bulk does not appear** = as designed.")
print("  however, since the fibre is a bead chain, the contact is evaluated as sphere x sphere. A true sphere x cylinder is an elliptical contact, so")
Rs=10.0; Rc=10.0
rt_ss=Rs*Rc/(Rs+Rc); rt_a=Rs*Rc/(Rs+Rc); rt_b=Rs
print("    sphere x sphere   : intersection-circle area ~ r~      = %.2f um"%rt_ss)
print("    sphere x cylinder : ellipse area   ~ sqrt(r~_a*r~_b) = sqrt(%.1fx%.1f) = %.2f um"%(rt_a,rt_b,np.sqrt(rt_a*rt_b)))
print("  -> the bead chain underestimates the sphere x fibre contact area by about %.0f%% (for equal radii)"%(100*(np.sqrt(rt_a*rt_b)/rt_ss-1)))
print()
print("=== (3) comparing the network under 3 bond-branch treatments (486 conditions) ===")
FAMS=['fiber_z','fiber_xy','jittered_ref']; D0=20.0
mats=list(itertools.product(MR.FILLERS_GRID,MR.MATRICES_GRID,MR.TREATMENTS_GRID))
rows=[]
for fam in FAMS:
    for seed in (0,1,2):
        ids,typ,pos,dia,mol,info=EN.make(fam,D0,seed)
        idx={int(v):kk for kk,v in enumerate(ids)}
        for f,m,t in mats:
            tb=MR.pair(f,t,'AlN',t,m)
            cfg=tsf.AnalysisConfig(gap_max_um=0.5*D0,eta_cut=1.0,slab_axis=2,slab_frac=0.06,viz=False)
            G=tsf.build_thermal_network(ids,typ,pos.copy(),dia,cfg,table=tb,mols=mol)
            src,snk=tsf.identify_slabs(ids,pos,cfg)
            if not src or not snk: continue
            Gt,_=tsf.attach_virtual_terminals(G,src,snk,cfg); sh=tsf.CellularSheaf(Gt)
            G0=tsf.effective_conductance(sh,0.0)
            if not (G0>0): continue
            be=[(u,v) for u,v,dd in Gt.edges(data=True) if u>=0 and v>=0 and dd.get('branch')=='bond']
            if not be:
                rows.append(dict(family=fam,seed=seed,filler=f,matrix=m,treat=t,G_cur=G0,G_1d=G0,G_iso=G0)); continue
            gmax=max(dd['g'] for _,_,dd in Gt.edges(data=True))
            sv=[(u,v,Gt[u][v]['g']) for u,v in be]
            for u,v in be:
                riu,rjv=dia[idx[u]]/2,dia[idx[v]]/2; Dd=Gt[u][v]['distance_um']
                ku=MR.KAPPA[f] if typ[idx[u]]==1 else MR.KAPPA['AlN']
                kv=MR.KAPPA[f] if typ[idx[v]]==1 else MR.KAPPA['AlN']
                Gt[u][v]['g']=1.0/one_d_R(riu,rjv,Dd,ku,kv)
            G1=tsf.effective_conductance(sh,0.0)
            for u,v in be: Gt[u][v]['g']=gmax*1e5
            G2=tsf.effective_conductance(sh,0.0)
            for u,v,g in sv: Gt[u][v]['g']=g
            rows.append(dict(family=fam,seed=seed,filler=f,matrix=m,treat=t,G_cur=G0,G_1d=G1,G_iso=G2))
d=pd.DataFrame(rows); d.to_csv('results_m7c_iso.csv',index=False)
d['r1d']=d.G_1d/d.G_cur; d['riso']=d.G_iso/d.G_cur
print('  family         1D/current   isothermal(short)/current   isothermal/1D')
for fam,g in d.groupby('family'):
    print('  %-13s %10.3f %16.3f %11.3f'%(fam,g.r1d.median(),g.riso.median(),(g.G_iso/g.G_1d).median()))
p=d.groupby(['family','seed','filler','matrix','treat'])[['G_cur','G_1d','G_iso']].first().reset_index()
z=p[p.family=='fiber_z'].set_index(['seed','filler','matrix','treat']); x=p[p.family=='fiber_xy'].set_index(['seed','filler','matrix','treat'])
print()
print('  anisotropy ratio fiber_z/fiber_xy : current %.3f / 1D %.3f / **isothermal %.3f**'%(
    (z.G_cur/x.G_cur).median(),(z.G_1d/x.G_1d).median(),(z.G_iso/x.G_iso).median()))
