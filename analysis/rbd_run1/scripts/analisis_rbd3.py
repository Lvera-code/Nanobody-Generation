#!/usr/bin/env python3
import pandas as pd, numpy as np, os, itertools, collections
pd.set_option('display.width',250); pd.set_option('display.max_columns',80)
SP=os.path.dirname(os.path.abspath(__file__))
OUT=open(os.path.join(SP,'RESULTADOS.txt'),'a')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); OUT.write(s+'\n')
d=pd.read_pickle(f'{SP}/d_annot.pkl')
t100=pd.read_csv(f'{SP}/csv/final_designs_metrics_100.csv')
FR={'7eow_caplacizumab':'WFRQAPGKGRELVA','7xl0_vobarilizumab':'WYRQAPGKGRELVA',
    '8coh_gefurulimab':'WFRQAPGQEREFVA','8z8v_ozoralizumab':'WVRQAPGKGLEWVS'}
def cl(s):
    h=[k for k,m in FR.items() if m in s]; return h[0] if len(h)==1 else 'NINGUNO'
t100['scaffold']=t100['designed_chain_sequence'].astype(str).apply(cl)

P('\n'+'='*100); P('6. TOP 100 Y TOP 10: CARACTERIZACION Y REDUNDANCIA'); P('='*100)
cols=['design_to_target_iptm','min_interaction_pae','design_ptm','complex_plddt','bb_rmsd',
      'bb_rmsd_design','bb_target_aligned_rmsd_design','bindsite_under_8rmsd','delta_sasa_refolded',
      'plip_hbonds_refolded','plip_saltbridge_refolded','num_design','liability_score','max_rank']
P('Resumen TOP-100:'); P(t100[cols].describe(percentiles=[.25,.5,.75]).round(4).T.to_string())
P('\nFlags booleanos TOP-100:')
for c in ['bb_rmsd<2.5','bb_target_aligned<2.5','bb_designability_rmsd_2','min_interaction_pae<5','design_ptm>80']:
    if c in t100: P(f'  {c:26s} {int(t100[c].sum()):3d} / 100')
P('\nComparativa TOP100 vs poblacion (media):')
for c in cols:
    if c in d.columns: P(f'  {c:32s} pobl {d[c].mean():10.4f}   top100 {t100[c].mean():10.4f}')

P('\n-- REDUNDANCIA DE SECUENCIA --')
def redund(df,name):
    seqs=df['designed_sequence'].astype(str).tolist()
    P(f'\n[{name}] n={len(seqs)}  unicas={len(set(seqs))}')
    cseq=df['designed_chain_sequence'].astype(str).tolist()
    P(f'  cadenas completas unicas: {len(set(cseq))}')
    L=collections.Counter(len(s) for s in seqs)
    P(f'  longitudes designed_sequence: min {min(map(len,seqs))} max {max(map(len,seqs))}')
    # identidad por pares sobre secuencias de igual longitud (ungapped)
    pares=0; sim=[]
    for a,b in itertools.combinations(range(len(seqs)),2):
        s1,s2=seqs[a],seqs[b]
        if len(s1)!=len(s2): continue
        pares+=1; sim.append(sum(x==y for x,y in zip(s1,s2))/len(s1))
    if sim:
        sim=np.array(sim)
        P(f'  pares de igual longitud: {pares} de {len(seqs)*(len(seqs)-1)//2}')
        P(f'  identidad media {sim.mean():.3f}  mediana {np.median(sim):.3f}  max {sim.max():.3f}')
        P(f'  pares con identidad >=0.9: {(sim>=0.9).sum()}   >=0.8: {(sim>=0.8).sum()}   >=0.7: {(sim>=0.7).sum()}')
    # clustering greedy 70% identidad (global, con normalizacion por longitud del mas corto)
    def ident(s1,s2):
        n=min(len(s1),len(s2))
        best=0
        for off in range(abs(len(s1)-len(s2))+1):
            a,b=(s1,s2[off:off+n]) if len(s2)>=len(s1) else (s1[off:off+n],s2)
            best=max(best,sum(x==y for x,y in zip(a,b))/n)
        return best
    reps=[]
    for s in seqs:
        if not any(ident(s,r)>=0.7 for r in reps): reps.append(s)
    P(f'  clusters greedy al 70% identidad: {len(reps)}')
    reps8=[]
    for s in seqs:
        if not any(ident(s,r)>=0.8 for r in reps8): reps8.append(s)
    P(f'  clusters greedy al 80% identidad: {len(reps8)}')
    # motivos CDR3-like: ultimos 12 aa antes de WGQG
    return seqs
redund(t100,'TOP-100')
# top10 intermedio
P('\n-- TOP 10 intermedio (intermediate_ranked_10_designs) --')
t10=d.nsmallest(10,'final_rank') if 'final_rank' in d else None
P('Nota: los 10 "intermediate_ranked" son los 10 de mayor calidad SIN el criterio de diversidad;')
P('los muestro como el top-10 por final_rank de all_designs_metrics.csv:')
t10=d.sort_values('final_rank').head(10)
P(t10[['id','final_rank','scaffold','design_to_target_iptm','min_interaction_pae','design_ptm',
       'bb_rmsd','bindsite_under_8rmsd','complex_plddt','num_design']].round(4).to_string(index=False))
redund(t10,'TOP-10')
P('\nSecuencias designed_sequence del TOP-10:')
for _,r in t10.iterrows(): P(f"  rank{int(r['final_rank']):4d} {r['scaffold'][:18]:18s} {r['designed_sequence']}")
P('\nSecuencias designed_sequence del TOP-100 (primeras 25):')
for i,(_,r) in enumerate(t100.head(25).iterrows()):
    P(f"  {i+1:3d} {r['scaffold'][:18]:18s} {r['designed_sequence']}")

P('\n'+'='*100); P('7. CORRELACIONES Y DOMINANCIA DEL RANKING'); P('='*100)
num=['design_to_target_iptm','design_iptm','design_iiptm','iptm','ptm','design_ptm','target_ptm',
     'min_interaction_pae','interaction_pae','min_design_to_target_pae','complex_plddt','complex_iplddt',
     'complex_ipde','bb_rmsd','bb_rmsd_design','bb_rmsd_target','bb_target_aligned_rmsd_design',
     'bindsite_under_6rmsd','bindsite_under_8rmsd','delta_sasa_refolded','plip_hbonds_refolded',
     'plip_saltbridge_refolded','design_hydrophobicity','liability_score','num_design','loop','helix','sheet']
cm=d[num].corr(method='spearman')
P('Correlacion de Spearman con design_to_target_iptm (ordenada):')
P(cm['design_to_target_iptm'].drop('design_to_target_iptm').sort_values(key=abs,ascending=False).round(3).to_string())
P('\nCorrelacion de Spearman con bb_rmsd (designabilidad; menor es mejor):')
P(cm['bb_rmsd'].drop('bb_rmsd').sort_values(key=abs,ascending=False).round(3).to_string())
P('\nCorrelacion de Spearman con bindsite_under_8rmsd:')
P(cm['bindsite_under_8rmsd'].drop('bindsite_under_8rmsd').sort_values(key=abs,ascending=False).round(3).to_string())
P('\n-- Dominancia del ranking final --')
rk=['rank_design_to_target_iptm','rank_design_ptm','rank_neg_min_design_to_target_pae',
    'rank_plip_hbonds_refolded','rank_plip_saltbridge_refolded','rank_delta_sasa_refolded']
P('El ranking de BoltzGen es max_rank = maximo de los 6 sub-rankings (peor posicion manda).')
P('Spearman de final_rank con cada sub-ranking (sobre los 407 que pasan filtros):')
sub=d[d['pass_filters']]
for r in rk+['max_rank']:
    P(f'  {r:38s} rho={sub["final_rank"].corr(sub[r],method="spearman"):.4f}')
P('\nFrecuencia con la que cada sub-ranking ES el maximo (el que define max_rank), n=407:')
arg=sub[rk].idxmax(axis=1).value_counts()
P(arg.to_string())
P('\nSobre los 9988:')
P(d[rk].idxmax(axis=1).value_counts().to_string())
P('\nSpearman de final_rank con las metricas crudas (407 que pasan filtros):')
for c in ['design_to_target_iptm','design_ptm','min_design_to_target_pae','plip_hbonds_refolded',
          'plip_saltbridge_refolded','delta_sasa_refolded','bb_rmsd','bindsite_under_8rmsd']:
    P(f'  {c:32s} rho={sub["final_rank"].corr(sub[c],method="spearman"):.4f}')
cm.to_csv(f'{SP}/correlaciones_spearman.csv')
OUT.flush()
