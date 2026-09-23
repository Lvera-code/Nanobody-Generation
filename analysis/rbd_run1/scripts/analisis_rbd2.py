#!/usr/bin/env python3
import pandas as pd, numpy as np, os
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 80)
SP = os.path.dirname(os.path.abspath(__file__))
OUT = open(os.path.join(SP, 'RESULTADOS.txt'), 'a')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); OUT.write(s+'\n')
d = pd.read_pickle(f'{SP}/d_annot.pkl')
t100 = pd.read_csv(f'{SP}/csv/final_designs_metrics_100.csv')

P('\n'+'='*100); P('4. DESGLOSE POR SCAFFOLD (las 4 plantillas de nanocuerpo)'); P('='*100)
P('AVISO: per_target_metrics_analyze.csv tiene UNA sola fila (target_id="design"):')
P('BoltzGen trata las 4 plantillas como una unica "diana". El desglose siguiente lo he')
P('reconstruido yo clasificando cada diseno por su framework FR2/FR3 (100% concordante con FR3).\n')
cols = ['design_to_target_iptm','min_interaction_pae','design_ptm','complex_plddt','complex_iplddt',
        'bb_rmsd','bb_rmsd_design','bb_target_aligned_rmsd_design','bindsite_under_8rmsd',
        'bindsite_under_6rmsd','delta_sasa_refolded','plip_hbonds_refolded','num_design','liability_score']
g = d.groupby('scaffold')
P('n por scaffold:'); P(g.size().to_string())
P('\nMEDIA por scaffold:'); P(g[cols].mean().round(4).to_string())
P('\nMEDIANA por scaffold:'); P(g[cols].median().round(4).to_string())
P('\nTasas de flags booleanos por scaffold (%):')
bl = ['bb_rmsd<2.5','bb_target_aligned<2.5','bb_designability_rmsd_2','bb_designability_rmsd_4',
      'min_interaction_pae<5','design_ptm>80','pass_filters']
P((g[bl].mean()*100).round(2).to_string())
P('\nRecuentos absolutos:'); P(g[bl].sum().to_string())
P('\nMaximos de design_to_target_iptm por scaffold:'); P(g['design_to_target_iptm'].max().round(4).to_string())
P('\nTop-100 final por scaffold:'); 
FR={'7eow_caplacizumab':'WFRQAPGKGRELVA','7xl0_vobarilizumab':'WYRQAPGKGRELVA',
    '8coh_gefurulimab':'WFRQAPGQEREFVA','8z8v_ozoralizumab':'WVRQAPGKGLEWVS'}
def cl(s):
    h=[k for k,m in FR.items() if m in s]; return h[0] if len(h)==1 else 'NINGUNO'
t100['scaffold']=t100['designed_chain_sequence'].astype(str).apply(cl)
P(t100['scaffold'].value_counts().to_string())
P('\nEnriquecimiento top100 vs base (%base -> %top100):')
for s in sorted(d['scaffold'].unique()):
    P(f'  {s:22s} {100*(d["scaffold"]==s).mean():5.2f}% -> {100*(t100["scaffold"]==s).mean():5.2f}%')
# test estadistico simple
from scipy import stats as st
P('\nKruskal-Wallis entre los 4 scaffolds:')
for c in ['design_to_target_iptm','min_interaction_pae','bb_rmsd','bb_target_aligned_rmsd_design','bindsite_under_8rmsd']:
    grp=[gg[c].values for _,gg in g]
    h,p = st.kruskal(*grp)
    P(f'  {c:32s} H={h:9.2f}  p={p:.3e}')

P('\n'+'='*100); P('5. HIT RATE CON CRITERIOS DE DISTINTA EXIGENCIA'); P('='*100)
P('Justificacion de umbrales: la literatura de diseno de binders (RFdiffusion/BindCraft/AF2)')
P('usa tipicamente ipTM>=0.8 y PAE_interaction<10 como filtro in-silico de "binder creible";')
P('BindCraft usa ipTM>=0.8. BoltzGen aporta ademas designabilidad (bb_rmsd<2.5 = el refold')
P('reproduce la pose disenada) y adherencia al sitio (bindsite_under_8rmsd).\n')
crit = {}
crit['A_literatura_estricto (ipTM>=0.8 & min_iPAE<10 & bb_rmsd<2.5 & bindsite8>=0.5)'] = (
    (d['design_to_target_iptm']>=0.8) & (d['min_interaction_pae']<10) & d['bb_rmsd<2.5'] & (d['bindsite_under_8rmsd']>=0.5))
crit['B_literatura_relajado (ipTM>=0.6 & min_iPAE<15 & bb_rmsd<2.5 & bindsite8>=0.5)'] = (
    (d['design_to_target_iptm']>=0.6) & (d['min_interaction_pae']<15) & d['bb_rmsd<2.5'] & (d['bindsite_under_8rmsd']>=0.5))
crit['C_pose_correcta (bb_rmsd<2.5 & bb_target_aligned<2.5 & bindsite8>=0.5)'] = (
    d['bb_rmsd<2.5'] & d['bb_target_aligned<2.5'] & (d['bindsite_under_8rmsd']>=0.5))
crit['D_designabilidad_sola (bb_rmsd<2.5)'] = d['bb_rmsd<2.5']
crit['E_pipeline_BoltzGen (pass_filters tal cual corrio)'] = d['pass_filters']
crit['F_pipeline+bindingsite (pass_filters & bindsite8>0)'] = d['pass_filters'] & (d['bindsite_under_8rmsd']>0)
crit['G_top_decil_confianza (ipTM>=0.3 & min_iPAE<15)'] = (d['design_to_target_iptm']>=0.3)&(d['min_interaction_pae']<15)
for k,v in crit.items():
    n=int(v.sum()); P(f'  {k:72s} n={n:5d}  {100*n/len(d):7.3f}%')
P('\nComponentes por separado (para ver que factor limita):')
comp = {'ipTM>=0.8':(d['design_to_target_iptm']>=0.8),'ipTM>=0.7':(d['design_to_target_iptm']>=0.7),
        'ipTM>=0.6':(d['design_to_target_iptm']>=0.6),'ipTM>=0.5':(d['design_to_target_iptm']>=0.5),
        'ipTM>=0.4':(d['design_to_target_iptm']>=0.4),'ipTM>=0.3':(d['design_to_target_iptm']>=0.3),
        'min_iPAE<5':(d['min_interaction_pae']<5),'min_iPAE<10':(d['min_interaction_pae']<10),
        'min_iPAE<15':(d['min_interaction_pae']<15),'bb_rmsd<2.5':d['bb_rmsd<2.5'],
        'bb_target_aligned<2.5':d['bb_target_aligned<2.5'],'bindsite8>=0.5':(d['bindsite_under_8rmsd']>=0.5),
        'bindsite8==1.0':(d['bindsite_under_8rmsd']>=0.999)}
for k,v in comp.items():
    n=int(v.sum()); P(f'  {k:26s} n={n:5d}  {100*n/len(d):7.3f}%')
P('\nCruce designabilidad x pose (tabla de contingencia bb_rmsd<2.5 vs bb_target_aligned<2.5):')
P(pd.crosstab(d['bb_rmsd<2.5'], d['bb_target_aligned<2.5']).to_string())
P('\nDe los 448 con bb_rmsd<2.5, distribucion de ipTM:')
sub=d[d['bb_rmsd<2.5']]
P(sub['design_to_target_iptm'].describe(percentiles=[.25,.5,.75,.9,.95,.99]).round(4).to_string())
P('\nDe los 171 con bb_target_aligned<2.5:')
sub2=d[d['bb_target_aligned<2.5']]
P(sub2[['design_to_target_iptm','min_interaction_pae','bb_rmsd','bindsite_under_8rmsd']].describe().round(4).to_string())
OUT.flush()
