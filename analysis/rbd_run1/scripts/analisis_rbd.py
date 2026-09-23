#!/usr/bin/env python3
"""Analisis exhaustivo campana BoltzGen RBD (eje 1 TFG). Solo lectura."""
import pandas as pd, numpy as np, json, os, sys
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 60)
SP = os.path.dirname(os.path.abspath(__file__))
OUT = open(os.path.join(SP, 'RESULTADOS.txt'), 'w')
def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s); OUT.write(s + '\n')

d  = pd.read_csv(f'{SP}/csv/all_designs_metrics.csv')
ag = pd.read_csv(f'{SP}/csv/aggregate_metrics_analyze.csv')
t100 = pd.read_csv(f'{SP}/csv/final_designs_metrics_100.csv')
P('='*100); P('0. INVENTARIO DE FILAS'); P('='*100)
P('all_designs_metrics.csv :', d.shape)
P('aggregate_metrics       :', ag.shape)
P('final_designs_100       :', t100.shape)

# --- duplicados: 9990 -> 9988
P('\n-- Los 2 duplicados eliminados (9990 -> 9988) --')
if 'designed_chain_sequence' in ag.columns:
    dup = ag[ag.duplicated('designed_chain_sequence', keep=False)].sort_values('designed_chain_sequence')
    P('filas en aggregate con designed_chain_sequence repetida:', len(dup))
    for s, g in dup.groupby('designed_chain_sequence'):
        P('  seq len', len(s), 'ids:', g['id'].tolist())
ids_ag = set(ag['id']); ids_d = set(d['id'])
P('ids en aggregate y no en all_designs:', sorted(ids_ag - ids_d))

# --- los 10 folds fallidos
allidx = set(range(10000))
present = set(int(i.split('_')[-1]) for i in ag['id'])
P('\n-- Disenos generados sin metricas (fold fallido) --')
P('faltantes:', sorted(allidx - present))
P('n faltantes:', len(allidx - present))

# ===================== SCAFFOLD =====================
P('\n' + '='*100); P('1. CLASIFICACION POR SCAFFOLD (framework FR2/FR3)'); P('='*100)
FR = {  # motivos FR2 discriminantes extraidos de los .cif de referencia
 '7eow_caplacizumab' : 'WFRQAPGKGRELVA',
 '7xl0_vobarilizumab': 'WYRQAPGKGRELVA',
 '8coh_gefurulimab'  : 'WFRQAPGQEREFVA',
 '8z8v_ozoralizumab' : 'WVRQAPGKGLEWVS',
}
FR3 = {
 '7eow_caplacizumab' : 'PDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYC',
 '7xl0_vobarilizumab': 'ADSVKGRFTISRDNAKNTLYLQMNSLRPEDTAVYYC',
 '8coh_gefurulimab'  : 'ADSVRGRFTNSRDNSKNTLYLQMNSLRAEDTAVYYC',
 '8z8v_ozoralizumab' : 'ADSVKGRFTISRDNAKTTLYLQMNSLRPEDTAVYYC',
}
def classify(seq, motifs):
    hits = [k for k, m in motifs.items() if m in seq]
    return hits[0] if len(hits) == 1 else ('AMBIGUO' if hits else 'NINGUNO')
for df, nm in ((d, 'all_designs'), (t100, 'top100')):
    df['scaffold']  = df['designed_chain_sequence'].astype(str).apply(lambda s: classify(s, FR))
    df['scaffold3'] = df['designed_chain_sequence'].astype(str).apply(lambda s: classify(s, FR3))
P('Clasificacion por FR2 (all_designs):'); P(d['scaffold'].value_counts().to_string())
P('\nConcordancia FR2 vs FR3:', (d['scaffold'] == d['scaffold3']).sum(), '/', len(d))
P('Discordantes/NINGUNO por FR3:'); P(d['scaffold3'].value_counts().to_string())
d.to_pickle(f'{SP}/d_annot.pkl')

# ===================== DISTRIBUCIONES =====================
P('\n' + '='*100); P('2. DISTRIBUCIONES GLOBALES'); P('='*100)
key = ['design_to_target_iptm','iptm','design_iptm','design_iiptm','design_residue_iptm','protein_iptm',
       'ptm','design_ptm','target_ptm','min_interaction_pae','interaction_pae','min_design_to_target_pae',
       'complex_plddt','complex_iplddt','complex_pde','complex_ipde',
       'design_ipsae_min','design_to_target_ipsae','target_to_design_ipsae',
       'bb_rmsd','bb_rmsd_design','bb_rmsd_target','bb_rmsd_design_target','bb_target_aligned_rmsd_design',
       'filter_rmsd','filter_rmsd_design',
       'delta_sasa_refolded','plip_hbonds_refolded','plip_saltbridge_refolded',
       'design_hydrophobicity','design_chain_hydrophobicity','liability_score','liability_num_violations',
       'num_design','num_prot_tokens','loop','helix','sheet']
qs = [0,.01,.05,.10,.25,.50,.75,.90,.95,.99,1.0]
desc = d[key].describe(percentiles=qs[1:-1]).T
P(desc.to_string())

P('\n-- Columnas booleanas (recuento True / %) --')
bools = [c for c in d.columns if d[c].dtype == bool]
for c in bools:
    P(f'  {c:38s} {d[c].sum():6d}  {100*d[c].mean():6.2f}%')

# ===================== EPITOPO =====================
P('\n' + '='*100); P('3. FIDELIDAD DE EPITOPO (bindsite_under_*)'); P('='*100)
P('SEMANTICA (boltzgen/task/analyze/analyze.py:832-843): para cada token del sitio de union')
P('se calcula la distancia center_coords al token DISEÑADO mas cercano; la columna es la')
P('FRACCION de residuos de epitopo con esa distancia < umbral (A). NO es un RMSD.')
P('Se mide sobre la estructura GENERADA (intermediate_designs_inverse_folded), no el refold.')
bs = [f'bindsite_under_{t}rmsd' for t in [3,4,5,6,7,8,9]]
P('\nEpitopo pedido: 16 residuos (chain Exp 148..163 = Spike P479..S494, motivo RBM/ACE2).')
P('Por eso la fraccion se mueve en multiplos de 1/16 = 0.0625.\n')
for c in bs:
    v = d[c]
    P(f'{c}: media {v.mean():.4f}  mediana {v.median():.4f}  max {v.max():.4f}  '
      f'n(>0) {int((v>0).sum())} ({100*(v>0).mean():.2f}%)')
P('\nDistribucion de valores discretos (n residuos de 16 en contacto):')
for c in bs:
    vc = (d[c]*16).round().astype(int).value_counts().sort_index()
    P(f'  {c}: ' + ', '.join(f'{k}res:{v}' for k, v in vc.items()))
P('\nCriterio canonico de BoltzGen (filter.py:255-260): filter_bindingsite usa')
P('bindsite_under_8rmsd > 0.0001 -> "al menos 1 residuo del epitopo a <8 A de un residuo disenado".')
P('En esta campana filter_bindingsite estaba en FALSE, asi que NO se aplico.')
n8 = (d['bindsite_under_8rmsd'] > 0).sum()
P(f'\nSi se hubiera aplicado: {n8} / {len(d)} = {100*n8/len(d):.2f}% pasarian.')
for thr in [1,2,4,8,16]:
    n = (d['bindsite_under_8rmsd']*16 >= thr-0.01).sum()
    P(f'  >= {thr:2d}/16 residuos de epitopo a <8 A : {n:5d} ({100*n/len(d):6.2f}%)')
P('\nnative_rmsd / native_rmsd_bb (no calculados, native=false):')
P('  native_rmsd unicos:', d['native_rmsd'].unique(), ' native_rmsd_bb unicos:', d['native_rmsd_bb'].unique())
P('  native_rmsd_refolded unicos:', d['native_rmsd_refolded'].unique()[:5],
  ' native_rmsd_bb_refolded unicos:', d['native_rmsd_bb_refolded'].unique()[:5])
OUT.flush()
