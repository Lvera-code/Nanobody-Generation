#!/usr/bin/env python3
import pandas as pd, numpy as np, os, datetime as dt
pd.set_option('display.width',250)
SP=os.path.dirname(os.path.abspath(__file__))
OUT=open(os.path.join(SP,'RESULTADOS.txt'),'a')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); OUT.write(s+'\n')
d=pd.read_pickle(f'{SP}/d_annot.pkl')
P('\n'+'='*100); P('8. COSTE REAL (tiempos extraidos de campaign_driver.log)'); P('='*100)
steps=[('1 design',290588.2),('2 inverse_folding',2509.3),('3 folding',417475.4),
       ('4 analysis',4768.4),('5 filtering',32.7)]
tot=sum(s for _,s in steps)
t0=dt.datetime(2026,9,15,12,38,29); t1=dt.datetime(2026,9,23,19,21,31)
wall=(t1-t0).total_seconds()
P(f'Inicio  : {t0}   Fin: {t1}   rc=0')
P(f'Wall clock total: {wall:.0f} s = {wall/3600:.2f} h = {wall/86400:.2f} dias')
P(f'Suma de los 5 pasos: {tot:.1f} s = {tot/3600:.2f} h  (cubre {100*tot/wall:.2f}% del wall clock)')
P(f'Overhead no contabilizado (configuracion/carga de modelos): {wall-tot:.0f} s = {(wall-tot)/60:.1f} min\n')
P(f'{"paso":22s} {"segundos":>12s} {"horas":>9s} {"% total":>9s} {"s/diseno":>10s} {"s/diseno/GPU":>14s}')
for n,s in steps:
    P(f'{n:22s} {s:12.1f} {s/3600:9.2f} {100*s/tot:8.2f}% {s/10000:10.2f} {s*4/10000:14.2f}')
P(f'{"TOTAL":22s} {tot:12.1f} {tot/3600:9.2f} {100:8.2f}% {tot/10000:10.2f} {tot*4/10000:14.2f}')
P(f'\nRendimiento agregado (4x Tesla T4):')
P(f'  {10000/(wall/3600):.2f} disenos/hora (wall clock, las 4 GPU juntas)')
P(f'  {10000/(tot/3600):.2f} disenos/hora (suma de pasos)')
P(f'  {10000/(wall/3600)/4:.2f} disenos/hora/GPU')
P(f'  {wall/10000:.1f} s/diseno (wall) ; {wall*4/10000:.1f} s-GPU/diseno')
P(f'  GPU-horas totales: {4*wall/3600:.1f}')
P('\nNOTA: el calculo s-GPU/diseno asume ocupacion exclusiva de las 4 T4. Carver NO es exclusiva')
P('(otros investigadores la usan en paralelo), asi que el coste en tiempo de pared puede incluir')
P('contencion no atribuible a esta campana. No he podido medirlo a posteriori.\n')
P('Coste por diseno "exitoso" segun cada criterio (GPU-horas por hit):')
gpuh=4*wall/3600
for k,n in [('A ipTM>=0.8 estricto',0),('B ipTM>=0.6 relajado',2),
            ('C pose correcta (bb_rmsd<2.5 & bb_target_aligned<2.5 & bindsite8>=.5)',169),
            ('D designabilidad bb_rmsd<2.5',448),('E pass_filters BoltzGen',407),
            ('top-100 entregado',100),('top-10 entregado',10)]:
    if n==0: P(f'  {k:70s} n=0   -> INFINITO (ningun hit)')
    else: P(f'  {k:70s} n={n:4d} -> {gpuh/n:8.1f} GPU-h/hit ; {wall/3600/n:7.2f} h-pared/hit')

P('\n'+'='*100); P('9. COMPARACION CON CORRIDAS ANTERIORES EN CARVER'); P('='*100)
import glob
P('Corridas encontradas (todas con all_designs_metrics.csv):')
info=[('rbd_campaign/run1 (ESTA)',f'{SP}/csv/all_designs_metrics.csv'),
      ('pfcsp_bench/run1 (piloto PfCSP 7RXP, 14-sep)',f'{SP}/csv/pfcsp_bench_all_designs_metrics.csv'),
      ('panel_tests/hard_targets/run1 (6M1U, 14-sep)',f'{SP}/csv/prev_panel_tests_hard_targets_run1.csv'),
      ('7RI1_cropped/run_A (HIV gp120/VHH-J3, 8-sep)',f'{SP}/csv/prev_7RI1_cropped_run_A.csv'),
      ('boltzgen_test/run (test instalacion, 4-sep)',f'{SP}/csv/prev_boltzgen_test_run.csv')]
rows=[]
for nm,p in info:
    if not os.path.exists(p): P(f'  {nm}: NO DESCARGADO'); continue
    x=pd.read_csv(p)
    r={'corrida':nm,'n':len(x)}
    for c in ['design_to_target_iptm','min_interaction_pae','bb_rmsd','design_ptm','bindsite_under_8rmsd']:
        r[c+'_med']=round(x[c].median(),4) if c in x.columns else np.nan
        r[c+'_max']=round(x[c].max(),4) if c in x.columns else np.nan
    for c in ['bb_rmsd<2.5','bb_target_aligned<2.5']:
        r[c]= f"{int(x[c].sum())}/{len(x)} ({100*x[c].mean():.1f}%)" if c in x.columns else 'NA'
    rows.append(r)
P(pd.DataFrame(rows).to_string(index=False))
P('\nCOMPARABILIDAD: solo pfcsp_bench y panel_tests/6M1U usan el mismo protocolo nanobody-anything;')
P('pero son pilotos de n muy pequeno (ver columna n), asi que no soportan comparacion estadistica.')
P('7RI1_cropped y boltzgen_test son de otra diana/otro momento del setup.')
P('El dato de referencia "0.41% de pose correcta" de la campana anterior NO esta en estos CSV:')
P('procede del panel de validacion previo, no de estos directorios. No lo he podido reproducir aqui.')
OUT.flush()
