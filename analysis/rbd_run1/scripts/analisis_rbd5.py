#!/usr/bin/env python3
import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
pd.set_option('display.width',250)
SP=os.path.dirname(os.path.abspath(__file__))
FIG=os.path.join(SP,'figuras'); os.makedirs(FIG,exist_ok=True)
OUT=open(os.path.join(SP,'RESULTADOS.txt'),'a')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); OUT.write(s+'\n')
d=pd.read_pickle(f'{SP}/d_annot.pkl')
t100=pd.read_csv(f'{SP}/csv/final_designs_metrics_100.csv')

P('\n'+'='*100); P('10. CONTROL: por que domina 8z8v/ozoralizumab?'); P('='*100)
P('Hipotesis: 8z8v es el scaffold con menos residuos rediseñados y CDR3 mas corto,')
P('lo que facilita que el refold reproduzca la pose (bb_rmsd<2.5).')
P('\nLongitud de region rediseñada (num_design) por scaffold:')
P(d.groupby('scaffold')['num_design'].describe().round(2).to_string())
P('\nTasa bb_rmsd<2.5 por scaffold DENTRO de bins de num_design (control):')
d['bin']=pd.cut(d['num_design'],[19,25,30,35,40,50])
tab=d.pivot_table(index='bin',columns='scaffold',values='bb_rmsd<2.5',aggfunc='mean',observed=True)*100
P(tab.round(2).to_string())
P('\nn por celda:')
P(d.pivot_table(index='bin',columns='scaffold',values='bb_rmsd<2.5',aggfunc='size',observed=True).to_string())
P('\nInterpretacion: la ventaja de 8z8v se mantiene dentro de cada bin de longitud,')
P('luego NO se explica solo por la longitud del diseño.')
P('\nSpearman num_design vs bb_rmsd por scaffold:')
for s,g in d.groupby('scaffold'):
    P(f'  {s:22s} rho={g["num_design"].corr(g["bb_rmsd"],method="spearman"):.3f}')
P('\ndesign_ptm por scaffold (confianza del pliegue del propio binder):')
P(d.groupby('scaffold')['design_ptm'].describe(percentiles=[.5,.9,.99]).round(4).to_string())

# ---------- FIGURAS ----------
plt.rcParams.update({'figure.dpi':130,'font.size':8,'axes.grid':True,'grid.alpha':.3})
SC=sorted(d['scaffold'].unique()); COL=dict(zip(SC,['#4C72B0','#DD8452','#55A868','#C44E52']))

fig,ax=plt.subplots(2,3,figsize=(13,7))
ax=ax.ravel()
ax[0].hist(d['design_to_target_iptm'],bins=80,color='#4C72B0'); ax[0].axvline(.8,c='r',ls='--',label='ipTM 0.8 (BindCraft)')
ax[0].axvline(.6,c='orange',ls='--',label='ipTM 0.6'); ax[0].set_title('ipTM diseño->diana'); ax[0].legend(fontsize=6); ax[0].set_yscale('log')
ax[1].hist(d['min_interaction_pae'],bins=80,color='#DD8452'); ax[1].axvline(10,c='r',ls='--',label='PAE 10')
ax[1].set_title('min_interaction_pae (A)'); ax[1].legend(fontsize=6)
ax[2].hist(d['bb_rmsd'],bins=80,color='#55A868'); ax[2].axvline(2.5,c='r',ls='--',label='2.5 A')
ax[2].set_title('bb_rmsd complejo (diseñado vs refold)'); ax[2].legend(fontsize=6)
ax[3].hist(d['design_ptm'],bins=80,color='#C44E52'); ax[3].set_title('design_ptm')
ax[4].hist(d['complex_plddt'],bins=80,color='#8172B3'); ax[4].set_title('complex_plddt')
ax[5].hist(d['bb_target_aligned_rmsd_design'],bins=80,color='#937860'); ax[5].axvline(2.5,c='r',ls='--')
ax[5].set_title('bb_target_aligned_rmsd_design (pose)'); ax[5].set_yscale('log')
fig.suptitle('Campaña RBD 7BEP / BoltzGen — distribuciones globales (n=9.988)',fontweight='bold')
fig.tight_layout(); fig.savefig(f'{FIG}/01_distribuciones_globales.png',bbox_inches='tight'); plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(12,4.2))
th=[3,4,5,6,7,8,9]
m=[d[f'bindsite_under_{t}rmsd'].mean() for t in th]
q1=[d[f'bindsite_under_{t}rmsd'].quantile(.25) for t in th]
q3=[d[f'bindsite_under_{t}rmsd'].quantile(.75) for t in th]
ax[0].plot(th,m,'o-',color='#4C72B0',label='media'); ax[0].fill_between(th,q1,q3,alpha=.25,color='#4C72B0',label='IQR')
ax[0].axvline(8,c='r',ls='--',label='umbral canonico BoltzGen (8 A)')
ax[0].set_xlabel('umbral de distancia centro-centro (A)'); ax[0].set_ylabel('fraccion del epitopo (16 res) en contacto')
ax[0].set_title('Adherencia al epitopo pedido (Spike 479-494)'); ax[0].legend(fontsize=7)
ax[1].hist(d['bindsite_under_8rmsd']*16,bins=np.arange(-.5,17.5,1),color='#55A868',edgecolor='w')
ax[1].set_xlabel('nº de residuos del epitopo (de 16) a <8 A de un residuo diseñado')
ax[1].set_ylabel('nº de diseños'); ax[1].set_title('bindsite_under_8rmsd — distribucion')
fig.tight_layout(); fig.savefig(f'{FIG}/02_epitopo.png',bbox_inches='tight'); plt.close(fig)

fig,ax=plt.subplots(1,4,figsize=(15,3.8))
for i,c in enumerate(['design_to_target_iptm','min_interaction_pae','bb_rmsd','design_ptm']):
    data=[d[d['scaffold']==s][c].values for s in SC]
    bp=ax[i].boxplot(data,labels=[s.split('_')[0] for s in SC],showfliers=False,patch_artist=True)
    for p,s in zip(bp['boxes'],SC): p.set_facecolor(COL[s])
    ax[i].set_title(c)
fig.suptitle('Por scaffold de nanocuerpo (n≈2.500 cada uno)',fontweight='bold')
fig.tight_layout(); fig.savefig(f'{FIG}/03_por_scaffold.png',bbox_inches='tight'); plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(11,4.5))
for s in SC:
    g=d[d['scaffold']==s]
    ax[0].scatter(g['bb_rmsd'],g['design_to_target_iptm'],s=2,alpha=.25,c=COL[s],label=s.split('_')[0])
ax[0].axvline(2.5,c='r',ls='--'); ax[0].axhline(.8,c='r',ls='--')
ax[0].set_xlabel('bb_rmsd (A)'); ax[0].set_ylabel('ipTM diseño->diana'); ax[0].legend(fontsize=7,markerscale=4)
ax[0].set_title('Designabilidad vs confianza de interfaz')
sc=ax[1].scatter(d['min_interaction_pae'],d['design_to_target_iptm'],s=2,alpha=.2,c=d['bindsite_under_8rmsd'],cmap='viridis')
ax[1].set_xlabel('min_interaction_pae (A)'); ax[1].set_ylabel('ipTM diseño->diana')
plt.colorbar(sc,ax=ax[1],label='bindsite_under_8rmsd')
ax[1].scatter(t100['min_interaction_pae'],t100['design_to_target_iptm'],s=14,facecolors='none',edgecolors='red',lw=.7,label='top-100')
ax[1].legend(fontsize=7); ax[1].set_title('PAE vs ipTM (top-100 marcado)')
fig.tight_layout(); fig.savefig(f'{FIG}/04_dispersion.png',bbox_inches='tight'); plt.close(fig)

# embudo
fig,ax=plt.subplots(figsize=(8,4.2))
lab=['generados','con metricas','sin duplicados','bb_rmsd<2.5','+rmsd_design','+GLY_frac\n(=pass_filters)','entregados\ntop-100','ipTM>=0.6']
val=[10000,9990,9988,448,420,407,100,2]
b=ax.bar(range(len(val)),val,color=['#999']*3+['#55A868']*3+['#4C72B0','#C44E52'])
ax.set_yscale('log'); ax.set_xticks(range(len(val))); ax.set_xticklabels(lab,fontsize=7)
for i,v in enumerate(val): ax.text(i,v*1.15,str(v),ha='center',fontsize=8,fontweight='bold')
ax.set_ylabel('nº de diseños (escala log)'); ax.set_title('Embudo de la campaña RBD',fontweight='bold')
fig.tight_layout(); fig.savefig(f'{FIG}/05_embudo.png',bbox_inches='tight'); plt.close(fig)

# correlaciones
cm=pd.read_csv(f'{SP}/correlaciones_spearman.csv',index_col=0)
fig,ax=plt.subplots(figsize=(10,8.5))
im=ax.imshow(cm.values,cmap='RdBu_r',vmin=-1,vmax=1)
ax.set_xticks(range(len(cm))); ax.set_xticklabels(cm.columns,rotation=90,fontsize=6)
ax.set_yticks(range(len(cm))); ax.set_yticklabels(cm.index,fontsize=6)
plt.colorbar(im,label='Spearman rho'); ax.set_title('Matriz de correlaciones (Spearman)',fontweight='bold')
ax.grid(False); fig.tight_layout(); fig.savefig(f'{FIG}/06_correlaciones.png',bbox_inches='tight'); plt.close(fig)
P('\nFiguras guardadas en',FIG)
for f in sorted(os.listdir(FIG)): P('  ',os.path.join(FIG,f))
OUT.flush()
