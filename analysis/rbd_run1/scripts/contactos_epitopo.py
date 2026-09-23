#!/usr/bin/env python3
"""Recalculo de contactos reales sobre los .cif del top-100 (refold y pre-refold).
Chain A = diana RBD (197 res, local 1..197 = Spike 332..528). Chain B = binder.
Epitopo pedido: local 148..163 = Spike 479..494."""
import numpy as np, os, glob, pandas as pd
SP=os.path.dirname(os.path.abspath(__file__))
OUT=open(os.path.join(SP,'RESULTADOS.txt'),'a')
def P(*a):
    s=' '.join(str(x) for x in a); print(s); OUT.write(s+'\n')
EPI=set(range(148,164)); CUT=5.0; OFF=331   # Spike = local + 331

def parse(path):
    A=[];B=[]
    for l in open(path):
        if not l.startswith(('ATOM','HETATM')): continue
        f=l.split()
        el=f[2]
        if el=='H': continue
        asym=f[6]; seqid=int(f[8]); xyz=(float(f[10]),float(f[11]),float(f[12]))
        (A if asym=='A' else B).append((seqid,xyz))
    return A,B

def contacts(path):
    A,B=parse(path)
    if not A or not B: return None
    ar=np.array([a[0] for a in A]); ax=np.array([a[1] for a in A])
    bx=np.array([b[1] for b in B])
    # min distance per target atom to any binder atom, chunked
    mind=np.full(len(ax),1e9)
    for i in range(0,len(bx),2000):
        dd=np.sqrt(((ax[:,None,:]-bx[None,i:i+2000,:])**2).sum(-1)).min(1)
        mind=np.minimum(mind,dd)
    res={}
    for r,m in zip(ar,mind):
        res[r]=min(res.get(r,1e9),m)
    cont=set(r for r,m in res.items() if m<CUT)
    return cont,res

for sub,lab in [('final_100_designs','REFOLD (Boltz-2) top-100'),
                ('final_100_designs/before_refolding','GENERADO (pre-refold) top-100')]:
    files=sorted(glob.glob(f'{SP}/estructuras/{sub}/rank*.cif'))
    P('\n'+'='*100); P(f'11. CONTACTOS REALES RECALCULADOS — {lab}  (heavy-atom < {CUT} A)'); P('='*100)
    P(f'n ficheros: {len(files)}')
    rows=[]
    for f in files:
        r=contacts(f)
        if r is None: continue
        cont,res=r
        inter=cont&EPI
        rows.append(dict(fichero=os.path.basename(f), n_contactos=len(cont),
            n_epitopo_contactado=len(inter), recall_epitopo=len(inter)/16,
            precision_epitopo=len(inter)/len(cont) if cont else 0.0,
            centro_contacto_local=float(np.mean(sorted(cont))) if cont else np.nan,
            min_dist_epitopo=min([res[r_] for r_ in EPI if r_ in res], default=np.nan)))
    df=pd.DataFrame(rows)
    df.to_csv(f'{SP}/contactos_{sub.replace("/","_")}.csv',index=False)
    P(df[['n_contactos','n_epitopo_contactado','recall_epitopo','precision_epitopo','min_dist_epitopo']]
      .describe(percentiles=[.25,.5,.75,.9]).round(3).to_string())
    P(f'\nDiseños con >=1 residuo del epitopo en contacto real (<{CUT} A): '
      f'{int((df.n_epitopo_contactado>=1).sum())}/{len(df)}')
    for k in [1,4,8,12,16]:
        P(f'  >= {k:2d}/16 residuos del epitopo contactados: {int((df.n_epitopo_contactado>=k).sum()):3d}/{len(df)}')
    P(f'\nPrecision (fraccion de la interfaz que cae DENTRO del epitopo pedido):')
    P(f'  media {df.precision_epitopo.mean():.3f}  mediana {df.precision_epitopo.median():.3f}  max {df.precision_epitopo.max():.3f}')
    P(f'  diseños con precision >=0.5: {int((df.precision_epitopo>=0.5).sum())}/{len(df)}')
    P(f'  diseños con precision ==0 (interfaz totalmente fuera del epitopo): {int((df.precision_epitopo==0).sum())}/{len(df)}')
    # residuos de diana mas contactados
    allc={}
    for f in files:
        r=contacts(f)
        if r is None: continue
        for res_ in r[0]: allc[res_]=allc.get(res_,0)+1
    top=sorted(allc.items(),key=lambda x:-x[1])[:25]
    P(f'\nResiduos de la diana mas contactados (local -> Spike, nº de diseños de {len(files)}):')
    P('  '+', '.join(f'{k}->S{k+OFF}:{v}' + ('*' if k in EPI else '') for k,v in top))
    P('  (* = dentro del epitopo pedido Spike 479-494)')
P('\nNOTA: recall_epitopo = fraccion de los 16 residuos del epitopo pedido que quedan a <5 A')
P('de algun atomo pesado del binder. precision_epitopo = fraccion de TODOS los residuos de la')
P('diana contactados que pertenecen al epitopo pedido (mide si ataca SOLO donde se le pide).')
OUT.flush()
