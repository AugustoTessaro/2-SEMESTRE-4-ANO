import numpy as np,os,json
def rd(p):
    L=open(p,encoding='latin1').read().split('\n')
    n=int(L[0]); names=[l.strip() for l in L[1:1+n]]
    v=np.array(' '.join(L[1+n:]).split(),float).reshape(-1,n)
    # remove duplicated time stamps (keep last = post-event)
    t=v[:,0]; keep=np.r_[t[1:]!=t[:-1],True]
    return names,v[keep]
base=os.path.expanduser('~/mnt/5.2/5.2.1')
TE={'E1_ABRE_LT9-14':1.0,'E2_FT_LT2-4_150ms':1.15}
BUS=[1,2,3,4,5,6,8,9,10,11,12,13,14]
LIM={'13-14':50,'6-13':50,'4-5':100,'9-14':50,'2-4':100,'9-10':50,'12-13':50,'1-2':100,'1-5':100}
out={}
for ev,te in TE.items():
  for c in ['ANTES','DEPOIS']:
    nm,d=rd(f'{base}/{c}/plt/{ev}.plt'); t=d[:,0]
    col=lambda tag,k:[i for i,n in enumerate(nm) if n.split()[0]==tag][k]
    # angle rel bus1
    gens=[2,3,6,8,12,14]; dl={g:d[:,col('DELT',i)] for i,g in enumerate(gens)}
    V={b:d[:,col('VOLT',b-1)] for b in range(1,15)}
    F={g:d[:,col('FMAQ',i)] for i,g in enumerate([1,2,3,6,8,12,14])}
    # oscillation period: from relative angle of most disturbed machine
    dev={g:np.abs(dl[g]-dl[g][0]).max() for g in gens}; gm=max(dev,key=dev.get)
    x=dl[gm][t>te]; tt=t[t>te]; xf=x-x[-1]
    # find local extrema
    ext=[i for i in range(1,len(x)-1) if (x[i]-x[i-1])*(x[i+1]-x[i])<0 and abs(xf[i])>0.05]
    T1=tt[ext[1]]-tt[ext[0]] if len(ext)>1 else None
    T= 2*T1 if T1 else 1.0
    t1=te+ (tt[ext[1]]-te if len(ext)>1 else 1.0)  # end of first swing: second extremum (return)
    r={'ger_mais_perturbado':gm,'ext_t':[round(tt[i],3) for i in ext[:6]],'T_osc':round(T,3),'fim_1a_osc':round(t1,3)}
    # synchronism
    r['ang_ini']={g:round(dl[g][0],2) for g in gens}
    r['ang_max']={g:round(dl[g][t>te].max(),2) for g in gens}
    r['ang_min']={g:round(dl[g][t>te].min(),2) for g in gens}
    r['ang_fin']={g:round(dl[g][-1],2) for g in gens}
    r['ang_pp_ult5s']={g:round(np.ptp(dl[g][t>=t[-1]-5]),3) for g in gens}
    r['max_abs_rel']=round(max(np.abs(dl[g]).max() for g in gens),2)
    # frequency
    r['f_min']={g:round(F[g][t>te].min(),4) for g in F}; r['f_max']={g:round(F[g][t>te].max(),4) for g in F}
    r['f_fin']={g:round(F[g][-1],4) for g in F}
    # voltages
    fault=(t>=1.0)&(t<te) if te>1.0 else np.zeros_like(t,bool)
    w1=(t>te)&(t<=t1); w2=t>t1; wd=t>=te+10
    rv={}
    for b in BUS:
        v=V[b]
        rv[b]=dict(Vini=round(v[0],4),Vfalta_min=round(v[fault].min(),4) if fault.any() else None,
          Vmin_1a=round(v[w1].min(),4),Vmin_dem=round(v[w2].min(),4),Vfin=round(v[-1],4),
          dV=round(v[-1]-v[0],4),pp10=round(np.ptp(v[wd]),5),Vmax_pos=round(v[t>te].max(),4))
    r['V']=rv
    fl={}
    for i,n in enumerate(nm):
      if n.startswith('FLXM'):
        w=n.split(); key=f'{w[1]}-{[x for x in w[2:] if x.isdigit()][0]}'
        fl[key]=dict(S0=round(d[0,i],1),Smax_pos=round(d[t>te+0.0,i].max(),1),Sfin=round(d[-1,i],1),lim=LIM[key],
                     pct_fin=round(100*d[-1,i]/LIM[key],1))
    r['fluxo']=fl
    out[f'{ev}|{c}']=r
json.dump(out,open(os.path.expanduser('~/work/crit.json'),'w'),indent=1,default=float)
for k,r in out.items():
    print('=====',k,'| ger',r['ger_mais_perturbado'],'ext',r['ext_t'],'fim1a',r['fim_1a_osc'],'T',r['T_osc'])
    print(' ang ini',r['ang_ini']); print(' ang max',r['ang_max']); print(' ang min',r['ang_min']); print(' ang fin',r['ang_fin']); print(' pp ult5s',r['ang_pp_ult5s'])
    print(' fmin',r['f_min']); print(' fmax',r['f_max']); print(' ffin',r['f_fin'])
    print(' bus  Vini  Vfalta  Vmin1a  VminDem  Vfin   dV    pp10   Vmaxpos')
    for b,v in r['V'].items(): print(f"  {b:2d} {v['Vini']:.4f} {v['Vfalta_min'] or 0:.4f} {v['Vmin_1a']:.4f} {v['Vmin_dem']:.4f} {v['Vfin']:.4f} {v['dV']:+.4f} {v['pp10']:.5f} {v['Vmax_pos']:.4f}")
    for k2,f in r['fluxo'].items(): print('  ',k2,f)
