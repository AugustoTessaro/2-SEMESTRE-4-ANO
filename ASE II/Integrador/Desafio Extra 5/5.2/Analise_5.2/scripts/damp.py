import numpy as np,os
exec(open(os.path.expanduser('~/work/crit.py')).read().split('base=')[0])
base=os.path.expanduser('~/mnt/5.2/5.2.1')
for ev,te in {'E1_ABRE_LT9-14':1.0,'E2_FT_LT2-4_150ms':1.15}.items():
  for c in ['ANTES','DEPOIS']:
    nm,d=rd(f'{base}/{c}/plt/{ev}.plt'); t=d[:,0]
    ii=[i for i,n in enumerate(nm) if n.startswith('DELT')][:6]; gens=[2,3,6,8,12,14]
    res=[]
    for k,i in enumerate(ii):
        x=d[:,i]; m=t>te; tt=t[m]; y=x[m]-x[-1]
        pk=[j for j in range(1,len(y)-1) if y[j]>y[j-1] and y[j]>=y[j+1] and y[j]>0]
        pk=[j for j in pk if y[j]>0.02*np.abs(y).max()]
        if len(pk)>=3:
            A=y[pk]; T=np.diff(tt[pk])[:3].mean(); delta=np.log(A[0]/A[2])/2
            zeta=delta/np.sqrt(4*np.pi**2+delta**2)
            # settling time: last time |y|>5% of max
            ts=tt[np.where(np.abs(y)>0.05*np.abs(y).max())[0][-1]]
            res.append((gens[k],round(1/T,2),round(zeta*100,1),round(ts-te,2),round(np.abs(y).max(),1)))
    pe=[i for i,n in enumerate(nm) if n.startswith('PELE')]; pm=[i for i,n in enumerate(nm) if n.startswith('PMEC')]
    print(ev,c,'(ger, f_modo Hz, zeta %, t_acomod 5% s, amplitude °):',res)
    print('   Pe total ini/fin', round(d[0,pe].sum(),1), round(d[-1,pe].sum(),1),' Pm', round(d[0,pm].sum(),1), round(d[-1,pm].sum(),1))
