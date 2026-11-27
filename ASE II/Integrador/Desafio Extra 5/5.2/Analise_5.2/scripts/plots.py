import numpy as np,os,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open(os.path.expanduser('~/work/crit.py')).read().split('base=')[0])
base=os.path.expanduser('~/mnt/5.2/5.2.1'); OUT=os.path.expanduser('~/mnt/5.2/Analise_5.2/figuras')
PAL=['#257de0','#f47230','#00b894','#f2a600','#ee7fac','#008a17','#6a4ce2','#ec4647']
TXT='#191919'; TX2='#666666'; GRID='#e6e6e6'
plt.rcParams.update({'font.size':8.5,'axes.edgecolor':'#bdbdbd','axes.labelcolor':TXT,'xtick.color':TX2,'ytick.color':TX2,
  'axes.grid':True,'grid.color':GRID,'grid.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
  'lines.linewidth':1.3,'legend.frameon':False,'figure.dpi':110,'savefig.dpi':200,'font.family':'DejaVu Sans'})
EV={'E1_ABRE_LT9-14':(1.0,'E1 – abertura da LT 9–14 em t = 1 s'),'E2_FT_LT2-4_150ms':(1.15,'E2 – curto monofásico na LT 2–4 (junto à barra 4), eliminado em 150 ms')}
D={}
for ev in EV:
  for c in ['ANTES','DEPOIS']:
    nm,d=rd(f'{base}/{c}/plt/{ev}.plt'); D[ev,c]=(nm,d)
def get(ev,c,tag,k):
    nm,d=D[ev,c]; ii=[i for i,n in enumerate(nm) if n.split()[0]==tag]; return d[:,0],d[:,ii[k]]
GEN_ALL=[1,2,3,6,8,12,14]; GEN_REL=[2,3,6,8,12,14]
CG={g:PAL[i] for i,g in enumerate(GEN_ALL)}   # cor fixa por máquina
def pair(ev,fname,title,ylabel,series,xlim=(0,10),hlines=(),shade=None,ylim=None):
    fig,ax=plt.subplots(1,2,figsize=(7.2,3.0),sharey=True)
    for a,c in zip(ax,['ANTES','DEPOIS']):
        for it in series:
            lab,tag,k,col=it[:4]; sc=it[4] if len(it)>4 else 1
            t,y=get(ev,c,tag,k); a.plot(t,y*sc,color=col,label=lab)
        for yv,ls,txt in hlines:
            a.axhline(yv,color='#8c8c8c',ls=ls,lw=0.9); a.text(xlim[1],yv,' '+txt,va='center',ha='left',fontsize=7,color=TX2,clip_on=False)
        if shade: a.axvspan(*shade,color='#f0f0f0',zorder=0,lw=0)
        a.set_title('Antes da conexão' if c=='ANTES' else 'Depois da conexão',fontsize=9,color=TXT,loc='left')
        a.set_xlim(*xlim); a.set_xlabel('Tempo (s)')
        if ylim: a.set_ylim(*ylim)
    ax[0].set_ylabel(ylabel)
    h,l=ax[0].get_legend_handles_labels()
    fig.legend(h,l,loc='lower center',ncol=min(len(l),7),fontsize=7.5,bbox_to_anchor=(0.5,-0.01))
    fig.suptitle(title,fontsize=9.5,color=TXT,x=0.01,ha='left')
    fig.tight_layout(rect=(0,0.08,0.97,0.95)); fig.savefig(f'{OUT}/{fname}.png',bbox_inches='tight'); plt.close(fig)
for ev,(te,tit) in EV.items():
    e=ev[:2]
    pair(ev,f'{e}_1_angulos','%s\nÂngulo do rotor em relação à máquina da barra 1'%tit,'δ − δ₁ (graus)',
         [(f'G{g}','DELT',i,CG[g]) for i,g in enumerate(GEN_REL)])
    pair(ev,f'{e}_2_frequencias','%s\nFrequência (velocidade) das máquinas'%tit,'Frequência (Hz)',
         [(f'G{g}','FMAQ',i,CG[g]) for i,g in enumerate(GEN_ALL)],xlim=(0,20),hlines=[(0.9917*60,':','ERAC 1º estágio')])
    FL=['13-14','6-13','4-5','9-14','2-4']
    nm,_=D[ev,'ANTES']; fk=[n for n in nm if n.startswith('FLXM')]
    keys=[f"{n.split()[1]}-{[x for x in n.split()[2:] if x.isdigit()][0]}" for n in fk]
    LIM={'13-14':50,'6-13':50,'4-5':100,'9-14':50,'2-4':100}
    pair(ev,f'{e}_5_fluxos','%s\nCarregamento das LTs (%% do limite)'%tit,'Carregamento (%)',
         [(f'LT {k} ({LIM[k]} MVA)','FLXM',keys.index(k),PAL[i],100/LIM[k]) for i,k in enumerate(FL)],hlines=[(100,'--','100 %')])

# ---- comparação direta antes x depois (variáveis-chave)
CA,CD='#257de0','#f47230'
def comp(ev,fname,title,panels,xlim=(0,10)):
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.0))
    for a,(ttl,ylab,tag,k,sc,hl) in zip(axs.flat,panels):
        for c,col,lab in [('ANTES',CA,'Antes'),('DEPOIS',CD,'Depois')]:
            t,y=get(ev,c,tag,k); a.plot(t,y*sc,color=col,label=lab)
        for yv,txt in hl: a.axhline(yv,color='#8c8c8c',ls='--',lw=0.9); a.text(xlim[1],yv,' '+txt,va='center',fontsize=7,color=TX2,clip_on=False)
        a.set_title(ttl,fontsize=8.5,loc='left',color=TXT); a.set_ylabel(ylab); a.set_xlim(*xlim)
    for a in axs[1]: a.set_xlabel('Tempo (s)')
    h,l=axs[0,0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=2,bbox_to_anchor=(0.5,-0.005))
    fig.suptitle(title,fontsize=9.5,x=0.01,ha='left',color=TXT)
    fig.tight_layout(rect=(0,0.04,0.96,0.95)); fig.savefig(f'{OUT}/{fname}.png',bbox_inches='tight'); plt.close(fig)
def fidx(ev,key):
    nm,_=D[ev,'ANTES']; fk=[n for n in nm if n.startswith('FLXM')]
    return [f"{n.split()[1]}-{[x for x in n.split()[2:] if x.isdigit()][0]}" for n in fk].index(key)
comp('E1_ABRE_LT9-14','E1_0_comparacao','E1 – abertura da LT 9–14: antes × depois da conexão',[
 ('Ângulo G14 − G1','graus','DELT',5,1,[]),('Frequência G14','Hz','FMAQ',6,1,[]),
 ('Tensão barra 13','pu','VOLT',12,1,[]),('Carregamento LT 13–14','% de 50 MVA','FLXM',fidx('E1_ABRE_LT9-14','13-14'),2,[(100,'100 %')])])
comp('E2_FT_LT2-4_150ms','E2_0_comparacao','E2 – curto monofásico na LT 2–4: antes × depois da conexão',[
 ('Ângulo G12 − G1','graus','DELT',4,1,[]),('Ângulo G14 − G1','graus','DELT',5,1,[]),
 ('Tensão barra 4','pu','VOLT',3,1,[(0.8,'0,80 pu')]),('Carregamento LT 4–5','% de 100 MVA','FLXM',fidx('E2_FT_LT2-4_150ms','4-5'),1,[(100,'100 %')])],xlim=(0,6))

def volt(ev,fname,title,buses,cols,te):
    fig,ax=plt.subplots(2,2,figsize=(7.2,5.2),sharey='row')
    for j,c in enumerate(['ANTES','DEPOIS']):
        for b,col in zip(buses,cols):
            t,y=get(ev,c,'VOLT',b-1)
            ax[0,j].plot(t,y,color=col,label=f'B{b}'); ax[1,j].plot(t,y,color=col)
        a=ax[0,j]; a.set_xlim(0.5,4.5); a.set_ylim(0.55,1.12)
        for yv,txt in [(0.80,'0,80 pu (b)'),(0.60,'0,60 pu (a)')]:
            a.axhline(yv,color='#8c8c8c',ls='--',lw=0.9); a.text(4.5,yv,' '+txt,va='center',fontsize=7,color=TX2,clip_on=False)
        a.set_title(('Antes' if c=='ANTES' else 'Depois')+' da conexão – zoom no transitório',fontsize=8.5,loc='left',color=TXT)
        a=ax[1,j]; a.set_xlim(0,20); a.set_ylim(0.93,1.08)
        a.axvspan(te+10,20,color='#f0f0f0',zorder=0,lw=0); a.text(te+10.3,0.935,'janela do critério (d)',fontsize=7,color=TX2)
        a.set_title(('Antes' if c=='ANTES' else 'Depois')+' – simulação completa (escala ampliada)',fontsize=8.5,loc='left',color=TXT); a.set_xlabel('Tempo (s)')
    ax[0,0].set_ylabel('Tensão (pu)'); ax[1,0].set_ylabel('Tensão (pu)')
    h,l=ax[0,0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=len(l),bbox_to_anchor=(0.5,-0.005),fontsize=7.5)
    fig.suptitle(title,fontsize=9.5,x=0.01,ha='left',color=TXT)
    fig.tight_layout(rect=(0,0.04,0.95,0.95)); fig.savefig(f'{OUT}/{fname}.png',bbox_inches='tight'); plt.close(fig)
for ev,(te,tit) in EV.items():
    e=ev[:2]
    volt(ev,f'{e}_3_tensao_geracao',tit+'\nTensão nas barras de geração',GEN_ALL,[CG[g] for g in GEN_ALL],te)
    volt(ev,f'{e}_4_tensao_carga',tit+'\nTensão nas barras de carga (sem geração)',[4,5,9,10,11,13],PAL[:6],te)
