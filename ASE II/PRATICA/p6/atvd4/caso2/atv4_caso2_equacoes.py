# -*- coding: utf-8 -*-
"""
Pratica 6 - ASE II - Atividade 4 - Caso 2
Falta trifasica no MEIO da linha 3-4 (representado pela barra 6); falta eliminada retirando a linha 3-4.
Calcula a rede reduzida (pre-falta, falta, pos-falta), imprime Pe1 e Pe2 (equacoes de oscilacao)
e estima o tempo critico por RK4 (referencia para o Anatem).
A linha 3-4 e dividida em duas metades (3-6 e 6-4), igual ao caso do Anarede/Anatem.
Dados dos geradores usados direto na base de 100 MVA, como no Grainger.
"""
import numpy as np
pol=lambda m,a: m*np.exp(1j*np.radians(a))
V={1:pol(1.030,8.88),2:pol(1.020,6.38),3:pol(1.0,0),4:pol(1.018,4.68),5:pol(1.011,2.27)}
Sg={1:3.5+0.712j,2:1.85+0.298j}; Sl={4:1.0+0.44j,5:0.5+0.16j}
Xd={1:0.067,2:0.1}; H={1:11.2,2:8.0}
E={k:V[k]+1j*Xd[k]*np.conj(Sg[k]/V[k]) for k in Sg}
Em=np.array([abs(E[1]),abs(E[2]),1.0]); d0=np.array([np.angle(E[1]),np.angle(E[2])])
base=[(1,4,0,0.022,0),(2,5,0,0.040,0),(3,5,0.008,0.047,0.098),(3,5,0.008,0.047,0.098),(4,5,0.018,0.110,0.226)]
l34=[(3,4,0.007,0.040,0.082)]
l34_split=[(3,6,0.0035,0.020,0.041),(6,4,0.0035,0.020,0.041)]
N=9 # nodes 1..6 buses, 7=E1',8=E2'
def ybus(rams,fault=None):
    Y=np.zeros((N,N),complex)
    def liga(a,b,y,bsh=0):
        Y[a,a]+=y+1j*bsh/2;Y[b,b]+=y+1j*bsh/2;Y[a,b]-=y;Y[b,a]-=y
    for a,b,r,x,B in rams: liga(a,b,1/complex(r,x),B)
    liga(7,1,1/(1j*Xd[1]));liga(8,2,1/(1j*Xd[2]))
    for b,s in Sl.items(): Y[b,b]+=np.conj(s)/abs(V[b])**2
    if fault: Y[fault,fault]+=1e10
    return Y
def kron(Y,elim):
    m=[7,8,3]
    return Y[np.ix_(m,m)]-Y[np.ix_(m,elim)]@np.linalg.solve(Y[np.ix_(elim,elim)],Y[np.ix_(elim,m)])
Yr={"pre (linha inteira)":kron(ybus(base+l34),[1,2,4,5]),
    "pre":kron(ybus(base+l34_split),[1,2,4,5,6]),
    "falta":kron(ybus(base+l34_split,6),[1,2,4,5,6]),
    "pos":kron(ybus(base),[1,2,4,5])}
def eqs(Y):
    out=[]
    for i in range(2):
        s=f"{Em[i]**2*Y[i,i].real:.4f}"
        for j in range(3):
            if j==i: continue
            C=Em[i]*Em[j]*abs(Y[i,j]); g=np.degrees(np.angle(Y[i,j]))-90
            if C<1e-6: continue
            arg=f"d{i+1}-d{j+1}" if j<2 else f"d{i+1}"
            s+=f" + {C:.4f} sen({arg} {'-' if g>=0 else '+'} {abs(g):.3f}°)"
        out.append(s)
    return out
WS=2*np.pi*60;PM=np.array([3.5,1.85]);HH=np.array([11.2,8.0])
def pe(d,Y):
    Ev=Em*np.exp(1j*np.array([d[0],d[1],0.]));return (Ev*np.conj(Y@Ev)).real[:2]
for k,Y in Yr.items():
    print("==",k); [print("  Pe%d ="%(i+1),e) for i,e in enumerate(eqs(Y))]; print("  Pe(d0)=",pe(d0,Y))
def f(x,Y): return np.concatenate([x[2:],WS/(2*HH)*(PM-pe(x[:2],Y))])
def rk(x,h,Y):
    k1=f(x,Y);k2=f(x+h/2*k1,Y);k3=f(x+h/2*k2,Y);k4=f(x+h*k3,Y);return x+h/6*(k1+2*k2+2*k3+k4)
def sim(tc,h=2e-4,T=3.0):
    x=np.array([*d0,0,0]);t=0;mx=np.degrees(x[:2])
    n=max(1,int(round(tc/h)));hh=tc/n
    for _ in range(n): x=rk(x,hh,Yr["falta"]);mx=np.maximum(mx,np.degrees(x[:2]))
    t=tc
    while t<T:
        x=rk(x,h,Yr["pos"]);t+=h;mx=np.maximum(mx,np.degrees(x[:2]))
        if mx.max()>360: return False,mx
    return True,mx
lo,hi=0,1.0
while hi-lo>1e-4:
    m=(lo+hi)/2
    if sim(m)[0]: lo=m
    else: hi=m
print("tcr caso 2 ~ entre %.4f e %.4f"%(lo,hi))
for tc in (0.15,0.2,round(lo-0.005,3),round(hi+0.005,3)): print(tc,sim(tc))
