# -*- coding: utf-8 -*-
"""
Prática 6 - ASE II - Atividade 4 - Caso 2 - alínea b)
Estabilidade multimáquinas: equações de oscilação resolvidas por Runge-Kutta
de 4ª ordem (RK4) e comparação com o Anatem (alínea a).

Caso 2: falta trifásica no MEIO da linha 3-4. A falta é eliminada retirando
a linha 3-4. Como no Anarede/Anatem, a linha 3-4 é dividida em duas metades
(3-6 e 6-4) e o ponto médio é a barra 6. Após a retirada da linha, a barra 6
fica isolada e sai da rede (equivale à opção ILHA do Anatem).

Sistema (Grainger, Ex. 16.9), base 100 MVA:
   G1 --(X'd1)-- [1] --T(0,022)-- [4] --- 6-4 --[6]-- 3-6 --- [3] barra infinita
                                    |                        ||
                                linha 4-5             linhas 3-5 (duas)
                                    |                        ||
   G2 --(X'd2)-- [2] --T(0,040)-- [5] ------------------------

Representação clássica:
  - cada gerador = E' constante atrás de X'd;
  - cargas = admitâncias constantes (calculadas com as tensões pré-falta);
  - a rede é reduzida (Kron) aos nós internos E1', E2' e à barra infinita.

Equação de oscilação de cada máquina i (δ em rad, Δω em rad/s elétricos):
   dδi/dt  = Δωi
   dΔωi/dt = (ωs / 2Hi) · (Pmi - Pei)
   Pei = Ei²·Gii + Σj Ei·Ej·|Yij|·cos(θij - δi + δj)

Como usar:  python atv4_caso2_rk4.py      (precisa de numpy e matplotlib)
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1) DADOS
# =============================================================================
F = 60.0
WS = 2 * np.pi * F

def polar(mod, ang_graus):
    return mod * np.exp(1j * np.radians(ang_graus))

# Fluxo de potência do enunciado (Tabela 2)
V = {1: polar(1.030, 8.88), 2: polar(1.020, 6.38), 3: polar(1.000, 0.0),
     4: polar(1.018, 4.68), 5: polar(1.011, 2.27)}
S_GER = {1: 3.50 + 0.712j, 2: 1.85 + 0.298j}       # geração (pu)
S_CARGA = {4: 1.00 + 0.44j, 5: 0.50 + 0.16j}       # cargas (pu)

# Geradores: usados direto na base de 100 MVA, como no Grainger
XD = {1: 0.067, 2: 0.10}
H = {1: 11.2, 2: 8.0}

# Ramos: (de, para, R, X, B total)
RAMOS = [(1, 4, 0.000, 0.022, 0.000),     # transformador 1-4
         (2, 5, 0.000, 0.040, 0.000),     # transformador 2-5
         (3, 6, 0.0035, 0.020, 0.041),    # metade 3-6 da linha 3-4 (retirada)
         (6, 4, 0.0035, 0.020, 0.041),    # metade 6-4 da linha 3-4 (retirada)
         (3, 5, 0.008, 0.047, 0.098),     # linha 3-5 circuito 1
         (3, 5, 0.008, 0.047, 0.098),     # linha 3-5 circuito 2
         (4, 5, 0.018, 0.110, 0.226)]     # linha 4-5

BARRA_FALTA = 6                       # ponto médio da linha 3-4
RAMOS_RETIRADOS = [(3, 6), (6, 4)]
T_FALTA = 0.1          # instante da falta (igual ao Anatem)
T_FINAL = 3.0
H_PASSO = 1e-4

PASTA = os.path.dirname(os.path.abspath(__file__))     # caso2/B (onde está este script)
PASTA_PLT = os.path.dirname(PASTA)                     # caso2   (onde estão os .plt do Anatem)

# =============================================================================
# 2) TENSÕES INTERNAS E'
# =============================================================================
E = {}
for g in (1, 2):
    I = np.conj(S_GER[g] / V[g])           # corrente do gerador
    E[g] = V[g] + 1j * XD[g] * I           # E' = Vt + jX'd·I
E_MOD = np.array([abs(E[1]), abs(E[2]), 1.0])          # |E1'|, |E2'|, |V3|
DELTA0 = np.array([np.angle(E[1]), np.angle(E[2])])    # ângulos iniciais (rad)
PM = np.array([S_GER[1].real, S_GER[2].real])          # Pm = Pe pré-falta
HH = np.array([H[1], H[2]])

# =============================================================================
# 3) YBUS AUMENTADA E REDUÇÃO DE KRON
#    nós: 1..6 = barras, 7 = E1', 8 = E2'
# =============================================================================
NO_E1, NO_E2 = 7, 8

def ybus(ramos, barra_curto=None):
    Y = np.zeros((9, 9), complex)
    def liga(a, b, y, b_sh=0.0):
        Y[a, a] += y + 1j * b_sh / 2
        Y[b, b] += y + 1j * b_sh / 2
        Y[a, b] -= y
        Y[b, a] -= y
    for a, b, r, x, bsh in ramos:
        liga(a, b, 1 / complex(r, x), bsh)
    liga(NO_E1, 1, 1 / (1j * XD[1]))
    liga(NO_E2, 2, 1 / (1j * XD[2]))
    for b, s in S_CARGA.items():                   # carga -> admitância constante
        Y[b, b] += np.conj(s) / abs(V[b]) ** 2
    if barra_curto is not None:                    # curto franco: admitância "infinita"
        Y[barra_curto, barra_curto] += 1e10
    return Y

def kron(Y, eliminar=(1, 2, 4, 5, 6)):
    """Reduz aos nós E1', E2' e barra infinita (3), nesta ordem."""
    m = [NO_E1, NO_E2, 3]
    e = list(eliminar)
    return Y[np.ix_(m, m)] - Y[np.ix_(m, e)] @ np.linalg.solve(Y[np.ix_(e, e)], Y[np.ix_(e, m)])

ramos_pos = [r for r in RAMOS if (r[0], r[1]) not in RAMOS_RETIRADOS]
YRED = {"pré-falta": kron(ybus(RAMOS)),
        "falta": kron(ybus(RAMOS, BARRA_FALTA)),
        "pós-falta": kron(ybus(ramos_pos), eliminar=(1, 2, 4, 5))}   # barra 6 isolada: fora

# =============================================================================
# 4) EQUAÇÕES DE OSCILAÇÃO (impressão no formato do livro)
# =============================================================================
def imprime_equacoes():
    print("\nTENSÕES INTERNAS")
    for g in (1, 2):
        print(f"  E{g}' = {abs(E[g]):.4f} ∠ {np.degrees(np.angle(E[g])):.2f}°")
    print("\nEQUAÇÕES DE OSCILAÇÃO  (δ em graus):  (H/180f)·d²δ/dt² = Pm − Pe")
    print(f"  G1: ({H[1]}/{180*F:.0f})·d²δ1/dt² = {PM[0]} − Pe1")
    print(f"  G2: ({H[2]}/{180*F:.0f})·d²δ2/dt² = {PM[1]} − Pe2")
    nomes = ["δ1", "δ2", "δ3"]
    for etapa, Y in YRED.items():
        print(f"\n  {etapa.upper()}")
        for i in range(2):
            termos = [f"{E_MOD[i]**2 * Y[i, i].real:.4f}"]
            for j in range(3):
                if j == i:
                    continue
                C = E_MOD[i] * E_MOD[j] * abs(Y[i, j])
                if C < 1e-6:
                    continue
                gama = np.degrees(np.angle(Y[i, j])) - 90
                arg = f"{nomes[i]}−{nomes[j]}" if j < 2 else nomes[i]
                termos.append(f"{C:.4f}·sen({arg} {'−' if gama >= 0 else '+'} {abs(gama):.3f}°)")
            print(f"    Pe{i+1} = " + " + ".join(termos))

# =============================================================================
# 5) RUNGE-KUTTA DE 4ª ORDEM
# =============================================================================
def pe(delta, Y):
    """Potência elétrica de cada gerador para os ângulos atuais."""
    Ev = E_MOD * np.exp(1j * np.array([delta[0], delta[1], 0.0]))
    return (Ev * np.conj(Y @ Ev)).real[:2]

def derivadas(x, Y):
    """x = [δ1, δ2, Δω1, Δω2]"""
    return np.concatenate([x[2:], WS / (2 * HH) * (PM - pe(x[:2], Y))])

def passo_rk4(x, h, Y):
    k1 = derivadas(x, Y)
    k2 = derivadas(x + h / 2 * k1, Y)
    k3 = derivadas(x + h / 2 * k2, Y)
    k4 = derivadas(x + h * k3, Y)
    return x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

def integra(x, t, t_fim, Y, h, ts, xs):
    n = max(1, int(np.ceil((t_fim - t) / h - 1e-9)))
    hh = (t_fim - t) / n                       # último passo cai exatamente em t_fim
    for _ in range(n):
        x = passo_rk4(x, hh, Y)
        t += hh
        ts.append(t)
        xs.append(x)
        if np.degrees(x[:2]).max() > 1000:     # mesmo limite do Anatem
            break
    return x, t

def simular(tc, h=H_PASSO):
    x = np.array([DELTA0[0], DELTA0[1], 0.0, 0.0])
    ts, xs = [0.0], [x]
    x, t = integra(x, 0.0, T_FALTA, YRED["pré-falta"], h, ts, xs)
    x, t = integra(x, t, T_FALTA + tc, YRED["falta"], h, ts, xs)
    if np.degrees(x[:2]).max() < 1000:
        x, t = integra(x, t, T_FINAL, YRED["pós-falta"], h, ts, xs)
    xs = np.array(xs)
    return np.array(ts), np.degrees(xs[:, 0]), np.degrees(xs[:, 1])

def estavel(tc, h=2e-4):
    """Mesmo critério do Anatem: instável se algum ângulo passar de 360°."""
    _, d1, d2 = simular(tc, h)
    return max(d1.max(), d2.max()) < 360

def tempo_critico(tol=1e-4):
    lo, hi = 0.0, 0.5
    while hi - lo > tol:
        meio = (lo + hi) / 2
        lo, hi = (meio, hi) if estavel(meio) else (lo, meio)
    return lo, hi

# =============================================================================
# 6) LEITURA DOS .PLT DO ANATEM
# =============================================================================
def ler_plt(caminho):
    with open(caminho, encoding="latin1") as f:
        linhas = f.read().split("\n")
    n = int(linhas[0])
    valores = np.array(" ".join(linhas[n + 1:]).split(), dtype=float)
    d = valores[: len(valores) // n * n].reshape(-1, n)
    return d[:, 0], d[:, 1], d[:, 2]                # t, DELT 1, DELT 2

def tc_do_nome(caminho):
    """atv4_caso2_0_32.plt -> RMCB = 0,32 s -> tc = 0,22 s"""
    nome = os.path.splitext(os.path.basename(caminho))[0]
    rmcb = float(nome.split("caso2_")[1].replace("_", "."))
    return round(rmcb - T_FALTA, 6)

# =============================================================================
# 7) PROGRAMA PRINCIPAL
# =============================================================================
if __name__ == "__main__":
    imprime_equacoes()

    lo, hi = tempo_critico()
    print(f"\nTEMPO CRÍTICO PELO RK4: entre {lo:.5f} s (estável) e {hi:.5f} s (instável)")

    arquivos = sorted(glob.glob(os.path.join(PASTA_PLT, "atv4_caso2_*.plt")), key=tc_do_nome)
    if not arquivos:
        print(f"\n(nenhum .plt do Anatem encontrado em {PASTA_PLT})")
        arquivos = []
    print("\nCOMPARAÇÃO RK4 x ANATEM")
    print(f"{'tc (s)':>7} {'δ1max RK4':>10} {'δ1max Anatem':>13} {'δ2max RK4':>10} "
          f"{'δ2max Anatem':>13} {'RK4':>9} {'Anatem':>9}")

    n = max(1, len(arquivos))
    fig, eixos = plt.subplots(n, 1, figsize=(10, 3.0 * n), sharex=True, squeeze=False)
    for ax, arq in zip(eixos[:, 0], arquivos):
        tc = tc_do_nome(arq)
        t, d1, d2 = simular(tc)
        ta, a1, a2 = ler_plt(arq)
        st = lambda a, b: "estável" if max(a.max(), b.max()) < 360 else "instável"
        print(f"{tc:7.3f} {min(d1.max(), 999):10.1f} {min(a1.max(), 999):13.1f} "
              f"{min(d2.max(), 999):10.1f} {min(a2.max(), 999):13.1f} "
              f"{st(d1, d2):>9} {st(a1, a2):>9}")
        ax.plot(t, d1, "b-", lw=1.6, label="δ1 RK4")
        ax.plot(ta, a1, "c--", lw=1.2, label="δ1 Anatem")
        ax.plot(t, d2, "r-", lw=1.6, label="δ2 RK4")
        ax.plot(ta, a2, "m--", lw=1.2, label="δ2 Anatem")
        ax.axvline(T_FALTA, color="gray", ls=":", lw=1)
        ax.axvline(T_FALTA + tc, color="gray", ls=":", lw=1)
        ax.set_ylim(-60, 400)
        ax.set_ylabel("δ (graus)")
        ax.set_title(f"Caso 2 – falta no meio da linha 3-4, linha retirada, tc = {tc:.3f} s ({st(d1, d2)})")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", ncol=2, fontsize=8)
    eixos[-1, 0].set_xlabel("tempo (s)")
    plt.tight_layout()
    plt.savefig(os.path.join(PASTA, "caso2_rk4_vs_anatem.png"), dpi=150)
    print("\nFigura salva em caso2_rk4_vs_anatem.png")
    plt.show()
