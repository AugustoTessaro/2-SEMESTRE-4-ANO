# -*- coding: utf-8 -*-
"""
Prática 6 - ASE II - Atividade 3
Solução da equação de oscilação pelo método de Runge-Kutta de 4ª ordem (RK4)
e comparação com o Anatem (Atividade 2) e com o Excel (Atividade 1).

Sistema da Figura 1 (máquina x barra infinita):

   E'∠δ --j0,2-- (1) --j0,1-- (3) ==== j0,4 (linha de cima) ==== (2) barra infinita
                                 \\=== j0,2 --(P)-- j0,2 ======//      V2 = 1,05∠0
   V1 = 1,0 (terminal do gerador)

Equação de oscilação (modelo clássico, sem amortecimento):

   dδ/dt = Δω
   dΔω/dt = (ωs / 2H) * (Pm - Pmax * sen δ)

   δ em rad, Δω = desvio de velocidade em rad/s elétricos, ωs = 2π·60.
   Pmax muda conforme a etapa: pré-falta, durante a falta e pós-falta.

Como usar:
   python atividade3_rk4.py
   (precisa de numpy e matplotlib:  pip install numpy matplotlib)
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1) DADOS DO SISTEMA (iguais aos do Excel e do Anatem)
# =============================================================================
F = 60.0                    # frequência (Hz)
WS = 2 * np.pi * F          # velocidade síncrona elétrica (rad/s)
V1 = 1.0                    # tensão no terminal do gerador (pu)
V2 = 1.05                   # tensão da barra infinita (pu)
XD = 0.2                    # reatância transitória X'd (pu)
XT = 0.1                    # transformador (pu)
XL = 0.4                    # cada linha de transmissão (pu)

T_FALTA = 0.1               # instante da falta (s) - igual ao Anatem
T_FINAL = 3.0               # tempo total simulado (s)
H_PASSO = 1e-4              # passo do RK4 (s)

# Casos a comparar: (Pm, H, letra, tc, arquivo .plt do Anatem)
PASTA = os.path.dirname(os.path.abspath(__file__))
ATV2 = os.path.join(PASTA, "..", "atvd2", "VariaçãoPM")
CASOS = [
    (0.7,  3, "a", 0.246,   os.path.join(ATV2, "H3 PM 0.7",  "atv2a_falta_barra3_caso_min.plt")),
    (1.0,  5, "a", 0.22456, os.path.join(ATV2, "H 5 PM 1",   "atv2a_falta_barra3.plt")),
    (1.2, 10, "a", 0.254,   os.path.join(ATV2, "H10 PM 1.2", "atv2a_falta_barra3.plt")),
]


# =============================================================================
# 2) CONDIÇÃO INICIAL E POTÊNCIAS MÁXIMAS
# =============================================================================
def condicao_inicial(Pm):
    """Calcula E' e δ0 a partir do fluxo de potência pré-falta."""
    X12 = XT + XL / 2                         # terminal -> barra infinita (2 linhas em paralelo)
    theta = np.arcsin(Pm * X12 / (V1 * V2))   # ângulo da tensão terminal
    Vt = V1 * np.exp(1j * theta)
    I = (Vt - V2) / (1j * X12)                # corrente do gerador
    E = Vt + 1j * XD * I                      # tensão interna E'
    return abs(E), np.angle(E)                # |E'| e δ0 (rad)


def potencias_maximas(E, letra):
    """Pmax = E'·V2 / X_transferência em cada etapa."""
    X_pre = XD + XT + XL / 2                  # duas linhas em paralelo
    if letra == "a":
        # a) falta trifásica na barra 3: nenhuma potência chega à barra infinita.
        #    Falta eliminada sem retirar linha -> pós-falta = pré-falta.
        X_falta = np.inf
        X_pos = X_pre
    else:
        # b) falta no meio da linha de baixo (ponto P).
        #    Estrela na barra 3: Xa = X'd+Xt (até E'), Xb = XL (linha de cima),
        #    Xc = XL/2 (até o ponto P em curto, ou seja, até a terra).
        #    Transformação estrela-triângulo: X_transf = (XaXb + XbXc + XcXa)/Xc
        Xa, Xb, Xc = XD + XT, XL, XL / 2
        X_falta = (Xa * Xb + Xb * Xc + Xc * Xa) / Xc
        X_pos = XD + XT + XL                  # linha de baixo retirada
    return E * V2 / X_pre, E * V2 / X_falta, E * V2 / X_pos


# =============================================================================
# 3) EQUAÇÃO DE OSCILAÇÃO E PASSO DO RUNGE-KUTTA
# =============================================================================
def derivadas(x, Pm, H, Pmax):
    """x = [δ, Δω]. Retorna [dδ/dt, dΔω/dt]."""
    delta, dw = x
    return np.array([dw, (WS / (2 * H)) * (Pm - Pmax * np.sin(delta))])


def passo_rk4(x, h, Pm, H, Pmax):
    """Um passo do método de Runge-Kutta de 4ª ordem."""
    k1 = derivadas(x, Pm, H, Pmax)
    k2 = derivadas(x + h / 2 * k1, Pm, H, Pmax)
    k3 = derivadas(x + h / 2 * k2, Pm, H, Pmax)
    k4 = derivadas(x + h * k3, Pm, H, Pmax)
    return x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def integrar(x, t, t_fim, Pm, H, Pmax, ts, xs):
    """Integra de t até t_fim com Pmax constante. O último passo é ajustado
    para cair exatamente em t_fim (o instante da eliminação da falta)."""
    n = max(1, int(np.ceil((t_fim - t) / H_PASSO - 1e-9)))
    h = (t_fim - t) / n
    for _ in range(n):
        x = passo_rk4(x, h, Pm, H, Pmax)
        t += h
        ts.append(t)
        xs.append(x)
        if x[0] > np.radians(1000):           # já perdeu o sincronismo há muito tempo
            break
    return x, t


def simular(Pm, H, letra, tc):
    """Simula: pré-falta (0 a 0,1 s) -> falta (duração tc) -> pós-falta."""
    E, d0 = condicao_inicial(Pm)
    P_pre, P_falta, P_pos = potencias_maximas(E, letra)
    x = np.array([d0, 0.0])
    ts, xs = [0.0], [x]
    x, t = integrar(x, 0.0, T_FALTA, Pm, H, P_pre, ts, xs)            # pré-falta
    x, t = integrar(x, t, T_FALTA + tc, Pm, H, P_falta, ts, xs)       # durante a falta
    x, t = integrar(x, t, T_FINAL, Pm, H, P_pos, ts, xs)              # pós-falta
    ts, xs = np.array(ts), np.array(xs)
    delta = xs[:, 0]
    Pmax_t = np.where(ts <= T_FALTA, P_pre, np.where(ts <= T_FALTA + tc, P_falta, P_pos))
    Pe = Pmax_t * np.sin(delta)
    return ts, np.degrees(delta), Pe


# =============================================================================
# 4) TEMPO CRÍTICO POR BISSECÇÃO (usando o próprio RK4)
# =============================================================================
def e_estavel(Pm, H, letra, tc):
    """Primeira oscilação: estável se a velocidade volta a zero (rotor freia)
    antes de δ passar do ponto de equilíbrio instável δu = 180° - δs."""
    E, d0 = condicao_inicial(Pm)
    P_pre, P_falta, P_pos = potencias_maximas(E, letra)
    du = np.pi - np.arcsin(Pm / P_pos)
    x = np.array([d0, 0.0])
    lixo_t, lixo_x = [], []
    x, t = integrar(x, 0.0, tc, Pm, H, P_falta, lixo_t, lixo_x)   # período de falta
    while t < 5.0:
        x = passo_rk4(x, H_PASSO, Pm, H, P_pos)
        t += H_PASSO
        if x[0] > du:
            return False
        if x[1] <= 0:
            return True
    return True


def tempo_critico(Pm, H, letra, tol=1e-7):
    lo, hi = 0.0, 2.0
    while hi - lo > tol:
        meio = (lo + hi) / 2
        if e_estavel(Pm, H, letra, meio):
            lo = meio
        else:
            hi = meio
    return (lo + hi) / 2


# =============================================================================
# 5) LEITURA DO .PLT DO ANATEM
# =============================================================================
def ler_plt(caminho):
    """Formato do .plt: 1ª linha = nº de colunas; depois o nome de cada coluna;
    depois os dados (tempo, DELT, PELE, ...)."""
    if not os.path.exists(caminho):
        return None
    with open(caminho, encoding="latin1") as f:
        linhas = f.read().split("\n")
    n = int(linhas[0])
    dados = np.array([[float(v) for v in l.split()] for l in linhas[n + 1:] if l.strip()])
    return dados[:, 0], dados[:, 1], dados[:, 2]       # t, δ (graus), Pe (MW)


# =============================================================================
# 6) PROGRAMA PRINCIPAL
# =============================================================================
if __name__ == "__main__":
    print("\nTEMPOS CRÍTICOS CALCULADOS PELO RK4")
    print(f"{'Pm':>5} {'H':>4} {'δ0 (°)':>8} {'tcr a) (s)':>11} {'tcr b) (s)':>11}")
    for Pm, H, *_ in CASOS:
        E, d0 = condicao_inicial(Pm)
        print(f"{Pm:5.1f} {H:4d} {np.degrees(d0):8.3f} "
              f"{tempo_critico(Pm, H, 'a'):11.5f} {tempo_critico(Pm, H, 'b'):11.5f}")

    print("\nCOMPARAÇÃO RK4 x ANATEM (mesmo tc usado no Anatem)")
    print(f"{'Caso':<14} {'tc (s)':>8} {'δ(tc) RK4':>10} {'δ(tc) Anatem':>13} "
          f"{'δmax RK4':>9} {'δmax Anatem':>12} {'RK4':>9} {'Anatem':>9}")

    fig, eixos = plt.subplots(len(CASOS), 1, figsize=(10, 3.2 * len(CASOS)), sharex=True)
    for ax, (Pm, H, letra, tc, plt_path) in zip(np.atleast_1d(eixos), CASOS):
        t, d, Pe = simular(Pm, H, letra, tc)
        i = np.searchsorted(t, T_FALTA + tc)
        status_rk = "estável" if d.max() < 180 else "instável"
        ax.plot(t, d, "b-", lw=1.6, label="RK4 (Python)")

        anatem = ler_plt(plt_path)
        if anatem is not None:
            ta, da, _ = anatem
            j = np.searchsorted(ta, T_FALTA + tc - 1e-9)
            status_an = "estável" if da.max() < 180 else "instável"
            ax.plot(ta, da, "r--", lw=1.2, label="Anatem")
            print(f"Pm={Pm} H={H:<4} {tc:8.5f} {d[i]:10.2f} {da[j]:13.2f} "
                  f"{min(d.max(), 999):9.1f} {min(da.max(), 999):12.1f} {status_rk:>9} {status_an:>9}")
        else:
            print(f"Pm={Pm} H={H:<4} {tc:8.5f} {d[i]:10.2f} {'(sem .plt)':>13}")

        ax.axvline(T_FALTA, color="gray", ls=":", lw=1)
        ax.axvline(T_FALTA + tc, color="gray", ls=":", lw=1)
        ax.set_ylim(min(-90, d.min() - 10), 360)
        ax.set_ylabel("δ (graus)")
        ax.set_title(f"Pm = {Pm} pu, H = {H} s, falta na barra 3, tc = {tc} s")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right")
    np.atleast_1d(eixos)[-1].set_xlabel("tempo (s)")
    plt.tight_layout()
    plt.savefig(os.path.join(PASTA, "rk4_vs_anatem.png"), dpi=150)
    print("\nFigura salva em rk4_vs_anatem.png")
    plt.show()
