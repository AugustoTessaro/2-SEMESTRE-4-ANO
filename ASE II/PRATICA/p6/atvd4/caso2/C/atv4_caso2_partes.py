# -*- coding: utf-8 -*-
"""
Prática 6 - ASE II - Atividade 4 - Caso 2 - alínea c)
Falta no meio da linha 3-4 (barra 6), linha 3-4 retirada.
Solução por partes (passo a passo) da curva de oscilação - multimáquinas.

Método (Grainger & Stevenson / slides 40-41):
  - Divide-se o tempo em intervalos Δt. Em cada intervalo, supõe-se que a
    potência de aceleração Pa = Pm - Pe é constante e igual ao valor calculado
    no INÍCIO do intervalo, e que a velocidade é constante no meio do intervalo.
  - Daí resulta, para cada máquina (δ em graus):

        Δδn = Δδn-1 + k · Pa,n-1        k = 180·f·Δt² / H
        δn  = δn-1 + Δδn

  - Descontinuidades (aplicação e eliminação da falta): se ocorrem num ponto
    do passo, usa-se a MÉDIA de Pa antes e depois da descontinuidade.
  - Pe de cada máquina vem da rede reduzida (pré-falta, falta, pós-falta),
    a mesma usada no RK4 da alínea b).

A rede e as funções comuns são importadas de ../B/atv4_caso2_rk4.py.
Estrutura de pastas esperada:
   caso2/             .plt do Anatem
   caso2/B/           atv4_caso2_rk4.py   (alínea b)
   caso2/C/           este script         (alínea c)
Como usar:  python atv4_caso2_partes.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

PASTA = os.path.dirname(os.path.abspath(__file__))     # caso2/C
PASTA_CASO = os.path.dirname(PASTA)                   # caso2 (onde estão os .plt)
sys.path.insert(0, os.path.join(PASTA_CASO, "B"))     # para importar o script da alínea b
import atv4_caso2_rk4 as base          # rede reduzida, E', Pm, H, RK4, leitor de .plt

F, PM, HH = base.F, base.PM, base.HH
Y_PRE, Y_FALTA, Y_POS = base.YRED["pré-falta"], base.YRED["falta"], base.YRED["pós-falta"]
T_FALTA, T_FINAL = base.T_FALTA, base.T_FINAL


def pa(delta_graus, Y):
    """Potência de aceleração Pa = Pm - Pe de cada máquina (pu)."""
    return PM - base.pe(np.radians(delta_graus), Y)


def solucao_por_partes(tc, dt, t_final=T_FINAL, tabela=False):
    """Retorna t, δ1, δ2 (graus). A falta ocorre em T_FALTA e é eliminada em T_FALTA+tc.
    Para usar a média nas descontinuidades, tc deve ser múltiplo de Δt."""
    for x, nome in ((T_FALTA, "instante da falta"), (tc, "tc")):
        if abs(x / dt - round(x / dt)) > 1e-6:
            raise ValueError(f"{nome} = {x} s não é múltiplo de Δt = {dt} s")
    k = 180 * F * dt ** 2 / HH                 # constante k de cada máquina
    n_falta = int(round(T_FALTA / dt))
    n_elim = n_falta + int(round(tc / dt))
    n_total = int(round(t_final / dt))

    delta = np.degrees(base.DELTA0)            # δ1, δ2 iniciais (graus)
    d_delta = np.zeros(2)                      # Δδ do intervalo anterior
    ts, d1, d2 = [0.0], [delta[0]], [delta[1]]
    linhas = []

    for n in range(n_total):
        t = n * dt
        # ---- potência de aceleração no ponto n ----
        if n < n_falta:
            Pa_n = pa(delta, Y_PRE)
            obs = "pré-falta"
        elif n == n_falta:                     # aplicação da falta: média
            Pa_n = (pa(delta, Y_PRE) + pa(delta, Y_FALTA)) / 2
            obs = "falta aplicada (média)"
        elif n < n_elim:
            Pa_n = pa(delta, Y_FALTA)
            obs = "durante a falta"
        elif n == n_elim:                      # eliminação da falta: média
            Pa_n = (pa(delta, Y_FALTA) + pa(delta, Y_POS)) / 2
            obs = "falta eliminada (média)"
        else:
            Pa_n = pa(delta, Y_POS)
            obs = "pós-falta"
        if tabela:
            linhas.append((t, *Pa_n, *delta, obs))

        # ---- passo do método ----
        d_delta = d_delta + k * Pa_n
        delta = delta + d_delta
        ts.append(t + dt)
        d1.append(delta[0])
        d2.append(delta[1])
        if delta.max() > 1000:                 # perdeu o sincronismo
            break

    if tabela:
        print(f"\nTABELA - SOLUÇÃO POR PARTES (Δt = {dt} s, tc = {tc} s; "
              f"k1 = {k[0]:.4f}, k2 = {k[1]:.4f})")
        print(f"{'t (s)':>7} {'Pa1':>9} {'Pa2':>9} {'δ1 (°)':>9} {'δ2 (°)':>9}   etapa")
        for t, p1, p2, a1, a2, obs in linhas:
            if t <= T_FALTA + tc + 0.5 + 1e-9:
                print(f"{t:7.3f} {p1:9.4f} {p2:9.4f} {a1:9.2f} {a2:9.2f}   {obs}")
    return np.array(ts), np.array(d1), np.array(d2)


def tc_critico_por_partes(dt, ini=0.200, fim=0.240):
    """Varre tc em múltiplos de Δt e devolve o último estável e o primeiro instável."""
    ultimo_estavel = None
    for m in range(int(round(ini / dt)), int(round(fim / dt)) + 1):
        tc = m * dt
        _, d1, d2 = solucao_por_partes(tc, dt)
        if max(d1.max(), d2.max()) < 360:
            ultimo_estavel = tc
        else:
            return ultimo_estavel, tc
    return ultimo_estavel, None


if __name__ == "__main__":
    # 1) Tabela no formato do livro, com Δt = 0,05 s e tc = 0,20 s
    solucao_por_partes(0.20, 0.05, tabela=True)

    # 2) Tempo crítico pelo método por partes, para alguns Δt
    print("\nTEMPO CRÍTICO PELO MÉTODO POR PARTES")
    for dt in (0.01, 0.005, 0.001, 0.0001):
        e, i = tc_critico_por_partes(dt)
        print(f"  Δt = {dt:<7} -> estável até tc = {e:.4f} s, instável a partir de {i:.4f} s")
    lo, hi = base.tempo_critico()
    print(f"  (RK4: tcr entre {lo:.5f} e {hi:.5f} s; Anatem: entre 0,22 e 0,23 s)")

    # 3) Gráficos: por partes (Δt = 0,05 e 0,01) x RK4 x Anatem
    casos = [0.200, 0.210]
    fig, eixos = plt.subplots(len(casos), 1, figsize=(10, 4 * len(casos)), sharex=True)
    for ax, tc in zip(eixos, casos):
        t, r1, r2 = base.simular(tc)
        ax.plot(t, r1, "b-", lw=1.4, label="δ1 RK4")
        ax.plot(t, r2, "r-", lw=1.4, label="δ2 RK4")
        arq = os.path.join(PASTA_CASO, "atv4_caso2_" + f"{T_FALTA + tc:.3f}".rstrip("0").replace(".", "_") + ".plt")
        if os.path.exists(arq):
            ta, a1, a2 = base.ler_plt(arq)
            ax.plot(ta, a1, "c--", lw=1.0, label="δ1 Anatem")
            ax.plot(ta, a2, "m--", lw=1.0, label="δ2 Anatem")
        passos = [d for d in (0.05, 0.02, 0.01, 0.005) if abs(tc / d - round(tc / d)) < 1e-6][:2]
        for dt, est in zip(passos, ("o", ".")):
            tp, p1, p2 = solucao_por_partes(tc, dt)
            cada = max(1, int(round(0.02 / dt)))     # desenha um ponto a cada 20 ms
            ax.plot(tp[::cada], p1[::cada], "k" + est, ms=4 if est == "o" else 3, label=f"δ1 por partes Δt={dt}")
            ax.plot(tp[::cada], p2[::cada], "g" + est, ms=4 if est == "o" else 3, label=f"δ2 por partes Δt={dt}")
        ax.axvline(T_FALTA, color="gray", ls=":", lw=1)
        ax.axvline(T_FALTA + tc, color="gray", ls=":", lw=1)
        ax.set_ylim(-60, 200)
        ax.set_ylabel("δ (graus)")
        ax.set_title(f"Caso 2 – solução por partes x RK4 x Anatem, tc = {tc:.3f} s")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", ncol=4, fontsize=7)
    eixos[-1].set_xlabel("tempo (s)")
    plt.tight_layout()
    plt.savefig(os.path.join(PASTA, "caso2_por_partes.png"), dpi=150)
    print("\nFigura salva em caso2_por_partes.png")
    plt.show()
