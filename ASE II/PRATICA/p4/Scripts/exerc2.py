import numpy as np
import cmath
from exerc1_direto import ler_arquivo_anarede, construir_zbarra_direta

def ret_para_polar(z):
    """Converte número complexo retangular para (módulo, ângulo em graus)."""
    mag = abs(z)
    ang = np.degrees(cmath.phase(z)) if mag > 1e-9 else 0.0
    return mag, ang

def calcular_faltas_trifasicas(Z_bus, barras, ramos, v_pre):
    """
    Recebe a matriz Z_bus e calcula I_falta, V_barras e I_ramos 
    para falta trifásica em cada barra do sistema.
    """
    idx = {b: i for i, b in enumerate(barras)}
    
    print("=" * 82)
    print(" ATIVIDADE 2: CÁLCULO DE FALTAS 3F EM CADA BARRA USANDO A MATRIZ Z_BUS")
    print("=" * 82)

    # 1. Tabela Resumo das Correntes de Falta
    print("\n1) RESUMO DAS CORRENTES DE CURTO-CIRCUITO (I_falta = V_k^0 / Z_kk):")
    print("-" * 76)
    print(f"{'Barra em Falta':^16} | {'Tensão Pré-Falta V0':^20} | {'Z_th = Z_kk (pu)':^16} | {'Corrente I_falta (pu)':^20}")
    print("-" * 76)
    
    for b_k in barras:
        k = idx[b_k]
        vk0 = v_pre[b_k]
        z_kk = Z_bus[k, k]
        i_falta = vk0 / z_kk
        
        v_m, v_a = ret_para_polar(vk0)
        i_m, i_a = ret_para_polar(i_falta)
        print(f"    Barra {b_k:<6} |  {v_m:7.4f} /_{v_a:6.2f}° pu |   {z_kk.imag:9.5f}j     |  {i_m:8.4f} /_{i_a:6.2f}° pu")
    print("-" * 76)

    # 2. Detalhamento de V_barras e I_ramos para falta em cada barra k
    print("\n2) DETALHAMENTO DE V_BARRAS E I_RAMOS PARA CADA FALTA:")
    for b_k in barras:
        k = idx[b_k]
        vk0 = v_pre[b_k]
        z_kk = Z_bus[k, k]
        
        # Corrente de falta na barra k: I_f = V_k^0 / Z_kk
        i_falta = vk0 / z_kk
        if_m, if_a = ret_para_polar(i_falta)

        print("\n" + "=" * 82)
        print(f">>> FALTA TRIFÁSICA NA BARRA {b_k}   (I_falta = {if_m:.4f} /_{if_a:.2f}° pu)")
        print("=" * 82)

        # Tensões pós-falta nas barras: V_j = V_j^0 - Z_jk * I_falta
        v_pos = {}
        print("  [Tensões Pós-Falta nas Barras (V_barras)]:")
        for b_j in barras:
            j = idx[b_j]
            vj0 = v_pre[b_j]
            v_pos[b_j] = vj0 - Z_bus[j, k] * i_falta
            vm, va = ret_para_polar(v_pos[b_j])
            fim_linha = "\n" if b_j % 3 == 0 else "   | "
            print(f"    V_{b_j} = {vm:7.4f} /_{va:6.2f}° pu", end=fim_linha)

        # Correntes nos ramos: I_ij = (V_i - V_j) / Z_b
        print("\n  [Correntes nos Ramos (I_ramos = (V_de - V_para) / Z_ramo)]:")
        print(f"    {'Ramo (De -> Para)':<24} | {'Z_ramo (pu)':<12} | {'I_ramo (pu)':<22}")
        print("    " + "-" * 64)
        for de, para, zb, desc in ramos:
            i_ij = (v_pos[de] - v_pos[para]) / zb
            im, ia = ret_para_polar(i_ij)
            print(f"    {desc:<24} |  {zb.imag:8.4f}j   | {im:8.4f} /_{ia:6.2f}° pu")


if __name__ == "__main__":
    barras, shunts, ramos = ler_arquivo_anarede("sistema6_barras.txt")
    Z_bus = np.load("zbus.npy")

    v_pre_falta = {
        1: cmath.rect(1, np.radians(0.0)),
        2: cmath.rect(1, np.radians(0.0)),
        3: cmath.rect(1, np.radians(0.0)),
        4: cmath.rect(1, np.radians(0.0)),
        5: cmath.rect(1, np.radians(0.0)),
        6: cmath.rect(1, np.radians(0.0)),
    }

    calcular_faltas_trifasicas(Z_bus, barras, ramos, v_pre_falta)