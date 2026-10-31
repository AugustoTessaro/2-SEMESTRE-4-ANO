import os
import numpy as np
import cmath
from exerc1_direto import ler_arquivo_anarede, construir_zbarra_direta
from exerc2 import ret_para_polar, calcular_faltas_trifasicas

def adicionar_carga_caso3(Z_atual, idx_barra, z_carga, nome_barra):
    """
    Aplica o Caso 3 (Redução de Kron) para conectar a impedância de carga 
    entre uma barra existente k e a referência (terra).
    """
    Z_nova = Z_atual.copy()
    k = idx_barra
    z_kk = Z_nova[k, k]
    col_k = Z_nova[:, k].reshape(-1, 1)
    lin_k = Z_nova[k, :].reshape(1, -1)
    
    Z_nova = Z_nova - (col_k @ lin_k) / (z_kk + z_carga)
    print(f" -> [Caso 3] Carga adicionada na Barra {nome_barra} | Z_carga = {z_carga.real:.5f} + {z_carga.imag:.5f}j pu")
    return Z_nova


def comparar_atividades_2_e_3(Z_sem_carga, Z_com_carga, barras, v_pre_2, v_pre_3):
    """Imprime uma tabela lado a lado comparando Z_kk e I_falta das Atividades 2 e 3."""
    idx = {b: i for i, b in enumerate(barras)}
    
    print("\n" + "=" * 95)
    print(" COMPARATIVO DIRETO: ATIVIDADE 2 (SEM CARGA) vs ATIVIDADE 3 (COM CARGAS NAS BARRAS 2 E 6)")
    print("=" * 95)
    print(f"{'Barra':^7} | {'|Z_kk| Ativ 2':^14} | {'|Z_kk| Ativ 3':^14} | {'I_falta Ativ 2 (pu)':^21} | {'I_falta Ativ 3 (pu)':^21} | {'Variação |I_f|'}")
    print("-" * 95)

    for b_k in barras:
        k = idx[b_k]
        # Atividade 2
        z2 = Z_sem_carga[k, k]
        if2 = v_pre_2[b_k] / z2
        if2_m, if2_a = ret_para_polar(if2)

        # Atividade 3
        z3 = Z_com_carga[k, k]
        if3 = v_pre_3[b_k] / z3
        if3_m, if3_a = ret_para_polar(if3)

        diff_perc = ((if3_m - if2_m) / if2_m) * 100.0

        print(f"  B{b_k:<4} |  {abs(z2):10.5f} pu |  {abs(z3):10.5f} pu | {if2_m:8.4f} /_{if2_a:6.2f}° | {if3_m:8.4f} /_{if3_a:6.2f}° |  {diff_perc:+6.2f} %")
    print("-" * 95)


if __name__ == "__main__":
    S_BASE = 100.0  # MVA
    barras, shunts, ramos = ler_arquivo_anarede("sistema6_barras.txt")
    idx_map = {b: i for i, b in enumerate(barras)}

    # 1. Carrega a Z_bus original (sem cargas) do Exercício 1
    if os.path.exists("zbus.npy"):
        Z_bus_sem_carga = np.load("zbus.npy")
    else:
        Z_bus_sem_carga = construir_zbarra_direta(barras, shunts, ramos, mostrar_passos=False)
        np.save("zbus.npy", Z_bus_sem_carga)

    print("=== EXERCÍCIO 3: INCLUSÃO DE CARGAS (200 MW e 75 Mvar nas Barras 2 e 6) ===")
    
    # 2. Calcula a impedância equivalente da carga: Z_L = |V|^2 / (P - jQ)
    P_pu = 200.0 / S_BASE  # 2.00 pu
    Q_pu = 75.0 / S_BASE   # 0.75 pu
    Z_carga = 1.0 / complex(P_pu, -Q_pu)  # Considerando V_nom = 1.0 pu

    # 3. Atualiza a matriz Z_bus adicionando as cargas nas Barras 2 e 6 (Caso 3)
    Z_bus_com_carga = Z_bus_sem_carga.copy()
    Z_bus_com_carga = adicionar_carga_caso3(Z_bus_com_carga, idx_map[2], Z_carga, nome_barra=2)
    Z_bus_com_carga = adicionar_carga_caso3(Z_bus_com_carga, idx_map[6], Z_carga, nome_barra=6)

    print("\nNova Matriz Z_barra com Cargas (Módulo |Z| em pu):")
    np.set_printoptions(precision=5, suppress=True, linewidth=120)
    print(np.abs(Z_bus_com_carga))

    # 4. Definição das tensões pré-falta
    # Atividade 2: 1.0 < 0° pu
    v_pre_ativ2 = {b: complex(1.0, 0.0) for b in barras}

    # Atividade 3: Se usar V_pre = 1.0 < 0° pu (ou substitua pelos valores do fluxo com carga do PW se preferir)
    v_pre_ativ3 = {
            1: cmath.rect(1.020, np.radians(0.0)),
            2: cmath.rect(0.935, np.radians(0.3582)),
            3: cmath.rect(1.020, np.radians(-0.0112)),
            4: cmath.rect(1.008, np.radians(0.4922)),
            5: cmath.rect(0.998, np.radians(0.4843)),
            6: cmath.rect(0.943, np.radians(0.3186)),
        }

    # 5. Executa os cálculos detalhados da Atividade 3 (I_falta, V_barras e I_ramos)
    calcular_faltas_trifasicas(Z_bus_com_carga, barras, ramos, v_pre_ativ3)

    # 6. Imprime a tabela comparativa entre Atividade 2 e Atividade 3
    comparar_atividades_2_e_3(Z_bus_sem_carga, Z_bus_com_carga, barras, v_pre_ativ2, v_pre_ativ3)

