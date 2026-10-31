import numpy as np
import cmath
from exerc1_direto import ler_arquivo_anarede, construir_zbarra_direta
from exerc2 import ret_para_polar

def preparar_sistema_falta_meio_linha(caminho_arquivo, barra_ficticia=7):
    """
    Lê o arquivo original e divide o Circuito 2 da linha 6-2 ao meio (50%),
    inserindo a barra fictícia 7 para representar o ponto de falta.
    """
    barras, shunts_orig, ramos_orig = ler_arquivo_anarede(
        caminho_arquivo, incluir_b_linha=False, incluir_shunt_barra=True
    )
    
    # Adiciona a barra fictícia 7 na lista de barras
    barras.append(barra_ficticia)
    barras = sorted(barras)

    novos_ramos = []
    novos_shunts = list(shunts_orig)

    # Vamos reconstruir os ramos e shunts de linha dividindo o Circuito 2 de 6-2
    with open(caminho_arquivo, 'r', encoding='utf-8') as f:
        secao = None
        for linha in f:
            l = linha.strip()
            if not l or l.startswith('(') or l.startswith('#'): continue
            if l.upper().startswith('DLIN'): secao = 'DLIN'; continue
            if l.startswith('99999') or l.upper().startswith('FIM'): secao = None; continue

            if secao == 'DLIN':
                p = l.split()
                de, para, circ = int(p[0]), int(p[1]), int(p[2])
                r, x, b_tot = float(p[3]), float(p[4]), float(p[5])
                z_serie = complex(r, x)

                # Verifica se é o Circuito 2 entre as barras 6 e 2 (ou 2 e 6)
                if ({de, para} == {2, 6}) and circ == 2:
                    # Divide impedância série e susceptância B ao meio (50%)
                    z_meio = z_serie / 2.0
                    b_meio = b_tot / 2.0

                    novos_ramos.append((de, barra_ficticia, z_meio, f"Meia-Linha {de}-{barra_ficticia} (50% C2)"))
                    novos_ramos.append((barra_ficticia, para, z_meio, f"Meia-Linha {barra_ficticia}-{para} (50% C2)"))

                    if b_meio != 0.0:
                        # Modelo Pi de cada metade: b_meio / 2 em cada extremidade do sub-trecho
                        z_b2_meio = 1.0 / complex(0.0, b_meio / 2.0)
                        novos_shunts.append((de, z_b2_meio, f"B/2 Sub-trecho {de}-{barra_ficticia} na Barra {de}"))
                        novos_shunts.append((barra_ficticia, z_b2_meio, f"B/2 Sub-trecho {de}-{barra_ficticia} na Barra {barra_ficticia}"))
                        novos_shunts.append((barra_ficticia, z_b2_meio, f"B/2 Sub-trecho {barra_ficticia}-{para} na Barra {barra_ficticia}"))
                        novos_shunts.append((para, z_b2_meio, f"B/2 Sub-trecho {barra_ficticia}-{para} na Barra {para}"))
                else:
                    # Ramos normais intactos
                    novos_ramos.append((de, para, z_serie, f"Ramo {de}-{para} (Circ {circ})"))
                    if b_tot != 0.0:
                        z_b2 = 1.0 / complex(0.0, b_tot / 2.0)
                        novos_shunts.append((de, z_b2, f"B/2 Linha {de}-{para} C{circ} na Barra {de}"))
                        novos_shunts.append((para, z_b2, f"B/2 Linha {de}-{para} C{circ} na Barra {para}"))

    return barras, novos_shunts, novos_ramos


def calcular_falta_barra_ficticia(Z_bus, barras, ramos, barra_falta=7):
    idx = {b: i for i, b in enumerate(barras)}
    k = idx[barra_falta]

    # Tensão pré-falta plana 1.0 < 0° pu (conforme orientação do professor para rede sem carga)
    v_pre = {b: complex(1.0, 0.0) for b in barras}
    
    vk0 = v_pre[barra_falta]
    z_kk = Z_bus[k, k]
    i_falta = vk0 / z_kk
    if_m, if_a = ret_para_polar(i_falta)

    print("\n" + "=" * 84)
    print(f" ATIVIDADE 5: FALTA 3F NO MEIO DA LINHA 6-2 (BARRA FICTÍCIA {barra_falta})")
    print("=" * 84)
    print(f" -> Impedância de Thévenin no meio da linha (Z_{barra_falta}{barra_falta}): {z_kk.imag:.5f}j pu")
    print(f" -> Corrente de Falta no meio da linha (I_falta)     : {if_m:.4f} /_{if_a:.2f}° pu")
    print("-" * 84)

    # 1. Tensões Pós-Falta nas Barras
    v_pos = {}
    print("\n1) TENSÕES PÓS-FALTA NAS BARRAS DO SISTEMA (V_barras):")
    for b_j in barras:
        j = idx[b_j]
        v_pos[b_j] = v_pre[b_j] - Z_bus[j, k] * i_falta
        vm, va = ret_para_polar(v_pos[b_j])
        rotulo = f"Barra {b_j} (Meio da Linha)" if b_j == barra_falta else f"Barra {b_j}"
        print(f"   {rotulo:<24}: V_{b_j} = {vm:7.4f} /_{va:6.2f}° pu")

    # 2. Correntes nos Ramos
    print("\n2) CORRENTES NOS RAMOS DO SISTEMA (I_ramos = (V_de - V_para) / Z_ramo):")
    print(f"   {'Ramo / Trecho':<32} | {'Z_ramo (pu)':<12} | {'I_ramo (pu)':<22}")
    print("   " + "-" * 72)
    for de, para, zb, desc in ramos:
        i_ij = (v_pos[de] - v_pos[para]) / zb
        im, ia = ret_para_polar(i_ij)
        print(f"   {desc:<32} |  {zb.imag:8.4f}j   | {im:8.4f} /_{ia:6.2f}° pu")


if __name__ == "__main__":
    print("=== CONSTRUINDO Z_BUS (7x7) COM A BARRA FICTÍCIA 7 NO MEIO DA LINHA 6-2 ===")
    barras, shunts, ramos = preparar_sistema_falta_meio_linha("sistema6_barras.txt", barra_ficticia=7)
    
    # Constrói a nova Z_bus 7x7 pelo método direto
    Z_bus_7x7 = construir_zbarra_direta(barras, shunts, ramos, mostrar_passos=True)

    # Calcula a falta na Barra 7 (Meio da linha 6-2)
    calcular_falta_barra_ficticia(Z_bus_7x7, barras, ramos, barra_falta=7)