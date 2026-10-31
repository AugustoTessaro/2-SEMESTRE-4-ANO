import numpy as np

def ler_arquivo_anarede(caminho_arquivo, s_base=100.0, incluir_b_linha=True, incluir_shunt_barra=True):
    shunts = []
    ramos = []
    barras = set()

    secao = None
    with open(caminho_arquivo, 'r', encoding='utf-8') as f:
        for linha in f:
            linha_limpa = linha.strip()
            if not linha_limpa or linha_limpa.startswith('(') or linha_limpa.startswith('#'):
                continue
            
            if linha_limpa.upper().startswith('DBAR'):
                secao = 'DBAR'
                continue
            elif linha_limpa.upper().startswith('DLIN'):
                secao = 'DLIN'
                continue
            elif linha_limpa.startswith('99999') or linha_limpa.upper().startswith('FIM'):
                secao = None
                continue

            partes = linha_limpa.split()
            
            if secao == 'DBAR':
                b = int(partes[0])
                r_g = float(partes[1])
                x_g = float(partes[2])
                q_sh = float(partes[3])
                barras.add(b)
                
                if r_g != 0.0 or x_g != 0.0:
                    shunts.append((b, complex(r_g, x_g), f"Gerador Barra {b}"))
                
                if incluir_shunt_barra and q_sh != 0.0:
                    y_sh = complex(0.0, q_sh / s_base)
                    shunts.append((b, 1.0 / y_sh, f"Shunt Barra {b} ({q_sh} Mvar)"))

            elif secao == 'DLIN':
                de = int(partes[0])
                para = int(partes[1])
                circ = int(partes[2])
                r = float(partes[3])
                x = float(partes[4])
                b_tot = float(partes[5])
                barras.update([de, para])
                
                ramos.append((de, para, complex(r, x), f"Ramo {de}-{para} (Circ {circ})"))
                
                if incluir_b_linha and b_tot != 0.0:
                    z_b2 = 1.0 / complex(0.0, b_tot / 2.0)
                    shunts.append((de, z_b2, f"B/2 Linha {de}-{para} C{circ} na Barra {de}"))
                    shunts.append((para, z_b2, f"B/2 Linha {de}-{para} C{circ} na Barra {para}"))

    return sorted(list(barras)), shunts, ramos


def imprimir_matriz_passo(Z, barras, conectadas):
    """Imprime a matriz jX formatada com índices das barras."""
    header = "       " + " ".join([f"  B{b}    " for b in barras])
    print(header)
    print("      " + "-" * (len(barras) * 9))
    for i, b_lin in enumerate(barras):
        valores = " ".join([f"{Z[i, j].imag:8.5f}" for j in range(len(barras))])
        status = "*" if b_lin in conectadas else " "
        print(f"B{b_lin} {status} | {valores}")
    print(f"Barras ativas na rede (*): {sorted(list(conectadas))}\n" + "=" * 68)


def construir_zbarra_direta(barras, shunts, ramos, mostrar_passos=True):
    n = len(barras)
    idx = {b: i for i, b in enumerate(barras)}
    Z = np.zeros((n, n), dtype=complex)
    conectadas = set()

    fila = [('shunt', b, None, zb, desc) for b, zb, desc in shunts if "Gerador" in desc] + \
           [('serie', d, p, zb, desc) for d, p, zb, desc in ramos] + \
           [('shunt', b, None, zb, desc) for b, zb, desc in shunts if "Gerador" not in desc]

    passo = 1
    while fila:
        for i, (tipo, b1, b2, zb, desc) in enumerate(fila):
            if tipo == 'shunt':
                k = idx[b1]
                if b1 not in conectadas:
                    # Caso 1: Barra nova para a referência
                    Z[k, k] = zb
                    conectadas.add(b1)
                    print(f"\n>>> Passo {passo:02d} [Caso 1 - Nova p/ Ref]: {desc} | Zb = {zb.imag:+.4f}j")
                    print(f"    Ação: Z[{b1},{b1}] = {zb.imag:.5f}j")
                else:
                    # Caso 3: Barra existente para a referência (Redução de Kron)
                    z_kk = Z[k, k]
                    denom = z_kk + zb
                    col_k = Z[:, k].reshape(-1, 1)
                    lin_k = Z[k, :].reshape(1, -1)
                    Z = Z - (col_k @ lin_k) / denom
                    print(f"\n>>> Passo {passo:02d} [Caso 3 - Exist p/ Ref]: {desc} | Zb = {zb.imag:+.4f}j")
                    print(f"    Ação: Redução de Kron na Barra {b1} | Denom (Z_{b1}{b1} + Zb) = {denom.imag:.5f}j")
                
                if mostrar_passos:
                    imprimir_matriz_passo(Z, barras, conectadas)
                fila.pop(i)
                passo += 1
                break

            elif tipo == 'serie':
                c1, c2 = (b1 in conectadas), (b2 in conectadas)
                if c1 != c2:
                    # Caso 2: Barra existente para barra nova
                    b_ex, b_nv = (b1, b2) if c1 else (b2, b1)
                    k, p = idx[b_ex], idx[b_nv]
                    Z[:, p] = Z[:, k]
                    Z[p, :] = Z[k, :]
                    Z[p, p] = Z[k, k] + zb
                    conectadas.add(b_nv)
                    print(f"\n>>> Passo {passo:02d} [Caso 2 - Exist p/ Nova]: {desc} | Zb = {zb.imag:+.4f}j")
                    print(f"    Ação: Copia linha/col da Barra {b_ex} para Barra {b_nv} | Z[{b_nv},{b_nv}] = Z[{b_ex},{b_ex}] + Zb = {Z[p, p].imag:.5f}j")
                    
                    if mostrar_passos:
                        imprimir_matriz_passo(Z, barras, conectadas)
                    fila.pop(i)
                    passo += 1
                    break

                elif c1 and c2:
                    # Caso 4: Entre duas barras existentes (Redução de Kron)
                    j, k = idx[b1], idx[b2]
                    diff_col = (Z[:, j] - Z[:, k]).reshape(-1, 1)
                    diff_lin = (Z[j, :] - Z[k, :]).reshape(1, -1)
                    z_bb = Z[j, j] + Z[k, k] - 2 * Z[j, k] + zb
                    Z = Z - (diff_col @ diff_lin) / z_bb
                    print(f"\n>>> Passo {passo:02d} [Caso 4 - Entre Exist]: {desc} | Zb = {zb.imag:+.4f}j")
                    print(f"    Ação: Redução de Kron entre {b1} e {b2} | Z_bb = Z_{b1}{b1} + Z_{b2}{b2} - 2*Z_{b1}{b2} + Zb = {z_bb.imag:.5f}j")
                    
                    if mostrar_passos:
                        imprimir_matriz_passo(Z, barras, conectadas)
                    fila.pop(i)
                    passo += 1
                    break
    return Z

if __name__ == "__main__":
    print("=== CONSTRUÇÃO DIRETA DA Z_BARRA PASSO A PASSO ===")
    barras, shunts, ramos = ler_arquivo_anarede("sistema6_barras.txt", incluir_b_linha=True, incluir_shunt_barra=True)
    Z_bus = construir_zbarra_direta(barras, shunts, ramos, mostrar_passos=True)
    np.save("zbus.npy", Z_bus)