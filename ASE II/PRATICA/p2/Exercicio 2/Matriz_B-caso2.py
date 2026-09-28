import numpy as np
from pypower.api import runpf, ppoption
from pypower.makeYbus import makeYbus
from pypower.ext2int import ext2int
from numpy.linalg import inv

def sistema_4barras():
    ppc = {"version": "2", "baseMVA": 100.0}
    ppc["bus"] = np.array([
        [1, 3,   0.00,    0.00, 0.0, 0.0, 1, 1.0000, 0.0, 230.0, 1, 1.1, 0.9],
        [2, 2,   0.00,    0.00, 0.0, 0.0, 1, 1.0000, 0.0, 230.0, 1, 1.1, 0.9],
        [3, 2, 220.00,  136.34, 0.0, 0.0, 1, 0.9605, 0.0, 230.0, 1, 1.1, 0.9],
        [4, 1, 280.00,  173.52, 0.0, 0.0, 1, 1.0000, 0.0, 230.0, 1, 1.1, 0.9]
    ])
    ppc["gen"] = np.array([
        [1,   0.0, 0.0, 999.0, -999.0, 1.0000, 100.0, 1, 450.0, 200.0],
        [2, 166.0, 0.0, 999.0, -999.0, 1.0000, 100.0, 1, 350.0, 150.0],
        [3, 111.0, 0.0, 999.0, -999.0, 0.9605, 100.0, 1, 235.0, 100.0]
    ])
    ppc["branch"] = np.array([
        [1, 4, 0.00744, 0.0372, 0.0775, 500, 500, 500, 0, 0, 1, -360, 360],
        [1, 3, 0.01008, 0.0504, 0.1025, 500, 500, 500, 0, 0, 1, -360, 360],
        [2, 3, 0.00744, 0.0372, 0.0775, 500, 500, 500, 0, 0, 1, -360, 360],
        [2, 4, 0.01272, 0.0636, 0.1275, 500, 500, 500, 0, 0, 1, -360, 360]
    ])
    return ppc

# =============================================================================
# FUNÇÃO PARA CALCULAR A MATRIZ B (Baseada no seu código original)
# =============================================================================
def calcular_matriz_B(ppc, resultados_pf):
    baseMVA = ppc['baseMVA']
    barras = resultados_pf['bus']
    geradores = resultados_pf['gen']
    num_barras = len(barras)
    num_gens = len(geradores)

    V = np.zeros(num_barras, dtype=complex)
    Pload, Qload = np.zeros(num_barras), np.zeros(num_barras)
    Pgen, Qgen = np.zeros(num_barras), np.zeros(num_barras)
    gen_buses = []

    for i in range(num_barras):
        V[i] = barras[i, 7] * np.exp(1j * np.deg2rad(barras[i, 8]))
        Pload[i] = barras[i, 2] / baseMVA
        Qload[i] = barras[i, 3] / baseMVA

    for i in range(num_gens):
        bus_idx = np.where(barras[:, 0] == geradores[i, 0])[0][0]
        gen_buses.append(bus_idx)
        Pgen[bus_idx] = geradores[i, 1] / baseMVA
        Qgen[bus_idx] = geradores[i, 2] / baseMVA

    ref_bus = 0 # Barra Slack

    ppc_int = ext2int(ppc)
    Ybus_esparsa, _, _ = makeYbus(ppc_int["baseMVA"], ppc_int["bus"], ppc_int["branch"])
    Zbus = inv(Ybus_esparsa.toarray())
    Rbus = Zbus.real

    Sload = Pload + 1j * Qload
    Iload = np.zeros_like(V, dtype=complex)
    mask_load = (np.abs(Sload) > 0)
    Iload[mask_load] = -np.conj(Sload[mask_load]) / np.conj(V[mask_load])

    ID = np.sum(Iload)
    d = np.zeros_like(V, dtype=complex)
    if np.abs(ID) > 0:
        d[mask_load] = Iload[mask_load] / ID
    load_buses = [k for k in range(num_barras) if mask_load[k]]

    den = sum(d[l] * Zbus[ref_bus, l] for l in load_buses)
    t = {g: Zbus[ref_bus, g] / den for g in gen_buses}
    t_ref = Zbus[ref_bus, ref_bus] / den

    K = len(gen_buses)
    C = np.zeros((num_barras, K + 1), dtype=complex)
    for m in range(num_barras):
        dm = d[m]
        for col, g in enumerate(gen_buses):
            C[m, col] += -dm * t[g]
        C[m, K] += -dm * t_ref

    for col, g in enumerate(gen_buses):
        C[g, col] += 1.0

    CRC = C.conj().T @ Rbus @ C

    alpha = np.zeros(K, dtype=complex)
    for i, g in enumerate(gen_buses):
        s = Qgen[g] / Pgen[g]
        alpha[i] = (1 - 1j * s) / np.conj(V[g])
    In0 = -V[ref_bus] / Zbus[ref_bus, ref_bus]

    Tdiag = np.concatenate([alpha, [In0]])
    Talpha = np.diag(Tdiag)
    Ts = Talpha.conj().T @ CRC @ Talpha

    Matriz_B = ((Ts + Ts.conj().T) / 2).real
    return Matriz_B, Pgen, Pload

# =============================================================================
# EXECUÇÃO DOS CASOS
# =============================================================================
np.set_printoptions(precision=6, suppress=True)
opcoes = ppoption(PF_ALG=1, OUT_ALL=0)

# --- CASO BASE ---
ppc_base = sistema_4barras()
res_base, _ = runpf(ppc_base, opcoes)
B_base, Pg_base, Pload_base = calcular_matriz_B(ppc_base, res_base)

print("="*50)
print("CASO BASE")
print("="*50)
print("Matriz B Base:\n", B_base)


# --- CASO 1: AUMENTO DE 5% (+5%) ---
ppc_mais = sistema_4barras()
ppc_mais["gen"][1, 1] *= 1.05
ppc_mais["gen"][2, 1] *= 1.05
res_mais, _ = runpf(ppc_mais, opcoes)
B_mais, Pg_mais, _ = calcular_matriz_B(ppc_mais, res_mais)

print("\n" + "="*50)
print("CASO 1: AUMENTO DE 5% (Barras 2 e 3)")
print("="*50)
print("Matriz B Exata para o Caso 1 (Veja como mudou!):\n", B_mais)

# Comparação do erro (Usando a B_Base para estimar, igual se faz em despacho)
pbar_mais = np.append(Pg_mais[[0,1,2]], 1.0) # Pegando apenas as barras geradoras + termo constante
Perdas_Est_Mais = float(pbar_mais @ B_base @ pbar_mais)
Perdas_Reais_Mais = np.sum(Pg_mais) - np.sum(Pload_base)
print(f"\nPerdas Reais (Fluxo):   {Perdas_Reais_Mais * 100:.4f} MW")
print(f"Perdas Estimadas (Kron): {Perdas_Est_Mais * 100:.4f} MW")


# --- CASO 2: REDUÇÃO DE 5% (-5%) ---
ppc_menos = sistema_4barras()
ppc_menos["gen"][1, 1] *= 0.95
ppc_menos["gen"][2, 1] *= 0.95
res_menos, _ = runpf(ppc_menos, opcoes)
B_menos, Pg_menos, _ = calcular_matriz_B(ppc_menos, res_menos)

print("\n" + "="*50)
print("CASO 2: REDUÇÃO DE 5% (Barras 2 e 3)")
print("="*50)
print("Matriz B Exata para o Caso 2 (Diferente de novo!):\n", B_menos)

# Comparação do erro (Usando a B_Base para estimar)
pbar_menos = np.append(Pg_menos[[0,1,2]], 1.0)
Perdas_Est_Menos = float(pbar_menos @ B_base @ pbar_menos)
Perdas_Reais_Menos = np.sum(Pg_menos) - np.sum(Pload_base)
print(f"\nPerdas Reais (Fluxo):   {Perdas_Reais_Menos * 100:.4f} MW")
print(f"Perdas Estimadas (Kron): {Perdas_Est_Menos * 100:.4f} MW")


# --- CASO EXTRA: AUMENTO DE 100% NA CARGA (Quebra real da aproximação) ---
ppc_carga = sistema_4barras()
# Dobrando a carga (Ativa e Reativa)
ppc_carga["bus"][2, 2] *= 2.0  # Pload barra 3
ppc_carga["bus"][2, 3] *= 2.0  # Qload barra 3
ppc_carga["bus"][3, 2] *= 2.0  # Pload barra 4
ppc_carga["bus"][3, 3] *= 2.0  # Qload barra 4

res_carga, _ = runpf(ppc_carga, opcoes)
_, Pg_carga, Pload_carga = calcular_matriz_B(ppc_carga, res_carga)

print("\n" + "="*50)
print("CASO EXTRA: AUMENTO DE 100% NA CARGA")
print("="*50)

# 1. Perdas pelo FLUXO DE CARGA (Novo estado estressado)
Perdas_Reais_Carga = np.sum(Pg_carga) - np.sum(Pload_carga)

# 2. Perdas pelo KRON (Tentando usar a B_Base do sistema original leve)
pbar_carga = np.append(Pg_carga[[0,1,2]], 1.0)
Perdas_Est_Carga = float(pbar_carga @ B_base @ pbar_carga)

erro_carga_perc = abs(Perdas_Est_Carga - Perdas_Reais_Carga) / Perdas_Reais_Carga * 100

print(f"Perdas Reais (Fluxo):    {Perdas_Reais_Carga * 100:.4f} MW")
print(f"Perdas Estimadas (Kron): {Perdas_Est_Carga * 100:.4f} MW")
print(f"Erro da aproximação:     {erro_carga_perc:.4f} %")
