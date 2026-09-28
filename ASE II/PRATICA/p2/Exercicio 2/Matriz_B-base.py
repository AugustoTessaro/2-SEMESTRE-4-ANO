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

# Executar Fluxo de Potência
ppc = sistema_4barras()
opcoes = ppoption(PF_ALG=1, OUT_ALL=0)
resultados, _ = runpf(ppc, opcoes)

baseMVA = ppc['baseMVA']
barras = resultados['bus']
geradores = resultados['gen']
num_barras = len(barras)
num_gens = len(geradores)

# 1. OBTER TENSÕES E POTÊNCIAS DO SISTEMA
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

ref_bus = 0 # Barra Slack (Índice 0)

print(Pgen)
print(Qgen)
print(V)


# 2. MATRIZ DE RESISTÊNCIAS NODAIS (R_BARRA)
ppc_int = ext2int(ppc)
Ybus_esparsa, _, _ = makeYbus(ppc_int["baseMVA"], ppc_int["bus"], ppc_int["branch"])
Zbus = inv(Ybus_esparsa.toarray())
Rbus = Zbus.real

# 3. FATORES DE DISTRIBUIÇÃO E TRANSLAÇÃO (Equações de matrizB.py)
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

# 4. CONSTRUIR MATRIZ DE TRANSFORMAÇÃO (C)
K = len(gen_buses)
C = np.zeros((num_barras, K + 1), dtype=complex)
for m in range(num_barras):
    dm = d[m]
    for col, g in enumerate(gen_buses):
        C[m, col] += -dm * t[g]
    C[m, K] += -dm * t_ref

for col, g in enumerate(gen_buses):
    C[g, col] += 1.0

# 5. MATRIZ HERMITIANA EXPANDIDA (T_alpha)
CRC = C.conj().T @ Rbus @ C

alpha = np.zeros(K, dtype=complex)
for i, g in enumerate(gen_buses):
    s = Qgen[g] / Pgen[g]
    alpha[i] = (1 - 1j * s) / np.conj(V[g])
In0 = -V[ref_bus] / Zbus[ref_bus, ref_bus]

Tdiag = np.concatenate([alpha, [In0]])
Talpha = np.diag(Tdiag)
Ts = Talpha.conj().T @ CRC @ Talpha

# 6. EXTRAÇÃO DOS COEFICIENTES
Matriz_B_Total = ((Ts + Ts.conj().T) / 2).real

B_ij = Matriz_B_Total[0:num_gens, 0:num_gens]
B_i0 = 2 * Matriz_B_Total[0:num_gens, -1]
B_00 = Matriz_B_Total[-1, -1]

# Validação
Pg_vec = np.array([Pgen[g] for g in gen_buses])
pbar = np.append(Pg_vec, 1.0)
Perdas_Polinomio = float(pbar @ Matriz_B_Total @ pbar)
Perdas_Fluxo = np.sum(Pg_vec) - np.sum(Pload)

np.set_printoptions(precision=8, suppress=True)
print("\nMatriz B (forma aumentada [Pg...,1]):\n", Matriz_B_Total)
print("\n--- Coeficientes Quadráticos (B_ij) ---")
print(B_ij)
print("\n--- Coeficientes Lineares (B_i0) ---")
print(B_i0)
print(f"\n--- Termo Constante (B_00) ---")
print(f"{B_00:.8f}")

print(f"\nPerdas pelo Polinômio de Kron: {Perdas_Polinomio * 100:.4f} MW")
print(f"Perdas exatas do Fluxo de Carga: {Perdas_Fluxo * 100:.4f} MW")