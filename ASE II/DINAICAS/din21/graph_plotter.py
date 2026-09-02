import matplotlib.pyplot as plt
import numpy as np

# Configuração do tamanho do gráfico
plt.figure(figsize=(12, 7))

# --- Curvas das Unidades Individuais ---

# Unidade 1
P1 = np.array([100, 420])
L1 = 13 + 0.02 * P1
plt.plot(P1, L1, label='Unidade 1', color='blue', marker='o')

# Unidade 2
P2 = np.array([120, 500])
L2 = 9 + 0.048 * P2
plt.plot(P2, L2, label='Unidade 2', color='red', marker='o')

# Unidade 3
P3 = np.array([50, 300])
L3 = 11.5 + 0.036 * P3
plt.plot(P3, L3, label='Unidade 3', color='green', marker='o')

# --- Curva do Sistema (Agregada) ---
# Pontos de quebra calculados:
P_sys = np.array([270, 310.56, 322.22, 953.33, 997.08, 1220])
L_sys = np.array([13.30, 14.76, 15.00, 21.40, 22.30, 33.00])
plt.plot(P_sys, L_sys, label='Sistema (Agregado)', color='purple', linewidth=2, linestyle='--', marker='s')

# --- Estilização do Gráfico ---
plt.title('Curvas de Custo Marginal - Unidades vs Sistema', fontsize=14)
plt.xlabel('Potência P (MW)', fontsize=12)
plt.ylabel('Custo Marginal $\lambda$ ($/MWh)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=11)

# Exibe o gráfico
plt.show()