import numpy as np
import matplotlib.pyplot as plt

# 1. Dados das curvas individuais (apenas dentro dos limites operacionais)
P1 = np.linspace(200, 450, 100)
lambda1 = 0.008 * P1 + 5.2

P2 = np.linspace(150, 350, 100)
lambda2 = 0.012 * P2 + 5.0

P3 = np.linspace(100, 235, 100)
lambda3 = 0.016 * P3 + 5.0

# 2. Curva agregada do sistema
# O lambda válido do sistema vai de 6.6 (mínimo do G3) até 9.2 (máximo do G2)
lambdas_sys = np.linspace(6.6, 9.2, 500)

# Calculamos as potências individuais para o sistema (com clip para travar nos limites)
P1_sys = np.clip((lambdas_sys - 5.2) / 0.008, 200, 450)
P2_sys = np.clip((lambdas_sys - 5.0) / 0.012, 150, 350)
P3_sys = np.clip((lambdas_sys - 5.0) / 0.016, 100, 235)

# Somamos as potências para criar a curva preta
P_Total = P1_sys + P2_sys + P3_sys

# 3. Configuração do Gráfico
plt.figure(figsize=(12, 7))

# Plotando os geradores (linhas sólidas flutuando apenas na área de operação)
plt.plot(P1, lambda1, label='G1 ($b=5.2, a=0.008$)', color='blue', linewidth=2)
plt.plot(P2, lambda2, label='G2 ($b=5.0, a=0.012$)', color='orange', linewidth=2)
plt.plot(P3, lambda3, label='G3 ($b=5.0, a=0.016$)', color='green', linewidth=2)

# Plotando a curva do sistema
plt.plot(P_Total, lambdas_sys, label='Sistema (G1 + G2 + G3)', color='black', linewidth=2.5, linestyle='--')

# Linhas de referência do Caso 4
plt.axvline(x=1020, color='red', linestyle=':', alpha=0.6, label='Demanda Caso 4 (1020 MW)')
plt.axhline(y=9.02, color='red', linestyle=':', alpha=0.6, label='$\lambda$ Redespacho = 9.02')

# Formatação visual
plt.title('Despacho Econômico: Curvas Estritas aos Limites Operacionais', fontsize=14)
plt.xlabel('Potência / Demanda (MW)', fontsize=12)
plt.ylabel('Custo Marginal $\lambda$ ($/MWh)', fontsize=12)
plt.legend(loc='lower right', fontsize=11)
plt.grid(True, linestyle=':', alpha=0.7)

# Ajuste visual do eixo X
plt.xlim(0, 1100)

plt.show()