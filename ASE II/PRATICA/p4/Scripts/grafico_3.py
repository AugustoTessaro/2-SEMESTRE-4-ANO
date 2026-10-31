import numpy as np
import matplotlib.pyplot as plt

# Dados extraídos exatamente do seu resultado_3.txt
barras = ['Barra 1', 'Barra 2\n(Carga)', 'Barra 3', 'Barra 4', 'Barra 5', 'Barra 6\n(Carga)']
x = np.arange(len(barras))

# 1. Impedâncias |Z_kk| (pu)
z_ativ2 = np.array([0.02811, 0.06470, 0.01831, 0.02376, 0.02977, 0.08377])
z_ativ3 = np.array([0.02751, 0.05940, 0.01791, 0.02295, 0.02855, 0.07604])

# 2. Tensões Pré-Falta |V_0| (pu)
v0_ativ2 = np.array([1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000])
v0_ativ3 = np.array([1.0200, 0.9350, 1.0200, 1.0080, 0.9980, 0.9430])

# 3. Módulo da Corrente de Falta |I_falta| (pu)
if_ativ2 = np.array([35.5709, 15.4561, 54.6145, 42.0907, 33.5950, 11.9376])
if_ativ3 = np.array([37.0753, 15.7398, 56.9368, 43.9262, 34.9583, 12.4015])

# 4. Ângulo da Corrente de Falta (graus)
ang_ativ2 = np.array([-90.00, -90.00, -90.00, -90.00, -90.00, -90.00])
ang_ativ3 = np.array([-87.97, -80.11, -87.99, -86.21, -85.43, -78.81])

# Variações percentuais da Ativ 2 para a Ativ 3
var_z  = ((z_ativ3 - z_ativ2) / z_ativ2) * 100
var_v0 = ((v0_ativ3 - v0_ativ2) / v0_ativ2) * 100
var_if = ((if_ativ3 - if_ativ2) / if_ativ2) * 100

# Configuração da Figura (3 painéis lado a lado em formato widescreen 16:9)
plt.rcParams['font.sans-serif'] = 'Segoe UI'
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.8))
fig.suptitle('Atividade 3: Impacto das Cargas (200 MW + 75 Mvar nas Barras 2 e 6) no Curto-Circuito', 
             fontsize=15, fontweight='bold', y=0.98)

largura = 0.36

# ==============================================================================
# PAINEL 1: Módulo da Corrente de Falta |I_falta| (Ativ 2 vs Ativ 3)
# ==============================================================================
b1 = ax1.bar(x - largura/2, if_ativ2, largura, label='Ativ. 2 (Sem Carga / V0=1.0)', color='#7f8c8d', alpha=0.85)
b2 = ax1.bar(x + largura/2, if_ativ3, largura, label='Ativ. 3 (Com Cargas B2 e B6)', color='#2980b9')

ax1.set_title('1) Módulo da Corrente de Falta |I_falta|', fontweight='bold', fontsize=12)
ax1.set_ylabel('Corrente de Falta (pu)', fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(barras, fontweight='bold')
ax1.set_ylim(0, 66)
ax1.grid(axis='y', linestyle='--', alpha=0.5)
ax1.legend(loc='upper right', fontsize=9.5)

# Rótulos com a variação percentual no topo de cada par de barras
for i in range(len(barras)):
    altura_max = max(if_ativ2[i], if_ativ3[i])
    ax1.annotate(f'+{var_if[i]:.2f}%',
                 xy=(x[i] + largura/2, altura_max),
                 xytext=(0, 5), textcoords="offset points",
                 ha='center', va='bottom', fontsize=9.5, fontweight='bold', color='#1b4f72')

# ==============================================================================
# PAINEL 2: Causa e Efeito (Variação % de V0, |Z_kk| e |I_falta|)
# ==============================================================================
w = 0.25
ax2.bar(x - w, var_v0, w, label='Variação Tensão Pré-Falta (ΔV0 %)', color='#e67e22')
ax2.bar(x,     var_z,  w, label='Variação Impedância (Δ|Z_kk| %)', color='#c0392b')
ax2.bar(x + w, var_if, w, label='Variação Resultante (Δ|I_falta| %)', color='#27ae60')

ax2.axhline(0, color='black', linewidth=1)
ax2.set_title('2) Causa e Efeito: Por que B2 e B6 subiram menos?', fontweight='bold', fontsize=12)
ax2.set_ylabel('Variação Percentual (%)', fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(barras, fontweight='bold')
ax2.set_ylim(-11.5, 7.5)
ax2.grid(axis='y', linestyle='--', alpha=0.5)
ax2.legend(loc='lower left', fontsize=9)

# Adiciona os valores nas barras de variação de B2 e B6 para destaque
for i in [1, 5]:  # Índices da Barra 2 e Barra 6
    ax2.text(x[i] - w, var_v0[i] - 1.2, f'{var_v0[i]:.1f}%', ha='center', fontsize=8.5, fontweight='bold', color='#d35400')
    ax2.text(x[i],     var_z[i] - 1.2,  f'{var_z[i]:.1f}%',  ha='center', fontsize=8.5, fontweight='bold', color='#922b21')
    ax2.text(x[i] + w, var_if[i] + 0.4, f'+{var_if[i]:.1f}%', ha='center', fontsize=8.5, fontweight='bold', color='#1e8449')

# ==============================================================================
# PAINEL 3: Deslocamento do Ângulo da Corrente de Falta (Efeito Resistivo 200 MW)
# ==============================================================================
ax3.plot(x, ang_ativ2, marker='o', linestyle='--', color='#7f8c8d', linewidth=2, markersize=7, label='Ativ. 2 (-90° Puramente Indutivo)')
ax3.plot(x, ang_ativ3, marker='s', linestyle='-',  color='#8e44ad', linewidth=2.5, markersize=8, label='Ativ. 3 (Com Carga Resistiva-Indutiva)')

ax3.fill_between(x, ang_ativ2, ang_ativ3, color='#8e44ad', alpha=0.12)
ax3.set_title('3) Ângulo da Corrente de Falta (Efeito dos 200 MW)', fontweight='bold', fontsize=12)
ax3.set_ylabel('Ângulo de I_falta (graus)', fontweight='bold')
ax3.set_xticks(x)
ax3.set_xticklabels(barras, fontweight='bold')
ax3.set_ylim(-92, -75)
ax3.grid(True, linestyle='--', alpha=0.5)
ax3.legend(loc='upper left', fontsize=9.5)

for i in range(len(barras)):
    ax3.annotate(f'{ang_ativ3[i]:.1f}°',
                 xy=(x[i], ang_ativ3[i]),
                 xytext=(0, 8), textcoords="offset points",
                 ha='center', fontsize=9.5, fontweight='bold', color='#6c3483')

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('comparativo_exerc3.png', dpi=300)
print(">>> Gráfico salvo com sucesso como 'comparativo_exerc3.png'!")
plt.show()