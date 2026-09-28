import numpy as np
import matplotlib.pyplot as plt

# Variação do índice de modulação
mu_vals = np.linspace(0, 1, 100)
eta_vals = (mu_vals**2) / (2 + mu_vals**2) * 100

# Pontos específicos
mu_pts = np.array([0.2, 0.5, 1.0])
eta_pts = (mu_pts**2) / (2 + mu_pts**2) * 100

plt.figure(figsize=(7, 4))
plt.plot(mu_vals, eta_vals, 'b-', linewidth=2, label='Eficiência \u03b7(\u03bc)')
plt.scatter(mu_pts, eta_pts, color='red', zorder=5)

for mu_i, eta_i in zip(mu_pts, eta_pts):
    plt.annotate(f'\u03bc={mu_i}: {eta_i:.2f}%', (mu_i, eta_i), 
                 textcoords="offset points", xytext=(-20,10), ha='center')

plt.title('Eficiência da Modulação AM vs. Índice de Modulação (\u03bc)')
plt.xlabel('Índice de Modulação (\u03bc)')
plt.ylabel('Eficiência \u03b7 (%)')
plt.grid(True)
plt.legend()
plt.show()