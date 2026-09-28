import numpy as np
import matplotlib.pyplot as plt

# Parâmetros de simulação
fc = 1000      # Frequência da portadora (Hz)
fm = 100       # Frequência da mensagem (Hz)
Ac = 1.0       # Amplitude da portadora
ka = 0.8       # Sensibilidade de amplitude
Am = 1.0       # Amplitude da mensagem
fs = 100000    # Frequência de amostragem (Hz)
t = np.linspace(0, 0.005, int(fs * 0.005))

# Sinais
m = Am * np.cos(2 * np.pi * fm * t)
c = Ac * np.cos(2 * np.pi * fc * t)
s = Ac * (1 + ka * m) * np.cos(2 * np.pi * fc * t)
envelope = Ac * (1 + ka * m)

# Gráficos
plt.figure(figsize=(10, 6))

plt.subplot(3, 1, 1)
plt.plot(t * 1e3, m, 'g')
plt.title('Mensagem m(t)')
plt.ylabel('Amplitude')

plt.subplot(3, 1, 2)
plt.plot(t * 1e3, c, 'b')
plt.title('Portadora c(t)')
plt.ylabel('Amplitude')

plt.subplot(3, 1, 3)
plt.plot(t * 1e3, s, 'r', label='Sinal AM s(t)')
plt.plot(t * 1e3, envelope, 'k--', label='Envelope')
plt.plot(t * 1e3, -envelope, 'k--')
plt.title('Sinal Modulado AM s(t) com Envelope Sobreposto')
plt.xlabel('Tempo (ms)')
plt.ylabel('Amplitude')
plt.legend()

plt.tight_layout()
plt.show()