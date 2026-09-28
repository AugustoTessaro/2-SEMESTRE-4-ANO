import numpy as np
import matplotlib.pyplot as plt

# Parâmetros
fc = 1000      # Frequência da portadora (Hz)
fm = 100       # Frequência da mensagem (Hz)
Ac = 1.0
mu = 0.8
fs = 20000
N = 8192
t = np.arange(N) / fs

# Sinal AM
m = np.cos(2 * np.pi * fm * t)
s = Ac * (1 + mu * m) * np.cos(2 * np.pi * fc * t)

# FFT
S = np.fft.fft(s) / N
freqs = np.fft.fftfreq(N, 1/fs)

# Apenas frequências positivas
pos_mask = freqs >= 0
freqs = freqs[pos_mask]
S_mag = np.abs(S[pos_mask]) * 2

plt.figure(figsize=(8, 4))
plt.plot(freqs, S_mag, 'b')
plt.xlim(fc - 3*fm, fc + 3*fm)
plt.title('Espetro de Frequência do Sinal AM (FFT)')
plt.xlabel('Frequência (Hz)')
plt.ylabel('Magnitude')
plt.grid(True)
plt.show()