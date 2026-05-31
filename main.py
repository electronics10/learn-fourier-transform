import numpy as np
import matplotlib.pyplot as plt

T = 2.6
fs = 200
t = np.arange(0, T, 1/fs)
N = len(t)

f0 = 5
x = np.sin(2*np.pi*f0*t)

X = np.fft.fft(x)
freq = np.fft.fftfreq(N, 1/fs)


fig, ax = plt.subplots(2,1, figsize=(9,4))
ax[0].plot(t,x)
ax[0].set_title(f"sin(2pi {f0}t), fs={fs}, T={T}, N={N}")
ax[1].plot(freq, np.abs(X))
ax[1].set_title("Absolute value of the spectrum")

fig.tight_layout()
plt.show()
