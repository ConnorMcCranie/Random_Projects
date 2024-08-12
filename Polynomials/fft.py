from scipy.fft import fft, ifft
import numpy as np
import matplotlib.pyplot as plt

def scaled_pow(N,p):
    x = np.array([n**p for n in range(N)])
    s = sum(x)
    return x/s


fig, ax = plt.subplots()
ax.plot([n for n in range(21)],fft(scaled_pow(21,3)).real)
ax.plot([n for n in range(21)],fft(scaled_pow(21,3)).imag)
plt.show()