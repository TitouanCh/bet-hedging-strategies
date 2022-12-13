import numpy as np
import matplotlib.pyplot as plt
from scipy.special  import erfinv

f, (ax1, ax2) = plt.subplots(1, 2)
t = np.linspace(0, 2.5, 60)

lognormal_quantile = lambda m, v, p : np.exp(m + np.sqrt(2 * v) * erfinv(2 * p - 1))
lognormal_density = lambda m, v, x : (1/(x * np.sqrt(v) * np.sqrt(2 * np.pi)) * np.exp(-np.power(np.log(x) - m, 2)/(2 * v)))

ax1.plot(t, lognormal_density(0, 1, t), c="darkgrey")
ax1.set_title("Densité d'une loi Lognormale de paramètres (0, 1)")
ax2.plot(t, lognormal_quantile(0, 1, t), c="darkgrey")
ax2.set_title("Fonction quantile d'une loi Lognormale de paramètres (0, 1)")

plt.show()