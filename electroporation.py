import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

C = 0.789322
k = 1.10092*(10**9)
Nep0 = 1


def nep(n0, t, Nep0=1):
    return (n0-(n0-Nep0*C)*np.exp(-k*C*t/n0))/C


n0 = np.linspace(0.000000001, 5, 100)
t = np.linspace(0.000000001, 0.0000001, 100)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
t, N0 = np.meshgrid(t, n0)
Nept = nep(N0, t, Nep0)

ax.set_xlabel('N0')
ax.set_ylabel('time in sec')
ax.set_zlabel('Nep0')
surf = ax.plot_surface(N0, t, Nept, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()
