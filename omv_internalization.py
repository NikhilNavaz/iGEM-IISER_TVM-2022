#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

k_d = 22*(10)**(-12)  # https://pubmed.ncbi.nlm.nih.gov/20696930/
k_a = 4.54*(10)**(10)  # https://pubmed.ncbi.nlm.nih.gov/20696930/
k_i = 0.394  # https://www.frontiersin.org/articles/10.3389/fmicb.2020.00057/full

# Da means C
# Db means Cbs
# Dc means Cb
# Dd means Ci

# Da = 100
# dDb/dt = -k_a*Db + k_d*Dc
# dDc/dt = k_a*Db*Da - k_d*Dc - k_i*Dc
# dDd/dt = k_i*Dc
# Solving the system of ODEs


def odes(x, t):
    Db = x[0]
    Dc = x[1]
    Dd = x[2]
    dDbdt = -k_a*Db + k_d*Dc
    dDcdt = k_a*Db*10**(-4) - k_d*Dc - k_i*Dc
    dDddt = k_i*Dc
    return [dDbdt, dDcdt, dDddt]


x0 = [10**(-4), 10**(-4), 10**(-4)]

t = np.linspace(0, 0.1, 1000)
x = odeint(odes, x0, t)

Db = x[:, 0]
Dc = x[:, 1]
Dd = x[:, 2]

plt.plot(t, Dd*10**3, label='Ci')
plt.xlabel('Concentration outside (uM)')
plt.ylabel('OMV concentration internalized (uM)')
plt.legend()
plt.show()
