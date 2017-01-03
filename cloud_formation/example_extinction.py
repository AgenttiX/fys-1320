# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017
# Strongly based on:
# Tampere University of Technology course FYS-1320 ESIMERKKI_EKSTINKTIOTEHOKKUUS.m by Joni Kalliokoski 2014-09-09

import toolbox
import numpy as np
import matplotlib.pyplot as plt

m = 1.33 + 0.001j
lambda_var = 635
dp = np.logspace(1,4,1000)

Qext = np.zeros(dp.size)

for k in range(dp.size):
    result = toolbox.mie(m, np.pi*dp[k] / lambda_var)
    Qext[k] = result[3]

plt.plot(dp, Qext)
plt.xscale("log")
plt.yscale("log")
plt.show()
