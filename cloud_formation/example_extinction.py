# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017
# Strongly based on:
# Tampere University of Technology course FYS-1320 ESIMERKKI_EKSTINKTIOTEHOKKUUS.m by Joni Kalliokoski 2014-09-09

# Our work is licensed with Creative Commons Attribution 4.0 International
# https://creativecommons.org/licenses/by/4.0/
# However, the license of the example code our work is based on is unclear and thereby so is the license
# for those parts of this code that are based on it

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
