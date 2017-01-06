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

"""
This code calculates extinction (absorption) of cloud gas by function of water droplet particle size

"""

N = 10000*10e6      # assuming number concentration of 10000/cm^3 --> 10^9/m^3
L = 1               # length of a pipe, assuming "about a meter" is a meter (Opintomoniste)


m = 1.33 + 0.001j
lambda_var = 635
dp = np.logspace(1,4,1000)

Qext = np.zeros(dp.size)

for k in range(dp.size):
    result = toolbox.mie(m, np.pi*dp[k] / lambda_var)
    Qext[k] = result[3]

#plt.plot(dp, Qext)
#plt.xscale("log")
#plt.yscale("log")
#plt.show()


sigma_ext=(np.pi*N*(dp*1e-9)**2*Qext)/4   # dp in nanometers
ext = 1- np.exp(-sigma_ext*L)

plt.plot(dp, ext)
plt.show()