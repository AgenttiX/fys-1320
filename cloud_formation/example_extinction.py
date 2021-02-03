# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017
# Strongly based on:
# Tampere University of Technology course FYS-1320 ESIMERKKI_EKSTINKTIOTEHOKKUUS.m by Joni Kalliokoski 2014-09-09

# Our work is licensed with the MIT license.
# However, the license of the example code our work is based on is unclear and thereby so is the license
# for those parts of this code that are based on it

import toolbox
import numpy as np
import matplotlib.pyplot as plt

m = 1.33 + 0.001j   # Refractive index for water
lambda_var = 635    # Laser wavelength (nm)

# Vector for particle sizes (Numpy array)
dp = np.logspace(1, 4, 1000)

Qext = np.zeros(dp.size)

# Fill the Qext vector with data
for k in range(dp.size):
    # Compute the Mie parameters using the toolbox
    result = toolbox.mie(m, np.pi*dp[k] / lambda_var)
    # The toolbox.mie returns a vector, in which the fourth item is the Qext. Other items aren't needed.
    Qext[k] = result[3]

# Draw a figure with logarithmic axes
plt.plot(dp, Qext)
plt.xscale("log")
plt.yscale("log")
plt.show()
