# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017
# Strongly based on Tampere University of Technology course FYS-1320 ESIMERKKI_dp.m by Joni Kalliokoski 2014-08-14

import toolbox
# import pyqtgraph as pg
import matplotlib.pyplot as plt

T = 296.15          # Temperature (K)
D = 0.282e-4        # Diffusion coefficient for water (m^2/s)
M = 18.016e-3       # Molar mass of water (kg/mol)
L = 2260e3          # Evaporation energy of water (J/kg)
ka = 0.0257         # Thermal conductivity of air (W/(m*K))
rho = 1000          # Air density (kg/m^3)
gamma = 27.8e-3     # Surface tension of water (N/m)
Ntot = 1e10         # particle density (#/m^3)
tmax = 2            # The maximum time up which to compute to (s)

# Initial particle diameter
dp0 = 100e-9
# Initial partial pressure for water
p0 = 1.2*toolbox.water_pvap(T)

# Solve the particle growth
t, dp, pw = toolbox.solve_growth(T, D, M, L, ka, rho, gamma, Ntot, tmax, dp0, p0)

"""
# Debug printing
print(t)
print(dp)
print(pw)
print(t.dtype)
print(dp.dtype)
print(pw.dtype)
print(t.size)
print(dp.size)
print(pw.size)
"""

# Draw figures
plt.figure(1)
plt.subplot(121)
plt.plot(t, dp)
plt.subplot(122)
plt.plot(t, pw)
plt.show()

"""
# PyQtGraph version of figures

app = pg.mkQApp()
win = pg.GraphicsWindow(title="example_dp")
win.setWindowTitle("example_dp")
pg.setConfigOptions(antialias=True)

p1 = win.addPlot(title="d_p (m)")
p1.plot(t, dp)
p2 = win.addPlot(title="P_w (Pa)")
p2.plot(t, pw)

app.exec_()
"""
