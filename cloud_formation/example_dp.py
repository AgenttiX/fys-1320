# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017
# Strongly based on Tampere University of Technology course FYS-1320 ESIMERKKI_dp.m by Joni Kalliokoski 2014-08-14

import toolbox
# import pyqtgraph as pg
import matplotlib.pyplot as plt

T = 296.15
D = 0.282e-4
M = 18.016e-3
L = 2260e3
ka = 0.0257
rho = 1000
gamma = 27.8e-3
Ntot = 1e10
tmax = 2

dp0 = 100e-9
p0 = 1.2*toolbox.water_pvap(T)

t, dp, pw = toolbox.solve_growth(T, D, M, L, ka, rho, gamma, Ntot, tmax, dp0, p0)

"""
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

"""
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

plt.figure(1)
plt.subplot(121)
plt.plot(t, dp)
plt.subplot(122)
plt.plot(t, pw)
plt.show()
