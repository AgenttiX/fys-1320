# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

# pylint: disable=wrong-import-order
import toolbox

import time

import numpy as np
import pyqtgraph as pg

# This is a value computed from our student numbers
studnum_a = 3

# Constants
temp = 296.15               # 23 deg C, from instructions
diff = 0.282e-4             # Diffusion coefficient for water (m^2/s)
r = 8.3144621               # (Pa*m^3)/(mol*K), from example code
surface_tension = 27.8e-3   # (N/m), From example code
m_mol = 18.016e-3           # Molar mass of water (kg/mol)
rho_wat = 998.20            # Density of water (kg/m^3)
evap_E = 2260e3             # Evaporation energy of water (J/kg)
thermal_con_air = 0.0257    # Thermal conductivity of air (W/(m*K))
rho_air = 1000              # Air density (kg/m^3)
heat_capacity_ratio = 1.4   # For (dry) air, from https://en.wikipedia.org/wiki/Heat_capacity_ratio


# Initial partial pressure for water
# p0 = 1.2*toolbox.water_pvap(temp)

# Molar volume of water
v_mol = m_mol / rho_wat

# Particle initialisation
#particle_sizes = np.array([5, 10, 20, 50, 100])*studnum_a*1e-9
#print("Particle diameters (m)")
#print(particle_sizes)

# AM_ch17.pdf table 17-1 data for water
water_a = 10.23
water_b = 1750
water_c = 38

# Set background white and foreground black
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


# Initialise figure environment
app = pg.mkQApp()
win = pg.GraphicsWindow(title="Adiabatic process")
# win.setWindowTitle("Cloud formation")
pg.setConfigOptions(antialias=True)


T_i = 296.15
p_i = 100000

p_f_vec = np.linspace(95000, 100000,1000)

T_f_vec = toolbox.final_temp(T_i, p_f_vec, p_i, heat_capacity_ratio)
p_sf_vec = toolbox.saturation_pressure(T_f_vec, water_a, water_b, water_c)

P_df_vec = toolbox.saturation_pressure(T_i, water_a, water_b, water_c) * (p_f_vec/p_i)

plot_p_s_by_p_f = win.addPlot()
plot_p_s_by_p_f.addLegend()
# Norm<span id="vertical-align: sub;">Sub</span>')
plot_p_s_by_p_f.plot(p_f_vec, p_sf_vec, pen=pg.mkPen((255, 0, 0), width=1.5), name='Kylläinen höyrynpaine p_s')
# Norm<span id="vertical-align: sub;">Sub</span>')
# ylä ja alaindeksien säätäminen ei onnistunut
plot_p_s_by_p_f.plot(p_f_vec, P_df_vec, pen=pg.mkPen((0, 255, 0), width=1.5), name='Höyryn osapaine p_d')

plot_p_s_by_p_f.setLabel("left", "pₛ", "Pa")
plot_p_s_by_p_f.setLabel("bottom", "P_f", "Pa")


# PyQtPlot main loop for graphs
app.exec_()

# This prevents a segmentation fault
# Don't ask us why
time.sleep(1)
