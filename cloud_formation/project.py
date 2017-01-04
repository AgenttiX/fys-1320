# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

"""
On the slide 7 of the first lecture it was mentioned that alternative programming environments can be used.
We thereby presume that the MATLAB specific instructions of this particular project are simply to prevent students from
making the project unnecessarily difficult for themselves.
However, we found out that there is a MATLAB API for Python and since both of us are familiar with Python and
passionately hate MATLAB, we chose to utilise Python instead.
"""

import toolbox
import numpy as np
import pyqtgraph as pg
import time
import scipy.interpolate

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
particle_sizes = np.array([5, 10, 20, 50, 100])*studnum_a*1e-9
print("Particle diameters (m)")
print(particle_sizes)

# AM_ch17.pdf table 17-1 data for water
water_a = 10.23
water_b = 1750
water_c = 38


# Initialise figure environment
app = pg.mkQApp()
win = pg.GraphicsWindow(title="Cloud formation")
# win.setWindowTitle("Cloud formation")
pg.setConfigOptions(antialias=True)


# 5)
print("\n----- 5 -----")

print()
saturation_ratios = toolbox.saturation_ratio(v_mol, surface_tension, r, temp, particle_sizes)
print("Saturation (=pressure) ratios")
print(saturation_ratios)

print()
print("Relative humidity (%)")
print(saturation_ratios*100)

print()
print("Saturation partial pressure of water at room temperature (Pa)")
print(toolbox.water_pvap(temp))

print()
part_pressures_on_particles = saturation_ratios*toolbox.water_pvap(temp)
print("Partial pressures on particle surfaces (Pa)")
print(part_pressures_on_particles)


def pressure_difference(nu, gamma, r, t, d_p, p_s):
    return (np.exp((4 *nu*gamma)/(r*t*d_p)) - 1)*p_s

print()
pressure_differences = pressure_difference(v_mol, surface_tension, r, temp, particle_sizes, toolbox.water_pvap(temp))
print("Pressure differences (Pa)")
print(pressure_differences)

# 7)
print("\n----- 7 -----")

final_pressure = 99000  # Pa


def particle_diameter_withvalues(p_i):
    return toolbox.minimum_particle_diameter_2(p_i, final_pressure, temp, heat_capacity_ratio, water_a, water_b, water_c, m_mol, surface_tension, rho_wat)

# Draw a figure
init_pressure_vec = np.arange(99050, 100500)
min_dp_vec = particle_diameter_withvalues(init_pressure_vec)

plot_dp_by_pi = win.addPlot(title="d_p (m) py P_i (Pa)")
plot_dp_by_pi.plot(init_pressure_vec, min_dp_vec)

# Compute intersections by interpolation
# init_pressures_for_particles = scipy.optimize.fixed_point(particle_diameter_withvalues, particle_sizes)

# The derivative of the function must be positive, thereby np.flipud() is necessary
init_pressures_for_particles = np.interp(particle_sizes, np.flipud(min_dp_vec), np.flipud(init_pressure_vec), left=-1, right=-2)

print()
print("Initial pressures for our particles (Pa)")
print(init_pressures_for_particles)

print()
print("Pressure ratios for our particles (Pa)")
print(final_pressure / init_pressures_for_particles)

print()
print("Pressure differences for our particles (Pa)")
print(init_pressures_for_particles - final_pressure)


# 8)
# particle_sizes_for8 = np.array([5, 10, 20])*studnum_a*1e-9

tmax = 4                    # The maximum time up which to compute to (s)
particle_dens = 1e4         # particle density (#/m^3)

# Compute particle growth
t_5, dp_5, pw_5 = toolbox.solve_growth(temp, diff, m_mol, evap_E, thermal_con_air, rho_air, surface_tension, particle_dens, tmax, 20*studnum_a*1e-9, part_pressures_on_particles[0])
t_10, dp_10, pw_10 = toolbox.solve_growth(temp, diff, m_mol, evap_E, thermal_con_air, rho_air, surface_tension, particle_dens, tmax, 20*studnum_a*1e-9, part_pressures_on_particles[1])
t_20, dp_20, pw_20 = toolbox.solve_growth(temp, diff, m_mol, evap_E, thermal_con_air, rho_air, surface_tension, particle_dens, tmax, 20*studnum_a*1e-9, part_pressures_on_particles[2])


# Particle diameter
plot_dp = win.addPlot(title="d_p (m) by t (s)")
plot_dp.plot(t_5, dp_5, pen=(255,0,0))
plot_dp.plot(t_10, dp_10,pen=(0,255,0))
plot_dp.plot(t_20, dp_20, pen=(0,0,255))

# Partial pressure of water
plot_pw = win.addPlot(title="p_w (Pa) by t(s)")
plot_pw.plot(t_5, pw_5, pen=(255,0,0))
plot_pw.plot(t_10, pw_10, pen=(0,255,0))
plot_pw.plot(t_20, pw_20, pen=(0,0,255))

app.exec_()

# This prevents a segmentation fault
# Don't ask us why
time.sleep(1)
