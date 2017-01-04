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
particle_dens = 1e4         # particle density (#/m^3)
tmax = 4                    # The maximum time up which to compute to (s)

# Initial partial pressure for water
p0 = 1.2*toolbox.water_pvap(temp)

# Molar volume of water
v_mol = m_mol / rho_wat

# Particle initialisation
particle_sizes = np.array([5, 10, 20, 50, 100])*studnum_a*1e-9
print("Particle diameters (m)")
print(particle_sizes)


# ----- USELESS -----

# Equation 17-1
def p_s(t, a, b, c):
    """

    :param t: temperature (K)
    :param a: constant from table 17-1
    :param b: constant from table 17-1
    :param c: constant from table 17-1
    :return: vapor pressure at the saturation or equilibrium condition
    """
    return 10**(a - b / (t - c))

"""
# Equation 17-2
def p_s_eq17_2(t, b, c):
    return 10**((-52.3*b)/t + c)
"""

# AM_ch17.pdf table 17-1 data for water
water_a = 10.23
water_b = 1750
water_c = 38

# Equation 17-3
def p_d(p_s, nu, gamma, r, t, d_p):
    """
    vapor pressure at the saturation on equilibrium condition for particles of this size
    :param p_s: vapor pressure at the saturation on equilibrium condition
    :param nu: molar volume of the liquid
    :param gamma: surface tension
    :param r: gas constant
    :param t: temperature
    :param d_p: particle diameter
    :return:
    """
    return p_s * np.exp((4*nu*gamma)/(r*t*d_p))

# -----


# 5)

def saturation_ratio(nu, gamma, r, t, d_p):
    """
    Compute the saturation ratio S_R = p_d/p_s for particles of given diameter
    :param nu: molar volume of the liquid
    :param gamma: surface tension
    :param r: gas constant
    :param t: temperature
    :param d_p: particle diameter
    :return: saturation ratio
    """
    return np.exp((4 * nu * gamma) / (r * t * d_p))

saturation_ratios = saturation_ratio(v_mol, surface_tension, r, temp, particle_sizes)
print("Saturation (=pressure) ratios")
print(saturation_ratios)

print("Saturation partial pressure of water at room temperature (Pa)")
print(toolbox.water_pvap(temp))

part_pressures_on_particles = saturation_ratios*toolbox.water_pvap(temp)
print("Partial pressures on particle surfaces (Pa)")
print(part_pressures_on_particles)


# 8)
particle_sizes_for8 = np.array([5, 10, 20])*studnum_a*1e-9

# Compute particle growth
t_5, dp_5, pw_5 = toolbox.solve_growth(temp, diff, m_mol, evap_E, thermal_con_air, rho_air, surface_tension, particle_dens, tmax, particle_sizes_for8[0], p0)
t_10, dp_10, pw_10 = toolbox.solve_growth(temp, diff, m_mol, evap_E, thermal_con_air, rho_air, surface_tension, particle_dens, tmax, particle_sizes_for8[1], p0)
t_20, dp_20, pw_20 = toolbox.solve_growth(temp, diff, m_mol, evap_E, thermal_con_air, rho_air, surface_tension, particle_dens, tmax, particle_sizes_for8[2], p0)

# Draw figures
app = pg.mkQApp()
win = pg.GraphicsWindow(title="example_dp")
win.setWindowTitle("example_dp")
pg.setConfigOptions(antialias=True)

p1 = win.addPlot(title="d_p (m)")
p1.plot(t_5, dp_5, pen=(255,0,0))
p1.plot(t_10, dp_10,pen=(0,255,0))
p1.plot(t_20, dp_20, pen=(0,0,255))
p2 = win.addPlot(title="P_w (Pa)")
p2.plot(t_5, pw_5, pen=(255,0,0))
p2.plot(t_10, pw_10, pen=(0,255,0))
p2.plot(t_20, pw_20, pen=(0,0,255))

app.exec_()
