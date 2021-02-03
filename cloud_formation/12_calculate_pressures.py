# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

"""
This file computes required initial pressures for twelve particles on range 6 - 200nm.
"""

# pylint: disable=wrong-import-order
import toolbox

import numpy as np


# Constants
temp = 296.15               # (K), 23 deg C, from instructions
surface_tension = 72.8e-3   # (N/m), From example code                      FIXED
m_mol = 18.016e-3           # Molar mass of water (kg/mol)
rho_wat = 998.20            # Density of water (kg/m^3)
heat_capacity_ratio = 1.4   # For (dry) air, from https://en.wikipedia.org/wiki/Heat_capacity_ratio

# AM_ch17.pdf table 17-1 data for water
water_a = 10.23
water_b = 1750
water_c = 38

final_pressure = 99000  # Pa


# 12 Particle sizes to compute initial pressures
particle_sizes = np.linspace(6e-9, 200e-9, 12)
print("Particle diameters (nm)")
print(np.round(particle_sizes*1e9, decimals=1))
print("")


def particle_diameter_withvalues(p_i):
    return toolbox.minimum_particle_diameter_2(p_i, final_pressure, temp, heat_capacity_ratio,
                                               water_a, water_b, water_c, m_mol, surface_tension, rho_wat)


init_pressure_vec = np.arange(99050, 109000)
min_dp_vec = particle_diameter_withvalues(init_pressure_vec)


# Uses interpolation to calculate zero-point for function particle_diameter_withvalues()
init_pressures_for_particles = np.interp(particle_sizes, np.flipud(min_dp_vec),
                                         np.flipud(init_pressure_vec), left=-1, right=-2)

print("Initial pressures for our 12 particles (Pa)")
print(np.round(init_pressures_for_particles))

print("Initial pressure for 5nm")
print(np.interp(
    5e-9, np.flipud(min_dp_vec),
    np.flipud(init_pressure_vec), left=-1, right=-2
))
