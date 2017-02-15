# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

# This program serves as a collection of supplementary functions for data_analysis.py

import toolbox  # To prevent Matlab from loading, comment this line and uncomment the line in simulate_extinction()
from math import factorial
import numpy as np


# From http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
# This function is thereby NOT covered by our licensing
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.

    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


def remove_noice(data, noice_reduction_number):
    if(noice_reduction_number == 0):
        return data
    if(noice_reduction_number < 0):
        raise ValueError("Noice reduction can't be negative")
    else:
        return savitzky_golay(data, (noice_reduction_number*50 + 1), 3) # window size 251, polynomial order 3


def flip_and_normalize(data):
    index = find_drop_index(data)
    zero_level = np.mean(data[index-10000:index-6000])
    return 1-data/zero_level


def find_drop_index(data):
    """
    Finds the index of beginning of the largest decline in vector. Intended for data=measurement.p_diff, in which the biggest
    difference drop indicates the moment of pressure release.
    :param data: type of array, works best for measurement.p_diff
    :return: index of the beginning of the pressure drop
    """
    # This funtion finds the coarse index of biggest drop first, and then it fine tunes the index to the beginning of drop.

    # Finding the biggest decline with coarse method, this is done by iterating differences on some interval.
    index_skip = 200        # every single index is not tested for performance reasons,
    index_spacing = 500     # the time interval between which difference is tested
    border_margin = np.int(index_spacing / index_skip) + 1  # for cropping off the ends of an array to avoid out_of_range
    greatest_difference = -1
    index_of_greatest_diff = -1

    # i is sparse index, meaning for every i, there's {index_skip=200} normal indexes
    for i in range(border_margin, (np.int(data.size/index_skip)-border_margin)):
        index = i * index_skip
        # finds the difference of data around the index.
        difference = data[index - index_spacing] - data[index + index_spacing]
        if difference > greatest_difference:
            greatest_difference = difference
            index_of_greatest_diff = index

    # Fine tuning of index to the beginning of decline:
    # Minimum and maximum of data before the drop occurs. Deviation is caused by random noise.
    minium_before_drop = np.amin(data[index_of_greatest_diff-4000: index_of_greatest_diff-2000])
    maxium_before_drop = np.amax(data[index_of_greatest_diff-4000: index_of_greatest_diff-2000])
    threshold = minium_before_drop - (maxium_before_drop-minium_before_drop)

    fine_tuned_index_of_drop = -1
    for i in range(index_of_greatest_diff-2000, index_of_greatest_diff+1000):
        if data[i] < threshold:
            fine_tuned_index_of_drop = i
            break
    if fine_tuned_index_of_drop == -1:
        raise UserWarning("Something went wrong in fine-tuning the index of pressure drop!")
    fine_tuned_index_of_drop -= 150     # -150 is approximation that sets index somewhat just before the drop starts

    return fine_tuned_index_of_drop


def get_pressure_change(measurement):
    """
    Returns initial and final pressures.

    :param measurement: must be instance of class Measurement
    :return: initial_pressure, final_pressure
    """
    drop = find_drop_index(measurement.p_diff)

    final_pressure = np.mean(measurement.p_abs[drop + 6000: drop + 16000])
    drop_height = np.mean(measurement.p_diff[drop - 6000: drop - 4000]) - np.mean(measurement.p_diff[drop + 4000: drop + 6000])

    initial_pressure = final_pressure + drop_height

    return initial_pressure*1000, final_pressure*1000   # Note that pressures in arrays are in kPa


def simulate_extinction(particle_size, p_i, p_f, particle_dens, tmax=10, saturation=1.0):
    """
    Simulates particle growth extinction with given parameters. Growing particles are all the same size
    :param particle_size:   size of particle
    :param p_i:             initial pressure
    :param p_f:             final pressure
    :param particle_dens:   the number of particles in unit volume  (#/m^3)
    :param tmax:            the maximum time up which to compute to (s)
    :param saturation       proportional fraction of partial water pressure from its maximum value, a.k.a vapor quality
    :return:
    """
    # To avoid Matlab from loading until executing this function, uncomment
    # the next line and comment the line in beginning of file
    # import toolbox

    # Constants
    temp_i = 296.15  # (K), 23 deg C, from instructions
    diff = 0.282e-4  # Diffusion coefficient for water (m^2/s)
    surface_tension = 72.8e-3  # (N/m), From example code
    m_mol = 18.016e-3  # Molar mass of water (kg/mol)
    rho_wat = 998.20  # Density of water (kg/m^3)
    evap_e = 2260e3  # Evaporation energy of water (J/kg)
    thermal_con_air = 0.0257  # Thermal conductivity of air (W/(m*K))
    heat_capacity_ratio = 1.4  # For (dry) air, from https://en.wikipedia.org/wiki/Heat_capacity_ratio
    m = 1.33 + 0.001  # Refractive index of water
    wavelength = 635  # Wavelength of our laser
    length = 1  # (m)

    partial_pressure = (p_f/p_i) * toolbox.water_pvap(temp_i) * saturation   # Partial water pressure after adiabatic expansion
    temp_f = toolbox.final_temp(temp_i, p_f, p_i, heat_capacity_ratio)  # Temperature after adiabatic expansion

    t, dp, pw = toolbox.solve_growth(temp_f, diff, m_mol, evap_e, thermal_con_air, rho_wat, surface_tension,
                                     particle_dens, tmax, particle_size, partial_pressure)
    q_ext = toolbox.q_ext(dp, m, wavelength)
    sigma_ext = toolbox.extinction_factor(particle_dens, dp, q_ext)
    ext = toolbox.extinction(sigma_ext, length)

    return ext, t


def minimum_particle_diameter(p_i, p_f, saturation=1.0):
    """
    Returns the smallest growing particle size.
    :param p_i:
    :param p_f:
    :param saturation:
    :return:
    """
    temp_i = 296.15  # (K), 23 deg C, from instructions
    heat_capacity_ratio = 1.4  # For (dry) air, from https://en.wikipedia.org/wiki/Heat_capacity_ratio
    water_a = 10.23
    water_b = 1750
    water_c = 38
    m_mol = 18.016e-3  # Molar mass of water (kg/mol)
    surface_tension = 72.8e-3  # (N/m), From example code
    rho_wat = 998.20  # Density of water (kg/m^3)

    temp_f = toolbox.final_temp(temp_i, p_f, p_i, heat_capacity_ratio)  # Temperature after adiabatic expansion

    smallest_growing_particle = toolbox.minimum_particle_diameter_2(p_i, p_f, temp_f, heat_capacity_ratio,
                                                                    water_a, water_b, water_c, m_mol, surface_tension,
                                                                    rho_wat, saturation)

    return smallest_growing_particle


def extinction_factor(extinction_fraction):
    """
    Calculates the extinction factor from equation:   extinction_fraction = exp( -extinction_factor * L)
    Symbol: sigma_ext
    :param extinction_fraction:     (1-I/I0), falls in range 0-1
    :return:
    """
    l = 1           # length of tube (m)

    return -np.log(1-extinction_fraction) / l


def particle_count(sigma_ext, p_size, q_ext):
    """
    Returns the particle concentration in #/m^3
    Symbol N
    :param sigma_ext:   extinction factor (1/m)
    :param p_size:      particle size (m)
    :param q_ext:       extinction efficiency "Ekstinktiotehokkuus"
    :return:
    """
    return (sigma_ext * 4) / (np.pi * p_size**2 * q_ext)


def particle_count_2(extinction_fraction):
    """
    Returns particle count. Combines two functions for better code appereance.
    :param extinction_fraction: (1-I/I0), falls in range 0-1
    :return: N
    """
    p_size = 1.76e-6  # Particle size at the beginning of first wrinkle. Project_mod.py's 5th plot.
    q_ext = 2.8  # extinction efficiency for above particle size. Project_mod.py's 4th plot.

    sigma_ext = extinction_factor(extinction_fraction)
    return particle_count(sigma_ext, p_size, q_ext)
