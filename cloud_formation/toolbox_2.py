# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

# Toolbox_2 was created to move functions out from data_analysis.py

import toolbox      # To avoid matlab from loading, comment this and uncomment line in simulate_extinction()
from math import factorial
import numpy as np

# From http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
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
    border_margin = np.int(index_spacing / index_skip) + 1 # for cropping off the ends of an array to avoid out_of_range
    greatest_difference = -1
    index_of_greatest_diff = -1
    i_of_greatest = -1      # i is sparce index, meaning for every i, there's {index_skip=200} normal indexes

    for i in range(border_margin, (np.int(data.size/index_skip)-border_margin)):
        index = i * index_skip
        # finds the difference of data around the index.
        difference = data[index - index_spacing] - data[index + index_spacing]
        if (difference > greatest_difference):
            greatest_difference = difference
            index_of_greatest_diff = index


    # Fine tuning of index to the beginning of decline:

    # Minium and maxium of data before the drop occurs. Deviation is caused by random noice.
    minium_before_drop = np.amin(data[index_of_greatest_diff-4000 : index_of_greatest_diff-2000])
    maxium_before_drop = np.amax(data[index_of_greatest_diff-4000 : index_of_greatest_diff-2000])
    threshold = minium_before_drop - (maxium_before_drop-minium_before_drop)

    fine_tuned_index_of_drop = -1
    for i in range(index_of_greatest_diff-2000, index_of_greatest_diff+1000):
        if (data[i] < threshold):
            fine_tuned_index_of_drop = i
            break
    if (fine_tuned_index_of_drop == -1):
        raise "Something went wrong in fine-tuning the index of pressure drop!"
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
    drop_heigth = np.mean(measurement.p_diff[drop - 6000 : drop - 4000]) - np.mean(measurement.p_diff[drop + 4000 : drop + 6000])

    initial_pressure = final_pressure + drop_heigth

    return initial_pressure*1000, final_pressure*1000   # Note that pressures in arrays are in kPa

def simulate_extinction(particle_size, p_i, p_f, particle_dens, tmax = 10):
    """
    Simulates particle growth extinction with given parameters. Growing particles are all the same size
    :param particle_size:   size of particle
    :param p_i:             initial pressure
    :param p_f:             final pressure
    :param particle_dens:   the number of particles in unit volume  (#/m^3)
    :param tmax:            the maximum time up which to compute to (s)
    :return:
    """
    # import toolbox # Uncomment this and comment the line in beginning of file to avoid matalb from loading

    # Constants
    temp_i = 296.15  # (K), 23 deg C, from instructions
    diff = 0.282e-4  # Diffusion coefficient for water (m^2/s)
    gas_const = 8.3144621  # (Pa*m^3)/(mol*K), from example code
    surface_tension = 72.8e-3  # (N/m), From example code                      FIXED
    m_mol = 18.016e-3  # Molar mass of water (kg/mol)
    rho_wat = 998.20  # Density of water (kg/m^3)
    evap_E = 2260e3  # Evaporation energy of water (J/kg)
    thermal_con_air = 0.0257  # Thermal conductivity of air (W/(m*K))
    heat_capacity_ratio = 1.4  # For (dry) air, from https://en.wikipedia.org/wiki/Heat_capacity_ratio
    m = 1.33 + 0.001  # Refractive index of water
    water_a = 10.23
    water_b = 1750
    water_c = 38
    wavelength = 635  # Wavelength of our laser
    length = 1  # (m)

    partial_pressure = (p_f/p_i) * toolbox.water_pvap(temp_i)   # Partial water pressure after adiabatic expansion
    temp_f = toolbox.final_temp(temp_i, p_f, p_i, heat_capacity_ratio)  # Temperature after adiabatic expansion

    t, dp, pw = toolbox.solve_growth(temp_f, diff, m_mol, evap_E, thermal_con_air, rho_wat, surface_tension,
                                           particle_dens, tmax, particle_size, partial_pressure)
    q_ext = toolbox.q_ext(dp, m, wavelength)
    sigma_ext = toolbox.extinction_factor(particle_dens, dp, q_ext)
    ext = toolbox.extinction(sigma_ext, length)

    smallest_growing_particle =  toolbox.minimum_particle_diameter_2(p_i, p_f, temp_f, heat_capacity_ratio,
                                           water_a, water_b, water_c, m_mol, surface_tension, rho_wat)
    return ext, t, smallest_growing_particle