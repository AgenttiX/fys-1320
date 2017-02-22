# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

# This file contains functions for project.py

import numpy as np
import scipy.integrate as integrate
import inspect


# Natural constants
m_e = 9.10939e-34
e = 1.60218e-19
c = 2.99792458e8

# Other constants

# Most common energy for Am241 is unknown, for now let's use one from wikipedia: https://en.wikipedia.org/wiki/Americium-241
E_1f = 59.5409e3 * e


def E_2f_fun(theta, E=E_1f):
    """
    Photons energy after compoton scatterin, (using energy E_1f)
    :param theta: angle for scattered photon
    :return:
    """
    return E / (1 + (E / (m_e * c ** 2)) * (1 - np.cos(theta)))

def E_2e_fun(theta, E=E_1f):
    """
    Electron energy after compoton scatterin, (using energy E_1f)
    :param theta: angle for scattered photon
    :return:
    """
    return E / (((m_e * c ** 2) / E) * (1 / (1 - np.cos(theta))) + 1)

def gauss(x,mu,var):
    """
    Gauss distribution value at x
    :param x:
    :param mu: expected value
    :param var: variance, (sigma^2)
    :return:
    """
    return 1/(np.sqrt(2*var*np.pi)) * np.exp(- (x-mu)**2/(2*var))

def gauss_cumulative(x,mu,var):
    """
        Gauss cumulative distribution value at x
        :param x:
        :param mu: expected value
        :param var: variance, (sigma^2)
        :return:
        """
    # Gauss distribution function:
    disrt_fun = lambda x_1: gauss(x_1, mu, var)
    # Cumulative function for cell element by numeric integration:
    cumul_fun = lambda x_2: integrate.quad(disrt_fun, -np.inf, x_2)[0] # [0] integral, [1] absolute error in integral
    # Vectorized function:
    cumul_vec_fun = np.vectorize(cumul_fun)

    cumul_vec = cumul_vec_fun(x)
    return cumul_vec

def cauchy(x,x0,gamma):
    """
    Cauchy distribution aka. Lorentz distribution
    :param x:
    :param x0: location parameter (not the expected value, as it is undefined)
    :param gamma: half-width at half maxium (not variance, as it is undefined)
    :return:
    """
    return 1/( np.pi*gamma*(1 + ((x-x0)/gamma)**2) )

def cauchy_cumulative(x,x0,gamma):
    """
    Cauchy distribution aka. Lorentz distribution
    :param x:
    :param x0: location parameter (not the expected value, as it is undefined)
    :param gamma: half-width at half maxium (not variance, as it is undefined)
    :return:
    """
    return 1/np.pi * np.arctan((x-x0)/gamma) + 0.5
