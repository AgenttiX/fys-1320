# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

# This file contains functions for project.py

import numpy as np
import scipy.integrate as integrate


# Natural constants
m_e = 9.10939e-31
e = 1.60218e-19
c = 2.99792458e8
h = 6.6260755e-34

# Other constants

# Most common energy for Am241 is unknown, for now let's use one from Wikipedia:
# https://en.wikipedia.org/wiki/Americium-241
e_1f = 59.5409e3 * e


def e_2f_fun(theta, E_init=e_1f):
    """
    Photon energy after Compton scattering, (using energy e_1f)
    :param theta: angle for scattered photon
    :param E_init: initial photon energy
    :return:
    """
    return E_init / (1 + (E_init / (m_e * c ** 2)) * (1 - np.cos(theta)))


def e_2e_fun(theta, e_init=e_1f):
    """
    Electron energy after Compton scattering, (using energy e_1f)
    :param theta: angle for scattered photon
    :param e_init: initial photon energy
    :return:
    """
    return e_init / (((m_e * c ** 2) / e_init) * (1 / (1 - np.cos(theta))) + 1)


def klein_nishina(theta, E_init=e_1f):
    """

    :param theta:
    :param E_init:
    :return:
    """
    alfa = 1/137
    r_c = h/(2*np.pi) * m_e * c

    def P(theta):
        return 1 / (1 + (E_init / (m_e * c ** 2)) * (1 - np.cos(theta)))
    d_sigma_d_omega = alfa**2 * r_c**2 * P(theta)**2 * (P(theta) + 1/P(theta) -(np.sin(theta))**2) / 2

    return d_sigma_d_omega


def gauss(x, mu, var, a=1):
    """
    Gauss distribution value at x
    :param x:
    :param mu: expected value
    :param var: variance, (sigma^2)
    :param a: coefficient in cases total area != 1
    :return:
    """
    return a/(np.sqrt(2*var*np.pi)) * np.exp(- (x-mu)**2/(2*var))


def gauss_cumulative(x, mu, var, a=1):
    """
        Gauss cumulative distribution value at x
        :param x:
        :param mu: expected value
        :param var: variance, (sigma^2)
        :param a: coefficient in cases x_inf != 1
        :return:
        """
    # Gauss distribution function:
    disrt_fun = lambda x_1: gauss(x_1, mu, var, a)
    # Cumulative function for cell element by numeric integration:
    cumul_fun = lambda x_2: integrate.quad(disrt_fun, -np.inf, x_2)[0]  # [0] integral, [1] absolute error in integral
    # Vectorized function:
    cumul_vec_fun = np.vectorize(cumul_fun)

    cumul_vec = cumul_vec_fun(x)
    return cumul_vec


def cauchy(x, x0, gamma, a=1):
    """
    Cauchy distribution aka. Lorentz distribution
    :param x:
    :param x0: location parameter (not the expected value, as it is undefined)
    :param gamma: half-width at half maxium (not variance, as it is undefined)
    :param a: coefficient in cases total area != 1
    :return:
    """
    return a / (np.pi*gamma*(1 + ((x-x0)/gamma)**2))


def cauchy_cumulative(x, x0, gamma, a=1):
    """
    Cauchy distribution aka. Lorentz distribution
    :param x:
    :param x0: location parameter (not the expected value, as it is undefined)
    :param gamma: half-width at half maxium (not variance, as it is undefined)
    :param a: coefficient in cases x_inf != 1
    :return:
    """
    return a/np.pi * np.arctan((x-x0)/gamma) + 0.5


def pseudo_voigt(x, x0, var, gamma, a=1):
    """
    Pseudo-Voigt profile, is accurate to 1 % to regular Voigt
    :param x:
    :param x0: location parameter
    :param var: Gauss parameter: variance
    :param gamma: Caychy parameter: half-width at half maxium
    :param a: coefficient in cases when distribution is not normalized
    :return:
    """
    # if RunTimeWarnings are thrown, its because curvefitting does try negative numbers. Setting bounds >0 makes situation
    # worse.
    sigma = np.sqrt(var)
    f_G = 2*sigma*np.sqrt(2*np.log(2))
    f_L = 2*gamma

    f = (f_G**5 + 2.69269*f_G**4*f_L + 2.42843*f_G**3*f_L**2 + 4.47163*f_G**2*f_L**3 + 0.07842*f_G*f_L**4 + f_L**5)**(1/5)

    n = 1.36603*(f_L/f) - 0.47719*(f_L/f)**2 + 0.11116*(f_L/f)**3

    V_p = n*cauchy(x,x0,gamma,a) + (1-n)*gauss(x,x0,var,a)

    return V_p


def voigt_cauchy_percentage(var, gamma):
    """
    Calculates the the amount of Cauchy in Voigt. This number not real thing, but more just intuitive evaluation
    :param var: Gauss parameter: variance
    :param gamma: Caychy parameter: half-width at half maxium
    :return: fraction of Cauchy (0-1)
    """
    sigma = np.sqrt(var)
    f_G = 2 * sigma * np.sqrt(2 * np.log(2))
    f_L = 2 * gamma

    f = (f_G ** 5 + 2.69269 * f_G ** 4 * f_L + 2.42843 * f_G ** 3 * f_L ** 2 + 4.47163 * f_G ** 2 * f_L ** 3 + 0.07842 * f_G * f_L ** 4 + f_L ** 5) ** (1 / 5)

    n = 1.36603 * (f_L / f) - 0.47719 * (f_L / f) ** 2 + 0.11116 * (f_L / f) ** 3
    return n
