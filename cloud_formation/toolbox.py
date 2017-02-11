# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017
# Strongly based on Tampere University of Technology course FYS-1320 example code

# Our work is licensed with Creative Commons Attribution 4.0 International
# https://creativecommons.org/licenses/by/4.0/
# However, the license of the example code our work is based on is unclear and thereby so is the license
# for those parts of this code that are based on it

# This toolbox utilises MATLAB Engine API for Python, since we couldn't translate some of the necessary code to Python
# https://se.mathworks.com/help/matlab/matlab-engine-for-python.html
# For Ubuntu the installation would be something like
# cd /usr/local/MATLAB/R2016b/extern/engines/python
# sudo python3 setup.py install

# It also requires the example code files to be in the same folder as this file

# Must be imported first to prevent errors
import matlab.engine

import numpy as np
import scipy.optimize
import scipy.special

# Start MATLAB
matlabeng = matlab.engine.start_matlab()

# Gas constant R
gas_const = 8.3144621


# ----- Based on Tampere University of Technology course FYS-1320 example code -----
# Not all of the functions have been tested

def dp_growth_rate(dp, t, rho, gamma, m, d, p, td):
    """
    Based on dp_growth_rate.m by Joni Kalliokoski 2014-08-13
    :param dp: particle diameter (m)
    :param t: environment temperature (K)
    :param rho: water density (kg / m^3)
    :param gamma: water surface tension (N/m)
    :param m: water molar mass (kg/mol)
    :param d: water diffusion coefficient (m^2 / s)
    :param p: environment water partial pressure (Pa)
    :param td: particle surface temperature (K)
    :return: speed of particle diameter change (m/s)
    """
    pd = water_pvap(td) * kelvin_ratio(dp, td, rho, gamma, m)
    lambda_var = free_path(t, d, m)

    return ((4 * d * m) / (gas_const * rho * dp)) * ((p / t) - (pd / td)) * fuchs_sutugin(dp, lambda_var)


def droplet_temp(t, d, m, l, ka, p, dp, rho, gamma):
    """
    Based on droplet_temp.m by Joni Kalliokoski 2014-08-13
    :param t: environment temperature (K)
    :param d: water diffusion coefficient (m^2 / s)
    :param m: water molar mass (kg/mol)
    :param l: water evaporation energy (J/kg)
    :param ka: thermal conductivity of air (W/(m*K))
    :param p: environment water partial pressure (Pa)
    :param dp: particle diamter (m)
    :param rho: water density (kg/m^3)
    :param gamma: water surface tension (N/m)
    :return: particle surface temperature (K)
    """
    lambda_var = free_path(t, d, m)
    # x = sympy.var("x")
    # return sympy.solvers.solve(x - T - ((D*M*L)/(R*ka))* ( (p/T) - \
    #       (water_pvap( x)*kelvin_ratio( dp,x,rho,gamma,M )/x) ) * FuchsSutugin( dp, lambda_var ), x)

    def droplet_temp_support(x):
        x - t - ((d * m * l) / (gas_const * ka)) * ((p / t) - (water_pvap(x) * kelvin_ratio(dp, x, rho, gamma, m) / x))\
                * fuchs_sutugin(dp, lambda_var)

    return scipy.optimize.fsolve(droplet_temp_support, t)


def dy_growth(t, y, temp, diff, m_mol, evap_e, ka, rho, gamma, ntot, dp0):
    """
    Based on dyKasvu.m by Joni Kalliokoski 2014-08-13
    :param t: time (s)
    :param y: y[0] = particle diameter (m), y[1] = environment water partial pressure (Pa)
    :param temp: environment temperature (K)
    :param diff: water diffusion coefficient (m^2/s)
    :param m_mol: water molar mass (kg/mol)
    :param evap_e: water evaporation energy (J/kg)
    :param ka: thermal conductivity of air (W/(m*K))
    :param rho: water density (kg/m^3)
    :param gamma: water surface tension (N/m)
    :param ntot: particle count in a volume (#/cm^3)
    :param dp0: initial particle diameter (m)
    :return: speeds of change: [0] = particle diameter (m), [1] = environment water partial pressure (Pa)

    It appears that the parameter t isn't used, even in the example code
    """
    dy = np.zeros(2)
    td = droplet_temp(temp, diff, m_mol, evap_e, ka, y[1], y[0], rho, gamma)
    dy[0] = dp_growth_rate(y[0], temp, rho, gamma, m_mol, diff, y[1], td)
    dy[1] = p_growth_rate(temp, rho, ntot, y[0], dy[0], m_mol)

    if (dy[0] < 0) and (y[0] <= dp0):
        dy[0] = 0
        dy[1] = 0

    return dy


def free_path(temp, diff, m_mol):
    """
    Calculates the free path for a molecule
    Based on free_path.m by Joni Kalliokoski 2014-08-13
    :param temp: temperature (K)
    :param diff: diffusion coefficient (m^2/s)
    :param m_mol: molar mass (kg/mol)
    :return: free path (m)
    """
    return 3*diff*np.sqrt((np.pi*m_mol)/(8*gas_const*temp))


def fuchs_sutugin(dp, lambda_var):
    """
    Calculates the Fuchs-Sutugin correction factor
    Based on FuchsSutugin.m by Joni Kalliokoski 2014-08-13
    :param dp:
    :param lambda_var: free path (m)
    :return: Fuchs-Sutugin correction factor
    """
    kn = 2*lambda_var / dp
    return (1 + kn) / (1 + 1.71*kn + 1.33*(kn**2))


def kelvin_ratio(dp, temp, rho, gamma, m_mol):
    """
    Calculates the kelvin ratio
    Based on kelvin_ratio.m by Joni Kalliokoski 2014-08-13
    :param dp: particle diameter (m)
    :param temp: temperature (K)
    :param rho: density (kg/m^3)
    :param gamma: surface tension (N/m)
    :param m_mol: molar mass (kg/mol)
    :return:
    """
    return np.exp((4 * m_mol * gamma)/(rho * gas_const * temp * dp))


def mie(m, x):
    """
    Based on Mie.m provided by TUT FYS-1320
    :param m:
    :param x:
    :return:
    """

    # The Mie.m fits this description quite nicely, so we didn't bother translating it to Python
    # https://xkcd.com/1513/

    """
    # Avoid singularity at x=0
    if x == 0:
        return [np.real(m), np.imag(m), 0, 0, 0, 0, 0, 0, 1.5]
    elif x > 0:
        nmax = round(2+x+4*(x**(1/3)))
        n1 = nmax -1

        n = np.arange(1, nmax+1)
        cn = 2*n+1
        c1n = n*(n+2) / (n+1)
        c2n = cn / n / (n+1)

        x2 = x*x

        #  TODO
    """

    return matlabeng.Mie(complex(m), float(x))[0]


def mie_abcd(m, x):
    """
    Based on mie_abcd.m provided by TUT FYS-1320
    :param m:
    :param x:
    :return:
    """

    # The mie_abcd.m fits this description quite nicely, so we didn't bother translating it to Python
    # https://xkcd.com/1513/

    """
    nmax = round(2+x+4*(x**(1/3)))
    n = np.arange(1, nmax+1)

    nu = (n+0.5)
    z = m*x
    m2 = m*m

    sqx = np.sqrt(0.5*np.pi / x)
    sqz = np.sqrt(0.5*np.pi / z)

    bx = scipy.special.jv(nu, x) * sqx
    bz = scipy.special.jv(nu, z) * sqz

    yx = scipy.special.yv(nu, x) * sqx
    hx = bx + 1j*yx

    # TODO
    """

    return matlabeng.mie_abcd(complex(m), float(x))


def p_growth_rate(temp, rho, ntot, dp, ddpdt, m_mol):
    """
    Calculates the speed of vapor pressure change
    Based on p_growth_rate.m by Joni Kalliokoski 2014-08-13
    :param temp: environment temperature (K)
    :param rho: water density (kg/m^3)
    :param ntot: particle count in a volume (#/cm^3)
    :param dp: particle diameter (m)
    :param ddpdt: speed of particle diamter change (m/s)
    :param m_mol: molar mass of water (kg/mol)
    :return: speed of partial pressure change (Pa/s)
    """
    return -(gas_const * temp * rho * ntot * np.pi * (dp**2) * ddpdt) / (2 * m_mol)


def solve_growth(t, d, m, l, ka, rho, gamma, ntot, tmax, dp0, p0):
    """
    Based on RatkaiseKasvu.m by Joni Kalliokoski 2014-08-13
    :param t: environment temperature (K)
    :param d: diffusion coefficient of water (m^2/s)
    :param m: molar mass of water (kg/mol)
    :param l: evaporation energy of water (J/kg)
    :param ka: thermal conductivity of air (W/(m*K))
    :param rho: density of water (kg/m^3)
    :param gamma: surface tension of water (N/m)
    :param ntot: particle count in a volume (#/m^3)
    :param tmax: the maximum time to compute to
    :param dp0: initial particle diameter (m)
    :param p0: initial partial pressure of water (Pa)
    :return: t (time np.array), dp (particle diameter np.array), pw (partial pressure of water np.array)
    """

    result = matlabeng.RatkaiseKasvu(float(t), float(d), float(m), float(l), float(ka), float(rho), float(gamma),
                                     float(ntot), float(tmax), float(dp0), float(p0), nargout=3)
    # array = np.array(list)
    # array = np.transpose(array)
    # return array
    t = np.array(result[0])
    dp = np.array(result[1])
    pw = np.array(result[2])

    t = t.transpose()[0]
    dp = dp.transpose()[0]
    pw = pw.transpose()[0]

    # print(t)
    # print(dp)
    # print(pw)

    return t, dp, pw


def water_pvap(temp):
    """
    Computes the saturated pressure at given temperature. Exactly the same as saturation_pressure(), but
    only with pre-defined constants for water.
    Based on water_pvap.m by Joni Kalliokoski 2014-08-13
    :param temp: temperature (K)
    :return: vapor pressure of water (Pa)
    """
    # This is the same function as on the page 4 of project instructions
    return np.exp(77.34491296-7235.424651 / temp - 8.2*np.log(temp) + 0.0057113*temp)


# ----- Our own functions -----

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


def saturation_pressure(t, a, b, c):
    """
    Computes the saturated pressure. Exactly the same as water_pvap(), but allows
    calculation for other liquids as well.
    :param t: temperature (K)
    :param a: constant from table 17-1
    :param b: constant from table 17-1
    :param c: constant from table 17-1
    :return: vapor pressure at the saturation or equilibrium condition
    """
    return 10**(a - b / (t - c))


def final_temp(t_i, p_f, p_i, gamma):
    """
    Computes the final temperature of adiabatic expansion
    :param t_i: initial temperature (K)
    :param p_f: total final pressure
    :param p_i: total initial pressure
    :param gamma: heat capacity ratio ("adiabaattivakio")
    :return: final temperature (K)
    """
    return t_i * (p_f / p_i)**((gamma-1)/gamma)


def minimum_particle_diameter(m, gamma, rho, t_f, s_r):
    """

    :param m: molar mass
    :param gamma: surface tension
    :param rho: density
    :param t_f: (final) temperature
    :param s_r: saturation ratio
    :return:
    """
    return (4 * m * gamma) / (rho * gas_const * t_f * np.log(s_r))


def minimum_particle_diameter_2(p_i, p_f, t_i, heat_capacity_ratio, a, b, c, m_mol, surface_tension, rho, satur_frac=1.0):
    """
    Calculates the minium growing particle size at adiabatic expansion.
    :param p_i: initial pressure
    :param p_f: final pressure
    :param t_i: initial temperature
    :param heat_capacity_ratio:
    :param a: contant for water in AM_ch17 table 17-1
    :param b: contant for water in AM_ch17 table 17-1
    :param c: contant for water in AM_ch17 table 17-1
    :param m_mol:  molar mass o (kg/mol)
    :param surface_tension:
    :param rho: density (kg/m^3)
    :param satur_frac: saturation fraction at initial stage, by default gas in fully saturated
    :return:
    """

    t_f = final_temp(t_i, p_f, p_i, heat_capacity_ratio)
    init_partial_press = saturation_pressure(t_i, a, b, c) * satur_frac
    final_partial_press = (p_f/p_i) * init_partial_press
    final_satur_press = saturation_pressure(t_f, a, b, c)
    """
    print("Final temperature:", t_f)
    print("Initial saturation pressure:", var_init_satur_press)
    print("Final saturation pressure:", var_final_satur_press)
    """

    return (4 * m_mol * surface_tension) / (rho * gas_const * t_f * np.log((final_partial_press) /
                                                                           (final_satur_press)))


def q_ext(dp, m, wavelength):
    """
    Ekstinktiotehokkuus
    :param dp: particle diameter (m)
    :param m:
    :param wavelength:
    :return:
    """
    vec_q_ext = np.zeros(dp.size)

    for k in range(dp.size):
        # 1e9 converts m to nm
        q_ext_value = mie(m, np.pi*dp[k]*1e9 / wavelength)
        vec_q_ext[k] = q_ext_value[3]

    return vec_q_ext


def extinction_factor(n, dp, var_q_ext):
    """

    :param n: particles in m^3
    :param dp:
    :param var_q_ext:
    :return:
    """
    return (np.pi*n*(dp**2)*var_q_ext)/4


def extinction(sigma_ext, length):
    return 1 - np.exp(-sigma_ext*length)
