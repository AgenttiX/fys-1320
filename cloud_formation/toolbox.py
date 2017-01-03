# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017
# Strongly based on Tampere University of Technology course FYS-1320 example code

# Must be imported first to prevent errors
import matlab.engine

import numpy as np
import sympy
import sympy.solvers
import scipy.optimize
import scipy.special

# Start MATLAB
matlabeng = matlab.engine.start_matlab()

R = 8.3144621


def dp_growth_rate(dp, T, rho, gamma, M, D, p, Td):
    """
    Based on dp_growth_rate.m by Joni Kalliokoski 2014-08-13
    :param dp: particle diameter (m)
    :param T: environment temperature (K)
    :param rho: water density (kg / m^3)
    :param gamma: water surface tension (N/m)
    :param M: water molar mass (kg/mol)
    :param D: water diffusion coefficient (m^2 / s)
    :param p: environment water partial pressure (Pa)
    :param Td: particle surface temperature (K)
    :return: speed of particle diameter change (m/s)
    """
    pd = water_pvap(Td) * kelvin_ratio(dp, Td, rho, gamma, M)
    lambda_var = free_path(T, D, M)

    return ((4*D*M)/(R*rho*dp))*((p/T)-(pd/Td))*fuchs_sutugin(dp, lambda_var)


def droplet_temp(T, D, M, L, ka, p, dp, rho, gamma):
    """
    Based on droplet_temp.m by Joni Kalliokoski 2014-08-13
    :param T: environment temperature (K)
    :param D: water diffusion coefficient (m^2 / s)
    :param M: water molar mass (kg/mol)
    :param L: water evaporation energy (J/kg)
    :param ka: thermal conductivity of air (W/(m*K))
    :param p: environment water partial pressure (Pa)
    :param dp: particle diamter (m)
    :param rho: water density (kg/m^3)
    :param gamma: water surface tension (N/m)
    :return: particle surface temperature (K)
    """
    lambda_var = free_path(T, D, M)
    # x = sympy.var("x")
    # return sympy.solvers.solve(x - T - ((D*M*L)/(R*ka))* ( (p/T)- (water_pvap( x)*kelvin_ratio( dp,x,rho,gamma,M )/x) ) * FuchsSutugin( dp, lambda_var ), x)

    def droplet_temp_support(x):
        x - T - ((D * M * L) / (R * ka)) * ((p / T) - (water_pvap(x) * kelvin_ratio(dp, x, rho, gamma, M) / x)) * fuchs_sutugin(dp, lambda_var)

    return scipy.optimize.fsolve(droplet_temp_support, T)


def dy_growth(t, y, T, D, M, L, ka, rho, gamma, Ntot, dp0):
    """
    Based on dyKasvu.m by Joni Kalliokoski 2014-08-13
    :param t: time (s)
    :param y: y[0] = particle diameter (m), y[1] = environment water partial pressure (Pa)
    :param T: environment temperature (K)
    :param D: water diffusion coefficient (m^2/s)
    :param M: water molar mass (kg/mol)
    :param L: water evaporation energy (J/kg)
    :param ka: thermal conductivity of air (W/(m*K))
    :param rho: water density (kg/m^3)
    :param gamma: water surface tension (N/m)
    :param Ntot: particle count in a volume (#/cm^3)
    :param dp0: initial particle diameter (m)
    :return: speeds of change: [0] = particle diameter (m), [1] = environment water partial pressure (Pa)
    """
    dy = np.zeros(2)
    Td = droplet_temp(T, D, M, L, ka, y[1], y[0], rho, gamma)
    dy[0] = dp_growth_rate(y[0], T, rho, gamma, M, D, y[1], Td)
    dy[1] = p_growth_rate(T, rho, Ntot, y[0], dy[0], M)

    if (dy[0] < 0) and (y[0] <= dp0):
        dy[0] = 0
        dy[1] = 0

    return dy


def free_path(T, D, M):
    """
    Calculates the free path for a molecule
    Based on free_path.m by Joni Kalliokoski 2014-08-13
    :param T: temperature (K)
    :param D: diffusion coefficient (m^2/s)
    :param M: molar mass (kg/mol)
    :return: free path (m)
    """
    return 3*D*np.sqrt((np.pi*M)/(8*R*T))


def fuchs_sutugin(dp, lambda_var):
    """
    Calculates the Fuchs-Sutugin correction factor
    Based on FuchsSutugin.m by Joni Kalliokoski 2014-08-13
    :param dp:
    :param lambda_var: free path (m)
    :return: Fuchs-Sutugin correction factor
    """
    Kn = 2*lambda_var / dp
    return (1 + Kn) / (1 + 1.71*Kn + 1.33*(Kn**2))


def kelvin_ratio(dp, T, rho, gamma, M):
    """
    Calculates the kelvin ratio
    Based on kelvin_ratio.m by Joni Kalliokoski 2014-08-13
    :param dp: particle diameter (m)
    :param T: temperature (K)
    :param rho: density (kg/m^3)
    :param gamma: surface tension (N/m)
    :param M: molar mass (kg/mol)
    :return:
    """
    kr = np.exp((4*M*gamma)/(rho*R*T*dp))


def mie(m, x):
    """
    Based on Mie.m provided by TUT FYS-1320
    :param m:
    :param x:
    :return:
    """

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

    return matlabeng.Mie(m, x)


def mie_abcd(m, x):
    """
    Based on mie_abcd.m provided by TUT FYS-1320
    :param m:
    :param x:
    :return:
    """

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

    return matlabeng.mie_abcd(m, x)


def p_growth_rate(T, rho, Ntot, dp, ddpdt, M):
    """
    Calculates the speed of vapor pressure change
    Based on p_growth_rate.m by Joni Kalliokoski 2014-08-13
    :param T: environment temperature (K)
    :param rho: water density (kg/m^3)
    :param Ntot: particle count in a volume (#/cm^3)
    :param dp: particle diameter (m)
    :param ddpdt: speed of particle diamter change (m/s)
    :param M: molar mass of water (kg/mol)
    :return: speed of partial pressure change (Pa/s)
    """
    return -(R*T*rho*Ntot*np.pi*(dp**2)*ddpdt) / (2*M)


def solve_growth(T, D, M, L, ka, rho, gamma, Ntot, tmax, dp0, p0):
    """
    Based on RatkaiseKasvu.m by Joni Kalliokoski 2014-08-13
    :param T: environment temperature (K)
    :param D: diffusion coefficient of water (m^2/s)
    :param M: molar mass of water (kg/mol)
    :param L: evaporation energy of water (J/kg)
    :param ka: thermal conductivity of air (W/(m*K))
    :param rho: density of water (kg/m^3)
    :param gamma: surface tension of water (N/m)
    :param Ntot: particle count in a volume (#/cm^3)
    :param tmax: the maximum time to compute to
    :param dp0: initial particle diameter (m)
    :param p0: initial partial pressure of water (Pa)
    :return:
    """

    result = matlabeng.RatkaiseKasvu(float(T), float(D), float(M), float(L), float(ka), float(rho), float(gamma), float(Ntot), float(tmax), float(dp0), float(p0), nargout=3)
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


def water_pvap(T):
    """
    Based on water_pvap.m by Joni Kalliokoski 2014-08-13
    :param T: temperature (K)
    :return: vapor pressure of water (Pa)
    """
    return np.exp(77.34491296-7235.424651/T-8.2*np.log(T)+0.0057113*T)
