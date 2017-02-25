# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

# This is the primary project file for pre-report in energy spectrum analysis

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import pyqtgraph as pg
import numpy as np
import toolbox



def set_plot_energy(plot_energy):
    """
    Plots photon and electron energy in compton scattering by function of scattered photon angle.
    :param plot_energy: pyqtgraph plot object, created by win.addPlot()
    :return:
    """
    plot_energy.addLegend()
    plot_energy.setLabel("left", "E", "eV")
    plot_energy.setLabel("bottom", "θ", "rad")

    curve_energy_photon = plot_energy.plot(pen=pg.mkPen((100, 255, 150)), name="fotoni")
    curve_energy_electron = plot_energy.plot(pen=pg.mkPen((100, 150, 255)), name="elektroni")

    theta = np.linspace(-np.pi,np.pi,1000)

    E_2f = toolbox.E_2f_fun(theta)
    E_2e = toolbox.E_2e_fun(theta)

    curve_energy_photon.setData(theta, E_2f/1.60218e-19)
    curve_energy_electron.setData(theta, E_2e/1.60218e-19)


def set_plot_probability(plot_probability, mu = 0, x0 = 0, var = 3, gamma = 1.38, x = np.linspace(-5,5,100)):
    """
    Plots Gauss and Cauchy probability distribution.
    :param plot_probability: pyqtgraph plot object, created by win.addPlot()
    :param mu: expected value for Gauss
    :param x0: center point for Couchy
    :param var: variance for Gauss
    :param gamma: deviation for Couchy
    :return:
    """
    plot_probability.addLegend()
    plot_probability.setLabel("left", "y")
    plot_probability.setLabel("bottom", "x")

    curve_probability_gauss = plot_probability.plot(pen=pg.mkPen((100, 150, 255)), name="Gauss")
    curve_probability_cauchy = plot_probability.plot(pen=pg.mkPen((100, 255, 150)), name="Cauchy")

    probability_gauss = toolbox.gauss(x,mu,var)
    probability_cachy = toolbox.cauchy(x,x0,gamma)

    curve_probability_gauss.setData(x, probability_gauss)
    curve_probability_cauchy.setData(x, probability_cachy)


def set_plot_cumulative(plot_cumulative, mu=0, x0=0, var=3, gamma=1.38, x = np.linspace(-5,5,100)):
    """
    Plots Gauss and Cauchy cumulative distribution.
    :param plot_cumulative: pyqtgraph plot object, created by win.addPlot()
    :param mu: expected value for Gauss
    :param x0: center point for Couchy
    :param var: variance for Gauss
    :param gamma: deviation for Couchy
    :return:
    """
    plot_cumulative.addLegend()
    plot_cumulative.setLabel("left", "y")
    plot_cumulative.setLabel("bottom", "x")

    curve_cumulative_gauss = plot_cumulative.plot(pen=pg.mkPen((100, 150, 255)), name="Gauss")
    curve_cumulative_cauchy = plot_cumulative.plot(pen=pg.mkPen((100, 255, 150)), name="Cauchy")

    cumulative_gauss = toolbox.gauss_cumulative(x, mu, var)
    cumulative_cauchy = toolbox.cauchy_cumulative(x, x0, gamma)

    curve_cumulative_gauss.setData(x, cumulative_gauss)
    curve_cumulative_cauchy.setData(x, cumulative_cauchy)



def main():
    pg.setConfigOptions(antialias=True, background="w", foreground="k")
    qapp = pg.mkQApp()

    win = pg.GraphicsWindow(title="Pre-report")

    plot_energy = win.addPlot()#title="Energy by angle")
    set_plot_energy(plot_energy)

    mu, x0 = 0, 0
    var, gamma = 0.394, 0.5
    x = np.linspace(-2.5, 2.5, 1000)

    plot_probability = win.addPlot()#title="Probability distribution")
    set_plot_probability(plot_probability, mu, x0, var, gamma, x)

    win.nextRow()

    plot_cumulative = win.addPlot()#title="Cumulative distribution")
    set_plot_cumulative(plot_cumulative, mu, x0, var, gamma, x)

    win.resize(950,800)

    # PyQtGraph main loop
    qapp.exec_()



main()