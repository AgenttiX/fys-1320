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
    plot_energy.setLabel("left", "E", "J")
    plot_energy.setLabel("bottom", "theta", "rad")

    curve_energy_photon = plot_energy.plot(pen=pg.mkPen((100, 255, 150)), name="photon")
    curve_energy_electron = plot_energy.plot(pen=pg.mkPen((100, 150, 255)), name="electron")

    theta = np.linspace(-np.pi,np.pi,1000)

    E_2f = toolbox.E_2f_fun(theta)
    E_2e = toolbox.E_2e_fun(theta)

    curve_energy_photon.setData(theta, E_2f)
    curve_energy_electron.setData(theta, E_2e)


def set_plot_gauss(plot_gauss):
    """
    Plots Gauss probability distribution and cumulative distribution.
    :param plot_gauss: pyqtgraph plot object, created by win.addPlot()
    :return:
    """
    plot_gauss.addLegend()
    plot_gauss.setLabel("left", "y")
    plot_gauss.setLabel("bottom", "x")

    curve_probability = plot_gauss.plot(pen=pg.mkPen((100, 255, 150)), name="Probability distribution")
    curve_cumulative = plot_gauss.plot(pen=pg.mkPen((100, 150, 255)), name="Cumulative distribution")

    x = np.linspace(-5,5,100)
    mu = 0
    var = 3

    prob = toolbox.gauss(x,mu,var)
    cumul = toolbox.gauss_cumulative(x, mu, var)

    curve_probability.setData(x, prob)
    curve_cumulative.setData(x, cumul)


def set_plot_cauchy(plot_cauchy):
    """
    Plots Cauchy probability distribution and cumulative distribution.
    :param plot_cauchy: pyqtgraph plot object, created by win.addPlot()
    :return:
    """

    plot_cauchy.addLegend()
    plot_cauchy.setLabel("left", "y")
    plot_cauchy.setLabel("bottom", "x")

    curve_probability = plot_cauchy.plot(pen=pg.mkPen((100, 255, 150)), name="Probability distribution")
    curve_cumulative = plot_cauchy.plot(pen=pg.mkPen((100, 150, 255)), name="Cumulative distribution")

    x = np.linspace(-5,5,100)
    x0 = 0
    gamma = 3

    prob = toolbox.cauchy(x,x0,gamma)
    cumul = toolbox.cauchy_cumulative(x, x0, gamma)

    curve_probability.setData(x, prob)
    curve_cumulative.setData(x, cumul)



def main():
    pg.setConfigOptions(antialias=True, background="w", foreground="k")
    qapp = pg.mkQApp()

    win = pg.GraphicsWindow(title="Pre-report")

    plot_energy = win.addPlot(title="Energy by angle")
    set_plot_energy(plot_energy)

    plot_gauss = win.addPlot(title="Gauss distribution")
    set_plot_gauss(plot_gauss)

    win.nextRow()

    plot_cauchy = win.addPlot(title="Cauchy distribution")
    set_plot_cauchy(plot_cauchy)

    # PyQtGraph main loop
    qapp.exec_()



main()