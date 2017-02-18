# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

## @package data_analysis
# This program is for analysing data created by logger_software v2 of the cloud formation experiment
# of the Tampere University of Technology student laboratory of physics


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

# toolbox_2 is a set of supplementary functions that have been developed together with this program
import toolbox_2

# npTDMS is a library for reading TDMS files created by LabVIEW
# https://pypi.python.org/pypi/npTDMS/
import nptdms

# PyQtGraph can be found from Ubuntu repositories as python3-pyqtgraph
# http://www.pyqtgraph.org/
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

import numpy as np


class Measurement:
    """ This class represents a single measurement and its data """
    def __init__(self, path):
        # Read the TDMS file
        tdms_file = nptdms.TdmsFile(path)
        # The group and channel names have been reverse engineered using a hex editor
        tdms_p_diff = tdms_file.object("GroupName", "Paine_ero_kPa")
        tdms_p_abs = tdms_file.object("GroupName", "Absoluuttipaine_kPa")
        tdms_ext = tdms_file.object("GroupName", "Light extinction_V")

        # Extract data from the objects
        # The constants have been reverse engineered from convert_TDMS_to_ASCII.exe with file_analysis.py
        self.p_diff = tdms_p_diff.data * 1.35322218e-03 + 2.20953580e+01
        self.p_abs = tdms_p_abs.data * 1.35312729e-02 + 2.35936196e+02
        self.ext = tdms_ext.data * 6.36900701e-04 + 1.06080729e+01
        self.time = tdms_p_diff.time_track()

        if not (self.p_diff.size == self.p_abs.size == self.ext.size == self.time.size):
            print("WARNING: this file is malformed:", path)


class Main:
    def __init__(self):
        # pg.setConfigOptions(antialias=True, background="w", foreground="k")
        qapp = pg.mkQApp()

        # File parameters
        meas_first = 1
        meas_last_pressure = 17
        meas_last_mixture = 8
        self.meas_selected_number = 1
        self.meas_selected_series = 1
        self.selected_data = 1
        self.noice_reduction_number = 0
        self.simulate_bool = False
        self.particle_size_number = 50
        self.particle_density_number = 10
        self.saturation_percentage = 100
        self.first_update = True

        self.measurement = None
        self.data = None
        self.index_of_drop = None

        self.smallest_particles = np.load("smallest_particles.npy")
        self.number_counts = np.load("number_counts.npy")
        self.number_counts_2 = np.load("number_counts_2.npy")

        self.particle_distribution_x = np.zeros(10)  # particle size (the plot in window)
        self.particle_distribution_y = np.zeros(10)  # number count for that size

        # Read data
        self.meas_pressure = self.read_to_list("data/by_pressure/", meas_first, meas_last_pressure)
        self.meas_mixture = self.read_to_list("data/by_mixture/", meas_first, meas_last_mixture)

        # Control window
        widget2 = QtGui.QWidget()
        widget2.setWindowTitle("Cloud formation data analysis")

        layout = QtGui.QGridLayout()
        widget2.setLayout(layout)

        # Text labels for the control window
        labels = ["Measurement number", "Measurement series", "Data type (1:p_diff, 2:p_abs, 3:ext)",
                  "Noice reduction (0 = off)", "", "Particle size (nm) (useless)",
                  "Particle density (n*1e10 1/m^3)", "Initial saturation"]

        for i, text in enumerate(labels):
            label = QtGui.QLabel()
            label.setText(text)
            layout.addWidget(label, i, 0)

        # Inputs for the control window
        self.__input_meas = pg.SpinBox(value=self.meas_selected_number, int=True, minStep=1, step=1)
        self.__input_meas.editingFinished.connect(self.update_change)
        layout.addWidget(self.__input_meas, 0, 1)

        self.__input_series = pg.SpinBox(value=self.meas_selected_series, int=True, minStep=1, step=1)
        self.__input_series.editingFinished.connect(self.update_change)
        layout.addWidget(self.__input_series, 1, 1)

        self.__input_select = pg.SpinBox(value=self.selected_data, int=True, minStep=1, step=1, min=1)
        self.__input_select.editingFinished.connect(self.update_change)
        layout.addWidget(self.__input_select, 2, 1)

        self.__input_noice = pg.SpinBox(value=self.noice_reduction_number, int=True, minStep=1, step=1, min=0)
        self.__input_noice.editingFinished.connect(self.update_change)
        layout.addWidget(self.__input_noice, 3, 1)

        self.__updateButton = QtGui.QPushButton("Update")
        self.__updateButton.clicked.connect(self.update_change)
        layout.addWidget(self.__updateButton, 4, 1)

        self.__input_p_size = pg.SpinBox(value=self.particle_size_number, int=True, minStep=1, step=1, min=0)
        self.__input_p_size.editingFinished.connect(self.update_change)
        layout.addWidget(self.__input_p_size, 5, 1)

        self.__input_p_density = pg.SpinBox(value=self.particle_density_number, int=True, minStep=1, step=1, min=0)
        self.__input_p_density.editingFinished.connect(self.update_change)
        layout.addWidget(self.__input_p_density, 6, 1)

        self.__input_saturation_percentage = pg.SpinBox(value=self.saturation_percentage, int=True, minStep=1, step=1, min=0, max=100)
        self.__input_saturation_percentage.editingFinished.connect(self.update_change)
        layout.addWidget(self.__input_saturation_percentage, 7, 1)

        self.__simulateButton = QtGui.QPushButton("Simulate")
        self.__simulateButton.clicked.connect(self.simulate_button_clicked)
        layout.addWidget(self.__simulateButton, 8, 1)

        widget2.show()

        # Graph window
        win = pg.GraphicsWindow(title="Cloud formation data analysis")

        self.plot_select = win.addPlot()#title="Region selection")
        self.plot_zoom = win.addPlot()#title="Zoom on selected region")
        #self.plot_simulate = win.addPlot()#title="Simulation on particle growth")
        win.nextRow()
        self.plot_distribution = win.addPlot()#title="Concentration distribtion")
        self.plot_rotatometer = win.addPlot()#title="Rotatometer_fractions")

        self.plot_distribution.addLegend()
        self.plot_rotatometer.addLegend()

        self.linear_region = pg.LinearRegionItem([100000, 110000])
        self.linear_region.setZValue(-10)
        self.plot_select.addItem(self.linear_region)

        self.line = pg.InfiniteLine(angle=90, movable=True)
        self.plot_zoom.addItem(self.line)

        self.set_labels()

        self.curve_select = None
        self.curve_zoom = None
        self.curve_simulate = None
        self.curve_distribution = None
        self.curve_distribution_cumulative = None
        self.curve_rotatometer = None
        self.curve_rotatometer_fit = None
        self.update_change()

        self.linear_region.sigRegionChanged.connect(self.update_zoom_plot)
        self.plot_zoom.sigXRangeChanged.connect(self.update_zoom_region)
        #self.plot_zoom.sigRangeChanged.connect(self.update_simulate_plot)

        self.line.sigPositionChangeFinished.connect(self.line_moved)

        win.resize(940,840)
        self.update_zoom_plot()

        # PyQtGraph main loop
        qapp.exec_()

    def update_change(self):
        """
        Updates the whole view based on the values of the configuration window.
        :return: -
        """
        self.meas_selected_number = self.__input_meas.value()
        self.meas_selected_series = self.__input_series.value()
        self.selected_data = self.__input_select.value()
        self.noice_reduction_number = self.__input_noice.value()
        self.particle_size_number = self.__input_p_size.value()
        self.particle_density_number = self.__input_p_density.value()
        self.saturation_percentage = self.__input_saturation_percentage.value()

        # print("Selected series:", self.meas_selected_series)
        # print("Selected measurement:", self.meas_selected_number)

        if not 1 <= self.meas_selected_number <= 17:  # In series 2 there's measurement 0, but that's a copy of the 8th
            raise ValueError
        if not 0 <= self.noice_reduction_number:
            raise ValueError
        if not 0 < self.particle_size_number:
            raise ValueError
        if not 0 <= self.particle_density_number:
            raise ValueError
        if not 0 <= self.saturation_percentage <= 100:
            raise ValueError

        meas_index = self.meas_selected_number - 1
        if self.meas_selected_series == 1:
            self.measurement = self.meas_pressure[meas_index]
        elif self.meas_selected_series == 2:
            self.measurement = self.meas_mixture[meas_index]
        else:
            raise ValueError

        if self.selected_data == 1:
            self.data = toolbox_2.remove_noise(self.measurement.p_diff, self.noice_reduction_number)
        elif self.selected_data == 2:
            self.data = toolbox_2.remove_noise(self.measurement.p_abs, self.noice_reduction_number)
        elif self.selected_data == 3:
            self.data = toolbox_2.remove_noise(toolbox_2.flip_and_normalize(self.measurement.ext), self.noice_reduction_number)
            self.plot_zoom.setYRange(-0.1, 1, padding=0)
        elif self.selected_data == 4:
            self.data = toolbox_2.remove_noise(self.measurement.ext, self.noice_reduction_number)
        else:
            raise ValueError

        self.index_of_drop = toolbox_2.find_drop_index(self.measurement.p_diff)
        time = self.measurement.time - self.measurement.time[self.index_of_drop]   # time vector starts from the beginning of pressure drop

        self.update_distribution()

        if self.first_update:
            self.curve_select = self.plot_select.plot(time, self.data)
            self.curve_zoom = self.plot_zoom.plot(self.measurement.time, self.data)
            self.curve_simulate = self.plot_zoom.plot(pen=pg.mkPen((100, 255, 200)))

            self.curve_distribution_cumulative = self.plot_distribution.plot(pen=pg.mkPen((100, 200, 255)), name="´   Kumulatiivinen pitoisuus", symbolBrush=(80, 160, 201), symbolPen='w')
            self.curve_distribution = self.plot_distribution.plot(name="Pitoisuusjakauma", symbolBrush=(50, 50, 255), symbolPen='w')
            self.curve_rotatometer = self.plot_rotatometer.plot(name="Mitattu pitoisuus", symbolBrush=(50, 50, 255), symbolPen='w')
            self.curve_rotatometer_fit = self.plot_rotatometer.plot(pen=pg.mkPen((100, 255, 200)), name="´  Ideaalinen pitoisuus")

            self.first_update = False
        else:
            self.curve_select.setData(time, self.data)
            self.curve_zoom.setData(time, self.data)

            self.curve_distribution.setData(self.particle_distribution_x, self.particle_distribution_y*1e-10)
            self.curve_distribution_cumulative.setData(self.smallest_particles, self.number_counts*1e-10)
            self.curve_rotatometer.setData(np.array([4, 6, 8, 10, 12, 14, 16, 18]), self.number_counts_2*1e-10)
            x = np.linspace(3.5, 20, 100)
            self.curve_rotatometer_fit.setData(x, self.number_counts_2[0]*4*(1/x)*1e-10)

        if self.simulate_bool:
            self.simulation()

        self.set_labels()

        # set the graphs to the point of pressure drop, units are in seconds
        self.plot_select.setXRange(-2, 4, padding=0)
        self.plot_zoom.setXRange(0, 0.3, padding=0)
        self.line.setX(0.1)

        self.update_zoom_region()
        # self.update_zoom_plot()

    def update_zoom_plot(self):
        """
        Updates the zoomed plot to match the region set in the primary graph
        :return: -
        """
        self.plot_zoom.setXRange(*self.linear_region.getRegion(), padding=0)

    def update_zoom_region(self):
        """
        Updates the zoom region in the primary graph to match the area of the zoomed graph
        :return: -
        """
        self.linear_region.setRegion(self.plot_zoom.getViewBox().viewRange()[0])

    def update_simulate_plot(self):
        """
        Updates the simulation plot view to match that of the zoomed signal
        :return:
        """
        a = self.plot_zoom.getViewBox().viewRange()
        self.plot_simulate.setXRange(a[0][0], a[0][1])
        self.plot_simulate.setYRange(a[1][0], a[1][1])

    def simulate_button_clicked(self):
        """
        Enable simulating and update the graphs.
        This function is called by the simulation button.
        :return: -
        """
        self.simulate_bool = True
        self.update_change()

    def line_moved(self):
        """
        This function is called when the line is moved in the zoomed graph.
        :return: -
        """

        # The line is supposed to be moved by hand to the beginning of first wrinkle.
        # The optimal spot is local maximum (not always visible)
        ext_index = self.index_of_drop + int(self.line.value() * 10000)
        ext_value = self.data[ext_index]

        p_i, p_f = toolbox_2.get_pressure_change(self.measurement)
        smallest_growing_particle = toolbox_2.minimum_particle_diameter(p_i, p_f, self.saturation_percentage / 100)

        n = toolbox_2.particle_count_2(ext_value)

        # measurement series 1
        if self.selected_data == 3 and 7 <= self.meas_selected_number <= 17 and self.meas_selected_series == 1:
            index = self.meas_selected_number - 7   # Assumes that first measurement is number 7
            self.smallest_particles[index] = smallest_growing_particle
            self.number_counts[index] = n

            self.update_distribution()
            # Update plot
            self.curve_distribution.setData(self.particle_distribution_x, self.particle_distribution_y*1e-10)
            self.curve_distribution_cumulative.setData(self.smallest_particles, self.number_counts*1e-10)

        # measurement series 2
        elif self.selected_data == 3 and self.meas_selected_series == 2:
            index = self.meas_selected_number - 1       # begins from 1, 0th measurement is just copy of 8th
            self.number_counts_2[index] = n

            self.curve_rotatometer.setData(np.array([4, 6, 8, 10, 12, 14, 16, 18]), self.number_counts_2*1e-10)
            x = np.linspace(3.5, 20, 100)
            self.curve_rotatometer_fit.setData(x, self.number_counts_2[0] * 4 * (1 / x) *1e-10)

        print(smallest_growing_particle, p_i - p_f)

    @staticmethod
    def read_to_list(folder, start, stop):
        """
        Read the measurements of a folder to Measurement objects in a list.
        :param folder: string of the folder to be scanned
        :param start: start index of the files
        :param stop: stop index of the files
        :return: list of Measurement objects
        """
        measurements = []
        for i in range(start, stop+1):
            measurements.append(Measurement(folder + str(i) + ".tdms"))
        return measurements

    def simulation(self):
        """
        Computes and displays the simulation of extinction as a function of time.
        The simulation is done using mathematical tools provided by toolbox_2.py.
        The toolbox_2.py utilises toolbox.py which in turn uses Matlab to run
        code provided by Tampere University of Technology.
        :return: -
        """

        t_max = 3
        if self.meas_selected_series == 1:
            particle_density_number = self.particle_density_number
        else: # series 2:
            factors = 4/np.array([4, 6, 8, 10, 12, 14, 16, 18])
            factor = factors[(self.meas_selected_number-1)]
            particle_density_number = self.particle_density_number * factor

        p_i, p_f = toolbox_2.get_pressure_change(self.measurement)
        size, time2 = toolbox_2.simulate_extinction(self.particle_size_number * 1e-9,
                                                    p_i, p_f,
                                                    particle_density_number * 1e10,
                                                    t_max, self.saturation_percentage / 100)
        smallest_growing_particle = toolbox_2.minimum_particle_diameter(p_i, p_f, self.saturation_percentage / 100)
        # short print:
        # print("M:", self.meas_selected_number, ", ", round((p_i - p_f) / 1000, 3), "kPa", ", ", self.saturation_percentage, "%", ", ", round(smallest_growing_particle * 1e9, 2), "nm", ", ", sep="")

        if smallest_growing_particle > 0:
            print("M:", self.meas_selected_number, " S:", self.meas_selected_series, " D:", self.selected_data,
                  ", smallest growing particle for pressure change (", round(p_i / 1000, 2), "-",
                  round(p_f / 1000, 2), " = ", round((p_i - p_f) / 1000, 2), "kPa) in ", self.saturation_percentage,
                  "% humidity is ", round(smallest_growing_particle * 1e9, 2), "nm", sep="")
        else:
            print("M:", self.meas_selected_number, " S:", self.meas_selected_series, " D:", self.selected_data,
                  ", no particle will grow in ", "(", round(p_i / 1000, 2), "-", round(p_f / 1000, 2), " = ",
                  round((p_i - p_f) / 1000, 2), "kPa)", " pressure change and ", self.saturation_percentage,
                  "% humidity ", sep="")

        self.curve_simulate.setData(time2+0.05, size)
        self.simulate_bool = False

    def update_distribution(self):
        #  sort arrays by minimum particle diameter
        inds = self.smallest_particles.argsort()
        sorted_smallest_particles = self.smallest_particles[inds]
        sorted_number_counts = self.number_counts[inds]
        for i in range(10):  # 11 measurements, iterates 10 in between gaps

            d_first = sorted_smallest_particles[i]
            d_second = sorted_smallest_particles[i+1]
            dn = sorted_number_counts[i+1] - sorted_number_counts[i]
            dd = d_second - d_first

            # takes in account that the intervals are different size, so the value is arbitarily number count per >>> NANOMETER <<<
            self.particle_distribution_y[i] = (dn / (dd * 1e9))
            self.particle_distribution_x[i] = (d_first + dd / 2)

            # if (i in range(10)):
            #    print(inds[i]+7, " dn", "%.2e"%dn, " dd", "%.2e"%dd, " d_first", "%.2e"%d_first, " N", "%.2e"%sorted_number_counts[i - 1])

    def save_numpy_array(self):
        """
        Saves the changes made to the result arrays.
        :return: -
        """
        np.save("smallest_particles.npy", self.smallest_particles)
        np.save("number_counts.npy", self.number_counts)
        np.save("number_counts_2.npy", self.number_counts_2)

    def set_labels(self):
        """
        Sets the graph labels to match the current data.
        :return: -
        """

        if 1 <= self.selected_data <= 2:
            self.plot_select.setLabel("left", "P (kPa)")
            self.plot_select.setLabel("bottom", "t", "s")
            self.plot_zoom.setLabel("left", "P (kPa)")
            self.plot_zoom.setLabel("bottom", "t", "s")

        elif self.selected_data == 3:
            self.plot_select.setLabel("left", "ext", "")
            self.plot_select.setLabel("bottom", "t", "s")
            self.plot_zoom.setLabel("left", "ext", "")
            self.plot_zoom.setLabel("bottom", "t", "s")

        elif self.selected_data == 4:
            self.plot_select.setLabel("left", "U", "V")
            self.plot_select.setLabel("bottom", "t", "s")
            self.plot_zoom.setLabel("left", "U", "V")
            self.plot_zoom.setLabel("bottom", "t", "s")

        #self.plot_simulate.setLabel("left", "ext", "")
        #self.plot_simulate.setLabel("bottom", "t", "s")

        self.plot_distribution.setLabel("left", "N ×10¹⁰")
        self.plot_distribution.setLabel("bottom", "d_p", "m")
        self.plot_distribution.showGrid(y=True)

        self.plot_rotatometer.setLabel("left", "N ×10¹⁰")
        self.plot_rotatometer.setLabel("bottom", "laimennusvirtaus")
        self.plot_rotatometer.showGrid(y=True)


def main():
    Main()


main()
