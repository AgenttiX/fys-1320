# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

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


# npTDMS is a library for reading TDMS files created by LabVIEW
# https://pypi.python.org/pypi/npTDMS/
import toolbox_2
import nptdms
# PyQtGraph can be found from Ubuntu repositories as python3-pyqtgraph
# http://www.pyqtgraph.org/
import pyqtgraph as pg
import numpy as np

from pyqtgraph.Qt import QtGui



class Measurement:
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

        # Read data
        self.meas_pressure = self.read_to_list("data/by_pressure/", meas_first, meas_last_pressure)
        self.meas_mixture = self.read_to_list("data/by_mixture/", meas_first, meas_last_mixture)

        # Control window
        widget2 = QtGui.QWidget()
        widget2.setWindowTitle("Cloud formation data analysis")

        layout = QtGui.QGridLayout()
        widget2.setLayout(layout)

        labels = ["Measurement number", "Measurement series", "Data type (1:p_diff, 2:p_abs, 3:ext)",
                  "Noice reduction (0 = off)", "", "Particle size (nm) (useless)",
                  "Particle density (n*1e-10 1/m^3)"]

        for i, text in enumerate(labels):
            label = QtGui.QLabel()
            label.setText(text)
            layout.addWidget(label, i, 0)

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

        self.plot_select = win.addPlot(title="Region selection")
        self.plot_zoom = win.addPlot(title="Zoom on selected region")
        self.plot_simulate = win.addPlot(title="Simulation on particle growth")

        self.linear_region = pg.LinearRegionItem([100000,110000])
        self.linear_region.setZValue(-10)
        self.plot_select.addItem(self.linear_region)

        self.curve_select = None
        self.curve_zoom = None
        self.curve_simulate = None
        self.update_change()

        self.linear_region.sigRegionChanged.connect(self.update_zoom_plot)
        self.plot_zoom.sigXRangeChanged.connect(self.update_zoom_region)
        self.plot_zoom.sigRangeChanged.connect(self.update_simulate_plot)

        win.resize(1100,600)
        self.update_zoom_plot()

        # PyQtGraph main loop
        qapp.exec_()

    def update_change(self):
        self.meas_selected_number = self.__input_meas.value()
        self.meas_selected_series = self.__input_series.value()
        self.selected_data = self.__input_select.value()
        self.noice_reduction_number = self.__input_noice.value()
        self.particle_size_number = self.__input_p_size.value()
        self.particle_density_number = self.__input_p_density.value()
        self.saturation_percentage = self.__input_saturation_percentage.value()

        # print("Selected series:", self.meas_selected_series)
        # print("Selected measurement:", self.meas_selected_number)

        if not 1 <= self.meas_selected_number <= 17:
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
            measurement = self.meas_pressure[meas_index]
        elif self.meas_selected_series == 2:
            measurement = self.meas_mixture[meas_index]
        else:
            raise ValueError

        if self.selected_data == 1:
            data = toolbox_2.remove_noice(measurement.p_diff, self.noice_reduction_number)
        elif self.selected_data == 2:
            data = toolbox_2.remove_noice(measurement.p_abs, self.noice_reduction_number)
        elif self.selected_data == 3:
            data = toolbox_2.remove_noice(toolbox_2.flip_and_normalize(measurement.ext), self.noice_reduction_number)
        else:
            raise ValueError

        index_of_drop = toolbox_2.find_drop_index(measurement.p_diff)
        time = measurement.time - measurement.time[index_of_drop]   # time vector starts from the beginning of pressure drop

        if self.first_update:
            self.curve_select = self.plot_select.plot(time, data)
            self.curve_zoom = self.plot_zoom.plot(measurement.time, data)
            self.curve_simulate = self.plot_simulate.plot()
            self.first_update = False
        else:
            self.curve_select.setData(time, data)
            self.curve_zoom.setData(time, data)

        # Data simulation
        if self.simulate_bool:
            t_max = 3
            p_i, p_f = toolbox_2.get_pressure_change(measurement)
            size, time2, smallest_growing_particle = toolbox_2.simulate_extinction(self.particle_size_number * 1e-9,
                                                                                   p_i, p_f,
                                                                                   self.particle_density_number * 1e10,
                                                                                   t_max, self.saturation_percentage/100)
            # short print:
            # print("M:", self.meas_selected_number, ", ", round((p_i - p_f) / 1000, 3), "kPa", ", ", self.saturation_percentage, "%", ", ", round(smallest_growing_particle * 1e9, 2), "nm", ", ", sep="")

            if (smallest_growing_particle > 0):
                print("M:", self.meas_selected_number, " S:", self.meas_selected_series, " D:", self.selected_data,
                      ", smallest growing particle for pressure change (", round(p_i / 1000, 2), "-",
                      round(p_f / 1000, 2), " = ", round( (p_i - p_f) / 1000, 2), "kPa) in ", self.saturation_percentage,
                      "% humidity is ", round(smallest_growing_particle * 1e9, 2), "nm", sep="")
            else:
                print("M:", self.meas_selected_number, " S:", self.meas_selected_series, " D:", self.selected_data,
                      ", no particle will grow in ", "(", round(p_i / 1000, 2), "-", round(p_f / 1000, 2), " = ",
                      round((p_i - p_f) / 1000, 2), "kPa)", " pressure change and ", self.saturation_percentage, "% humidity ", sep="")

            self.curve_simulate.setData(time2, size)
            self.simulate_bool = False

        # set the graphs to the point of pressure drop, units are in seconds
        self.plot_select.setXRange(-2, 4, padding=0)
        self.plot_zoom.setXRange(0, 0.5, padding=0)

        self.update_zoom_region()
        #self.update_zoom_plot()

    def update_zoom_plot(self):
        self.plot_zoom.setXRange(*self.linear_region.getRegion(), padding=0)

    def update_zoom_region(self):
        self.linear_region.setRegion(self.plot_zoom.getViewBox().viewRange()[0])

    def update_simulate_plot(self):
        a = self.plot_zoom.getViewBox().viewRange()
        self.plot_simulate.setXRange(a[0][0], a[0][1])
        self.plot_simulate.setYRange(a[1][0], a[1][1])

    def simulate_button_clicked(self):
        self.simulate_bool = True
        self.update_change()

    def read_to_list(self, folder, start, stop):
        measurements = []
        for i in range(start, stop+1):
            measurements.append(Measurement(folder + str(i) + ".tdms"))
        return measurements






def main():
    Main()


main()
