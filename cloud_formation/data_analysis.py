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
import nptdms

# PyQtGraph can be found from Ubuntu repositories as python3-pyqtgraph
# http://www.pyqtgraph.org/
import pyqtgraph as pg
import numpy as np

from pyqtgraph.Qt import QtGui
from math import factorial


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
                  "Noice reduction (0 = off)"]

        for i, text in enumerate(labels):
            label = QtGui.QLabel()
            label.setText(text)
            layout.addWidget(label, i, 0)

        self.__input_meas = pg.SpinBox(value=self.meas_selected_number, int=True, minStep=1, step=1)
        layout.addWidget(self.__input_meas, 0, 1)

        self.__input_series = pg.SpinBox(value=self.meas_selected_series, int=True, minStep=1, step=1)
        layout.addWidget(self.__input_series, 1, 1)

        self.__input_select = pg.SpinBox(value=self.selected_data, int=True, minStep=1, step=1, min=1)
        layout.addWidget(self.__input_select, 2, 1)

        self.__input_noice = pg.SpinBox(value=self.noice_reduction_number, int=True, minStep=1, step=1, min=0)
        layout.addWidget(self.__input_noice, 3, 1)

        self.__updateButton = QtGui.QPushButton("Update")
        layout.addWidget(self.__updateButton, 4, 1)

        self.__updateButton.clicked.connect(self.update_change)

        widget2.show()

        # Graph window
        win = pg.GraphicsWindow(title="Cloud formation data analysis")

        self.plot_select = win.addPlot(title="Region selection")
        self.plot_zoom = win.addPlot(title="Zoom on selected region")

        self.linear_region = pg.LinearRegionItem([100000,110000])
        self.linear_region.setZValue(-10)
        self.plot_select.addItem(self.linear_region)

        self.curve_select = None
        self.curve_zoom = None
        self.update_change()

        self.linear_region.sigRegionChanged.connect(self.update_zoom_plot)
        self.plot_zoom.sigXRangeChanged.connect(self.update_zoom_region)
        self.update_zoom_plot()

        # PyQtGraph main loop
        qapp.exec_()

    def update_change(self):
        self.meas_selected_number = self.__input_meas.value()
        self.meas_selected_series = self.__input_series.value()
        self.selected_data = self.__input_select.value()
        self.noice_reduction_number = self.__input_noice.value()

        # print("Selected series:", self.meas_selected_series)
        # print("Selected measurement:", self.meas_selected_number)
        meas_index = self.meas_selected_number - 1
        if self.meas_selected_series == 1:
            measurement = self.meas_pressure[meas_index]
        elif self.meas_selected_series == 2:
            measurement = self.meas_mixture[meas_index]
        else:
            raise ValueError

        if self.selected_data == 1:
            data = self.remove_noice(measurement.p_diff)
        elif self.selected_data == 2:
            data = self.remove_noice(measurement.p_abs)
        elif self.selected_data == 3:
            data = self.remove_noice(measurement.ext)

        else:
            raise ValueError

        if self.first_update:
            self.curve_select = self.plot_select.plot(data)
            self.curve_zoom = self.plot_zoom.plot(data)
            self.first_update = False
        else:
            self.curve_select.setData(data)
            self.curve_zoom.setData(data)

        # sets the graphs for the drop-point
        index_of_drop = self.find_drop_index(measurement.p_diff)
        self.plot_select.setXRange(index_of_drop - 10000, index_of_drop + 10000, padding=0)
        self.plot_zoom.setXRange(index_of_drop - 1000, index_of_drop + 3000, padding=0)

        self.update_zoom_region()
        #self.update_zoom_plot()

    def update_zoom_plot(self):
        self.plot_zoom.setXRange(*self.linear_region.getRegion(), padding=0)

    def update_zoom_region(self):
        self.linear_region.setRegion(self.plot_zoom.getViewBox().viewRange()[0])

    def read_to_list(self, folder, start, stop):
        measurements = []
        for i in range(start, stop+1):
            measurements.append(Measurement(folder + str(i) + ".tdms"))
        return measurements

    # From http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
    def savitzky_golay(self, y, window_size, order, deriv=0, rate=1):
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

    def remove_noice(self, data):
        if(self.noice_reduction_number == 0):
            return data
        if(self.noice_reduction_number < 0):
            raise ValueError("Noice reduction can't be negative")
        else:
            return self.savitzky_golay(data, (self.noice_reduction_number*50 + 1), 3) # window size 251, polynomial order 3

    def find_drop_index(self, data):
        # This function will find the index of biggest drop on data vector. Intended for data=measurement.p_diff.
        # Finds the biggest difference drop, or in other words the most negative coarse derivative.
        index_skip = 200        # every single index is not tested for performance reasons
        index_spacing = 500    # the time interval between which difference is tested
        border_margin = np.int(index_spacing / index_skip) + 1 # for cropping off the ends of vector

        greatest_difference = -1
        index_of_greatest = -1

        for i in range(border_margin, (np.int(data.size/index_skip)-border_margin)):
            index = i * index_skip

            # find the difference of data around the index.
            difference = data[index - index_spacing] - data[index + index_spacing]

            if (difference > greatest_difference):
                greatest_difference = difference
                index_of_greatest = index

        return index_of_greatest


def main():
    Main()


main()
