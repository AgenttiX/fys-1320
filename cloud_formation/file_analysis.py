# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

## @package file_analysis
# This program is for reverse engineering the TDMS -> TXT conversion process done by
# the convert_TDMS_to_ASCII.exe of Tampere University of Technology student laboratory of physics

# npTDMS is a library for reading TDMS files created by LabVIEW
# https://pypi.python.org/pypi/npTDMS/
import nptdms
import numpy as np
import pyqtgraph as pg


class TDMS_Measurement:
    """ This class represents a single measurement loaded from a TDMS file """
    def __init__(self, path):
        # Read the TDMS file
        tdms_file = nptdms.TdmsFile(path)
        # The group and channel names have been reverse engineered using a hex editor
        tdms_p_diff = tdms_file.object("GroupName", "Paine_ero_kPa")
        tdms_p_abs = tdms_file.object("GroupName", "Absoluuttipaine_kPa")
        tdms_ext = tdms_file.object("GroupName", "Light extinction_V")

        # Extract data from the objects
        self.p_diff = tdms_p_diff.data
        self.p_abs = tdms_p_abs.data
        self.ext = tdms_ext.data
        self.time = tdms_p_diff.time_track()


class TXT_Measurement:
    """ This class represents a single measurement loaded from a TXT file """
    def __init__(self, path):
        data = np.loadtxt(path)
        self.p_diff = data[0]
        self.p_abs = data[1]
        self.ext = data[2]


class Main:
    def __init__(self):

        folder = "data/by_pressure/"

        self.show_contents(folder + "1.tdms")

        start = 1
        stop = 3

        # File loading
        tdms_measurements = []
        txt_measurements = []
        for i in range(start, stop + 1):
            tdms_measurements.append(TDMS_Measurement(folder + str(i) + ".tdms"))
            txt_measurements.append(TXT_Measurement(folder + str(i) + ".txt"))

        # Data extraction
        meas_number = 0

        tdms_p_diff = tdms_measurements[meas_number].p_diff
        txt_p_diff = txt_measurements[meas_number].p_diff

        tdms_p_abs = tdms_measurements[meas_number].p_abs
        txt_p_abs = txt_measurements[meas_number].p_abs

        tdms_ext = tdms_measurements[meas_number].ext
        txt_ext = txt_measurements[meas_number].ext

        print("Pressure difference")
        print(tdms_p_diff)
        print(txt_p_diff)

        print("Absolute pressure")
        print(tdms_p_abs)
        print(txt_p_abs)

        print("Extinction")
        print(tdms_ext)
        print(txt_ext)

        # Computation
        fit_p_diff = np.polyfit(tdms_p_diff, txt_p_diff, 1)
        fit_p_abs = np.polyfit(tdms_p_abs, txt_p_abs, 1)
        fit_ext = np.polyfit(tdms_ext, txt_ext, 1)

        print("Linear fits")
        print(fit_p_diff)
        print(fit_p_abs)
        print(fit_ext)

        err_p_diff = tdms_p_diff*fit_p_diff[0] + fit_p_diff[1] - txt_p_diff
        err_p_abs = tdms_p_abs*fit_p_abs[0] + fit_p_abs[1] - txt_p_abs
        err_ext = tdms_ext*fit_ext[0] + fit_ext[1] - txt_ext

        print("Maximum errors")
        print(err_p_diff.max())
        print(err_p_abs.max())
        print(err_ext.max())

        # Plotting
        qapp = pg.mkQApp()

        win = pg.GraphicsWindow(title="Cloud formation file analysis")

        plot_p_diff = win.addPlot(title="Pressure difference")
        plot_p_diff.plot(tdms_p_diff, txt_p_diff, pen=None, symbol="s")

        plot_p_abs = win.addPlot(title="Absolute pressure")
        plot_p_abs.plot(tdms_p_abs, txt_p_abs, pen=None, symbol="s")

        plot_ext = win.addPlot(title="Extinction")
        plot_ext.plot(tdms_ext, txt_ext, pen=None, symbol="s")

        # Main loop
        qapp.exec_()

    @staticmethod
    def show_contents(path):
        """
        Prints the contents of a TDMS file (channel display seems not to work)
        :param path: path to file (str)
        :return: -
        """
        tdms_file = nptdms.TdmsFile("data/by_pressure/1.tdms")
        groups = tdms_file.groups()
        for group in groups:
            print("Group:", group)
            print("Channels:")
            tdms_file.group_channels(group)


def main():
    Main()

main()
