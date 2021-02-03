# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

## @package data_analysis
# This program is for analysing data gathered from spectrometer in Compton
# energy spectrum analyis -physics experiment.

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import simps
import toolbox

import matplotlib.pyplot as plt


# Constants
_eV = 1.60218e-19
_keV = 1.60218e-16


class Measurement:
    """ This class represents a single measurement and its data """
    def __init__(self, path_to_file):
        a, b, c = self.read_data(path_to_file)
        self.count = a
        self.channel = b
        self.energy = (c[0]*b + c[1]) * _keV

    @staticmethod
    def read_data(path_to_file):
        """
        Reads .txt data files
        :param path_to_file:
        :return:    numpy_array count_array         # counts
                    numpy_array channel_array       # respective channels
                    list        coefficients        # a,b (a*x + b)-coefficients
                    list<tuple> calibration         # calibration channel-value -pairs (values in keV)
        """
        if path_to_file[-4:] != ".txt":
            raise "file must be .txt"
        file = open(path_to_file)
        line = ""
        coefficients = []

        # Skips to correct line, and determines wheter it is calibration data or calibrated data
        calibration = False
        for x in range(10):
            if (file.readline() == "Channel\tValue\n"):
                calibration = False
                break
            if (file.readline() =="Channel\tData\n"):
                calibration = True
                break
            if (x==9):
                raise "malformed txt-file"

        # Calibrated data:
        if (calibration == False):
            # Reads channel-value calibration-pairs (unused)
            calibration = []
            for x in range(3):
                channel, value = file.readline().split()
                calibration.append((float(channel), float(value)))
            file.readline()     # skip the "A  B  C" -line

            # Reads A, B, C -coefficients
            c = file.readline().split()
            coefficients = [float(c[1]), float(c[0])]
            for x in range(4): # Skipping to correct line
                line = file.readline()

        # Calibration data:
        if (calibration == True):
            line = file.readline()
            # in case of calibration data coefficients are set to match others
            coefficients = [0.0223487, -0.0131686]

        # Reads the vectors
        channel_list = []   # "Channel"
        count_list = []     # "Data"
        while (line != ""):
            channel, count = line.split()
            channel_list.append(int(channel))
            count_list.append(int(count))
            line = file.readline()
        file.close()

        channel_array = np.array(channel_list)
        count_array= np.array(count_list)

        return count_array, channel_array, coefficients









class Main:
    def __init__(self):
        pg.setConfigOptions(antialias=True, background="w", foreground="k")
        qapp = pg.mkQApp()

        self.angles = [75,80]+list(range(105, 160, 5))
        self.measurement_list = self.read_to_list(self.angles)
        self.measurement = self.measurement_list[2] # initial value is set to measurement at 105 degrees
        self.count = self.measurement.count
        self.energy = self.measurement.energy
        self.energy_range_indx = [np.argmin(np.abs(self.energy-38*_keV)), np.argmin(np.abs(self.energy-58.5*_keV))]



        self.selected_angle = 105
        self.correction = False


        # Control window
        widget2 = QtGui.QWidget()
        widget2.setWindowTitle("Energy spectrum data analysis")
        layout = QtGui.QGridLayout()
        widget2.setLayout(layout)

        labels = ["Angle", "Count correction", "Gauss fit", "Cauchy fit", "Pseudo-Voigt fit"]
        for i, text in enumerate(labels):
            label = QtGui.QLabel()
            label.setText(text)
            layout.addWidget(label, i, 0)

        self.__input_angle = pg.SpinBox(value=self.selected_angle, int=True, minStep=5, step=5)
        self.__input_angle.editingFinished.connect(self.update_change)
        layout.addWidget(self.__input_angle, 0, 1)
        self.__correctionCheckbox = QtGui.QCheckBox()
        layout.addWidget(self.__correctionCheckbox,1,1)
        self.__fitGaussCheckbox = QtGui.QCheckBox()
        layout.addWidget(self.__fitGaussCheckbox, 2, 1)
        self.__fitCauchyCheckbox = QtGui.QCheckBox()
        layout.addWidget(self.__fitCauchyCheckbox, 3, 1)
        self.__fitPseudoVoigtCheckbox = QtGui.QCheckBox()
        layout.addWidget(self.__fitPseudoVoigtCheckbox, 4, 1)
        self.__updateButton = QtGui.QPushButton("Update")
        layout.addWidget(self.__updateButton, 5, 1)
        self.__updateButton.clicked.connect(self.update_change)

        widget2.show()


        # Graph window
        win = pg.GraphicsWindow(title="Energy spectrum data analysis")
        self.plot_count = win.addPlot()#title="Energy spectrum")
        self.plot_diff = win.addPlot()#title="Difference to expected")
        self.plot_klein_nishina = win.addPlot()#title="Klein–Nishina")
        win.nextRow()
        self.plot_cauchy_in_voigt = win.addPlot()#title="Cauchy in Voigt")

        self.plot_int = win.addPlot()
        self.plot_conv = win.addPlot()
        self.plot_conv.setLogMode(x=True, y=None)

        self.set_labels()
        self.curve_count = self.plot_count.plot(self.measurement.energy/_eV, self.measurement.count,
                                                pen=pg.mkPen((70, 170, 255)), name="Mittaus")
        self.curve_gauss = self.plot_count.plot(pen=pg.mkPen((70, 255, 170), width=2), name="Gauss sovite")
        self.curve_cauchy = self.plot_count.plot(pen=pg.mkPen((170, 255, 70), width=2), name="Lorentz")
        self.curve_pseudovoigt = self.plot_count.plot(pen=pg.mkPen((255, 170, 70), width=2), name="Pseudo-Voigt")

        self.curve_diff_measured = self.plot_diff.plot(pen=pg.mkPen((70, 170, 255)), name="Mitattu energiamaksimi")
        self.curve_diff_expected = self.plot_diff.plot(pen=pg.mkPen((70, 255, 170)), name="Laskettu energiamaksimi")

        self.curve_cross_section = self.plot_klein_nishina.plot(pen=pg.mkPen((70, 170, 255)), name="Siroamistodennäköisyys maksimista")
        self.curve_cauchy_fraction = self.plot_cauchy_in_voigt.plot(pen=pg.mkPen((70, 170, 255)), name="Cauchyn osuus Voigt:sta")

        # Integrals
        self.curve_int_cauchy = self.plot_int.plot(pen=(85, 210, 45), name="Lorentz")
        self.curve_int_raw = self.plot_int.plot(pen=(70, 170, 255), name="Mittausdata")

        self.curve_conv_trap = self.plot_conv.plot(pen=(85, 210, 45), name="Puolisuunnikas")
        self.curve_conv_simp = self.plot_conv.plot(pen=(70, 170, 255), name="Simpson")

        win.resize(1600, 900)
        # draws measured and compton-calculated figures
        self.compare_to_excepted()
        # draws Klein-Nishina cross sections and voigt cauchy fractions
        self.additional_analyzing()


        # PyQtGraph main loop
        qapp.exec_()

    @staticmethod
    def read_to_list(angles):
        """
        Reads all single measurements to one list, notice that it also reads
        calibration- and block-measurements at angles 75 and 80 respectively.
        :param angles:
        :return: measurements-list
        """
        measurements = []
        for angle in angles:
            if (angle!=75 and angle!=80):
                measurements.append(Measurement("data/" + str(angle) + ".txt"))
            elif (angle==75):
                measurements.append(Measurement("data/calibration.txt"))
            elif (angle==80):
                measurements.append(Measurement("data/block-80.txt"))
        return measurements

    def set_labels(self):
        """
        Sets labels for all plots in window
        :return:
        """
        self.plot_count.addLegend()
        self.plot_count.setLabel("left", "fotonien lukumäärä")
        self.plot_count.setLabel("bottom", "E", "eV")

        self.plot_diff.addLegend()
        self.plot_diff.setLabel("left", "E", "eV")
        self.plot_diff.setLabel("bottom", "kulma (°) ")

        # self.plot_klein_nishina.addLegend()
        self.plot_klein_nishina.setLabel("left", "Todennäköisyys (%)")
        self.plot_klein_nishina.setLabel("bottom", "kulma (°)")

        # self.plot_cauchy_in_voigt.addLegend()
        self.plot_cauchy_in_voigt.setLabel("left", "Osuus (%)")
        self.plot_cauchy_in_voigt.setLabel("bottom", "kulma (°)")

        self.plot_int.addLegend()
        self.plot_int.setLabel("left", "pinta-ala")
        self.plot_int.setLabel("bottom", "mittaus")

        self.plot_conv.addLegend()
        self.plot_conv.setLabel("left", "pinta-ala")
        self.plot_conv.setLabel("bottom", "osavälejä")

    def update_change(self):
        """
        Updates plots.
        Gets signaled from update-button and change of spinbox value.
        :return:
        """
        self.selected_angle = self.__input_angle.value()
        self.correction = self.__correctionCheckbox.checkState()

        # Energy spectrum
        if self.selected_angle in self.angles:
            self.measurement = self.measurement_list[self.angles.index(self.selected_angle)]
            self.energy = self.measurement.energy
            if (self.correction):
                self.count = self.measurement.count / self.count_correction(self.energy)
            else:
                self.count = self.measurement.count
        else:
            raise ValueError

        # [a:b] is the range between which the mean square distribution fit is made.
        a = self.energy_range_indx[0]
        b = self.energy_range_indx[1]

        print("----- Data analysis -----")
        # print(self.measurement.count)
        # print(self.measurement.energy)
        # print(self.measurement.channel)

        energy_vec = self.energy[a:b] * self.count[a:b]
        # print(energy_vec)
        print("Raw E sum:", np.sum(energy_vec))
        print("Raw E trapz:", np.trapz(energy_vec))
        print("Raw E simps:", simps(energy_vec))
        print("Raw sum:", np.sum(self.count[a:b]))

        # Gauss
        if self.__fitGaussCheckbox.checkState():
            coeff = self.fit_gauss(self.energy[a:b], self.count[a:b])
            # print(coeff)
            gauss = toolbox.gauss(self.energy[a:b], coeff[0], coeff[1], coeff[2])
            # print(gauss)
            self.curve_gauss.setData(self.energy[a:b]/_eV, gauss)

            gauss_energy = self.energy[a:b] * gauss
            print("Gauss E sum:", np.sum(gauss_energy))
            print("Gauss E trapz:", np.trapz(gauss_energy))
            print("Gauss E simps:", simps(gauss_energy))

            print("Gauss sum:", np.sum(gauss))
            print("Gauss trapz:", np.trapz(y=gauss, x=self.energy[a:b]))
            print("Gauss simps:", simps(y=gauss, x=self.energy[a:b]))

            print("Gauss int:", coeff[2])
            # print("Gauss cumulative:", toolbox.gauss_cumulative(self.energy[b], coeff[0], coeff[1], coeff[2]))
            # self.curve_test.setData(gauss)

        else:
            self.curve_gauss.setData([0], [0])

        # Cauchy
        if self.__fitCauchyCheckbox.checkState():
            coeff = self.fit_cauchy(self.energy[a:b], self.count[a:b])
            cauchy = toolbox.cauchy(self.energy[a:b], coeff[0], coeff[1], coeff[2])
            print(self.energy[a])
            print(self.energy[b])
            print(coeff)
            self.curve_cauchy.setData(self.energy[a:b]/_eV, cauchy)

            cauchy_energy = self.energy[a:b] * cauchy
            print("Cauchy E sum:", np.sum(cauchy_energy))
            print("Cauchy E trapz:", np.trapz(cauchy_energy))
            print("Cauchy E simps:", simps(cauchy_energy))

            print("Cauchy sum:", np.sum(cauchy))
            print("Cauchy trapz:", np.trapz(y=cauchy, x=self.energy[a:b]))
            print("Cauchy simps:", simps(y=cauchy, x=self.energy[a:b]))

            print("Cauchy total int:", coeff[2])
            print("Cauchy int:", toolbox.cauchy_cumulative(self.energy[b], coeff[0], coeff[1], coeff[2]) - toolbox.cauchy_cumulative(self.energy[a], coeff[0], coeff[1], coeff[2]))
            # self.curve_cauchy.setData(self.energy[a:b]/_eV, toolboxASDFASDFSA)
        else:
            self.curve_cauchy.setData([0], [0])

        # Pseudo Voigt
        if self.__fitPseudoVoigtCheckbox.checkState():
            coeff = self.fit_pseudovoigt(self.energy[a:b], self.count[a:b])
            pseudo_voigt = toolbox.pseudo_voigt(self.energy[a:b], coeff[0], coeff[1], coeff[2], coeff[3])
            self.curve_pseudovoigt.setData(self.energy[a:b]/_eV, pseudo_voigt)

            pseudo_voigt_energy = self.energy[a:b] * pseudo_voigt
            print("Pseudo Voigt E sum:", np.sum(pseudo_voigt_energy))
            print("Pseudo Voigt E trapz:", np.trapz(y=pseudo_voigt_energy, x=self.energy[a:b]))
            print("Pseudo Voigt E simps:", simps(y=pseudo_voigt_energy, x=self.energy[a:b]))

            print("Pseudo Voigt sum:", np.sum(pseudo_voigt))
            print("Pseudo Voigt trapz:", np.trapz(y=pseudo_voigt, x=self.energy[a:b]))
            print("Pseudo Voigt simps:", simps(y=pseudo_voigt, x=self.energy[a:b]))

            print("Pseudo Voigt int:", coeff[3])
        else:
            self.curve_pseudovoigt.setData([0], [0])

        self.curve_count.setData(self.energy/_eV, self.count)

        print("-----")

    @staticmethod
    def count_correction(energy):
        """
        Returns efficiency factors for energies measured by XR100CR x-ray detector.
        Assumin 500e-6m thick Si-detector and only concidering range of 3-70 keV
        https://www1.aps.anl.gov/files/download/DET/Detector-Pool/Spectroscopic-Detectors/Amptek_CdZnTe/XR100CR_PX2_rev080713.pdf
        :param energy: energy-vector
        :return: (inverse) coefficient vector
        """

        # Photon energies (x), and corresponding measurement efficiencies (y)
        x = np.array([-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]) * _keV # (energy[0] < 0)
        y = np.array([1, 1, 1, 0.4, 0.12, 0.05, 0.023, 0.012, 0.01, 0.01, 0.01, 0.01])

        factor = np.zeros(len(energy))

        for i,en in enumerate(energy):
            for j in range(len(x)-1):
                if (x[j] <= en and en < x[j+1]):
                    # linear interpolation
                    factor[i] = (y[j+1]-y[j])/(x[j+1]-x[j]) * (en-x[j]) + y[j]
                    break
        return factor

    @staticmethod
    def fit_gauss(energy,counts):
        """
        Fits non-linear least squares gauss-fit for scattered peak.
        :param energy: energy-vector
        :param counts: count-vector
        :return: parameters for gauss() -function
        """
        # gauss(x, mu, var)
        guess = [52*_keV, _keV**2, 1e-12]
        coeff, cov = curve_fit(toolbox.gauss, energy, counts, p0 = guess)
        return coeff

    @staticmethod
    def fit_cauchy(energy,counts):
        """
        Fits non-linear least squares cauchy-fit for scattered peak.
        :param energy: energy-vector
        :param counts: count-vector
        :return: parameters for cauchy() -function
        """
        # cauchy(x, x0, gamma)
        guess = [52*_keV, _keV, 1e-12]
        coeff, cov = curve_fit(toolbox.cauchy, energy, counts, p0 = guess)

        return coeff

    @staticmethod
    def fit_pseudovoigt(energy,counts):
        """
        Fits Pseudo-Voigt profile to the given data. Pseudo-Voigt is approximation with linear combination of Gauss and
        Cauchy distributions. This method does not use convolution, but is accurate within 1%.
        https://en.wikipedia.org/wiki/Voigt_profile
        :param energy:
        :param counts:
        :return:
        """
        # pseudo_voigt(x, x0, var, gamma, a=1)
        guess = [52 * _keV, _keV**2, _keV, 1e-12]
        coeff, cov = curve_fit(toolbox.pseudo_voigt, energy, counts, p0=guess)
        return coeff

    def compare_to_excepted(self):
        """
        Plots measured energies and theoretically predicted energies
        :return:
        """

        # [a:b] is the range between which the mean square distribution fit is made.
        a = self.energy_range_indx[0]
        b = self.energy_range_indx[1]

        measured = []
        excepted = []

        for i, meas in enumerate(self.measurement_list):
            energy = meas.energy
            count = meas.count / self.count_correction(self.energy)
            # using Cauchy distribution instead of Voigt, because it seems to be more accurate
            coeff = self.fit_cauchy(energy[a:b], count[a:b])

            measured_energy_maxium_eV = coeff[0] / _eV
            excepted_energy_maxium_eV = toolbox.e_2f_fun(self.angles[i]*2*np.pi/360) / _eV

            measured.append(measured_energy_maxium_eV)
            excepted.append(excepted_energy_maxium_eV)


        # (measurement data and block-test at 80 degrees is omitted)
        self.curve_diff_measured.setData(self.angles[2:], measured[2:])
        self.curve_diff_expected.setData(self.angles[2:], excepted[2:])



    def interpolate(self,vec_x,vec_y):
        vec_size = 10*len(vec_x)
        new_vec_x = np.arange(vec_x[0],vec_x[-1],float((vec_x[-1]-vec_x[0])/vec_size))
        new_vec_y = np.zeros(vec_size)

        for i, x in enumerate(new_vec_x):
            for j in range(len(vec_x)-1):
                if (vec_x[j] <= x and x < vec_x[j+1]):
                    new_vec_y[i] = (vec_y[j+1]-vec_y[j])/(vec_x[j+1]-vec_x[j]) * (x-vec_x[j]) + vec_y[j]
                    break

        return new_vec_x, new_vec_y

    def additional_analyzing(self):
        """
        Plots Klein-Nishina cross sections and Voigt Cauchy-fractions
        :return:
        """

        # [a:b] is the range between which the mean square distribution fit is made.
        a = self.energy_range_indx[0]
        b = self.energy_range_indx[1]

        cross_sections = []
        cauchy_fractions = []

        for i, meas in enumerate(self.measurement_list):
            energy = meas.energy
            count = meas.count / self.count_correction(self.energy)
            coeff = self.fit_pseudovoigt(energy[a:b], count[a:b])

            cross_section_arbitary = toolbox.klein_nishina(self.angles[i] * 2 * np.pi / 360)/(toolbox.klein_nishina(0)) * 100
            cauchy_fraction = toolbox.voigt_cauchy_percentage(coeff[1],coeff[2])*100

            cross_sections.append(cross_section_arbitary)
            cauchy_fractions.append(cauchy_fraction)

        # (measurement data and block-test at 80 degrees is omitted)
        self.curve_cross_section.setData(self.angles[2:], cross_sections[2:])
        self.curve_cauchy_fraction.setData(self.angles[2:], cauchy_fractions[2:])

        # ASDFASDF
        analyticals = []
        numerics = []
        print("----- Cauchy integrals -----")
        for i, meas in enumerate(self.measurement_list):
            count = meas.count / self.count_correction(meas.energy)
            coeff = self.fit_cauchy(meas.energy[a:b], count[a:b])

            # cauchy = toolbox.cauchy(self.energy[a:b], coeff[0], coeff[1], coeff[2])

            cauchy = toolbox.cauchy_cumulative(meas.energy[b], coeff[0], coeff[1], coeff[2]) - toolbox.cauchy_cumulative(
                      meas.energy[a], coeff[0], coeff[1], coeff[2])

            raw = np.trapz(y=count[a:b], x=meas.energy[a:b])

            print("Raw", raw)
            print("Cauchy", cauchy)
            analyticals.append(cauchy)
            numerics.append(raw)

        print("-----")

        interpolated_angle, interpolated_raw = self.interpolate(self.angles[2:],numerics[2:])
        interpolated_angle, interpolated_cauchy = self.interpolate(self.angles[2:], analyticals[2:])

        c, d = np.argmin(np.abs(interpolated_angle-117)) , np.argmin(np.abs(interpolated_angle-141))

        area_raw = np.trapz(y=interpolated_raw, x=interpolated_angle)
        area_cauchy = np.trapz(y=interpolated_cauchy, x=interpolated_angle)
        area_raw_117_141 = np.trapz(y=interpolated_raw[c:d], x=interpolated_angle[c:d])
        area_cauchy_117_141 = np.trapz(y=interpolated_cauchy[c:d], x=interpolated_angle[c:d])

        print("Raw: area in 117-141 / total area: ", area_raw_117_141, "/", area_raw, "=", area_raw_117_141/area_raw)
        print("Cauchy: area in 117-141 / total area: ", area_cauchy_117_141, "/", area_cauchy, "=", area_cauchy_117_141/area_cauchy)





        self.curve_int_cauchy.setData(self.angles[2:], analyticals[2:])
        self.curve_int_raw.setData(self.angles[2:], numerics[2:])

        meas = self.measurement_list[2]
        count = meas.count / self.count_correction(meas.energy)
        print("Exact:", analyticals[2])

        a = self.energy_range_indx[0]
        b = self.energy_range_indx[1]
        print(a, b)

        # stepcounts = [10]
        stepcounts = np.logspace(1, 3, num=1000)

        trapzs = []
        simpss = []

        cauchy_params = self.fit_cauchy(meas.energy[a:b], count[a:b])
        # OK
        # print(cauchy_params)

        a = 6.08859964548e-15
        b = 9.3720685487e-15

        for i_stepcount, stepcount in enumerate(stepcounts):
            x = np.arange(start=float(a), stop=float(b), step=float((b - a) / stepcount))

            cauchy_vec = toolbox.cauchy(x, cauchy_params[0], cauchy_params[1], cauchy_params[2])

            trapz = np.trapz(y=cauchy_vec, x=x)
            trapzs.append(trapz)
            simp = simps(y=cauchy_vec, x=x)
            simpss.append(simp)

            #print(trapz)
            #print(simp)

        self.curve_conv_trap.setData(stepcounts, trapzs)
        self.curve_conv_simp.setData(stepcounts, simpss)

        print("-----")


def main():
    Main()

main()
