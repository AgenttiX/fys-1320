# Mika "AgenttiX" Mäki & Alpi Tolvanen, 2017

## @package data_analysis
# This program is for analysing data gathered from spectrometer in Compton
# energy spectrum analyis -physics experiment.


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
from pyqtgraph.Qt import QtGui
import numpy as np
from scipy.optimize import curve_fit
import toolbox


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


    def read_data(self, path_to_file):
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

        self.set_labels()
        self.curve_count = self.plot_count.plot(self.measurement.energy/_eV, self.measurement.count,
                                                pen=pg.mkPen((70, 170, 255)), name="Mittaus")
        self.curve_gauss = self.plot_count.plot(pen=pg.mkPen((70, 255, 170), width=2), name="Gauss sovite")
        self.curve_cauchy = self.plot_count.plot(pen=pg.mkPen((170, 255, 70), width=2), name="Lorentz")
        self.curve_pseudovoigt = self.plot_count.plot(pen=pg.mkPen((255, 170, 70), width=2), name="Pseudo-Voigt")

        self.curve_diff_measured = self.plot_diff.plot(pen=pg.mkPen((70, 170, 255)), name="Mitattu energiamaksimi")
        self.curve_diff_expected = self.plot_diff.plot(pen=pg.mkPen((70, 255, 170)), name="Laskettu energiamaksimi")

        self.curve_cross_section = self.plot_klein_nishina.plot(pen=pg.mkPen((70, 170, 255)), name="Siroamistodennäköisyys maksimista")
        self.curve_cauchy_fraction = self.plot_cauchy_in_voigt.plot(pen=pg.mkPen(70, 255, 170), name="Cauchyn osuus Voigt:sta")


        win.resize(1600, 900)
        # draws measured and compton-calculated figures
        self.compare_to_excepted()
        # draws Klein-Nishina cross sections and voigt cauchy fractions
        self.additional_analyzing()


        # PyQtGraph main loop
        qapp.exec_()


    def read_to_list(self,angles):
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

        #self.plot_klein_nishina.addLegend()
        self.plot_klein_nishina.setLabel("left", "Todennäköisyys (%)")
        self.plot_klein_nishina.setLabel("bottom", "kulma (°)")

        #self.plot_cauchy_in_voigt.addLegend()
        self.plot_cauchy_in_voigt.setLabel("left", "Osuus (%)")
        self.plot_cauchy_in_voigt.setLabel("bottom", "kulma (°)")


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

        # Gauss
        if (self.__fitGaussCheckbox.checkState()):
            coeff = self.fit_gauss(self.energy[a:b],self.count[a:b])
            self.curve_gauss.setData(self.energy[a:b]/_eV, toolbox.gauss(self.energy[a:b],coeff[0],coeff[1], coeff[2]))
        else:
            self.curve_gauss.setData([0],[0])
        # Cauchy
        if (self.__fitCauchyCheckbox.checkState()):
            coeff = self.fit_cauchy(self.energy[a:b], self.count[a:b])
            self.curve_cauchy.setData(self.energy[a:b]/_eV, toolbox.cauchy(self.energy[a:b], coeff[0], coeff[1], coeff[2]))
        else:
            self.curve_cauchy.setData([0],[0])
        # Pseudo Voigt
        if (self.__fitPseudoVoigtCheckbox.checkState()):
            coeff = self.fit_pseudovoigt(self.energy[a:b], self.count[a:b])
            self.curve_pseudovoigt.setData(self.energy[a:b]/_eV, toolbox.pseudo_voigt(self.energy[a:b], coeff[0], coeff[1], coeff[2], coeff[3]))
            n = toolbox.voigt_cauchy_percentage(coeff[1], coeff[2])
            #print("toolbox: % cauchy/gauss in voigt", n, "/", 1 - n)

        else:
            self.curve_pseudovoigt.setData([0],[0])



        self.curve_count.setData(self.energy/_eV, self.count)

    def count_correction(self,energy):
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

    def fit_gauss(self,energy,counts):
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

    def fit_cauchy(self,energy,counts):
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

    def fit_pseudovoigt(self,energy,counts):
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

def main():
    Main()


main()