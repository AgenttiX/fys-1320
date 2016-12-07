# FYS-1320 Gas Analysis
# This software is written to analyse gas samples measured with
# the NDIR setup of TUT physics student laboratory
# It is based on code which we wrote on Octave, and that code in
# turn is somewhat based on the example code of the course.

# Copyright Mika "AgenttiX" Mäki and Alpi Tolvanen, 2016


# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import numpy as np
# import numpy.linalg
import csv
import scipy.io

from pyqtgraph.Qt import QtGui
import pyqtgraph as pg

# import matplotlib.pyplot as plt
# import mpl_toolkits.mplot3d


def readfile(path, filename):
    """ Parses a file created by the measurement software (which is based on LabVIEW)
    :param path: str
    :param filename: str
    :return: np.array
    """
    rawfile = open(path + filename)
    file = csv.reader(rawfile, delimiter="\t")
    datalist = []
    for row in file:
        datalist.append(row)

    data = np.zeros((len(datalist), 4))

    for i, row in enumerate(datalist):
        for j, cell in enumerate(row):
            data[i, j] = float(cell.replace(",", "."))

    return data


class Pipe:
    """
    This class represents a single sample pipe with all the data collected and computed from it
    """
    def __init__(self, path, filenames):
        """
        :param path: string
        :param filenames: np.array
        """
        self.datapoints = 60
        filters = filenames[:, 1].size

        # Create matrixes for data
        # Two columns: measurement and refernce
        # A row for each measurement
        # Each cell contains a vector of the measurement data
        self.data = np.zeros((filters, 2, self.datapoints))
        self.data2 = np.zeros((filters, 2, self.datapoints))
        self.temp_data = np.zeros((filters, 2, self.datapoints))

        # Fill matrixes with data
        for filter_i in range(filters):
            measurement = readfile(path, filenames[filter_i, 0])
            reference = readfile(path, filenames[filter_i, 1])

            # First data channel (the one actually connected to the light sensor
            self.data[filter_i, 0, :] = measurement[:self.datapoints, 1]
            self.data[filter_i, 1, :] = reference[:self.datapoints, 1]

            # Second data channel (so far unconnected)
            self.data2[filter_i, 0, :] = measurement[:self.datapoints, 2]
            self.data2[filter_i, 1, :] = reference[:self.datapoints, 2]

            # Temperature sensor resistance
            self.temp_data[filter_i, 0, :] = measurement[:self.datapoints, 3]
            self.temp_data[filter_i, 1, :] = reference[:self.datapoints, 3]

        # Compute properties from data
        self.mean = np.zeros((filters, 2))
        self.mean2 = np.zeros((filters, 2))
        self.temp_mean = np.zeros((filters, 2))

        # avgavgerr is in Finnish keskiarvon keskivirhe
        self.avgavgerr = np.zeros((filters, 2))
        self.temp_avgavgerr = np.zeros((filters, 2))

        # Fill the property matrixes with data
        for i in range(filters):
            for j in range(2):
                self.mean[i, j] = np.mean(self.data[i, j, :])
                self.mean2[i, j] = np.mean(self.data2[i, j, :])
                self.temp_mean[i, j] = np.mean(self.temp_data[i, j, :])

                self.avgavgerr[i, j] = np.std(self.data[i, j, :]) / np.sqrt(self.datapoints)
                self.temp_avgavgerr[i, j] = np.std(self.temp_data[i, j, :]) / np.sqrt(self.datapoints)

        self.avgavgerr_max = np.amax(self.avgavgerr)
        self.temp_avgavgerr_max = np.amax(self.temp_avgavgerr)

        self.avgavgerr_rel = self.avgavgerr / self.mean
        self.avgavgerr_max_rel = np.amax(self.avgavgerr_rel)

        self.temp_avgavgerr_rel = self.temp_avgavgerr / self.temp_mean
        self.temp_avgavgerr_rel_max = np.amax(self.temp_avgavgerr)

        self.I = self.mean[:, 0]
        self.I_norm = self.I * (self.mean[:, 1] / self.mean[0, 1])


class MatrixComputer:
    def __init__(self):
        self.gasnames = ["H20", "CO2", "CO", "CH4", "C2H2"]

        # Gas matrix creation
        # The file "kaasut_carbon.mat" is provided in the course materials
        abs_dict = scipy.io.loadmat("vectors/kaasut_carbon.mat")
        abs_mat = np.array([abs_dict["H2O"][:, 0], abs_dict["CO2"][:, 0], abs_dict["CO"][:, 0], abs_dict["CH4"][:, 0], abs_dict["C2H2"][:, 0]])
        abs_mat *= 1e-4
        lambda_vec = abs_dict["lambda"][:, 0]

        # dlambda
        dlambda_temp = np.diff(lambda_vec)
        self.__dlambda = np.zeros(dlambda_temp.size + 1)
        self.__dlambda[0] = abs(dlambda_temp[0])
        self.__dlambda[1:] = abs(dlambda_temp)

        # Filter matrix creation
        filternames = ["2710", "3060", "3220", "3977", "4270", "4740", "5756"]
        self.filtercount = len(filternames)

        # Filter matrix preparation
        # The .mat files are provided in the course materials
        for i, filtername in enumerate(filternames):
            filter_dict = scipy.io.loadmat("vectors/" + filtername + ".mat")
            filter_vec = filter_dict["NB" + filtername]
            if i == 0:
                self.__filter_mat = np.zeros((len(filternames), filter_vec.size))
            self.__filter_mat[i, :] = filter_vec[:, 0]

        # It was necessary to change 1e-2 to 1e-6 from Octave, reason unknown.
        A = np.array([
            np.sum(np.array([abs_mat[0, :]] * self.filtercount) * self.__filter_mat * 1e-6 * np.array([[self.__dlambda] * self.filtercount]), axis=2)[0],
            np.sum(np.array([abs_mat[1, :]] * self.filtercount) * self.__filter_mat * 1e-6 * np.array([[self.__dlambda] * self.filtercount]), axis=2)[0],
            np.sum(np.array([abs_mat[2, :]] * self.filtercount) * self.__filter_mat * 1e-6 * np.array([[self.__dlambda] * self.filtercount]), axis=2)[0],
            np.sum(np.array([abs_mat[3, :]] * self.filtercount) * self.__filter_mat * 1e-6 * np.array([[self.__dlambda] * self.filtercount]), axis=2)[0],
            np.sum(np.array([abs_mat[4, :]] * self.filtercount) * self.__filter_mat * 1e-6 * np.array([[self.__dlambda] * self.filtercount]), axis=2)[0],
        ])

        p = 1.013e5
        k = 1.3806505e-23
        T = 298
        L = 0.2
        self.__A = np.transpose(A * p / (k * T) * L)

    def compute_pipe(self, sample, reference):
        """
        :param sample: Pipe
        :param reference: Pipe
        :return:
        """

        # It was necessary to change 1e-2 to 1e-6 from Octave, reason unknown.
        b_vec = (1 - sample.I_norm / reference.I_norm) * np.sum(self.__filter_mat * np.array([[self.__dlambda] * self.filtercount]) * 1e-6, axis=2)[0]

        return np.linalg.lstsq(self.__A, b_vec)


def main():
    # File names for pipes
    # We had three samples and two measurement sessions
    # At the first session we measured the sample pipe of the laboratory, and
    # at the second session we collected samples from Hiukkanen quild room
    # and near the main doors of TUT Tietotalo.

    M1_N2_filenames = np.array([
        ["N2_2710-050.txt", "N2_2710-050-reference.txt"],
        ["N2_3060-060.txt", "N2_3060-060_reference.txt"],
        ["N2_3220-060.txt", "N2_3220-060_reference.txt"],
        ["N2_3977-065.txt", "N2_3977-065.txt"],
        ["N2_4270-070.txt", "N2_4270-070_reference.txt"],
        ["N2_4740-140.txt", "N2_4740-140_reference.txt"],
        ["N2_5756-080.txt", "N2_5756-080_reference.txt"]]
    )

    M1_sample_filenames = np.array([
        ["Sample_2710-050.txt", "Sample_2710-050_reference.txt"],
        ["Sample_3060-060.txt", "Sample_3060-060_reference.txt"],
        ["Sample_3220-060.txt", "Sample_3220-060_reference.txt"],
        ["Sample_3977-060.txt", "Sample_3977-060.txt"],
        ["Sample_4270-070.txt", "Sample_4270-070_reference.txt"],
        ["Sample_4740-140_2.txt", "Sample_4740-140_2_reference.txt"],
        ["Sample_5756-080.txt", "Sample_5756-080_reference.txt"]
    ])

    M2_N2_filenames = np.array([
        ["M2_N2_2710-050.txt", "M2_N2_3977-065.txt"],
        ["M2_N2_3060-060.txt", "M2_N2_3977-065.txt"],
        ["M2_N2_3220-060.txt", "M2_N2_3977-065.txt"],
        ["M2_N2_3977-065.txt", "M2_N2_3977-065.txt"],
        ["M2_N2_4270-070.txt", "M2_N2_3977-065.txt"],
        ["M2_N2_4740-140.txt", "M2_N2_3977-065.txt"],
        ["M2_N2_5756-080.txt", "M2_N2_3977-065.txt"]
    ])

    M2_Hiu_filenames = np.array([
        ["M2_Hiu_2710-050.txt", "M2_Hiu_3977-065.txt"],
        ["M2_Hiu_3060-060.txt", "M2_Hiu_3977-065.txt"],
        ["M2_Hiu_3220-060.txt", "M2_Hiu_3977-065.txt"],
        ["M2_Hiu_3977-065.txt", "M2_Hiu_3977-065.txt"],
        ["M2_Hiu_4270-070.txt", "M2_Hiu_3977-065.txt"],
        ["M2_Hiu_4740-140.txt", "M2_Hiu_3977-065.txt"],
        ["M2_Hiu_5756-080.txt", "M2_Hiu_3977-065.txt"]
    ])

    M2_out_filenames = np.array([
        ["M2_out_2710-050.txt", "M2_out_3977-065.txt"],
        ["M2_out_3060-060.txt", "M2_out_3977-065.txt"],
        ["M2_out_3220-060.txt", "M2_out_3977-065.txt"],
        ["M2_out_3977-065.txt", "M2_out_3977-065.txt"],
        ["M2_out_4270-070.txt", "M2_out_3977-065.txt"],
        ["M2_out_4740-140.txt", "M2_out_3977-065.txt"],
        ["M2_out_5756-080.txt", "M2_out_3977-065.txt"]
    ])

    # Pipe creation
    M1_N2 = Pipe("data/", M1_N2_filenames)
    M1_sample = Pipe("data/", M1_sample_filenames)

    M2_N2 = Pipe("data2/", M2_N2_filenames)
    M2_Hiu = Pipe("data2/", M2_Hiu_filenames)
    M2_out = Pipe("data2/", M2_out_filenames)

    # Matrix computation
    computer = MatrixComputer()

    sample_results = computer.compute_pipe(M1_sample, M1_N2)
    Hiu_results = computer.compute_pipe(M2_Hiu, M2_N2)
    out_results = computer.compute_pipe(M2_out, M2_N2)

    c_mat = np.array([
        sample_results[0],
        Hiu_results[0],
        out_results[0]
    ])

    # Values to ppm
    c_mat *= 1e6

    # Output printing
    print("     sample             Hiukkanen           out")
    for i, row in enumerate(np.transpose(c_mat)):
        print(computer.gasnames[i], row)

    # Additional data output
    residual_mat = np.array([
        sample_results[1],
        Hiu_results[1],
        out_results[1]
    ])
    residual_mat *= 1e-6
    print()
    print(residual_mat)

    rank_mat = np.array([
        sample_results[2],
        Hiu_results[2],
        out_results[2]
    ])
    print()
    print(rank_mat)

    singular_mat = np.array([
        sample_results[3],
        Hiu_results[3],
        out_results[3]
    ])
    print()
    print(singular_mat)

    # Pyqtgraph preparation
    app = pg.mkQApp()
    pg.setConfigOptions(antialias=True)

    win = pg.GraphicsWindow(title="N2 reference changes")
    plot = win.addPlot(title="N2 reference changes")
    for i in range(7):
        plot.plot(M1_N2.data[i, 1, :])

    plot.plot(np.mean(M1_N2.data[:, 1, :], axis=0), pen=(255, 0, 0))

    win2 = QtGui.QMainWindow()
    win2.resize(1280, 720)
    win2.setWindowTitle("Measurement 1 references")
    imv = pg.ImageView()
    win2.setCentralWidget(imv)

    ref_image = np.zeros((60, 15))

    # M1_N2.data[0,1,0] = 0
    # M1_sample.data[0,1,0] = 0

    # x = time
    # y = measurement (top first)
    ref_image[:, :7] = np.transpose(M1_N2.data[:, 1, :])

    # Is already zeros
    # ref_image[:,7] = np.zeros(60)

    ref_image[:, 8:] = np.transpose(M1_sample.data[:, 1, :])
    imv.setImage(ref_image)

    win2.show()

    # Pyqtgraph loop
    app.exec_()

main()
