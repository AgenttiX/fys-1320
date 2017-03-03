# fys-1320
Tampere University of Technology FYS-1320 Methods in Physics

The course is composed of three experiments, each of which has a pre-report and a seminar presentation. To prevent plagiarism, our measurement data, pre-reports and seminar presentations aren't publicly available.

- [Course description in Finnish](http://www.tut.fi/opinto-opas/wwwoppaat/opas2016-2017/perus/laitokset/Fysiikka/FYS-1320.html)
- [Course description in English](http://www.tut.fi/opinto-opas/wwwoppaat/opas2016-2017/kv/laitokset/Fysiikka/FYS-1326.html)
- [Moodle page of the course](https://moodle2.tut.fi/course/view.php?id=7583)


## Gas analysis
This is an [NDIR spectroscopy](https://en.wikipedia.org/wiki/Nondispersive_infrared_sensor) experiment to measure gas concentrations based on the characteristic infrared absorption of the gas molecules.

Files
- **gas_analysis.py**
  - For easy computing of the gas concentrations from the files provided by the measurement software
- **laske_konsentraatiot.m**
  - The script we used to create most of the graphs of our pre-report and seminar presentation. The code is written specifically for Octave and contains some clever abuse of its plotting system. It's dependent on the two other Octave files.
- lue_putki.m
  - An Octave function for reading the data of the measurements of a single sample pipe
- pilkut.m
  - An Octave function for dealing with the measurement files that use comma as their decimal separator

## Cloud formation
This is an experiment to measure cloud formation on condensation particles using laser absorption, [Beer-Lambert law](https://en.wikipedia.org/wiki/Beer%E2%80%93Lambert_law) and [Mie theory](https://en.wikipedia.org/wiki/Mie_scattering).

We were unable to translate the growth simulation and Mie theory code provided by the course to Python due to the lack of comments and other documentation, and thereby our toolbox.py is dependent on Matlab API for Python and the course example code.

Files
- 12_calculate_pressures.py
  - For calculating the necessary pressures for the particle growth with different pressure changes, as described in the task 12 of the instructions
- 9_extinction.py
  - For calculating extinction as described in the task 9 of the instructions
- Doxyfile
  - An automatically generated file for generating documentation of our code using Doxygen
- **data_analysis.py**
  - Our software for analysing data created by the measurement system, and thereby the primary source of the graphs in our presentation. This program is dependent on toolbox_2.py.
- example_dp.py
  - Python translation of ESIMERKKI_dp.m provided by the course
- example_extinction.py
  - Python translation of ESIMERKKI_EKSTINKTIOTEHOKKUUS.m provided by the course
- example_interp.py
  - Python translation of ESIMERKKI_interp.m provided by the course
- file_analysis.py
  - For reverse engineering the conversion coefficients of the convert_TDMS_to_ASCII.exe of the measurement computer
- matlab_test.py
  - For testing whether the Matlab API for Python has been set up correctly
- number_counts.py
  - A numpy array required by data_analysis.py for storing the results
- number_counts_2.py
  - A numpy array required by data_analysis.py for storing the results
- project.py
  - Old version of the primary program used to generate data and graphs for the pre-report
- **project_mod.py**
  - The primary program for generating data and graphs for the pre-report
- smallest_particles.npy
  - Required by the data_analysis.py
- **toolbox.py**
  - A library for project.py, project_mod.py and data_analysis.py. Requires Matlab API for Python and the course example code files.
- **toolbox_2.py**
  - A library for data_analysis.py. Requires toolbox.py.

## Energy spectrum analysis
An experiment about [Compton scattering](https://en.wikipedia.org/wiki/Compton_scattering).

Files
- project.py
  - The program for generating the graphs of our pre-report
- toolbox.py
  - A library for project.py

Work in progress.
