# FYS-1320 Methods in Physics

This repository contains our analysis code for the
[Tampere University of Technology](https://en.wikipedia.org/wiki/Tampere_University_of_Technology)
course
[FYS-1320 Methods in Physics](https://www.tuni.fi/studentsguide/curriculum/course-units/tut-cu-g-36085).

The course is composed of three experiments, each of which has a pre-report and a seminar presentation.

We'd be delighted to see our code being used in the course materials. Please contact us if licensing changes are necessary.

Mika Mäki & Alpi Tolvanen, 2016-2017


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
This is an experiment about [Compton scattering](https://en.wikipedia.org/wiki/Compton_scattering) of gamma rays emitted by [Americium-241](https://en.wikipedia.org/wiki/Americium-241). The statistical distributions used for data analysis are the [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution), the [Cauchy distribution](https://en.wikipedia.org/wiki/Cauchy_distribution) and the [Voigt profile](https://en.wikipedia.org/wiki/Voigt_profile). The [Monte-Carlo method](https://en.wikipedia.org/wiki/Monte_Carlo_method) and [Klein–Nishina formula](https://en.wikipedia.org/wiki/Klein%E2%80%93Nishina_formula) are also used.

Files
- convergence.py
  - Creates a 3D graph of Monte-Carlo integral error for visualising convergence
- **data_analysis.py**
  - Analysis software for data generated by the measurement system
- monte_carlo.py
  - Generates a 2D graph of a Monte-Carlo integral
- **project.py**
  - The program for generating the graphs of our pre-report
- test_1.html
  - An example of a 3D graph generated using convergence.py. Can be opened using a web browser
- toolbox.py
  - A library for data_analysis.py and project.py. Contains functions for Compton scattering and statistical distributions.
