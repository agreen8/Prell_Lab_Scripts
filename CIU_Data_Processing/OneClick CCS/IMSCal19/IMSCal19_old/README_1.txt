# Instructions

#### About this program ######
Uses a purist Bayesian approach for calibration because calibration functions are highly non-linear and not easily inverted. 

Unknown parameters -> a (alpha scaling factor), g (gamma scaling factor), c (radial exponent coefficient), and t0 (time offset). 

These parameters have an assigned plausible range and prior probability distribution. These are fixed internally.

An important goal of this program is correct propagation of uncertainties. A global uncertainty as a fixed percentage is given for combined given drift times and uncertainties in reference CCS measurements. 

For calibration, a reference file (calibrants) is required to perform calibration. An input file (ions to be calibrated) supplied will have predictions in the output file. Reference file and input file have exactly the same format.

############

##### Reference and Input File format ##########

The format is one line per ion with space separated values:

id mass charge ccs drift

where id is identity of ion, mass is mass in daltons, charge is charge in integer, ccs is reference ccs value in sq. Angstrom units, and drift is experimental drift time in ms.

For input file, where ions do not have CCS, you can put any number in ccs field. However, ccs field has to be filled! Put 0 or 999 or any other number if that's easy.

##################################################


############## Output File Format #########################

The output format is a text file containing three sections in csv format with headers:

[PARAMETERS]
1. g, dg, a, da, c, dc, t0, dt0

These are the inferred values of the parameters g, a, c, and t0 with associated uncertainties. This information is mainly for diagnostic purposes: these values are not used directly in the calibration.

[DIAGNOSTICS]
2. ID, Mass, Z, Mobility, Alpha, Gamma, Model Velocity, Exp. Velocity, Error%

This is diagnostics concerning the reference ions (calibrants) including the model average ion velocities at the above values of a, g,
c, and t0, and the residual percentage velocity errors. 

[CALIBRATED DATA]
3. ID, Mass, Z, Drift, CCS, CCS Std. Dev.
This reports the calibrated data from input file with calculated CCS with uncertainties. Note that the CCS values and its uncertainties represent the mean and standard deviations of samples from a posterior probability distribution over the unknown parameters, and cannot be obtained directly from the parameter estimates in section 1 [PARAMETERS].

##########################################################


############ How to Run the Application ###################

Optional values and defaults in parenthesis

/bin/TWaveCalibrate.exe

    - ref   :reference data file
    (-input)    :data to be calibrated
    -velocity   :wave velocity (m/s)
    -voltage    :wave height (V)
    -pressure   :(mbar)
    -accuracy   : (%)
    (-temp)     : cell temperature (300K)
    (-length)   :(0.254m)
    (-lambda)   :wavelength(0.012m)
    (-t0)       :fixed time offset (s)
    (-a)        :fixed vel. relaxation parameter
    (-c)        :fixed raidal parameter


You can also batch run many calibrations at once by using TWaveCalibrate.bat file. 

Set up calibration syntax in .bat file for multiple calibrations at once!

Template files are provided for each file indicated above!

####################################################