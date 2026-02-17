# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 14:44:38 2025

@author: Austin
"""

# %%


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# from scipy.optimize import minimize, curve_fit
import glob
import subprocess
import time
from CIUSuite3 import CIU2_Main, Raw_Data_Import 

# from CIUSuite3 import CIU_analysis_obj
from scipy import interpolate, integrate
# import pandas as pd
import os
from pathlib import Path
# import itertools
import scipy.signal as sp
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
from matplotlib.backends.backend_pdf import PdfPages
# import copy as cx
import warnings
from scipy.optimize import curve_fit
warnings.filterwarnings("ignore")
# import scipy.signal
# import scipy.interpolate
import scipy.stats
# import pygubu
import logging
logger = logging.getLogger('main')


# %%
total_times = []

#User inputs for poly-alanine extraction
extractor_path = r"TWIMExtract\jars\TWIMExtract.jar"
# polya_raw_data_location = r"C:\Users\Austin\Documents\Grad School\One_click_CCS_calibration\Polya_raw_data\250528_POLYA_AA_A.raw"
# polya_raw_data_location = r"C:\Users\Austin\Documents\Grad School\CIU comparison paper\Synapt raw data\251201 BSA\PolyA\251201_POLYA_FULL_AA.raw"
# polya_raw_data_location = r"C:\Users\Austin\Documents\Grad School\One_click_CCS_calibration\Polya_raw_data\251001_POLA_FULL_AA.raw"
polya_raw_data_location = r"C:\Users\Austin\Documents\Grad School\One_click_CCS_calibration\Polya_raw_data\251201_POLYA_FULL_AB.raw"
polya_raw_data_location = r"D:\260206 New BSA\260206_POLYA.raw"
polya_rangefiles_location = r"Polya_range_files\*.txt" 
polya_out_folder = '\\' + polya_raw_data_location.split('\\')[-1][:-4] #dont change this one
polya_outputfile_location = os.getcwd() + '\\Polya_raw_data' + polya_out_folder
if not os.path.exists(polya_outputfile_location):        
    os.makedirs(polya_outputfile_location)
polya_expected_peaks = 1
polya_num_cores = 1     #currently, this needs to be one, any higher and it breaks
polya_dt_tolerance = 4

#User inputs for creating the IMSCal_19 calibration input file
calibration_file_name = 'One_Click_CCS_test.dat'
cal_velocity = 500
cal_voltage = 30
cal_pressure = 1.77
cal_accuracy = 2

#User inputs for creating CIUSuite3 input files
# CUISuite3_input_name = 'BSA_one_click_CCS_input.csv'
# raw_CIU_data = glob.glob(r"C:\Users\Austin\Documents\Grad School\One_click_CCS_calibration\CIU_raw_data\*.raw")
# output_csvs_location = os.getcwd() + '\\Extracted_csvs'
# multiplex = False

CUISuite3_input_name = '260206_BSA_17_CCS_extraction_HALF1.csv'
raw_CIU_data = glob.glob(r"D:\260206 New BSA\half\*.raw")
output_csvs_location = os.getcwd() + '\\BSA_17_CCS_new'
multiplex = True
if not os.path.exists(polya_outputfile_location):        
    os.makedirs(polya_outputfile_location)
trips = ['A']
proteins = ['BSA']
charges = ['17']
# range_starts = [4400, 4120, 3875]
# range_ends = [4550, 4230, 3975]
range_starts = [3875]
range_ends = [3975]
masses = [66465]

#Base TWIMExtract parameters that don't need changing
polya_extraction_mode = 1
polya_rulemode_boolean = False
polya_combine_boolean = True
polya_ms_boolean = True
polya_check_list = glob.glob(polya_outputfile_location + '\\*.csv')

###################################################################################33

print('Starting the poly-alanine extraction:')
start_polya = time.time()
for polya_range in glob.glob(polya_rangefiles_location):
    
    short_range = polya_range.split('\\')[-1]
    cmd_arg = f'java -jar "{extractor_path}" -i "{polya_raw_data_location}" -o "{polya_outputfile_location}" -m {polya_extraction_mode} -r "{polya_range}" -rulemode {polya_rulemode_boolean} -combinemode {polya_combine_boolean} -ms {polya_ms_boolean}'
    
    completed_proc = subprocess.run(cmd_arg)
    if not completed_proc.returncode == 0:
        # process finished successfully
        logger.error('Error in extraction for file {} with range file {}. Data NOT extracted. Check that this is a Waters raw data file and that appropriate range values were provided.'.format(polya_raw_data_location, polya_range))

    else:
        print(f'Done with {short_range}', flush=True)
end_polya = time.time()
polya_time = (end_polya - start_polya)
total_times.append(polya_time)
print(f'Finished the poly-alanine extraction, took {polya_time} seconds', flush=True)




# %%

# PolyA Gaussian Fitting

polya_distributions = glob.glob(polya_outputfile_location + '\\*_raw.csv')
print()
start_gauss = time.time()
print('Starting the PolyA Gaussian Fitting')

CIUsuite_class = CIU2_Main.CIUSuite2()
CIU2_Main.CIUSuite2.load_raw_files(CIUsuite_class, polya_distributions)

print('Loaded Files')


CIU2_Main.CIUSuite2.on_button_gaussfit_clicked(CIUsuite_class, polya_expected_peaks, polya_num_cores, polya_dt_tolerance)

end_gauss = time.time()
gauss_time = end_gauss - start_gauss
total_times.append(gauss_time)
print()
print(f'FInished the poly-alanine gaussian fitting, took {gauss_time} seconds', flush=True)

# %%

print('Creating the IMSCal19 input file for CIUSuite3')
start_cal = time.time()

polya_gaussian_file = r"Polya_raw_data" + polya_out_folder + '\\All_gaussians.csv'
base_calibration_file = 'base_calibration_file.txt'



with open(polya_gaussian_file, 'r') as R:                
    gaussians_file = [x.split() for x in R.readlines()]

in_centroid_data = False
data_count = 0
centroid_dict = {}
for line in gaussians_file:
    if len(line) == 0:
        in_centroid_data = True
        continue
    if in_centroid_data:
        if 'File,CV,' in line[0]:
            
            data_count += 1
            continue
        if data_count == 1:
            line = line[0].split(',')
            
            polya = line[0].split('#')[1][:-4]
            centroid = line[-1]
            centroid_dict[polya] = centroid

        else:
            continue


with open(base_calibration_file, 'r') as R:
    cal_base = [x.split() for x in R.readlines()]
    
W =  open(calibration_file_name, 'w')

for line in cal_base:
    if '#' in line:
        delim = ' '
        
        if 'velocity' in line:
            line[-1] = str(cal_velocity)
        if 'voltage' in line:
            line[-1] = str(cal_voltage)
        if 'pressure' in line:
            line[-1] = str(cal_pressure)
        if 'accuracy' in line:
            line[-1] = str(cal_accuracy)
        
    else:
        delim = '\t'

        poly = line[0]

        if 'z' in poly:
            charge = poly.split('z')[1][0]
        else:
            charge = 1

        num = poly.split('_')[1]

        if len(num) == 1:
            numa = '0' + num
        else:
            numa = num
            
        polya_key = f'A{numa}_{charge}+'
        

        line[-1] = str(centroid_dict[polya_key])

        

    write_line = delim.join(line) + '\n'
    W.write(write_line)
W.close()

end_cal = time.time()
cal_time = end_cal - start_cal
total_times.append(cal_time)
print(f'Created calibration file {calibration_file_name} in {cal_time} seconds')

# %%

print()
print('Creating the CIUSuite3 extraction .csv')

start_suite_setup = time.time()




base_setup = r"base_CIUSuite_file.csv"
calibration_file_path = os.path.abspath(calibration_file_name)


with open(base_setup, 'r') as R:
    setup_base = [x.split(',') for x in R.readlines()]

suite_line = setup_base[1]
delim = ','
S = open(CUISuite3_input_name, 'w')
S.write(delim.join(setup_base[0]))

    
    
if multiplex:
    # for i in range(len(trips)):
    for protein, charge, range_start, range_end, mass in zip(proteins, charges, range_starts, range_ends, masses):
        for raw in raw_CIU_data:
            short_raw = raw.split('\\')[-1]
        
            current_trip = short_raw.split('_')[-1][0]
            current_volt_step = short_raw.split('_')[-1][1]
            
            if current_volt_step == 'A' and current_trip in trips:
                # print()
                new_setup_line = suite_line
                
                trip_name = f'{protein}{charge}{current_trip}'
                
                new_setup_line[1] = trip_name
                new_setup_line[2] = str(range_start)
                new_setup_line[3] = str(range_end)
                new_setup_line[8] = str(mass)
                new_setup_line[9] = str(charge)
                new_setup_line[10] = output_csvs_location
                new_setup_line[11] = calibration_file_path
                new_setup_line[12] = output_csvs_location
                new_setup_line[14] = str(cal_velocity)
                new_setup_line[15] = str(cal_voltage)
                new_setup_line[16] = str(cal_pressure)
                
                
                
                write_line = delim.join(new_setup_line) + '\n'
                
                S.write(write_line)
            
            
            S.write(raw + ',,,,,,,,,,,,,,,,,' + '\n')
    S.close()
    
if not multiplex:
    for raw in raw_CIU_data:
        short_raw = raw.split('\\')[-1]
        if 'BSA' in short_raw:
            protein = 'BSA'
            mass = 66879
        
            if '_15_' in short_raw:
                charge = 15
                range_start = 4400
                range_end = 4550
            if '_16_' in short_raw:
                charge = 16
                range_start = 4120
                range_end = 4230
            if '_17_' in short_raw:
                charge = 17
                range_start = 3875
                range_end = 3975
        
        current_trip = short_raw.split('_')[-1][0]
        current_volt_step = short_raw.split('_')[-1][1]
        
        if current_volt_step == 'A' and current_trip in trips:
            # print()
            new_setup_line = suite_line
            
            trip_name = f'{protein}{charge}{current_trip}'
            
            new_setup_line[1] = trip_name
            new_setup_line[2] = str(range_start)
            new_setup_line[3] = str(range_end)
            new_setup_line[8] = str(mass)
            new_setup_line[9] = str(charge)
            new_setup_line[10] = output_csvs_location
            new_setup_line[11] = calibration_file_path
            new_setup_line[12] = output_csvs_location
            new_setup_line[14] = str(cal_velocity)
            new_setup_line[15] = str(cal_voltage)
            new_setup_line[16] = str(cal_pressure)
            
            
            
            write_line = delim.join(new_setup_line) + '\n'
            
            S.write(write_line)
        
        
        S.write(raw + ',,,,,,,,,,,,,,,,,' + '\n')
        
        
        # print(current_trip)






    S.close()

end_suite_setup = time.time()
suite_setup_time = end_suite_setup - start_suite_setup
total_times.append(suite_setup_time)
print(f'Done with CIUSuite3 setup in {suite_setup_time} seconds')


# %%


print()
print('Running CIUSuite3 CCS calibrated extraction', flush=True)
start_extract = time.time()
cwd = os.getcwd()
mycwd = os.getcwd()
suite_file = os.path.abspath(CUISuite3_input_name)

root_dir = os.path.dirname(__file__)
hard_twimextract_path = os.path.join(root_dir, 'TWIMExtract', 'jars', 'TWIMExtract.jar')

CCS_class = Raw_Data_Import.WatersImportTypeUI(hard_twimextract_path)

Raw_Data_Import.WatersImportTypeUI.on_button_ccsextraction_clicked(CCS_class, suite_file)
os.chdir(cwd)
end_extract = time.time()
extract_time = end_extract - start_extract
total_times.append(extract_time)
print(f'Done with CIUSuite3 extraction in {extract_time} seconds')









# %%


ccs_calibrated_files = glob.glob(output_csvs_location + '\\*ccs_raw.csv')
print()
start_plots = time.time()
print('Starting the plotting and feature detection of the CCS calibrated data')

CIUsuite_class = CIU2_Main.CIUSuite2()
CIU2_Main.CIUSuite2.load_raw_files(CIUsuite_class, ccs_calibrated_files)

print('Loaded Files')

CIU2_Main.CIUSuite2.on_button_oldplot_clicked(CIUsuite_class)

min_feature_length = 10
min_feature_width_tolerance = 200

CIU2_Main.CIUSuite2.on_button_feature_detect_clicked(CIUsuite_class, min_feature_length, min_feature_width_tolerance)

end_plots = time.time()

plot_time = end_plots - start_plots
total_times.append(plot_time)
print(f'Done with plots and feature detection in {plot_time} seconds')


# %%

total_time = np.sum(total_times)

print(f'Done with One Click CCS calibration in {total_time:0.2f} seconds')




















































