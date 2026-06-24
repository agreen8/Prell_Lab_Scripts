# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 07:58:19 2025

@author: Austin
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import minimize, curve_fit
import glob
import io 
import os
from pathlib import Path
import shutil
import re
import json
import copy


# file_header = 'agilent_right_bsa_1V'
# exp_folder = 'agilent_right_bsa_1V_expdata'
# file_header = 'agilent_right_bsa_3V'
# exp_folder = 'agilent_right_bsa_3V_expdata'
file_header = 'no_twave_synapt_right_bsa_3V'
exp_folder = 'no_twave_synapt_right_bsa_3V_expdata'
# file_header = 'nistmab_fab'
# exp_folder = 'nistmab_fab_expdata'

protein = 'BSA'

sb = open('base_sbatch.sbatch') 
# send_runs = io.open('send_runs.sbatch', 'a', newline='\n')
base = []




dir_path = str(Path().absolute())
fit_path = dir_path + f'\\{file_header}_fits'
in_path = dir_path + f'\\{file_header}_inputs'
run_path = dir_path + f'\\{file_header}_sbatchs'
db_path = dir_path + f'\\{file_header}_databases'
plot_path = dir_path + f'\\{file_header}_fit_plots'

if not os.path.exists(fit_path):
    os.makedirs(fit_path)
if not os.path.exists(in_path):
    os.makedirs(in_path)
if not os.path.exists(run_path):
    os.makedirs(run_path)
if not os.path.exists(db_path):
    os.makedirs(db_path)
if not os.path.exists(plot_path):
    os.makedirs(plot_path)


for line in sb:
    base.append(f'{line[:-1]}')
    # send_runs.write(line[:-1] + '\n')    
sb.close()

F = io.open(f'{file_header}_fits\send_{file_header}_fits.sbatch', 'w', newline='\n')
    
for line in base:
    if r'job-name' in line:
        line = line.replace('STX11M', f'send')
    if r'output' in line:
        line = line.replace('ionspa_out', f'send_OUT')
    if r'error' in line:
        line = line.replace('ionspa_err', f'send_ERR')
    if r'--ntasks-per-node=28':
        line = line.replace('28', '1')
    if r'--partition' in line:
        line = line.replace('compute', 'compute')
        
    # V.write(line + '\n')
    F.write(line + '\n')
F.close()

S = io.open(f'{file_header}_sbatchs\send_{file_header}_runs.sbatch', 'w', newline='\n')
for line in base:
    if r'job-name' in line:
        line = line.replace('STX11M', f'send')
    if r'output' in line:
        line = line.replace('ionspa_out', f'send_OUT')
    if r'error' in line:
        line = line.replace('ionspa_err', f'send_ERR')
    if r'--ntasks-per-node=28':
        line = line.replace('28', '1')
    if r'--partition' in line:
        line = line.replace('compute', 'compute')
        
    # V.write(line + '\n')
    S.write(line + '\n')
S.close()

with open(r'Input Files\\Myo_8_high_CIU_a.txt') as d:
    hcstruct = json.load(d)
d.close() 

F = io.open(f'{file_header}_fits\send_{file_header}_fits.sbatch', 'a', newline='\n')
S = io.open(f'{file_header}_sbatchs\send_{file_header}_runs.sbatch', 'a', newline='\n')

# for i, expdata in enumerate(glob.glob(r'agilent_right_bsa_1V_expdata\*')):
# for i, expdata in enumerate(glob.glob(r'agilent_right_bsa_3V_expdata\*')):
for i, expdata in enumerate(glob.glob(r'synapt_right_bsa_3V_expdata\*')):
# for i, expdata in enumerate(glob.glob(r'nistmab_expdata\*')):
# for i, expdata in enumerate(glob.glob(r'nistmab_fab_expdata\*')):

    csv_name = expdata.split('\\')[-1]
    # print(csv_name)
    # continue
    if 'high' in csv_name:
        pressure = 'high'
        cell = 'Synapt_CIU_high.txt'
    elif 'mid' in csv_name:
        pressure = 'mid'
        cell = 'Synapt_CIU_mid.txt'
    elif 'low' in csv_name:
        pressure = 'low'
        cell = 'Synapt_CIU_low.txt'
    elif '6560c':
        pressure = ''
        cell = 'a_CIU_Cell.json'
    if '6560c' not in csv_name:
        pressure = 'high'
        # cell = 'Synapt_CIU_high.txt'
        cell = 'Synapt_new_CIU_high.txt'
    
    if 'CIU' in csv_name:
        ciu = 'CIU'
    # else:
    #     ciu = ''
    if '_A.' in csv_name:
        trip = 'a'
    if '_B.' in csv_name:
        trip = 'b'
    if '_C.' in csv_name:
        trip = 'c'
    if '_D.' in csv_name:
        trip = 'd'

    if 'ConA' in csv_name:
        ion = 'ConA_20_ion.txt'
        db = r'ConA_20_CIU_1_0_A'
    
    
    if '6560c' in csv_name:
        if 'BSA_15' in csv_name:
            if 'CIU_1' in csv_name:
                ion = '6560c_BSA_15_ion_1.txt'
                db = '6560c_BSA_15_1'
            if 'CIU_2' in csv_name:
                ion = '6560c_BSA_15_ion_2.txt'
                db = '6560c_BSA_15_2'
            if 'CIU_3' in csv_name:
                ion = '6560c_BSA_15_ion_3.txt'
                db = '6560c_BSA_15_3'
                
        if 'BSA_16' in csv_name:
            if 'CIU_1' in csv_name:
                ion = '6560c_BSA_16_ion_1.txt'
                db = '6560c_BSA_16_1'
            if 'CIU_2' in csv_name:
                ion = '6560c_BSA_16_ion_2.txt'
                db = '6560c_BSA_16_2'
            if 'CIU_3' in csv_name:
                ion = '6560c_BSA_16_ion_3.txt'
                db = '6560c_BSA_16_3'
            
        if 'BSA_17' in csv_name:
            if 'CIU_1' in csv_name:
                ion = '6560c_BSA_17_ion_1.txt'
                db = '6560c_BSA_17_1'
            if 'CIU_2' in csv_name:
                ion = '6560c_BSA_17_ion_2.txt'
                db = '6560c_BSA_17_2'
            if 'CIU_3' in csv_name:
                ion = '6560c_BSA_17_ion_3.txt'
                db = '6560c_BSA_17_3'
            if 'CIU_4' in csv_name:
                ion = '6560c_BSA_17_ion_4.txt'
                db = '6560c_BSA_17_4'
            
    
    if '6560c' not in csv_name:
        if 'BSA_15' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'BSA_15_ion_1.txt'
                db = 'BSA_15_1'
            if 'CIU_2' in csv_name:
                ion = 'BSA_15_ion_2.txt'
                db = 'BSA_15_2'
            if 'CIU_3' in csv_name:
                ion = 'BSA_15_ion_3.txt'
                db = 'BSA_15_3'
                
        if 'BSA_16' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'BSA_16_ion_1.txt'
                db = 'BSA_16_1'
            if 'CIU_2' in csv_name:
                ion = 'BSA_16_ion_2.txt'
                db = 'BSA_16_2'
            if 'CIU_3' in csv_name:
                ion = 'BSA_16_ion_3.txt'
                db = 'BSA_16_3'
            
        if 'BSA_17' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'BSA_17_ion_1.txt'
                db = 'BSA_17_1'
            if 'CIU_2' in csv_name:
                ion = 'BSA_17_ion_2.txt'
                db = 'BSA_17_2'
            if 'CIU_3' in csv_name:
                ion = 'BSA_17_ion_3.txt'
                db = 'BSA_17_3'
            if 'CIU_4' in csv_name:
                ion = 'BSA_17_ion_4.txt'
                db = 'BSA_17_4'
       
        if 'NistmAb_13' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_13_ion_1.txt'
                db = 'NistmAb_13_1'
        if 'NistmAb_13' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_13_ion_2.txt'
                db = 'NistmAb_13_2'    
        if 'NistmAb_14' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_14_ion_1.txt'
                db = 'NistmAb_14_1'
        if 'NistmAb_14' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_14_ion_2.txt'
                db = 'NistmAb_14_2'    
        if 'NistmAb_19' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_19_ion_1.txt'
                db = 'NistmAb_19_1'
        if 'NistmAb_19' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_19_ion_2.txt'
                db = 'NistmAb_19_2'    
        if 'NistmAb_20' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_20_ion_1.txt'
                db = 'NistmAb_20_1'
        if 'NistmAb_13' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_20_ion_2.txt'
                db = 'NistmAb_20_2'    
        
        if 'NistmAb_24' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_24_ion_1.txt'
                db = 'NistmAb_24_1'
        if 'NistmAb_24' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_24_ion_2.txt'
                db = 'NistmAb_24_2'
        
        if 'NistmAb_25' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_25_ion_1.txt'
                db = 'NistmAb_25_1'
        if 'NistmAb_25' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_25_ion_2.txt'
                db = 'NistmAb_25_2'
                
        if 'NistmAb_26' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_26_ion_1.txt'
                db = 'NistmAb_26_1'
        if 'NistmAb_26' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_26_ion_2.txt'
                db = 'NistmAb_26_2'
                
        if 'NistmAb_27' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_27_ion_1.txt'
                db = 'NistmAb_27_1'
        if 'NistmAb_27' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_27_ion_2.txt'
                db = 'NistmAb_27_2'
                
        if 'NistmAb_28' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_28_ion_1.txt'
                db = 'NistmAb_28_1'
        if 'NistmAb_28' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_28_ion_2.txt'
                db = 'NistmAb_28_2'
                
        if 'NistmAb_29' in csv_name:
            if 'CIU_1' in csv_name:
                ion = 'NistmAb_29_ion_1.txt'
                db = 'NistmAb_29_1'
        if 'NistmAb_29' in csv_name:
            if 'CIU_2' in csv_name:
                ion = 'NistmAb_29_ion_2.txt'
                db = 'NistmAb_29_2'
        
        
                
    idic = copy.deepcopy(hcstruct)
    
    idic['cellfilename'] = f'{file_header}_inputs/{cell}'
    idic['exp']['folder'] = f'{exp_folder}'
    idic['exp']['file'] = csv_name
    idic['ionfilename'] = f'{file_header}_inputs/{ion}'
    idic['dhds_start'] = [80, -100]
    idic['tmax'] = 0.005
    
    if '6560c' in csv_name:
        idic['Nrepeats'] = 50
    if '6560c' not in csv_name:
        idic['Nrepeats'] = 50
    
    iname = csv_name[:-5] + f'{trip}.txt'
    input_name = in_path + f'\\{csv_name[:-5]}{trip}.txt'
    
    with open(input_name, 'w') as ijson_file:
        json.dump(idic, ijson_file, indent=4)
    ijson_file.close()
    
    # print(f'Expirement file is: {csv_name}')
    # print(f'Created input file: {iname}')
    # print()

    '''
    Now creating the fit sbatches
    '''
    
    fit_name = csv_name[:-5] + f'{trip}_fit.sbatch'
    
    K = io.open(f'{fit_path}\\{fit_name}', 'w', newline='\n')
    for line in base:
        if r'job-name' in line:
            line = line.replace('STX11M', f'{csv_name[:-5]}{trip}')
        if r'output' in line:
            line = line.replace('ionspa_out', f'{csv_name[:-5]}{trip}_OUT')
        if r'error' in line:
            line = line.replace('ionspa_err', f'{csv_name[:-5]}{trip}_ERR')
        if r'--ntasks-per-node=28':
            line = line.replace('28', '1')
        if r'--partition' in line:
            line = line.replace('compute', 'compute')
        K.write(line + '\n')
    K.write('cd ../' + '\n')
    K.close()
    
    K = io.open(f'{fit_path}\\{fit_name}', 'a', newline='\n')
    
    # fit_line = f'python spa_db_fit.py --db {file_header}_databases/{db}.sqlite3 {file_header}_inputs/{iname} --outfile {file_header}_fits.txt --plotfile None --plotmap None --tpoints 1000 --weighting combined --updatedb False'
    fit_line = f'python spa_db_fit.py --db {file_header}_databases/{db}.sqlite3 {file_header}_inputs/{iname} --outfile {file_header}_fits.txt --plotfile {file_header}_fit_plots/{iname[:-4]}.png --plotmap none --tpoints 1000 --weighting combined --updatedb False'
    # fit_line = f'python spa_db_new_fit.py --db {file_header}_databases/{db}.sqlite3 {file_header}_inputs/{iname} --outfile {file_header}_fits.txt --plotfile {file_header}_fit_plots/{iname[:-4]}.png --plotmap none --tpoints 1000 --weighting combined --updatedb False'
    
    fit_line = f'python spa_avg_fit.py --db {file_header}_databases/{db}.sqlite3 {file_header}_inputs/{iname} --outfile {file_header}_avg_fits.txt --plotfile {file_header}_fit_plots/{iname[:-4]}.png --plotmap none --tpoints 1000 --weighting combined'
    print(fit_line)
    
    K.write(fit_line + '\n')
    K.close()
    
    send_fit_line = f'sbatch {fit_name}'
    F.write(f'{send_fit_line}\n')
    
    
    
    
    # print(f'Created fit sbatch with name: {fit_name}')
    
    
    '''
    Now creating the send sbatches
    '''
    
    if trip != 'c':
        # print()
        continue

    run_name = csv_name[:-5] + f'{trip}.sbatch'
    
    R = io.open(f'{run_path}\\{run_name}', 'w', newline='\n')
    for line in base:
        if r'job-name' in line:
            line = line.replace('STX11M', f'{csv_name[:-5]}{trip}')
        if r'output' in line:
            line = line.replace('ionspa_out', f'{csv_name[:-5]}{trip}_OUT')
        if r'error' in line:
            line = line.replace('ionspa_err', f'{csv_name[:-5]}{trip}_ERR')
        if r'--ntasks-per-node=28':
            if '6560c' in csv_name:
                line = line.replace('28', '51')
            if '6560c' not in csv_name:
                line = line.replace('28', '51')
        if r'--partition' in line:
            line = line.replace('compute', 'compute')
        R.write(line + '\n')
    R.write('cd ../' + '\n')
    R.close()
    
    R = io.open(f'{run_path}\\{run_name}', 'a', newline='\n')
    
    if '6560c' in csv_name:
        send_line = f'mpiexec -np 51 python spa_mpirun.py --db {file_header}_databases/{db}.sqlite3 {file_header}_inputs/{iname}'
    if '6560c' not in csv_name:
        send_line = f'mpiexec -np 51 python spa_mpirun.py --db {file_header}_databases/{db}.sqlite3 {file_header}_inputs/{iname}'
    R.write(send_line + '\n')
    R.close()
    
    send_run_line = f'sbatch {run_name}'
    S.write(f'{send_run_line}\n')
    
    
    # print(f'Created run sbatch with name: {run_name}')
    # print()
    
    
    # break
    
F.close()
S.close()




















































































