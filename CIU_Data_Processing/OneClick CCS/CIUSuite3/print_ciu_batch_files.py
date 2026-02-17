# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 11:31:47 2025

@author: Austin
"""
import string
import glob
alph = list(string.ascii_uppercase)
# path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Raw Insturment Data\Stx CIU experiments\*\*\*.raw'
# path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Raw Insturment Data\Evan Myo CIU\*\*\*.raw'
# path = r'D:\UBQ for sams paper\*\*\*.raw'
# path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\UBQ for sams paper\250709\11+\*\*.raw'
# path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\UBQ for sams paper\7+\*\*.raw'
path = r"D:\Dhiman Myoglobin\Raw Data\*.raw"
num_files = len(glob.glob(path))


current_trip = 'A'
for path in glob.glob(path):
    if '241118' in path:
        continue
    trip = path[-6]
    
    if trip != current_trip:
        charge = path.split('\\')[-3][:-1]
        pressure = path.split('\\')[-2]
        if charge == '11':
            mz_low = 3504
            mz_high = 3604
        if charge == '12':
            mz_low = 3201
            mz_high = 3301
        if charge == '13':
            mz_low = 2958
            mz_high = 3058
        # print(f'*\tstx{charge}_{pressure}_{trip}\t{mz_low}\t{mz_high}\t0\t100\t1\t200\t39104\t{charge}\tC:\\Users\\Austin\\Documents\\Grad School\\Sams IonSPA paper data\\CIU Raw DT CSVs\tC:\\path_to_where_to_put_calibratedfiles\t0.245\t300\t30\t4.2\t298')
        print()
        current_trip = trip
        # break
    
    print(path)
    
    
    # break
    # print(i)
