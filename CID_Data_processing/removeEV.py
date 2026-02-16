# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 11:24:50 2023

@author: Austin
"""


import glob

# "C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\UBQ for sams paper\250709\11+\high\250709_Ubq_11_CIU_high_AA.raw"
# for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Raw Insturment Data\Evan Myo CIU\8+\high\240918_Myo_8_high_CIU_AA.raw\_extern.inf'):
# for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\UBQ for sams paper\7+\*\*.raw\_extern.inf'):
# for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\UBQ for sams paper\250701\*\*.raw\_extern.inf'):
# for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\UBQ for sams paper\250709\11+\*\*.raw\_extern.inf'):
# for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Bradykinin\1+\*\_extern.inf'):
for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Bradykinin\250825\1+\mid\*\*_extern.inf'):
    print(name)
    filelist = []
    
    with open(name,'r') as R:
        for line in R:
            if 'Trap Collision Energy (eV)' in line:
                line = line.replace('Trap Collision Energy (eV)','Trap Collision Energy     ')
            if 'Transfer Collision Energy (eV)' in line:
                line = line.replace('Transfer Collision Energy (eV)','Transfer Collision Energy     ')
            filelist.append(line)
    R.close()
    W = open(name,'w')
    for line in filelist:
        W.write(line)
    W.close()














