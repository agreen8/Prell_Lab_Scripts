# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 11:24:50 2023

@author: Austin
"""


import glob



# for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Raw Insturment Data\Stx CIU experiments\*\*\*\_extern.inf'):
# for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Raw Insturment Data\Evan Myo CIU\*\*\*\_extern.inf'):
# for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\UBQ for sams paper\250709\11+\*\*\_extern.inf'):
for name in glob.glob(r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\UBQ for sams paper\7+\*\*\_extern.inf'):
    print(name)
    filelist = []
    
    with open(name,'r') as R:
        for line in R:
            if 'Trap Collision Energy     ' in line.split('\t')[0]:
                # print(line.split('\t'))
                
                line = line.replace('Trap Collision Energy     ','Trap Collision Energy (eV)')
            if 'Transfer Collision Energy     ' == line.split('\t')[0]:
                line = line.replace('Transfer Collision Energy     ','Transfer Collision Energy (eV)')
            filelist.append(line)
    R.close()
    W = open(name,'w')
    for line in filelist:
        W.write(line)
    W.close()














