# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 13:33:54 2024

@author: Austin
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

# trial_name_path = r'6545xt_Argon_Myo_9_*_*.csv' #for stx 13+
# trial_name_path = r'6545xt_Argon_Myo_9_high_with_neutral_*.csv' #for stx 12+
# trial_name_path = r'6545xt_Argon_Myo_10_high_*.csv' #for stx 12+
# trial_name_path = r'Myo_9_low_*.csv' #for stx 12+
# trial_name_path = r'Stx_13_high_*.csv' #for stx 12+
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sam IonSPA paper inputs\Full redone Breakdown CID curves\Stx_12_low_A.csv' #for stx 12+
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sam IonSPA paper inputs\Full redone Breakdown CID curves\Stx_12_CIU_high_A.csv' #for stx 12+
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sam IonSPA paper inputs\Full redone Breakdown CID curves\Myo_9_low_A.csv' #for stx 12+
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sam IonSPA paper inputs\Full redone Breakdown CID curves\Myo_9_low_CIU_A.csv' #for stx 12+

# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Full redone Breakdown CID curves\*.csv'

""
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Breakdown Curves\Myo_8_low_CIU_A*.csv'
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Breakdown Curves\6545xt_Myo_8_low_A*.csv'
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Breakdown Curves\Myo_*CIU*.csv'
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Breakdown Curves\Stx_*CIU*.csv'
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Breakdown Curves\6545xt_Stx_13_high_A*.csv'
# 
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Breakdown Curves\6545xt_Argon_Myo_9_low_*.csv'
# trial_name_path = r'C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\Breakdown Curves\Myo_9_low_CIU_A*.csv'
# trial_name_path = r"C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\ubq_bdcs_background_subtracted\Ubq_11_CIU_high_*.csv"
# trial_name_path = r"C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\ubq_bdcs_background_subtracted_old_11+\Ubq_11_CIU_high_*.csv"
# trial_name_path = r"C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\ubq_bdcs_250628\Ubq_11_CIU_high_*.csv"
# trial_name_path = r"C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\ubq_bdcs_fixed\Ubq_7_CIU_high_*.csv"
# trial_name_path = r"C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\ubq_bdcs_sliced\Ubq_7_CIU_high_*.csv"
trial_name_path = r"C:\Users\Austin\Documents\Grad School\Sams IonSPA paper data\bradykinin_everything_bdcs_left\BK_1_CIU_mid_C*.csv"
# plt.figure()
for trialname in glob.glob(trial_name_path):
    # print(trialname)
    # if '125' in trialname:
    #     continue
    # if 'CIU' not in trialname:
    #     continue
    data = np.loadtxt(trialname, delimiter=',')
    Vins = data.T[0]
    fracs = data.T[1]
    
    ibelow = np.where(fracs <= 0.5)[0][0]
    iabove = ibelow - 1
    slope = (fracs[ibelow] - fracs[iabove])/(Vins[ibelow] - Vins[iabove]) 
    intercept = fracs[iabove] - slope*Vins[iabove]
    CID50 = (0.5 - intercept)/slope
    # if '_A.csv' in trialname:
    #     Vins = Vins + 9
    print(CID50)
    tname = trialname.split('\\')[-1]
    sV = Vins[0]
    eV = Vins[-1]
    print(tname)
    # print(f'{tname}: {sV}-{eV}')
    
    # # if 'CIU' in tname:
    # if 'Myo' in tname:
        
    #     if '8' in tname:
    #         Vins = Vins*8
    #         # if 'high' in tname:
    #         #     if '_A.csv' in tname:
    #         #         Vins = Vins - 18
    #         #     else:
    #         #         Vins = Vins - 9 
    #         # if 'mid' in tname:
    #         #     Vins = Vins - 0 
    #         # if 'low' in tname:
    #         #     Vins = Vins - 9 
                
    #     if '9' in tname:
    #         Vins = Vins*9
    #         # if 'high' in tname:
    #         #     Vins = Vins - 4
    #         # if 'mid' in tname:
    #         #     Vins = Vins - 3 
    #         # if 'low' in tname:
    #         #     Vins = Vins - 4
                
    #     if '10' in tname:
    #         Vins = Vins*10
    #         # if 'high' in tname:
    #         #     Vins = Vins - 4
    #         # if 'mid' in tname:
    #         #     Vins = Vins - 4
    #         # if 'low' in tname:
    #         #     Vins = Vins - 4 
    
    # else:
    #     if '11' in tname:
    #         Vins = Vins*11
    #         # if 'high' in tname:
    #         #     if '_A.csv' in tname:
    #         #         Vins = Vins - 18
    #         #     else:
    #         #         Vins = Vins - 9 
    #         # if 'mid' in tname:
    #         #     Vins = Vins - 0 
    #         # if 'low' in tname:
    #         #     Vins = Vins - 9 
                
    #     if '12' in tname:
    #         Vins = Vins*12
    #         # if 'high' in tname:
    #         #     Vins = Vins - 4
    #         # if 'mid' in tname:
    #         #     Vins = Vins - 3 
    #         # if 'low' in tname:
    #         #     Vins = Vins - 4
                
    #     if '13' in tname:
    #         Vins = Vins*13
    #         # if 'high' in tname:
    #         #     Vins = Vins - 4
    #         # if 'mid' in tname:
    #         #     Vins = Vins - 4
    #         # if 'low' in tname:
    #         #     Vins = Vins - 4 
    # # tname = trialname.split('\\')[-1]
    # # nsV = Vins[0]
    # # neV = Vins[-1]
    
    if 'low' in tname:
        if 'CIU' in tname:    
            c = 'c'
        if 'CIU' not in tname:
            c = 'g'
    if 'mid' in tname:
        if 'CIU' in tname:    
            c = 'darkorange'
        if 'CIU' not in tname:
            c = 'b'
    if 'high' in tname:
        if 'CIU' in tname:    
            c = 'm'
        if 'CIU' not in tname:
            c = 'r'
    

    # print(c)
    # print(f'{tname}: {sV}-{eV}\t {nsV}-{neV}')
    plt.ylim(-0.05, 1.05)
    plt.plot(Vins,fracs, linewidth=1, color='r', markersize=12, marker='o')
    # plt.plot(Vins[2:],fracs, linewidth=1, color='r', markersize=12, marker='o')
    plt.plot()
    # plt.title('Breakdown curves for Argon Gas', fontsize=55)
    plt.xlabel('collision voltage',fontsize=55)
    plt.ylabel('relative abundance',fontsize=55)
    plt.minorticks_on()
    plt.xticks(fontsize=42)
    plt.yticks(fontsize=42)
    plt.legend()
    plt.show
    temp_var=np.column_stack((Vins, fracs))
    # temp_var= zip(voltlist,precursor)
    # temp_var=np.append(temp_var, np.array(charged_no_heme), axis=0)
    # temp_var=np.append(temp_var, neutral_no_heme)
    
    # save_name = f'Breakdown Curves\{tname}'

    # np.savetxt(save_name, temp_var, delimiter=',')
    
# plt.savefig('whatever_you_want.pdf', dpi=100)
# plt.show()



















































































