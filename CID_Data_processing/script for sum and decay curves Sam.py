# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 11:28:16 2023

@author: prelllab
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os

trialpathname_a='gas on V1.25 mg'
# trialpathname= r'Synapt G2-si/240205/13+/mid/*midA*.txt' #insert path to your data files
# trialpathname= r'Synapt G2-si/13+/high/*highA*.txt' #insert path to your data files

# trialpathname= r'240627_Argon_Myo9_nESI/240627_Myo9_low_C/*_C.csv' #insert path to your data files
# trialpathname= r'240814_Argon_Myo10_nESI/240812_Myo10_high_A/Spec*.csv' #insert path to your data files
# trialpathname= r'240814_Argon_Myo8_nESI/240815_Myo8_low_A/Spec*.csv' #insert path to your data files
# trialpathname= r'240815_Nitrogen_Myo8_nESI/240815_Myo8_N_high_C/Spec*.csv' #insert path to your data files
# trialpathname= r'240815_Nitrogen_Myo10_nESI/240815_Myo10_N_low_C/Spec*.csv' #insert path to your data files
# trialpathname= r'Agilent 6545xt/240130_Stx_13_high_A/*_A.csv' #insert path to your data files
# trialpathname= r'Sam Data/Agilent Mg high 1/*.csv' #insert path to your data files
trialpathname= r'Sam Data/Synapt Low C/*.txt' #insert path to your data files
voltstep=2#these should be updated tl reflect the experimental conditions
voltstart=2

globtest=sorted(glob.glob(trialpathname,recursive=False), key=os.path.getmtime)
# plt.figure()

precursor=[]
product = []
neutral_loss = []
heme = []

matrix=[]

monomer_4=[]
monomer_5=[]
monomer_3=[]
monomer_6=[]
monomer_7=[]

tetramer_6=[]
tetramer_7=[]
tetramer_8=[]
tetramer_9=[]
tetramer_10=[]

##this is for myoglobin trials, and creating a list of used voltages
# voltlist=np.random.rand(len(globtest))
# for i in range(len(voltlist)):
#       voltlist[i]=i*voltstep+voltstart
# ends here
#39=L

[41.55,42.91,42.99,42.78,42.89]
[42.7,42.8,42.7]



voltlist=[6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60]
# voltlist = np.arange(voltstart, voltstep*len(globtest)+1, voltstep)




setlength=int(len(globtest))
for i in range((setlength)):
    temp_ms=np.loadtxt(globtest[i])
    # temp_ms=np.loadtxt(globtest[i],delimiter=',')                   #Use this for Agilent Data
    
    # plt.legend()    
    # plt.xlim(2000,2500)
    # plt.xlabel('m/z')
    # plt.ylabel('abundance')
    # plt.title(trialpathname)
    matrix=np.zeros([8000])
    mzs = np.arange(0,8000,1)
    for j in range(len(temp_ms)):
        mz=int(temp_ms[j,0])
        matrix[mz]=temp_ms[j,1]
    # plt.plot(mzs,matrix)
    
    # if i == 19:
    #     plt.figure()
    #     plt.plot(temp_ms[:,0],temp_ms[:,1], label=str(15+i*5))
    #     plt.plot(mzs+1, matrix)
    #     plt.title('8+ data')
    # precursor_sum=sum(matrix[3540:3560]) #11+ precursor
    # precursor_sum=sum(matrix[3249:3269])   #12+ precursor
    # precursor_sum=sum(matrix[2998:3018])   #13+ precursor
    

    
    # begin_pre = np.where(temp_ms.T[0] <= 1957)[0][-1]       #9+
    # end_pre = np.where(temp_ms.T[0] >= 1963.5)[0][0]
    
    # begin_prod = np.where(temp_ms.T[0] <= 2124)[0][-1]
    # end_prod = np.where(temp_ms.T[0] >= 2132)[0][0]
    
    # begin_neut = np.where(temp_ms.T[0] <= 1880)[0][-1]
    # end_neut = np.where(temp_ms.T[0] >= 1900)[0][0]
    
    begin_pre = np.where(temp_ms.T[0] <= 1950)[0][-1]       #9+ Sam Data
    end_pre = np.where(temp_ms.T[0] >= 1958.5)[0][0]
    
    begin_prod = np.where(temp_ms.T[0] <= 2118)[0][-1]
    end_prod = np.where(temp_ms.T[0] >= 2125.5)[0][0]
    
    begin_neut = np.where(temp_ms.T[0] <= 1880)[0][-1]
    end_neut = np.where(temp_ms.T[0] >= 1900)[0][0]
    
    # begin_mon_3
    # end_mon_3
    
    
    # begin_pre = np.where(temp_ms.T[0] <= 1761)[0][-1]       #10+
    # end_pre = np.where(temp_ms.T[0] >= 1768)[0][0]
    
    # begin_prod = np.where(temp_ms.T[0] <= 1889)[0][-1]
    # end_prod = np.where(temp_ms.T[0] >= 1895)[0][0]
    
    # begin_neut = np.where(temp_ms.T[0] <= 1697)[0][-1]
    # end_neut = np.where(temp_ms.T[0] >= 1707)[0][0]
    
    # begin_pre = np.where(temp_ms.T[0] <= 2202.5)[0][-1]       #8+
    # end_pre = np.where(temp_ms.T[0] >= 2209)[0][0]
    
    # begin_prod = np.where(temp_ms.T[0] <= 2428)[0][-1]
    # end_prod = np.where(temp_ms.T[0] >= 2436)[0][0]
    
    # begin_neut = np.where(temp_ms.T[0] <= 2124)[0][-1]
    # end_neut = np.where(temp_ms.T[0] >= 2132)[0][0]
    
    
    begin_heme = np.where(temp_ms.T[0] <= 617)[0][-1]
    end_heme = np.where(temp_ms.T[0] >= 621)[0][0]
    
    # begin_heme = np.where(temp_ms.T[0] <= 613)[0][-1]           #Sam Data
    # end_heme = np.where(temp_ms.T[0] >= 621)[0][0]
    
    baseline = np.mean(matrix[700:1500])
    
    abundance = temp_ms.T[1] - baseline
    abundance = np.where(abundance < 0, 0, abundance)
    
    # precursor_sum=sum(matrix[1957:1960])# - baseline)   #Myo 9+ precursor
    precursor_sum = sum(temp_ms.T[1][begin_pre:end_pre] - baseline)
    precursor.append(precursor_sum)
    
    # monomer_3_sum = sum(temp.ms.T[1][begin_monomer_3:end_monomer_3])
    # monomer_3.append(monomer_3)
    
    # product_sum = sum(matrix[2125:2128])# - baseline)
    product_sum = sum(temp_ms.T[1][begin_prod:end_prod] - baseline)
    product.append(product_sum)
    
    neutral_sum = sum(temp_ms.T[1][begin_neut:end_neut] - baseline)
    neutral_loss.append(neutral_sum)
    
    heme_sum = sum(temp_ms.T[1][begin_heme:end_heme] - baseline)
    heme.append(heme_sum)
    
    
    '''
    # if i == 1:
    #     # print()
    #     plt.plot(temp_ms[:,0],temp_ms[:,1], label=str(15+i*5))
    #     plt.plot([temp_ms[:,0][0], temp_ms[:,0][-1]], [baseline, baseline])
    #     plt.show()
    monomer_3_sum = sum(matrix[2597:2617])
    monomer_4_sum = sum(matrix[1946:1966])
    monomer_5_sum = sum(matrix[1555:1575])
    monomer_6_sum = sum(matrix[1294:1314])
    monomer_7_sum = sum(matrix[1107:1127])
    
    tetramer_6_sum = sum(matrix[5202:5222])
    tetramer_7_sum = sum(matrix[4459:4479])
    tetramer_8_sum = sum(matrix[3906:3916])
    tetramer_9_sum = sum(matrix[3466:3486])
    tetramer_10_sum = sum(matrix[3114:3134])
    
    monomer_3.append(monomer_3_sum)
    monomer_4.append(monomer_4_sum)
    monomer_5.append(monomer_5_sum)
    monomer_6.append(monomer_6_sum)
    monomer_7.append(monomer_7_sum)
    tetramer_6.append(tetramer_6_sum)
    tetramer_7.append(tetramer_7_sum)
    tetramer_8.append(tetramer_8_sum)
    tetramer_9.append(tetramer_9_sum)
    tetramer_10.append(tetramer_10_sum)
    precursor.append(precursor_sum)
    '''
    
    
#to plot non-normalized data comment out the following for loop (lines 44 through 49) selet all and ctrl 1
totalions=np.zeros([len(precursor)])    
# for i in range(len(precursor)):
#     totalions[i]= precursor[i] + monomer_4[i] + monomer_5[i] + monomer_6[i] + tetramer_10[i] + tetramer_9[i] + tetramer_8[i]
    
    # precursor[i] = precursor[i]/
    # charged_no_heme[i]=charged_no_heme[i]/totalions[i] #(take the hashtag off for normalized plot)
    # neutral_no_heme[i]=neutral_no_heme[i]/totalions[i] #(take the hashtag off for normalized plot)
    # mg_with_heme[i]=mg_with_heme[i]/totalions[i] #(take the hashtag off for normalized plot)

monomer_3 = np.array(monomer_3)
monomer_4 = np.array(monomer_4)
monomer_5 = np.array(monomer_5)
monomer_6 = np.array(monomer_6)
monomer_7 = np.array(monomer_7)

tetramer_6 = np.array(tetramer_6)
tetramer_7 = np.array(tetramer_7)
tetramer_8 = np.array(tetramer_8)
tetramer_9 = np.array(tetramer_9)
tetramer_10 = np.array(tetramer_10)

totalions = precursor + tetramer_6 + tetramer_7 + tetramer_8 #for 11+
# totalions = precursor + monomer_3 + monomer_4 + monomer_5 + monomer_6 + tetramer_6 + tetramer_7 + tetramer_8 + tetramer_9 #for 12+
# totalions = precursor + 2*(monomer_3 + monomer_4 + monomer_5 + monomer_6 + monomer_7 + tetramer_6)# + tetramer_7 + tetramer_8 + tetramer_9 + tetramer_10
precursor = np.array(precursor)
product = np.array(product)
neutral_loss = np.array(neutral_loss)
heme = np.array(heme)
totalions = precursor + product + heme


n_precursor = precursor/totalions
n_product = product/totalions
n_neutral_loss = neutral_loss/totalions
n_heme = heme/totalions

tot_prod = precursor + product
tot_heme = precursor + heme*30

just_prod = precursor/tot_prod
just_heme = precursor/tot_heme

'''
precursor = precursor/totalions
monomer_3 = monomer_3/totalions
monomer_4 = monomer_4/totalions
monomer_5 = monomer_5/totalions
monomer_6 = monomer_6/totalions
monomer_7 = monomer_7/totalions
tetramer_6 = tetramer_6/totalions
tetramer_7 = tetramer_7/totalions
tetramer_8 = tetramer_8/totalions
tetramer_9 = tetramer_9/totalions
tetramer_10 = tetramer_10/totalions
'''

# fig=plt.figure
# ax = plt.subplot(111)
# plt.figure()
# plt.title('low pressure normalized',fontsize=36) #(take the hashtag off for normalized plot)
# plt.xlim([0,55])
# plt.plot(voltlist, charged_no_heme, label='charged heme loss', marker='o')
# plt.plot(voltlist, neutral_no_heme, label='neutral heme loss', marker='o')
# plt.plot(voltlist, mg_with_heme, label='precursor', marker='o')
# #plt.plot(voltlist, totalions, label='total ions')
# ax.minorticks_on()
# plt.xlabel('Wave Voltage (velocity 500m/s)',fontsize=20)
# ax.tick_params(axis='x', labelsize=18)
# ax.tick_params(axis='y', labelsize=18)
# plt.ylabel('relative abundance',fontsize=20)
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
# # Put a legend to the right of the current axis
# #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=20)
# plt.rcParams["figure.figsize"] = (20,12)
# plt.savefig('test2png.png', dpi=100)
# plt.show()

# halfmax = max(mg_with_heme)/2
# percapp = []
# for i in range(len(voltlist)): #creates a nested array of the decay curves for each trial
#     percapp.append(float(mg_with_heme[i]))
# for i in range(len(voltlist)): #This loop calculates the CID50 for each trial. It takes the max percentage in the data set, divides it
#     if mg_with_heme[i] <= halfmax: #by two and approximates where on the curve it should be. 
#         slope = (mg_with_heme[i-1]-mg_with_heme[i])/(voltlist[i-1]-voltlist[i])
#         intercept = mg_with_heme[i] - slope*voltlist[i]
#         CID50 = (halfmax - intercept)/slope
#         print(trialpathname+'CID is '+ str(CID50))



fig=plt.figure
ax = plt.subplot(111)
xlim=([0,55])
plt.plot(voltlist, n_precursor, label='both', marker='o',markersize=12, linewidth=4, color='r')
plt.plot(voltlist, just_prod, label='just product', marker='o',markersize=12, linewidth=4, color='k')
plt.plot(voltlist, just_heme, label='just heme', marker='o',markersize=12, linewidth=4, color='b')
# plt.plot(voltlist, product, label='charged heme loss', marker='o',markersize=12, linewidth=4)
# plt.plot(voltlist, neutral_loss, label='neutral heme loss', marker='o',markersize=12, linewidth=4)
# plt.plot(voltlist, heme, label='heme', marker='o',markersize=12, linewidth=4)
# plt.plot(voltlist, monomer_3, label='charged heme loss', marker='o',markersize=12, linewidth=4)
# plt.plot(voltlist, monomer_4, label='charged heme loss', marker='o',markersize=12, linewidth=4)
# plt.plot(voltlist, monomer_5, label='charged heme loss', marker='o',markersize=12, linewidth=4)
# plt.plot(voltlist, neutral_no_heme, label='neutral heme loss', marker='o',markersize=12, linewidth=4)
# plt.plot(voltlist, mg_with_heme, label='precursor', marker='o',markersize=12, linewidth=4)
#plt.plot(voltlist, totalions, label='total ions')
plt.title('Myo breakdown curve', fontsize=34)
# ax.minorticks_on()
plt.xlabel('m/z',fontsize=26)
# ax.tick_params(axis='x', labelsize=24)
ax.tick_params(axis='y', labelsize=24)
plt.ylabel('ion count',fontsize=26)
# plt.ylim(-0.1,1.1)
box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=20)
plt.rcParams["figure.figsize"] = (20,12)
# plt.savefig('Mg_far_A_1.pdf', dpi=100)


temp_var=np.column_stack((voltlist, precursor))
# temp_var= zip(voltlist,precursor)
# temp_var=np.append(temp_var, np.array(charged_no_heme), axis=0)
# temp_var=np.append(temp_var, neutral_no_heme)


# np.savetxt(
#     '6545xt_Myo_10_low_C'+'.csv', temp_var, delimiter=',')
plt.show()