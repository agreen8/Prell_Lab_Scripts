# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 08:23:02 2026

@author: Austin
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# from scipy.optimize import minimize, curve_fit
# import glob
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
import yaml


dir_path = str(Path().absolute())
out_path = dir_path + r'\260220_Everything_right_Breakdown_curves'

if not os.path.exists(out_path):        #Creates Folder for Breakdown curve output if one isn't already created
    os.makedirs(out_path)


#these set the parameters of the breakdown curve generator
with open(f"{dir_path}\\breakdown_curve_settings.yaml", 'r') as file:
    Breakdown_Curve_Settings =  yaml.safe_load(file)


protein_code = "6560c_BSA" #Eventually, we can have it manually get the code from the filename
# protein_code = "ConA" #Eventually, we can have it manually get the code from the filename

protein_dict = Breakdown_Curve_Settings[f"{protein_code}"]
vin_interp = protein_dict["vin_interp"]
dt_interp_factor = protein_dict["dt_interp_factor"]
savgol_smooth_window = protein_dict["savgol_smooth_window"]
savgol_smooth_order = protein_dict["savgol_smooth_order"]
choice_val = protein_dict["choice_val"]
low_dt_bound = protein_dict["low_dt_bound"]
high_dt_bound = protein_dict["high_dt_bound"]
other_dt_bound = protein_dict["other_dt_bound"]
amp_low = protein_dict["amp_low"]
amp_high = protein_dict["amp_high"]
std_low = protein_dict["std_low"]
std_high = protein_dict["std_high"] #set to 0.9ish to underaccount/fit both the gaussian and the noise, 1.1/1.2 for the reverse
after_vin = protein_dict["after_vin"] #set to cut off breakdowncurve if fragmentation is seen-- tune for first 1-2 transitions
plot_low = protein_dict["plot_low"]
plot_high = protein_dict["plot_high"]
gaussian_file = protein_dict["gaussian_file"]
feature_file = protein_dict["feature_file"]
subtract = True #always true
ccs = True

# these next two functions are straight from CIUSuite for help plotting the fingerprints
def sgolay2d(z, window_size, order, derivative=None):
    """
    ADAPTED FROM THE SCIPY COOKBOOK, at http://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html
    accessed 2/14/2018. Performs a 2D Savitzky-Golay smooth on the the provided 2D array z using parameters
    window_size and polynomial order.
    :param z: 2D numpy array of data to smooth
    :param window_size: filter size (int), must be odd to use for smoothing
    :param order: polynomial order (int) to use
    :param derivative: optional (string), values = row, col, both, or None
    :return: smoothed 2D numpy array
    """

    # number of terms in the polynomial expression
    n_terms = (order + 1) * (order + 2) / 2.0

    if window_size % 2 == 0:
        logger.error('2D SG smooth window_size ({}) must be odd'.format(window_size))
        raise ValueError('window_size must be odd')

    if window_size**2 < n_terms:
        logger.error('2D SG smooth order is too high for the window size of {}'.format(window_size))
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial.
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ...
    # this line gives a list of two item tuple. Each tuple contains
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [(k-n, n) for k in range(order + 1) for n in range(k+1)]

    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat(ind, window_size)
    dy = np.tile(ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    a_mat = np.empty((window_size ** 2, len(exps)))
    for i, exp in enumerate(exps):
        a_mat[:, i] = (dx ** exp[0]) * (dy ** exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    z_mat = np.zeros(new_shape)
    # top band
    band = z[0, :]
    z_mat[:half_size, half_size:-half_size] = band - np.abs(np.flipud(z[1:half_size+1, :]) - band)
    # bottom band
    band = z[-1, :]
    z_mat[-half_size:, half_size:-half_size] = band + np.abs(np.flipud(z[-half_size-1:-1, :]) - band)
    # left band
    band = np.tile(z[:, 0].reshape(-1, 1), [1, half_size])
    z_mat[half_size:-half_size, :half_size] = band - np.abs(np.fliplr(z[:, 1:half_size+1]) - band)
    # right band
    band = np.tile(z[:, -1].reshape(-1, 1), [1, half_size])
    z_mat[half_size:-half_size, -half_size:] = band + np.abs(np.fliplr(z[:, -half_size-1:-1]) - band)
    # central band
    z_mat[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0, 0]
    z_mat[:half_size, :half_size] = band - np.abs(np.flipud(np.fliplr(z[1:half_size+1, 1:half_size+1])) - band)
    # bottom right corner
    band = z[-1, -1]
    z_mat[-half_size:, -half_size:] = band + np.abs(np.flipud(np.fliplr(z[-half_size-1:-1, -half_size-1:-1])) - band)
    # top right corner
    band = z_mat[half_size, -half_size:]
    z_mat[:half_size, -half_size:] = band - np.abs(np.flipud(z_mat[half_size+1:2*half_size+1, -half_size:]) - band)
    # bottom left corner
    band = z_mat[-half_size:, half_size].reshape(-1, 1)
    z_mat[-half_size:, :half_size] = band - np.abs(np.fliplr(z_mat[-half_size:, half_size+1:2*half_size+1]) - band)

    # solve system and convolve
    if derivative is None:
        m = np.linalg.pinv(a_mat)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(z_mat, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(a_mat)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(z_mat, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(a_mat)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(z_mat, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(a_mat)[1].reshape((window_size, -1))
        r = np.linalg.pinv(a_mat)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(z_mat, -r, mode='valid'), scipy.signal.fftconvolve(z_mat, -c, mode='valid')

def get_contour_levels(ciu_data, merge_cutoff=10, num_contours=100):
    """
    Generates contours for CIU plots with the bottom 10% of data merged into a single contour for
    easier post-processing. Detects min/max values to ensure that contours match the data.
    :param ciu_data: analysis_obj.ciu_data - 2D numpy array of DT/CV ciu data
    :param merge_cutoff: Percent (int) value below which to merge all contours together. Default 10% for data scaled to 100
    :param num_contours: approximate number of contour levels to generate. Default 100
    :return: list of ints (contour levels). Pass returned list directly to pyplot.contourf as 'levels' arg
    """
    possible_steps = np.asarray([0.001, 0.01, 0.1, 1, 10, 100])

    max_val = int(round(np.max(ciu_data) * 100)) + 1  # +1 and -1 for max/min vals are to prevent white spots from appearing after rounding
    min_val = int(round(np.min(ciu_data) * 100)) - 1
    
    step = (max_val - min_val) / float(num_contours)      # round magnitude step to get close to 100 contours/bins total
    step = possible_steps[(np.abs(possible_steps - step)).argmin()]

    # ensure the min value in the data is above the merge cutoff
    if merge_cutoff < min_val:
        merge_cutoff = min_val
        levels = [x for x in np.arange(merge_cutoff, max_val, step)]
    else:
        levels = [x for x in np.arange(merge_cutoff, max_val, step)]
        levels.insert(0, min_val)
    levels = [x / 100.0 for x in levels]  # convert from integers (percent) back to float (relative intensity)
    return levels

"""Parsing Features File"""
feature_ranges = {}     #these next couple of variable are just to create them to analyze the gaussian and feature file


# feature_dts = []
file_names = []

peaks_by_file_dict = {}

file_being_processed = ""


gaussian_data = {}


#reads in the CIUSuite all gaussians .csv output file
with open(gaussian_file, 'r') as R:                
    gaussians_file = [x.split() for x in R.readlines()]


#reads in the CIUSuite all feature .csv output file
with open(feature_file, 'r') as B:                
    feature_file = [x.split() for x in B.readlines()]

collect = False
gd = []
file_to_process = ''

# reads the guassian CIUSuite data file and gets all the raw .csv file names to batch process
# will also read the feature file and get the number of featers and their median driftimes and create the general dictionary structure to store everyting

for line in gaussians_file:
    
    if len(line) != 0:
        if '#' in line[0]:
            # print(line)
            if 'DT' in line[1] or r'BSA' in line[1]:       
                if collect:
                    gaussian_data[file_to_process]['Gaussian Data'] = gd
                    gd = []
                file_to_process = line[1][:-5]
                print(file_to_process)
                old_file_to_process = file_to_process
                gaussian_data[file_to_process] = {}
                gaussian_data[file_to_process]['Data File'] = f'{file_to_process}_raw.csv'
                for fline in feature_file:
                    if file_to_process in fline[0]:
                        gaussian_data[file_to_process]['Number of Features'] = len(fline[0].split(',')) - 1
                        # print(fline[0].split(','))
                        for i in range(len(fline[0].split(',')) - 1):
                            fl = fline[0].split(',')
                            gaussian_data[file_to_process][f'Feature {i+1}'] = {}
                            gaussian_data[file_to_process][f'Feature {i+1}']['Median Centroid'] = float(fl[i+1])
                        break
                
                        
                # print(file_to_process)
                file_names.append(file_to_process)
                
                
                
    
        else:
            collect = True
            if collect:
                gd.append(line)
        

    else:
        if line == []:
            print(line)
            if collect:
                gaussian_data[file_to_process]['Gaussian Data'] = gd
                gd = []
            break
# fit class to fit the CIU data with gaussians
class fitClass:
    def __init__(self):
        pass
    def gauss_func(self, x, amp, sigma):
        
        return amp*np.exp(-(x - self.centroid)**2/(2*sigma**2))
# usefull function to plot the guassian, just so we don't have to keep creating a guassian everytime we want to plot
def plot_gauss_func(x, amp, centroid, sigma):
    
    return amp*np.exp(-(x - centroid)**2/(2*sigma**2))


# begins the batch processing
for key1 in gaussian_data.keys():
    print(key1)
    # if '24b' in key1:
        # break
    Vins = []    
    centroids = []
    fwhms = []
    
    num_feats = gaussian_data[key1]['Number of Features']
    
    gaussian_data[key1]['Total Abundance'] = []
    
    feat_range = False
    for fline in feature_file:
        if len(fline) == 0:
            continue
        if 'CV' in fline:           #this bit gets the feature range as CIUSuite was able to detect it, Essentialy it looks for key words in each line in the feature file
            feat_range = True       #to look for and then will grab the feature range
        if feat_range:
            if key1 in fline[0]:
                split = fline[0].split(',')
                for i in range(num_feats):
                    feat_beg = float(split[i+1].split('-')[0])
                    feat_end = float(split[i+1].split('-')[1])
                    gaussian_data[key1][f'Feature {i+1}']['Feature Begin'] = feat_beg
                    gaussian_data[key1][f'Feature {i+1}']['Feature End'] = feat_end

    for i in range(num_feats):
        cent = gaussian_data[key1][f'Feature {i+1}']['Median Centroid']
        gaussian_data[key1][f'Feature {i+1}']['Centroids'] = []
        gaussian_data[key1][f'Feature {i+1}']['FWHMs'] = []
        gaussian_data[key1][f'Feature {i+1}']['BDC'] = []
        gaussian_data[key1][f'Feature {i+1}']['AUC'] = []
        gaussian_data[key1][f'Feature {i+1}']['Non gaussian'] = np.array([0])
        gaussian_data[key1][f'Feature {i+1}']['GAUSSIAN_AUC'] = []
        gaussian_data[key1][f'Feature {i+1}']['LEFTOVER_AUC'] = []
        gaussian_data[key1][f'Feature {i+1}']['RIGHT_AUC'] = []
        centroids.append(cent)
    centroids = np.array(centroids) #so this is a bit of legacy code, where i tried to get the exact centroids that CIUSuite would fit but I've switched it to always use the average centroid that CIUSuite detects

    for row in gaussian_data[key1]['Gaussian Data']: #Gets the experimental Vins
        # print(row)
        split_row = row[0].split(',')
        # print(split_row)
        Vins.append(float(split_row[0]))
    
    gaussian_data[key1]['Vins'] = Vins

    for row in gaussian_data[key1]['Gaussian Data']:    #Starting to analyze the Gaussian data to get the features centroid and FWHM
        split_row = row[0].split(',')
        num_feat_vin = int(len(split_row)/4)
        vin = float(split_row[0])
        temp_centroids = []
        temp_fwhms = []
        temp_features = []
        
        array_feats = np.arange(0, num_feats, 1)
        
        
        
        for j in range(num_feat_vin):
            if j+1 > num_feats:         #I can't remember why this is here, probably to catch run away code or not look for 1 more feature than necessary
                continue
            exp_centroid = float(split_row[j*4 + 2])
            exp_fwhm = float(split_row[j*4 + 3])
            index = abs(centroids - exp_centroid).argmin()  #Finds the feature who's median drift time is closest to the experimental value
            
            # if vin >= 160 and index <= num_feats - 2:   #this doesn't catch all edge cases and may need to be adjusted on some occasions
                # continue
            
            if 'SF6' in key1 or '_N2_ ' in key1 or ccs:         #just catches a couple of edgecases for the 6560c BSA data
                if exp_centroid > 12000 or exp_centroid < 3000:
                    continue
                if exp_fwhm < 100:
                    continue
                
            if abs(exp_centroid - centroids[index])/centroids[index] > 0.025: #if the experiemtnal centroid is over 5% larger than the median value just use the median
                exp_centroid = centroids[index]
                continue
            temp_centroids.append(exp_centroid)
            temp_fwhms.append(exp_fwhm)
            temp_features.append(index)
           
        for cent, fwhm, indx in zip(temp_centroids, temp_fwhms, temp_features):     #now assigns the exp centroids into the correct feature in the overall dictionary
            
            gaussian_data[key1][f'Feature {indx+1}']['Centroids'].append(cent)
            gaussian_data[key1][f'Feature {indx+1}']['FWHMs'].append(fwhm)
            array_feats = np.delete(array_feats, np.where(array_feats == indx))
        
        for leftover in array_feats:                                                #for a specific Vin, if a exp drift time wasn't assigned to a feature, that feature then gets assigned it's median drift time
            gaussian_data[key1][f'Feature {leftover+1}']['Centroids'].append(centroids[leftover])
            gaussian_data[key1][f'Feature {leftover+1}']['FWHMs'].append(0) #I don't like this being 0, need to think of something better if there isn't an exp FWHM value

    for i in range(num_feats):  #ok this sets all those zeros to the median FWHM
        
        fwhms = np.array(gaussian_data[key1][f'Feature {i+1}']['FWHMs'])
        nonzero_fwhms = fwhms[fwhms != 0]
        median_fwhm = np.median(nonzero_fwhms)
        fwhms[fwhms == 0] = median_fwhm
        gaussian_data[key1][f'Feature {i+1}']['FWHMs'] = fwhms
        gaussian_data[key1][f'Feature {i+1}']['Median FWHM'] = median_fwhm

    if ccs and 'N2' not in key1:    
        abund_data = np.loadtxt(f'{key1}ccs_raw.csv', dtype='str', delimiter=',')      #loads the _raw.csv data file and prepares it for analysis
    elif 'N2' in key1:
        abund_data = np.loadtxt(f'{key1}_raw.csv', dtype='str', delimiter=',')      #loads the _raw.csv data file and prepares it for analysis
    else:
        abund_data = np.loadtxt(f'{key1}_raw.csv', dtype='str', delimiter=',')      #loads the _raw.csv data file and prepares it for analysis
    abund_data[0][0] = 0.0
    abund_data = abund_data.astype('float')
    stripped_abund = abund_data[1:,1:]
    exp_vins = abund_data[0][1:]
    exp_drift_times = abund_data[1:,0]
    
    gauss_drift_times = np.linspace(exp_drift_times.min(),exp_drift_times.max(), len(exp_drift_times)*dt_interp_factor)    #interpolates the drift time axis (x-axis)
    
    # f = interpolate.interp2d(exp_vins, exp_drift_times, stripped_abund, kind='cubic')
    # f = interpolate.RegularGridInterpolator((exp_vins, exp_drift_times), stripped_abund.T)
    f = interpolate.RectBivariateSpline(exp_vins, exp_drift_times, stripped_abund.T)

    new_Vins = np.linspace(min(exp_vins), max(exp_vins), vin_interp)                            #This is the interpolation of Vins, since most CIU experiments are from 5 to 200 V in 5 V steps, 40 is no interpolation, 79 is a factor of 2, ie 2.5 V step sizes. For the purposes of CIU in IonSPA I wouldn't go smaller than a 2.5 V step size
    new_drifttimes = np.linspace(exp_drift_times.min(),exp_drift_times.max(), len(exp_drift_times)*dt_interp_factor)       #This is the interpolation of the drifttimes, I've found higher numbers to be better, may need tuning for different situations
    new_xx, new_yy = np.meshgrid(new_Vins, new_drifttimes)
    new_abundances = f(new_Vins, new_drifttimes)
    
    smoothed_new_abundances = sgolay2d(new_abundances, savgol_smooth_window, savgol_smooth_order)       #smooths the data using the functions in CIUSUite, so use the same smoothing settings as you did with CIUSuite
    
    plot_new_abundances = sgolay2d(stripped_abund, 5, 2)
    
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink']
    
    
    pp = PdfPages(f"{key1}_plots_EVERYTHING_TO_THE_RIGHT_UE5.pdf")  #creates the PDF that will store all the gaussian fits and the BDCS
    
    inst = fitClass()

    num_feats = num_feats - 1
    
    # Now onto the Fitting
    for v in range(len(new_Vins)):
        interp_vin = v*(new_Vins[1] - new_Vins[0]) + exp_vins[0]
        # if round(interp_vin) != 180:
        #     continue
        print(f'Vin: {interp_vin}')
        
        smoothed_abund = smoothed_new_abundances[v]
        smoothed_abund = np.where(smoothed_abund <= 0, 1, smoothed_abund)   # setting all abundances to be positive
        smoothed_centroid = new_drifttimes[np.where(smoothed_abund == max(smoothed_abund))[0][0]]
        fit_abund = []
        fit_drift = []
        
        
        fig1 = plt.figure(layout='constrained')
        # fig2, ax2 = plt.subplots()
        if num_feats > 2 and num_feats < 7:
            num_rows = 2
        elif num_feats >= 7:
            num_rows = 3
        else:
            num_rows = 1
        
        num_colls = num_feats - num_rows - 1        
        while num_rows*num_colls < num_feats:
            num_colls += 1
        
        
        for i in range(num_feats-1, -1, -1):    #this fits the second to last feature first, and then descends to the first feature (in general I've found this reverse order to be better)
            loc = int(f'{num_rows}{num_colls}{i+1}')
            
            # interp_centroids = np.interp(new_Vins, exp_vins, gaussian_data[key1][f'Feature {i+1}']['Centroids'])
            # interp_fwhms = np.interp(new_Vins, exp_vins, gaussian_data[key1][f'Feature {i+1}']['FWHMs'])
            ax = fig1.add_subplot(loc)
            
            fwhm = gaussian_data[key1][f'Feature {i+1}']['Median FWHM']
            # fwhm = interp_fwhms[v]
            std = fwhm / (2 * np.sqrt(2*np.log(2)))
            centroid = gaussian_data[key1][f'Feature {i+1}']['Median Centroid']
            # centroid = interp_centroids[v]
            
            #these are for figuring out the bounds for the gaussian fitting
            if i == 0: # if it is the first feature use the lower bound as a percentage, the upper is a percentage of the drift time distance to the centroid of the next feature
                lower_dt = centroid*low_dt_bound
                higher_dt = abs(centroid - gaussian_data[key1][f'Feature {i+2}']['Median Centroid'])*other_dt_bound
            elif i == num_feats-1: # otherwise if it's the last feature then do the reverse 
                higher_dt = centroid*high_dt_bound
                lower_dt = abs(centroid - gaussian_data[key1][f'Feature {i}']['Median Centroid'])*other_dt_bound
            else:
                lower_dt = abs(centroid - gaussian_data[key1][f'Feature {i}']['Median Centroid'])*other_dt_bound
                higher_dt = abs(centroid - gaussian_data[key1][f'Feature {i+2}']['Median Centroid'])*other_dt_bound
            
            drift_index = abs(new_drifttimes - centroid).argmin()
            amplitude = smoothed_abund[drift_index]
            
        
            # if i == num_feats-1:
                # area_h_dt = centroid + 3*std
            # else:
            area_h_dt = (centroid + gaussian_data[key1][f'Feature {i+2}']['Median Centroid'])/2
            if i == 0:
                area_l_dt = centroid - 4*std
            else:
                area_l_dt = (centroid + gaussian_data[key1][f'Feature {i}']['Median Centroid'])/2

            lower_dt_index = np.argmin(np.abs((new_drifttimes - area_l_dt)))
            higher_dt_index = np.argmin(np.abs((new_drifttimes - area_h_dt))) 
            
            append_abund = smoothed_abund[lower_dt_index:higher_dt_index]
            append_drift = new_drifttimes[lower_dt_index:higher_dt_index]
            
            #reverse them for appending 
            fit_abund = fit_abund[::-1]
            fit_drift = fit_drift[::-1]
            
            for a in append_abund[::-1]:
                fit_abund.append(a)
            for d in append_drift[::-1]:
                fit_drift.append(d)
            
            #and puts them right for fitting
            fit_abund = fit_abund[::-1]
            fit_drift = fit_drift[::-1]
            
            
            # fit_abund = deque(list(fit_abund))
            # fit_drift = deque(list(fit_drift))
            # fit_abund.appendleft(append_abund)
            # fit_drift.appendleft(append_drift)
            # fit_abund = np.array(fit_abund[0]).T
            # fit_drift = np.array(fit_drift[0]).T
            
            
            #fit the smoothed abundances to a gaussian at roughly the feature location
            popt, pcov = curve_fit(plot_gauss_func, fit_drift, fit_abund, bounds=([amplitude*amp_low,centroid - lower_dt,std*std_low], [amplitude*amp_high,centroid + higher_dt,std*std_high]), maxfev=25000)
            plot_fit_abund = fit_abund
            plot_fit_drift = fit_drift
            fit_amp, fit_centroid, fit_sigma = popt
            gaussian = plot_gauss_func(new_drifttimes, fit_amp, fit_centroid, fit_sigma)
            gauss_sub = plot_gauss_func(fit_drift, fit_amp, fit_centroid, fit_sigma)
            fit_abund = np.array(fit_abund) - gauss_sub
            fit_abund = np.where(fit_abund <= 0, 1, fit_abund)   # setting all abundances to be positive
            fit_abund = list(fit_abund)
            
            
            # leftover_abund = np.array(fit_abund)
            leftover_abund = smoothed_abund - gaussian
            leftover_abund = np.where(leftover_abund <= 0, 1, leftover_abund)   # setting all abundances to be positive
            centroid_index = np.argmin(np.abs((fit_drift - (fit_centroid + 0*fit_sigma))))
            right_leftover_index = np.argmin(np.abs((new_drifttimes - (fit_centroid + 3*fit_sigma))))
            left_leftover_index = np.argmin(np.abs((fit_drift - (fit_centroid - 1*fit_sigma))))
            # if i == num_feats - 1:
                # right_leftover_index = left_leftover_index
            population_a = fit_abund[left_leftover_index:centroid_index]
            population_b = leftover_abund[right_leftover_index:]
            # break
            gauss_auc = integrate.trapezoid(gaussian, gauss_drift_times)      #integrating under the curve for the abundance data
            popa_auc = integrate.trapezoid(population_a, fit_drift[left_leftover_index:centroid_index])
            popb_auc = integrate.trapezoid(population_b, new_drifttimes[right_leftover_index:])
            
            prec_auc = gauss_auc + popa_auc
            
            
            total_auc = gauss_auc + popa_auc + popb_auc
            # fracr = (gauss_auc + popa_auc) / total_auc
            
            gaussian_data[key1][f'Feature {i+1}']['GAUSSIAN_AUC'].append(gauss_auc)
            gaussian_data[key1][f'Feature {i+1}']['LEFTOVER_AUC'].append(popa_auc)
            gaussian_data[key1][f'Feature {i+1}']['RIGHT_AUC'].append(popb_auc)
            gaussian_data[key1][f'Feature {i+1}']['AUC'].append(prec_auc)
            # gaussian_data[key1][f'Feature {i+1}']['BDC'].append(fracr)
            
            
            
            ax.fill_between(gauss_drift_times, gaussian, 0, color=colors[i])
            ax.fill_between(fit_drift[left_leftover_index:centroid_index], population_a, 0, color='r')
            ax.fill_between(new_drifttimes[right_leftover_index:], population_b, 0, color='k')
            ax.plot(new_drifttimes, smoothed_abund, color='k', ls='-')
            ax.plot(plot_fit_drift, plot_fit_abund, color='white', ls=':')
            ax.set_xlim(plot_low,plot_high)
            ax.tick_params(axis='both', which='major', labelsize=5)
            # break
        fig1.suptitle(f'Vin = {interp_vin}')
        fig1.tight_layout()
        pp.savefig(fig1)
        plt.close(fig1)
        # break
    fig,ax = plt.subplots()
    
    total_prec_auc = 0
    for i in range(num_feats):
        
        total_prec_auc = total_prec_auc + np.array(gaussian_data[key1][f'Feature {i+1}']['AUC'])  # Get the Area under the curve for each feature, (explain our successive feature assumption)
        bdc = total_prec_auc / (total_prec_auc + gaussian_data[key1][f'Feature {i+1}']['RIGHT_AUC'])
        
        
        
        # bdc = gaussian_data[key1][f'Feature {i+1}']['BDC'] #grab the bdc
        bdc = np.interp(exp_vins, new_Vins, bdc)        # and undo the Vin interpolation we did at the start so we are now only using the experimental Vins
        bdc_min_index = np.argmin(bdc)
        bdc_max_index = np.argmax(bdc)
        if bdc_min_index != len(bdc) - 1:
            bdc[bdc_min_index:] = bdc[bdc_min_index]
        if bdc_max_index != 0:
            bdc[:bdc_max_index] = bdc[bdc_max_index]
        gaussian_data[key1][f'Feature {i+1}']['BDC'] = bdc
        ax.plot(exp_vins, bdc, color=colors[i])
        ax.set_title('BDCs', fontsize=25)
        
        # These are just to help naming the breakdown curves so they are nice and neat titles
        if 'Myo' in key1:
            protein = 'Myo'
            feature = i + 1
            if 'myo7' in key1:
                charge = 7
            if '_8_' in key1:
                charge = 8
            if '_9_' in key1:
                charge = 9
            if 'myo10' in key1:
                charge = 10
            if 'AA' in key1:
                trip = 'A'
            if 'BA' in key1:
                trip = 'B'
            if 'CA' in key1:
                trip = 'C'
            if 'WM' not in key1:
                buffer = 'AA'
            if 'WM' in key1:
                buffer = 'WM'
            if 'high' in key1:
                pressure = 'high'
            if 'mid' in key1:
                pressure = 'mid'
            if 'low' in key1:
                pressure = 'low'
            
        
        
        if 'BSA' in key1:
            protein = 'BSA'
            buffer = ''
            if 'z17' in key1 or '_17_' in key1 or '3900' in key1:
                charge = 17
            if 'z16' in key1 or '_16_' in key1:
                charge = 16
            if 'z15' in key1 or '_15_' in key1:
                charge = 15
            if r'SF6' in key1 or 'N2' in key1 or '_seg' in key1:
                # charge = 17
                protein = '6560c_BSA'
            feature = i + 1
            
            if key1[-1] == 'A':
                trip = 'A'
            if key1[-1] == 'B':
                trip = 'B'
            if key1[-1] == 'C':
                trip = 'C'
            if key1[-1] == 'D':
                trip = 'D'
            
            
            
            if 'BSA15A' in key1:
                trip = 'A'
            if 'BSA15D' in key1:
                trip = 'B'
            if 'BSA15G' in key1:
                trip = 'C'
            
            if 'BSA16A' in key1:
                trip = 'A'
            if 'BSA16D' in key1:
                trip = 'B'
            if 'BSA16G' in key1:
                trip = 'C'
            
            
            if 'BSA17D' in key1:
                trip = 'A'
            if 'BSA17G' in key1:
                trip = 'B'
            if 'BSA17J' in key1:
                trip = 'C'
            
            if '_01_' in key1:
                trip = 'A'
            if '_02_' in key1:
                trip = 'B'
            if '_03_' in key1:
                trip = 'C'
            # if '_06_' in key1:
                # trip = 'D'
            
            # if '_A' in key1 in key1 or 'BSA17D' in key1: 
            # # if 'BSA17A' in key1:
            #     trip = 'A'
            # if '_B' in key1 or '_04' in key1 or 'BSA17G' in key1:
            # # if 'BSA17E' in key1:
            #     trip = 'B'
            # if '_C' in key1 or '_05' in key1 or 'BSA17J' in key1:
            # # if 'BSA17I' in key1:
            #     trip = 'C'
            # if '_3_mz' in key1 or '_06' in key1 or '_D' in key1:
            #     trip = 'D'
            # if '_4_mz' in key1:
            #     trip = 'B'
            # if '_5_mz' in key1:
            #     trip = 'C'
            # if '_6_mz' in key1:
            #     trip = 'D'
            # if '_7_mz' in key1:
            #     trip = 'E'
        
        if 'PK' in key1:
            protein = 'PK'
            charge = 33
            feature = i + 1
            buffer = ''
            if 'A.txt' in key1:
                trip = 'A'
            if 'B.txt' in key1:
                trip = 'B'
            if 'C.txt' in key1:
                trip = 'C'
        
        if 'CONA' in key1:
            protein = "ConA"
            feature = i+1
            if '20' in key1:
                charge = 20
            if '12' in key1:
                charge = 12
            if '11' in key1:
                charge = 11
            
            if r'a.txt' in key1 or r'A.txt' in key1:
                trip = 'A'
            elif r'I.txt' in key1 or r'i.txt' in key1:
                trip = 'B'
            elif r'c.txt' in key1 or r'C.txt' in key1:
                trip = 'C'

            else:
                trip = 'A'
                
            bound = ''
            
            
            
        if 'NISTmAb' in key1:
            protein = 'NistmAb'
            feature = i + 1
            buffer = ''
            if r'#24' in key1:
                charge = 24
            if r'_25_' in key1:
                charge = 25
            if r'#26' in key1:
                charge = 26
            if r'#27' in key1:
                charge = 27
            if r'#28' in key1:
                charge = 28 
            
            
            if r'a.txt' in key1 or r'A.txt' in key1:
                trip = 'A'
            if r'b.txt' in key1 or r'B.txt' in key1:
                trip = 'B'
            if r'c.txt' in key1 or r'C.txt' in key1:
                trip = 'C'
                
        # print(trip)
        outname = f'{out_path}/{protein}_{charge}_CIU_{feature}_{trip}.csv'
        
        output_array = np.array([exp_vins, bdc]).T
        
        np.savetxt(outname, output_array, delimiter=',')  #saves the breakdown curves in a folder called 'Breakdown_Curves' 
    
    pp.savefig(fig)
    plt.close(fig)
    
     
    # Creates the BDCs on fingerprint plots
    
    new_abundances[np.isnan(new_abundances)] = 0
    
    levels = get_contour_levels(stripped_abund/stripped_abund.max(axis=0))
    
    
    ciufig, ax1 = plt.subplots()
    ax1.contourf(exp_vins, exp_drift_times, plot_new_abundances/plot_new_abundances.max(axis=0), levels=levels, cmap='jet')
    ax2 = ax1.twinx()
    # ax1.set_ylim(0,12)
    # ax2.set_xlim(5,200)
    # ax2.set_ylim(0,1)
    
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.tick_params(axis='y', which='major', pad=10)
    ax1.tick_params(axis='x', which='major', pad=15)
    # ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    # ax1.minorticks_on()
    ax1.tick_params(axis='y', which='minor', length=15, width=2)
    ax1.tick_params(axis='x', which='minor', length=0, width=0)
    
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax2.tick_params(axis='y', which='major', pad=10)
    # ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    # ax1.minorticks_on()
    ax2.tick_params(axis='y', which='minor', length=15, width=2)
    ax2.tick_params(axis='x', which='minor', length=0, width=0)
    
    ax1.set_xlabel(r'Collision Voltage (V)',fontsize=10)
    ax1.set_ylabel(r'Drift time (ms)',fontsize=10)
    # ax2.set_ylabel(r'Fraction remaining',fontsize=55)
    # ax1.set_title(r'CIU and CID for BSA 17+ ion', fontsize=10)
    # plt.savefig('myo8_ciu_cid.pdf')
    # plt.figure()
    
    for i in range(num_feats):
        # 
        # Vins = gaussian_data[key1]['Vins']
        bdc = gaussian_data[key1][f'Feature {i+1}']['BDC']
        ax2.plot(exp_vins, bdc, marker='o',markersize=10, linewidth=3, color='w')
        
        ax1.set_ylim(plot_low, plot_high)
        ax2.set_ylim(0,1)
        
    pp.savefig(ciufig)
    pp.close()
    # break























































































