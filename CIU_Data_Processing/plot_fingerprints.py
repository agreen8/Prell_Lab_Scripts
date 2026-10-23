# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 13:50:01 2026

@author: Austin
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# from scipy.optimize import minimize, curve_fit
import glob
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

def get_contour_levels(ciu_data, merge_cutoff=5, num_contours=1000):
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


bsa_15_6560c_1V_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260220 6560C 1V BSA\15+\BSA_CIU_1V_01_z15_seg_raw.csv"
bsa_16_6560c_1V_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260220 6560C 1V BSA\16+\BSA_CIU_1V_02_z16_seg_raw.csv"
bsa_17_6560c_1V_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260220 6560C 1V BSA\17+\BSA_CIU_1V_01_z17_seg_raw.csv"

bsa_15_6560c_ave_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\250323 6560C 10V BSA\15+\20250323_BSA_CIU_N2_15_A_Avg_raw.csv"
bsa_15_6560c_10V_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\250323 6560C 10V BSA\15+\20250323_BSA_CIU_N2_15_A_raw.csv"

bsa_16_6560c_ave_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\250323 6560C 10V BSA\16+\20250323_BSA_CIU_N2_16_A_Avg_raw.csv"
bsa_16_6560c_10V_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\250323 6560C 10V BSA\16+\20250323_BSA_CIU_N2_16_A_raw.csv"

bsa_17_6560c_ave_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\250323 6560C 10V BSA\17+\20250323_BSA_CIU_N2_03_mz3900-3950_seg_Avg_raw.csv"
bsa_17_6560c_10V_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\250323 6560C 10V BSA\17+\20250323_BSA_CIU_N2_03_mz3900-3950_seg_raw.csv"

bsa_15_synapt_old_ave_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260216 Synapt BSA\15+\DT_260216_BSA_15_3V_AA.raw_#BSA15A-_Avg_raw.csv"
bsa_15_synapt_old_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260216 Synapt BSA\15+\DT_260216_BSA_15_3V_AA.raw_#BSA15A-ccs_raw.csv"

bsa_16_synapt_old_ave_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260216 Synapt BSA\16+\DT_260218_BSA_16_3V_AA.raw_#BSA16A-_Avg_raw.csv"
bsa_16_synapt_old_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260216 Synapt BSA\16+\DT_260218_BSA_16_3V_AA.raw_#BSA16A-ccs_raw.csv"

bsa_17_synapt_old_ave_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260216 Synapt BSA\17+\DT_260216_BSA_17_3V_DA.raw_#BSA17D-_Avg_raw.csv"
bsa_17_synapt_old_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260216 Synapt BSA\17+\DT_260216_BSA_17_3V_DA.raw_#BSA17D-ccs_raw.csv"

bsa_15_synapt_n_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260525 Synapt BSA\15+\DT_260525_BSA_POS_3V_AA.raw_#BSA15A-ccs_raw.csv"
bsa_16_synapt_n_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260525 Synapt BSA\16+\DT_260525_BSA_POS_3V_AA.raw_#BSA16A-_Avg_raw.csv"
bsa_17_synapt_n_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\260525 Synapt BSA\17+\DT_260525_BSA_POS_3V_AA.raw_#BSA17A-ccs_raw.csv"


bsa_15_synapt_new_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\261008 Synapt new BSA\15+\DT_261008_BSA_New_3V_AA.raw_#BSA15A-ccs_raw.csv"
bsa_16_synapt_new_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\261008 Synapt new BSA\16+\DT_261008_BSA_New_3V_AA.raw_#BSA16A-ccs_raw.csv"
bsa_17_synapt_new_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\261008 Synapt new BSA\17+\DT_261008_BSA_New_3V_AA.raw_#BSA17A-ccs_raw.csv"

bsa_15_synapt_new_new_csv = r"C:\Users\Austin\Documents\Grad School\Austin Green BSA CIU Paper OFFLINE\CIU_BDC_generation\261013 Synapt new BSA\15+\DT_261013_BSA_New_3V_AA.raw_#BSA15A-ccs_raw.csv"

abund_data = np.loadtxt(bsa_17_synapt_old_ave_csv, dtype='str', delimiter=',')    
abund_data[0][0] = 0.0
abund_data = abund_data.astype('float')
abund_data = abund_data[::1, ::1]


stripped_abund = abund_data[1:,1:]
exp_vins = abund_data[0][1:]
exp_drift_times = abund_data[1:,0]


dx = exp_drift_times[1] - exp_drift_times[0]

plot_new_abundances = sgolay2d(stripped_abund, 5, 2)
plot_new_abundances = stripped_abund
levels = get_contour_levels(plot_new_abundances/plot_new_abundances.max(axis=0))

# test = plot_new_abundances/max(integrate.trapezoid(plot_new_abundances, axis=0, dx=dx))
# test2 = test/np.max(test)
# levels = get_contour_levels(test2)

ciufig, ax1 = plt.subplots()
ax1.contourf(exp_vins, exp_drift_times, plot_new_abundances/plot_new_abundances.max(axis=0), levels=levels, cmap='jet')
# ax1.contourf(exp_vins, exp_drift_times, test2, levels=levels, cmap='jet')
# ax2 = ax1.twinx()


ax1.set_ylim(4000, 7500)
ax1.tick_params(axis='both', which='major', labelsize=30, length=15)
























































































