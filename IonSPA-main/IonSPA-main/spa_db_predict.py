####################################
#
# spa_predict.py
#
# Original work by Ken Newton, et.at. at University of Oregon
# Copyright (c) 2024 James S. Prell, all rights reserved
#
# non-MPI script to run ionspa for given ion with given cell over a range of
# input energies and, with given deltaH and deltaS, predict the fragmentation curve.
# (derived from exp_fit.py)
#
# Use spa_heat_impulsive2 model for ion heating, running a collection of ions
# Ion time-temperature arrays are interpolated to fixed time points, averaged
# and used with the fracremains function and deltaH, deltaS pair to pridict
# an experimental fragmentation curve.
#
# The ion, cell, and experiment data files can be specified via a JSON
# formatted input file.
#
# Remove all MPI code from exp_predict and read saved data from the DB file
#
####################################


import numpy as np
import pandas as pd
import json
from os import path
from ionspa import loadjfile
import spa_db_fit as spafit
import sqlite3
import matplotlib.pyplot as plt
import argparse

plt.interactive(False)
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
'''
This program run without the mpi system so will normally be called as:
python spa_db_predict.py --db {database_num} {input_file_name} --dH {dH_value} --dS {dS_value}

DB table creation:
CREATE TABLE IF NOT EXISTS "cells" ("cid" integer primary key,"type" text,"L" float,"T" float,"rho" float,"Mgamu" float,"z0" float,"otherjson" text);
CREATE TABLE IF NOT EXISTS "ions" ("iid" integer primary key, "name" text, "mass" float, "charge" float, "CCS" float, "num_atoms" integer, "T0" float, "hcprofname" text);
CREATE TABLE IF NOT EXISTS "runs" ("rid" integer primary key, "cid" integer, "iid" integer, "Vin" integer, "tmax" float);
#CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float);
CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float, "collnum" float, "frac" float, "dfrac" float, "bfrag" float, "ma" float);

'''



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---------------------- Common code -----------------------------

gdbname = 'spa_db.sqlite3'


def main(celld, iond, dH, dS, Vins, fracs, title=''):
    '''Show predicted fraction remaining vs. Vin for given inputs.
    Plot results with exp data if present.'''

    print(f'pspa_db_predict:\n  celld: {celld}\n  iond: {iond}\n  dH: {dH}  dS: {dS}\n  Vins: {Vins}\n  fracs: {fracs}')

    # Now just get the time/temperature data for each Vin, and calculate the frac_remaining vs Vin.
    vtTd = {}   # {5.0: [[times...], [temps...]], 30.0: [[times...], [temps...]], ...}
    cid = spafit.find_cid(gdbname, celld)
    iid = spafit.find_iid(gdbname, iond)

    # read data from DB and build the tTd structure
    tTd = {}
    for indx, Vin in enumerate(Vins):
        rid, tmax = spafit.find_rid(gdbname, cid, iid, Vin)
        dd = spafit.get_rundat(gdbname, rid, 'Temp')
        times, Temps = np.array(spafit.average_data(dd, 'Temp'))
        tTd[Vin] = [times, Temps, tmax]

    rs = []
    sumsq = 0.0
    print(f'Calculating remaining fraction with dH: {dH:0.2f}, dS: {dS:0.2f}')
    fitter = spafit.fitclass(Vins, tTd, fracs)
    sumsq = fitter.diffsq([dH,dS], printvals=True, savers=True)
    print(f'diffsq returned sumsq: {sumsq:0.2e}')
    rs = fitter.rs

#    for indx, Vin in enumerate(Vins):
#        times, Temps, tmax = tTd[Vin]
#        r = spah.fracremains(times, Temps, dH*1000.0, dS, tmax)
#        rs.append(r)
#        print(f'   Vin: {Vin:0.1f}  tmax: {tmax*1000:0.2f} ms   exp: {fracs[indx]:0.2f} --> calcr: {r:0.2f}')
#        sumsq += (fracs[indx] - r)**2

    plt.plot(Vins, fracs, 'x')
    plt.plot(Vins, rs)
    plt.xlabel('Vin [V]')
    plt.ylabel('frac remaining')
    plt.title(title + f'  R2: {sumsq:0.02e}')
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='python spa_db_predict.py',
        description='''Given a .json file (as from spa_mpirun), and either the dhds_start from
        that file or specified dh and ds, plot the fraction remaining vs Vin.''')
    parser.add_argument('jsonfile', nargs='+')
    parser.add_argument('--dH') #   , nargs=1)
    parser.add_argument('--dS') #   , nargs=1)
    parser.add_argument('--db', default = 'spa_db.sqlite3')
    args = parser.parse_args()
    gdbname = path.splitext(args.db)[0] + '.sqlite3'

#    print(args)
    jfile = args.jsonfile
    paramd = loadjfile(jfile)
    # paramd = {}
    # for jfile in args.jsonfile:
    #     with open(jfile, 'r') as f:
    #         paramd.update(json.load(f))
    # cellfname = paramd.get("cellfilename")
    # if cellfname:
    #     with open(cellfname, 'r') as fc:
    #         fcdict = json.load(fc)
    #         paramd.update(fcdict)
    # ionfname = paramd.get("ionfilename")
    # if ionfname:
    #     with open(ionfname, 'r') as fi:
    #         fidict = json.load(fi)
    #         paramd.update(fidict)
    
    expd = paramd['exp']
    celld = paramd['cell']
    cellVoffset = celld.get('Voffset', 0.0)     # read the offset voltage to be added to Vins, default is 0.0
    iond = paramd['ion']
    newtemp = celld.get('T0offset', 0.0) + iond.get('T0')
    iond.update({'T0' : newtemp})
    dH = float( args.dH or paramd['dhds_start'][0] )
    dS = float( args.dS or paramd['dhds_start'][1] )

#    print(jfile, dH, dS)

    Vins = expd.get('Vins', None)
    fracs = None
    if not Vins:
        expdatpath = path.join(expd['folder'], expd['file'])
        gexp = pd.read_table(expdatpath, delimiter=',', header=None, names=['Vin', 'efrac'])
        #print(gexp)
        Vins = list(gexp.Vin.values)
        fracs = list(gexp.efrac.values)

    Vins = [vin+cellVoffset for vin in Vins]
    
    title = f'{jfile} dH={dH}, dS={dS}'
    main(celld, iond, dH, dS, Vins, fracs, title)
