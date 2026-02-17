####################################
#
# spa_db_plot.py
#
# Original work by Ken Newton, et.at. at University of Oregon
# Copyright (c) 2024 James S. Prell, all rights reserved
#
# script to plot spa_heat data saved in the database format
#
# The ion, cell, and experiment data files can be specified via a JSON
# formatted input file.
#
# Will use the same json input file from the spa_mpirun program, adding a command line option
# to choose which plots to do.
#
####################################


import numpy as np
import pandas as pd
from scipy.optimize import minimize
from mpi4py import MPI
import json, sys
from os import path
from ionspa import const, ionclass, makecell, loadjfile
from time import sleep
import spa_db_fit as spafit
import sqlite3
import argparse
from spa_avg import main as avg_main
import matplotlib.pyplot as plt
plt.interactive(False)
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
'''
This program runs without the mpi system so will normally be called as:
python spa_db_plot.py --db {database_name} {input_file_name} --tpoints 1000 --xaxis {xaxis} --yaxis {yaxis} --vin {Vin}

DB table creation:
CREATE TABLE IF NOT EXISTS "cells" ("cid" integer primary key,"type" text,"L" float,"T" float,"rho" float,"Mgamu" float,"z0" float,"otherjson" text);
CREATE TABLE IF NOT EXISTS "ions" ("iid" integer primary key, "name" text, "mass" float, "charge" float, "CCS" float, "num_atoms" integer, "T0" float, "hcprofname" text);
CREATE TABLE IF NOT EXISTS "runs" ("rid" integer primary key, "cid" integer, "iid" integer, "Vin" integer, "tmax" float);
#CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float);
CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float, "collnum" float, "frac" float, "dfrac" float, "bfrag" float, "ma" float);
'''

def find_cid(dbname, celld):
    '''Search the database cells table to return the cid with all params matching the celld dict.'''
    cid = None

    od = celld.copy()   # copy the celld dict and delete normal required fields to be left with 'other' fields
    del(od['L'])
    del(od['type'])
    del(od['T'])
    try:
        del(od['rho'])
    except:
        del(od['pressure1'])
    del(od['Mgamu'])
    del(od['z0'])
    otherjson = json.dumps(od, separators=(',',':'),sort_keys=True)         # compact form with sorted keys
    # do the search
    tnames = 'cid,type,L,T,pressure1,rho,Mgamu,z0,otherjson'
    query = f"SELECT {tnames} from cells where type=\'{celld['type']}\' and T={celld['T']};"
    with sqlite3.connect(dbname) as db:
        cur = db.execute(query)
        fdat = cur.fetchall()
        for itemlist in fdat:
            itemd = dict(zip(tnames.split(','),itemlist))
            # each itemd matches cell type and Temperature, but also needs to match other fields
            # so check for match for L, rho, Mgamu, z0, otherjson
            try:
                if itemd['L'] == celld['L'] and itemd['rho'] == celld['rho'] and itemd['Mgamu'] == celld['Mgamu'] and itemd['z0'] == celld['z0'] and itemd['otherjson'] == otherjson:
                    print(f'found {itemd}')
                    cid = itemd['cid']
            except:
                if itemd['L'] == celld['L'] and itemd['pressure1'] == celld['pressure1'] and itemd['Mgamu'] == celld['Mgamu'] and itemd['z0'] == celld['z0'] and itemd['otherjson'] == otherjson:
                    print(f'found {itemd}')
                    cid = itemd['cid']
    return cid

def find_iid(dbname, iond):
    '''Search the database ions table to return the iid with all params matching the iond dict.'''
    iid = None

    # do the search
    tnames = 'iid,name,mass,charge,CCS,num_atoms,T0,hcprofname,dH,dS'
    query = f"SELECT {tnames} from ions where name=\'{iond['name']}\' and hcprofname=\'{iond['hcprofname']}\';"
    with sqlite3.connect(dbname) as db:
        cur = db.execute(query)
        fdat = cur.fetchall()
        for itemlist in fdat:
            itemd = dict(zip(tnames.split(','),itemlist))
            # each itemd matches ion name and hcprofname, but also needs to match other fields
            # so check for match for mass, charge, CCS, num_atoms, T0
            if itemd['mass'] == iond['mass'] and itemd['charge'] == iond['charge'] and itemd['CCS'] == iond['CCS'] and itemd['num_atoms'] == iond['num_atoms'] and itemd['T0'] == iond['T0']:
                print(f'found {itemd}')
                iid = itemd['iid']
    return iid

def find_rid(dbname, cid, iid, Vin):
    '''Search the database runs table to return the rid associated with cid and iid and Vin.
    returns rid,tmax for the (last) matching run.
    '''
    rid = None
    tmax = 0.0
    tnames = 'rid,cid,iid,Vin,tmax'
    query = f'SELECT {tnames} from runs where cid={cid} and iid={iid} and Vin={Vin}'
#    print(query)
    with sqlite3.connect(dbname) as db:
        cur = db.execute(query)
        fdat = cur.fetchall()
        for itemlist in fdat:
            itemd = dict(zip(tnames.split(','),itemlist))
            rid = itemd['rid']
            tmax = itemd['tmax']
    return rid, tmax

# CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float);

def get_rundat(dbname, rid, colname='Temp', xname='t', tpoints=None):
    if xname == 'z':
        query = f'SELECT irid, z, {colname} FROM trajs where rid={rid}'
    else:
        query = f'SELECT irid, {xname}, {colname}, z FROM trajs where rid={rid}'

    with sqlite3.connect(dbname) as db:
        df = pd.read_sql_query(query, db)

    dd = {}
    minirid = df.irid.min()
    maxirid = df.irid.max()
    for irid in range(minirid,maxirid+1):
        ta = df[xname].loc[df['irid'] == irid]
        va = df[colname].loc[df['irid'] == irid]
        zs = df['z'].loc[df['irid'] == irid]
        ta = np.array(ta)
        va = np.array(va)
        zs = np.array(zs)
        maxzarg = zs.argmax()
        if tpoints:
            tta = np.linspace(ta[0], ta[maxzarg], tpoints)
            vva = np.interp(tta, ta, va)
            dd[irid] = np.array((tta,vva))
        if not tpoints:
            dd[irid] = np.array((ta,va))
    return dd


def average_data(dd, colname='Temp', plotfile='screen', logvalues=None):
    k0 = list(dd.keys())[0]
    nions = len(dd.keys())
    d = dd[k0]
    t = d[0]    # might actually be xname (z or collnum)
    v = np.zeros_like(d[1])
    ionN = 0
    times = []
    for d in dd.values():
        ionN += 1
        v += d[1]
        if ionN < 301:
            plt.plot(d[0], d[1], lw=0.5)
            times.append(d[0][-1])    
            print(times)
    v = v / nions
    print(f'# of data points: {len(t)}')
    end_time = np.mean(times)
    ave_t = np.linspace(0, end_time, len(v))
    plt.plot(ave_t, v, 'k-', lw=1.5, label='Monte Carlo Avg')

    if logvalues:
        for key, logdat in logvalues.items():
            avg_ts = logdat.t
            avg_Tis = logdat.Ti
            plt.plot(avg_ts, avg_Tis, lw=1.5, label=key)

    title = 'data plot for all runs with selected Vin'
    plt.title(title)
    if xname == 't':
        plt.xlabel('t [s]')
    else:
        plt.xlabel(xname)
    plt.ylabel(colname)
    if plotfile == 'screen':
        plt.legend()
        plt.show()
    else:
        plt.savefig(plotfile)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---------------------- Common code -----------------------------

gdbname = 'spa_db.sqlite3'

# celld and iond not used, read from exp_input.json filename supplied on command line
celld = json.loads('''{
        "type": "aCIUcell",
        "L": 0.1,
        "Mgamu": 28,
        "T": 298.0,
        "rho": 2.42e-5,
        "z0": 0.0
    }''')
iond = json.loads('''{
       "name": "CC1222",
        "mass": 1.222,
        "charge": 1,
        "CCS": 2.8125,
        "num_atoms": 91,
        "T0": 298,
        "hcprofname": "tunemix1222"
    }''')
Vin = 425.0     # override this with value from command line


def main(Vin, colname='Temp', xname='t', plotfile='screen', tpoints=None, logvalues=None):

    cid = spafit.find_cid(gdbname, celld)
    iid = spafit.find_iid(gdbname, iond)
    rid,tmax = find_rid(gdbname, cid, iid, Vin)
    if rid is None:
        print(f'db records not found for Vin {Vin}')
        exit(1)
    
    print(f'cid: {cid}    iid: {iid}    rid: {rid}   tmax: {tmax}')
    
    dat = get_rundat(gdbname, rid, colname, xname, tpoints)
    average_data(dat, colname, plotfile, logvalues)
    return dat


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='python spa_db_plot.py',
        description='''Given a .json file (as from spa_mpirun), plot
        the variable vname vs time for the givin vin voltage.
        the --db argument can specify an alternate database file name.''')
    parser.add_argument('jsonfile', nargs='+')
    parser.add_argument('--yaxis', choices=['Temp', 'u', 'z', 'KE', 'collnum', 'frac', 'dfrac', 'bfrag', 'ma'], default='Temp') #   , nargs=1)
    parser.add_argument('--vin', default=10) #   , nargs=1)
    parser.add_argument('--db', default = 'spa_db.sqlite3')
    parser.add_argument('--xaxis', choices=['t', 'z', 'collnum'], default='t')
    parser.add_argument('--plotfile', default='screen', help='if specified, output plot to this filename.')
    parser.add_argument('--tpoints', default='none', help='number of time points to interpolate, default is no interpolation')
    args = parser.parse_args()
    gdbname = path.splitext(args.db)[0] + '.sqlite3'
#    print(args)
    paramd = loadjfile(args.jsonfile)
    paramd.update(colname=args.yaxis, Vin=args.vin)
    xname = args.xaxis
    plotfile = args.plotfile
#    for jname in args.jsonfile:
#        with open(jname, 'r') as f:
#            paramd.update(json.load(f))
#    cellfname = paramd.get("cellfilename")
#    if cellfname:
#        with open(cellfname, 'r') as fc:
#            fcdict = json.load(fc)
#            paramd.update(fcdict)
#    ionfname = paramd.get("ionfilename")
#    if ionfname:
#        with open(ionfname, 'r') as fi:
#            fidict = json.load(fi)
#            paramd.update(fidict)

    print(paramd, '\n')

    colname = paramd.get('colname')
    sVin = paramd.get('Vin')
    tpointsarg = args.tpoints
    if tpointsarg == 'none':
        tpoints = None
    else:
        tpoints = int(tpointsarg)

    celld = paramd['cell']      # don't need celld.Voffset here since plot is for specific Vin from command line
    iond = paramd['ion']
    newtemp = celld.get('T0offset', 0.0) + iond.get('T0')
    iond.update({'T0' : newtemp})
    Vin = float(sVin)
    logvalues = {}
#    logvalues['IICT'], extra = avg_main(cell=celld, ion=iond, injectV=Vin-5.5, option='IICT', plotfile='none')
#    logvalues['ELASTIC'], extra = avg_main(cell=celld, ion=iond, injectV=Vin-5.5, option='ELASTIC', plotfile='none')
#    logvalues['ELASTICHO'], extra = avg_main(cell=celld, ion=iond, injectV=Vin-5.5, option='ELASTICHO', plotfile='none')
#    logvalues['INELTH'], extra = avg_main(cell=celld, ion=iond, injectV=Vin-5.5, option='INELTH', plotfile='none')

#    print(f'celld: {celld}\niond: {iond}\ncolname: {colname}\nVin: {Vin}')
    dat = main(Vin, colname, xname, plotfile, tpoints, logvalues)
