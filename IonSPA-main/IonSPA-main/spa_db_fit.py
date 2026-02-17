####################################
#
# spa_fit.py
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
# The dH and dS values will be adjusted to get a best fit of the fragmentation curve.
#
# Fits are done using several steps to speed the fit times.  Assumptions made about fit function shape tested for this process.
#
####################################


import numpy as np
import pandas as pd
from scipy.optimize import minimize, curve_fit, fsolve
import json
from os import path
from ionspa import loadjfile, fracloss, lossrate, fracremains
from time import time
import random
import math
import sqlite3
import matplotlib.pyplot as plt
import argparse
from random import random
import matplotlib
from scipy import stats,integrate
import ionspa
import warnings
warnings.filterwarnings("error")
const = ionspa.const
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

plt.interactive(False)

'''
This program runs without the mpi system so will normally be called as:
python spa_db_predict.py --db {database_name} {input_file_name} --{options}

DB table creation done in spa_mpirun.py:
CREATE TABLE IF NOT EXISTS "cells" ("cid" integer primary key,"type" text,"L" float,"T" float,"rho" float,"Mgamu" float,"z0" float,"otherjson" text);
CREATE TABLE IF NOT EXISTS "ions" ("iid" integer primary key, "name" text, "mass" float, "charge" float, "CCS" float, "num_atoms" integer, "T0" float, "hcprofname" text);
CREATE TABLE IF NOT EXISTS "runs" ("rid" integer primary key, "cid" integer, "iid" integer, "Vin" integer, "tmax" float);
#CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float);
#CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float, "collnum" float, "frac" float, "dfrac" float, "bfrag" float
CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float, "collnum" float, "frac" float, "dfrac" float, "bfrag" float, "ma" float);

'''

### DB find functions mostly copied from spa_db_plot.py

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
    dbcellkeys = tnames.split(',')
    with sqlite3.connect(dbname) as db:
        cur = db.execute(query)
        fdat = cur.fetchall()
        for itemlist in fdat:
            itemd = dict(zip(dbcellkeys,itemlist))
            # each itemd matches cell type and Temperature, but also needs to match other fields
            # so check for match for L, rho, Mgamu, z0, otherjson
            cmatch = True
            for ckey in dbcellkeys[1:-1]:
                fitem = itemd.get(ckey)
                fcitem = celld.get(ckey)
                if fitem != fcitem:
                    cmatch = False
#                    print(f'not matching for {ckey}:')
#                    print(f'  {fitem}')
#                    print(f'  {fcitem}\n')
            fitem = itemd.get('otherjson')
            if cmatch and fitem != otherjson:
                cmatch = False
#                print(f'not matching for otherjson')
#                print(f'  {fitem}')
#                print(f'  {otherjson}\n')            
            if cmatch:
                print(f'found {itemd}\n')
                cid = itemd['cid']
            else:
                ecid = itemd['cid']
#                print(f'unmatched: for cid {ecid}\n  {itemd}\n {celld}\n')
#            if itemd['L'] == celld['L'] and itemd['rho'] == celld['rho'] and itemd['Mgamu'] == celld['Mgamu'] and itemd['z0'] == celld['z0'] and itemd['otherjson'] == otherjson:
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
            match_mass =  (abs(itemd['mass'] - iond['mass']) < 0.001)
            match_charge =  (abs(itemd['charge'] - iond['charge']) < 1)
            match_ccs = (abs(itemd['CCS'] - iond['CCS']) < 0.001)
            match_numatoms = (itemd['num_atoms'] == iond['num_atoms'])
            match_T0 =  (abs(itemd['T0'] - iond['T0']) < 0.05)
            ionfound = match_mass and match_charge and match_ccs and match_numatoms and match_T0
            if ionfound:
                print(f'found {itemd}')
                iid = itemd['iid']
                return iid
            else:
                print(f'unmatched:  {itemd}')
                print(f'match to: {iond}')
                print(f'matches: mass:{match_mass}, charge:{match_charge}, ccs:{match_ccs}, num_atoms:{match_numatoms}, T0:{match_T0}')
    return iid


def find_rid(dbname, cid, iid, Vin):
    '''Search the database runs table to return the rid associated with cid and iid and Vin.
    returns rid,tmax for the (last) matching run.
    '''
    rid = None
    tmax = 0.0
    tnames = 'rid,cid,iid,Vin,tmax'
    query = f'SELECT {tnames} from runs where cid={cid} and iid={iid} and Vin={Vin}'
    with sqlite3.connect(dbname) as db:
        cur = db.execute(query)
        fdat = cur.fetchall()
        for itemlist in fdat:
            itemd = dict(zip(tnames.split(','),itemlist))
            rid = itemd['rid']
            tmax = itemd['tmax']
    return rid, tmax


def get_fracdat(dbname, rid):
    '''Return a dataframe with data relevant to recomputing the frac, dfrac and bfrag columns. Used when updating the database, currently not working and unsupported'''
    with sqlite3.connect(dbname) as db:
        query = f'SELECT irid, t, z, Temp, frac, dfrac, bfrag FROM trajs where rid={rid}'
        df = pd.read_sql_query(query, db)
    return df


def get_rundat(dbname, rid, colname='Temp'):
    '''Return a dict extracted from the DB for the given runid with keys as
    the individual irid and values as a pair with time and the value from
    the column specified.'''
    with sqlite3.connect(dbname) as db:
        query = f'SELECT irid, t, {colname} FROM trajs where rid={rid}'
        df = pd.read_sql_query(query, db)

    dd = {}
    maxirid = df.irid.max()
    for irid in range(1,maxirid+1):
        ta = df.t.loc[df['irid'] == irid]
        va = df[colname].loc[df['irid'] == irid]
        dd[irid] = np.array((ta,va))
    return dd


def average_data(dd, colname='Temp'):
    '''dd is a dict with keys the run id for each ion and values that are
    2-tuples with time and a value specified earlier in the get_rundat function.
    With dd, compute the averate value for all ions at each time.
    Note: the colname param is not used at this time.
    The column averaged determined by the data passed to the function.
    The averaging is the average value of the selected parameter at each time.'''
    k0 = list(dd.keys())[0]
    nions = len(dd.keys())
    d = dd[k0]
    t = d[0]
    v = np.zeros_like(d[1])
    ionN = 0
    for d in dd.values():
        ionN += 1
        dsub = d[1][:len(v)]
#        print(v.shape, d[1].shape)
        v += dsub
    v = v / nions
    return t, v


# new functions to compute average fraction remaining of all individual trajectories for this run (rid) at the matching Vin
def get_tzT(dbname, rid, tpoints=0, Tscale=1.0):
    '''Return data needed for fragmentation calculation.
    Will need one of these data structs for each Vin (each rid)'''
    df = None
    print('Opening the database', flush=True)
    with sqlite3.connect(dbname) as db:
        query = f'SELECT irid, t, z, Temp FROM trajs WHERE rid={rid}'
        df = pd.read_sql_query(query, db)
    print('Got the data', flush=True)
    # print(df['Temp'])
    # print()
    # now reorganize for easy use
    # df['Temp'] = df['Temp']
    
    # print(df['Temp'])
    maxirid = df.irid.max()
    # print(df)
    iirids = []
    ttas = []
    vvas = []
    vva95s = []
    zzas = []
    if tpoints != 0:
        for irid in range(1,maxirid+1):
            ta = df['t'].loc[df['irid'] == irid]
            #new_ta = [0, ta[0]-0.0000001+ 0.001278]
            #for t in ta:
            #    new_ta.append(t + 0.001278)
            va = df['Temp'].loc[df['irid'] == irid]
            #new_va = [400, 400]
            #for v in va:
            #    new_va.append(v)
            va95 = df['Temp'].loc[df['irid'] == irid]*0.95
            #new_va95 = [400, 400]
            #for v in va95:
            #    new_va95.append(v)
            za = df['z'].loc[df['irid'] == irid]
            #new_za = [0,0.000000000001]
            #for z in za:
            #    new_za.append(z)
            #ta = np.array(new_ta)
            #va = np.array(new_va)
            #va95 = np.array(new_va95)
            #za = np.array(new_za)
            ta = np.array(ta)
            va = np.array(va)
            va95 = np.array(va95)
            za = np.array(za)
            maxzarg = za.argmax()
            tta = np.linspace(ta[0], ta[maxzarg], tpoints)
            # print(len(tta))
            vva = np.interp(tta, ta, va)
            vva95 = np.interp(tta, ta, va95)
            zza = np.interp(tta, ta, za)
            iirid = np.full(len(tta), irid)
            # print(len(iirid))
            # iirids.append(ii for ii in iirid)
            for ii, tt, vv, zz, vv95 in zip(iirid, tta, vva, zza, vva95): 
                iirids.append(ii)
                ttas.append(tt)
                vvas.append(vv)
                vva95s.append(vv95)
                zzas.append(zz)
        dfv = np.array([iirids, ttas, zzas, vvas]).T
        dfv95 = np.array([iirids, ttas, zzas, vva95s]).T
    else:
        df95 = df.copy()
        df95['Temp'] = df95['Temp']*0.95
        dfv = df.values
        dfv95 = df95.values
    # print(dfv)
    return dfv, dfv95


def computefracs(dfv, dH, dS):
    '''Given an np.array dfv with columns(irid, t, z, T) and a specified dH and dS, 
    Compute the average fraction remaining at the end of the ion trajectories.
    '''
    # dfv is a np.array with columns irid, t, z, Temp
    # each irid is a different ion trajectory
    # will use this by finding dt array for each irid, and use Temp vs. dt to get ionloss until point where ion reaches max(z)
    iridset = set(dfv[:,0])
    Ntrajs = len(iridset)
    lastfraca = np.zeros(Ntrajs)
    ionnum = 0
    for irid in iridset:
        t = dfv[:,1][dfv[:,0]==irid]
        dt = t[1:] - t[:-1]
        z = dfv[:,2][dfv[:,0]==irid][1:]
        T = dfv[:,3][dfv[:,0]==irid][1:]
        # maxz = max(z)
        maxzarg = z.argmax()
        dfrac = lossrate(T, dH, dS)            # full array of loss rates (close to 0)
        dfrac[maxzarg:] = 0.0                       # declare zero loss after ion reaches the end
        dfracremains = np.exp(-dfrac*dt)            # with small losses, usually close to 1
        xa = np.arange(0, len(dfracremains))
        #plt.plot(xa, dfracremains)
        #plt.show()
        lastfraca[ionnum] = dfracremains.prod()
        ionnum = ionnum + 1
    quartile = int(max(iridset)/4)                          #sub-selects the middle 50% frac remains
    lastfraca = np.sort(lastfraca)                          #seems to improve fit results
    fracremains = np.mean(lastfraca[quartile:-quartile])
    return fracremains

def get_CID50(Vins, fracs):
    '''Finds the Vin where the fraction remaining is 0.5 by linear interpolation'''
    fracs = np.array(fracs)
    try:
        ibelow = np.where(fracs <= 0.5)[0][0]
    except:
        ibelow = np.where(fracs == fracs[-1])[0][0]
    iabove = ibelow - 1
    slope = (fracs[ibelow] - fracs[iabove])/(Vins[ibelow] - Vins[iabove]) 
    intercept = fracs[iabove] - slope*Vins[iabove]
    CID50 = (0.5 - intercept)/slope
    
    cindex = (np.abs(fracs - 0.5)).argmin()
    
    return CID50, cindex

def get_lifetime(Vins, fracs):
    '''Finds the Vin where the fraction remaining is 1/e by linear interpolation'''
    fracs = np.array(fracs)
    try:
        ibelow = np.where(fracs <= 1/np.exp(1))[0][0]
    except:
        ibelow = np.where(fracs == fracs[-1])[0][0]
    iabove = ibelow - 1
    slope = (fracs[ibelow] - fracs[iabove])/(Vins[ibelow] - Vins[iabove]) 
    intercept = fracs[iabove] - slope*Vins[iabove]
    CIDlife = (1/np.exp(1) - intercept)/slope
    
    cindex = (np.abs(fracs - 1/np.exp(1))).argmin()
    
    return CIDlife, cindex

def compute_max_T(dfv):
    '''Computes the average max Temperature of the IonSPA trajectories'''
    iridset = set(dfv[:,0])
    # print(max(iridset))

    max_Ts = []
    for irid in iridset:
        t = dfv[:,1][dfv[:,0]==irid]
        dt = t[1:] - t[:-1]
        z = dfv[:,2][dfv[:,0]==irid][1:]
        T = dfv[:,3][dfv[:,0]==irid][1:]
        max_T = max(T)
        max_Ts.append(max_T)
    
    quartile = int(max(iridset)/4)                  #sub-selects the middle 50% temperatures to match data used in fit
    max_Ts = np.sort(max_Ts)                        
    max_T = np.mean(max_Ts[quartile:-quartile])
    return max_T

def compute_rate_integral(dfv, dH, dS):
    '''Given an np.array dfv with columns(irid, t, z, T) and a specified dH and dS, 
    Compute the average fraction remaining at the end of the ion trajectories.
    '''
    # dfv is a np.array with columns irid, t, z, Temp
    # each irid is a different ion trajectory
    # will use this by finding dt array for each irid, and use Temp vs. dt to get ionloss until point where ion reaches max(z)
    iridset = set(dfv[:,0])
    Ntrajs = len(iridset)
    int_results = np.zeros(Ntrajs)
    ionnum = 0
    max_Ts = []
    max_ts = []
    mean_Ts = []
    for irid in iridset:
        t = dfv[:,1][dfv[:,0]==irid]
        dt = t[1:] - t[:-1]
        z = dfv[:,2][dfv[:,0]==irid][1:]
        T = dfv[:,3][dfv[:,0]==irid][:]
        # print(t)
        # maxz = max(z)
        maxzarg = z.argmax()
        dfrac = lossrate(T, dH, dS)            # full array of loss rates (close to 0)
        try:
            bi = np.where(dfrac > 1)[0][0]
            ei = np.where(dfrac > 1)[0][-1]
            integral = integrate.trapz(dfrac[bi:ei], t[bi:ei])
            t_max = t[ei] - t[bi]
            # print('HERE', tmax)
        except:
            t_max = max(t) 
            integral = integrate.trapz(dfrac, t)
        max_ts.append(t_max)
        # print(integral)
        max_T = max(T)
        max_Ts.append(max_T)
        int_results[ionnum] = integral
        mean_T = np.mean(T)
        mean_Ts.append(mean_T)
        ionnum = ionnum + 1
    quartile = int(max(iridset)/4)                  #sub-selects the middle 50% temperatures to match data used in fit
    max_Ts = np.sort(max_Ts)           
    # print(dfv)
    # print(max_Ts)             
    max_ts = np.sort(max_ts)
    mean_Ts = np.sort(mean_Ts)
    max_T = np.mean(max_Ts[quartile:-quartile])
    max_t = np.mean(max_ts[quartile:-quartile])
    kint = np.mean(int_results[quartile:-quartile])
    # mean_T = np.mean(mean_Ts[quartile:-quartile])
    return kint, max_T, max_t

def effective_temp(T, *args):
    '''Given the average maximum flight time of the trajectories, dH, dS, 
    and the integral under the Temperature vs time curve, calculates the
    constant Temperature that will give the same kinetic results as the
    time varying one'''
    
    dH, dS, k_int, t_max = args
    dHscale = 1000.0*dH/const.R
    dSscale = dS/const.R
    
    func = const.kBoverh*T*np.exp(dSscale - dHscale/T) - k_int/t_max
    
    
    return func

def error_calc(kint, Teff, dG, teff, expfrac, fitfrac, weights, Rcid50, dS):
    
    relTerr = 0.05
    kb = 8.314/1000 #boltzmann constant in kj/mol*K
    
    kint = kint/teff
    
    # print(f'k_cid50: {kint}')
    # print(f'T_eff_cid50: {Teff}')
    # print(f'dG: {dG}')
    # print(f'teff: {teff}')
    # print(f'Rcid50: {Rcid50}')
    
    kbt_dG = kb*Teff/dG
    kbT = kb*Teff
    # print(f'kbT/dG: {kbt_dG}')
    # print(f'kbT: {kbT}')
    dR2 = sum((expfrac - fitfrac)**2 * weights**2)/(sum(weights) - 1)
    top_sum = sum((expfrac - fitfrac)**2 * weights**2)
    bot_sum = (sum(weights) - 1)
    # print(f'dR2: {dR2}')
    # print(f'top sum: {top_sum}')
    # print(f'bot sum: {bot_sum}')
    
    
    rel_dk_squared = (1/kint**2) * (1/teff**2) * 1/(Rcid50**2) * sum((expfrac - fitfrac)**2 * weights**2)/(sum(weights) - 1)
    dk = kint*np.sqrt(rel_dk_squared)
    # print(f'rel_dk_squared: {rel_dk_squared}')
    # print(f'dk: {dk}')
    
    rel_ddG_squared = (relTerr)**2 * (1 + 2*kb*Teff/dG + (kb*Teff/dG)**2) + rel_dk_squared*(kb*Teff/dG)**2
    T_err = np.sqrt((relTerr)**2 * (1 + 2*kb*Teff/dG + (kb*Teff/dG)**2))
    k_err = np.sqrt(rel_dk_squared*(kb*Teff/dG)**2)
    # print(f'rel_ddG_squared: {rel_ddG_squared}')
    # print(f'dG error from T: {T_err}')
    # print(f'dG error from k: {k_err}')
    
    ddG = np.sqrt(rel_ddG_squared) * dG
    
    dT = Teff * relTerr
    # print(f'dT: {dT}')
    ddG2 = ddG**2
    st = (dS*dT/1000)**2
    
    # print(f'ddG2: {ddG2}')
    # print(f'Second Term: {st}')
    # print(f'dS: {dS}')
    #ddH = np.sqrt(((ddG**2 - (dS*dT/1000)**2)/2))
    ddH = 0
    ddS = 1000 * ddH / Teff
    
    # print(f'ddG: {ddG}')
    # print(f'ddH: {ddH}')
    # print(f'ddS: {ddS}')
     
    outputs = (kint, Rcid50, kbt_dG, dR2, top_sum, bot_sum, rel_dk_squared, dk, rel_ddG_squared, T_err, k_err, ddG, ddH, ddS)
    
    return outputs
    


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---------------------- Common code -----------------------------

gdbname = 'spa_db.sqlite3'
usenewfrac = False

def doupdate(celld, iond, newdH, newdS):
    '''
    Recalculate the frac remaining for each ion and store back into the
    database.
    This code is probably not complete.
    '''
    updatecount = 0

    cid = find_cid(gdbname, celld)
    iid = find_iid(gdbname, iond)
    if not (cid and iid):
        print(f'cid {cid}    iid {iid}  error')
        exit(1)

    # get list of Vins
    query = f"SELECT rid, Vin from runs where cid={cid} and iid={iid};"
    fdat = None
    with sqlite3.connect(gdbname) as db:
        cur = db.execute(query)
        fdat = cur.fetchall()
        # now extract Vin list from fdat (associated with rid?)
    testcount = 0
    for rid, vin in fdat:      # this might work
        testcount += 1
    
        df = get_fracdat(gdbname, rid)
        # df has columns: irid, t, Temp, frac, dfrac, bfrag
        # we could use a dt column as t[1:] - t[:-1]
        # then call fracloss(deltat, T, deltaH, deltaS) for each row
        # insert appropriate results in frac, dfrac, bfrag

        # chatgpt suggested approach:
        '''
        df = df.sort_values(by=['rid', 'irid', 't'])
        dtfill = df.t.iloc[1] = df.t.iloc[0]
        df['dt'] = df.groupby(['rid', 'irid'])['t'].diff().fillna(dtfill)
        # with current code, dt is the same throughout a run, so fill 0 vales with the first value
        # so remove the fillna above, get the first value, and use fillna on that

        lastdfrac = 0
        lastdt = 0
        def calcdfrac(Temp, dt):
#            xdt = dt or lastdt     # don't know if resetting dt might be needed, possibly when new irid value is encountered
            thisloss = fracloss(dt, Temp, newdH, newdS)
            xdfrac = thisloss / dt
            return xdfrac
        
        df['dfracnew'] = df.apply(lambda x calcdfrac(x['Temp'], x['dt']), axis=1)
        df.drop('dt', axis=1)
        print(df)
        '''

        # complicated way to make a column with dt values
        # may not be needed since dt will be constant through each run
        min_irid = df.irid.min()
        max_irid = df.irid.max()
        print(f'for rid {rid}, irid range is {min_irid} to {max_irid}')
#        dfview = df['where irid is {iirid}']                   
        dfview = df[df.irid ==3]
        print(dfview)
        return
    
        xdt = np.zeros(len(df.t))
        xdt[1:] = df.t[1:].values - df.t[:-1].values
        xdt[0] = xdt[1]
        df['dt'] = xdt

        xfrac = np.zeros_like(xdt)
        xfrac[0] = 1
        xdfrac = np.zeros_like(xdt)
        xbfrag = np.zeros_like(xdt)
        for i in range(1, len(xdt)):
            xT = df.loc[i].Temp
            dt = xdt[i]
            thisloss = fracloss(dt, xT, newdH, newdS)
            xdfrac[i] = thisloss / xdt[i]    # a dissociation rate units 1/s
            cdfrac = 1-np.exp(-thisloss)
            fracfrac = random.uniform(0, 1)
            xbfrag[i] = xbfrag[i-1] or (fracfrac < cdfrac)
            xfrac[i] = xfrac[i-1] * (1-thisloss)
        df['frac'] = xfrac
        df['dfrac'] = xdfrac
        df['bfrag'] = xbfrag
        lfrac = df.iloc[-1].frac
        ldfrac = df.iloc[-1].dfrac
        print(f'{vin:4.0f}   {lfrac:0.3f}   {ldfrac:0.4f}')
#        query = f''UPDATE trajs set (frac,dfrac,bfrag)=(?,?,?) where rid=? and irid=? and t=?'
#        with sqlite3.connect(gdbname) as db:
#            db.executemany(query, (df.frac, df.dfrac, df.bfrag, df.rid, df.irid, df.t))

    return testcount


class quadratic3point():
    '''Do a 3-point quadratic fit given 3 x and 3 y values.
    x-points must be evenly spaced and ordered low to high.
    Also return x value for a minimum (max) and x values for y=const.
    Fit coeficients are from f(x) = ax^2 + bx + c
    '''
    def __init__(self, x3, y3):
        self.x3 = x3
        self.y3 = y3
        self.xmiddle = self.x3[1]    # keep middle x for reference
        self.d = x3[1] - x3[0]      # x spacing, assumed equal to x3[2]-x3[1]
        self.a = (y3[0] - 2*y3[1] + y3[2]) / (2 * self.d**2)
        self.b = (y3[2] - y3[0]) / (2 * self.d)
        self.c = y3[1]
        xshift = -self.b / (2 * self.a)     # minimum is shifted from middle pt by this much
        self.xc = self.xmiddle + xshift 
        self.yminc = self.a*xshift**2 + self.b * xshift + self.c
    
    def coef(self):
        '''return the coefficients.'''
        return self.a, self.b, self.c
    
    def eval(self, x):
        xshift = x - self.xmiddle
        return self.a*xshift**2 + self.b * xshift + self.c
    
    def xmin(self):
        '''return x value of the function minimum.'''
        return self.xc
    
    def ymin(self):
        return self.yminc
    
    def xcross(self, ymult):
        '''Solve the quadratic for x points where function crosses ymin * ymult
        ymin*ymult = ax^2 + bx + c
            or
        0 = ax^2 + bx + (c-ymin*ymult)
            or
        0 = ax^2 + bx + coff            with coff = (c-ymin*ymult)
        xcross = (-b +/- sqrt(b^2 - 4*a*coff)) / 2a
        The values returned are offsets from self.xmiddle
        '''
        coff = self.c - self.yminc * ymult
        xc1 = (-self.b - math.sqrt(self.b**2 - 4 * self.a * coff)) / (2 * self.a)
        xc2 = (-self.b + math.sqrt(self.b**2 - 4 * self.a * coff)) / (2 * self.a)
        return (xc1, xc2)


class fitclass():
    '''A fitclass holds any information needed to do a fit from a spa_heat model to experimental data.
    Specifically, it needs to store appropriate values for:
        Vins, tTd, rfrac arrays
        Can create a weights list by analyzing the experimental data and assigning weights.
        And can support computing the average squared difference between model and experiment.
    '''
    def __init__(self, Vins, tTd, expfracs):
        self.Vins = Vins
        self.efracs = np.array(expfracs)
        self.tTd = tTd
        self.weights = [1 for v in Vins]
        self.wtsum = np.sum(self.weights)
        self.hasprinted = False
        self.lastfittime = None
        self.opt = None
        self.Hpts = None
        self.Spts = None
        self.fit_fracs = np.array(expfracs)
        self.rescaled_fit = np.array(expfracs)
        self.rescaled_fracs = self.efracs
        self.popt = None
        self.pcov = None
        self.b = 25

    def setHSpts(self, Hpts, Spts):
        self.Hpts = Hpts
        self.Spts = Spts
        
    def sigmoid_fitting(self):
        '''Fit a generic sigmoid curve to the experimental data.
        This is done so we can use the baseline offset and maximum L to rescale the experimental data.
        The rescaling is intended to get a better estimate of the 100% parent ion abundance at low energies
        as well as account for non-zero baseline (chemical?) noise at high energies.
        '''
        def sigmoid (v, L, v0, k, b):
          # L: the maximum value of the function
          # v0: the value of v at the midpoint of the sigmoid
          # k: the steepness of the curve
          # b: the offset from zero
          return L / (1 + np.exp (-k * (v - v0))) + b
        
        xdata = self.Vins
        ydata = self.efracs
        try:                #if the sigmoid fit fails it will use exp data and CID50 for rescale/weighting rather than the sigmoid fit parameters
            self.popt, self.pcov = curve_fit(sigmoid, xdata, ydata, p0=[-1,30,.5,1])
            print(f'sigmoid popt: {self.popt}\n        pcov: {self.pcov}')
            
            self.fit_fracs = sigmoid(xdata, *self.popt)
            
            self.fit_fracs = np.where(self.fit_fracs <= 0, 0, self.fit_fracs)   # forces the sigmoid curve points > 0
            
            # now rescale the fit points to range from 0 to 1
            self.rescaled_fit = (self.fit_fracs - min(self.fit_fracs)) /(max(self.fit_fracs) - min(self.fit_fracs))
            self.rescaled_fit = np.where(self.rescaled_fit <= 0, 0, self.rescaled_fit)  # and ensure > 0
            self.b = self.popt[1]  #midpoint of gaussian
        except:
            print('Sigmoid fit failed. Using Experiment Data for rescaling and CID50 for combined weighting')
            self.rescaled_fit = (self.expfracs - min(self.expfracs)) /(max(self.expfracs) - min(self.expfracs))
            self.rescaled_fit = np.where(self.rescaled_fit <= 0, 0, self.rescaled_fit)  # and ensure > 0
            ibelow = np.where(self.rescaled_fit <= 0.5)[0][0]
            iabove = ibelow - 1
            slope = (self.rescaled_fit[ibelow] - self.rescaled_fit[iabove])/(Vins[ibelow] - Vins[iabove]) 
            intercept = self.rescaled_fit[iabove] - slope*Vins[iabove]
            CID50 = (0.5 - intercept)/slope
            self.b = CID50
            
            
            
    def rescalefracs(self):
        '''first pass as rescaling, can use it if you want to'''
        
        self.rescaled_fracs = (self.efracs - min(self.fit_fracs)) /(max(self.efracs) - min(self.fit_fracs))
        self.rescaled_fracs = np.where(self.rescaled_fracs <= 0 , 0, self.rescaled_fracs)
        # self.efracs = rescaled_fracs
        # self.rescaled_fracs = self.efracs

    def formula_weights(self, weighting='combined'):        
        '''There are three different weighting options: none, square, combined.
        none: all Vins are weighted at 100%.
        square: uses weights formula, Vins with weights over 10% get rounded to 100%, otherwise 0%,
                then ensures that there are at least 10 non-zero weight values,
        combined: creates a 'square' weighting for roughly points between 95-5 percent,
                  falling off gaussian centered at b (CID50/Sigmoid fit center) w/ 1 std being
                  3 times the difference between b and the first Vin where the weights fall under 95%.'''
        
        # assume efracs range is from 0 to 1
        # we want 0->0, 0.5->1, 1.0->0
        # weights = 1 - 4*(self.rescaled_fit - (1/2))**2  
        weights = 1 - 16*(self.rescaled_fit - (1/2))**4  
        
        
        
        if weighting == 'none':
            weights = np.ones(len(self.Vins))
            print(f'CID50 (roughly): {self.b}')
        if weighting == 'square':
            cutoff = 0.1
            weights = np.where(weights >= cutoff, 1, weights) # This creates a 'square' weighting, it tends to be pretty strict
            weights = np.where(weights < cutoff, 0, weights)
            
            ilist = np.where(weights >= 0.5)[0]
            num_weights = len(ilist)
            while num_weights <= 10:
                weights[ilist[0] - 1] = 1
                weights[ilist[-1] + 1] = 1
                ilist = np.where(weights >= 0.5)[0]
                num_weights = len(ilist)
            print(f'CID50 (roughly): {self.b}')
        if weighting == 'combined':
            
            percent = 0.95
            
            
            c = (self.b - self.Vins[np.where(self.rescaled_fit <= percent)[0][0]] + 1)*3
            gaussian = np.exp(-(1/2)*((self.Vins - self.b)/c)**2)
            # plt.plot(self.Vins, gaussian)
            cutoff = 0.05
            weights = np.where(weights >= cutoff, 1, gaussian) # This creates a 'square' weighting for roughly points between 90-10 percent, 
            print(f'CID50 (roughly): {self.b}')
            
        
        
      
        self.weights = weights
        
        return self.weights

    def diffsq(self, paramv, dH=None, printvals=False, savers=False): #this is the diffsq that does the fitting for experiment data
        '''Uses Vins, times, Temps, rfrac arrays from class context.
        Also uses weights for each point (Vin value).
        Calculates fracremains for each Vin in Vins list.
        Computes the weighted sum of square difference between experimental and computed fractions.
        Returns weighted average.
        '''
        if savers: rs = []
        if not dH:
            ddH = paramv[0] * 1000      # convert to J/mol from kJ/mol
            ddHns = paramv[0]
            ddS = paramv[1]
        else:
            ddH = dH * 1000
            ddHns = dH
            ddS = paramv[0]
        
        sumsq = 0.0
        sumwt = 0.0
        efracs = self.rescaled_fracs    # for local use for the fit
        # new_fracs = list(np.array(fracs) - min(fracs))/(max(fracs) - min(fracs))  #First pass as rescaling, can use it if you want to
        for Vin,efrac,weight in zip(self.Vins, efracs, self.weights):
            if weight <= 0.01:          #if the weight is less than 1%, skip the calculation to speed up fitting times
                if savers:
                    rs.append(efrac)
                continue
            if usenewfrac:
                dfv = self.tTd[Vin]
                finalr = computefracs(dfv, ddHns, ddS)
#                print(f'computefracs: dH,dS ({ddHns:0.2f}, {ddS:0.2f}, frac {finalr:0.3f})')
            else:
                times, Temps, tmax = self.tTd[Vin]
                finalr = fracremains(times, Temps, ddH, ddS, maxtime=tmax)
#                print(f'oldfrac: dH,dS ({ddHns:0.2f}, {ddS:0.2f}, frac {finalr:0.3f})')

            sumsq += weight*(efrac - finalr)**2
            sumwt += weight
            if savers: rs.append(finalr)
            if printvals and not self.hasprinted:
                print(f'   Vin: {Vin:4.1f}  tmax: {tmax*1000:4.2f} ms   exp: {efrac:0.2f} --> calcr: {finalr:0.2f}')
        self.hasprinted = True

        if savers:
            self.rs = rs            # save for plotting any time we do this
        
        return sumsq / sumwt

    def setfitdirec(self, d_array):
        self.d_array = d_array

    def fit(self, dH, dS):
        '''Do curve fit with starting point [dH, dS].'''
        self.lastfittime = time()
        param0 = [dH, dS]                    ############################ pick a better starting value
        #res = minimize(self.diffsq, param0, method='Nelder-Mead', tol=2e-2)
        darray = self.d_array.copy()
        res = minimize(self.diffsq, param0, method='Powell', tol=2e-5, options={'direc':darray})
        self.lastfittime = time() - self.lastfittime
        return res

    def fitH(self, dH, dS):
        '''Do curve fit with starting point [dH, dS]. Only vary dH for a 1-D fit.'''
        self.lastfittime = time()
        param0 = [dS]                    ############################ pick a better starting value
        res = minimize(self.diffsq, param0, args=(dH,), method='Nelder-Mead', tol=0.01)
        self.lastfittime = time() - self.lastfittime
        return res

    def get3lowpoints(self, axisvec, step, dH0, dS0):
        '''Start with 2 points, add extra points until middle one is a minimum.'''
        dHa = [dH0]
        dSa = [dS0]
        fa = [self.diffsq([dHa[-1], dSa[-1]])]
        imin = 0
        dHa.append(dH0 - axisvec[0]*step)
        dSa.append(dS0 - axisvec[1]*step)
        fa.append(self.diffsq([dHa[-1], dSa[-1]]))
        if fa[0] > fa[1]:
            imin = 1
            famin = fa[1]
        else:
            imin = 0
            famin = fa[0]

#        print('start set of 2 with imin {imin}')
#        for h,s,f in zip(dHa,dSa,fa):
#            print(f'   [{h:0.3f}, {s:0.3f}]: {100*f:0.3f}')

        while len(fa) < 3:
            if imin == 0:
                dH = dHa[0] + axisvec[0]*step
                dS = dSa[0] + axisvec[1]*step
                fav = self.diffsq([dH, dS])
#                print(f'eval0 [{dH:0.3f},{dS:0.3f}]: {fav*100:0.3f}')
                dHa.insert(0, dH)
                dSa.insert(0, dS)
                fa.insert(0, fav)
                if fav < famin:     # the three are ordered from low to high, drop the last and try another step
                    dHa.pop(-1)
                    dSa.pop(-1)
                    fa.pop(-1)
                    famin = fav
                    # imin is still 0
                else:
                    imin = 1        # we now have 3 points with min in the center
            else:   # imin must be 1 here, insert at end
                dH = dHa[1] - axisvec[0]*step
                dS = dSa[1] - axisvec[1]*step
                fav = self.diffsq([dH, dS])
#                print(f'eval1 [{dH:0.3f},{dS:0.3f}]: {fav*100:0.3f}')
                dHa.insert(2, dH)
                dSa.insert(2, dS)
                fa.insert(2, fav)
                if fav < famin:     # the three are ordered from high to low, drop the first and try another step
                    dHa.pop(0)
                    dSa.pop(0)
                    fa.pop(0)
                    famin = fav
                else: 
                    imin = 1        # we now have 3 points with min in center
        return dHa, dSa, fa

    def fits(self, dHstart, dSstart, Tstart=500, step=10):
        '''Do fit assuming well-behaved eliptical region using multiple lines across the region.
        Note: while this works for the minor direction scans, the major scan seems to cross the minimum 
        too quickly as if it is not properly aligned.'''
        self.lastfittime = time()
        dH = dHstart
        dS = dSstart
        Tfit = Tstart
        kTfit = Tfit / 1000     # convenient to use kK (1000 K) for temp scale
        minorlen = step
        majorlen = minorlen / kTfit
        self.dH = dH
        self.dS = dS

        # debugging plots
#        fHS = plt.figure()
#        fHS.add_subplot(111)
#        axHS = fHS.get_axes()[0]
#        ff = plt.figure()
#        ff.add_subplot(111)
#        axff = ff.get_axes()[0]
        px = []
        py = []
        pf = []

        print(f'Tfit starts as {Tfit:0.2f} {kTfit:0.5f}\n')

        # start with dHstart, dSstart, evaluate 5 points on line perpendicular to major axis defined by Tstart
        major = np.array([kTfit, 1]) / math.sqrt(1 + kTfit**2)    # unit vector along assumed major axis
        minor = np.array([-major[1], major[0]])                 # unit vector along minor axis
#        print(f'major: [{major[0]:0.3f}, {major[1]:0.3f}]    minor: [{minor[0]:0.3f}, {minor[1]:0.3f}]')


        maxratio = 1.1
        # first scan with assumed Tfit line
        # repeat this with different steps until ratio between the max(fa) and min(fa) < 2
        dH, dS = dHstart, dSstart
        while True:
            dHa, dSa, fa = self.get3lowpoints(minor, step, dH, dS)
            px.extend(dHa); py.extend(dSa); pf.extend(fa)       # debugging plots
            faratio = max(fa) / min(fa)
            dH, dS, fav = dHa[1], dSa[1], fa[1]
            if faratio < maxratio:
                break
            else:
                step = step / 2

        # now get quad approx to best point
        quadobj = quadratic3point(dHa, fa)
        dH = quadobj.xc
        dS = dSa[1] +  (dH-quadobj.xmiddle)*minor[1]/minor[0]
        # save scan0 point
        dHscan0, dSscan0 = dH, dS
        favd0 = self.diffsq([dHscan0,dSscan0])
        px.append(dHscan0); py.append(dSscan0); pf.append(favd0)       # debugging plots
        print(f'1st scan: [{dHscan0:0.3f}, {dSscan0:0.3f}] -> {fav*100:0.3f}, {favd0*100:0.3f}')



        # now shift down the major axis and do another scan
        dH = dHscan0 - major[0]*step*10
        dS = dSscan0 - major[1]*step*10
#        print(f'moved scan from [{dHscan0:0.3f}, {dSscan0:0.3f}] to [{dH:0.3f}, {dS:0.3f}]')
        while True:
            dHa, dSa, fa = self.get3lowpoints(minor, step, dH, dS)
            px.extend(dHa); py.extend(dSa); pf.extend(fa)       # debugging plots
            faratio = max(fa) / min(fa)
#            print(f'found set of 3 with center minimum, step {step:0.2f} ratio {faratio:0.2f}')
#            for h,s,f in zip(dHa,dSa,fa):
#                print(f'   [{h:0.3f}, {s:0.3f}]: {100*f:0.3f}')
            dH, dS, fav = dHa[1], dSa[1], fa[1]
            if faratio < maxratio:
                break
            else:
                step = step / 2

        quadobj = quadratic3point(dHa, fa)
        dH = quadobj.xc
        dS = dSa[1] +  (dH-quadobj.xmiddle)*minor[1]/minor[0]
        # save scan1 point
        dHscan1, dSscan1 = dH, dS
        favd1 = self.diffsq([dHscan1,dSscan1])
        px.append(dHscan1); py.append(dSscan1); pf.append(favd1)       # debugging plots
        # now with 2 scans at starting Tfit minor axis, we can compute a corrected Tfit
        kTfit = (dHscan1 - dHscan0) / (dSscan1 - dSscan0)
        Tfit = 1000 * kTfit
        self.Tfit = Tfit

        print(f'2nd scan: [{dHscan1:0.3f}, {dSscan1:0.3f}] -> {fav*100:0.3f}, {favd1*100:0.3f}   Tfit={Tfit:0.2f}')



        major = np.array([kTfit, 1]) / math.sqrt(1 + kTfit**2)    # unit vector along assumed major axis
        minor = np.array([-major[1], major[0]])                 # unit vector along minor axis
#        print(f'major: [{major[0]:0.3f}, {major[1]:0.3f}]    minor: [{minor[0]:0.3f}, {minor[1]:0.3f}]')

        # and now we can scan along major axis
        # now shift down the new major axis and do another scan along major axis this time
        if favd1 < favd0:
            dH = dHscan1 - major[0]*step*10
            dS = dSscan1 - major[1]*step*10
        else:
            dH = dHscan1 + major[0]*step*5
            dS = dSscan1 + major[1]*step*5
#        print(f'scanning major axis from [{dH:0.3f}, {dS:0.3f}]')
        while True:
            dHa, dSa, fa = self.get3lowpoints(major, step*10, dH, dS)
            px.extend(dHa); py.extend(dSa); pf.extend(fa)       # debugging plots
            faratio = max(fa) / min(fa)
            dH, dS, fav = dHa[1], dSa[1], fa[1]
            if faratio < maxratio:
                break
            else:
                step = step / 2


        quadobj = quadratic3point(dHa, fa)
        dH = quadobj.xc
        dS = dSa[1] +  (dH-quadobj.xmiddle)*minor[1]/minor[0]
        # save scan0m point
        dHscan0m, dSscan0m = dH, dS
        favd = self.diffsq([dHscan0m,dSscan0m])
        px.append(dHscan0m); py.append(dSscan0m); pf.append(favd)       # debugging plots
        kTfit = (dHscan0m - dHscan1) / (dSscan0m - dSscan1)
        Tfit = 1000 * kTfit
        self.Tfit = Tfit
        print(f'major scan: [{dHscan0m:0.3f}, {dSscan0m:0.3f}] -> {fav*100:0.3f}, {favd*100:0.3f}   Tfit={Tfit:0.2f}')

        majorlenv = quadobj.xcross(3);  majorlen = abs(majorlenv[1]-majorlenv[0])
        majorlen10v = quadobj.xcross(10);   majorlen10 = abs(majorlen10v[1]-majorlen10v[0])
        print(f'major len estimates at 3sd and 10sd: {majorlen:0.2f}, {majorlen10:0.2f}')
 
        major = np.array([kTfit, 1]) / math.sqrt(1 + kTfit**2)    # unit vector along assumed major axis
        minor = np.array([-major[1], major[0]])                 # unit vector along minor axis
 


        dH, dS = dHscan0m, dSscan0m
#        print(f'scanning minor axis from [{dH:0.3f}, {dS:0.3f}]')
        while True:
            dHa, dSa, fa = self.get3lowpoints(minor, step*2, dH, dS)
            px.extend(dHa); py.extend(dSa); pf.extend(fa)       # debugging plots
            faratio = max(fa) / min(fa)
            dH, dS, fav = dHa[1], dSa[1], fa[1]
            if faratio < maxratio:
                break
            else:
                step = step / 2
        quadobj = quadratic3point(dHa, fa)
        dH = quadobj.xc
        dS = dSa[1] +  (dH-quadobj.xmiddle)*minor[1]/minor[0]
        minorlenv = quadobj.xcross(3);  minorlen = abs(minorlenv[1]-minorlenv[0])
        minorlen10v = quadobj.xcross(10);   minorlen10 = abs(minorlen10v[1]-minorlen10v[0])
        print(f'minor len estimates at 3sd and 10sd: {minorlen:0.2f}, {minorlen10:0.2f}')


       # debugging plots
#        axHS.plot(px, py, '.-')
#        axff.plot(px, pf, '.-')
#        plt.show()

        self.Tfit = Tfit
        self.majorlen = majorlen
        self.minorlen = minorlen
        self.majorlenx = majorlen10
        self.minorlenx = minorlen10
        self.lastfittime = time() - self.lastfittime
        print(f'fits() completed in {self.lastfittime:0.2f} s\n')
        return dH, dS, Tfit, majorlen, minorlen

fitter = None

def main(celld, iond, dH, dS, Vins, fracs, title='', outfile=None, plotfile=None, weighting='combined', tpoints=0, updatedb=False, Tfactor=1.0):
    '''Show predicted fraction remaining vs. Vin for given inputs.
    Plot results with exp data if present.'''
    global fitter

    print(f'main(celld, iond, Vins)\n  celld: {celld}\n  iond: {iond}\n  dH: {dH}  dS: {dS}\n  Vins: {Vins}\n  fracs: {fracs}')

    # Now just get the time/temperature data for each Vin, and calculate the frac_remaining vs Vin.
    vtTd = {}   # {5.0: [[times...], [temps...]], 30.0: [[times...], [temps...]], ...}
    cid = find_cid(gdbname, celld)
    iid = find_iid(gdbname, iond)
    if not (cid and iid):
        print(f'cid {cid}    iid {iid}  error')
        exit(1)
    
    # Tfactors = [0.95, 1.0]
    # for Tscale in Tfactors:
        
    # read data from DB and build the tTd structure
    tTd = {}
    tTd95 = {}
    for indx, Vin in enumerate(Vins):
        rid, tmax = find_rid(gdbname, cid, iid, Vin)

        # if usenewfrac:
        dfv, dfv95 = get_tzT(gdbname, rid, tpoints)
        tTd[Vin] = dfv
        tTd95[Vin] = dfv95
        # else:
        #     dd = get_rundat(gdbname, rid, 'Temp')
        #     times, Temps = np.array(average_data(dd, 'Temp'))
        #     tTd[Vin] = [times, Temps, tmax]
            
        # max_T = compute_max_T(dfv)
            
    
    for i, dic in enumerate([tTd95, tTd]):
        # print(dic)
        # construct the fitter object 
        # do any rescaling and weights calculations and get the initial average squared difference
        fitter = fitclass(Vins, dic, fracs)
        fitter.sigmoid_fitting()
        fitter.rescalefracs()
        fitter.formula_weights(weighting=weighting)
        avesqdiff = fitter.diffsq([dH, dS])
        print(f'initial R2 with class: {avesqdiff:0.3f}')
    
        # try the fits function:
        dH, dS, Tfit, majorlen, minorlen = fitter.fits(dH, dS, 420)
        kTfit = Tfit / 1000
        # now starting with best fit from fits()
        Hpts = []
        Spts = []
    
        direc1 = np.array([kTfit, 1]) / math.sqrt(1 + kTfit**2)
        direc2 = np.array([-direc1[1], direc1[0]])
        direc = np.array([direc1*majorlen, direc2*minorlen])
        fitter.setfitdirec(direc)
    
        # do the curve fit, (fitter calls the ionspa fracremains many times as needed)
        res = fitter.fit(dH, dS)
        dH, dS = res.x
        Hpts.append(dH);    Spts.append(dS)
    
        print(f'class minimize:      time:{fitter.lastfittime:0.2f} s     dH:{dH:0.3f}, dS:{dS:0.3f}     R2:{res.fun:0.2e}')
        print(res)
    
    
        # redo a sweep with the best fit values for dH, dS
        avesq = fitter.diffsq([dH, dS], savers=True)
        rs = fitter.rs
        fitter.setHSpts(Hpts, Spts)
    
        Tfit = fitter.Tfit
        dGfit = dH - Tfit * dS/1000
        majorlen = fitter.majorlen
        minorlen = fitter.minorlen
        majorlenx = fitter.majorlenx
        minorlenx = fitter.minorlenx
    
        kints = []
        max_Ts_Vins = []
        max_ts_Vins = []
        Teffs = []
        for indx, Vin in enumerate(Vins):
            dfv = dic[Vin]
            kint, max_T, max_t = compute_rate_integral(dfv, dH, dS)
            
            args = (dH, dS, kint, max_t)
            Teff = fsolve(effective_temp, max_T, args=args, full_output=True)[0][0]
            
            kints.append(round(kint,2))
            max_Ts_Vins.append(round(max_T,2))
            max_ts_Vins.append(max_t)
            Teffs.append(round(Teff,2))
            
        T_harm = stats.hmean(max_Ts_Vins)               #Computes several different temperature values to use for dG calculations
                                                        #T_harm is harmonic average temperature of all Vins
        CID50, indx50 = get_CID50(Vins, fitter.rescaled_fracs)

        
        CIDlife, indxlife = get_lifetime(Vins, fitter.rescaled_fracs)

        
        T_cid50_eff = Teffs[indx50]
        
        dG = dH - dS * T_cid50_eff/1000
        avesq = 1 - avesq
        
        #indx 50 is the Vin index that is immediately below the CID50 point
        
        old_T_cid50_eff = Teffs[indx50]
        T_cid50_eff = np.interp(fitter.b, Vins, Teffs)
        
        print()
        print(f'Old way Teff: {old_T_cid50_eff:3.2f} K at Vin: {Vins[indx50]} V')
        print(f'New way Teff: {T_cid50_eff:3.2f} K at Vin: {fitter.b:3.2f} V', flush=True)
        print()
        T_life = max_Ts_Vins[indxlife]
        T_cid50 = max_Ts_Vins[indx50]
        T_eff_harm = stats.hmean(Teffs)
        
        teff = max_ts_Vins[indx50]
        if i == 1:
            print('Finished with the first fit', flush=True)
            dH100 = dH
            dS100 = dS
            dG100 = dG
        if i == 0:
            dH95 = dH
            dS95 = dS
            dG95 = dG
            print('Finished with the assumed Temperature error fit', flush=True)
        
        
    errdH = abs(dH100 - dH95)
    errdS = abs(dS100 - dS95)
    errdG = abs(dG100 - dG95)
        
    
    kint = kints[indx50]
    Rcid50 = rs[indx50]
    outputs = error_calc(kint, T_cid50_eff, dG100, teff, fitter.rescaled_fracs, rs, fitter.weights, Rcid50, dS100)
    
    kint, Rcid50, kbt_dG, dR2, top_sum, bot_sum, rel_dk_squared, dk, rel_ddG_squared, T_err, k_err, ddG, ddH, ddS = outputs
    # dG_harm = dH - T_harm * dS/1000
    # dG_cid50 = dH - T_cid50 * dS/1000
    # dG_life = dH - T_life * dS/1000
    
    updateargs = (max_ts_Vins, max_Ts_Vins, Teffs, fracs, rs, fitter.rescaled_fracs, fitter.weights, kints, 
                  round(dH100,2), round(dS100,2), round(dG100,2), avesq, round(T_cid50_eff,2), round(T_harm,2), round(T_cid50,2), round(T_life,2), round(T_eff_harm,2), round(Tfit,2), round(CID50,2))
    
#    if updatedb:
#        update_db(gdbname, cid, iid, Vins, *updateargs)    
#        print('Updated databases with fit results')

    
    rate_constant = lossrate(T_cid50_eff, dH100, dS100)

    print(dH95)
    print(dS95)
    print(dG95)
    outhead = 'infile, dH, errdH, dS, errdS, dG, errdG, R2, Teff, rate_constant, Rcid50, CID50, effective_time'
    outline = f'{title}, {dH100:0.2f}, {errdH:0.4g}, {dS100:0.2f}, {errdS:0.4g}, {dG100:0.2f}, {errdG:0.4g}, {avesq:0.2e}, {T_cid50_eff:0.2f}, {rate_constant:0.4g}, {Rcid50:0.4g}, {fitter.b:0.2g}, {teff:0.4g}'

    if outfile:
        with open(outfile, 'a') as fh:
            fh.write(outhead + '\n')
            fh.write(outline + '\n')

    if plotfile:
        plt.plot(Vins, rs, label='fit')
        plt.plot(Vins, fitter.weights, label='weights')
        plt.plot(Vins, fracs, 'x', label='exp fracs')
        plt.plot(Vins, fitter.rescaled_fracs, 'o', label='rescaled fracs')
        plt.xlabel('Vin [V]')
        plt.ylabel('frac remaining')
        title += f'dH={dH100:0.2f},  dS={dS100:0.2f}   R2:{avesq:0.3f}'
        plt.title(title)
        plt.legend()
        if plotfile == 'screen':
            plt.show()
        else:
            plt.savefig(plotfile)

    print(outhead)
    print(outline)

    return dH100, dS100, tTd


def plotfitregion(newdH, newdS, tTd, title='', rdH=None, rdS=None, plotmap='screen', maprange=None, doscan=False):
    '''Plot contours of residual^2 in the region around the best fit.
    Range is from the minimum value to 10x the minimum.'''

    global fitter   # so we can use diffsq() and retrieve sample points
    Tfit = fitter.Tfit
    kTfit = Tfit / 1000

    majorlen = fitter.majorlen
    minorlen = fitter.minorlen
    majorlenx = fitter.majorlenx
    minorlenx = fitter.minorlenx
    dH, dS = newdH, newdS

    major = np.array([kTfit, 1]) / math.sqrt(1 + kTfit**2)
    minor = np.array([-major[1], major[0]])
    direc = np.array([major*majorlen, minor*minorlen])


#    dHlist = [newdH + delta for delta in np.arange(-35, 35, 0.3)]
#    dSlist = [newdS + delta for delta in np.arange(-50, 50, 0.3)]
    if doscan:
        xds = newdS
        Slist = [xds + delta for delta in np.arange(-20, 20, 0.3)]

        resS = np.full_like(Slist, np.nan)
        for i, dS in enumerate(Slist):
            r = fitter.diffsq([newdH,dS], savers=False)
            if r < 0.01:
                resS[i] = r
        plt.figure()
        plt.plot(Slist, resS)
        plt.show()
        print('scan plot done')
    else:
        pass

    # approximate ellipse around residual ellipse to limit number of points calculated below
    ##### Better plan is to sample a few lines across the ellipse and plot the ellipses directly from the functional form.
    ##### This should work well since the residuals make a very well-defined elliptical shape.

#    v0 = fitter.d_array[0]
#    v1 = fitter.d_array[1]
#    size = 50
#    major = v1 * size
#    minor = v0 * size

    ctr = [newdH, newdS]
    a = majorlen        # wild guess to keep ellipse minor axis a better size. major/minor may be misnamed.
    b = minorlen
    theta = math.atan2(major[1], major[0])
#    def inellipse(point):
#        '''Calculate if point is within ellipse defined by center and major/minor axis vectors.'''
#        x = point[0] - ctr[0]
#        y = point[1] - ctr[1]
#        x_rot = x * math.cos(theta) + y * math.sin(theta)
#        y_rot = -x * math.sin(theta) + y * math.cos(theta)
#        return (x_rot/a)**2 + (y_rot/b)**2 <= 1

    angles = np.linspace(0,2*math.pi,73)
    pdH = np.zeros_like(angles)
    pdS = np.zeros_like(angles)
    pxdH = np.zeros_like(angles)
    pxdS = np.zeros_like(angles)


    ptheta = theta          # redefine theta from above
    print(f'theta: {theta*180/3.14:0.1f}   ptheta: {ptheta*180/3.14:0.1f}')
    for i,angle in enumerate(angles):
        # 3R ellipse
        ddH = math.cos(angle)*majorlen/2
        ddS = math.sin(angle)*minorlen/2
        pdH[i] = newdH + ddH * math.cos(ptheta) - ddS * math.sin(ptheta)
        pdS[i] = newdS + ddH * math.sin(ptheta) + ddS * math.cos(ptheta)
        # 10R ellipse
        ddH = math.cos(angle)*majorlenx/2
        ddS = math.sin(angle)*minorlenx/2
        pxdH[i] = newdH + ddH * math.cos(ptheta) - ddS * math.sin(ptheta)
        pxdS[i] = newdS + ddH * math.sin(ptheta) + ddS * math.cos(ptheta)
    

    ##### the ellipse values above can be used to find actual values for major/minor axes by doing 1D scans across the 
    ##### approximate major/minor axes, then updating the vectors, and theta, a, b values.
    ##### First scan across center along the minor axis direction to verify the lowest point.
    #####     This can return a good approximation for the minor axis minR2*2, minR2*3, minR2*10 lengths.  
    ##### Then offset by dH +/- 20 and corresponding dS range and do another scan to get an accurate major axis vector.
    ##### Then scan along the major axis to get minR2*2, minR2*3, minR2*10 lengths.  
    ##### Now simply plot those ellipses directly and avoid the thousands of calls to diffsq().
    
    # map response for the full space
#    pointstoplot = 0
#    resA = np.ones((len(dHlist), len(dSlist)))
#    totpoints = len(dHlist) * len(dSlist)
#    dHplotted = []
#    dSplotted = []
#    for i, dH in enumerate(dHlist):
#        for j, dS in enumerate(dSlist):
#            pt = [dH,dS]
#            if inellipse(pt):
#                resA[i,j] = fitter.diffsq(pt, savers=False)
#                dHplotted.append(dH)
#                dSplotted.append(dS)
#                pointstoplot += 1
#    print(f'plot count: {pointstoplot} out of {totpoints}')

#    resAmin = resA.min()
#    resAc = resA # - resAmin
#    resAc = np.clip(resAc, resAmin, 10*resAmin)

    if plotmap:
        plt.figure()
#        levels = np.arange(1, 10, 1.0) * resAmin
#        plt.plot(dHplotted, dSplotted, ',', label='calc locations')
#        plt.contour(dHlist, dSlist, resAc.T, 10)    # label doesn't work here!
        plt.plot([newdH], [newdS], 'o', label='best fit point')
        if rdH and rdS:
            plt.plot([rdH], [rdS], 'x', label='reported BIRD point')
        plt.plot(pdH, pdS, '-r', label='ellipse3R')
        plt.plot(pxdH, pxdS, '-b', label='ellipse10R')
        if maprange and len(maprange)==4:
            plt.xlim(maprange[0], maprange[1])
            plt.ylim(maprange[2], maprange[3])
        plt.xlabel('dH')
        plt.ylabel('dS')
        plt.title(title)
        plt.legend()
        if plotmap == 'screen':
            plt.show()
        else:
            plt.savefig(plotmap)


if __name__ == '__main__':
    fitter = None
    parser = argparse.ArgumentParser(
        prog='python spa_db_fit.py',
        description='''Given a .json file (as from spa_mpirun), and either the dhds_start from
        that file or specified dH and dS, get a best fit to the data and plot the fraction remaining vs Vin.
        Optionally, write results to an output file and/or update the fracs values in the DB.''')
    parser.add_argument('jsonfile', nargs='+', help='json input files to specify ion, cell and other parameters.')
    parser.add_argument('--dH', help='overrides the initial dH value from input file') #   , nargs=1)
    parser.add_argument('--dS', help='overrides the initial dS value from input file') #   , nargs=1)
    parser.add_argument('--plotfile', default='screen', help='normally show interactive plot, otherwise to file')
    parser.add_argument('--plotmap', default=None, help='also plot residuals around best point')   # add this option to show plot of fit vs dH,dS
    parser.add_argument('--db', default = 'spa_db.sqlite3', help='use this for the database file')
    # parser.add_argument('--updatefracs', action = 'store_true', help='update the fracs values in the DB')
    parser.add_argument('--outfile', help='if given, write fit values to this file.')
    parser.add_argument('--maprange', nargs=4, type=float, default=None)
    parser.add_argument('--newfrac', default=True, action='store_false', help='use new fracremains version.')
    parser.add_argument('--weighting', default='combined', help='choose what type a fit weights to use: none, square, combined')
    parser.add_argument('--tpoints', default=0, help='number of time points to interpolate, default is no interpolation')
    parser.add_argument('--updatedb', default=False, help='update the db with the results of the fit')
    args = parser.parse_args()
    gdbname = path.splitext(args.db)[0] + '.sqlite3'
#    print(args)

    print(f'args taken from {args.jsonfile} files in order.')
    jfile = args.jsonfile
    paramd = loadjfile(jfile)

    expd = paramd['exp']
    celld = paramd['cell']
    cellVoffset = celld.get('Voffset', 0.0)     # read the offset voltage to be added to Vins, default is 0.0
    iond = paramd['ion']
    newtemp = celld.get('T0offset', 0.0) + iond.get('T0')
    iond.update({'T0' : newtemp})
    dH = float( args.dH or paramd['dhds_start'][0] )
    dS = float( args.dS or paramd['dhds_start'][1] )
    plotfile = args.plotfile
    plotmap = args.plotmap
    maprange = args.maprange
    print(f'plotfile: {plotfile}')
    print(f'plotmap: {plotmap}   range: {maprange}')
    # updatefracs = args.updatefracs
    usenewfrac = args.newfrac
    tpoints = int(args.tpoints)
    updatedb = args.updatedb
    print(f'newfrac: {usenewfrac}     tpoints: {tpoints}')
    # print(f'updatefracs: {updatefracs}')
    outfile = args.outfile      # may be None if not specified
        # the main function might choose to write a simple results output to this output file
        # format is TBD
    weighting = args.weighting
#    print(jfile, dH, dS)

    Vins = expd.get('Vins', None)
    fracs = None
    if not Vins:
        print(expd)
        expdatpath = path.join(expd['folder'], expd['file'])
        print(expdatpath)
        gexp = pd.read_table(expdatpath, delimiter=',', header=None, names=['Vin', 'efrac'])
        Vins = gexp.Vin
        fracs = list(gexp.efrac.values)
    Vins = [vin + cellVoffset for vin in Vins]
    print(Vins)
    print(fracs)
    title = f'{jfile} '
    newdH, newdS, tTd = main(celld, iond, dH, dS, Vins, fracs, title, outfile, plotfile, weighting, tpoints, updatedb)
    rdH = iond.get('refdH', None)
    rdS = iond.get('refdS', None)

    if plotmap:
        plotfitregion(newdH, newdS, tTd, title, rdH, rdS, plotmap, maprange)

    # if updatefracs:
        # ucount = doupdate(celld, iond, newdH, newdS)
        # print(f'{ucount} records updated. Implementation is not done yet.')
