####################################
#
# spa_db_refrag.py
#
# Original work by Ken Newton, et.at. at University of Oregon
# Copyright (c) 2024 James S. Prell, all rights reserved
#
# non-MPI script to run ionspa system for given ion with given cell over a range of
# input energies and, with given deltaH and deltaS, predict the fragmentation curve.
#
# with a specified dH and dS, recalculate thee columns in the database table for
# frac, dfrac, bfrag using the ion temperatures and times from an earlier run.
# This allows one to do a standalone recalculation of just the fragmentation
# estimates without rerunning all the ion trajectories.
#
####################################

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from os import path
from ionspa import fracloss, fracremains, lossrate, loadjfile, printparams
import spa_db_fit as spafit
from time import sleep, time
import random
import sqlite3
import argparse
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})


'''
This program runs without the mpi system so will normally be called as:
python spa_db_refrag.py --db {database_name} {input_file_name} --xaxis {xaxis_option} --yaxis {yaxis_option} --dH {dH_value} --dS {dS_value}

DB table creation:
CREATE TABLE IF NOT EXISTS "cells" ("cid" integer primary key,"type" text,"L" float,"T" float,"rho" float,"Mgamu" float,"z0" float,"otherjson" text);
CREATE TABLE IF NOT EXISTS "ions" ("iid" integer primary key, "name" text, "mass" float, "charge" float, "CCS" float, "num_atoms" integer, "T0" float, "hcprofname" text);
CREATE TABLE IF NOT EXISTS "runs" ("rid" integer primary key, "cid" integer, "iid" integer, "Vin" integer, "tmax" float);
#CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float);
#CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float, "collnum" float, "frac" float, "dfrac" float, "bfrag" float
CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float, "collnum" float, "frac" float, "dfrac" float, "bfrag" float, "ma" float);
'''


gdbname = 'spa_db.sqlite3'
gdH = 80
gdS = -120
gxlabel = 't'
gylabel = 'frac'

def get_iridminmax(dbname, rid):
    '''Return minimum and maximum irid values for a given rid in the trajs table.'''
    with sqlite3.connect(dbname) as db:
        query = f'SELECT min(irid), max(irid) from trajs where rid = {rid}'
        cur = db.execute(query)
        fdat = cur.fetchall()
        return fdat[0]

df = None
lastrid = -1
def get_fracdat(dbname, rid, irid):     # rid different for each vin, irid for each ion in the group
    '''Return a dataframe with data relevant to recomputing the frac, dfrac and bfrag columns.'''
    global df, lastrid
    with sqlite3.connect(dbname) as db:
#        query = f'SELECT irid, t, z, Temp, frac, dfrac, bfrag FROM trajs where rid={rid} and irid = {irid}'
        if rid != lastrid:
            query = f'SELECT irid, t, z, Temp FROM trajs where rid={rid}'
            lastrid = rid
            df = pd.read_sql_query(query, db)
        # OK but returns nan in rows with where clause is false
    dfv = df[df.irid==irid].values
    return dfv


def doupdate(celld, iond, newdH, newdS):
    '''
    Recalculate the frac remaining for each ion
    ---- not implemented: store back into the database.

    This is really used to get plots showing predicted fragmentation
    '''
    updatecount = 0

    print(f'recalculating fragmentation for dH {newdH}, dS {newdS}')

    cid = spafit.find_cid(gdbname, celld)
    iid = spafit.find_iid(gdbname, iond)
    if not cid:
        print(f'cid {cid} not found')
    if not iid:
        print(f'iid {iid} not found')
    if not (cid and iid):
        print(f'Error in locating run data for cell and ion.')
        exit(1)

    # get list of Vins
    query = f"SELECT rid, Vin, tmax from runs where cid={cid} and iid={iid};"
    fdat = None
    with sqlite3.connect(gdbname) as db:
        cur = db.execute(query)
        fdat = cur.fetchall()
        # fdat is a list of 2-tuples, each pair is rid, Vin
    # now extract Vin list from fdat (associated with rid?)

    Vins = np.array([vin for rid,vin,tmax in fdat])
    tmaxs = np.array([tmax for rid,vin,tmax in fdat])
    cfracs = np.zeros_like(tmaxs)   # will hold final parent fraction at end of each trajectory
    gfracs = np.zeros_like(tmaxs)

    tcount = 0
    for rid, vin, tmax in fdat:     # loop over vins
        min_irid, max_irid = get_iridminmax(gdbname, rid)

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

#        print(f'for rid {rid} Vin = {vin}, irid range is {min_irid} to {max_irid}')
        # sums for this vin
        slfrac = 0
        sfrac = None
        sdfrac = None
        sbfrag = None
        sTemp = None
        sz = None

#        starttime = time()

        for ionnum in range(min_irid, max_irid+1):
            npa = get_fracdat(gdbname, rid, ionnum)
            # columns are irid, t, z, Temp
            cdct = dict(irid=0,  t=1, z=2, Temp=3)
            dft = npa[:,cdct['t']]
            dfT = npa[:,cdct['Temp']]
            dfz = npa[:,cdct['z']]
            clength = len(dft) - 1          # we will discard one value since we need delta t at each step
            nafrac = np.ones(clength)
        #    nadfrac = np.zeros(clength)
        #    nabfrag = np.zeros(clength)
            if sfrac is None:               # delayed since we didn't know the ion count before doing get_fracdat
                sfrac = np.zeros(clength)
                sdfrac = np.zeros(clength)
                sbfrag = np.zeros(clength)
                sTemp = np.zeros(clength)
                sz = np.zeros(clength)

#            tss = ''
#            for i in range(1,10):
#                tss = tss + f' {(dft[i]-dft[i-1])*1e6:0.2f}'
#            print(f'{tss}    {ion_max_t:0.2e}   {ts:0.2e} s    {ion_max_z:0.2f} m')

            lastzindx = dfz.argmax()                        # locate the index where ion reaches the end
            dt = dft[1:] - dft[:-1]
            nadfrac = lossrate(dfT[:-1], newdH, newdS)      # full array of loss rates (close to 0)
            nadfrac[lastzindx:] = 0.0                       # declare zero loss after ion reaches the end
            dfracremains = np.exp(-nadfrac*dt)              # with small losses, usually close to 1
            cdfracs = 1-np.exp(-nadfrac*dt)                 # full array of fractional losses
            rfrac = np.random.random(len(dfT)-1)
            nabfrag = rfrac < cdfracs                      # array with true values at times when fragmentation starts
            nabfrag[nabfrag.argmax():] = 1    # once fragmented, remains so
            lastfrac = dfracremains.prod()
            Nvals = len(dfracremains)
            for i in range(1, Nvals):
                nafrac[i] = nafrac[i-1] * dfracremains[i]

#            print(f'ion: {ionnum:3}   fetchtime: {fetchtime:4.2f}    calctime: {calctime:4.2f}   frac[-1]: {nafrac[-1]:6.4f}')
            slfrac = slfrac + lastfrac
            sfrac = sfrac + nafrac
            sdfrac = sdfrac + nadfrac
            sbfrag = sbfrag + nabfrag
            sTemp = sTemp + dfT[:-1]
            sz = sz + dfz[:-1]

            xlabel = gxlabel     #'t' or = 'z'
            if xlabel == 't':
                xdata = dft[:-1]
                xadata = xdata
            elif xlabel == 'z':
                xdata = dfz[:-1]
                xadata = sz / max_irid
            elif xlabel == 'Temp':
                xdata = dfT[:-1]
                xadata = sTemp / max_irid

            ylabel = gylabel    # or 'dfracs' or 'fragmenteds' or 'Temp'
            if ylabel == 'frac':
                ydata = nafrac
            elif ylabel == 'dfrac':
                ydata = nadfrac
            elif ylabel == 'bfrag':
                ydata = nabfrag
            elif ylabel == 'Temp':
                ydata = dfT[:-1]
            elif ylabel == 'z':
                ydata = dfz[:-1]

            title = f'Vin = {vin}'
            # try plotting dfview.Temp, fracs, fragmenteds, dfracs
            plt.plot(xdata, ydata, lw=0.5)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)

        if ylabel == 'frac':
            yadata = sfrac / max_irid
        elif ylabel == 'dfrac':
            yadata = sdfrac / max_irid
        elif ylabel == 'bfrag':
            yadata = sbfrag / max_irid
        elif ylabel == 'Temp':
            yadata = sTemp / max_irid
        elif ylabel == 'z':
            yadata = sz / max_irid

        cfracs[tcount] = sfrac[-1] / max_irid
        aveTemp = sTemp / max_irid
#        ctime = time() -starttime

        gfrac = fracremains(dft, aveTemp, newdH*1000, newdS, tmax)
        gfracs[tcount] = gfrac
        print(f'Vin: {vin}   cfrac: {slfrac/max_irid:0.3f}   gfrac: {gfrac:0.3f}    diff:{cfracs[tcount]-gfrac:0.3e}')
#        vintime = time() - starttime

        plt.plot(xadata, yadata, lw=2, label='Avg.')
        plt.legend()
        plt.show()

        tcount += 1
#        print(f'ctime: {(ctime):0.3g}   gtime: {(vintime-ctime):0.3g}   vintime: {vintime:0.3g}')

        # finally, we would like to see the frac remaining vs Vin based on average of individual fracs or based on averaged Temps.
        # 
        # now we've already computed average Temp vs time, so can call the spa_heat frag routine
        # times is time arrays for all Vins, Temps is temperature arrays for all Vins ddH, ddS are the specified dH,dS
        # tmax needed to define end time for each run (max for all Vins)
#        finalr = fracremains(times, Temps, ddH, ddS, maxtime=tmax)
        # this gives final frac remaining at all Vins
        # to get the comparable info from individual runs, we just gather the final frac value for each Vin above

        # Another item is to do averages in space rather than time.  But when ion trajectories have some backtracking, this
        # folds the temperature vs time plot. So this isn't well-defined.
    plt.plot(Vins, cfracs, label='average frac from trajs')
    plt.plot(Vins, gfracs, label='frac from average T')
    plt.title('compare frac vs Vin calculations')
    plt.xlabel('Vin')
    plt.ylabel('frac remaining')
    plt.legend()
    plt.show()
    return
    
""" Next steps:
        - Might rename program to some other name reflecting its use.
        - also get the position data from the db and plot variables vs. position in cell.
        - Do "standard" processing with averaged temperature vs. time and single run for fragmentation.
        - compare this to computing fragmentation for each ion and then averaging for fragment fractions.
        - and can compare averaging for temperature or fragmentation amount in space vs. in time.
    if False:
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

    return tcount
    """

def main():
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='python spa_db_refrag.py',
        description='''Given a .json file (as from spa_mpirun), and specified dH and dS, 
        recalculate the values for frac, dfrac, bfrag for each trajectory time for each 
        run for the range of Vin voltages.''')
    parser.add_argument('jsonfile', nargs='+', help='a json input file to specify ion, cell and other parameters.')
    parser.add_argument('--dH', help='overrides the initial dH value from input file') #   , nargs=1)
    parser.add_argument('--dS', help='overrides the initial dS value from input file') #   , nargs=1)
    parser.add_argument('--db', default = 'spa_db.sqlite3', help='use this for the database file')
    parser.add_argument('--xaxis', default = 't', help='x axis choice. May be t or z')
    parser.add_argument('--yaxis', default = 'frac', help='y axis choice may be frac, dfrac, bfrag, Temp, z')

    args = parser.parse_args()
    gdbname = path.splitext(args.db)[0] + '.sqlite3'
#    print(args)

    gxlabel = args.xaxis
    gylabel = args.yaxis

    paramd = loadjfile(args.jsonfile)
#    printparams(paramd)

    expd = paramd['exp']
    celld = paramd['cell']
    # might need this offset for update feature not yet implemented
    cellVoffset = celld.get('Voffset', 0.0)     # read the offset voltage to be added to Vins, default is 0.0
    iond = paramd['ion']
    newtemp = celld.get('T0offset', 0.0) + iond.get('T0')
    iond.update({'T0' : newtemp})

    dH = float( args.dH or gdH )
    dS = float( args.dS or gdS )

    cnt = doupdate(celld, iond, dH, dS)
    print(f'update count: {cnt}')
