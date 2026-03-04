####################################
#
# spa_mpirun.py
#
# Original work by Ken Newton, et.at. at University of Oregon
# Copyright (c) 2024 James S. Prell, all rights reserved
#
# MPI script to run ionspa for given ion with given cell over a range of
# input energies.
# (developed after test_tune_MPI_B.py)
#
# Use spa_heat.py model for ion heating, running a collection of ions
# Ion time-temperature arrays are interpolated to fixed time points, averaged
# and saved to a database file.
#
# The ion, cell, and experiment data files can be specified via a JSON
# formatted input file.
#
####################################


import numpy as np
import pandas as pd
from mpi4py import MPI
import json, sys
from os import path
from ionspa import const, ionclass, makecell, loadjfile, modelclass
#import spa_heat as spah
#from spa_heat import const, ionclass, makecell
#from spa_paramd import loadjfile
from time import sleep
import sqlite3
import argparse

'''
MPI ion_spa xeff tune strategy:
- The Rank0 host will direct the full run and do all work other than the actual spa_heat model runs.
- The workers will only be used to do the time-consuming spa_heat model runs.
  - The workers need to know the cell parameters, ion parameters, and model parameters.
  - The cell, ion and general model params are read from input files and shared by all.
  - The Vin model parameters may be different for some runs so must be passed from the host.
  - So host needs to send a message with Vin. When worker receives this, it starts a run. (timestep and tmax also sent per run)
  - The worker runs the run_model() function, returing the tT array to the host.
  - To simplify the program logic, all workers at any one time will run with the same Vin, so will just
    be doing the "repeat" runs.  Anticipated number of workers is in the range of 2 to 100 (likely 1 less than cpu core count).  So the
    value for Nrepeat can be set to an integer multiple of the number of workers. (If not, only a little efficiency is lost.)
  - Since all workers get the same parameters at any one time, this can be done with a broadcast kind of message send.
  - A further simplification would be to have all workers return temperatures at the identical sequence of times.
    - This might be done by having the workers do the job of interpolating to a fixed time sequence, starting
      at time 0, with a specified end time and step time. (Add these params to the host-worker message.)
    - This adds a requirement to get a reasonable estimate of total time for each condition -- mainly transit time vs Vin.
      - Or we can just overestimate, but also have each worker-run return the actual max time before extrapolation.
        (more complex return structure)
      For all runs where the ion completes in less than the given end time, the worker would just insert constant temperature
      values (equal to the final temperature). This is the "normal" behavior for the np.interp() function.
  - Dividing the summed array by the number of workers yields an average.  Average arrays from multiple repetitions of the
    model runs gives the full average temperature.
  - Once the host has the time vs temperature array averaged from all the desired Nrepeats and also loops over all Vin values,
    it can do the remaining work within the single host process.  For this program, that includes cycling through the list of
    Vin values and storing into a combined vtTd array. This array is written to the database file.
  - Later the spa_db_fit program is used to read the db file and to run the optimization to match the experimental data
    until an optimized dH, dS pair is found.
'''


# MPI definitions
def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def db_insert_celld(dbname, celld):
    '''Return a cell id for the cell specified by the cell dict celld.
    If db is new and table not present, create cells table.
    If the cell is already present, find it and return the cellid.
    If the cell is not yet present, insert a new cell item from celld and return the new cellid.'''
    cellid = 0
    with sqlite3.connect(dbname) as db:
        db.execute('''CREATE TABLE IF NOT EXISTS "cells" ("cid" integer primary key,"type" text,"L" float,"T" float,"pressure1" float, "rho" float,"Mgamu" float,"z0" float,"otherjson" text);''')

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
        otherjson = json.dumps(od, separators=(',', ':'), sort_keys=True)         # compact form with sorted keys
        tnames = 'cid,type,L,T,pressure1,rho,Mgamu,z0,otherjson'
        query = f"SELECT {tnames} from cells where type=\'{celld['type']}\' and T={celld['T']};"
#        print(query)
        cur = db.execute(query)  # query returns a cursor object
        fdat = cur.fetchall()    # which can fetch the data as a list of items, each item being a tuple
        
        for itemlist in fdat:
            itemd = dict(zip(tnames.split(','), itemlist))
            # each itemd matches cell type and Temperature, but also needs to match other fields
            # so check for match for L, rho, Mgamu, z0, otherjson
            try:
                if itemd['L'] == celld['L'] and itemd['rho'] == celld['rho'] and itemd['Mgamu'] == celld['Mgamu'] and itemd['z0'] == celld['z0'] and itemd['otherjson'] == otherjson:
#                    print(f'found {itemd}')
                    cellid = itemd['cid']
            except:
                if itemd['L'] == celld['L'] and itemd['pressure1'] == celld['pressure1'] and itemd['Mgamu'] == celld['Mgamu'] and itemd['z0'] == celld['z0'] and itemd['otherjson'] == otherjson:
#                    print(f'found {itemd}')
                    cellid = itemd['cid']

        if cellid != 0:
            return cellid
        else:
            # add record to database
            try:
                db.execute('''INSERT INTO cells(type, L, T, rho, Mgamu, z0, otherjson) VALUES(?,?,?,?,?,?,?);''',
                           [ celld['type'], celld['L'], celld['T'], celld['rho'], celld['Mgamu'], celld['z0'], otherjson ] )
            except:
                db.execute('''INSERT INTO cells(type, L, T, pressure1, Mgamu, z0, otherjson) VALUES(?,?,?,?,?,?,?);''',
                           [ celld['type'], celld['L'], celld['T'], celld['pressure1'], celld['Mgamu'], celld['z0'], otherjson ] )
            cur = db.execute('select last_insert_rowid();')
            ids = cur.fetchone()        # returns (id, )
            cellid = ids[0]
            # print(f'new cid is {cellid}')
    return cellid


def db_insert_iond(dbname, iond):
    '''Return a ion id for the ion specified by the ion dict iond.
    If db is new and table not present, create ions table.
    If the ion is already present, find it and return the ionid.
    If the ion is not yet present, insert a new ion item from iond and return the new ionid.'''
    ionid = 0
    with sqlite3.connect(dbname) as db:
        db.execute('CREATE TABLE IF NOT EXISTS "ions" ("iid" integer primary key, "name" text, "mass" float, "charge" float, "CCS" float, "num_atoms" integer, "T0" float, "hcprofname" text, "dH" float, "dS" float);')

        tnames = 'iid,name,mass,charge,CCS,num_atoms,T0,hcprofname'
        query = f"SELECT {tnames} from ions where name=\'{iond['name']}\' and hcprofname=\'{iond['hcprofname']}\';"
#        print(query)
        cur = db.execute(query)  # query returns a cursor object
        fdat = cur.fetchall()    # which can fetch the data as a list of items, each item being a tuple
        for itemlist in fdat:
            itemd = dict(zip(tnames.split(','), itemlist))
            # each itemd matches ion name and hcprofname, but also needs to match other fields
            # so check for match for mass, charge, CCS, num_atoms, T0
            if itemd['mass'] == iond['mass'] and itemd['charge'] == iond['charge'] and itemd['CCS'] == iond['CCS'] and itemd['num_atoms'] == iond['num_atoms'] and itemd['T0'] == iond['T0']:
#                print(f'found {itemd}')
                ionid = itemd['iid']

        if ionid != 0:
            return ionid
        else:
            # only insert new ion into table if it is now yet present in table.
            db.execute('''INSERT INTO ions(name,mass,charge,CCS,num_atoms,T0,hcprofname,dH,dS)
                VALUES(?,?,?,?,?,?,?,?,?);''',
                [ iond['name'], iond['mass'], iond['charge'], iond['CCS'],
                iond['num_atoms'], iond['T0'], iond['hcprofname'],
                iond.get('dH', 80), iond.get('dS', -120) ])
            cur = db.execute('select last_insert_rowid();')
            ids = cur.fetchone()        # returns (id, )
            ionid = ids[0]
#            print(f'new iid is {ionid}')
    return ionid


def db_insert_VinRuns(gdbname, cellid, ionid, Vins, Nrepeats):
    '''Return the highest run id for the run specified by cellid, ionid and Vins list.
    If db is new and table not present, create runs table.
    If the runs are present, find it and return the (runid,tmax, newNrepeats).
    If the ion is not yet present, insert a new ion item from iond and return the new (runid, 0.0, Nrepeats).
    newNrepeats is the number of additional runs needed to get to a total of Nrepeats runs for each rid
    and if it is <= 0, the runs will not be done.'''
    runids = [ [0, 0.0] for Vin in Vins]
    runidsd = dict(zip(Vins, runids))
    with sqlite3.connect(gdbname) as db:
        db.execute('CREATE TABLE IF NOT EXISTS "runs" ("rid" integer primary key, "cid" integer, "iid" integer, "Vin" integer, "tmax" float)')
        # and create the trajectories table if needed
        db.execute('CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float, "collnum" float, "frac" float, "dfrac" float, "bfrag" float, "ma" float)')

        tnames = 'rid,Vin,mirid,tmax'
        query = f'select runs.rid,Vin,max(irid) as mirid, tmax from runs join trajs on runs.rid=trajs.rid where cid={cellid} and iid={ionid} group by Vin'
        cur = db.execute(query)
        fdat = cur.fetchall()
        runs = {}
        for itemlist in fdat:
            rid,Vin,mirid,tmax = itemlist
            # itemd = dict(zip(tnames.split(','), itemlist))   # itemd is like dict(rid, cid, iid, Vin, tmax)
            runs[Vin] = [rid, tmax, Nrepeats-mirid]
        # print(f'runs from last round is {runs}')
        # now ddat has list of Vins already in db with given cid, iid
        for thisVin in runidsd:
            if thisVin in runs:
                # print(f'already present {runidsd[thisVin]} -> {runs[thisVin]}')
                runidsd[thisVin] = runs[thisVin]
            else:
                db.execute('insert into runs(cid,iid,Vin,tmax) Values(?,?,?,?);', [cellid,ionid,thisVin,0.0])
                runid = db.execute('select last_insert_rowid();').fetchone()[0]
                runidsd[thisVin] = [runid, 0.0, Nrepeats]
    return runidsd


def db_insert_rundata(gdbname, trajs, rid, rtmax):
    '''Insert trajectory information into trajectory table in db.'''
    with sqlite3.connect(gdbname) as db:
        query = f'select max(irid) from trajs where rid={rid}'
        ret = db.execute(query)
        irid = ret.fetchone()[0] or 0

        for traj in trajs:  # multiple runs for this rid combination
            t,T,u,z,KE,collnum,frac,dfrac,bfrag,ma = traj
            nvals = len(t)
            irid += 1   # note next ion at this Vin, rid combination
            query = f'update runs set tmax={rtmax} where rid={rid}'
    #        print(query)
            ret = db.execute(query)
            # print(ret.rowcount)   # I think we expect rowcount to be 1 here
            for i in range(nvals):
                if t[i] <= rtmax:
                    ret = db.execute('insert into trajs(rid,irid,t,Temp,u,z,KE,collnum,frac,dfrac,bfrag,ma) Values(?,?,?,?,?,?,?,?,?,?,?,?);',
                                     [rid, irid, t[i], T[i], u[i], z[i], KE[i], collnum[i], frac[i], dfrac[i], float(bfrag[i]), ma[i]] )
                    # print('  ', ret)


def controller(celld, cellid, iond, ionid, runids):
    '''Main program for rank 0.
     Control the sequence of the model runs, sending work to distributed
     workers as needed.
    '''
    import matplotlib.pyplot as plt

    # wait for all workers to report as READY
    num_workers = size - 1
    num_ready = 0
    workers = dict()
    for wr in range(1, comm.size):
        workers[wr] = dict(state=tags.NONE)
    # print(workers, flush=True)
    while num_ready < num_workers:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        # print('recvd tag {} from worker {}'.format(tag, source), flush=True)
        workers[source]['state'] = tag
        if tag == tags.READY:
            num_ready += 1
    # print(workers, num_ready, flush=True)

    Vins = list(runids.keys())

    Nrepeats = gNrepeats
    # tstep = 6e-8
    tstep = gtstep  # 1.8e-7
    tmax = gtmax
    tt = np.linspace(0.0, tmax, int(tmax/tstep)+1)
    # tt = np.linspace(0.0, tmax, 1001)
# 
    vtTd = {}
    for Vin in Vins:
        rid, dtmax, Nrepeats = runids[Vin]    # the runid and old tmax from the database (or tmax=0.0)
        if Nrepeats <= 0:
            continue
        rt_tmax, trajs = distribute_model(workers, Nrepeats, Vin, tstep, tmax)
        rtmax = max(dtmax, rt_tmax)
        runids[Vin][1] = rtmax     # update tmax for this Vin value
        db_insert_rundata(gdbname, trajs, rid, rtmax)   # update rtmax in runs table, add rows to trajs table
        vtTd[Vin] = [tt, trajs, rt_tmax]
        if gplot_tT:
            for traj in trajs:
                plt.plot(traj[0], traj[1])
            plt.show()

    print('Done with work, ask all workers to EXIT\n', flush=True)
    for wr, worker in workers.items():
        # print('sending EXIT to worker {}'.format(wr))
        comm.send(None, dest=wr, tag=tags.EXIT)

    print('sleeping for 2 sec')
    sleep(2)


def distribute_model(workers, Nrepeats, Vin, tstep, tmax):
    # The remainder is done for each Vin
    # returning the associated average(maxrts) and temperature array.

    nTasks = Nrepeats
    tt = np.linspace(0.0, tmax, int(tmax/tstep)+1)
    aveT = np.zeros_like(tt)
    tsk = dict(Vin=Vin, tstep=tstep, tmax=tmax)
    tasks = [tsk]*Nrepeats

    # send initial tasks for this group
    for wr, worker in workers.items():
        try:
            tsk = tasks.pop()
            comm.send(tsk, dest=wr, tag=tags.START)
            # print("Sending task '{}' to worker {}".format(tsk, wr), flush=True)
        except: # no more tasks
            pass

    # now, we can wait for responses and send more tasks till done
    NdatRcvd = 0
    num_workers = size - 1
    closed_workers = 0
    maxrts = []
    alltrajs = []
    # print("Master starting with {} workers".format(num_workers), flush=True)
    while (closed_workers < num_workers) and (NdatRcvd < Nrepeats) :         # if num_workers is zero, don't run
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.DONE:
            rtmax, traj = data
            alltrajs.append(traj)   # full list of all trajectories from all runs with these conditions
            maxrts.append(rtmax)
            # print(maxrts)
            # tts = traj[0]     # elements of the traj array are set by what the worker sets.  Examine that code below for details.
            Temps = traj[1]
            # Uvals = traj[2]
            # zvals = traj[3]
            # kevals = traj[4]
            # collnums = traj[5]
            # fracs = traj[6]
            # dfracs = traj[7]
            # bfrags = traj[8]
            # mas = traj[9]
            # maxT = Temps.max()
            # lastT = Temps[-1]
            # print("xeff {:0.2f}, Vin {:0.1f}, tave {:0.2f}, maxT {:0.1f}, lastT {:0.1f} from worker {}".format(xeff, Vin, rtmax*1000, maxT, lastT, source), flush=True)
            # aveT = np.zeros_like(traj[0])
            aveT = aveT + Temps / float(nTasks)
            NdatRcvd += 1
            try:
                tsk = tasks.pop()
                comm.send(tsk, dest=source, tag=tags.START)
                # print("Sending task '{}' to worker {}".format(tsk, source), flush=True)
            except: # no more tasks
                pass
        elif tag == tags.EXIT:
            # print("Worker {} exited.".format(source), flush=True)
            closed_workers += 1

    # print("Master finishing", flush=True)
    maxrts = np.array(maxrts)
    rt_tmax = np.max(maxrts)
    
    Tmax = max(aveT)
    print(f'{Nrepeats} runs complete for Vin={Vin} with rt_tmax={rt_tmax*1000:0.2f} ms and Tmax= {Tmax:0.1f}', flush=True)
    return rt_tmax, alltrajs


def worker():
    '''Do the requested computational work as requested by the controller.
    '''

    tmpVin = 5        # temporary values to be replaced later
    mmodel = modelclass(gcell, gion, tmpVin, ma_kg_mol=gpsMkg)
    if rank == 1:
        print(repr(mmodel), flush=True)

    def run_model(Vin, tt, tmax):
        '''With given Vin input energy, fly ion through cell, returning complete history with t, T, u, z, KE.
        This is intended to be used to replace workrun() below and doesn't require modifications to spa_heat_impulsive.
        '''
        if mmodel.cell.type == 'A_CIU':
            mmodel.cell.setCIUV(Vin)
            z0 = mmodel.cell.z0
            mmodel.start(z=z0, dVin=0)
        else:
            mmodel.start(z=0, dVin=Vin)

        t = [mmodel.t]        # put starting values into lists
        T = [mmodel.temp]
        u = [mmodel.u/const.qe]
        position = [mmodel.z]
        KE = [Vin]
        collnum = [mmodel.collisioncount]
        frac = [mmodel.frac]
        dfrac = [mmodel.dfrac]
        bfrag = [mmodel.fragmented]
        ma = [mmodel.ma]

        time_cut = []
        Temp_cut = []
        KE_cut = []
        while mmodel.z <= mmodel.cell.L and mmodel.tcount < gcountmax and mmodel.t <= tmax:
            mmodel.steptime()
            
            if (mmodel.KE) / (mmodel.ion.Z) < 0.5 and abs(mmodel.temp - mmodel.Tg(mmodel.z))/mmodel.Tg(mmodel.z) < 0.025:
                time_cut.append(mmodel.t)
                Temp_cut.append(mmodel.temp)
                KE_cut.append(mmodel.KE/const.qe)
                
                if (time_cut[-1] - time_cut[0]) > 0.0005:
                    #print('broken here')
                    #print(time_cut[-1])
                    break
            else:
                time_cut = []
                Temp_cut = []
                KE_cut = []
            
            if mmodel.collided or mmodel.tcount < 4:      # only add to array if ion just collided (or first few points)
                t.append(mmodel.t)
                T.append(mmodel.temp)
                u.append(mmodel.u/const.qe)
                position.append(mmodel.z)
                KE.append(mmodel.KE / const.qe)
                collnum.append(mmodel.collisioncount)
                frac.append(mmodel.frac)
                dfrac.append(mmodel.dfrac)
                bfrag.append(mmodel.fragmented)
                ma.append(mmodel.ma)

        trajinfo = np.array([t, T, u, position, KE, collnum, frac, dfrac, bfrag, ma])   # sets contents of the trajinfo array
        rtmax = max(t)      # this is time just before z exceeds cell.L
        # would we want to return binned data or raw?  Currently raw data is returned.
        # np.interp can be used interpolate each column to a linear timeline one column at a time
        # alternately, we might wish to do this interpolation to a linear position array.
        # most of the behavior is fairly smooth, so around 200 time points and 200 position points
        # might have enough resolution to characterize the data.  (Agilent cell is 180 mm long)
        # Note that positions don't necessarily increase monotonically in time (ions can travel backwards)
        # tt = np.linspace(0.0, rtmax, 1001)
        lentt = len(tt)
        traj = np.resize(tt, (10,lentt))                # make an array to hold data.  Following lines define the contents
        traj[1] = np.interp(tt, trajinfo[0], trajinfo[1])   # T
        traj[2] = np.interp(tt, trajinfo[0], trajinfo[2])   # u
        traj[3] = np.interp(tt, trajinfo[0], trajinfo[3])   # position
        traj[4] = np.interp(tt, trajinfo[0], trajinfo[4])   # KE
        traj[5] = np.interp(tt, trajinfo[0], trajinfo[5])   # collisioncount
        traj[6] = np.interp(tt, trajinfo[0], trajinfo[6])   # frac
        traj[7] = np.interp(tt, trajinfo[0], trajinfo[7])   # dfrac
        traj[8] = np.interp(tt, trajinfo[0], trajinfo[8])   # bfrag
        traj[9] = np.interp(tt, trajinfo[0], trajinfo[9])   # ma)
        return(rtmax, traj)

    # Worker processes execute code below
    # name = MPI.Get_processor_name()
    # print("I am a worker with rank {} on {}.".format(rank, name), flush=True)
    while True:
        comm.send(None, dest=0, tag=tags.READY)
        task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()

        if tag == tags.START:
            # hope to get a task dict with run information
#            taskd = dict(Vin=5, xeff=0.9, tstep=1e-6, tmax=1e-3)
            taskd = task
            Vin = taskd['Vin']
            tstep = taskd['tstep']
            tmax = taskd['tmax']
            tt = np.linspace(0.0, tmax, int(tmax/tstep)+1)   
            # tt = np.linspace(0.0, tmax, 1001)

            # Do the work here
            rtmax, tT = run_model(Vin, tt, tmax)
            result = (rtmax, tT)
            # result = "task '{}' in rank({})".format(task, rank)
            comm.send(result, dest=0, tag=tags.DONE)
        elif tag == tags.EXIT:
            break
    comm.send(None, dest=0, tag=tags.EXIT)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---------------------- Common code for all ranks -----------------------------

# ---------------- Common model setup -----------------

#########################################
#
# First define ions and store into dict ions
# Then define the Waters cell object wcell
# Then read in exp data files and store into dict expdats
#
#########################################

# Define ions. Once you know the correct AA sequence, you can get number of atoms from https://web.expasy.org/protparam/ easily.
gion = None
gcell = None
gexp = None


# global settings for the run
gplot_vf = False
gplot_tT = False
gNrepeats = 28  # using multiple of 7 bc we have 7 worker threads on lab computer
gtstep = 1.8e-7
gtmax = gtstep * 7500

# ---------------- Common MPI setup -----------------

# Define MPI message tags
# We'll need distinct tags for each different message sent between Rank0 host and RankN workers.
tags = enum('NONE', 'READY', 'DONE', 'EXIT', 'START')

# Initializations and preliminaries
# Both the Rank0 host and workers will get their own copy of these variables/objects
comm = MPI.COMM_WORLD   # get MPI communicator object
size = comm.size        # total number of processes
rank = comm.rank        # rank of this process
status = MPI.Status()   # get MPI status object

# defaults for any parameter not supplied
paramd = dict(cell=None, ion=None, exp=None)
iond = dict(name="CC322", mass=0.322, charge=1, CCS=1.5367, num_atoms=37, T0=298.0, hcprofname="tunemix322")

parser = argparse.ArgumentParser(
    prog='python spa_mpirun.py',
    description='''Given a .json file and an optional db name.''')
parser.add_argument('jsonfile', nargs='+')
parser.add_argument('--db', default="spa_db.sqlite3") #   , nargs=1)
parser.add_argument('--psM', default=64)
args = parser.parse_args()
gdbname = path.splitext(args.db)[0] + '.sqlite3'
gpsMkg = float(args.psM) * 0.001

gbasename = None
paramd = loadjfile(args.jsonfile)


gcell = makecell(paramd.get("cell"))
ion = paramd.get('ion')
newtemp = paramd.get("cell").get('T0offset', 0.0) + ion.get('T0')
ion.update({'T0' : newtemp})
gion = ionclass(ion)
gplot_vf = paramd.get("plotvf") or False
gplot_tT = paramd.get("plottT") or False
gNrepeats = paramd.get("Nrepeats") or 28
g_dhds_start = paramd.get("dhds_start") or [100, -100]
gtstep = paramd.get("tstep") or gtstep
gtmax = paramd.get("tmax") or gtmax
gcountmax = paramd.get("maxsteps") or 1e7

if size <= 1:
    print()
    print('Run this with the command mpiexec python spa_mpirun.py exp_fit_xxx.json')
    print(' - optional arguments are --db spa_db_test --psM 68  to specify db file and pseudoatom mass.')
    print('This version is very verbose with diagnostic text.')
    print('It calculates nTasks = 140 runs of the WAVE model and averages the T vs t data.')
    exit(1)
if rank == 0:
    print(args)
    print('running from input files', args.jsonfile, 'with gplot_vf', gplot_vf, 'and gplot_tT', gplot_tT, 'and dhds', g_dhds_start)
    print(paramd, '\n')
    print(f'pseudo-atom mass in kg/mol is {gpsMkg:0.4f}')
    print(f'Will use sqlite3 DB {gdbname} to store results')
    print('gion', gion)
    print('gcell', gcell)
    print()

    celld = paramd.get("cell")
    cellVoffset = celld.get('Voffset', 0.0)     # read the offset voltage to be added to Vins, default is 0.0
    cellid = db_insert_celld(gdbname, celld)
    print(f'cellid is {cellid}')

    iond = paramd.get('ion')
    ionid = db_insert_iond(gdbname, iond)
    print(f'ionid is {ionid}')

    expd = paramd.get("exp") or dict(folder="expdata", file="CC322.csv")
    Vins = expd.get('Vins', None)       # if "Vins" defined, use that list, otherwise lookup exp data Vins
    if not Vins:
        gexp = pd.read_table(path.join(expd['folder'], expd['file']), delimiter=',', header=None, names=['Vin', 'efrac'])
        Vins = gexp.Vin 

    Vins = [vin+cellVoffset for vin in Vins]
    
    runids = db_insert_VinRuns(gdbname, cellid, ionid, Vins, gNrepeats)
    
    
    
    print(f'runids is {runids}', flush=True)
    controller(celld, cellid, iond, ionid, runids)
else:
    worker()