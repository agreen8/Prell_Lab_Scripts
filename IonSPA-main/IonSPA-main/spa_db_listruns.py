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


#import numpy as np
import pandas as pd
from scipy.optimize import minimize
#from mpi4py import MPI
import sys
#from os import path
#import spa_heat_impulsive2_QMhc as spah
#from spa_heat_impulsive2_QMhc import const, ionclass, makecell
#from time import sleep
import sqlite3
#import matplotlib.pyplot as plt
#plt.interactive(False)

'''
This program rusn without the mpi system so will normally be called as:
python spa_db_listruns.py {database_name}

DB table creation:
CREATE TABLE IF NOT EXISTS "cells" ("cid" integer primary key,"type" text,"L" float,"T" float,"rho" float,"Mgamu" float,"z0" float,"otherjson" text);
CREATE TABLE IF NOT EXISTS "ions" ("iid" integer primary key, "name" text, "mass" float, "charge" float, "CCS" float, "num_atoms" integer, "T0" float, "hcprofname" text);
CREATE TABLE IF NOT EXISTS "runs" ("rid" integer primary key, "cid" integer, "iid" integer, "Vin" integer, "tmax" float);
#CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float);
CREATE TABLE IF NOT EXISTS "trajs" ("rid" integer, "irid" integer, "t" float, "Temp" float, "u" float, "z" float, "KE" float, "collnum" float, "frac" float, "dfrac" float, "bfrag" float, "ma" float);

'''

gdbname = 'spa_db.sqlite3'

if __name__ == '__main__':
    args = sys.argv
    if len(args) > 1:
        dbname = args[1]
    else:
        dbname = gdbname

    with sqlite3.connect(dbname) as db:
        cnames = 'cid type L T pressure1 rho Mgamu z0 otherjson'.split()
        csel = ', '.join(cnames)
        query = f'SELECT {csel} from cells'
        df = pd.read_sql_query(query, db)
        print ('cells table:')
        print(df)
        print()

        inames = 'iid name mass charge CCS num_atoms T0 hcprofname dH dS'.split()
        isel = ', '.join(inames)
        query = f'SELECT {isel} FROM ions'
        df = pd.read_sql_query(query, db)
        print('ions table:')
        print(df)
        print()

        query = 'select rid, Vin, cells.cid, cells.type, ions.iid, ions.name, tmax from runs join cells on cells.cid=runs.cid join ions on ions.iid=runs.iid'
        df = pd.read_sql_query(query, db)
        print('runs table:')
        print(df.to_string())
        print()

        query = 'select count() as recs from trajs'
        df = pd.read_sql_query(query, db)
        print(f"total ion trajectory points stored: {(df['recs'].values)[0]}")
#        print((df['recs'].values)[0])
