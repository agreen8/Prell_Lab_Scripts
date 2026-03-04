# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

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



gdbname = 'synapt_right_bsa_3V_databases/BSA_17_1.sqlite3'

dbname = gdbname

with sqlite3.connect(dbname) as db:
    cnames = 'cid type L T pressure1 rho Mgamu z0 otherjson'.split()
    csel = ', '.join(cnames)
    query = f'SELECT {csel} from cells'
    cell_df = pd.read_sql_query(query, db)
    print ('cells table:')
    print(cell_df)
    print()

    inames = 'iid name mass charge CCS num_atoms T0 hcprofname dH dS'.split()
    isel = ', '.join(inames)
    query = f'SELECT {isel} FROM ions'
    ion_df = pd.read_sql_query(query, db)
    print('ions table:')
    print(ion_df)
    print()

    query = 'select rid, Vin, cells.cid, cells.type, ions.iid, ions.name, tmax from runs join cells on cells.cid=runs.cid join ions on ions.iid=runs.iid'
    df = pd.read_sql_query(query, db)
    print('runs table:')
    print(df.to_string())
    print()


    cid = cell_df['cid'][0]
    iid = ion_df['iid'][0]
    Vin = 83.0
    tnames = 'rid,cid,iid,Vin,tmax'
    query = f'SELECT {tnames} from runs where cid={cid} and iid={iid} and tmax!=0.0'
# with sqlite3.connect(dbname) as db:
    cur = db.execute(query)
    fdat = cur.fetchall()

    # for itemlist in fdat:
    itemd = dict(zip(tnames.split(','),fdat))

    rids = [v[0] for v in fdat]
    Vins = [v[3] for v in fdat]

    query = f'SELECT rid, irid, t, z, Temp FROM trajs'#' WHERE rid={rids[0]}'
    trajs_df = pd.read_sql_query(query, db)

db.close()
db = None
for rid in rids:
    
    sub_df = trajs_df[trajs_df['rid'] == rid]
    
    ts = sub_df['t']
    Temp = sub_df['Temp']
    plt.plot(ts, Temp)

























































































