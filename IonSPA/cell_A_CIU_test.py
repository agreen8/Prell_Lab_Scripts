# cell_A_CIU_test.py
# test program for Agilent in-source CIU cell
# define cell and ion and model and run
# plot results of run to verify operation of the model

#import numpy as np
#import pandas as pd
#from mpi4py import MPI
#import json, sys
#from os import path
import ionspa

from ionspa import const, ionclass, makecell, loadjfile, modelclass
print(f'in test script, ionspadir = {const.ionspadir}')
#import spa_heat as spah
#from time import sleep
#import sqlite3
#import argparse
import matplotlib.pyplot as plt

celld =  {    
    "type":"aCIUcell",
    "L":0.0200000, 
     "Mgamu":28, 
     "T":298, 
     "z0":-0.09000, 
     "CIUV":115
    }
iond = {
    	"name": "Stx13",
    	"mass": 39.111,
    	"charge": 13,
    	"CCS": 31.4,
    	"num_atoms": 5495,
    	"T0": 300,
        "hcprofname": "peptide",
        "refdH": 196,
        "refdS": 184
    }

ion = ionclass(iond)
print(ion)

gvz = []
vz = []
gT = []
grho = []
P = []
T = []
z = []
t = []
KE = []
E = []

for zstart in [-0.1805, -0.18, -0.17, -0.15, -0.13, -0.11, -0.09, -0.07, -0.05, -0.03, -0.01, -0.005, -0.001, 0, 0.001]:
    celld["z0"] = zstart
    cell = makecell(celld)
    print(cell)

    m = modelclass(cell, ion, 0)
    m.cell.setCIUV(111)
    print(f"cell vgz(0): {m.cell.vgz(0)}")
    z0 = m.cell.z0
    m.start(z=z0, dVin=0)
    m.vz = 0.95 * m.cell.vgz(m.cell.z0)

    print(repr(m))
    print('\n')

    print(m.headerstr())
    print(m)

    gvz.append(m.cell.vgz(z0))
    vz.append(m.vz)
    gT.append(m.cell.T(z0))
    grho.append(m.cell.rho(z0))
    P.append(m.cell.rho(m.z) * 1e-5 * m.cell.T(m.z) * const.kB / cell.Mg)
    T.append(m.temp)
    z.append(m.z)
    t.append(m.t)
    KE.append(m.KE/const.qe)
    E.append(m.cell.Ez(z0, 0)/1000)

    while m.z < m.L and m.collisioncount < 30000 and m.tcount < 100000000:
        m.steptime()

        if m.collided and (m.collisioncount % 100 == 1):
            #print(m)
            gvz.append(m.cell.vgz(m.z))
            vz.append(m.vz)
            gT.append(m.cell.T(m.z))
            grho.append(m.cell.rho(m.z))
            P.append(m.cell.rho(m.z) * 1e-5 * m.cell.T(m.z) * const.kB / cell.Mg)
            T.append(m.temp)
            z.append(m.z)
            t.append(m.t)
            KE.append(m.KE/const.qe)
            E.append(m.cell.Ez(m.z, m.t)/1000)

            if m.collisioncount % 1000 == 1:
                zz = m.z
                print(f'{1e6*m.t:9.5f}  {1000*zz:7.5f}   {P[-1]:5.3f}  {m.collisioncount:7d}   {m.cell.vgz(zz):6.2f}  {m.vz:6.2f}    {m.cell.T(zz):6.2f}  {m.temp:6.2f}    {0.001*m.cell.Ez(zz,0):6.3f}   {m.KE/const.qe:5.2f}')

#print(f'z min {min(z)}    z max {max(z)}')

plt.figure()
plt.plot(z, T, label='T')
plt.plot(z, gT, label='gT')
plt.legend()
plt.title("T vs z")
plt.show()

plt.figure()
plt.plot(z, gvz, label='gvz')
plt.plot(z, vz, label='vz')
plt.plot(z, gT, label='gT')
plt.plot(z, T, label='T')
plt.plot(z, KE, label='KE')
plt.plot(z, E, label='E')
plt.xlabel('z')
plt.legend()
plt.show()
