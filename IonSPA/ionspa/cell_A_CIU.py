# Copyright (c) 2024 James S. Prell, all rights reserved

import numpy as np
import pandas as pd
from collections import defaultdict
from .cell import cellclass
from .spa_const import const
from os import path
import matplotlib.pyplot as plt

'''
Working notes:

This is a class to represent the Agilent CIU in-source activation device using data modeled with a CFD program and an electric field program.
Output data from these programs will be stored in files:
- CIUb_vTrho.txt  z, vz, T, rho columns
- CIUb_zVsc.txt   z, Vbase, dV/V columns
    Voltage for any applied CIU voltage is Vbase(z) + Vciu * dV/V(z)
    Electric field calculated from voltages.
    Length of cell can be specified less than max z in files.
    Ion flight should start at z=0.

This model is still in development and testing, specifically with regard to the modeled physical parameters.
'''


def loadfile(fname='CIUb_zVsc.txt'):
    '''Loads CIU model file with voltage as function of z position.
    The first voltage column is the base voltage, the voltage with
    a zero applied CIU potential. The second column is the voltage
    scale factor to be multiplied by the applied CIU potential.
    Returns a dataframe with columns z, Vbase, Vscale.
    Alternately can load the vTrho file. In this case the columns
    of the returned dataframe are z, v, T, rho.
    '''
    pname = path.join(const.ionspadir, fname)
    with open(pname, 'r') as f:
        line = f.readline()
    vtlist = line.split()
    #print(vtlist)
    # vlist = [float(vt) for vt in vtlist]

    vdat = np.loadtxt(pname, skiprows=1)
    vf = pd.DataFrame(vdat, columns=vtlist)
    return vf


class CIUcell(cellclass):
    '''Derived cell class for the Agilent CIU in-source activation device.'''


    def __init__(self, **args):
        pdict = defaultdict(lambda:None)    # return None for any key not defined in the dict
        pdict.update(args)
        super(CIUcell, self).__init__(L=pdict['L'], Mgamu=pdict['Mgamu'], T=pdict['T'], rho=pdict['rho'], z0=pdict['z0'])
        self.type = 'A_CIU'
        # ffile=pdict['ffile']
        self.CIUV = pdict.get('CIUV', 0)
        dfV = loadfile(r'CIUb_zVsc.txt')
        self.vbase = dfV.base
        self.dV_V = dfV.dV_V
        # dfC = loadfile(r'CIUb_vTrho.txt')
        dfC = loadfile(r'CIU_vTrho_new_6560c_1mmleft.txt')
        # capEndZ = 0.8124495909
        capEndZ = 0.0
        self._zElec = dfV.zz
        self._zflow = dfC.z - capEndZ
        #print(f'in cellA: maxzE:{self._zElec.max()}    maxzflo:{self._zflow.max()}')
        self.L = min(self.L, self._zElec.max(), self._zflow.max())
        self._Tz = dfC.temperature
        self._rhoz = dfC.rho
        self._vz = dfC.velocity
        self.setCIUV(self.CIUV)
#        self._Vz = dfV.base + dfV.dV_V * self.CIUV
        self.Ez = self._Ez
        self.T = self._T
        self.rho = self._rho
        self.vgz = self._vgz

    def setCIUV(self, ciuV):
        self.CIUV = ciuV
        self._Vz = self.vbase + self.dV_V * self.CIUV

    def _T(self, z):
        '''Return temperature at z location, interpolating from file data.'''
        T = np.interp(z, self._zflow, self._Tz)
        return T

    def _rho(self, z):
        '''Return temperature at z location, interpolating from file data.'''
        rho = np.interp(z, self._zflow, self._rhoz)
        return rho

    def _vgz(self, z):
        vgz = np.interp(z, self._zflow, self._vz)
        return vgz

    def _Ez(self, z, t):
        '''Return Ez at z location, interpolating from file data.'''
        dz = 0.2e-3     # 0.2 mm
        V0 = np.interp(z, self._zElec, self._Vz)
        V1 = np.interp(z+dz, self._zElec, self._Vz)
        Ez = (V0-V1)/dz
        return Ez

    def __str__(self):
        outs = super(CIUcell, self).__str__() + '\n'
        outs += 'CIU: {}'.format('CIUfile')
        return outs


cell = None


def testCIU(Tion=298, CIUV=10):
    expl = '''Test run with CIU cell.  The initial ion temperature is 298 K.
    The internal ion temperature would be expected to also settle close to
    that same temperature.
    The ion is myoglobin.
    '''
    from spa_heat import ionclass, modelclass

    global cell
    print(expl)
    cell = CIUcell(Mgamu=28, z0=0, CIUV=115)
    #iond = dict(name='myoglobin', mass=17.585, charge=9, CCS=20.2, num_atoms=2411, T0=Tion,polarizability=0)
    iond = dict(name='IS322', mass=0.322, charge=1, CCS=1.5367, num_atoms=37, T0=298, polarizability=0, hcprofname='tunemix322')
    #iond = dict(name='IS622', mass=0.622, charge=1, CCS=2.0267, num_atoms=55, T0=298, polarizability=0)
    #iond = dict(name='IS922', mass=0.922, charge=1, CCS=2.4305, num_atoms=73, T0=298, polarizability=0)
    #iond = dict(name='IS1222', mass=1.222, charge=1, CCS=2.8125, num_atoms=91, T0=298, polarizability=0)

    ion = ionclass(iond)
    m = modelclass(cell, ion, 0)
    print(repr(m))
    mm, tT = m.run_model(0, plotresult=True, verbose=3)

def plotvals():
    import matplotlib.pyplot as plt
    celld0 =  {    
    "type":"aCIUcell",
    "L":-0.0100000, 
     "Mgamu":28, 
     "T":298, 
     "z0":-0.010800, 
     "CIUV":0
    }
    celld1 =  {    
    "type":"aCIUcell",
    "L":-0.0100000, 
     "Mgamu":28, 
     "T":298, 
     "z0":-0.010800, 
     "CIUV":100
    }
    cell0 = makecell(celld0)
    cell1 = makecell(celld1)

    zz = [-0.180 + i*0.001 for i in range(200)]
    T0 = [cell0.T(z) for z in zz]
    T1 = [cell1.T(z) for z in zz]
    plt.plot(zz, T0, label='T0')
    plt.plot(zz, T1, label='T1')
    plt.xlabel('z')
    plt.ylabel('T')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    testCIU()
