# Copyright (c) 2024 James S. Prell, all rights reserved

from .spa_const import const
from collections import defaultdict
from pathlib import Path
import random, math
import numpy as np
from os import path

class cellclass():
    '''This is a new generic cell parent class to define a cell interface that can be used
    in the modelclass for doing the collision heating calculations.
    Options:
    1. I could make this work as a cell with constant density, temperature and field. (simple CID cell).
    2. Or, could force every cell type to be defined separately.
    I think option 1 makes sense to me.
    '''
    min_dT = 0.1        # minimum change in T to trigger velocity distribution recalculation

    def __init__(self, **args):
        # args is a dict with any named parameters supplied by caller
        # a simple cell with constant fields can initialize as:
        #       cell = cellclass(L=0.5, Mgamu=28, T=298, rho=2.42e-5, E=5, z0=0)
        #       cell = cellclass(L=0.5, Mg=4.65e-26, T=298, pressure=2.133e-5, E=5, z0=0)
        # override in child class for more complex cell definitions
        self.type = ''  # base type is empty string
        pdict = defaultdict(lambda:None)    # return None for any key not defined in the dict
        pdict.update(args)
        self.args = args
        self.L = pdict.get('L', 0.5)              # length of cell; units m
        gas = pdict.get('gas', None)
        Mgamu = pdict.get('Mgamu', 28)            # mass of the cell gas in amu units
        Mg = pdict.get('Mg', None)                # mass of the cell gas in kg
        self.Mg = Mg or Mgamu*const.mamu          # mass of the cell gas in kg units
        glist = 'ar n2 nn he'.split()
        try:
            gas = gas.lower()
            if gas not in glist:
                gas = None      # not a recognized gas
        except:
            gas = None
        # now gas is either None or one of the standard ones in lowercase.
        if gas == 'ar':
            Mgamu = 39.95
        elif gas == 'n2':
            Mgamu = 28.0134
        elif gas == 'nn':
            Mgamu = 28.0134
        elif gas == 'he':
            Mgamu = 4.002602
        elif gas == 'xe':
            Mgamu = 131.293
        else:
            # try to guess gas from the mass.
            m_int = int(0.5 + self.Mg / const.mamu)
            if m_int == 40:
                gas = 'ar'
                Mgamu = 39.95
            elif m_int == 28:
                gas = 'n2'
                Mgamu = 28.0134
            elif m_int == 4:
                gas = 'he'
                Mgamu = 4.002602
            elif m_int == 131 or m_int == 132:
                gas = 'xe'
                Mgamu = 131.293
            else:       # final catch. No known gas, mass not on list, default to 'ar'
                gas = 'ar'
                Mgamu = 39.95
        self.gas = gas
        self.Mg = Mg or Mgamu*const.mamu
        T = pdict.get('T', 298.0)
        P = pdict.get('pressure', None)           # the formula for rho needs pressure in Bar
        rho = pdict.get('rho', None)
        test = P or rho
        if test is None:    # both not supplied
            rho = 2.42e-5           # use this if neither pressure or rho is given
            P = rho * 1e-5 * T * const.kB / self.Mg
        elif rho is None:   # pressure given rho absent
            rho = P * 1e5 * self.Mg / const.kB / T      # assumes constant T
        else:               # rho given, pressure absent
            out = f'pressure {P} or rho {rho} must be specified.'
            print(out)
            P = rho * 1e-5 * T * const.kB / self.Mg
        z0 = pdict.get('z0',0.0)
        E = pdict.get('E', None) or pdict.get('dV', 7) / (self.L - z0) # priority is to specify E directly, otherwise use dV/L, with fallback

        self.P = P
        self.T = lambda X: T                 # a function of position vector and time; units K
        self.rho = lambda X: rho        # a function of position vector; units kg/m^3
        self.z0 = z0                    # cell start z location (x,y assumed = 0.0, zfinal assumed= z0+L)
#        self.V = lambda X, t: 0         # a function of position vector and time; units V
        self.Ez = lambda X, t: E        # a function of position vector and time derived from V; units V/m
        self.vgz = lambda X: 0          # default is no gas flow, so constant 0

        # internal values for doing velocity distribution
        self._Tglast = -1                # use to save most recent gas temperature used
        self._Tilast = -1               # use to save most recent ion temperature used
        self._Nvel = 2500
        self._vrange = np.arange(0,self._Nvel)
        self._gprobs = None
        self._aprobs = None

    def __str__(self):
        z0 = self.z0
        outf = 'L={} m; z0={} m; T(z0)={} K; rho(z0)={:0.3e} kg/m^3; Ez(z0)={:0.2f} V/m; Mg={:0.3e} kg; P={} bar; gas={}'
        return outf.format(self.L, z0, self.T(z0), self.rho(z0), self.Ez(z0,0), self.Mg, self.P, self.gas)

    def __repr__(self):
        return str(self.args)
        
    def _veldist(self, T, m):
        '''return an array representing the probability distribution of velocities
        for mass mgkg and pseudo-atom mass ma temperatures Tg and Ti.

        New code wants two different velocity distributions.  But the place this is used can
        only use one distribution at a time.
        But since we are trying to save work when dT from last call is small, we need to
        distinguish between the two distributions.

        The first distribution is for the gas in the cell and uses the gas mass self.Mg
        The second is for the ion pseudo-atom and uses the ion effective mass ma.
        '''

        # ------------------------------------------------------------------------
        # Boltzmann distribution of gas velocities for gas with mass mg, at
        # temperature, Tg
        # ------------------------------------------------------------------------
        kb = const.kB
        pi = np.pi
        v = self._vrange    # saved velocity range array

        pg = ((m/(2*pi*kb*T))**(3/2))*4*pi*v**2*np.exp(-v**2*m/(2*kb*T))
        probs = np.array(pg)
        probs /= probs.sum()

        return probs

    def randomgasvel(self, z, Ti, ma):
        '''Return z,x,y components of velocity for a background gas molecule with
        gas mass at temperature Tg [K]. Sampled from a Boltzman distribution
        of velocities.
        Also return z,x,y components of velocity for the pseudoatom at the 
        ion temperature Ti and pseudoatom mass ma.
        The vectors are uniformly distributed with a Maxwel-Boltzmann distrubution.
        As such, these are correct in the frame of the moving gas.'''
        Nvel = 2500
        Tg = self.T(z)
        # select gas unit vector uniformly from a cylindrical districution, convert to rectangular coordinates
        zg = random.uniform(-1, 1)
        phi = random.uniform(0, 2*np.pi)
        xg = math.cos(phi) * math.sin(math.acos(zg))
        yg = math.sin(phi) * math.sin(math.acos(zg))

        dTg = abs(self._Tglast - Tg)
        if dTg > self.min_dT:
            self._gprobs = self._veldist(Tg, self.Mg)
            self._Tglast = Tg   # save current temperature for next call
        gprobs =  self._gprobs

        dTi = abs(self._Tilast - Ti)
        if dTi > self.min_dT:
            self._aprobs = self._veldist(Ti, ma)
            self._Tilast = Ti   # save current temperature for next call
        aprobs =  self._aprobs

        # sample velocities from 0 to 2000 m/s from a Boltzmann distribution
        vg = np.random.choice(Nvel, 1, p=gprobs)[0]
        va = np.random.choice(Nvel, 1, p=aprobs)[0]
#        vg, va = np.random.choice(Nvel, 1, p=self._veldist(Tg, Ti, ma))[0]          # use single scalar value

        # scale gas unit vector by velocity
        zgas = vg*zg
        xgas = vg*xg
        ygas = vg*yg

        za = random.uniform(-1, 1)
        phi = random.uniform(0,2*np.pi)
        xa = math.cos(phi)*math.sin(math.acos(za))
        ya = math.sin(phi)*math.sin(math.acos(za))

        # scale gas unit vector by velocity
        zatom = va*za
        xatom = va*xa
        yatom = va*ya

        return zgas, xgas, ygas, zatom, xatom, yatom

    def pseudoatom_mass(self,temp,veli):
        '''Return parameters for calculating pseudo-atom mass [units are Daltons]
        pa, pb = cell.pseudomass_params()
        ma = (pa -pb*vn) / 1000 / const.Av  # /1000 to get mass in kg/mol, /Av to get mass in kg
        
        Currently, the code is set up to use the 'full' fit of the MD pseudo-atom mass. To use the 'select' fit, uncomment lines 154 and 170, and comment
        out lines 156 and 172. If you want to use the 'piecewise' fit, comment out lines 156, 172, 154, and 170 and uncomment lines 158, 159, 174, 175,
        184-197, and comment out line 199
        '''
#        m_int = int(self.Mg / const.mamu)
        a, b, c = 0.0013, -1.1859, 1297.5 # default PA parameter values (same as argon gas parameters)
        d, e, f = -2e-05, 0.0324, 76.606
        # or build into the dict passed in to the cell __init__ function, possibly from the json file
#        if m_int == 40:
        if self.gas == 'ar':
            '''
            a, b, c = 0.012406, -0.011740, 92.511274    #argon piecewise fit, this is the 'upper' fit
            d, e, f = -0.022176, 0.056141, 39.375040    #this is the lower fit# the given values are for Argon
            '''
            #AMBER PARAMETERS
            # a, b, c = 0.0013, -1.1859, 1297.5
            # d, e, f = 1e-05, -0.0043, 85.667
        
            #NEW AMBER PARAMETERS
            a, b, c = -0.0042, 6.2829, 338.07
            d, e, f = 6e-06, -0.0169, 67.152
            
            #OPLS PARAMETERS
            #a, b, c = 0.003, -4.443, 2899.4 
            #d, e, f = -1e-05, 0.0257, 66.676
            
#        elif m_int == 28:
        elif self.gas == 'n2' or self.gas == 'nn':
            
            '''
            # a, b, c = -0.000876, 0.001141, 65.523096    #nitrogen piecewise fit, this is the 'upper' fit
            # d, e, f = -0.030848, 0.63956, 9.972156      #this is the lower fit
            
            # Diatomic PA piecewise fit numbers
            a, b, c = -0.002426, 0.004578, 46.641539    #nitrogen piecewise fit, this is the 'upper' fit
            d, e, f = -0.041589, 0.059209, 8.778090      #this is the lower fit
            '''
            
            # Unified PA fit numbers
            # a, b, c = -0.0019, 2.7834, 836.36 
            # d, e, f = -2e-05, 0.0311, 59.267 
            
            #NEW AMBER Unified PA fit numbers
            a, b, c = -3e-17, 4e-14, 2599.9
            d, e, f = -5e-06, -0.0048, 54.362
            
            # Diatomic PA fit numbers
            # a, b, c = -0.0006, 1.9499, 529.2 
            # d, e, f = -7e-05, 0.0851, 32.986
            
            #NEW AMBER Diatomic PA fit numbers
            # a, b, c = -3e-17, 4e-14, 2599.9
            # d, e, f = 3e-07, -0.0113, 46.57

#        elif m_int == 4:
        elif self.gas == 'he':
            # a, b, c = -3e-17, 4e-14, 2599.9
            # d, e, f = -4e-06, -0.0055, 19.479
            
            # NEW AMBER HELIUM
            a, b, c = 0.0109, -14.88, 6066.4
            d, e, f = 0.001, -0.3666, 34.867
        
        elif self.gas == 'xe':
            # a, b, c = 0, 0, 750.5
            # d, e, f = -9e-5, 0.1816, 122.97
            
            # NEW OPLS XENON
            a, b, c = 0.0009, -0.7966, 1836.8
            d, e, f = 3e-06, -0.0279, 155.8
            
        '''
        direction = np.cross([a,b,-1],[d,e,-1])
        
        y0 = -(a*(f - c)/(a - d) + c)/(a*(e - b)/(a - d) + b)
        x0 = ((e - b)*y0 + (f-c))/(a - d)
        
        t = (temp - x0)/direction[0]
        
        vel_check = y0 + t*direction[1]
        
        if veli >= vel_check:
            ma = a*temp + b*veli + c
            # print('here')
        if veli < vel_check:
            ma = d*temp + e*veli + f
        
        #ma = a*temp + b*veli + c
        '''
        
        coeff2_T = a*temp**2 + b*temp + c 
        ma_T = d*temp**2 + e*temp + f
        coeff1_T = np.exp(1)*(ma_T - (self.Mg/const.mamu))/coeff2_T
        ma = coeff1_T*veli*np.exp(-(veli/coeff2_T)) + (self.Mg/const.mamu)
        
        return ma #this will be in daltons
         
class fullAcell(cellclass):
    '''Cell class for the three hexapoles of an agilent 6545xt'''
    
    def __init__(self, **args):
        self.type = 'fullAcell'
        pdict = defaultdict(lambda:None)    # return None for any key not defined in the dict
        pdict.update(args)
        super(fullAcell,self).__init__(L=pdict['L'], Mgamu=pdict['Mgamu'], T=pdict['T'], rho=pdict['rho'], z0=pdict['z0'])
        # self.rho0 = pdict.get('rho', 4.18e-5)
        # self.rho1 = pdict.get('rho1', 4.1e-7)
        # self.rho2 = pdict.get('rho2', 4.1e-7)
        T = pdict['T']
        self.rho0 = pdict.get('rho', None)   
        self.rho1 = pdict.get('rho2', None)  
        self.rho2 = pdict.get('rho3', None)  
        self.P0 = pdict.get('pressure1', None)
        self.P1 = pdict.get('pressure2', None)
        self.P2 = pdict.get('pressure3', None)
        self.rho = self.rho0
        test = self.P0 or self.rho0
        if test is None:    # both not supplied
            self.rho0 = 2.42e-5           # use this if neither pressure or rho is given
            self.P0 = self.rho0 * 1e-5 * T * const.kB / self.Mg
            print(f'YOU DID NOT SPECIFY A PRESSURE OR RHO, USING DEFAULT VALUES OF P: {self.P0:0.2f} rho: {self.rho0:0.2f}')
            print(f'PRESSURE OR RHO MUST BE SPECIFIED')
        elif self.rho0 is None:   # pressure given rho absent
            self.rho0 = self.P0 * 1e5 * self.Mg / const.kB / T      # assumes constant T
            self.rho1 = self.P1 * 1e5 * self.Mg / const.kB / T      # assumes constant T
            self.rho2 = self.P2 * 1e5 * self.Mg / const.kB / T      # assumes constant T
        else:               # rho given, pressure absent
            self.P0 = self.rho0 * 1e-5 * T * const.kB / self.Mg
            self.P1 = self.rho1 * 1e-5 * T * const.kB / self.Mg
            self.P2 = self.rho1 * 1e-5 * T * const.kB / self.Mg
        
        
        
        self.Ez = self.wEz
        self.Ez1 = pdict.get('E1', 38.89)
        self.Ez2 = pdict.get('E2', 28.57)
        self.Ez3 = pdict.get('E3', 0)
        self.rho = self._rho
        
        self.h3apt = 2.5
        self.zlocs = np.array([0.182, 0.2155])
        self.dVs = np.array([1, 1])
        self.Emags = 1.5*0.896*self.dVs / (self.h3apt/1000)
        self.c = 0.005
        
        self.zshift = self.h3apt
        self.zbarrier = 218
        self.zpdrop = self.zbarrier - self.zshift
        
        
        
    def _rho(self, z):
        
        if z < 0.182:
            rhoz = self.rho0
        
        if z >= 0.182 and z < 0.2155:
            dz = z - .182
            rhoz = self.rho1 + (self.rho0 - self.rho1)*(0.0025/(0.002 + (z - 0.1825)))**2
        if z >= 0.2155:
            dz = z - self.zpdrop
            rhoz = self.rho2 + (self.rho1 - self.rho2)*(0.0025/(0.002 + (0.2155 - dz)))**2
            
        
        return rhoz
    
    def wEz(self, z, t):

        
        
        if z < 0.182:
            e = self.Ez1
        if z >= 0.182 and z < 0.2155:
            e = self.Ez2
        if z >= 0.2155:
            e = self.Ez3
            
        for zloc, Emag in zip(self.zlocs, self.Emags):
            if Emag == 0:
                continue
            # print(zloc)
            de = 4.6*(self.c**2)*Emag/(1000*(z - zloc)**2 + self.c**2)
            e += de
        return e
        
    

class WAVEcell(cellclass):
    '''Derived cell class for the Waters WAVE cell type.  Consider moving to separate module.'''
    '''myoglobin: mi = 2.020000e-17 kg, z = 8 or 9+,
        pressure = 2.133e-05 bar, wave velocity = 300 m/s, wave height = 2 V
        DEPRECIATED, PLEASE USE WAVEcell4'''
    def __init__(self, **args): # L=0.125, wl = 0.0121, wh = 2, vw = 300, Mgamu=28, T=298, rho=4.18e-5, z0=0):
        self.type = 'WAVE'
        pdict = defaultdict(lambda:None)    # return None for any key not defined in the dict
        pdict.update(args)
        super(WAVEcell,self).__init__(L=pdict['L'], Mgamu=pdict['Mgamu'], T=pdict['T'], rho=pdict['rho'], z0=pdict['z0'])
        self._vfactor = 0.896
        self.wl = pdict.get('wl', 0.0121)
        self.wh = pdict.get('wh', 2.0)
        self.vw = pdict.get('vw', 300.0)
        self.k = 2 * math.pi / self.wl
        self.w = self.k * self.vw
        self.V = self.wV
        self.Ez = self.wEz

    def __str__(self):
        outs = super(WAVEcell, self).__str__() + '\n'
        outs += 'WAVE: wl={}, wh={}, vw={}'.format(self.wl, self.wh, self.vw)
        return outs

    def wV(self, z, t):    #  z, vw, wl, ts=0, t=0, wh=0
        '''goal here is for this to be the Voltage for TWIMs version.
        As is, is voltage at the cell rings, but needs to be modified
        to include decrease in amplitude at cell center.
        Paper says the center voltage is 0.896 times edge voltage.'''
        Vt = self._vfactor * self.wh / 2 * np.sin(self.k * (z+self.z0) - self.w * t) #note that Mortensen et al. JASMS 2017, 1282 says this isn't very accurate
        return Vt

    def wEz(self, z, t):
        #e = self._vfactor * self.k * self.wh / 2 * np.cos(self.k * (z+self.z0) - self.w * t)
        e = self.wh*(213.7*np.cos(self.k * (z+self.z0) - self.w * t)-25.1*np.cos(3*self.k * (z+self.z0) - 3*self.w * t)+2.0*np.cos(5*self.k * (z+self.z0) - 5*self.w * t)) #from SIMION via Mortensen et al. JASMS 2017, 1282
        return e


class WAVEcell2(cellclass):
    '''Derived cell class for the Waters WAVE cell type.  Consider moving to separate module.
       pressure = 2.133e-05 bar, wave velocity = 300 m/s, wave height = 2 V
       Modified version to model pressure changes through cell.
       DEPRECIATED, PLEASE USE WAVEcell4'''
    def __init__(self, **args): # L=0.125, wl = 0.0121, wh = 2, vw = 300, Mgamu=28, T=298, rho=4.18e-5, z0=0):
        self.type = 'WAVE2'
        pdict = defaultdict(lambda:None)    # return None for any key not defined in the dict
        pdict.update(args)
        self.rho0 = pdict.get('rho', 4.18e-5)
        self.rho1 = pdict.get('rho1', 4.1e-7)
        self.zopen = pdict.get('zopen', 0.09)
        super(WAVEcell2,self).__init__(L=pdict['L'], Mgamu=pdict['Mgamu'], T=pdict['T'], rho=pdict['rho'], z0=pdict['z0'])
        self._vfactor = 0.896
        self.wl = pdict.get('wl', 0.0121)
        self.wh = pdict.get('wh', 2.0)
        self.vw = pdict.get('vw', 300.0)
        self.k = 2 * math.pi / self.wl
        self.w = self.k * self.vw
        self.V = self.wV
        self.Ez = self.wEz
        self.rho = self._rho        # a function of position vector; units kg/m^3
        # self.savez = 0.0


    def __str__(self):
        outs = super(WAVEcell2, self).__str__() + '\n'
        outs += 'WAVE: wl={}, wh={}, vw={}'.format(self.wl, self.wh, self.vw)
        return outs

    def wV(self, z, t):    #  z, vw, wl, ts=0, t=0, wh=0
        '''goal here is for this to be the Voltage for TWIMs version.
        As is, is voltage at the cell rings, but needs to be modified
        to include decrease in amplitude at cell center.
        Paper says the center voltage is 0.896 times edge voltage.'''
        Vt = self._vfactor * self.wh / 2 * np.sin(self.k * (z+self.z0) - self.w * t) #note that Mortensen et al. JASMS 2017, 1282 says this isn't very accurate
        return Vt

    def wEz(self, z, t):
        #e = self._vfactor * self.k * self.wh / 2 * np.cos(self.k * (z+self.z0) - self.w * t)
        e = self.wh*(213.7*np.cos(self.k * (z+self.z0) - self.w * t)-25.1*np.cos(3*self.k * (z+self.z0) - 3*self.w * t)+2.0*np.cos(5*self.k * (z+self.z0) - 5*self.w * t)) #from SIMION via Mortensen et al. JASMS 2017, 1282
        return e

    def _rho(self, z):
        '''Return temperature at z location, interpolating from file data.'''
        dz = z - self.zopen
        if dz <= 0:
            rhoz = self.rho0
        else:
            # start with simple quadratic decrease in density following division to open part of cell
            rhoz = self.rho1 + (self.rho0 - self.rho1) * (0.003/(0.003+dz))**2
        
        # if z > (self.savez + 0.001): 
        #     print(f'        {dz}   {self.rho0}   {self.rho1}       {rhoz}')
        #     self.savez = z
        # else:
        #     pass
        return rhoz

class WAVEcell3(cellclass):
    '''Derived cell class for the Waters WAVE cell type. 
       Modified version to model pressure changes through extended cell, 
       including the trap cell, ion guides, low pressure IMS cell, 
       and final transfer cell.
       Low pressure through ion guides and IMS cell should result in low
       collision frequency, but require sufficient model time (and tmax for recording)
       to account for full transit through the system.
       DEPRECIATED, PLEASE USE WAVEcell4'''
    def __init__(self, **args):
        '''Arguments should specify basic cell arguments:
            L:  total cell length,
            Mgamu: gas mass in amu,
            T: gas temperature in cell, 
            rho: gas density in cell [kg/m^3], 
            z0: usually 0.0 -- ion position at cell opening
        WAVE arguments:
            wl: wave length [m], 
            wh: wave height [V], 
            vw: wave velocity [m/s]
        and extended arguments:
            zopen: z position where gas density drops (location of trap exit aperture),
            zclose: z position where gas density returns to initial pressure (aperture location)
            zV1: z position for first voltage change    - 4.5 cm
            dV1: amount of first voltage change         - -1 V
            zV2: z position for 2nd voltage change      - same as zopen, 9 cm
            dv2: amount of 2nd change                   - -2 V
            zV3: z position for 3rd voltage change      - 13 cm
            dv3: amount of 3rd change                   - -1 V
            zV4: z position for 4th voltage change      - 38 cm
            dv4: amount of 4th change                   - 0 V
            zV5: z position for 5th voltage change      - 42 cm -- same as zclose
            dv5: amount of 5th change                   - +3 V
            zwaveoff: position where travelling wave stops
            zwaveon: position where travelling wave restarts
            '''
        self.type = 'WAVE3'
        pdict = defaultdict(lambda:None)    # return None for any key not defined in the dict
        pdict.update(args)
        self.rho0 = pdict.get('rho', 3.213e-5)   # This is pressure of trap and transfer regions
        self.rho1 = pdict.get('rho1', 3.213e-7)   # This is base pressure for ion guides and IMS cell
        self.zopen = pdict.get('zopen', 0.09)   # location of trap exit aperture
        self.zclose = pdict.get('zclose', 0.42)
        super(WAVEcell3,self).__init__(L=pdict['L'], Mgamu=pdict['Mgamu'], T=pdict['T'], rho=pdict['rho'], z0=pdict['z0'])
        self.wl = pdict.get('wl', 0.0121)
        self.wh = pdict.get('wh', 4.0)
        self.vw = pdict.get('vw', 311.0)
        self.k = 2 * math.pi / self.wl
        self.zV1 = pdict.get('zV1', 0.045)
        self.dV1 = pdict.get('dV1', -1)
        self.zV2 = pdict.get('zV2', 0.09)
        self.dV2 = pdict.get('dV2', -2)
        self.zV3 = pdict.get('zV3', 0.13)
        self.dV3 = pdict.get('dV3', -1)
        self.zV4 = pdict.get('zV4', 0.38)
        self.dV4 = pdict.get('dV4', 0)
        self.zV5 = pdict.get('zV5', 0.42)
        self.dV5 = pdict.get('dV5', 3)
        self.zwaveoff = pdict.get('zwaveoff', 0.13)
        self.zwaveon = pdict.get('zwaveon', 0.38)
            # cell end at 51 cm, voltage drop 5 V -- but no further collisions so not in model
        self._vfactor = 0.896
        self.d0 = 0.003          # approx the diameter of the aperture between high and low pressure zones
        self.zlocs = [self.zV1, self.zV2, self.zV3, self.zV4, self.zV5]
        self.dVs = [self.dV1, self.dV2, self.dV3, self.dV4, self.dV5]
        self.w = self.k * self.vw
        self.period = self.w/self.vw
        self.t0rand = self.period*random.random()    # a random phase adjustment for the rf in the cell
#        print(f'random time shift: {1e6 * self.t0rand:0.2f} us')
        self.V = self.wV
        self.Ez = self.wEz
        self.rho = self._rho        # a function of position vector; units kg/m^3
        # self.savez = 0.0


    def __str__(self):
        outs = super(WAVEcell3, self).__str__() + '\n'
        outs += f'WAVE: wl={self.wl}, wh={self.wh}, vw={self.vw}, zopen/close=[{self.zopen}, {self.zclose}]\n'
        outs += f'     [(zV,dV)= [({self.zV1},{self.dV1}), ({self.zV2},{self.dV2}), ({self.zV3},{self.dV3}), ({self.zV4},{self.dV4}), ({self.zV5},{self.dV5})]'
        return outs

    def wV(self, z, t):    #  z, vw, wl, ts=0, t=0, wh=0
        '''goal here is for this to be the Voltage for TWIMs version.
        As is, is voltage at the cell rings, but needs to be modified
        to include decrease in amplitude at cell center.
        Paper says the center voltage is 0.896 times edge voltage.'''
        Vt = self._vfactor * self.wh / 2 * np.sin(self.k * (z+self.z0) - self.w * t) #note that Mortensen et al. JASMS 2017, 1282 says this isn't very accurate
        return Vt

    def wEz(self, z, t):
        '''For WAVE3, need to add dc fields in defined dV transition areas.
        UNTESTED!!!
        Plan:
            at each location where dVn not zero, electric field adds dVn/self.d0 to the wave field
            Will need to check sign of this to ensure ions accelerate in correct direction!
        '''
        dzlocs = np.array([abs(z-zloc) for zloc in self.zlocs])
        mindz = dzlocs.min()
        minidx = dzlocs.argmin() # which one is least?
        if mindz < self.d0:                         # position relative to dV transition is +/- d0
            de = -0.5 * self.dVs[minidx] / self.d0  # add const electric field dV/(2*d0)
        else:
            de = 0
        tt = t + self.t0rand    # shift time for random phase of entrance relative to wave rf
        zz = z + self.z0

        #e = self._vfactor * self.k * self.wh / 2 * np.cos(self.k * (z+self.z0) - self.w * t)
        if zz < self.zwaveoff or zz > self.zwaveon:
            e = self.wh*(213.7*np.cos(self.k * zz - self.w * tt)-25.1*np.cos(3*self.k * zz - 3*self.w * tt)+2.0*np.cos(5*self.k * zz - 5*self.w * tt)) #from SIMION via Mortensen et al. JASMS 2017, 1282
        else:
            e = 0
        e += de
        return e

    def _rho(self, z):
        '''Return temperature at z location, interpolating from file data.'''
        if (z <= self.zopen - self.d0/3) or (z > self.zclose + self.d0/3):
            rhoz = self.rho0
            return rhoz
        if z < self.L / 2:
            dz = z + self.d0/3 - self.zopen
        else:
            dz = self.zclose - z + self.d0/3
            # start with simple quadratic decrease in density following division to open part of cell
        rhoz = self.rho1 + (self.rho0 - self.rho1) * (self.d0/(self.d0+dz))**2
        return rhoz
    
class WAVEcell4(cellclass):
    '''Derived cell class for the Waters WAVE cell type. 
       Modified version to model pressure changes through extended cell, 
       including the trap cell, ion guides, low pressure IMS cell, 
       and final transfer cell.
       Low pressure through ion guides and IMS cell should result in low
       collision frequency, but require sufficient model time (and tmax for recording)
       to account for full transit through the system.'''
    def __init__(self, **args):
        '''Arguments should specify basic cell arguments:
            L:  total cell length,
            Mgamu: gas mass in amu,
            T: gas temperature in cell, 
            rho: gas density in cell [kg/m^3], 
            z0: usually 0.0 -- ion position at cell opening
        WAVE arguments:
            wl: wave length [m], 
            wh1: wave height [V],           Parameters for trap Twave
            vw1: wave velocity [m/s]
            wh2: wave height [V],           Parameters for transfer Twave
            vw2: wave velocity [m/s]
        and extended arguments:
            zopen: z position where gas density drops (location of trap exit aperture),
            zclose: z position where gas density returns to initial pressure (aperture location)
            zV1: z position for first voltage change    - 4.5 cm
            dV1: amount of first voltage change         - -1 V
            zV2: z position for 2nd voltage change      - same as zopen, 9 cm
            dv2: amount of 2nd change                   - -2 V
            zV3: z position for 3rd voltage change      - 13 cm
            dv3: amount of 3rd change                   - -1 V
            zV4: z position for 4th voltage change      - 38 cm
            dv4: amount of 4th change                   - 0 V
            zV5: z position for 5th voltage change      - 42 cm -- same as zclose
            dv5: amount of 5th change                   - -2 V
            zV6: z position for 6th voltage change      - 46.5 cm -- same as zclose
            dv6: amount of 6th change                   - -1 V
            zwaveoff: position where travelling wave stops
            zwaveon: position where travelling wave restarts
            '''
        self.type = 'WAVE4'
        pdict = defaultdict(lambda:None)    # return None for any key not defined in the dict
        pdict.update(args)
        self.rho0 = pdict.get('rho', None)   # This is pressure of trap and transfer regions
        self.rho1 = pdict.get('rho1', None)   # This is base pressure for ion guides and IMS cell
        self.P0 = pdict.get('pressure1', None)
        self.P1 = pdict.get('pressure2', None)
        super(WAVEcell4,self).__init__(L=pdict['L'], Mgamu=pdict['Mgamu'], T=pdict['T'], rho=pdict['rho0'], z0=pdict['z0'])
        T = pdict['T']
        test = self.P0 or self.rho0
        if test is None:    # both not supplied
            self.rho0 = 2.42e-5           # use this if neither pressure or rho is given
            self.P0 = self.rho0 * 1e-5 * T * const.kB / self.Mg
            print(f'YOU DID NOT SPECIFY A PRESSURE OR RHO, USING DEFAULT VALUES OF P: {self.P0:0.2f} rho: {self.rho0:0.2f}')
            print(f'PRESSURE OR RHO MUST BE SPECIFIED')
        elif self.rho0 is None:   # pressure given rho absent
            self.rho0 = self.P0 * 1e5 * self.Mg / const.kB / T      # assumes constant T
            self.rho1 = self.P1 * 1e5 * self.Mg / const.kB / T      # assumes constant T
        else:               # rho given, pressure absent
            self.P0 = self.rho0 * 1e-5 * T * const.kB / self.Mg
            self.P1 = self.rho1 * 1e-5 * T * const.kB / self.Mg
        
        self.zopen = pdict.get('zopen', 0.09)   # location of trap exit aperture
        self.zclose = pdict.get('zclose', 0.42)
        self.wl = pdict.get('wl', 0.0121)
        self.wh1 = pdict.get('wh1', 4.0)
        self.vw1 = pdict.get('vw1', 311.0)
        self.wh2 = pdict.get('wh2', 0.2)
        self.vw2 = pdict.get('vw2', 247.0)
        self.k = 2 * math.pi / self.wl
        self.zV1 = pdict.get('zV1', 0.045)
        self.dV1 = pdict.get('dV1', -1)
        self.zV2 = pdict.get('zV2', 0.09)
        self.dV2 = pdict.get('dV2', -2)
        self.zV3 = pdict.get('zV3', 0.13)
        self.dV3 = pdict.get('dV3', 0)
        self.zV4 = pdict.get('zV4', 0.38)
        self.dV4 = pdict.get('dV4', 0)
        self.zV5 = pdict.get('zV5', 0.42)
        self.dV5 = pdict.get('dV5', -1)
        self.zV6 = pdict.get('zV6', .465)
        self.dV6 = pdict.get('dV6', -2)
        self.zV7 = pdict.get('zV7', .49)
        self.dV7 = pdict.get('dV7', 0)
        self.path_bool = False
        self.zwaveoff = pdict.get('zwaveoff', 0.13)
        self.zwaveon = pdict.get('zwaveon', 0.38)
            # cell end at 51 cm, voltage drop 5 V -- but no further collisions so not in model
        self._vfactor = 0.896
        self.d0 = 0.005       # approx the diameter of the aperture between high and low pressure zones
        self.c = 0.0007
        self.zlocs = np.array([self.zV1, self.zV2, self.zV3, self.zV4, self.zV5, self.zV6, self.zV7])
        self.dVs = np.array([self.dV1, self.dV2, self.dV3, self.dV4, self.dV5, self.dV6, self.dV7])
        self.Emags = -1.5*0.896*self.dVs / self.d0
        # self.t0rand = 10*(random.random()-0.5) * self.d0 / self.vw1   # a random phase adjustment for the rf in the cell
        # print(f'random time shift: {1e6 * self.t0rand:0.2f} us')
        self.w = self.k * self.vw1
        self.w1 = self.k * self.vw1
        self.w2 = self.k * self.vw2
        self.period1 = self.wl/self.vw1
        self.period2 = self.wl/self.vw2
        self.t0rand1 = self.period1*random.random() #generates a random time shift for wave1 phase in trap cell
        # self.t0rand1 = 10*(random.random()-0.5) * self.d0 / self.vw1
        # self.t0rand2 = self.period2*random.random() #generates a random time shift for wave2 phase in transfer cell
        # self.t0rand2 = 10*(random.random()-0.5) * self.d0 / self.vw2
        self.V = self.wV
        self.Ez = self.wEz
        self.rho = self._rho        # a function of position vector; units kg/m^3
        # self.savez = 0.0
        self.simion_array = pdict.get('Simion Potential Arrays', None)
#        self.array_file_path = Path(f'ionspa/{self.simion_array}')
        self.array_file_path = path.join(const.ionspadir, self.simion_array)
        if path.exists(self.array_file_path):
            self.path_bool = True
#        if self.array_file_path.is_file():
            self.Efield = np.loadtxt(self.array_file_path, delimiter=',')
            self.cell_zs = self.Efield[0]
            self.cell_efield = self.Efield[1]

        # print(self.Emags)
    def __str__(self):
        outs = super(WAVEcell4, self).__str__() + '\n'
        outs += f'WAVE: wl={self.wl}, wh={self.wh1}, vw={self.vw1}, zopen/close=[{self.zopen}, {self.zclose}]\n'
        outs += f'     [(zV,dV)= [({self.zV1},{self.dV1}), ({self.zV2},{self.dV2}), ({self.zV3},{self.dV3}), ({self.zV4},{self.dV4}), ({self.zV5},{self.dV5})]'
        return outs

    def wV(self, z, t):    #  z, vw, wl, ts=0, t=0, wh=0
        '''goal here is for this to be the Voltage for TWIMs version.
        As is, is voltage at the cell rings, but needs to be modified
        to include decrease in amplitude at cell center.
        Paper says the center voltage is 0.896 times edge voltage.'''
        Vt = self._vfactor * self.wh1 / 2 * np.sin(self.k * (z+self.z0) - self.w * t) #note that Mortensen et al. JASMS 2017, 1282 says this isn't very accurate
        return Vt

    def wEz(self, z, t):
        '''For WAVE3, need to add dc fields in defined dV transition areas.
        UNTESTED!!!
        Plan:
            at each location where dVn not zero, electric field adds dVn/self.d0 to the wave field
            Will need to check sign of this to ensure ions accelerate in correct direction!
        '''
        
        tt1 = t + self.t0rand1    # shift time for random phase of entrance relative to wave rf for trap
        tt2 = t #+ self.t0rand2    # shift time for random phase of entrance relative to wave rf for transfer
        zz = z + self.z0
        
        '''
        dzlocs = np.array([abs(z-zloc) for zloc in self.zlocs])
        mindz = dzlocs.min()
        minidx = dzlocs.argmin() # which one is least?
        if mindz < self.d0:                         # position relative to dV transition is +/- d0
            de1 = -0.5 * self.dVs[minidx] / self.d0  # add const electric field dV/(2*d0)
        else:
            de1 = 0
        '''
        de1 =0
        if self.path_bool:
            de1 += np.interp(zz, self.cell_zs, self.cell_efield)
            # print('Using SIMION potential arrays')
        else:
            for zloc, Emag in zip(self.zlocs, self.Emags):
                if Emag == 0:
                    continue
                de1 += (self.c**2)*Emag/((zz - zloc)**2 + self.c**2)
            # print('Using Emag formula')
        
        
        

        #e = self._vfactor * self.k * self.wh / 2 * np.cos(self.k * (z+self.z0) - self.w * t)
        if zz < self.zwaveoff:
            e = self.wh1*(213.7*np.cos(self.k * zz - self.w1 * tt1)-25.1*np.cos(3*self.k * zz - 3*self.w1 * tt1)+2.0*np.cos(5*self.k * zz - 5*self.w1 * tt1)) #from SIMION via Mortensen et al. JASMS 2017, 1282
            # print(e)
        elif zz > self.zwaveon:
            e = self.wh2*(213.7*np.cos(self.k * zz - self.w2 * tt2)-25.1*np.cos(3*self.k * zz - 3*self.w2 * tt2)+2.0*np.cos(5*self.k * zz - 5*self.w2 * tt2)) #from SIMION via Mortensen et al. JASMS 2017, 1282
        else:
            e = 0
            
        e += de1# + de2
        # print(de2)
        # if e!t(e)
        return e

    def _rho(self, z):
        '''Return temperature at z location, interpolating from file data.'''
        if (z <= self.zopen - self.d0/5) or (z > self.zclose + self.d0/5):
            rhoz = self.rho0
            return rhoz
        if z < self.L / 2:
            dz = z + self.d0/5 - self.zopen
        else:
            dz = self.zclose - z + self.d0/5
            # start with simple quadratic decrease in density following division to open part of cell
        rhoz = self.rho1 + (self.rho0 - self.rho1) * (self.d0/(self.d0+dz))**2
        return rhoz

