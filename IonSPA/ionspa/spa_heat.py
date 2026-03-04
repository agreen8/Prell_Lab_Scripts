# -*- coding: utf-8 -*-
####################################
#
# spa_heat as part of the IonSPA (ionspa) system
#
# Original work by Jim Prell, Sam Shepherd, et.at. at University of Oregon
# Copyright (c) 2024 James S. Prell, all rights reserved
#
# Modifications by Ken Newton, University of Oregon
#   cell, ion, and model are implemented as classes
#
# This script can be run standalone for testing as:
#     python -m ionspa.spa_heat
####################################

import sys, os
import numpy as np
import scipy.special as special
import random
import math

# make other file functions and classes appear as localss to this file
import ionspa
const = ionspa.const
hcprofileclass = ionspa.hcprofileclass
cellclass = ionspa.cellclass
CIUcell = ionspa.CIUcell
WAVEcell = ionspa.WAVEcell
makecell = ionspa.makecell


# def makecell(celldict):
#     '''Given a celldict with either a base "cell" type or a Waters cell "wcell" type, return a cell object.'''
#     if celldict is None:
#         cell = cellclass()
#     elif celldict['type'] == "cell":
#         cell = cellclass(**celldict)
#     elif celldict['type'] == "wcell":
#         cell = WAVEcell(**celldict)
#     elif celldict['type'] == "wcell2":
#         cell = WAVEcell2(**celldict)
#     elif celldict['type'] == "wcell3":
#         cell = WAVEcell3(**celldict)
#     elif celldict['type'] == "wcell4":
#         cell = WAVEcell4(**celldict)
#     elif celldict['type'] == "aCIUcell":
#         cell = CIUcell(**celldict)
#     return cell


def fracremains(times, Temps, deltaH, deltaS, maxtime=None):
    '''Return remaining fraction of ions.
    times is an ion time array in s  (computed for ion, cell, Vin, etc.)
    Temps is a corresponding ion temperature array in K
    deltaH is the modeled deltaH in J/mol
    deltaS is the modeled deltaS in J/mol/K
    if maxtime is not None, only return fraction remaining until maxtime
    '''
    kBoverh = 20836810339.4  # kB/h in units of 1/(K*s)
    idgasR = const.R  # in J/mol*K

    t = 0
    Ntimes = len(times)
    r = np.zeros(Ntimes)
    r[0] = 1.0
    finalr = 1.0
    deltat = times[1]-times[0]
    while r[t] >= 2e-6 and r[t] <= 1.0 and t < Ntimes-1:
        # can change power of temperature for non-Arrhenius, pure Arrhenius would just be T in prefactor
        r[t+1] = r[t]*(1.0 - kBoverh*deltat*Temps[t]*np.exp(deltaS/idgasR)*np.exp(-deltaH/(idgasR*Temps[t])))
        if r[t+1] <= r[t] and r[t+1] >= 0.0:
            finalr = r[t+1]
            t += 1
        else:
            # print(deltaS, deltaH/Temps[t], t, Temps[t], r[t], r[t+1])
            if t < Ntimes/10:
                finalr = 1e-6
            break
        if maxtime is not None:
            if times[t] > maxtime:
                break  # don't go past maxtime
    return finalr


# default values for deltaH, deltaS
deltaH = 80
deltaS = -120


def fracloss(deltat, T, deltaH, deltaS):
    '''Compute fractional loss (fragmentation) in a single time step (deltat)
    given current temperature T and assumed deltaH, deltaS.
    Note this function uses deltaH in kJ/mol. and deltaS in J/mol/K
    The function should work for either scalar or arrays for deltat and T.
    '''

#    kBoverh = const.kBoverh     # kB/h in units of 1/(K*s)
#    idgasR = const.R            # in J/mol*K
    dH = deltaH * 1000.0        # deltaH assumed in kJ/mol, convert to J/mol
    dHscale = dH/const.R
    dSscale = deltaS/const.R
    dr = const.kBoverh*deltat*T*np.exp(dSscale)*np.exp(-dHscale/T)
    # note for large time steps, this returns fractional loss that can be large (>1)
    # in that case, it is important to treat the loss as
    # loss = 1-exp(-dr) to model the continual loss rate with large time steps
    return dr


def lossrate(T, deltaH, deltaS):
    '''Just compute a loss rate at this temperature given dH, dS. Same as fracloss but without deltat.
    So it is a dissociation rate in units 1/s
    Like fracloss, should work with T as np.array
    '''

    dHscale = 1000.0*deltaH/const.R
    dSscale = deltaS/const.R
    dr = const.kBoverh*T*np.exp(dSscale -dHscale/T)
    return dr


class ionclass():
    '''The ionclass should encapsulate all the variables directly tied to the ion properties. This includes a name
    for the ion and mass, charge, collision cross section, number of atoms, and initial internal temperature.
    We choose to design the system with flight and energy transfer functions in a separate object.'''
    def __init__(self,
                 iondict=dict(name='hp1821Ag+213e', mass=1.822, charge=1, CCS=3.5,
                              num_atoms=126, T0=298, hcprofname='tunemix', refdH=deltaH, refdS=deltaS)):
        #     inputs = [dict(name='hp1821Ag+213e', mass=1.822, charge=1, CCS=3.5,
        #                     num_atoms=126, T0=298,hcprofname='tunemix')]
        self.name = iondict['name']
        self.mi = iondict['mass']/const.Av   # mass in kg (input in kg/mol)
        self.Z = iondict['charge']*const.qe  # charge in C (input in num charges)
        self.ccs = iondict['CCS']*1e-18      # ccs in m^2 (input was in nm^2)
        self.natoms = iondict['num_atoms']
        self.T0 = iondict['T0']              # units: K
        self.hcprofname = iondict['hcprofname']
        self.hcprofileobject = hcprofileclass()
        self.dof = 3 * self.natoms - 6
        self.hcprofileobject.buildfunctions(self.hcprofname, self.dof)
        self.dH = iondict.get('refdH', deltaH)
        self.dS = iondict.get('refdS', deltaS)

    def asdict(self):
        return dict(name=self.name, mass=self.mi*const.Av, charge=self.Z, CCS=self.ccs*1e18, num_atoms=self.natoms,
                    T0=self.T0, hcprofname=self.hcprofname)

    def __repr__(self):
        return str(self.asdict())

    def UfromT(self, temperature):
        return self.hcprofileobject.UfromT(temperature)

    def TfromU(self, U):
        return self.hcprofileobject.TfromU(U)


class modelclass():
    '''A modelclass holds an ionclass object and cellclass object. In addition,
    it defines the starting Vin (as V, multiply by charge for actual start energy).
    Functions for the user include start(), and steptime(), which can reinitialize
    or single-step the model and runmodel() to do a full run.
    __str__() provides a single-line output string for run state, with header() the corresponding header line.
    __repr__() gives a more verbose summary of the model variables.

    Recorded variables for each timestep can be extracted between calls to steptime() if desired.
    '''
    def __init__(self, cellobj, ionobj, Vin, ma_kg_mol=0.032):
        self.cell = cellobj
        self.ion = ionobj
        self.Vin = Vin
        self.KE = Vin * self.ion.Z      # Z is charge in C, KE is in J
        self.ma_kg_mol = ma_kg_mol
        self.ma = self.ma_kg_mol/const.Av

        self.start(z=0, dVin=Vin)
        self.collided = False
        self.coeffrest = 1.0
        self.frac = 1   # start with all ions unfragmented
        self.fragmented = False     # boolean to indicate if the ion would have fragmented at current time
        self.dfrac = 0.0            # probability of fragmentation at this timestep

    def headerstr(self):
        '''A single line header string for use with __str__()'''
        return 'tcount ccount t z ts vz mfp u0 u umax T0 T Tmax KE u/kT ma  f'

    def __str__(self):
        '''A single line model state output.'''
        outs = '{:5} {:4}  {:8.2e}  {:7.5f}  {:7.1e}  {:7.3f}  {:8.2e}  {:7.3f} {:7.3f} {:7.3f}'\
            '  {:5.1f} {:5.1f} {:5.1f}  {:5.3f}  {:6.1f}  {:5.2f}  {:5.3f}'.format(
                self.tcount,
                self.collisioncount,
                self.t,
                self.z,
                self.ts,
                self.vz,
                self.mfp(self.z),
                self.u0/const.qe,
                self.u/const.qe,
                self.umax/const.qe,
                self.T0,
                self.temp,
                self.Tmax,
                self.KE/const.qe,
                self.u/const.kB/self.temp,
                self.ma*const.Av*1000,
                self.frac
                )
        return outs

    def __repr__(self):
        '''A multi-line model summary.'''
        outs = '\n'
        outs += 'cell = ' + str(self.cell) + '\n'
        outs += 'ion = ' + self.ion.__repr__() + '\n'
        outs += self.headerstr() + '\n'
        outs += self.__str__() + '\n'
        outs += f'   pseudoAtom = {self.ma*const.Av} kg/mol\n'
        outs += f'   Vin        = {self.Vin} [V]\n'
        outs += f'   L          = {self.L} [m] \n'
        outs += f'   coeffrest  = {self.coeffrest} \n'
        outs += f'   mgkg       = {self.mgkg:0.3e} [kg]\n'
        outs += f'   name       = {self.name}\n'
        outs += f'   mi         = {self.mi:0.3e} [kg]\n'
        outs += f'   mu         = {self.mu:0.3e} [kg]\n'
        outs += f'   ccs        = {self.ccs:0.3e} [m^2]\n'
        outs += f'   natoms     = {self.natoms}\n'
        outs += f'   dof        = {self.dof}\n'
        outs += f'   u/kT       = {self.u/const.kB/self.temp:0.2f} [~dof]\n'
        outs += f'   hc         = {self.hc:0.3e}    [J/K ?]\n'
        outs += f'   z          = {self.z:0.1f} [m]\n'
        outs += f'   Tg(z)      = {self.Tg(self.z):0.1f} [K]\n'
        outs += f'   n(z)       = {self.n(self.z):0.3e} [#/m^3]\n'
        outs += f'   mfp(z)     = {self.mfp(self.z):0.6f} [m]\n'
        outs += f'   ncoll0     = {self.ncoll0:0.1f} \n'
        outs += f'   driftfield = {self.driftfield:0.2f} [V/m]\n'
        outs += f'   mobility   = {self.mobility:0.2f} [m/s / (V/m)]\n'
        outs += f'   driftvel   = {self.driftvel:0.2f} [m/s]\n'
        outs += f'   thermvel   = {self.thermvel:0.2f} [m/s]\n'
        outs += f'   vxyz       = {self.vx:0.2f}, {self.vy:0.2f}, {self.vz:0.2f} [m/s]\n'

        return outs

    def start(self, z=0.0, dVin=0):
        '''Initial setup of the model variables before starting steptime().
        This could be used to reinitialize the model before starting a new run.'''
        self.natoms = self.ion.natoms
        self.dof = 3 * self.natoms - 6          # number of degrees of vibrational freedom in the ion
        self.collided = False

        # init starting variable values.
        self.T0 = self.ion.T0
        self.umax = 0
        self.Tmax = 0
        self.z = z
        self.vz = 0
        self.t = 0
        self.ts = 0

        self.zlen = self.cell.L       # run from z=0 to zlen as minimum defined field length

        # cell variables -- direct translate from main() above
        self.L = self.cell.L
        self.mgkg = self.cell.Mg
        # self.pa, self.pb, self.pc = self.cell.pseudomass_params()

        # self.pressure = self.cell.Pg  # only used for number density and mobility scaling
        self.Tg = self.cell.T                                       # Tg(z)

        # assign variables for given ion (protein)
        self.name = self.ion.name
        self.mi = self.ion.mi       # mass in kg
        self.Z = self.ion.Z         # charge in C
        self.ccs = self.ion.ccs     # ccs in m^2

        self.hc = self.ion.UfromT(self.T0)*self.dof/const.Av
        # 32 Da pseudoatom mass is from Uggerud/Derrick, good for small molecules,
        # self.mi/self.natoms*(self.dof*const.kB/self.hc)
        self.mu = (self.mi*self.mgkg)/(self.mi+self.mgkg)      # reduced mass of protein and gas
        self.n = lambda z: self.cell.rho(z)/self.cell.Mg  # now can do n(z)
        self.mfp = lambda z: 1/(self.ccs * self.cell.rho(z)/self.cell.Mg)

        self.ncoll0 = self.L/self.mfp(self.z)      # predicted number of collisions
        self.driftfield = self.cell.Ez(self.z, self.t)
        self.mobility = 3*self.Z/(16*self.n(self.z))*math.sqrt(
            2*np.pi/(self.mu*const.kB*self.cell.T(self.z)))*1/self.ccs  # mobility, K, at Tg in m^2/V*s
        self.driftvel = self.driftfield*self.mobility
        self.thermvel = math.sqrt(8*const.kB*self.Tg(self.z)/(np.pi*self.mi))  # mean Maxwell-Boltzmann speed distribution

        # explore possible time steps
        # for initial data set timesteps may be from 1e-11 to 1e-8 s
        vel_currentmax = max(self.thermvel, self.cell.vgz(self.z))
        self.ts = 0.05 * self.mfp(self.z) / vel_currentmax      # thermvel at entrance point
        self.ts = min(1e-8, self.ts)
        self.vx = 0
        self.vy = 0
        self.KE = dVin * self.Z                                 # start param dVin in V to KE in J
        self.vz = math.sqrt(2 * self.KE / self.mi)# self.cell.vgz(self.z) + math.sqrt(2 * self.KE / self.mi)
        self.vn = math.sqrt(self.vz**2+self.vx**2+self.vy**2)

        self.t = 0
        self.thermolimitT = 1e5
        self.uatthermolim = 1e10
        self.uatTg = 0

        self.temp = self.T0  # in K
        self.Tmax = self.temp
        self.u0 = self.ion.UfromT(self.T0)
        self.u = self.u0
        self.umax = self.u
        self.uatTg = self.ion.UfromT(self.Tg(z))

        self.tcount = 0
        self.collisioncount = 0
        self.frac = 1
        self.fragmented = False
        self.dfrac = 0.0

    def hcT(self, temperature):
        '''An internal convenience function to calculate heat capacity.'''
        hc = self.ion.UfromT(temperature) * self.dof/const.Av
        return hc

    def steptime(self):
        '''Step the model through a single time step.'''
        self.collided = False       # did it collide in this time step?
        gasflow_vz = self.cell.vgz(self.z)  # 0 by default
        vn_mfp = np.sqrt(self.vx**2 + self.vy**2 + (self.vz - gasflow_vz)**2)
        self.ma = self.cell.pseudoatom_mass(self.temp,vn_mfp)/1000/const.Av #added this to update as a variable pseudo-atom mass. could only fail at unrealistically large velocity for these cells
        # note that vaz, vax, vay are in the rest frame of the ion
        vgz, vgx, vgy, vaz, vax, vay = self.cell.randomgasvel(self.z, self.temp, self.ma)
        
        vgz += gasflow_vz           # now vgx,vgy,vgz are gas velocity in lab frame, including bulk gas flow

        # need these first
        self.vrel = np.sqrt((self.vx-vgx)**2 + (self.vy-vgy)**2 + (self.vz-vgz)**2)
        self.vn = np.sqrt(self.vx**2 + self.vy**2 + self.vz**2)
        # self.ma = self.cell.pseudoatom_mass(self.temp,self.vn)/1000/const.Av #added this to update as a variable pseudo-atom mass. could only fail at unrealistically large velocity for these cells
        vgmode = np.sqrt(2*const.kB*self.Tg(self.z)/self.mgkg)
        vgmean = np.sqrt(4/np.pi)*vgmode
        #vgmode += gasflow_vz
        #vgmean += gasflow_vz
        
        s = vn_mfp/vgmode
        meanvrel = vgmean*( (s+1/(2*s)) * np.sqrt(np.pi)/2*special.erf(s) + 1/2*np.exp(-s**2) )
        mfp = vn_mfp/(self.n(self.z) * self.ccs * meanvrel)        # rigorously correct effective mfp
        self.ts = min(1e-8, mfp/vn_mfp*0.05)       # will be 10 ns or smaller as needed.

        # calculate the distance that the ion travels in the time step, ts
        dtr = vn_mfp * self.ts   # total distance traveled in time step
        dz = self.vz * self.ts    # distance traveled in z direction during time step

        # compute collison rate and determine collision probability
        cdf = 1-np.exp(-dtr/mfp)
        coll = random.uniform(0, 1)

        # calculate acceleration from electric field
        az = self.Z/self.mi * self.cell.Ez(self.z, self.t)
        # changevel = 0
        if coll < cdf:
            # model a collision
            self.collisioncount += 1
            self.collided = True

            # COLLISION

            # compute random line of centers vector at collision, (x_ion-math.sin(math.acos(zcoll))x_gas)
            # zcoll = math.cos(math.asin(np.sqrt(0.999999999*random.uniform(0, 1))))  # zcoll is defined parallel to vrel
            # now we just pick a random x- and y- coordinate in the lab frame, consistent with the length of zcoll, and rotate it around the cross product of vrel and the lab z-axis
            # phicoll = random.uniform(0,2*np.pi)
            # msmasz = math.sin(math.acos(zcoll))
            # xcoll = math.cos(phicoll) * msmasz
            # ycoll = math.sin(phicoll) * msmasz
            # now we need a unit vector for the axis of rotation, which is just vrel x lab-z-axis (normalized)
            # thetavrel = math.acos((self.vz-vgz)/self.vrel)
            # vrelxylen = np.sqrt((self.vy-vgy)**2+(self.vx-vgx)**2)
            # axrotx = (self.vy-vgy)/vrelxylen
            # axroty = -(self.vx-vgx)/vrelxylen
            # note that axrotz would be 0 by design, and now (axrotx,axroty,0) is a unit vector
            # now we rotate (xcoll,ycoll,zcoll) about this axis clockwise by the angle (theta) between vrel and the lab z-axis
            # to do this, we use Rodrigues's formula: vrot = v cos (-theta) + k x v (sin (-theta)) + k (k dot v) (1-cos (-theta))
            # cthetavrel = math.cos(-thetavrel)
            # sthetavrel = math.sin(-thetavrel)
            # kdotv = axrotx*xcoll + axroty*ycoll #unit rotation axis dotted with vrel
#            xcollrot = xcoll * cthetavrel + sthetavrel * zcoll * axroty + axrotx * kdotv * (1-cthetavrel)
#            ycollrot = ycoll * cthetavrel - sthetavrel * zcoll * axrotx + axroty * kdotv * (1-cthetavrel)
#            zcollrot = zcoll * cthetavrel + sthetavrel * (axrotx*ycoll - axroty*xcoll)
            # since xcoll, ycoll, zcoll is a unit vector, so is its rotated version

            # where did the collision occur?
            y = coll*((1/mfp) - (1/mfp)*np.exp(-dtr/mfp)) + (1/mfp)*np.exp(-dtr/mfp)   # JP: this is the only part of the code I don't recall how to derive
            dcoll = -mfp*np.log(mfp*y)  # how far in absolute distance to the collision, always +
            dzcoll = dcoll*dz/dtr       # what that translates to in terms of z direction, can be + or -, depending on sign of dz

            # self.z += dzcoll            # this can move d backwards, if vz is negative, as it should

            # tnewstep = abs(dcoll/self.vn)
            tnewstep = abs(dcoll/vn_mfp)

            self.t += tnewstep

            # calculate change in KE in z direction due to cell potential
            # vznew = self.vz + az*tnewstep
            self.z += self.vz * tnewstep + 0.5 * az * tnewstep*tnewstep     # use old vz here and acceleration and keep in LAB frame #incrementing position twice
            # vzflow = vznew
            self.vn = np.sqrt(self.vz**2+self.vx**2+self.vy**2)

            # now calculate post-collision lab-frame velocity components of the ion
#            vrelx = self.vx-vgx
#            vrely = self.vy-vgy
#            vrelz = self.vz-vgz

            # c = 1/2+1/2*np.sqrt(1-(1-self.coeffrest**2)/nreldotncoll2)
            c = 1.0
            # self.ma = self.cell.pseudoatom_mass(self.temp,abs(self.vn - gasflow_vz))/1000/const.Av #added this to update as a variable pseudo-atom mass. could only fail at unrealistically large velocity for these cells
            dvelscale = 2*c*self.mgkg/(self.mgkg+self.ma)
            vznew = self.vz + self.ma/self.mi*(vaz-dvelscale*(vaz-(vgz-self.vz)))  # note that vzflow here is includes any field acceleration
            vxnew = self.vx + self.ma/self.mi*(vax-dvelscale*(vax-(vgx-self.vx)))
            vynew = self.vy + self.ma/self.mi*(vay-dvelscale*(vay-(vgy-self.vy)))
            
            # changevel += self.ma/self.mi*(vaz-dvelscale*(vaz-(vgz-self.vz)))
            
            # print(changevel)
            # print(fmt.format(ua/const.qe, ub/const.qe, uc/const.qe, ud/const.qe, ue/const.qe, ut/const.qe))
            # This now directly calculates total internal energy at this timestep

            # mg = self.cell.Mg

            self.u = self.u + 1/2*self.ma*(
                    -4*self.mgkg*c/(self.ma+self.mgkg)*(vax*(vax-(vgx-self.vx))+vay*(vay-(vgy-self.vy))+vaz*(vaz-(vgz-self.vz)))
                    + 4*self.mgkg**2*c**2/(self.ma+self.mgkg)**2*((vax-(vgx-self.vx))**2+(vay-(vgy-self.vy))**2+(vaz-(vgz-self.vz))**2))

            self.vx = vxnew
            self.vy = vynew
            vzflow = vznew
            vnew = np.sqrt(vxnew**2+vynew**2+vznew**2)

            if self.u > self.umax:
                self.umax = self.u
            self.temp = self.ion.TfromU(self.u)
            self.Tmax = self.ion.TfromU(self.umax)
            self.hc = self.hcT(self.temp)  # defined in J/mol*K
            self.KE = 0.5*self.mi*vnew**2

        else:
            # model movement with no collision
            # no collision, update velocities due to cell potential and advance
            self.t += self.ts
            self.z += self.vz*self.ts + 0.5*az*self.ts**2         # movement includes gasflow offset velocity

            vznew = self.vz + az*self.ts
            # print(az)
            vzflow = vznew
            self.vn = np.sqrt(self.vz**2+self.vx**2+self.vy**2)

            self.KE = 0.5*self.mi*self.vn**2
            self.temp = self.ion.TfromU(self.u)
            self.hc = self.hcT(self.temp)

        thisloss = fracloss(self.ts, self.temp, self.ion.dH, self.ion.dS)
        self.dfrac = thisloss / self.ts    # a dissociation rate units 1/s

        cdfrac = 1-np.exp(-thisloss)
        fracfrac = random.uniform(0, 1)
        self.fragmented = self.fragmented or (fracfrac < cdfrac)
        self.frac = self.frac * (1-thisloss)
        #print(f'dfrac: {self.dfrac}   cdfrac: {cdfrac}   fracfrac: {fracfrac}   fragd: {self.fragmented}   frac: {self.frac}')

#        print(f'{self.t*1e6:5.2f} us  {self.z*1e3:5.1f} mm   T: {self.temp:6.1f}   thisloss: {thisloss:0.2e}   frac: {1-self.frac:0.4f}')
        self.tcount += 1
        # print(vzflow)
        self.vz = vzflow

    def run_model(self, Vin, plotresult=False, verbose=False):
        '''With default input parameters, does a full run, returning
        a reference to the model and a combined tT (time, Temperature) array.
        '''
        if plotresult:
            import matplotlib.pyplot as plt

        self.start(z=0, dVin=Vin)

        if verbose > 1:
            print(self.headerstr())
            print(self)

        t = [self.t]        # put starting values into lists
        T = [self.temp]
        f = [self.frac]
        ma_list=[self.ma]
#        u = [self.u/const.qe]
        position = [self.z]
#        KE = [self.KE]
        df = [self.dfrac]
        frag = [self.fragmented]

    #    print(m.headerstr()); print(m)
        while self.z < self.cell.L and self.z >= 0 and self.tcount < 1e7:
            self.steptime()
            if verbose > 1:
                if self.tcount < 2 or (self.collided and self.collisioncount <= 3):
                    print(self)
            if verbose > 2:
                if self.tcount % 5000 == 0:
                    print(self)
            if self.collided or self.tcount < 4:      # only add to array if ion just collided
                t.append(self.t)
                T.append(self.temp)
                f.append(self.frac)
                ma_list.append(self.ma)
#                u.append(self.u/const.qe)
                position.append(self.z)
#                KE.append(self.KE)
                df.append(self.dfrac)
                frag.append(self.fragmented)

        t.append(self.t)        # add final step to history also
        T.append(self.temp)
        f.append(self.frac)
        position.append(self.z)
        df.append(self.dfrac)
        frag.append(self.fragmented)

        tT = np.array([t, T])
        avevdrift = self.cell.L/t[-1]
        exptemp = self.Tg(self.z)+1/3*self.mgkg*avevdrift**2/const.kB
        avgfinaltemp = sum(T[-1000:])/1000

        if verbose:
            print(self)     # a final line for the single line prints
            print('{} timesteps, {} collisions in {:0.2f} ms     umax = {:0.3f} eV    Tmax = {:0.2f} K   T = {:0.2f}'.format(
                self.tcount, self.collisioncount, self.t*1000, self.umax/const.qe, self.Tmax, self.temp))
            print(f'average vdrift = {avevdrift:0.3f} m/s, expected ion temperature due to drift = {exptemp:0.1f}, avg ion temp for last 1000 collisions = {avgfinaltemp:0.1f}')

        if plotresult:
            plt.figure(); plt.plot(t, T); plt.title('T vs t'); plt.show()
            # plt.figure(); plt.plot(t, u); plt.title('u vs t'); plt.show()
            plt.figure(); plt.plot(t, position); plt.title('z vs t'); plt.show()
            # plt.figure(); plt.plot(t,KE); plt.title('KE vs t'); plt.show()
            plt.figure(); plt.plot(t, f); plt.title('f vs t')
            plt.plot(t, df)
            plt.plot(t,frag)
            plt.plot(t, [Tx/1000 for Tx in T])
            plt.show()

        return(self, tT)


def testWAVE(Tion=298, Vin=1, **args):
    '''A test function for the WAVE cell.
    Can test either the simple WAVE cell or the pressure-varying cell.
    Change the args dict to modify cell parameters.'''
    rho = args.get('rho', 4.179e-5)
    Tcell = args.get('T', 298)

    gamu = 40   # collision gas mass
    p = rho / (gamu/1000) * const.R * Tcell
    pmtorr = p * 7.5006167382113

    cellL = 0.51
    # cellrho = rho

    # cell = WAVEcell(L=0.125, Mgamu=40, T=298, rho=4.19e-5) #normally rho is around 4e-5 for Trap
    # cell = WAVEcell(L=0.254, wh=16, vw=500, Mgamu=28, T=298, rho=2.82/0.5*0.000565) #use these for TWIMS, rho of 0.000565 is 0.5 mbar N2
    args['L'] = cellL
    cell = makecell(args)  # normally rho is around 4e-5 for Trap

#    iond = dict(name='CC322', mass=0.322, charge=1, CCS=1.5367, num_atoms=37, T0=298, hcprofname='tunemix322')
#    iond = dict(name='cc1222', mass=1.222, charge=1, CCS=2.8125, num_atoms=91, T0=Tion, hcprofname='tunemix1222')
#    iond = dict(name='myoglobin', mass=17.585, charge=9, CCS=20.2, num_atoms=2411, T0=Tion, hcprofname='peptide')
    iond = {
    	"name": "Stx11",
    	"mass": 39.111,
    	"charge": 11,
    	"CCS": 31.4,
    	"num_atoms": 5495,
    	"T0": 300,
        "hcprofname": "peptide",
        "refdH": 267,
        "refdS": 326
    }
    ion = ionclass(iond)

    m = modelclass(cell, ion, Vin)

    expl = f'''Test run with WAVE cell ({cellL} meter)
    and input energy {Vin:0.1f}.  The initial ion temperature is {Tion} K and the cell
    temperature is {Tcell} K at a pressure of {pmtorr:0.2f} mTorr.
    The internal ion temperature would be expected to also settle
    a little higher than the initial temperature.
    Expected drift velocity is {m.driftvel:0.3f} m/s.
    The ion is myoglobin.
    '''
    print(expl)

    print(repr(m))
#    mm, tT = m.run_model(Vin, plotresult=True, verbose=3)
    m.start(z=0, dVin=Vin)

    print(m.headerstr())
    print(m)

    t = [m.t]        # put starting values into lists
    T = [m.temp]
    position = [m.z]
    KE = [m.KE/ion.Z]
    rhos = [m.cell.rho(m.z)]
    Es = [m.cell.Ez(m.z, m.t)]


#    print(m.headerstr()); print(m)
    while m.z < m.cell.L and m.z >= 0 and m.tcount < 5e5:
        m.steptime()
        if m.tcount < 2 or (m.collided and m.collisioncount <= 3):
            print(m)
        if m.tcount % 5000 == 0:
                print(m)
        if m.collided or m.tcount < 4 or m.tcount % 100 == 0:      # only add to array if ion just collided
            t.append(m.t)
            T.append(m.temp)
            position.append(m.z)
            KE.append(m.KE/ion.Z)
            rhos.append(m.cell.rho(m.z))
            Es.append(m.cell.Ez(m.z, m.t))


    t.append(m.t)        # add final step to history also
    T.append(m.temp)
    position.append(m.z)
    KE.append(m.KE/ion.Z)
    rhos.append(m.cell.rho(m.z))
    Es.append(m.cell.Ez(m.z, m.t))

    tT = np.array([t, T])
    avevdrift = m.cell.L/t[-1]
    exptemp = m.Tg(m.z)+1/3*m.mgkg*avevdrift**2/const.kB
    avgfinaltemp = sum(T[-1000:])/1000

    print(m)     # a final line for the single line prints
    print('{} timesteps, {} collisions in {:0.2f} ms     umax = {:0.3f} eV    Tmax = {:0.2f} K   T = {:0.2f}'.format(
            m.tcount, m.collisioncount, m.t*1000, m.umax/const.qe, m.Tmax, m.temp))
    print(f'average vdrift = {avevdrift:0.3f} m/s, expected ion temperature due to drift = {exptemp:0.1f}, avg ion temp for last 1000 collisions = {avgfinaltemp:0.1f}')

    import matplotlib.pyplot as plt
    plt.figure(); plt.plot(t, T); plt.title('T vs t'); plt.show()
    plt.figure(); plt.plot(t, position); plt.title('z vs t'); plt.show()
    plt.figure(); plt.plot(t,KE); plt.title('KE/Z vs t'); plt.show()
    plt.figure(); plt.plot(position,KE); plt.title('KE/Z vs z'); plt.show()
    plt.figure(); plt.plot(position,rhos); plt.title('rho vs z'); plt.show()
    plt.figure(); plt.plot(position,Es); plt.title('E vs z'); plt.show()




def testCID(Tion=298, Tcell=298, Vin=1, L=0.18, dV=7, rho=2.42e-5):
    '''A test function for the simple collision cell.'''
    gamu = 28   # collision gas mass
    p = rho / (gamu/1000) * const.R * Tcell                 # p in Pa
    pmtorr = p * 7.5006167382113    # from pressure in Pa to mTorr

    cellL = L
    cellV = dV
    cellrho = rho

    cell = cellclass(L=cellL, Mgamu=gamu, T=Tcell, rho=cellrho, E=(cellV/cellL), z0=0)
    # cell = WAVEcell(L=0.125, Mgamu=40, T=298, rho=4.19e-5) #normally rho is around 4e-5 for Trap
    # cell = WAVEcell(L=0.254, wh=16, vw=500, Mgamu=28, T=298, rho=2.82/0.5*0.000565) #use these for TWIMS, rho of 0.000565 is 0.5 mbar N2
#    iond = dict(name='cc1222', mass=1.222, charge=1, CCS=2.8125, num_atoms=91, T0=Tion, hcprofname='tunemix1222')
    iond = dict(name='myoglobin', mass=17.585, charge=9, CCS=20.2, num_atoms=2411, T0=Tion, hcprofname='peptide', dH=80, dS=-100)
    ion = ionclass(iond)
#    ion = ionclass(dict(name='CC322', mass=0.322, charge=1, CCS=1.5367, num_atoms=37, T0=298, hcprofname='tunemix322'))
    m = modelclass(cell, ion, Vin)

    expl = f'''Test run with CID cell ({L} meter) with axial field of {dV/L:0.2f}
    and input energy {Vin:0.1f}.  The initial ion temperature is {Tion} K and the cell
    temperature is {Tcell} K at a pressure of {pmtorr:0.2f} mTorr.
    The internal ion temperature would be expected to also settle
    a little higher than the initial temperature.
    Expected drift velocity is {m.driftvel:0.3f} m/s.
    The ion is myoglobin.
    '''
    print(expl)

    print(repr(m))
    mm, tT = m.run_model(Vin, plotresult=True, verbose=3)


if __name__ == '__main__':
    '''Test code to validate basic function when run standalone.'''

#    testCID(Tion=300, Tcell=300, Vin=65.0, rho=2.42e-5)  # myoglobin with standard CC at 16.1 mTorr
#    testCID(Tion=300, Tcell=300, Vin=0.1, L=0.02, dV=0.7, rho=2.42e-4)  # high pressure test run testing low field limit
#    testCID(Tion=300, Tcell=300, Vin=0.1, rho=2.42e-9)  # very low (near zero) pressure for no collisions
    wd = {
    "type": "wcell3",
    "L": 0.51,
    "wl": 0.0121,
    "wh": 4.0,
    "vw": 311.0,
    "Mgamu": 40,
    "T": 298.0,
    "rho": 3.213e-5,
    "rho1": 3.213e-7,
    "zopen": 0.09,
    "zclose": 0.42,
    "zV1": 0.045,
    "dV1": -1,
    "zV2": 0.09,
    "dV2": -2,
    "zV3": 0.13,
    "dV3": -1,
    "zV4": 0.38,
    "dV4": 0,
    "zV5": 0.42,
    "dV5": 1.0,               # 3 from instrument actuals, needs to be lower to improve transmission
    "zwaveoff": 0.13,
    "zwaveon": 0.38,
    "z0": 0.0,
    "T0offset": 0.0
    }
    testWAVE(Tion=300, Vin=1.0, **wd)

    exit(0)     # skip all code after this

    # test comparison code currently disabled.
    import sys
    args = sys.argv
    Vin = 5
    runWaters = False

    if 'CID' in args:
        runWaters = False
    else:
        runWaters = True

    if len(args) > 1:
        for arg in args:  # this uses the last number in command line as Vin
            try:
                Vin = float(arg)
            except:
                pass

    # two example cell types
    cellCID = cellclass(L=0.18, Mgamu=28, T=298, rho=2.42e-5, E=(7/0.18), z0=0)       # default (Agilent 6560) CID cell L = 0.18
    # cellCID = cellclass(L=0.18, Mgamu=28, T=298, rho=2.42e-5, E=(7/0.18), z0=0)       # default (all other Agilent QTOFs) CID cell L = 0.18
#    cellWAVE = WAVEcell()           # default Waters WAVE cell
    cellWAVE = WAVEcell(Mgamu=40, T=298, rho=4.18e-5)

    runWaters = False

    if runWaters:
        cell = cellWAVE
        cell1 = WAVEcell()
    else:
        cell = cellCID

    # short sample with this ion type
    uion = dict(name='myoglobin', mass=17.585, charge=9, CCS=20.2, num_atoms=2411, T0=298, hcprofname='peptide')
    ion = ionclass(uion)

    m = modelclass(cell, ion, Vin)
    m1 = modelclass(cell1, ion, Vin)

    print(m)
    mm, tT = m.run_model(Vin, plotresult=False, verbose=True)

    print(m1)
    mm1, tT = m1.run_model(Vin, plotresult=False, verbose=True)
