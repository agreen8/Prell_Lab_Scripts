# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 12:00:53 2022

@author: jprell
Copyright (c) 2024 James S. Prell, all rights reserved
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy as sp
#from scipy.optimize import curve_fit
from ionspa import loadjfile, makecell, ionclass, const
#import cell         # contains cell support code including pseudoatom gas parameters
import argparse

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

#data from Bush CCS database for native-like proteins in N2 buffer gas
Bushzs = [6,7,7,8,11,12,13,14,15,16,15,16,17,14,15,16,17,19,20,21,22,22,23,24,25,26,23,24,25,26,31,32,33,34,35,31,32,33,34,37,38,39,40,41,65,66,67,68,69,70,71]
Bushms = [12,12,18,18,37,37,37,56,56,56,64,64,64,66,66,66,66,103,103,103,103,125,125,125,125,125,143,143,143,143,237,237,237,237,237,250,250,250,250,336,336,336,336,336,801,801,801,801,801,801,801]
BushN2CCSs = [14.9,15.9,19.5,20.3,32.3,33.1,34.3,38.4,38.5,38.8,41.5,41.5,41.6,44.9,44.9,44.7,44.9,60.6,60.8,60.9,60.5,76.3,76.0,74.6,73.1,72.8,74.2,74.5,74.4,75.0,111,111,110,110,110,111,112,112,111,134,134,134,134,135,218,220,220,219,219,218,219]

kB = const.kB       # Boltzman
qe = const.qe       # electron charge
NA = const.Av       # Avogadro
mamu = const.mamu   # mass in kg for unit mass

#Ar = {"ma0K":68.56,"maslope":0.0222,"mg":39.948} #units are Da, Da/K, Da, values from AMBER99 GROMACS collision simulations
#N2 = {"ma0K":62.486,"maslope":0.0169,"mg":28.02}
#He = {"ma0K":32.592,"maslope":-0.007,"mg":4.0026}

defcell = {'type': 'cell', 'L': 0.18, 'Mgamu': 28, 'T': 298.15, 'E': 0.03889, 'rho': 2.712e-05, 'z0': 0.0, 'Voffset': 0.0, 'T0offset': 100}
defion = {'name': 'cytc7', 'mass': 17.594, 'charge': 7, 'CCS': 15.9, 'num_atoms': 1694, 'T0': 298.15, 'hcprofname': 'peptide', 'refdH': 95, 'refdS': -100} 
defV = 15

vgmode = 0
vgmean = 0

class vlogger():
    '''My logger class.
    Initialize as mylog=vlogger(x=1, y=2)
    Append values with mylog.append(x=2, y=3) # same names
    Access logged lists as mylog.x or mylog.y
    -- general assistance to store variables in a loop for plotting, etc.
    '''
    def __init__(self, **parms):
        # make named list for each parm
        # and initialize with first value
        self.pdict = {}
        for pname in parms:
            if parms[pname] is None:
                self.pdict[pname] = []      # allow initialization with empty list
            else:
                self.pdict[pname] = [ parms[pname] ]
            setattr(self, pname, self.pdict[pname])

    def append(self, **parms):
        '''Append parms values to lists'''
        for pname in parms:
            if pname in self.pdict:
                self.pdict[pname].append( parms[pname] )

    def vnames(self):
        '''Return list of parameter names'''
        return [name for name in self.pdict]


def main(**args):
    """Main calculation function.
    Uses args['cell'], args['ion'], injectV, option, sizefactor, Ngfactor.
    May use cellfilename or ionfilename or possibly exp/Vin params in future.
    May add other params as needed for plot selection, etc.
    Returns a logger object with values at each collision.
    """
    print('in main():')

    for arg in args:
        print(f'{arg}: {args[arg]}')
    print()

    # extract input parameters
    celld = args.get('cell', defcell)
    iond = args.get('ion', defion)
    injectV = args.get('injectV', defV)
    if not type(injectV) == list: injectV = [injectV]   # use this as a list even for single value
    option = args.get('option', 'IICT')
    sizefactor = args.get('sizefactor', 1.0)
    Ngfactor = args.get('Ngfactor', 1.0)
    xaxis = args.get('xaxis', 'coll')
    yaxis = args.get('yaxis', 'Ti')
    plotfile = args.get('plotfile', 'screen')

    # make cell and ion objects so we can use methods from those classes
    ccell = makecell(celld)
    cion = ionclass(iond)

    # Get values we use from the cell definition
    cellVoff = celld.get('Voffset', 0)

    rho_g = celld.get('rho', 2.712e-5)
    Lcell = celld.get('L', 0.18)            # length of cell, m, Agilent 45XT and 6560 default is 0.18 m
    Ecell = celld.get('E', 38.89)
    Tg = ccell.T(0)         # initial gas temperature, usually constant
    mgamu = ccell.Mg / mamu # gas mass in Daltons
    mg = ccell.Mg

    print('celld:', celld)
    print(f'cell class:\n  {ccell}')
    print(f'Cellparms\n   Mgamu: {mgamu}\n   Lcell: {Lcell}\n   Ecell: {Ecell}\n   Tcell: {Tg}\n   rhog: {rho_g}\n   Voff: {cellVoff}\n')

    # Get values we use from the ion definition
    Mion = iond.get('mass', 17.594)      # 17.594 kg/mol for myo9, 
    z = cion.Z/qe                        # ion charge state
    numatoms = cion.natoms               # number atoms in the ion
    dof = numatoms*3 - 6                 # degrees of freedom
    ccs = cion.ccs                       # collision cross section in m^2
    Tion0 = cion.T0                         # initial ion temperature
    hcprofname = cion.hcprofname
    dof = 3 * numatoms - 6
#    cion.hcprofileobject.buildfunctions(iond['hcprofname'], dof=dof)

    print('iond: ', iond)
    print(f'ion class:\n  {cion}')
    print(f'ion parms\n   Miamu: {Mion*1000}\n   z: {z}\n   natoms: {numatoms}\n   CCS: {ccs/1e-18}\n   Ti: {Tion0}\n   hctype: {hcprofname}')

    injectV = [iV + cellVoff for iV in injectV]
    print(f'injectV offset by {cellVoff} to {injectV} using cell[Voffset]\n')

    print(f'Calculation option is {option}')
    print(f'Size scaling is {sizefactor}')
    print(f'Density scaling is {Ngfactor}\n')

    # Some initial calculations not expected to change during a run
    vgmode = np.sqrt(2*kB*Tg/mg)
    vgmean = np.sqrt(4/np.pi)*vgmode
    print(f'vgmode: {vgmode:0.1f} m/s      vgmean: {vgmean:0.1f} m/s')

    KE = z*qe*injectV[0]       # ion kinetic energy, in J
    size = sizefactor
    mi = Mion/NA*size
    Ng = Ngfactor * rho_g/mg      # gas number density, 1/m^3, 2.712e-5 = Collision Cell normal pressure


    def meanfreet(ng,kinen): #see SIMION documentation from David Manura for derivation
        '''In addition to input params, uses mi, vgmode(Tg,mg), vgmean(vgmode), ccs.
        Defined inside main() to pick up those values. If Tg changes with position, 
        would need to also recalculate vgmode, vgmean and meanfreet'''
        vn = np.sqrt(2*kinen/mi)
        s = vn/vgmode
        meanvrel = vgmean*((s+1/(2*s))*np.sqrt(np.pi)/2*sp.special.erf(s)+1/2*np.exp(-s**2))
        freet = 1/(ng * ccs * meanvrel)
        return freet

    print(f'Initial meanfreet: {meanfreet(Ng, KE):0.2e}\n')

    Vlogger = vlogger(injV=None, KE=None, Ui=None, T=None)

    if plotfile != 'none':
        fig, ax = plt.subplots(figsize=(7, 5))
    postdictlist = []
    for injectv in injectV:
        postd = dict(injectv=injectv)
        Ti0 = Tion0
        ui0 = cion.UfromT(Ti0)
        
        KE = z*qe*injectv
        #KE = 3/2*kB*(Ti*mg+Tg*(ma0K+maslope*Ti))/((ma0K+maslope*Ti)+mg-(ma0K+maslope*Ti)*mg/mi)

        ui = ui0
        Ti = Ti0

        veli0 = np.sqrt(2*KE/mi)  # initial ion speed in m/s
        ma0 = ccell.pseudoatom_mass(Ti0,veli0)/NA/1000
        chi0 = 4*ma0*mg/(ma0+mg)**2

        time = 0
        collision = 0
        zpos = 0.0

        Vlogger.append(injV=injectv, KE=KE/qe, Ui=ui0/qe, T=Ti)
        tlogger = vlogger(coll=0, t=time, zp=0.0, KE=KE/qe, ui=ui0/qe, Ti=Ti0, chi=chi0)

        while (zpos < Lcell):           # and collision < 2E7):
            Ecell = ccell.Ez(zpos, t=0)
            Ng = Ngfactor * ccell.rho(zpos)/mg      # allow for density change within the cell
            collision += 1
            deltat = meanfreet(Ng,KE)
            time += deltat
            zpos = zpos + np.sqrt(2*KE/mi)*deltat
            veli = np.sqrt(2*KE/mi)
            ma = ccell.pseudoatom_mass(Ti,veli)/NA/1000
            # ma = mg #un-comment if you want to override ma and force chi = 1
            chi = 4*ma*mg/(ma+mg)**2
            #chi = 4*(68.56 + 0.0222*Ti)/1000/NA*mg/((68.56 + 0.0222*Ti)/1000/NA+mg)**2
            #chi=chi0

            if option == 'xa': #omit option a!
                # old IICT model (had typo in KE formula)
                KEn = KE + z*qe*Ecell*np.sqrt(KE*2/mi)*deltat + chi/(2*mi)*(3*kB*Ti*mg+3*kB*Tg*ma) + chi/mi*(-(ma+mg)+mg*ma/mi)*KE #likely incorrect (from old derivation)
                uin = ui + chi/2*3*kB*(Tg-Ti)+chi*mg/mi*KE #IICM                
            elif option == 'IICT': #call this 'Improved Impulsive Collision Theory' ('IICT' for short)
                # correct IICT model        
                KEn = KE + z*qe*Ecell*np.sqrt(KE*2/mi)*deltat + 3/2*kB*Ti*ma/mi*(1-chi) + 3/2*ma/mi*chi*kB*Tg + chi/mi*(ma*mg/mi-ma-mg)*KE #correct version for IICM
                uin = ui + chi/2*3*kB*(Tg-Ti)+chi*mg/mi*KE #IICM
            elif option == 'ELASTIC': #break these into two options, with the same uin, call the first one 'Elastic' and the second 'Elastic Head-On' ('Elastic HO' for short)
                # fully elastic model        
                KEn = KE + z*qe*Ecell*np.sqrt(KE*2/mi)*deltat - 1/3*4*mi*mg/(mi+mg)**2*KE + 1/2*4*mi*mg/(mi+mg)**2*kB*Tg #fully elastic model with random impact parameter 
                #KEn = KE + z*qe*Ecell*np.sqrt(KE*2/mi)*deltat - 4*mi*mg/(mi+mg)**2*KE + 3/2*4*mi*mg/(mi+mg)**2*kB*Tg #fully elastic model with only head-on collisions
                uin = ui0
            elif option == 'ELASTICHO': #break these into two options, with the same uin, call the first one 'Elastic' and the second 'Elastic Head-On' ('Elastic HO' for short)
                # fully elastic model        
                #KEn = KE + z*qe*Ecell*np.sqrt(KE*2/mi)*deltat - 1/3*4*mi*mg/(mi+mg)**2*KE + 1/2*4*mi*mg/(mi+mg)**2*kB*Tg #fully elastic model with random impact parameter 
                KEn = KE + z*qe*Ecell*np.sqrt(KE*2/mi)*deltat - 4*mi*mg/(mi+mg)**2*KE + 3/2*4*mi*mg/(mi+mg)**2*kB*Tg #fully elastic model with only head-on collisions
                uin = ui0
            elif option == 'INELTH': #call this 'Inelastic Thermalizing' ('Inel. Therm.' for short)
                # fully inelastic model
                hc = (cion.UfromT(Ti+1)-cion.UfromT(Ti)) # for calculations that invoke hc explicitly
                KEn = KE + z*qe*Ecell*np.sqrt(KE*2/mi)*deltat + (-1+mi**2/(mi+mg)**2+3/2*mg**2/(mi+mg)**2*kB/hc)*KE + 3/2*mi*mg/(mi+mg)**2*(1+3/2*kB/hc)*kB*Tg + 3/2*mg/(mi+mg)*kB*Ti # fully inelastic model with thermal re-emission
                uin = ui - 3/2*kB*Ti + mg/(mi+mg)*(1-3/2*kB/hc)*KE + 3/2*mi/(mi+mg)*(1-3/2*kB/hc)*kB*Tg # from fully inelastic model with thermal re-emission

            Tin = cion.TfromU(uin)
            KE = KEn
            ui = uin
            Ti = Tin
            KEeVz = KEn/qe/z
            KEsz = tlogger.KE[-1]/z
            if KEeVz < 1.0 and KEsz > 1.0:
                print(f'KE per charge dropped below 1.0 at {zpos}')

            tlogger.append(coll=collision, t=time, zp=zpos, KE=KE/qe, ui=ui/qe, Ti=Ti, chi=chi)
            ### End of collision step loop

        # this time increment is extraneous as it was already done and this is out of the loop    
#        time += deltat

#        uideltas = []
#        for u in range(1,len(tlogger.ui)):
#            uideltas.append(tlogger.ui[u]-tlogger.ui[0])

        ### section for a multipart plot
        # plt.subplot(1,4,1)
        #plt.plot(times,zs)
        #plt.title('kinetic energy (eV) vs approx distance (m)')
        #plt.xlim([-0.01,0.18])
        #plt.subplot(1,4,1)
        #plt.plot(collisions,uis)
        #plt.title('internal energy (eV) vs collision number')
        #plt.subplot(1,4,2)
        #plt.plot(collisions,KEs)
        #plt.title('kinetic energy (eV) vs collision number')
        #plt.xlim([-10,50000])
        # plt.subplot(1,4,2)
        #plt.plot(times,KEs)
        #plt.title('kinetic energy (eV) vs time (s)')
        #plt.xlim([-0.000005,0.0007])
        # plt.subplot(1,4,3)
        #plt.plot(collisions,KEs)

        if False:
            print('plot options', xaxis, yaxis, plotfile)
            # possible xaxis, yaxis values might be coll, t, zp, KE, ui, Ti, chi
            # add other possible x values options here
            if xaxis == 't':
                xvals = tlogger.t
                pltxlabel = 't [s]'
            elif xaxis == 'coll':
                xvals = tlogger.coll
                pltxlabel = 'collisions'
            elif xaxis == 'zp':
                xvals = tlogger.zp
                pltxlabel = 'z position [m]'

            # add other possible y values options here
            if yaxis == 'Ti':
                yvals = tlogger.Ti
                pltylabel = 'T [K]'
            elif yaxis == 'coll':
                yvals = tlogger.coll
                pltylabel = 'collisions'
            elif yaxis == 'KE':
                yvals = tlogger.KE
                pltylabel = 'KE [eV]'
            elif yaxis == 'ui':
                yvals = tlogger.ui
                pltylabel = 'ui [eV]'
            elif yaxis == 'chi':
                yvals = tlogger.chi
                pltylabel = 'chi'
            elif yaxis == 'zp':
                yvals = tlogger.zp
                pltylabel = 'z position [m]'
            if plotfile == 'none':
                pass
            else:
                plt.plot(xvals, yvals)
                plt.suptitle(f'{cion.name}     option "{option}"')                # main title
                plt.title('{repr(ccell)}', fontsize=8)     # becomes subtitle
                plt.xlabel(pltxlabel)
                plt.ylabel(pltylabel)
                if plotfile == 'screen':
                    plt.show()
                else:
                    plt.savefig(plotfile)
        elif plotfile == 'none':
            pass
#            print(f'plotfile is {plotfile}, no plot done.')
        elif plotfile != 'none':
            ptitle = f"{iond['name']} with {option}"
            plabel = f'V={injectv:0.1f}'
            doplot(ax, tlogger, xaxis, yaxis, plabel, ptitle,  plotfile)


        #plt.xlim([-50,115000])
        #plt.ylim([-5,140])
        #plt.xlim([-0.00001,0.00035])
        #plt.title('internal energy (eV) vs time (s)')
        #plt.xlim([-0.000005,0.0007])
        # plt.subplot(1,4,4)
        #plt.plot(times,Tis)
        #plt.title('internal temperature (K) vs time (s)')
        #plt.xlim([-0.000005,0.0007])
        #plt.show()


        ## Post analysis section
        maxdeltat = max(tlogger.Ti)-Ti0
        maxdeltau = (max(tlogger.ui)-tlogger.ui[0])
        maxreldeltau = maxdeltau/tlogger.ui[0]
        driftvel = np.sqrt(2*(np.sum(tlogger.KE[-101:-1])/100*qe-3/2*kB*Tg)/mi)
        postd['maxdT'] = maxdeltat
        postd['maxdU'] = maxdeltau
        postd['maxreldU'] = maxreldeltau
        postd['vdrift'] = driftvel
#        print(f'maximum deltaT: {maxdeltat:0.02f} K, maximum deltaU: {maxdeltau:0.01f} eV, max deltaU relative to init U: {maxreldeltau:0.01f}, sim drift velocity: {driftvel:0.03f} m/s')
        
        #MSvdrift is the classic Mason-Schamp result
        MSvdrift = 3*z*qe/(16*Ng)*Ecell*np.sqrt(2*np.pi*(mi+mg)/(mi*mg)/(kB*Tg))/ccs
        #MSmodvdrift is the modified Mason-Schamp equation, from solving Eq 35 in Revercomb and Mason, Anal Chem, 1975
        b = 3*kB*Tg/mg
        c = -1/mg*(mg+mi)/(mg*mi)*(z*qe*Ecell/(Ng*ccs))**2
        MSmodvdrift = np.sqrt(-b/2 + 1/2*np.sqrt(b**2-4*c))
#        MSmodvdrifts.append(MSmodvdrift)
        EoverN = Ecell/Ng*1E21 #field strength in Townsend, 1E-21 V/m^2
        maxdelui = max(tlogger.ui)-tlogger.ui[0]
        maxeff = maxdelui/(z*injectv)

        postd['MSvdrift'] = MSvdrift
        postd['MSmodvdrift'] = MSmodvdrift
        postd['EoverN'] = EoverN
        postd['maxdelui'] = maxdelui
        postd['maxeff'] = maxeff
        print(f'E/N is {EoverN:0.2f} [Td]')
#        print(f'MSvdrift: {MSvdrift:0.3f}   MSmodvdrift: {MSmodvdrift:0.3f}   EoverN: {EoverN:0.2f} Td')
#        print(f'maximum delta ui: {maxdelui:0.2f}   maxefficiency = {maxeff:0.4f}')
        
#        thermalv = np.sqrt(3*kB*Tg/mi)
#        thermalvs = []
#        for i in range(len(thermalvs)):
#            thermalvs.append(thermalv)
        #plt.plot(Ngfactor,thermalvs)
    #    ratios = []
    #    for i in range(len(Ngfactor)):
    #        ratios.append(vdrifts[i]/MSmodvdrifts[i]) 
    #    #plt.plot(Ngfactor,ratios)
        
        C=5.0*dof # assumed constant heat capacity
        # parameters from analytical solution work-up
        #old/incorrect parameters (only h is wrong)
        #a=-3*chi0*kB*NA/(2*C)
        #b=chi0*mg/mi
        #c=3/2*chi0*kB*Tg
        #h=3/2*chi0*kB*NA/C*mg/mi
        #k=-chi0/mi*(ma+mg-mg*ma/mi)
        #l=3/2*chi0*kB*Tg*ma/mi
        
        #new/correct parameters
        a = -3*chi0*kB*NA/(2*C)
        b = chi0*mg/mi
        c = 3/2*chi0*kB*Tg
        h = 3/2*(ma/mi)*(1-chi0)*kB*NA/C
        k = -chi0/mi*(ma+mg-mg*ma/mi)
        l = 3/2*chi0*kB*Tg*ma/mi
        
        lamplus=1/2*(a+k+np.sqrt((a-k)**2+4*b*h)) # less negative rate constant w.r.t. collision number, cooling, very nearly = a under many conditions
        lamminus=1/2*(a+k-np.sqrt((a-k)**2+4*b*h)) # more negative rate constant w.r.t. collision number, heating, very nearly = k under many conditions
        cplus = (b*(tlogger.KE[0]-tlogger.KE[-1])+(a-lamminus)*(tlogger.ui[0]-tlogger.ui[-1]))/(lamplus-lamminus)*qe
        cminus = -(b*(tlogger.KE[0]-tlogger.KE[-1])+(a-lamplus)*(tlogger.ui[0]-tlogger.ui[-1]))/(lamplus-lamminus)*qe
        zmax = np.log(-cminus*lamminus/(cplus*lamplus))/(lamplus-lamminus)
#        print(f'zmax = {zmax}')
        umaxapp=(cplus*(-cminus*lamminus/(cplus*lamplus))**(lamplus/(lamplus-lamminus))+cminus*(-cminus*lamminus/(cplus*lamplus))**(lamminus/(lamplus-lamminus)))/qe+tlogger.ui[-1]
#        print(f'umaxapp: {umaxapp}')
        postd['zmax'] = zmax
        postd['umaxapp'] = umaxapp

#        uinfapp=C*Tg/(1-mg/mi)/NA # this will obviously be too high if it uses C at room temperature

#        uapps=[]
#        for Z in range(len(tlogger.coll)):
#            uapps.append((cplus*np.exp(lamplus*Z)+cminus*np.exp(lamminus*Z))/qe+tlogger.ui[0]-(cplus+cminus)/qe)

        #plt.subplot(1,4,1)
        #plt.plot(tlogger.z,uapps)
        fplus = (lamplus-a)*cplus/b
        fminus = (lamminus-a)*cminus/b
        KEapps=[]
        KEinfapp=3/2*kB*Tg*mi/(mi-mg)
        for Z in range(len(tlogger.coll)):
            KEapps.append((fplus*np.exp(lamplus*Z)+fminus*np.exp(lamminus*Z)+KEinfapp)/qe)
        #plt.subplot(1,4,2)
        #plt.plot(collisions,KEapps)
        dd = dict(a=a, b=b, c=c, h=h, k=k, l=l, lamplus=lamplus, lamminus = lamminus, cplus=cplus, cminus=cminus)
        postd.update(dd)
#        print(f'a={a}  b={b}  c={c}  h={h}  k={k}  l={l}  lamplus={lamplus}  lamminus={lamminus}  cplus={cplus}  cminus={cminus}')
        postdictlist.append(postd)
    # section to verify we captured the logging information
    # aVlog = np.array([Vlogger.pdict[n] for n in Vlogger.pdict])
    # print('aVlog:')
    # print(aVlog)
    # print()
    # atlog = np.array([tlogger.pdict[n] for n in tlogger.pdict])
    # print('atlog')
    # print(atlog)
    # print()

    #plot vertical lines showing predictions from analytical theory
    #plt.vlines(516.1465,uis[0],237.80,colors='c')
    #plt.vlines(515.129,uis[0],59.61,colors='b')

    #Optimizations of scaling laws from Bush CCS database (using N2 data)

    # def plaw(x,a,p):
    #     return a*x**p
    #zvsmfit = plaw
    #popt, pcov = curve_fit(zvsmfit, Bushms, Bushzs)
    #plt.plot(Bushms,Bushzs,marker='o',linestyle='None')
    #massrange = np.linspace(0,1000,100)
    #plt.plot(massrange,zvsmfit(massrange,*popt))

    #CCSvsmfit = plaw
    #popt, pcov = curve_fit(CCSvsmfit, Bushms, BushN2CCSs)
    #plt.plot(Bushms,BushN2CCSs,marker='o',linestyle='None')
    #massrange = np.linspace(0,1000,100)
    #plt.plot(massrange,CCSvsmfit(massrange,*popt))

    # with open(r"C:\Users\jprell\Desktop\myo9Tvst.csv","w",newline="") as csvfile:
    #     writer = csv.writer(csvfile)
    #     for value in range(len(times)):
    #         writer.writerow([times[value],Tis[value]])
    # csvfile.close()

    if plotfile == 'screen':
        plt.show()
    else:
        plt.savefig(plotfile)

    return tlogger, postdictlist

def doplot(ax, tlogger, xaxis, yaxis, plabel, ptitle, plotfile):
    if xaxis == 't':
        xvals = tlogger.t
        pltxlabel = 't [s]'
    elif xaxis == 'coll':
        xvals = tlogger.coll
        pltxlabel = 'collisions'
    elif xaxis == 'zp':
        xvals = tlogger.zp
        pltxlabel = 'z position [m]'

    # add other possible y values options here
    if yaxis == 'Ti':
        yvals = tlogger.Ti
        pltylabel = 'T [K]'
    elif yaxis == 'coll':
        yvals = tlogger.coll
        pltylabel = 'collisions'
    elif yaxis == 'KE':
        yvals = tlogger.KE
        pltylabel = 'KE [eV]'
    elif yaxis == 'ui':
        yvals = tlogger.ui
        pltylabel = 'ui [eV]'
    elif yaxis == 'chi':
        yvals = tlogger.chi
        pltylabel = 'chi'
    elif yaxis == 'zp':
        yvals = tlogger.zp
        pltylabel = 'z position [m]'

    ax.plot(xvals, yvals, label=plabel)
    ax.legend()
    ax.set_title(ptitle)
    ax.set_xlabel(pltxlabel)
    ax.set_ylabel(pltylabel)

#    plt.suptitle(f'{"cion.name"}     option {"option"}')                # main title
#    plt.title('{repr(ccell)}', fontsize=8)     # becomes subtitle
#    plt.xlabel(pltxlabel)
#    plt.ylabel(pltylabel)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='python spa_avg_new.py',
        description='''Use direct methods to compute ion internal energy and kinetic energy.
        If no arguments given, will do a default computation. Arguments may be used to specify
        cell, ion and choice of computation method or outputs.
        xaxis may be t, coll, zp
        yaxis may be Ti, ui, KE, coll, chi, zp''')
    parser.add_argument('jsonfiles', nargs='*', help='a json input files to specify ion, cell and other parameters.')
    parser.add_argument('--injectV', nargs = '+', type=float, default=defV)    # injection V, cell.Voffset adds to this
    parser.add_argument('--option', default='IICT')    # choose which calculation option to do. TBD
    parser.add_argument('--sizefactor', type=float, default=1.0)
    parser.add_argument('--Ngfactor', type=float, default=1.0)
    parser.add_argument('--xaxis', default= 'coll')
    parser.add_argument('--yaxis', default = 'Ti')
    parser.add_argument('--plotfile', default='screen')


#    parser.add_argument('--dH', help='overrides the initial dH value from input file') #   , nargs=1)
#    parser.add_argument('--dS', help='overrides the initial dS value from input file') #   , nargs=1)
#    parser.add_argument('--xlabel', default = 't', help='x axis choice. May be t or z')
#    parser.add_argument('--ylabel', default = 'frac', help='y axis choice may be fracs, dfracs, fragmenteds, Temp')

    args = parser.parse_args()
#    print(args)

    paramd = loadjfile(args.jsonfiles)
    paramd['injectV'] = args.injectV
    paramd['option'] = args.option
    paramd['sizefactor'] = args.sizefactor
    paramd['Ngfactor'] = args.Ngfactor
    paramd['xaxis'] = args.xaxis
    paramd['yaxis'] = args.yaxis
    paramd['plotfile'] = args.plotfile

#    print(paramd)

    tlogvalues, postd = main(**paramd)

#    print(f'maximum deltaT: {maxdeltat:0.02f} K, maximum deltaU: {maxdeltau:0.01f} eV, max deltaU relative to init U: {maxreldeltau:0.01f}, sim drift velocity: {driftvel:0.03f} m/s')
#    print(f'MSvdrift: {MSvdrift:0.3f}   MSmodvdrift: {MSmodvdrift:0.3f}   EoverN: {EoverN:0.2f} Td')
#    print(f'maximum delta ui: {maxdelui:0.2f}   maxefficiency = {maxeff:0.4f}')
#    print(f'zmax = {zmax}')
#    print(f'umaxapp: {umaxapp}')
#    print(f'a={a}  b={b}  c={c}  h={h}  k={k}  l={l}  lamplus={lamplus}  lamminus={lamminus}  cplus={cplus}  cminus={cminus}')
    print(postd)
