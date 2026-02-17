# -*- coding: utf-8 -*-
# hcprofiles.py
####################################
#
# hcprofiles
#
# Original work by Jim Prell, Sam Shepherd, Ken Newton, et.at. at University of Oregon
# Copyright (c) 2024 James S. Prell, all rights reserved
#
# class for loading, storing and access to heat capacity profiles
####################################

import numpy as np
from ionspa import const
import json
from matplotlib import pylab as plt
import sys
import os.path

# we need to be able to find the hcprofiles*.json files
fnloc = os.path.dirname(os.path.abspath(__file__))
#print(f'fnloc is {fnloc}, file is {__file__}')

class hcprofileclass():
    '''utility class to hold heat capacity curves.'''
    def __init__(self, filepath='hcprofiles2.json'):
        self.filepath = os.path.join(fnloc, filepath)
        self.proftype = None
        self.hcstruct = None    # holds the full hc data for all types
        self.Tsamples = None
        self.hcsamples = None
        self._readfile()
        self.UfromTsamples = None
    
    def __repr__(self):
        srep = str(dict(path=self.filepath, type=self.proftype))
        return srep

    def _readfile(self, profname='lipid'):
        with open(self.filepath, 'r') as f:
            self.hcstruct = json.load(f)
        self.proftype = profname
        if type(self.hcstruct[profname]) is dict:
#            print(f'hc new structure: {profname}')
            self.Tsamples = self.hcstruct[profname]['temperature']
            self.hcsamples = self.hcstruct[profname]['hc']
        else:
#            print(f'hc old structure: {profname}')
            self.hcsamples = self.hcstruct[profname]
            try:
                indx = 54
                print(f'deleting sample for T {self.Tsamples[indx]}')
                del self.hcsamples[indx]                      # remove the value for T = 298
            except:
                pass
            if len(self.hcsamples) < 100:
                self.Tsamples = list(np.arange(0,3001,50, dtype=float))
            else:
                self.Tsamples = list(np.arange(0,101,2, dtype=float))
                self.Tsamples.extend(np.arange(150,3001,50, dtype=float))

    def buildfunctions(self, profilename, dof):
        # self.Tsamples=np.arange(0,3001,50) #ADJUST TEMPERATURE RANGE AND STEP SIZE AS NEEDED TO MATCH LAWREN'S HEAT CAPACITY DATA
        if not self.hcstruct or profilename != self.proftype:
            self._readfile(profilename)
        # print(f'    lT = {len(self.Tsamples)}   lhc = {len(self.hcsamples)}')
        self.UfromTsamples = []
        Uint = 0
        self.UfromTsamples.append(Uint)

        for i in range(len(self.hcsamples)-1):
            Uint += dof/const.Av*(self.Tsamples[i+1]-self.Tsamples[i])*1/2*(self.hcsamples[i]+self.hcsamples[i+1])
            self.UfromTsamples.append(Uint)


    def TfromU(self,Ucurrent):
        T = np.interp(Ucurrent,self.UfromTsamples,self.Tsamples)
        return T

    def UfromT(self,temperature):
        U = np.interp(temperature,self.Tsamples,self.UfromTsamples)
        return U

def loadOldHCprof(fname):
    hcprof = None
    with open(fname, 'r') as f:
        hcprof = json.load(f)
    return hcprof

def loadNewHCprof():
    fname = 'hcprofiles2.json'
    hcprof = None
    with open(fname, 'r') as f:
        hcprof = json.load(f)
    return hcprof

def interpolateprof(mlo, hclo, mhi, hchi, m_interp):
    '''Given a low mass and matching hcprofile, 
    and a high mass and matching hcprofile,
    make a new profile for the mass m_interp.
    Return the new profile values.
    '''
    newhc = []
    for hlo, hhi in zip(hclo, hchi):
        sc = (m_interp-mlo) / (mhi-mlo)
        newh = hlo + (hhi-hlo) * sc
        newhc.append(newh)
    return newhc

def newtransformSave():
    hcold = hcprofileclass('hcprofiles_LRP.json')
    hckeys = ['lipid', 'oligonucleotide', 'drug', 'sugar', 'tunemix', 'tunemix322', 'tunemix622', 'tunemix922', 'tunemix1222', 'tunemix1522', 'tunemix2122', 'peptide']
    hcnew = dict()
    for key in hckeys:
        hcold._readfile(key)
        hcnew[key] = dict(temperature=hcold.Tsamples, hc=hcold.hcsamples)
        print(key, len(hcnew[key]['temperature']), len(hcnew[key]['hc']) )
    # Now interpolate the 622, 922, 1522 calibrant ions from the adjacent ones
    hcnew['tunemix622']['temperature'] = hcnew['tunemix322']['temperature']
    hcnew['tunemix622']['hc'] = interpolateprof(322.0, hcnew['tunemix322']['hc'], 1222.0, hcnew['tunemix1222']['hc'], 622.0)
    hcnew['tunemix922']['temperature'] = hcnew['tunemix322']['temperature']
    hcnew['tunemix922']['hc'] = interpolateprof(322.0, hcnew['tunemix322']['hc'], 1222.0, hcnew['tunemix1222']['hc'], 922.0)
    hcnew['tunemix1522']['temperature'] = hcnew['tunemix322']['temperature']
    hcnew['tunemix1522']['hc'] = interpolateprof(1222.0, hcnew['tunemix1222']['hc'], 2122.0, hcnew['tunemix2122']['hc'], 1522.0)

    jstr = json.dumps(hcnew, indent=3, separators=(', ', ': '))
    print(jstr)
    with open('hcprofiles2.json', 'w') as f:
        json.dump(hcnew, fp=f, indent=3, separators=(', ', ': '))


def transform(hcprofold):
    hckeys = ['lipid', 'oligonucleotide', 'drug', 'sugar', 'tunemix', 'tunemix322', 'tunemix622', 'tunemix922', 'tunemix1222', 'tunemix1522', 'tunemix2122', 'peptide']
    for key in hckeys:
        print(f'{key}: {len(hcprofold[key])} ')
    
    # Now interpolate the 622, 922, 1522 calibrant ions from the adjacent ones
    hcprofold['tunemix622'] = interpolateprof(322.0, hcprofold['tunemix322'], 1222.0, hcprofold['tunemix1222'], 622.0)
    hcprofold['tunemix922'] = interpolateprof(322.0, hcprofold['tunemix322'], 1222.0, hcprofold['tunemix1222'], 922.0)
    hcprofold['tunemix1522'] = interpolateprof(1222.0, hcprofold['tunemix1222'], 2122.0, hcprofold['tunemix2122'], 1522.0)
    hcnew = dict()
    for key, vlist in hcprofold.items():
        print(key,vlist)
        print()
#        Tsamples = None
#        hcnew[key] = dict(temperature=Tsamples, hc=vlist)
#    return hcnew

def update_to_new_LRP_values():
    '''Read old profile file, update to new format with temperature for each profile.'''
    fname = 'hcprofiles_LRP.json'
    hcprofold = loadOldHCprof(fname)
#    print(hcprofold)
#    print(json.dumps(hcprofold, separators=(', ', ':\n')))
#    print()
    hcnew = transform(hcprofold)
#    print(json.dumps(hcnew, separators=(', ', ':\n')))
#    print()
#    print(hcnew)
#    print()
    jstr = json.dumps(hcnew, indent=3, separators=(', ', ': '))
    print(jstr)
#    with open('hcprofiles2.json', 'w') as f:
#        json.dump(hcnew, fp=f, indent=3, separators=(', ', ': '))

#    fname = 'Full HC Database (IonSPA).xlsx'
#    df = pandas.read_excel(fname)
#    print(df)


if __name__ == '__main__':
#    newtransformSave()
#    sys.exit(0)

    hcold = hcprofileclass('hcprofiles2.json')
#    hcnew = hcprofileclass('hcprofiles2.json')
    hckeys = list(hcold.hcstruct.keys())
#    hckeys = ['lipid', 'oligonucleotide', 'drug', 'sugar', 'tunemix', 'tunemix322', 'tunemix622', 'tunemix922', 'tunemix1222', 'tunemix1522', 'tunemix2122', 'peptide']
    for key in hckeys:#['tunemix322']: #hckeys:#[0:1]:
        hcold._readfile(key)
        print(f'profile type: {hcold.proftype}')
        dof = 3 * 35 - 6        # tunemix322 has 37 atoms (not 13), others have other sizes
        dof = 33                # override with a bogus value for plotting
        hcold.buildfunctions(key, dof)
#        print(hcold.Tsamples)
#        print(len(hcold.Tsamples), len(hcold.hcsamples))
        ufromt = np.array(hcold.UfromTsamples) / 1.6021892e-19
        if len(hcold.hcsamples) > 0:
            plt.plot(hcold.Tsamples, hcold.hcsamples, label='Thc_'+key)   #, '-o')
            plt.plot(hcold.Tsamples, ufromt, label='Tu_'+key)
            plt.legend()
            plt.xlabel('T [K]')
            plt.ylabel('hc   or   U [eV]')

    plt.show()