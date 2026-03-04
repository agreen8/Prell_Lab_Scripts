# Copyright (c) 2024 James S. Prell, all rights reserved
from os import path

class const():
    '''Define useful constants for the ionspa scripts.
    Rbar, R, Av, mamu, kB, qe'''
    Rbar = 8.3144e-5        # pressure in bar
    R = 8.3145              # SI units  m3⋅Pa⋅K−1⋅mol−1
    Av = 6.02214076e23      # Avogadro number #/mol
#    mamu = 1.66053906660e-27   # recommended value
    mamu = 1.66053906717386e-27  # fudged by 3e-10 so mamug * Av = 1.0
    mamug = 1000*mamu       # amu unit mass in g
    kB = 1.380649e-23       # Boltzmann constant in J/K
    qe = 1.6021892e-19      # electron charge in Coulombs
    k = 8.978e9             # Coulomb constant in N m^2/C^2
    h = 6.62607015e-34      # Planck constant in J s
    kBoverh = kB/h          # convenient for fracloss calc
    ionspadir = path.dirname(__file__)