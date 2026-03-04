# Copyright (c) 2024 James S. Prell, all rights reserved

VERSION = "1.0.1f"

def get_version():
    return VERSION

__all__ = ["get_version"]


import sys
from .spa_const import const
#from .cell import cellclass, WAVEcell, WAVEcell2, WAVEcell3, WAVEcell4
from .cell import *
from .cell_A_CIU import *
from .hcprofiles import hcprofileclass
from .spa_paramd import loadjfile, printparams


# since makecell is used within spa_heat, this is defined
# before the spa_heat imports
def makecell(celldict=None):
    '''Given a celldict with either a base "cell" type or a Waters cell "wcell" type, return a cell object.'''
    if celldict is None or celldict.get('type') is None:
        cell = cellclass()
    elif celldict['type'] == "cell":
        cell = cellclass(**celldict)
    elif celldict['type'] == "wcell":
        cell = WAVEcell(**celldict)
    elif celldict['type'] == "wcell2":
        cell = WAVEcell2(**celldict)
    elif celldict['type'] == "wcell3":
        cell = WAVEcell3(**celldict)
    elif celldict['type'] == "wcell4":
        cell = WAVEcell4(**celldict)
    elif celldict['type'] == "aCIUcell":
        cell = CIUcell(**celldict)
    elif celldict['type'] == "fullAcell":
        cell = fullAcell(**celldict)
    return cell


from .spa_heat import fracloss, fracremains, ionclass, lossrate
from .spa_heat import modelclass
#from .spa_heat import runWaters, testCID, testWAVE

# do this to prevent circular import references.
del sys.modules['ionspa.spa_heat']
# In practice, IonSPA.spa_heat will still be in sys.modules after the delete
del sys.modules['ionspa.hcprofiles']


__all__ = ["get_version", 
#            "spa_heat", 
#            "const",
#            "cellclass",
#            "WAVEcell",
#            "WAVEcell2",
#            "WAVEcell3",
#            "WAVEcell4",
#            "spa_paramd",
#            "cell",
#            "cell_A_CIU",
#            "hcprofiles",
#            "makecell"
            ]  # Make these available from the package
