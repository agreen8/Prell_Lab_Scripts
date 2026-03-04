####################################
#
# spa_paramd.py
# Copyright (c) 2024 James S. Prell, all rights reserved
#
# Utility function to get all parameters from json file 
# The json file name is usually specified on the command line.
# Contents include sections to define the cell, ion and other run parameters
# This may be used by any of the main IonSPA programs or any other programs
# that wish to use the same definitions
#
####################################

import json
from pathlib import Path, PureWindowsPath


def loadjfile(jsonfile):
    '''Load the parameters from the jsonfile.
    Note that jsonfile is a list of file names.
    Each filename after the first (if included) will overwrite any values
    specified in the previous json files. This allows the option of specifying 
    cell or ion params in separate files without using the cellfilename or ionfilename
    parameters. 
    '''
    paramd = dict(jsonfile=jsonfile)

    for jfile in jsonfile:
        cellfname = None
        ionfname = None
        with open(jfile, 'r') as f:
            jdict = json.load(f)
            paramd.update(jdict)
#            print(f'loaded/updated: {jdict}')
            cellfname = jdict.get("cellfilename")
            if cellfname:
                if '\\' in cellfname:
                    from pathlib import Path, PureWindowsPath
                    cellfname = Path(PureWindowsPath(cellfname))
                    paramd['cellfilename'] = cellfname
                with open(cellfname, 'r') as fc:
                    fcdict = json.load(fc)
                    paramd.update(fcdict)
#                print('cell updated')
            ionfname = jdict.get("ionfilename")
            if ionfname:
                if '\\' in ionfname:
                    from pathlib import Path, PureWindowsPath
                    ionfname = Path(PureWindowsPath(ionfname))
                    paramd['ionfilename'] = ionfname
                with open(ionfname, 'r') as fi:
                    fidict = json.load(fi)
                    paramd.update(fidict)
#                print('ion updated')
#        print(f'from {jfile}: {paramd}\n')

    return paramd

#    print('printing paramd')
#    print(paramd)
#    print()

def printparams(paramd):
    print('paramd = \n{')
    for dkey,dvalue in paramd.items():
        if type(dvalue) is dict:
            print(f"    '{dkey}': " + "{")
#            print(dvalue)
            for skey,svalue in dvalue.items():
                print(f"        '{skey}': {svalue}")
            print('    }')
        else:
            print(f"    '{dkey}': {dvalue}")
    print('}')
 

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        prog='python spa_paramd.py',
        description='''Given a .json file (as from spa_mpirun), read the file and 
        referenced files and display the dict structure returned from the loadjfile() function.
        ''')
    parser.add_argument('jsonfile', nargs='+', help='json input files to specify ion, cell and other parameters.')
    args = parser.parse_args()
#    print(args)
#    print()

    jsonfile = args.jsonfile
    paramd = loadjfile(jsonfile)
    printparams(paramd)
    

"""
Example paramd object after loading.
The cell and ion dicts were loaded from cellfilename and ionfilename, respectively.

{
    'cellfilename': 'IonSPA_1.0_inputs/Synapt_trap3_high.txt', 
    'ionfilename': 'IonSPA_1.0_inputs/Stx_13_ion.txt', 
    'exp': {'folder': 'IonSPA_1.0_expdata', 'file': 'Stx_13_high_C.csv'}, 
    'tmax': 0.002, 
    'plotvf': 1, 
    'plottT': 0, 
    'Nrepeats': 50, 
    'dhds_start': [196, 184], 
    
    'cell': {
        'type': 'wcell3', 
        'L': 0.51, 
        'wl': 0.0121, 
        'wh': 4.0, 
        'vw': 311.0, 
        'Mgamu': 40, 
        'T': 298.0, 
        'rho': 5.9895e-05, 
        'rho1': 1.6063e-06, 
        'zopen': 0.09, 
        'zclose': 0.42, 
        'zV1': 0.045, 
        'dV1': -1, 
        'zV2': 0.09, 
        'dV2': -2, 
        'zV3': 0.13, 
        'dV3': -1, 
        'zV4': 0.38, 
        'dV4': -1, 
        'zV5': 0.42, 
        'dV5': -4, 
        'zwaveoff': 0.13, 
        'zwaveon': 0.38, 
        'z0': 0.0, 
        'T0offset': 0.0},
    'ion': {
        'name': 'Stx13', 
        'mass': 39.111, 
        'charge': 13, 
        'CCS': 31.4, 
        'num_atoms': 5495, 
        'T0': 300, 'hcprofname': 
        'peptide', 
        'refdH': 196, 
        'refdS': 184}
}
 """