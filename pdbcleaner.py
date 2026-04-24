'''Input a PDB file for cleaning or rechaining'''

import argparse
parser = argparse.ArgumentParser(
    prog='python pdbcleaner.py',
    description='''Given a PDB file, execute one of the following
    clean: remove all non-ATOM information and any previous charge information,
    and resolve any Alt-Loc notation by keeping the "A" atoms and deleting the "B" atoms.
    rechain: automatically set the chain IDs to A-B-C-D-etc. in sequence''')
parser.add_argument('infile')
parser.add_argument('outfile')
parser.add_argument('-mode', choices=['clean','rechain'], default='clean') #"clean" will not rechain by default
parser.add_argument('-keep', choices=['A','B'], default='A') #can optionally keep the B Alt Loc residues and delete the A residues instead
args = parser.parse_args()

with open(f'{args.infile}', 'r') as infile, open(f'{args.outfile}', 'w') as outfile:
    if args.mode == 'clean':
        if args.keep == 'A':
            idx = 1
            for line in infile:
                if line[:4] == 'ATOM':
                    if line[16] == 'A':
                        outfile.write(line[:11-len(str(idx))] + str(idx) + line[11:16] + ' ' + line[17:78] + '\n')
                        idx += 1
                    elif line[16] == 'B':
                        continue
                    else: 
                        outfile.write(line[:11-len(str(idx))] + str(idx) + line[11:78] + '\n')
                        idx += 1
        if args.keep == 'B':
            idx = 1
            for line in infile:
                if line[:4] == 'ATOM':
                    if line[16] == 'B':
                        outfile.write(line[:11-len(str(idx))] + str(idx) + line[11:16] + ' ' + line[17:78] + '\n')
                        idx += 1
                    elif line[16] == 'A':
                        continue
                    else: 
                        outfile.write(line[:11-len(str(idx))] + str(idx) + line[11:78] + '\n')
                        idx += 1

    #The "clean" loop looks at each line of the infile that starts with “ATOM” and copies it
    #to the outfile, trimming off anything after column 65 (i.e.,
    #it trims off any charges that were in the input file.

    #It also deletes any "B" Alt Loc atoms and removes the "A" from the others

    #It also manually increments the atom number by 1 for each out-written line
    #This accounts for the loss of atom numbers due to the deletion of the "B" Alt Loc atoms

    if args.mode == 'rechain':
        chain_ID = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'] #Max: 16 chains
        abc_idx = -1
        for line in infile:
            if line[:4] == 'ATOM':
                new_chain = line[:21] + chain_ID[abc_idx] + line[22:] #replace the chain ID character with new chain ID
                outfile.write(new_chain)
                if line[22:26] == '   1' and line[13:16] == 'N  ':
                    abc_idx += 1 #increment to next index of the chain ID list when encountering the N terminus of a new chain
            else:
                outfile.write(line) #preserves any lines that are not ATOM lines

    #This loop automatically updates the chain ID of a PDB file to A-B-C in sequence
    #Chain ID increments based on the N-terminus of each chain