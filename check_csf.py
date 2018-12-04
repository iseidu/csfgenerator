#!/usr/bin/env python
"""
Script to check multiplicities and electron counts in a CSF list

Usage
-----
check_csf.py [csf filename]
"""
import sys
import numpy as np


def main():
    """The main routine."""
    csf_fname = sys.argv[1]

    # get number of internals and externals
    with open(csf_fname, 'r') as f:
        line = f.readline().split()[1:]
    for i, el in enumerate(line):
        if ':' in el:
            break
    n_int = i
    n_ext = (len(line) - i) // 2

    # generate column list
    rng = list(range(1, n_int+1)) + list(range(n_int+2, n_int+2+2*n_ext, 2))

    # import csf file
    csf = np.genfromtxt(csf_fname, usecols=rng)

    # find number of unpaired electrons
    n1 = np.sum(csf == 1, axis=1)
    n2 = np.sum(csf == 2, axis=1)

    # find total number of electrons
    n3 = 2*np.sum(csf == 3, axis=1)
    ne = n1 + n2 + n3

    # find total spins
    s = np.abs(n2 - n1)

    # output indices (lines number) of each electron count
    print('Electron count')
    print('--------------')
    clist = np.unique(ne)
    for ct in clist:
        lind = np.where(ne == ct)[0] + 1
        nl = len(lind)
        print('{:d} with {:d} electrons'.format(nl, ct))
        print(lind)
        print('\n')

    # output indices (line number) of each multiplicity
    print('Multiplicity')
    print('------------')
    mult = ['singlets', 'doublets', 'triplets', 'quartets']
    slist = np.unique(s)
    for st in slist:
        lind = np.where(s == st)[0] + 1
        nl = len(lind)
        print('{:d} {:s}'.format(nl, mult[st]))
        print(lind)
        print('\n')


if __name__ == '__main__':
    main()
