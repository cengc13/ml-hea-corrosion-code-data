#!/usr/bin/env python3
import sys
from ase.io import read
import numpy as np

##### Write atoms objects to MTP configuration file

print("Usage: convert2cfg.py [input file] [Element string with , as delimiter] [output file]")
ind = ':'
if len(sys.argv) < 2:
    raise RuntimeError('ASE compatabile input file must be supplied.')
if len(sys.argv) < 3:
    raise RuntimeError('Element list is not provided.')
if len(sys.argv) < 4:
    output = 'train.cfg'
else:
    output = str(sys.argv[3])

filename = str(sys.argv[1])
ele_list = str(sys.argv[2]).split(',')
ele_dict = {ele.capitalize(): int(i) for i, ele in enumerate(ele_list)}
if len(sys.argv) == 5:
    ind = '-1:'
write_f, write_e = True, True

atoms = read(filename, index=ind)
f = open(output, 'w')
for image in atoms:
    try:
        e = image.get_potential_energy()
    except:
        write_e = False
    f.write('BEGIN_CFG\n')
    f.write(' Size\n')
    size = len(image)
    f.write(f'    {int(size)}\n')
    f.write(' Supercell\n')
    cell = image.get_cell()
    f.write("{0:<9}{1}      {2}      {3}\n".format(
        '', cell[0][0], cell[0][1], cell[0][2]))
    f.write("{0:<9}{1}      {2}      {3}\n".format(
        '', cell[1][0], cell[1][1], cell[1][2]))
    f.write("{0:<9}{1}      {2}      {3}\n".format(
        '', cell[2][0], cell[2][1], cell[2][2]))
    try:
        fs = image.get_forces()
    except:
        write_f = False
    if write_f:
        f.write(' AtomData:  id type       cartes_x      cartes_y'
                '      cartes_z           fx          fy          fz\n')
    else:
        f.write(' AtomData:  id type       cartes_x      cartes_y'
                '      cartes_z\n')
    pos = image.positions
    symbols = image.symbols
    for i in range(size):
        aid = int(i+1)
        atype = ele_dict[symbols[i]]
        x, y, z = pos[i]
        if write_f:
            f_x, f_y, f_z = fs[i]
            f.write('{0:>14}{1:>5}{2:>16.8f}{3:>16.8f}{4:>16.8f}{5:>12.6f}{6:>12.6f}{7:>12.6f}\n'.format(
                aid, atype, x, y, z, f_x, f_y, f_z))
        else:
            f.write('{0:>14}{1:>5}{2:>16.8f}{3:>16.8f}{4:>16.8f}\n'.format(aid, atype, x, y, z))
    if write_e:
        f.write(' Energy\n')
        f.write(f'{e:16.6f}\n')
    f.write('END_CFG\n')
    f.write('\n')
f.close()