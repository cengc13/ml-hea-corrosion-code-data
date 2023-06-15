#!/usr/bin/env  python3

from ase import Atoms
from ase.visualize import view
from ase.io import read
from ase.calculators.emt import EMT
from ase.constraints import StrainFilter
from ase.optimize import BFGS as qn
from gpaw import GPAW, FermiDirac, PW

import numpy as np
from itertools import combinations

def get_calc(txt, atoms):
    calc = GPAW(txt='out_%s.txt'%(txt),
              mode=PW(ecut=350),
              xc='PBE',
              kpts=np.int32(30//atoms.cell.lengths()) + 1,
              occupations=FermiDirac(0.1),
              spinpol=True,
              maxiter=1000,)
    return calc

base_file = 'sqs_quinary.traj'

def lat_const():
#    calc = EMT()
    metal = read(base_file)
    for atom in metal:
        if atom.symbol in magmom_dict.keys():
            atom.magmom = magmom_dict[atom.symbol]
    file = 'FeCrCoNiAl'
    calc = get_calc(txt=f'{file}', atoms=metal)
    metal.calc = calc
    sf = StrainFilter(metal)
    dyn = qn(sf, trajectory=f'{file}.traj')
    try:
        dyn.run(fmax=0.05)
    except:
        pass

    # rattled cell
    metal.rattle(stdev=0.05)
    calc = get_calc(txt=f'{file}_R', atoms=metal)
    metal.calc = calc
    try:
        metal.get_potential_energy()
        metal.get_forces()
        metal.write(f'{file}_R.traj')
    except:
        pass

if __name__ == "__main__":

    # elementlist=['Fe', 'Co', 'Ni', 'Cr', 'Al', 'Mo']
    elementlist=['Fe', 'Co', 'Ni', 'Cr', 'Al']
    # elementlist=['Ni','Al', 'Pt', 'Au', 'Pd']
    magmom_dict={'Fe':2.0, 'Co':1.5, 'Ni':1.0,}

    lat_const()
