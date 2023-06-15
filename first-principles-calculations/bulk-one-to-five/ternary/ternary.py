#!/usr/bin/env  python3

from ase import Atoms
from ase.visualize import view
from ase.io import Trajectory, write, read
from ase.calculators.emt import EMT
from ase.constraints import StrainFilter
from ase.optimize import BFGS as qn
# from ase.lattice.surface import fcc111, add_adsorbate
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
              maxiter=333,)
    return calc

base_file = 'sqs_FeCoNi.traj'
_metal = read(base_file)
symbols = _metal.symbols
fe_indices = np.arange(len(_metal))[np.array(symbols) == 'Fe']
co_indices = np.arange(len(_metal))[np.array(symbols) == 'Co']
ni_indices = np.arange(len(_metal))[np.array(symbols) == 'Ni']
# print(fe_indices, co_indices, ni_indices)

def lat_const(ele_0, ele_1, ele_2):
    # calc = EMT()
    metal = read(base_file)
    for i in fe_indices:
        metal[i].symbol = metal[i].symbol.replace('Fe', ele_0)
    for i in co_indices:
        metal[i].symbol =  metal[i].symbol.replace('Co', ele_1)
    for i in ni_indices:
        metal[i].symbol =  metal[i].symbol.replace('Ni', ele_2)
    for atom in metal:
        if atom.symbol in magmom_dict.keys():
            atom.magmom = magmom_dict[atom.symbol]
    calc = get_calc(txt=f'{ele_0}_{ele_1}_{ele_2}', atoms=metal)
    metal.calc = calc
    sf = StrainFilter(metal)
    dyn = qn(sf, trajectory=f'{ele_0}_{ele_1}_{ele_2}.traj')
    try:
        dyn.run(fmax=0.05, steps=5)
    except:
        pass

    # rattled cell
    metal.rattle(stdev=0.05)
    calc = get_calc(txt=f'{ele_0}_{ele_1}_{ele_2}_R', atoms=metal)
    metal.calc = calc
    try:
        metal.get_potential_energy()
        metal.get_forces()
        metal.write(f'{ele_0}_{ele_1}_{ele_2}_R.traj')
    except:
        pass

if __name__ == "__main__":

    # elementlist=['Fe', 'Co', 'Ni', 'Cr', 'Al', 'Mo']
    elementlist=['Fe', 'Co', 'Ni', 'Cr', 'Al']
    # elementlist=['Ni','Al', 'Pt']
    magmom_dict={'Fe':2.0, 'Co':1.5, 'Ni':1.0,}

    for ele_0, ele_1, ele_2 in combinations(elementlist, 3):
        lat_const(ele_0, ele_1, ele_2)