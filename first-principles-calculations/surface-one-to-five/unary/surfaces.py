#!/usr/bin/env  python3

from ase.build import fcc111
from ase.calculators.eam import EAM
from ase.optimize import BFGS as qn
from ase.constraints import FixAtoms
from gpaw import GPAW, FermiDirac, PW

import numpy as np
from itertools import combinations

def get_calc(txt, atoms):
    calc = GPAW(txt=f'out_{txt}.txt',
              mode=PW(ecut=350),
              xc='PBE',
              kpts=np.int32(30 // atoms.cell.lengths()),
              occupations=FermiDirac(0.1),
              spinpol=True,
              poissonsolver={'dipolelayer': 'xy'},
              symmetry={'point_group': False},
              maxiter=2000,)
    return calc

def run(elements, a):
    element = elements[0]
    metal = fcc111(symbol=element, size=(1, 1, 5), a=a, vacuum=10)
    for atom in metal:
        if atom.symbol in magmom_dict.keys():
            atom.magmom = magmom_dict[atom.symbol]
    formula = metal.get_chemical_formula()
#    calc = get_calc(txt=f'{formula}', atoms=metal)
    fixed_atoms = FixAtoms(indices=[atom.index for atom in metal if atom.tag == 5])
    metal.set_constraint(fixed_atoms)
    calc = EAM(potential='../../FeNiCrCoAl-heaweight.eam.alloy')
    metal.calc = calc
    dyn = qn(metal, trajectory=f'{formula}.traj')
    dyn.run(fmax=0.05)

if __name__ == "__main__":

    elementlist=['Fe', 'Co', 'Ni', 'Cr', 'Al']
    lc_dict = {
    'Fe': 2.439,
    'Co': 2.489,
    'Ni': 2.480,
    'Cr': 2.454,
    'Al': 2.862,
    'Mo': 2.722,
    }
    magmom_dict={'Fe':2.0, 'Co':1.5, 'Ni':1.0,}

    for eles in combinations(elementlist, 1):
        a = (2**0.5) * np.average([lc_dict[_] for _ in eles])
        run(eles, a=a)

