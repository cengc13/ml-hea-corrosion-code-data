#!/usr/bin/env  python3

from ase.lattice.compounds import L1_0
from ase.calculators.emt import EMT
from ase.constraints import StrainFilter
from ase.optimize import BFGS as qn
from gpaw import GPAW, FermiDirac, PW, Mixer

import numpy as np
from itertools import combinations

def get_calc(txt, atoms):
    calc = GPAW(txt='out_%s.txt'%(txt),
              mode=PW(ecut=350),
              xc='PBE',
              kpts=np.int32(30 // atoms.cell.lengths()) + 1,
              occupations=FermiDirac(0.1),
              spinpol=True,
              maxiter=1000,)
    return calc


def lat_const(ele_0, ele_1, a0, c0):
    # calc=EMT()
    metal = L1_0(symbol = (ele_0, ele_1),
               latticeconstant={'a': a0, 'c': c0,},
               size=(1, 1, 1))
    for atom in metal:
        if atom.symbol in magmom_dict.keys():
            atom.magmom = magmom_dict[atom.symbol]
    calc = get_calc(txt=f'{ele_0}-{ele_1}', atoms=metal)
    metal.calc = calc
    sf = StrainFilter(metal)
    dyn = qn(sf, trajectory=f'{ele_0}_{ele_1}.traj')
    try:
        dyn.run(fmax=0.05, steps=5)
    except:
        pass

    # cell dimension doubled
    metal = metal.repeat(2)
    metal.rattle(stdev=0.05)
    calc = get_calc(txt=f'{ele_0}-{ele_1}-D', atoms=metal)
    metal.calc = calc
    metal.get_potential_energy()
    metal.get_forces()
    metal.write(f'{ele_0}_{ele_1}_D.traj')

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

    for ele_0, ele_1 in combinations(elementlist, 2):
        a0 = lc_dict[ele_0] * 2**0.5
        c0 = lc_dict[ele_1] * 2 ** 0.5
        c0 = c0/2 + a0/2
        lat_const(ele_0, ele_1, a0=a0, c0=c0)

