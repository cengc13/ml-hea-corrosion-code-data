#!/usr/bin/env  python3

from ase.lattice.compounds import L1_2
from ase.calculators.emt import EMT
from ase.constraints import StrainFilter
from ase.optimize import QuasiNewton as qn
# from ase.lattice.surface import fcc111, add_adsorbate
from gpaw import GPAW, FermiDirac, PW

import numpy as np
from itertools import combinations, permutations

def get_calc(txt, atoms):
    calc = GPAW(txt='out_%s.txt'%(txt),
              mode=PW(ecut=350),
              xc='PBE',
              kpts=np.int32(30//atoms.cell.lengths()) + 1,
              occupations=FermiDirac(0.1),
              spinpol=True,
              maxiter=1000,)
    return calc


def lat_const(ele_0, ele_1, a0):
    # calc=EMT()
    metal = L1_2(symbol=(ele_0, ele_1),
               latticeconstant=a0,
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
    try:
        metal.get_potential_energy()
        metal.get_forces()
        metal.write(f'{ele_0}_{ele_1}_D.traj')
    except:
        pass


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

    for ele_0, ele_1 in permutations(elementlist, 2):
        print(ele_0, ele_1)
        a0 = lc_dict[ele_0] * 2**0.5
        c0 = lc_dict[ele_1] * 2 ** 0.5
        a0 = c0/2 + a0/2
        lat_const(ele_0, ele_1, a0=a0)
