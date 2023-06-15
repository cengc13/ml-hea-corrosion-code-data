#!/usr/bin/env python3
import numpy as np
from ase.build import bulk
from ase.io import Trajectory
from ase.optimize import BFGS
from ase.constraints import StrainFilter
from ase.calculators.emt import EMT
from gpaw import GPAW, FermiDirac, PW
from ase.visualize import view
from ase.spacegroup import crystal


def get_calc(txt, atoms):
    calc = GPAW(txt='out_%s.txt'%(txt),
              mode=PW(ecut=350),
              xc='PBE',
              kpts=np.int32(30//atoms.cell.lengths()) + 1,
              occupations=FermiDirac(0.1),
              spinpol=True,
              maxiter=333,)
    return calc

def lat_const(element='Pt', ):
    # calc = EMT()

    # primitive cell
    metal = bulk(name=element)
    # metal = crystal(element, [(0,0,0)], spacegroup=225,
    #                 cellpar=[a, a, a, 90, 90, 90])
    for atom in metal:
        if atom.symbol in magmom_dict.keys():
            atom.magmom = magmom_dict[atom.symbol]
    # print(reference_states[atomic_numbers[element]])
    calc = get_calc(txt=element, atoms=metal)
    metal.calc = calc
    sf = StrainFilter(metal)
    opt = BFGS(sf)
    traj = Trajectory('{}.traj'.format(element), 'w', metal)
    opt.attach(traj)
    opt.run(0.05)

    # cell dimension doubled
    metal = metal.repeat(3)
    metal.rattle(stdev=0.05)
    calc = get_calc(txt=element + 'double', atoms=metal)
    metal.calc = calc
    metal.get_potential_energy()
    metal.get_forces()
    metal.write('{}_double.traj'.format(element))


if __name__ == "__main__":
    elementlist=['Fe', 'Cr', 'Co', 'Ni', 'Al', 'Mo']
#    elementlist=['Pt']
    magmom_dict={'Fe':2.0, 'Co':1.5, 'Ni':1.0,}
    # latcnt_dict = {'Fe':3.467, 'Co':3.518, 'Ni':3.507, 'Cr': , 'Al': , 'Mo':}
    for element in elementlist:
        lat_const(element=element,)

