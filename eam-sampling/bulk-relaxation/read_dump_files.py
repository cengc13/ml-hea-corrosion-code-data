from ase.io.cfg import read_cfg
from ase.io import write
from ase.visualize import view
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np


def enforce_magnetic_moments(atoms, supplied_magmom_dict=None):
    magmom_dict={'Fe':2.0, 'Co':2.1, 'Ni':1.5, 'Ru':2.0, 'Rh':1.8}
    if supplied_magmom_dict is not None:
        magmom_dict.update(supplied_magmom_dict)
    for atom in atoms:
        if atom.symbol in magmom_dict.keys():
            atom.magmom = magmom_dict[atom.symbol]
    return atoms

def get_cfg_with_e_and_f(atoms):
    e = atoms.arrays['c_eng'].sum()
    fs = np.transpose([atoms.arrays['fx'], atoms.arrays['fy'], atoms.arrays['fz']])
    # print(fs.shape)

    _atoms = Atoms(numbers=atoms.arrays['numbers'],
                   positions=atoms.arrays['positions'],
                   pbc=atoms._pbc,
                   cell=atoms._cellobj
                   )

    calc = SinglePointCalculator(_atoms, energy=e, forces=fs)
    _atoms.calc = calc
    return _atoms

al_ratios = ['0', '5', '10', '15', '20']
n_opt_steps = [497, 320, 408, 401, 400] # You need to find the indices of relaxed dump.*.cfg files.
for al_ratio, n_opt_step in zip(al_ratios, n_opt_steps):
    all_atoms = []
    num = n_opt_step//10
    steps = 10*np.arange(0, num+1)
    steps = np.r_[steps, n_opt_step]
    for step in steps:
        atoms = read_cfg(f'dump-files/dump.al0{al_ratio}_{int(step)}.cfg')
        _atoms = get_cfg_with_e_and_f(atoms)
        _atoms = enforce_magnetic_moments(_atoms)
        all_atoms.append(_atoms)
    write(f'dump-files/al0{al_ratio}.traj', images=all_atoms)


# view(atoms)
