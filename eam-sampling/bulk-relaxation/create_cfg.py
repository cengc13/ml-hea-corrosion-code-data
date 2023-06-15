from ase.calculators.eam import EAM
from ase.build import bulk
from ase.visualize import view
import numpy as np
import random
random.seed(827)
from random import shuffle
from ase.constraints import StrainFilter
from ase.optimize import BFGS
from ase.io import write
from ase.io.lammpsdata import write_lammps_data


## High-entropy alloys
atoms = bulk('Al', crystalstructure='fcc', cubic=True)

atoms = atoms.repeat(3)
n_atoms = len(atoms)

al_ratios = np.arange(0, 21, step=5) / 100

symbols = ['Fe', 'Cr', 'Ni', 'Co']

for al_ratio in al_ratios:
    n_al_atoms = int(al_ratio * n_atoms)
    print(n_al_atoms)
    for i in range(len(atoms)-n_al_atoms):
        atoms[i].symbol = symbols[i % 4]
    for i in range(len(atoms)-n_al_atoms, len(atoms)):
        atoms[i].symbol = 'Al'
    shuffle(atoms.symbols)
    write_lammps_data(f'data.fcc.fecrnicoal{al_ratio:.2f}', atoms, specorder=['Fe', 'Ni', 'Cr', 'Co', 'Al'])


