#!/usr/bin/env python3
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --partition=short
#SBATCH --job-name=hea-fcc-vs-bcc
#SBATCH --mem=60GB
import pickle
import numpy as np
import random
random.seed(827)
from random import shuffle, sample

from ase.optimize import BFGS
from ase.spacegroup import crystal
from ase.io import read, Trajectory
from ase.constraints import UnitCellFilter, StrainFilter
from ase.lattice.compounds import L1_2
from ase.visualize import view

from mtp import MTP

elements = ['Fe', 'Ni', 'Cr', 'Co', 'Al']
mtp_file = '../../mtp-training/20.mtp/pot.mtp'
calc = MTP(mtp_file, unique_elements=elements)

a0, c_al, c_cr = [3.54987734e+00, 3.26672181e-03, 6.81928982e-05]
al_pers = np.arange(0, 26, step=5)
cr_pers = np.arange(10, 31, step=5)

symbols = ['Ni', 'Co', 'Fe']
for al_per in al_pers:
    for cr_per in cr_pers:
        ### B2 structures
        a = a0 + c_al * al_per + c_cr * cr_per
        n_Al = int(al_per/100*500)
        n_Cr = int(cr_per/100*500)
        L12 = L1_2(symbol=['Al', 'Pt'],
                   latticeconstant=a, size=(5, 5, 5))
        indices_pt = [atom.index for atom in L12 if atom.symbol == 'Pt']
        indices_al = [atom.index for atom in L12 if atom.symbol == 'Al']
        to_replace = sample(indices_al, len(indices_al)-n_Al)
        indices = np.r_[to_replace, indices_pt]
        for_cr = sample(list(indices), n_Cr)
        for i in for_cr:
            L12[int(i)].symbol = 'Cr'
        _indices = [int(_) for _ in indices if _ not in for_cr]
        shuffle(_indices)
        for i, ind in enumerate(_indices):
            L12[ind].symbol = symbols[i % 3]

        L12.calc = calc
        # L12.set_constraint([])
        dyn = BFGS(L12)
        dyn.run(0.05)
        ucf = StrainFilter(L12, mask=[1,1,1,0,0,0])
        dyn = BFGS(ucf)
        dyn.run(0.05)
        L12.set_constraint([])
        dyn = BFGS(L12)
        dyn.run(0.05)
        L12.write(f'L12_al_{al_per}_cr_{cr_per}.traj')
