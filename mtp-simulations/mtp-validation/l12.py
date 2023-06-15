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

al_pers = np.arange(0, 41, step=5)
lcs = [3.55 , 3.566, 3.586, 3.604, 3.621, 3.638, 3.655, 3.671, 3.687,  3.703, 3.727]
low, high = 0, len(al_pers)

symbols = ['Cr', 'Ni', 'Co', 'Fe']
for al_per, a in zip(al_pers[low:high], lcs[low:high]):
    ### L12 structures
    n_Al = int(al_per/100*500)
    L12 = L1_2(symbol=['Al', 'Pt'],
               latticeconstant=a, size=(5, 5, 5))
    indices_pt = [atom.index for atom in L12 if atom.symbol == 'Pt']
    indices_al = [atom.index for atom in L12 if atom.symbol == 'Al']
    if n_Al <= len(indices_al):
        to_replace = sample(indices_al, len(indices_al)-n_Al)
        if len(to_replace) < 1:
            indices = indices_pt
        else:
            indices = np.r_[to_replace, indices_pt]
        shuffle(indices)
        for i, ind in enumerate(indices):
            L12[ind].symbol = symbols[i % 4]
        # view(L12)
    else:
        for_al = sample(indices_pt, n_Al-len(indices_al))
        for_other = [_ for _ in indices_pt if _ not in for_al]
        for ind in for_al:
            L12[ind].symbol = 'Al'
        shuffle(for_other)
        for i, ind in enumerate(for_other):
            L12[ind].symbol = symbols[i % 4]
        # view(L12)
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
    L12.write(f'L12_al_{al_per}.traj')
