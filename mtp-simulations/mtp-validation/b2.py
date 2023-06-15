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
from ase.lattice.compounds import L1_0
from ase.visualize import view

from mtp import MTP

elements = ['Fe', 'Ni', 'Cr', 'Co', 'Al']
mtp_file = '../../mtp-training/20.mtp/pot.mtp'
calc = MTP(mtp_file, unique_elements=elements)

al_pers = np.arange(0, 41, step=5)
lcs = [3.55 , 3.566, 3.586, 3.604, 3.621, 3.638, 3.655, 3.671, 3.687,  3.703, 3.727]
# low, high = -3, -2
low, high = 0, len(al_pers)

symbols = ['Cr', 'Ni', 'Co', 'Fe']
for al_per, a in zip(al_pers[low:high], lcs[low:high]):
    ### B2 structures
    a *= (2/3)**0.5
    n_Al = int(al_per/100*(7**3)*2)
    b2 = crystal('Al', [(0,0,0)], spacegroup=229,
                  cellpar=[a, a, a, 90, 90, 90])
    for atom in b2:
        if atom.index % 2 == 0:
            atom.symbol = 'Pt'
    b2 = b2.repeat(7)
    indices_pt = [atom.index for atom in b2 if atom.symbol == 'Pt']
    indices_al = [atom.index for atom in b2 if atom.symbol == 'Al']
    to_replace = sample(indices_al, len(indices_al)-n_Al)
    if len(to_replace) < 1:
        indices = indices_pt
    else:
        indices = np.r_[to_replace, indices_pt]
    shuffle(indices)
    for i, ind in enumerate(indices):
        b2[ind].symbol = symbols[i % 4]
    # view(b2)
    b2.calc = calc
    dyn = BFGS(b2)
    dyn.run(0.05)
    # ucf = StrainFilter(b2, mask=[1,1,1,0,0,0])
    ucf = StrainFilter(b2, mask=[1,1,1,0,0,0])
    dyn = BFGS(ucf)
    dyn.run(0.05)
    b2.set_constraint([])
    dyn = BFGS(b2)
    dyn.run(0.05)
    b2.write(f'B2_al_{al_per}.traj')