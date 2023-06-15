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
from random import shuffle

from ase.optimize import BFGS
from ase.spacegroup import crystal
from ase.io import read, Trajectory
from ase.constraints import UnitCellFilter, StrainFilter

from mtp import MTP

elements = ['Fe', 'Ni', 'Cr', 'Co', 'Al']
mtp_file = '../../mtp-training/20.mtp/pot.mtp'
calc = MTP(mtp_file, unique_elements=elements)

al_pers = np.arange(0, 41, step=5)
lcs = [3.44, 3.51, 3.62, 3.49, 4.08]
a0 = np.average(lcs[:-1])
a_Al = lcs[-1]

symbols = ['Cr', 'Ni', 'Co', 'Fe']
for al_per in al_pers:
    ### FCC_A1 structures
    a = (a0 * (100 - al_per) + a_Al * al_per) / 100
    fcc = crystal('Pt', [(0,0,0)], spacegroup=225,
                  cellpar=[a, a, a, 90, 90, 90])
    fcc = fcc.repeat(rep=5)
    n_Al = int(al_per/100*len(fcc))
    for i in range(len(fcc)-n_Al):
        fcc[i].symbol = symbols[i % 4]
    for i in range(len(fcc)-n_Al, len(fcc)):
        fcc[i].symbol = 'Al'
    shuffle(fcc.symbols)
    fcc.calc = calc
    dyn = BFGS(fcc)
    dyn.run(0.05)
    ucf = StrainFilter(fcc, mask=[1,1,1,0,0,0])
    dyn = BFGS(ucf)
    dyn.run(0.05)
    fcc.set_constraint([])
    dyn = BFGS(fcc)
    dyn.run(0.05)
    fcc.write(f'fcc_al_{al_per}.traj')
