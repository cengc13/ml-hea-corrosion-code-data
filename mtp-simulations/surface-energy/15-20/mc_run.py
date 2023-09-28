#!/usr/bin/env python3
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --partition=short
#SBATCH --job-name=hea-surface-mcmc
#SBATCH --mem=60GB
import numpy as np
import random
from ase.io import read
import os
from ase.visualize import view
from mtp import MTP
import sys

path = '../'
sys.path.append(f'{path}')
from monte_carlo import main as CEMC

mtp_file = f'{path}/pot_aug.mtp'
elements = ['Fe', 'Ni', 'Cr', 'Co', 'Al']
calc = MTP(mtp_file, unique_elements=elements)
alcr = str(os.path.split(os.getcwd())[-1])
print(alcr)
al_per, cr_per = alcr.split('-')
atoms = read(f'{path}/al_{al_per}_cr_{cr_per}.traj')
#random.seed(527)
#random.shuffle(atoms.numbers)
# view(atoms)
atoms.calc = calc

db_file = f'mc_al_{al_per}_cr_{cr_per}.db'
saved_db_file = f'mc_al_{al_per}_cr_{cr_per}-2.db'
if os.path.exists(db_file):
    atoms = read(db_file)
    os.system(f'mv {db_file} {saved_db_file}')
ratio = CEMC(atoms, dbname=db_file, T=500, steps=80000)
print(ratio)
