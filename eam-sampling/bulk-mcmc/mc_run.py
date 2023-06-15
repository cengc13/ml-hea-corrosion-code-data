#!/usr/bin/env python3
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --partition=short
#SBATCH --job-name=mc-run
#SBATCH --mem=60GB
import numpy as np
import random
from ase.io import read
import os
from ase.visualize import view
from ase.calculators.eam import EAM

pwd = os.getcwd()
os.chdir('../')
path = os.getcwd()
# print(path)
import sys
sys.path.append(f'{path}/')
from monte_carlo import main as CEMC
os.chdir(pwd)

eam_file = f'{path}/FeNiCrCoAl-heaweight.eam.alloy'
calc = EAM(potential=eam_file)
al_ratio = str(os.path.split(os.getcwd())[-1])
atoms = read(f'{path}/bulk-relaxation/dump-files/al0{al_ratio}.traj')
random.seed(527)
random.shuffle(atoms.numbers)
# view(atoms)
atoms.calc = calc

db_file = f'mc_al{al_ratio}.db'
saved_db_file = f'mc_al{al_ratio}-2.db'
if os.path.exists(db_file):
    atoms = read(db_file)
    os.system(f'mv {db_file} {saved_db_file}')
ratio = CEMC(atoms, dbname=db_file, T=500, steps=20000)
print(ratio)
