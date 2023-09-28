#!/usr/bin/env python3
#SBATCH --account=ap31-condo
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=batch
#SBATCH --mem=8g
#SBATCH --mail-type=end
#SBATCH --job-name=bulk-mcmc
#SBATCH --mail-user=cheng_zeng1@brown.edu
import numpy as np
import random
from ase.io import read
import os
from ase.visualize import view
from ase.calculators.eam import EAM

pwd = os.getcwd()
os.chdir('../../')
path = os.getcwd()
print(path)
import sys
sys.path.append(f'{path}/scripts/')
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
