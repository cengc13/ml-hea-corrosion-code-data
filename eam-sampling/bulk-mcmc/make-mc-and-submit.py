#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
from ase.io import read
import numpy as np

def make_run_and_submit(traj_path, job):
    fullpath = os.path.join('{:s}'.format(traj_path), '{:s}'.format(job))
    if not os.path.exists(fullpath):
        print(fullpath + ': Creating and submitting.')
        os.makedirs(fullpath)
    shutil.copy('mc_run.py', os.path.join(fullpath, 'mc_run.py'))
    pwd = os.getcwd()
    os.chdir(fullpath)
    subprocess.run(['sbatch', 'mc_run.py'])
    # subprocess.run(['python3', 'mc_run.py'])
    os.chdir(pwd)

traj_path = './'
al_ratios = ['25', '30', '35', '40']
for job in al_ratios:
    job = str(job)
    make_run_and_submit(traj_path, job)
