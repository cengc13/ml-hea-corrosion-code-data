#!/usr/bin/env python3
import os, pickle
import numpy as np
from ase.io import Trajectory, read
from ase.calculators.singlepoint import SinglePointCalculator
from collections import Counter

traj = Trajectory('hea_sp_bulk_al_xper.traj', mode='w')
atoms = read('al_xper.traj', index=':')
indices = range(len(atoms))

pwd = os.getcwd()
dftjobpath = "dft-jobs"
os.chdir(dftjobpath)
fulldftpath = os.getcwd()
for index in indices:
    os.chdir('%i' % index)
    f = open('completed')
    completed = str(f.read())
    f.close()
    if completed == 'False':
        print(index)
        os.chdir(fulldftpath)
        continue
    atoms_  = atoms[index]
    with open("results", "rb") as f:
        results = pickle.load(f)
    energy = results['energy']
    forces = results['forces']
    sp = SinglePointCalculator(atoms_, energy=energy,
                               forces=forces)
    atoms_.set_calculator(sp)
    traj.write(atoms_)
    os.chdir(fulldftpath)
os.chdir(pwd)
