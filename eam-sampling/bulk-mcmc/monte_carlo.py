#!/usr/bin/env python3
import numpy as np
import random
from ase.units import kB
from ase.db import connect
from ase.neighborlist import NeighborList
from ase.calculators.singlepoint import SinglePointCalculator as SPC
from ase.calculators.eam import EAM
from ase.optimize import MDMin
from amp.utilities import Logger

import os

def main(atoms, dbname, T=500, steps=20000):
    log = Logger('mc.log', overwrite=True)
    pwd = os.getcwd()
    os.chdir('../../')
    path = os.getcwd()
    os.chdir(pwd)
    eam_file = f'{path}/FeNiCrCoAl-heaweight.eam.alloy'
    db = connect(dbname)
    # Setting up variables for grand canonical MC
    symbols = atoms.get_chemical_symbols()
    sym = list(set(symbols))
    chem_bins = {_: [] for _ in sym}

    for i, s in enumerate(symbols):
        chem_bins[s] += [i]
    # Ensure sym1 has the lower concentration
    if len(chem_bins[sym[0]]) < len(chem_bins[sym[1]]):
        sym.reverse()
    # Calculate the initial energy and store it
    atoms.calc = atoms.get_calculator()
    nrg = atoms.get_potential_energy()
    log('MC initialized...')
    log(f'Initial energy: {nrg:.2f}.')
    # Write the initial configuration
    if db.count() == 0:
        dummy = atoms.copy()
        dummy.set_calculator(SPC(atoms, energy=nrg))
        db.write(dummy)
        # Construct a Neighbors list
        r = 1.4
        nl = NeighborList([r]*len(atoms),
                          self_interaction=False,
                          bothways=True)
        # Perform MC steps
        attempt, success = 0, 0
        while attempt < steps:
            ind1 = None
            log(f'Attempt {attempt}')
            while ind1 is None:
                # First, choose a random index from sym[0]
                random.shuffle(chem_bins[sym[0]])
                ind0 = chem_bins[sym[0]][-1]
                # Calculate nearest neighbors
                nl.update(atoms)
                indices, _ = nl.get_neighbors(ind0)
                # Determine if sym2 neighbors exist and choose one
                sym1_neighbors = [i for i in indices if atoms[i].symbol == sym[1]]
                if sym1_neighbors:
                    ind1 = random.sample(sym1_neighbors, 1)[0]
                else:
                    continue

                # Create new atoms object to test
                new_atoms = atoms.copy()
                new_atoms[ind0].symbol, new_atoms[ind1].symbol = sym[1], sym[0]
                new_atoms.calc = EAM(potential=eam_file)
                new_nrg = new_atoms.get_potential_energy()
                log(f'Energy: {new_nrg:.2f}')
                # Determine if lower than previous energy
                if new_nrg < nrg:
                    log(f'Attempt {attempt} accepted due to lower energy.')
                    atoms = new_atoms
                    nrg = new_nrg
                    chem_bins[sym[1]][-1] = ind0
                    chem_bins[sym[0]][-1] = ind1

                    dummy = atoms.copy()
                    dummy.set_calculator(SPC(atoms, energy=nrg))
                    db.write(dummy)
                    success += 1
                elif np.exp(-(new_nrg - nrg) / (kB * T)) > np.random.rand():
                    log(f'Attempt {attempt} accepted due to probabilistic seleciton.')
                    atoms = new_atoms
                    nrg = new_nrg
                    chem_bins[sym[1]][-1] = ind0
                    chem_bins[sym[0]][-1] = ind1
                    dummy = atoms.copy()
                    dummy.set_calculator(SPC(atoms, energy=nrg))
                    db.write(dummy)
                    success += 1
                else:
                    log('Not accepted, restart from previous configuration.')
                attempt += 1
        return success/attempt
