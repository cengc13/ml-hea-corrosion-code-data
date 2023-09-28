#!/usr/bin/env python3
from testamp.nft.dft_jobs import relaxation_dft_jobs
from ase.io import read
from gpaw import GPAW, FermiDirac, PW, Mixer

def get_calc(text=None, atoms=None):
    calc = GPAW(txt=text,
                mode=PW(350),
                xc='PBE',
                kpts = (2, 2, 2),
                occupations=FermiDirac(0.1),
                spinpol=True,
                maxiter=666,
#                mixer=Mixer(0.02, 5, 100),
                convergence={'energy': 1e-4})
    return calc

calc = get_calc()

atoms = read('bulk-mcmc.traj', index=':')

relaxation_dft_jobs(atoms=atoms, parent_calc=calc,
               cores=64, memory='128GB')
