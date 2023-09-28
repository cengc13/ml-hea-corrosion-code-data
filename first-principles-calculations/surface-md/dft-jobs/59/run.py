import os
import json
import pickle
import numpy as np
from ase.io import read, Trajectory
from ase.parallel import paropen
from gpaw import GPAW, PW, FermiDirac, Mixer

with open("../gpaw_params.pkl", "rb") as f:
    gpaw_params = pickle.load(f)
gpaw_params.pop('txt', None)
def get_calc(text=None, atoms=None):
    kpts = [1, 1, 1]
    if atoms is not None:
        pbc = atoms.pbc
        cell = atoms.cell.lengths()
        for i, each_pbc in enumerate(pbc):
            kpts[i] = np.int32(30//cell[i]) if each_pbc else 1
    gpaw_params.update(dict(kpts=kpts))
    calc = GPAW(txt=text, **gpaw_params)
    return calc

index = int(os.path.split(os.getcwd())[-1])
traj = '../dft-images/%i.traj' % index
atoms = read(traj)

calc = get_calc(text='%s.txt' % str(index), atoms=atoms)
atoms.calc = calc

completed = True
results = {}
try:
    e = atoms.get_potential_energy()
    f = atoms.get_forces(apply_constraint=False)
    results.update(energy=e)
    results.update(forces=f)
except:
    completed = False

f = paropen('completed', 'w')
f.write(str(completed))
f.close()

with paropen('results', 'wb') as f:
    pickle.dump(results, f)
