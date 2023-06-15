import sys
from ase.io import read, write
from numpy.random import choice
from random import sample

filename = sys.argv[1]

_name = filename.split('.')[0]

atoms = read(filename, index=':')

print(len(atoms))

selected = sample(atoms, k=50)

write(f'{_name}_50.traj', images=selected)
