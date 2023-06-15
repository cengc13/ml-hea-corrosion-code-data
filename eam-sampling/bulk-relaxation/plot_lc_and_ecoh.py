from matplotlib import pyplot as plt
import numpy as np

# Al percentage in the 500-atom model systems
al_ratios = 0.2 * np.linspace(0, 1, num=5, endpoint=True)

# lattice constants in Angstrom
lcs = [3.538, 3.560, 3.580, 3.602, 3.626]

# Cohesive energies in eV/atom
ecohs = [-4.363, -4.307, -4.247, -4.194, -4.141]

fig, ax = plt.subplots(2, 1, sharex=True, figsize=(5, 6))

ax[0].plot(100*al_ratios, lcs, '--o', label='lattice constant', color='brown')
ax[0].set_ylabel(r'Lattice constant [$\AA$]')
ax[1].plot(100*al_ratios, ecohs, '-.*', label='Cohesive energy', color='navy')
ax[1].set_ylabel(r'Cohesive energy [eV]')
ax[1].set_xlabel('Al content [%]')
ax[0].legend()
ax[1].legend()
plt.xlim((0, 20))
plt.tight_layout()
plt.show()
