import numpy as np
from ase.eos import EquationOfState
from ase.units import kJ
import matplotlib.pyplot as plt

from ase.io import read, Trajectory
from mtp import MTP

elements = ['Fe', 'Ni', 'Cr', 'Co', 'Al']
mtp_file = '../fecrconial-mtp/20.mtp-run1/pot_aug.mtp'
calc = MTP(mtp_file, unique_elements=elements)

al_pers = np.arange(0, 41, step=5)
lcs = [3.44, 3.51, 3.62, 3.49, 4.08]
a0 = np.average(lcs[:-1])
a_Al = lcs[-1]

symbols = ['Cr', 'Ni', 'Co', 'Fe']
for al_per in al_pers:
    ag = read(f"fcc_al_{al_per}.traj")
    ag.calc = calc
    cell = ag.get_cell()
    traj = Trajectory(f'eos_al_{al_per}.traj', 'w')
    for x in np.linspace(0.92, 1.08, 10):
        ag.set_cell(cell * x, scale_atoms=True)
        ag.get_potential_energy()
        traj.write(ag)

    configs = read(f'eos_al_{al_per}.traj', index=':')
    # Extract volumes and energies:
    volumes = [ag.get_volume() for ag in configs]
    energies = [ag.get_potential_energy() for ag in configs]
    eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
    v0, e0, B = eos.fit()
    print('Al', al_per, round(B / kJ * 1.0e24, 2), 'GPa')
    print('Al', al_per, round(v0 ** (1/3)/5, 2), 'angstrom')
    eos.plot(f'al_{al_per}-eos.png')
    plt.clf()
