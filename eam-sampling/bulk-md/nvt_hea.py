#!/usr/bin/env python3
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --partition=short
#SBATCH --job-name=hea-fcc-vs-bcc
#SBATCH --mem=60GB
from ase import units
from ase.io import read
from ase.calculators.eam import EAM
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import (MaxwellBoltzmannDistribution,
                                         Stationary, ZeroRotation)

eam = EAM(potential='../FeNiCrCoAl-heaweight.eam.alloy')

T = 1700  # Kelvin

path = '../bulk-relaxation/dump-files'

al_ratios = ['0', '5', '10', '15', '20']

for al_ratio in al_ratios:
  atoms = read(f'{path}/al0{al_ratio}.traj')
  MaxwellBoltzmannDistribution(atoms, temperature_K=1200)
  atoms.calc = eam
  def printenergy(a=atoms):  # store a reference to atoms in the definition.
      """Function to print the potential, kinetic and total energy."""
      epot = a.get_potential_energy() / len(a)
      ekin = a.get_kinetic_energy() / len(a)
      print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
            'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
# We want to run MD with constant energy using the Langevin algorithm
# with a time step of 5 fs, the temperature T and the friction
# coefficient to 0.002 atomic units.
  dyn = Langevin(atoms, 5 * units.fs, temperature_K=T, friction=0.002)
  dyn.attach(printenergy, interval=50)

  # We also want to save the positions of all atoms after every 100th time step.
  traj = Trajectory(f'bulk_md_nvt_1700K_al{al_ratio}.traj', 'w', atoms)
  dyn.attach(traj.write, interval=50)

  # Now run the dynamics
  dyn.run(5000)