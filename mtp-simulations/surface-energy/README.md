- make-mc-and-submit.py: generate parallel multiple jobs on slurm

- monte_carlo.py: implement MCMC algorithm

- mc_run.py: define IO and MCMC simulation parameters, and initiate MCMC runs

- compute_surface_energy.py: calculate surface energies

- count_ratios.py: calculate the degree of surface segregation

- initial 20x20x5 FCC111 surface structure is saved in the trajectory file `al_[Al composition]_cr_[Cr composition].traj`.

- MCMC simulations were carried in the folder named `[Al composition]-[Cr composition]`. Since the full ase-db file is too large, only the last struture of MC trajectory with 80000 attempts was saved.
