- make-mc-and-submit.py: generate parallel multiple jobs on slurm

- monte_carlo.py: implement MCMC algorithm

- mc_run.py: define IO and MCMC simulation parameters, and initiate MCMC runs

- random_select_50.py: select images from the mc\*.db files

- Each MCMC job and its results were saved in the folder named with Al compositions. The MCMC trajectory is in the ase-db format and the randomly selected 50 structures for each composition were saved in the trajectory file.
