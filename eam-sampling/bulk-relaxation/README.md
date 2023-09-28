 - create_cfg.py: generate initial structure input for LAMMPS; generated input data is named as `data.fcc.fecrnicoal[al composition]`. 

 - opt_cell_and_positions.ipynb: structural optimization in LAMMPS. The relaxation files are saved in the dump-files folder.

 - hea.in: example lammps inputfile
 
 - read_dump_files.py: convert lammps relaxed output files to ASE trajectory files. The converted and selected structures for each Al composition (0% to 20% with an interval of 5%) are saved in the file `dump-files/al[Al composition].traj`. The combined file is named "al_xper.traj", which will be sent to DFT calculations.

 - plot_lc_and_ecoh.py: plot the lattice constants and cohesives energies
