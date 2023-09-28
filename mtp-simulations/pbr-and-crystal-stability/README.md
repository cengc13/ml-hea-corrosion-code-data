
- b2.py: generate b2 structures and optimize them

- l12.py: generate l12 structures and optimize them

- l12_vs_b2.py: compare the cohesive energies of B2 versus L12 for each pair of Al and Cr compositoins of AlxCry(FeNiCo)_{100-x-y}  

- pbr_al_cr.py: calculate the pbr_cr for oxidation of Cr from the high-entropy alloys AlxCry(FeNiCo)_{100-x-y}

- The trajectory files for a structure with certain Cr and Al compositions are saved in the format [crystal structure]_al_[Al composition]_cr_[Cr composition].traj

Note: When the cell and atomic positions were simultaneously optimized, ASE gave some errors. Instead, those two are optimized sequentially.
We first perform optimization of atomic positions, then cell geometry and lastly another round of position adjustment.
