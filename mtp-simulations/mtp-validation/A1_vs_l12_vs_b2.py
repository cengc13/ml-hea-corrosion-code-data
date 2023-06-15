import numpy as np
from ase.io import read
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib
font = {'size': 18}
matplotlib.rc('font', **font)

elements = ['Fe', 'Ni', 'Cr', 'Co', 'Al']
atomic_es = [-4.868341787389, -3.436339185232, -5.9059268684, -4.433012581685, -2.065514396292]
e_dict = {ele: atomic_e for ele, atomic_e in zip(elements, atomic_es)}

cohs_fcc, cohs_l12, cohs_b2 = [], [], []
cohs_min = []
al_pers = np.arange(0, 41, step=5)
for al_per in al_pers:
    fcc =read(f'fcc_al_{al_per}.traj')
    symbols_dict = Counter(fcc.symbols)
    e_fcc = fcc.get_potential_energy()
    coh_fcc = e_fcc - sum(value*e_dict[key] for key, value in symbols_dict.items())
    coh_fcc /= len(fcc)
    cohs_fcc.append(coh_fcc)

    l12 =read(f'L12_al_{al_per}.traj')
    symbols_dict = Counter(l12.symbols)
    e_l12 = l12.get_potential_energy()
    coh_l12 = e_l12 - sum(value*e_dict[key] for key, value in symbols_dict.items())
    coh_l12 /= len(l12)
    cohs_l12.append(coh_l12)

    b2 =read(f'B2_al_{al_per}.traj')
    symbols_dict = Counter(b2.symbols)
    e_b2 = b2.get_potential_energy()
    coh_b2 = e_b2 - sum(value*e_dict[key] for key, value in symbols_dict.items())
    coh_b2 /= len(b2)
    cohs_b2.append(coh_b2)
    coh_min = min(coh_fcc, coh_l12, coh_b2)
    cohs_min.append(coh_min)

plt.figure(figsize=(7, 5))
plt.plot(al_pers, cohs_fcc, 's', color='b', label='FCC_A1')
plt.plot(al_pers, cohs_l12, 'o', color='r', label=r'L1$_2$')
plt.plot(al_pers, cohs_b2, '^', color='g', label='B2')
plt.plot(al_pers, cohs_min, '--', color='k')
plt.xlabel('Al composition [%]')
plt.ylabel('Cohesive energy [eV/atom]')
plt.xticks([0, 10, 20, 30, 40])
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig('A1_vs_L12_vs_B2.pdf')
plt.show()