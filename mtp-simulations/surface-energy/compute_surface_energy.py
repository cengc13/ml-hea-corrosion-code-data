import numpy as np
from ase.io import read
import matplotlib
matplotlib.rc('font', **{'size':16})
from matplotlib import pyplot as plt
from math import exp


R = 8.314
T = 298
alpha = 0.5

# print(exp(alpha/R/T * (4000)))
al_pers = np.arange(0, 26, step=5)
cr_pers = np.arange(10, 31, step=5)

l12_path = '../pbr-and-crystal-stability'
toJpermpower2 = 1.6e-19 * 1e20
toJpermol = 1.6e-19 * 6.02e23
# print(toJpermol)
gammas = [] # surface energy
als, crs, gammas = [], [], []
for al_per in al_pers:
    for cr_per in cr_pers:
        als.append(al_per)
        crs.append(cr_per)
        file = f"{al_per}-{cr_per}/last.traj"
        final = read(file, index='-1')
        lengths = final.get_cell().lengths()
        # print(lengths)
        area = lengths[0] * lengths[1] * np.sin(60*np.pi/180)
        e_surf = final.get_potential_energy()
        bulk = read(f"{l12_path}/L12_al_{al_per}_cr_{cr_per}.traj")
        e_bulk = bulk.get_potential_energy()
        gamma = (e_surf - 4*e_bulk) / area / 2 * toJpermpower2
        # gamma = (e_surf - 4*e_bulk) / 400 / 2
        # gamma = (e_surf - e_bulk) / 2 /
        gammas.append(gamma)

X = np.array(als).reshape((len(al_pers), len(cr_pers)))
Y = np.array(crs).reshape((len(al_pers), len(cr_pers)))
Z = np.array(gammas).reshape((len(al_pers), len(cr_pers)))
fig,ax=plt.subplots(1,1, figsize=(5,4))
# cp = ax.imshow(Z, cmap='Reds')
cp = ax.contourf(X, Y, Z, cmap='jet')
fig.colorbar(cp, label=r'surface energy [J/m$^2$]',
# fig.colorbar(cp, label=r'surface energy [eV/atom]',
             location='right') # Add a colorbar to a plot
# ax.plot(8, 46, '^', markersize=10, color='k', label='Xie et al.')
# ax.text(5, 15, r'FCC(L1$_2$)', color='k', fontsize=24)
# ax.text(15, 26, r'BCC(B2)', color='w', fontsize=24)

ax.set_xticks(al_pers)
ax.set_yticks(cr_pers)
# ax.set_title(r'ORR mass activity')
ax.set_xlabel('Al composition [%]')
ax.set_ylabel('Cr composition [%]')
plt.tight_layout()
plt.savefig('surface_energy.pdf')
plt.show()