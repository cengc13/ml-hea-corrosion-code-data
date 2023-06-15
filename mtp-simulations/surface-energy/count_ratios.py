import numpy as np
from collections import Counter, OrderedDict
from ase.io import read
import matplotlib
matplotlib.rc('font', **{'size':16})
from matplotlib import pyplot as plt
import pandas as pd


al_pers = np.arange(0, 26, step=5)
cr_pers = np.arange(10, 31, step=5)

tags = [(1, 5), (2, 3, 4)]
surf_al_pers = []
surf_cr_pers = []
surf_fe_pers = []
surf_co_pers = []
surf_ni_pers = []
als_avg = []
crs_avg = []
als, crs = [], []

# cr_per = 30
for al_per in al_pers:
    for cr_per in cr_pers:
        als.append(al_per)
        crs.append(cr_per)
        file = f"al_{al_per}_cr_{cr_per}.traj"
        initial = read(file, index='0')
        final = read(f"{al_per}-{cr_per}/last.traj", index='-1')
        als_avg.append(al_per)
        avg = (100 - al_per - cr_per) / 3
        crs_avg.append(avg)
        for tag in tags:
            initial_syms = [atom.symbol for atom in initial if atom.tag in tag]
            final_syms = [atom.symbol for atom in final if atom.tag in tag]
            cm_initial = Counter(initial_syms)
            cm_final = Counter(final_syms)
            cm_diff = {key: cm_final[key] - cm_initial[key]
                       for key in ['Co', 'Cr', 'Fe', 'Ni', 'Al']}
            # print(cm_initial)
            # print(f"Al ratio {al_per}, layer {tag}:")
            # print(cr_per)
            surf_al_pers.append(cm_final['Al']/800 - al_per/100) if tag == (1, 5) else None
            surf_cr_pers.append(cm_final['Cr']/800 - cr_per/100) if tag == (1, 5) else None
            surf_fe_pers.append(cm_final['Fe']/800 - avg/100) if tag == (1, 5) else None
            surf_co_pers.append(cm_final['Co']/800 - avg/100) if tag == (1, 5) else None
            surf_ni_pers.append(cm_final['Ni']/800 - avg/100) if tag == (1, 5) else None
            # print(OrderedDict(sorted(cm_diff.items())))
data = {
    'al': surf_al_pers,
    'cr': surf_cr_pers,
    'fe': surf_fe_pers,
    'co': surf_co_pers,
    'ni': surf_ni_pers,
}

df = pd.DataFrame(data)

corr_matrix = df.corr()
print(corr_matrix)

X = np.array(als).reshape((len(al_pers), len(cr_pers)))
Y = np.array(crs).reshape((len(al_pers), len(cr_pers)))
Z = np.array(surf_cr_pers).reshape((len(al_pers), len(cr_pers)))
# Z = np.array(surf_ni_pers).reshape((len(al_pers), len(cr_pers)))
# Z = np.array(surf_fe_pers).reshape((len(al_pers), len(cr_pers)))
# Z = np.array(surf_co_pers).reshape((len(al_pers), len(cr_pers)))
# Z = np.array(surf_al_pers).reshape((len(al_pers), len(cr_pers)))
Z *= 100
fig,ax=plt.subplots(1,1, figsize=(5,4))
# cp = ax.imshow(Z, cmap='Reds')
# cp = ax.contourf(X, Y, Z, cmap='Blues')
cp = ax.contourf(X, Y, Z, cmap='Reds_r')
fig.colorbar(cp, label=r'$\Delta$x$_{\mathrm{Cr}}$ [%]',
# fig.colorbar(cp, label=r'$\Delta$ x$_{\mathrm{Al}}$ [%]',
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
plt.savefig('net_cr_gain.pdf')
# plt.savefig('net_ni_gain.png')
# plt.savefig('net_fe_gain.png')
# plt.savefig('net_co_gain.png')
# plt.savefig('net_al_gain.pdf')
plt.show()