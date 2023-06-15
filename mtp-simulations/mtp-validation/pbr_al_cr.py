import numpy as np
from scipy.stats import linregress
import matplotlib
matplotlib.rc('font', **{'size':16})
from matplotlib import pyplot as plt

w_al2o3 = 102 # g/mol
rho_al2o3 = 3.987 # g/cm^3
v_al2o3 = w_al2o3 / rho_al2o3
w_cr2o3 = 152 # g/mol
rho_cr2o3 = 5.22 # g/cm^3
v_cr2o3 = w_cr2o3 / rho_cr2o3
a0, c_al, c_cr = [3.54987734e+00, 3.26672181e-03, 6.81928982e-05]
R = 6.02e23

v0_cr = 2/4 * R * (3.62 * 1e-8)**3
pbr_cr = v_cr2o3 / v0_cr
print(pbr_cr)

pbrs_low, pbrs_high = [], []
# al_pers = np.arange(1, 26, step=1)
al_pers = np.arange(0, 26, step=1)
cr_pers = np.arange(10, 31, step=1)
# al_pers = [25, 50, 75]
al_per = 20
als, crs = [], []
for al_per in al_pers:
    for cr_per in cr_pers:
        als.append(al_per)
        crs.append(cr_per)
        a = a0 + c_al * al_per + c_cr * cr_per
        a0_cr = a0 + c_al * al_per
        a0_al = a0 + c_cr * cr_per
        # print(a)
        v_tot = 2/cr_per * 25 * R * (a * 1e-8)**3
        # v_tot = 2/al_per * 25 * R * (a * 1e-8)**3
        v_left = 2/cr_per * (100-cr_per)/4 * R * (a0_cr * 1e-8) ** 3
        # v_left = 2/al_per * (100-al_per)/4 * R * (a0_al * 1e-8) ** 3
        # print(v_tot, v_left)
        v_al = v_tot - v_left
        v_al_high =  2/4 * R * (a * 1e-8)**3
        pbr_low = v_cr2o3 / v_al
        pbr_high = v_cr2o3 / v_al_high
        # print(al_per, pbr, v_al)
        # pbr_low = v_al2o3 / v_al
        # pbr_high = v_al2o3 / v_al_high

        pbrs_low.append(pbr_low)
        pbrs_high.append(pbr_high)

X = np.array(als).reshape((len(al_pers), len(cr_pers)))
Y = np.array(crs).reshape((len(al_pers), len(cr_pers)))
Z = np.array(pbrs_high).reshape((len(al_pers), len(cr_pers)))
# Z = np.array(np.array(Z).reshape((no_X, no_Y)))
fig,ax=plt.subplots(1,1, figsize=(5, 4))
cp = ax.contourf(X, Y, Z, cmap='Reds')
fig.colorbar(cp, label=r'PBR$_{\mathrm{Cr}}$',
# fig.colorbar(cp, label=r'PBR$_{\mathrm{Al}}$',
             location='right') # Add a colorbar to a plot
# ax.plot(size_pt_comm, 0, 's', markersize=10, color='k', label='Pt comm')
# ax.plot(test_size, test_x, 'h', markersize=10, color='k', label='Li et al.')
# ax.plot(8.9, 15, 'o', markersize=10, color='k', label='Li et al.')
# ax.plot(8, 46, '^', markersize=10, color='k', label='Xie et al.')
# ax.plot(X[ind_max_Z], Y[ind_max_Z], '*', markersize=10, color='k',
#         label='Optimum')

# ax.set_title(r'ORR mass activity')
ax.set_xticks(np.arange(0, 26, step=5))
ax.set_yticks(np.arange(10, 31, step=5))
ax.set_xlabel('Al composition [%]')
ax.set_ylabel('Cr composition [%]')
plt.tight_layout()
# plt.savefig('pbr_al.png')
# plt.savefig('pbr_cr.pdf')
plt.show()
