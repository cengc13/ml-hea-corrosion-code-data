import matplotlib.pyplot as plt
from ase.io import read
import numpy as np

al_pers = np.arange(0, 51, step=5)
lcs = []
for al_per in al_pers:
    # fcc = read(f"fcc_al_{al_per}.traj")
    # fcc = read(f"B2_al_{al_per}.traj")
    fcc = read(f"L12_al_{al_per}.traj")
    lengths = fcc.get_cell().lengths()
    lc = lengths.mean()
    lcs.append(lc/5)
    print(al_per, np.std(lengths))
print(repr(np.round(lcs,3)))
plt.plot(al_pers, lcs, '--*', color='brown')
# plt.savefig('fcc_lcs.png')
# plt.savefig('B2_lcs.png')
plt.savefig('L12_lcs.png')
# plt.show()