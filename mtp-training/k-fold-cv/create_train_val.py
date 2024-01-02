import os
import random
import numpy as np
from ase.io import read, write
from ase.visualize import view
from sklearn.model_selection import KFold

images = read('training.traj', index=':')
# view(images[0])
random.shuffle(images)
# view(images[0])

kf = KFold(n_splits=5, random_state=None, shuffle=False)

for i, (train_index, val_index) in enumerate(kf.split(images)):
	# print(len(train_index), len(val_index))
	# print(train_index, val_index)
	train = [images[_] for _ in train_index]
	val = [images[_] for _ in val_index]
	# view(train)
	write(f'train_fold{i}.traj', train)
	os.system(command=f'../training-data/convert2cfg.py '
		f'train_fold{i}.traj Al,Cr,Fe,Co,Ni train_fold{i}.cfg')
	write(f'val_fold{i}.traj', val)
	os.system(command=f'../training-data/convert2cfg.py '
		f'val_fold{i}.traj Al,Cr,Fe,Co,Ni val_fold{i}.cfg')
