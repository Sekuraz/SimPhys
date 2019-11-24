#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import os

plt.rcParams['text.usetex'] = True

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/time_data_verlet_comparison.txt'
dir_path2 = os.path.dirname(os.path.realpath(__file__)) + '/time_data.txt'

data = np.genfromtxt(dir_path)
data2 = np.genfromtxt(dir_path2)

plt.plot(data[:,0], data[:,1], 'b', label=r'with verlet list')
plt.plot(data2[:,0], data2[:,1], 'r', label=r'without verlet list')
plt.xlabel(r'number of particles $N$')
plt.ylabel(r'runtime $t/\mathrm{s}$')
plt.legend()
plt.show()
