#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import os

plt.rcParams['text.usetex'] = True

dir_path = os.path.dirname(os.path.realpath(__file__)) + '/time_data_verlet.txt'

data = np.genfromtxt(dir_path)

plt.plot(data[:,0], data[:,1], 'b', label=r'measured data')
plt.xlabel(r'skin$/\sigma$')
plt.ylabel(r'runtime $t/\mathrm{s}$')
plt.legend()
plt.show()
