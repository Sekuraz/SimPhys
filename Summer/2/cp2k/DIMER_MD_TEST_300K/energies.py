import numpy as np
import os
import matplotlib.pyplot as plt

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/water_dimer_md-1.ener'
data= np.loadtxt(dir_file,unpack=False)

print(np.mean(data[:,4]))
