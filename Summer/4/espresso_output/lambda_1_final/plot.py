import matplotlib.pyplot as plt
import numpy as np
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/rod-energy_1.0.dat'
data= np.loadtxt(dir_file,unpack=False)

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

plt.plot(data[:,0], data[:,1])
plt.xlim((0.0, 2000.0))
plt.xlabel(r'Time $t$')
plt.ylabel(r'Coulomb energy')
plt.tight_layout()
plt.show()
