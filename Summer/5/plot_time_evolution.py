import numpy as np
import matplotlib.pyplot as plt
import os

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

for f in [0.001, 0.01, 0.1]:

    dir_file = os.path.dirname(os.path.realpath(__file__)) + '/velocity_' + str(f).replace('.','') + '.dat'

    data = np.loadtxt(dir_file,unpack=False)

    plt.plot(data[:,0], data[:,1])
    plt.xlim((0.0, 1000.0))
    plt.ylim((0.0, f * 120))
    plt.xlabel(r'$t$')
    plt.ylabel(r'$u_x\left(y = d_y/2,t\right)$')
    plt.tight_layout()
    plt.show()
