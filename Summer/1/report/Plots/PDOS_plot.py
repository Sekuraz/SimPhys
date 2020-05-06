import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/diamond_pdos_2.pdos'
#dir_file = dir_path + '/smeared4.dat'
data= np.loadtxt(dir_file,unpack=False,delimiter=',').T

plt.plot(27.2*(data[1,:]-0.457204), data[3,:],'k', label='$s$-orbitals')
plt.plot(27.2*(data[1,:]-0.457204), data[4,:],'b',label='$p$-orbitals')
plt.plot(27.2*(data[1,:]-0.457204), data[5,:],'r',label='$d$-orbitals')
#plt.plot(data[0,:], data[1,:]/np.amax(data[1,:]))
#plt.plot(data[1,:], data[5,:],label='Simulation results')

plt.xlim((-25.0,30.0))
plt.ylim((0.0,1.05))
plt.xlabel(r'Energy in eV')
plt.ylabel(r'PDOS in arbitrary units')
plt.legend()
plt.tight_layout()
plt.show()
