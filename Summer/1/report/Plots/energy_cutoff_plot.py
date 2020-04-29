import matplotlib.pyplot as plt
import numpy as np
import os

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/energy_cutoff.csv'
e_cutoff, e_total = np.loadtxt(dir_file, delimiter =',', unpack=False).T

plt.plot(e_cutoff, e_total, linestyle='none', marker='o')

plt.xlim((0.0,500.0))
plt.xlabel(r'Cut-Off energy $E_\mathrm{cut-off}$ in Ry')
plt.ylabel(r'Total energy $E_\mathrm{tot}$ in Ry')
plt.tight_layout()
plt.show()
