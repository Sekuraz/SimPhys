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
dir_file = dir_path + '/lattice_constants.csv'
a, e_total = np.loadtxt(dir_file, delimiter =',', unpack=False).T

#Fit
def fit(x, a, b, c):
    return a + b * (x - c) ** 2
param, cov = curve_fit(fit, a**3, e_total, maxfev = 1000000)
print(param)
print(cov)

V=np.linspace(48.0,53.0, 2000)

plt.plot(V, fit(V, *param), label=r'Fit')
plt.plot(a**3, e_total, linestyle='none', marker='o', label='Simulation results')
#plt.plot(a, e_total, linestyle='none', marker='o', label='Simulation results')

plt.xlim((48.7,53.0))
plt.xlabel(r'Volume $a^3$ in \r{A}$^3$')
plt.ylabel(r'Total energy $E_\mathrm{tot}$ in Ry')
plt.legend()
plt.tight_layout()
plt.show()
