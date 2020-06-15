import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/gyration.csv'
N, rg2, rg2_lj = np.loadtxt(dir_file, delimiter =',', unpack=False).T

#Fit
def fit(x, a, b):
    return a + b * x
param, cov = curve_fit(fit, np.log(N), np.log(np.sqrt(rg2)), maxfev = 1000000)
print(param)

param2, cov2 = curve_fit(fit, np.log(N), np.log(np.sqrt(rg2_lj)), maxfev = 1000000)
print(param2)

plt.plot(np.log(N), np.log(np.sqrt(rg2)), 'k', linestyle='none', marker='o', label=r'ideal chain')
plt.plot(np.log(N), np.log(np.sqrt(rg2_lj)), 'r', linestyle='none', marker='o', label=r'real chain')
plt.plot(np.log(N), fit(np.log(N), *param), 'k')
plt.plot(np.log(N), fit(np.log(N), *param2), 'r')

plt.xlabel(r'$\mathrm{log}\, N$')
plt.ylabel(r'$\mathrm{log}\, \sqrt{\langle R_\mathrm{g}^2\rangle}$')
plt.legend()
plt.tight_layout()
plt.show()
