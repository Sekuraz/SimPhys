import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

def fit(t, a, b):
    return a + b * t

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/MSD_TIP3P.csv'
data= np.loadtxt(dir_file,unpack=False).T

#t = np.linspace(0.0,500.0,1000)

#param, cov = curve_fit(fit, np.log(data[0,1:3]), np.log(data[1,1:3]), maxfev = 1000000)
#print(param)
#print(cov)

#plt.plot(np.log(data[0,:]), np.log(data[1,:]))
#plt.plot(np.log(data[0,:]), fit(np.log(data[0,:]),*param), label='Fit')
plt.plot(data[0,:], data[1,:])
#plt.plot(t, 6 * 5.0608 * 10 ** (-3) * t)

#plt.xlim((-25.0,30.0))
#plt.ylim((0.0,1.05))
plt.xlabel(r'time $t$ in ps')
#plt.xlabel(r'$\mathrm{log}\left(\frac{t}{1\mathrm{ps}}\right)$')
plt.ylabel(r'MSD $\langle\Delta r^2\left(t\right)\rangle$ in nm$^2$')
#plt.ylabel(r'$\mathrm{log}\left(\frac{\langle\Delta r^2\left(t\right)\rangle}{1\mathrm{nm}^2}\right)$')
#plt.legend()
plt.tight_layout()
plt.show()
