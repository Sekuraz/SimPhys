import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

def infrared_absorption(nu, dt):
    t = dt * np.arange(len(data[:,0]))
    integral = 0.0
    for i in range(3):
        integrand = np.multiply(np.gradient(data[:,i], dt), np.exp(-1j*nu*t))
        integral += np.abs(np.trapz(integrand, dx = dt)) ** 2
    integral /= t[len(t)-1]
    return integral

infrared_absorption = np.vectorize(infrared_absorption)

dt = 0.5 #(fs)

#Set plot size
width = 5.787
height = width*0.6
plt.rc('figure', figsize=(width,height))

#Use LaTeX for fonts
plt.rc('font',**{'family':'serif','serif':['Computer Modern']})
plt.rc('text', usetex=True)

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_file = dir_path + '/dipole_dimer.txt'
data= np.loadtxt(dir_file,unpack=False)

f = np.linspace(0.0, 1.0, 1000)

plt.plot(5309 * f, infrared_absorption(f, dt), linewidth=0.5)
plt.ylabel(r'absorption cross-section $\alpha$ in arbitrary units')
plt.xlabel(r'wavenumber in $\mathrm{cm}^{-1}$')
plt.tight_layout()
plt.show()
