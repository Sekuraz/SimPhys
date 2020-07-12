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


def analytical_profile(y, f, eta, d):
    ret = f * (y-0.5) * (d - (y-0.5)) / (2 * eta)
    if y < 0.5:
        ret = 0.0
    elif y > 30.5:
        ret = 0.0
    return ret

analytical_profile = np.vectorize(analytical_profile)

y = np.linspace(-0.5, 31.5, 1000)

for f in [0.001, 0.01, 0.1]:

    dir_file = os.path.dirname(os.path.realpath(__file__)) + '/velocity_' + str(f).replace('.','') + '.vtk'

    data = np.loadtxt(dir_file,unpack=False, skiprows=10)
    profile = data[0:1024:32,0]

    plt.plot(y, analytical_profile(y, f, 1.0, 30.0), label=r'Analytical result')
    plt.plot(profile, linestyle='none', marker='o', label=r'Lattice Boltzmann method')
    plt.legend()
    plt.xlim((-0.5, 31.5))
    plt.xlabel(r'$y$')
    plt.ylabel(r'$u_x(y)$')
    plt.tight_layout()
    plt.show()
